import streamlit as st
import pandas as pd
import os
import tempfile
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from io import StringIO
# 共通関数
from modules.common import find_tool_path, generate_alignment_html_from_df, run_command

def app_phylo():
    st.header("🌳 PhyloPipeline Pro")
    st.info("MAFFT → (trimAl) → 編集 → IQ-TREE")

    # パス設定 (デフォルト値を入れておく)
    mafft_def = find_tool_path("mafft") or "mafft"
    trimal_def = find_tool_path("trimal") or "trimal"
    iqtree_def = find_tool_path("iqtree") or find_tool_path("iqtree2") or "iqtree2"
    
    with st.expander("🔧 ツール詳細設定", expanded=False):
        t1, t2, t3 = st.tabs(["MAFFT", "trimAl", "IQ-TREE"])
        with t1:
            mafft_binary = st.text_input("MAFFT Path", value=mafft_def)
            mafft_algo = st.selectbox("Algorithm", ["--auto", "--linsi", "--fftnsi"])
            c1, c2 = st.columns(2)
            mafft_op = c1.text_input("Op", value="1.53")
            mafft_ep = c2.text_input("Ep", value="0.0")
        with t2:
            trimal_binary = st.text_input("trimAl Path", value=trimal_def)
            use_trimal = st.checkbox("trimAlを使用", value=False)
            trimal_method = st.selectbox("Method", ["automated1", "gappyout"])
        with t3:
            iqtree_binary = st.text_input("IQ-TREE Path", value=iqtree_def)
            bootstrap = st.number_input("Bootstrap", 1000, step=100)
            model_sel = st.text_input("Model", placeholder="e.g. GTR+G")

    if 'phylo_step' not in st.session_state: st.session_state.phylo_step = 1
    if 'phylo_aligned_df' not in st.session_state: st.session_state.phylo_aligned_df = None

    uploaded_file = st.file_uploader("FASTAアップロード", type=["fasta", "fas", "fa"], key="phylo_up")

    if uploaded_file:
        if st.session_state.get('current_phylo_file') != uploaded_file.name:
            st.session_state.phylo_step = 1
            st.session_state.current_phylo_file = uploaded_file.name
            st.session_state.phylo_aligned_df = None
            stringio = uploaded_file.getvalue().decode("utf-8")
            raw_seqs = list(SeqIO.parse(StringIO(stringio), "fasta"))
            data = [{"Include": True, "ID": s.id, "Sequence": str(s.seq)} for s in raw_seqs]
            st.session_state.phylo_initial_df = pd.DataFrame(data)

        # Step 1: アラインメント
        if st.session_state.phylo_step == 1:
            st.subheader("1. アラインメント")
            input_df = st.data_editor(st.session_state.phylo_initial_df, key="phylo_ed1", hide_index=True)
            
            if st.button("🚀 解析を実行", key="phylo_run"):
                sel = input_df[input_df["Include"]==True]
                if len(sel) < 2: st.error("最低2配列必要です"); st.stop()
                
                with tempfile.TemporaryDirectory() as td:
                    inp = os.path.join(td, "in.fa")
                    out_aln = os.path.join(td, "out.aln")
                    
                    # 入力ファイル作成
                    SeqIO.write([SeqRecord(Seq(r["Sequence"]), id=r["ID"], description="") for i,r in sel.iterrows()], inp, "fasta")
                    
                    # MAFFT実行
                    cmd_mafft = [mafft_binary, mafft_algo, "--op", mafft_op, "--ep", mafft_ep, inp]
                    
                    # 【デバッグ用】実行するコマンドを表示
                    st.info(f"実行コマンド: `{' '.join(cmd_mafft)}`")
                    
                    with st.spinner("Running MAFFT..."):
                        # common.pyで修正したシンプルな実行関数
                        res = run_command(cmd_mafft)
                        
                        # エラーチェック
                        if res.returncode != 0:
                            st.error("❌ MAFFT Error: アラインメントに失敗しました。")
                            # エラー詳細を表示
                            st.text_area("Error Log", res.stderr, height=100)
                            st.stop()
                        
                        # 成功したらファイルに書き込む (ここで書き込む)
                        with open(out_aln, "w") as f:
                            f.write(res.stdout)

                    # trimAl
                    final_aln_file = out_aln
                    if use_trimal:
                        out_trim = os.path.join(td, "out_trim.aln")
                        cmd_trimal = [trimal_binary, "-in", out_aln, "-out", out_trim, "-" + trimal_method]
                        st.info(f"trimAl実行: `{' '.join(cmd_trimal)}`")
                        
                        res_t = run_command(cmd_trimal)
                        if res_t.returncode != 0:
                            st.error("trimAl Error")
                            st.code(res_t.stderr)
                            st.stop()
                        final_aln_file = out_trim

                    # 結果読み込み
                    recs = list(SeqIO.parse(final_aln_file, "fasta"))
                    if not recs:
                        st.error("結果ファイルが空です。")
                        st.stop()
                        
                    st.session_state.phylo_aligned_df = pd.DataFrame([{"Include":True, "ID":s.id, "Sequence":str(s.seq)} for s in recs])
                    st.session_state.phylo_step = 2
                    st.rerun()

        # Step 2: 確認・編集
        elif st.session_state.phylo_step == 2:
            st.subheader("2. 確認・編集")
            if st.session_state.phylo_aligned_df is None or st.session_state.phylo_aligned_df.empty:
                 st.warning("データがありません。最初からやり直してください。")
                 if st.button("Reset"): st.session_state.phylo_step=1; st.rerun()
            else:
                edf = st.data_editor(st.session_state.phylo_aligned_df, key="phylo_ed2", hide_index=True)
                st.markdown(generate_alignment_html_from_df(edf), unsafe_allow_html=True)
                
                c1, c2 = st.columns(2)
                with c1:
                    if st.button("🔄 再アラインメント", key="phylo_realign"):
                        # ガード: カラムチェック
                        if "Include" not in edf.columns:
                            st.error("データ破損: Include列がありません")
                            st.stop()
                        sel = edf[edf["Include"]==True]
                        new_data = [{"Include":True,"ID":r["ID"],"Sequence":r["Sequence"].replace("-","")} for i,r in sel.iterrows()]
                        st.session_state.phylo_initial_df = pd.DataFrame(new_data)
                        st.session_state.phylo_step = 1
                        st.rerun()
                with c2:
                    if st.button("🌳 IQ-TREEを実行", key="phylo_iqtree"):
                        sel = edf[edf["Include"]==True]
                        if len(sel) < 3: st.error("Min 3 seqs required"); st.stop()

                        with tempfile.TemporaryDirectory() as td:
                            aln = os.path.join(td, "aln.fa")
                            SeqIO.write([SeqRecord(Seq(r["Sequence"]), id=r["ID"], description="") for i,r in sel.iterrows()], aln, "fasta")
                            
                            cmd = [iqtree_binary, "-s", aln, "-bb", str(bootstrap), "-pre", os.path.join(td,"out"), "-nt", "AUTO"]
                            if model_sel: cmd.extend(["-m", model_sel])
                            
                            st.info(f"IQ-TREE実行: `{' '.join(cmd)}`")
                            
                            with st.spinner("Running IQ-TREE..."):
                                res = run_command(cmd)
                                if res.returncode != 0:
                                    st.error("IQ-TREE Error")
                                    st.code(res.stderr)
                                    # 続行可能な場合もあるのでstopしない
                                
                                if os.path.exists(os.path.join(td,"out.treefile")):
                                    with open(os.path.join(td,"out.treefile"), "r") as f: st.session_state['phylo_tree'] = f.read()
                                if os.path.exists(os.path.join(td,"out.iqtree")):
                                    with open(os.path.join(td,"out.iqtree"), "r") as f: st.session_state['phylo_report'] = f.read()
                                
                                st.session_state['phylo_log'] = res.stdout
                                st.session_state.phylo_step = 3
                                st.rerun()

        # Step 3: 完了
        elif st.session_state.phylo_step == 3:
            st.subheader("3. 解析完了")
            c1, c2 = st.columns(2)
            if 'phylo_tree' in st.session_state: c1.download_button("📥 Treefile", st.session_state['phylo_tree'], "phylo.treefile")
            if 'phylo_report' in st.session_state: c2.download_button("📄 Report", st.session_state['phylo_report'], "report.iqtree")
            with st.expander("Log"): st.code(st.session_state.get('phylo_log'))
            if st.button("最初から", key="phylo_reset"): st.session_state.phylo_step=1; st.rerun()
