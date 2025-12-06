import streamlit as st
import pandas as pd
import os
import tempfile
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from io import StringIO
# 共通関数をインポート
from modules.common import find_tool_path, generate_alignment_html_from_df, run_command

def app_phylo():
    st.header("🌳 PhyloPipeline Pro")
    st.info("MAFFT → (trimAl) → 編集 → IQ-TREE")

    # ツールパス設定
    mafft_def = find_tool_path("mafft") or "mafft"
    trimal_def = find_tool_path("trimal") or "trimal"
    iqtree_def = find_tool_path("iqtree") or find_tool_path("iqtree2") or "iqtree2"
    
    with st.expander("🔧 ツール詳細設定 (Parameters)", expanded=False):
        tab1, tab2, tab3 = st.tabs(["MAFFT", "trimAl", "IQ-TREE"])
        
        with tab1:
            st.caption("MAFFT (Alignment)")
            c1, c2 = st.columns(2)
            mafft_binary = c1.text_input("MAFFT Path", value=mafft_def)
            mafft_algo = c2.selectbox("Algorithm", ["--auto", "--linsi", "--ginsi", "--einsi", "--fftnsi"], index=0)
            c3, c4 = st.columns(2)
            mafft_op = c3.text_input("Op", value="1.53")
            mafft_ep = c4.text_input("Ep", value="0.0")
            mafft_extra = st.text_input("MAFFT Extra", value="")

        with tab2:
            st.caption("trimAl")
            trimal_binary = st.text_input("trimAl Path", value=trimal_def)
            use_trimal = st.checkbox("trimAlを使用する", value=False)
            trimal_method = st.selectbox("Method", ["automated1", "gappyout", "strictplus", "strict"], index=0)
            trimal_gt = st.slider("Gap Threshold", 0.0, 1.0, 0.9)
            trimal_extra = st.text_input("trimAl Extra", value="")

        with tab3:
            st.caption("IQ-TREE")
            iqtree_binary = st.text_input("IQ-TREE Path", value=iqtree_def)
            c5, c6 = st.columns(2)
            bootstrap = c5.number_input("Bootstrap", value=1000, step=100)
            model_sel = c6.text_input("Model", placeholder="e.g. GTR+G")
            c7, c8 = st.columns(2)
            iqtree_nt = c7.text_input("Threads", value="AUTO")
            iqtree_extra = st.text_input("IQ-TREE Extra", value="-bnni")

    # ステート管理
    if 'phylo_step' not in st.session_state: st.session_state.phylo_step = 1
    if 'phylo_aligned_df' not in st.session_state: st.session_state.phylo_aligned_df = None

    uploaded_file = st.file_uploader("FASTAファイルをアップロード", type=["fasta", "fas", "fa"], key="phylo_up")

    if uploaded_file:
        if st.session_state.get('current_phylo_file') != uploaded_file.name:
            st.session_state.phylo_step = 1
            st.session_state.current_phylo_file = uploaded_file.name
            st.session_state.phylo_aligned_df = None
            stringio = uploaded_file.getvalue().decode("utf-8")
            raw_seqs = list(SeqIO.parse(StringIO(stringio), "fasta"))
            data = [{"Include": True, "ID": s.id, "Sequence": str(s.seq)} for s in raw_seqs]
            st.session_state.phylo_initial_df = pd.DataFrame(data)

        # Step 1: アラインメント実行
        if st.session_state.phylo_step == 1:
            st.subheader("1. アラインメント")
            if 'phylo_initial_df' in st.session_state:
                input_df = st.data_editor(st.session_state.phylo_initial_df, key="phylo_ed1", hide_index=True)
            else:
                st.error("データの読み込みに失敗しました。")
                st.stop()
            
            if st.button("🚀 解析を実行", key="phylo_run"):
                sel = input_df[input_df["Include"]==True]
                if len(sel) < 2: st.error("最低2配列必要です"); st.stop()
                
                with tempfile.TemporaryDirectory() as td:
                    inp = os.path.join(td, "in.fa")
                    out_aln = os.path.join(td, "out.aln")
                    SeqIO.write([SeqRecord(Seq(r["Sequence"]), id=r["ID"], description="") for i,r in sel.iterrows()], inp, "fasta")
                    
                    # MAFFT
                    cmd_mafft = [mafft_binary, mafft_algo, "--op", mafft_op, "--ep", mafft_ep]
                    if mafft_extra: cmd_mafft.extend(mafft_extra.split())
                    cmd_mafft.extend([inp])
                    
                    with st.spinner("Running MAFFT..."):
                        # common.pyの修正版run_commandを使用
                        with open(out_aln, "w") as f:
                            res = run_command(cmd_mafft, stdout=f)
                            
                        # エラーチェック
                        if res.returncode != 0:
                            st.error("MAFFT Error: アラインメントに失敗しました。")
                            st.code(res.stderr)
                            st.stop()

                    # trimAl
                    final_aln_file = out_aln
                    if use_trimal:
                        out_trim = os.path.join(td, "out_trim.aln")
                        cmd_trimal = [trimal_binary, "-in", out_aln, "-out", out_trim]
                        if trimal_method == "automated1": cmd_trimal.append("-automated1")
                        elif trimal_method == "gappyout": cmd_trimal.append("-gappyout")
                        elif trimal_method == "strict": cmd_trimal.append("-strict")
                        elif trimal_method == "strictplus": cmd_trimal.append("-strictplus")
                        
                        if trimal_extra: cmd_trimal.extend(trimal_extra.split())
                        elif "-automated1" not in cmd_trimal and "-gappyout" not in cmd_trimal: 
                            cmd_trimal.extend(["-gt", str(trimal_gt)])
                        
                        run_command(cmd_trimal)
                        final_aln_file = out_trim

                    # 結果読み込みと空チェック
                    recs = list(SeqIO.parse(final_aln_file, "fasta"))
                    if not recs:
                        st.error("アラインメント結果が空です。MAFFTまたはtrimAlの設定を確認してください。")
                        # 空でもDataFrameの構造だけは維持する
                        st.session_state.phylo_aligned_df = pd.DataFrame(columns=["Include", "ID", "Sequence"])
                    else:
                        st.session_state.phylo_aligned_df = pd.DataFrame([{"Include":True, "ID":s.id, "Sequence":str(s.seq)} for s in recs])
                        st.session_state.phylo_step = 2
                        st.rerun()

        # Step 2: 確認・編集
        elif st.session_state.phylo_step == 2:
            st.subheader("2. 確認・編集")
            
            # データフレームが存在し、かつ空でないか確認
            if st.session_state.phylo_aligned_df is None or st.session_state.phylo_aligned_df.empty:
                st.warning("表示するアラインメントデータがありません。Step 1に戻ってください。")
                if st.button("戻る"):
                    st.session_state.phylo_step = 1
                    st.rerun()
            else:
                edf = st.data_editor(st.session_state.phylo_aligned_df, key="phylo_ed2", hide_index=True)
                st.markdown(generate_alignment_html_from_df(edf), unsafe_allow_html=True)
                
                c1, c2 = st.columns(2)
                with c1:
                    if st.button("🔄 再アラインメント", key="phylo_realign"):
                        sel = edf[edf["Include"]==True]
                        new_data = [{"Include":True,"ID":r["ID"],"Sequence":r["Sequence"].replace("-","")} for i,r in sel.iterrows()]
                        st.session_state.phylo_initial_df = pd.DataFrame(new_data)
                        st.session_state.phylo_step = 1
                        st.rerun()
                with c2:
                    if st.button("🌳 IQ-TREEを実行", key="phylo_iqtree"):
                        # ここでKeyErrorが起きていたのでガード
                        if "Include" not in edf.columns:
                            st.error("データエラー: 'Include'列が見つかりません。")
                            st.stop()
                            
                        sel = edf[edf["Include"]==True]
                        if len(sel) < 3:
                            st.error("系統樹構築には最低3配列必要です。")
                            st.stop()

                        with tempfile.TemporaryDirectory() as td:
                            aln = os.path.join(td, "aln.fa")
                            SeqIO.write([SeqRecord(Seq(r["Sequence"]), id=r["ID"], description="") for i,r in sel.iterrows()], aln, "fasta")
                            
                            cmd = [iqtree_binary, "-s", aln, "-bb", str(bootstrap), "-pre", os.path.join(td,"out"), "-nt", "AUTO"]
                            if model_sel: cmd.extend(["-m", model_sel])
                            if iqtree_extra: cmd.extend(iqtree_extra.split())
                            
                            with st.spinner("Running IQ-TREE..."):
                                res = run_command(cmd)
                                if res.returncode != 0:
                                    st.error("IQ-TREE Error")
                                    st.code(res.stderr)
                                    # ログは出すが停止はしない（解析できる場合もあるため）
                                
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
            if 'phylo_tree' in st.session_state: c1.download_button("📥 系統樹 (.treefile)", st.session_state['phylo_tree'], "phylo.treefile")
            if 'phylo_report' in st.session_state: c2.download_button("📄 レポート (.iqtree)", st.session_state['phylo_report'], "report.iqtree")
            with st.expander("Log"): st.code(st.session_state.get('phylo_log'))
            if st.button("最初からやり直す", key="phylo_reset"): st.session_state.phylo_step=1; st.rerun()
