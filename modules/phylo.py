import streamlit as st
import pandas as pd
import os
import tempfile
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from io import StringIO

# モジュール分割した機能をインポート
from modules.common import find_tool_path, generate_alignment_html_from_df, run_command
from modules.phylo_logic import run_simple_asap_logic
from modules.phylo_editor import open_alignment_editor

def app_phylo():
    st.header("🌳 PhyloPipeline Pro")
    st.info("MAFFT → (trimAl) → 編集 → IQ-TREE + 種区分解析")

    # --- ツールパス設定 ---
    mafft_def = find_tool_path("mafft") or "mafft"
    trimal_def = find_tool_path("trimal") or "trimal"
    iqtree_def = find_tool_path("iqtree") or find_tool_path("iqtree2") or "iqtree2"
    
    # --- ツール詳細設定 ---
    with st.expander("🔧 ツール詳細設定 (Tool Settings)", expanded=False):
        c1, c2, c3 = st.columns(3)
        with c1:
            st.markdown("#### MAFFT")
            mafft_bin = st.text_input("Path", value=mafft_def, key="m_path")
            mafft_algo = st.selectbox("Algo", ["--auto", "--linsi", "--fftnsi"], key="m_algo")
            mafft_op = st.text_input("Op", value="1.53", key="m_op")
            mafft_ep = st.text_input("Ep", value="0.0", key="m_ep")
        with c2:
            st.markdown("#### trimAl")
            trimal_bin = st.text_input("Path", value=trimal_def, key="t_path")
            use_trimal = st.checkbox("Use trimAl", value=False, key="t_use")
            trimal_met = st.selectbox("Method", ["automated1", "gappyout"], key="t_met")
        with c3:
            st.markdown("#### IQ-TREE")
            iqtree_bin = st.text_input("Path", value=iqtree_def, key="i_path")
            boot = st.number_input("Bootstrap", 1000, step=100, key="i_boot")
            model_list = ["Auto (ModelFinder)", "GTR+G", "HKY+G", "TIM2+I+G", "GTR+I+G"]
            model_sel_ui = st.selectbox("Model", model_list, key="i_model_sel")
            model_str = "" if "Auto" in model_sel_ui else model_sel_ui

    # --- ステート初期化 ---
    if 'phylo_step' not in st.session_state: st.session_state.phylo_step = 1
    if 'phylo_aligned_df' not in st.session_state: st.session_state.phylo_aligned_df = None
    
    # --- ファイルアップロード ---
    uploaded_file = st.file_uploader("FASTAファイルをアップロード", type=["fasta", "fas", "fa"], key="phylo_up")

    if uploaded_file:
        # 新規ファイル読み込み処理
        if st.session_state.get('current_phylo_file') != uploaded_file.name:
            st.session_state.phylo_step = 1
            st.session_state.current_phylo_file = uploaded_file.name
            st.session_state.phylo_aligned_df = None
            
            # 文字コード対応
            file_bytes = uploaded_file.getvalue()
            decoded = None
            for enc in ['utf-8', 'shift_jis', 'cp932', 'latin-1']:
                try: decoded = file_bytes.decode(enc); break
                except: continue
            if decoded is None: st.error("Encode Error"); st.stop()

            try:
                raw_seqs = list(SeqIO.parse(StringIO(decoded), "fasta"))
                data = [{"Include": True, "ID": s.id, "Sequence": str(s.seq)} for s in raw_seqs]
                st.session_state.phylo_initial_df = pd.DataFrame(data)
            except Exception as e:
                st.error(f"Error parsing FASTA: {e}")
                st.stop()

        # === Step 1: アラインメント実行 ===
        if st.session_state.phylo_step == 1:
            st.subheader("1. アラインメント実行")
            if 'phylo_initial_df' in st.session_state:
                with st.expander("入力データ確認", expanded=True):
                    input_df = st.data_editor(st.session_state.phylo_initial_df, key="p_ed1", hide_index=True)
                
                if st.button("🚀 アラインメント開始 (MAFFT)", key="p_run", type="primary"):
                    sel = input_df[input_df["Include"]==True]
                    if len(sel) < 2: st.error("最低2配列必要です"); st.stop()
                    
                    with tempfile.TemporaryDirectory() as td:
                        inp = os.path.join(td, "in.fa")
                        out_aln = os.path.join(td, "out.aln")
                        SeqIO.write([SeqRecord(Seq(r["Sequence"]), id=r["ID"], description="") for i,r in sel.iterrows()], inp, "fasta")
                        
                        with st.spinner("Running MAFFT..."):
                            cmd = [mafft_bin, mafft_algo, "--op", mafft_op, "--ep", mafft_ep, inp]
                            # common.pyのrun_commandを使用（エラーハンドリング付き）
                            with open(out_aln, "w") as f: run_command(cmd, stdout=f)
                        
                        final_aln = out_aln
                        if use_trimal:
                            trim = os.path.join(td, "trim.fa")
                            cmd_t = [trimal_bin, "-in", out_aln, "-out", trim, "-" + trimal_met]
                            run_command(cmd_t)
                            final_aln = trim

                        recs = list(SeqIO.parse(final_aln, "fasta"))
                        if not recs:
                            st.error("アラインメント結果が空です。")
                        else:
                            st.session_state.phylo_aligned_df = pd.DataFrame([{"Include":True, "ID":s.id, "Sequence":str(s.seq)} for s in recs])
                            st.session_state.phylo_step = 2
                            st.rerun()
            else:
                st.error("データなし")

        # === Step 2: 確認・編集・解析 ===
        elif st.session_state.phylo_step == 2:
            st.subheader("2. アラインメント確認・解析")
            
            if st.session_state.phylo_aligned_df is None or st.session_state.phylo_aligned_df.empty:
                st.warning("データが空です。Step 1に戻ってください。")
                if st.button("戻る"): st.session_state.phylo_step = 1; st.rerun()
            else:
                # 簡易プレビュー
                st.markdown(generate_alignment_html_from_df(st.session_state.phylo_aligned_df), unsafe_allow_html=True)
                
                # ツールバー
                c_tools = st.columns([1, 1, 2])
                with c_tools[0]:
                    if st.button("🔍 エディタを開く", use_container_width=True):
                        open_alignment_editor(st.session_state.phylo_aligned_df) # モジュール呼び出し
                with c_tools[1]:
                    if st.button("🔄 再整列", use_container_width=True):
                        # ギャップ除去して再実行へ
                        sel = st.session_state.phylo_aligned_df[st.session_state.phylo_aligned_df["Include"]==True]
                        new_data = [{"Include":True, "ID":r["ID"], "Sequence":r["Sequence"].replace("-","")} for i,r in sel.iterrows()]
                        st.session_state.phylo_initial_df = pd.DataFrame(new_data)
                        st.session_state.phylo_step = 1
                        st.rerun()

                st.divider()
                c_iq, c_asap = st.columns(2)
                
                # IQ-TREE
                with c_iq:
                    st.markdown("### 🌳 系統樹構築")
                    if st.button("Run IQ-TREE", type="primary", use_container_width=True):
                        sel = st.session_state.phylo_aligned_df[st.session_state.phylo_aligned_df["Include"]==True]
                        with tempfile.TemporaryDirectory() as td:
                            aln = os.path.join(td, "aln.fa")
                            SeqIO.write([SeqRecord(Seq(r["Sequence"]), id=r["ID"], description="") for i,r in sel.iterrows()], aln, "fasta")
                            
                            cmd = [iqtree_bin, "-s", aln, "-bb", str(boot), "-pre", os.path.join(td,"out"), "-nt", "AUTO"]
                            if model_str: cmd.extend(["-m", model_str])
                            
                            with st.spinner("Running IQ-TREE..."):
                                res = run_command(cmd)
                                if os.path.exists(os.path.join(td,"out.treefile")):
                                    with open(os.path.join(td,"out.treefile")) as f: st.session_state.ptree = f.read()
                                if os.path.exists(os.path.join(td,"out.iqtree")):
                                    with open(os.path.join(td,"out.iqtree")) as f: st.session_state.preport = f.read()
                                st.session_state.plog = res.stdout
                                st.session_state.phylo_step = 3
                                st.rerun()

                # ASAP
                with c_asap:
                    st.markdown("### 🧬 種区分解析 (ASAP-like)")
                    asap_thresh = st.slider("Distance Threshold", 0.00, 0.10, 0.02, 0.005)
                    if st.button("Run Analysis", use_container_width=True):
                        sel = st.session_state.phylo_aligned_df[st.session_state.phylo_aligned_df["Include"]==True]
                        with tempfile.TemporaryDirectory() as td:
                            aln = os.path.join(td, "aln_asap.fa")
                            SeqIO.write([SeqRecord(Seq(r["Sequence"]), id=r["ID"], description="") for i,r in sel.iterrows()], aln, "fasta")
                            
                            # モジュール呼び出し
                            df_res, dist_mat = run_simple_asap_logic(aln, asap_thresh)
                            
                            if df_res is not None:
                                st.session_state.asap_res = df_res
                                st.session_state.asap_dist = dist_mat
                                st.session_state.phylo_step = 3
                                st.rerun()
                            else:
                                st.error(dist_mat)

        # === Step 3: 結果 ===
        elif st.session_state.phylo_step == 3:
            st.subheader("3. 解析結果")
            t1, t2 = st.tabs(["IQ-TREE", "ASAP"])
            
            with t1:
                if 'ptree' in st.session_state:
                    st.success("Finished!")
                    c1, c2 = st.columns(2)
                    c1.download_button("📥 Treefile", st.session_state.ptree, "phylo.treefile")
                    c2.download_button("📄 Report", st.session_state.preport, "report.iqtree")
                    with st.expander("Log"): st.code(st.session_state.get('plog'))
                else:
                    st.info("IQ-TREE results not available.")

            with t2:
                if 'asap_res' in st.session_state:
                    st.success("Finished!")
                    st.dataframe(st.session_state.asap_res, use_container_width=True)
                    if 'asap_dist' in st.session_state:
                        fig, ax = plt.subplots(figsize=(8, 6))
                        sns.heatmap(st.session_state.asap_dist, ax=ax, cmap="viridis")
                        st.pyplot(fig)
                    csv = st.session_state.asap_res.to_csv(index=False).encode('utf-8')
                    st.download_button("📥 Download CSV", csv, "species.csv")
                else:
                    st.info("ASAP results not available.")
            
            if st.button("最初からやり直す", key="p_rst"): st.session_state.phylo_step=1; st.rerun()
