import streamlit as st
import pandas as pd
import os
import tempfile
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import pdist, squareform
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from io import StringIO
from .common import find_tool_path, generate_alignment_html_from_df, run_command

# --- 種区分解析 (ASAP-like Distance Clustering) ---
def run_simple_asap(aligned_fasta_path, threshold=0.02):
    """
    簡易的な距離ベースの種区分解析 (Barcode Gap Analysis)
    """
    # 配列読み込み
    seqs = list(SeqIO.parse(aligned_fasta_path, "fasta"))
    if len(seqs) < 3: return None, "配列数が少なすぎます"
    
    ids = [s.id for s in seqs]
    
    # 配列をnumpy配列化
    max_len = max(len(s.seq) for s in seqs)
    matrix = []
    for s in seqs:
        # ギャップ含めて数値化
        seq_arr = np.array(list(str(s.seq).upper().ljust(max_len, '-')))
        matrix.append(seq_arr)
    matrix = np.array(matrix)
    
    # 距離行列計算 (p-distance)
    def p_dist(s1, s2):
        valid = (s1 != '-') & (s2 != '-') & (s1 != 'N') & (s2 != 'N')
        if np.sum(valid) == 0: return 0.0
        diff = (s1[valid] != s2[valid])
        return np.sum(diff) / np.sum(valid)

    n = len(seqs)
    dist_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            d = p_dist(matrix[i], matrix[j])
            dist_matrix[i, j] = dist_matrix[j, i] = d
    
    # クラスタリング (UPGMA)
    Z = linkage(squareform(dist_matrix), method='average')
    
    # 閾値でグループ分け
    clusters = fcluster(Z, t=threshold, criterion='distance')
    
    # 結果整形
    result_df = pd.DataFrame({"ID": ids, "Cluster": clusters})
    result_df = result_df.sort_values("Cluster")
    
    return result_df, dist_matrix

# --- メイン ---
def app_phylo():
    st.header("🌳 PhyloPipeline Pro")
    st.info("MAFFT → (trimAl) → 編集 → IQ-TREE + 種区分解析")

    # ツールパス探索
    mafft_def = find_tool_path("mafft") or "mafft"
    trimal_def = find_tool_path("trimal") or "trimal"
    iqtree_def = find_tool_path("iqtree") or find_tool_path("iqtree2") or "iqtree2"
    
    # --- ツール設定 ---
    with st.expander("🔧 ツール詳細設定 (Tool Settings)", expanded=False):
        c1, c2, c3 = st.columns(3)
        with c1:
            st.markdown("#### MAFFT")
            mafft_bin = st.text_input("Path", value=mafft_def, key="m_path")
            mafft_algo = st.selectbox("Algo", ["--auto", "--linsi", "--fftnsi"], key="m_algo")
            mafft_op = st.text_input("Op", value="1.53", key="m_op")
        with c2:
            st.markdown("#### trimAl")
            trimal_bin = st.text_input("Path", value=trimal_def, key="t_path")
            use_trimal = st.checkbox("Use trimAl", value=False, key="t_use")
            trimal_met = st.selectbox("Method", ["automated1", "gappyout"], key="t_met")
        with c3:
            st.markdown("#### IQ-TREE")
            iqtree_bin = st.text_input("Path", value=iqtree_def, key="i_path")
            boot = st.number_input("Bootstrap", 1000, step=100, key="i_boot")
            # モデル選択
            model_list = ["Auto (ModelFinder)", "GTR+G", "HKY+G", "TIM2+I+G", "GTR+I+G"]
            model_sel_ui = st.selectbox("Model", model_list, key="i_model_sel")
            model_str = "" if "Auto" in model_sel_ui else model_sel_ui

    # --- ステート管理 ---
    if 'phylo_step' not in st.session_state: st.session_state.phylo_step = 1
    if 'phylo_aligned_df' not in st.session_state: st.session_state.phylo_aligned_df = None
    if 'show_dots_mode' not in st.session_state: st.session_state.show_dots_mode = False

    uploaded_file = st.file_uploader("FASTAファイルをアップロード", type=["fasta", "fas", "fa"], key="phylo_up")

    # --- 編集用ダイアログ関数 ---
    @st.dialog("配列エディタ (Alignment Editor)", width="large")
    def open_editor():
        st.write("配列の選択・除外や内容の確認ができます。")
        show_dots = st.toggle("同一塩基をドット(.)で表示", value=st.session_state.show_dots_mode)
        st.session_state.show_dots_mode = show_dots
        
        edited = st.data_editor(
            st.session_state.phylo_aligned_df, 
            hide_index=True, 
            use_container_width=True,
            height=400
        )
        st.markdown("###### プレビュー")
        st.markdown(
            generate_alignment_html_from_df(edited, max_seqs=20, show_dots=show_dots), 
            unsafe_allow_html=True
        )
        if st.button("変更を保存して閉じる", type="primary"):
            st.session_state.phylo_aligned_df = edited
            st.rerun()

    # --- メインロジック ---
    if uploaded_file:
        if st.session_state.get('current_phylo_file') != uploaded_file.name:
            st.session_state.phylo_step = 1
            st.session_state.current_phylo_file = uploaded_file.name
            st.session_state.phylo_aligned_df = None
            
            # 【修正】文字コード自動判別
            file_bytes = uploaded_file.getvalue()
            decoded_string = None
            # UTF-8 -> Shift-JIS -> Latin-1 の順で試行
            for enc in ['utf-8', 'shift_jis', 'cp932', 'latin-1']:
                try:
                    decoded_string = file_bytes.decode(enc)
                    break
                except UnicodeDecodeError:
                    continue
            
            if decoded_string is None:
                st.error("ファイルの文字コードを判別できませんでした。UTF-8で保存し直してください。")
                st.stop()

            stringio = StringIO(decoded_string)
            try:
                raw_seqs = list(SeqIO.parse(stringio, "fasta"))
                if not raw_seqs:
                    st.warning("有効なシーケンスが見つかりませんでした。ファイル形式を確認してください。")
                    st.stop()
                data = [{"Include": True, "ID": s.id, "Sequence": str(s.seq)} for s in raw_seqs]
                st.session_state.phylo_initial_df = pd.DataFrame(data)
            except Exception as e:
                st.error(f"FASTA解析エラー: {e}")
                st.stop()

        # Step 1: アラインメント
        if st.session_state.phylo_step == 1:
            st.subheader("1. アラインメント実行")
            
            # データフレームが存在することを確認
            if 'phylo_initial_df' in st.session_state and not st.session_state.phylo_initial_df.empty:
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
                            with open(out_aln, "w") as f: run_command(cmd, stdout=f)
                        
                        final_aln = out_aln
                        if use_trimal:
                            trim = os.path.join(td, "trim.fa")
                            cmd_t = [trimal_bin, "-in", out_aln, "-out", trim, "-" + trimal_met]
                            run_command(cmd_t)
                            final_aln = trim

                        recs = list(SeqIO.parse(final_aln, "fasta"))
                        st.session_state.phylo_aligned_df = pd.DataFrame([{"Include":True, "ID":s.id, "Sequence":str(s.seq)} for s in recs])
                        st.session_state.phylo_step = 2
                        st.rerun()
            else:
                st.error("データが読み込まれていません。")

        # Step 2: 確認・編集・解析
        elif st.session_state.phylo_step == 2:
            st.subheader("2. アラインメント確認・解析")
            
            if st.session_state.phylo_aligned_df is None or st.session_state.phylo_aligned_df.empty:
                st.warning("データが空です。Step 1に戻ってください。")
                if st.button("戻る"): st.session_state.phylo_step = 1; st.rerun()
            else:
                st.markdown(
                    generate_alignment_html_from_df(
                        st.session_state.phylo_aligned_df, 
                        show_dots=st.session_state.show_dots_mode
                    ), 
                    unsafe_allow_html=True
                )
                
                col_tools = st.columns([1, 1, 1, 1])
                with col_tools[0]:
                    if st.button("🔍 エディタを開く (編集)", use_container_width=True):
                        open_editor()
                
                with col_tools[1]:
                    if st.button("🔄 再アラインメント", use_container_width=True):
                        sel = st.session_state.phylo_aligned_df[st.session_state.phylo_aligned_df["Include"]==True]
                        new_data = [{"Include":True, "ID":r["ID"], "Sequence":r["Sequence"].replace("-","")} for i,r in sel.iterrows()]
                        st.session_state.phylo_initial_df = pd.DataFrame(new_data)
                        st.session_state.phylo_step = 1
                        st.rerun()

                st.divider()
                c_iq, c_asap = st.columns(2)
                
                # --- IQ-TREE ---
                with c_iq:
                    st.markdown("### 🌳 系統樹構築 (IQ-TREE)")
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

                # --- ASAP ---
                with c_asap:
                    st.markdown("### 🧬 種区分解析 (ASAP-like)")
                    asap_thresh = st.slider("Distance Threshold (e.g. 0.02 = 2%)", 0.00, 0.10, 0.02, 0.005)
                    if st.button("Run Species Delimitation", use_container_width=True):
                        sel = st.session_state.phylo_aligned_df[st.session_state.phylo_aligned_df["Include"]==True]
                        with tempfile.TemporaryDirectory() as td:
                            aln = os.path.join(td, "aln_asap.fa")
                            SeqIO.write([SeqRecord(Seq(r["Sequence"]), id=r["ID"], description="") for i,r in sel.iterrows()], aln, "fasta")
                            
                            df_res, dist_mat = run_simple_asap(aln, asap_thresh)
                            
                            if df_res is not None:
                                st.session_state.asap_res = df_res
                                st.session_state.asap_dist = dist_mat
                                st.session_state.phylo_step = 3
                                st.rerun()
                            else:
                                st.error(dist_mat)

        # Step 3: 結果
        elif st.session_state.phylo_step == 3:
            st.subheader("3. 解析結果")
            t1, t2 = st.tabs(["IQ-TREE Results", "Species Delimitation"])
            
            with t1:
                if 'ptree' in st.session_state:
                    st.success("IQ-TREE Finished!")
                    c_d1, c_d2 = st.columns(2)
                    c_d1.download_button("📥 Treefile", st.session_state.ptree, "phylo.treefile")
                    c_d2.download_button("📄 Report", st.session_state.preport, "report.iqtree")
                    with st.expander("Show Log"): st.code(st.session_state.get('plog'))
                else:
                    st.info("No IQ-TREE results yet.")

            with t2:
                if 'asap_res' in st.session_state:
                    st.success("Delimitation Finished!")
                    st.dataframe(st.session_state.asap_res, use_container_width=True)
                    if 'asap_dist' in st.session_state:
                        st.write("Distance Matrix Heatmap:")
                        fig, ax = plt.subplots(figsize=(8, 6))
                        sns.heatmap(st.session_state.asap_dist, ax=ax, cmap="viridis")
                        st.pyplot(fig)
                    csv = st.session_state.asap_res.to_csv(index=False).encode('utf-8')
                    st.download_button("📥 Download Partition (CSV)", csv, "species_partition.csv")
                else:
                    st.info("No Delimitation results yet.")
            
            if st.button("最初からやり直す", key="p_rst"): st.session_state.phylo_step=1; st.rerun()
