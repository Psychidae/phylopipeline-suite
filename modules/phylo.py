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
from modules.phylo_logic import run_asap_scan, get_partition_by_threshold, generate_methods_log
from modules.phylo_editor import open_alignment_editor

def app_phylo():
    st.header("🌳 PhyloPipeline Pro")
    st.info("MAFFT → (trimAl) → 編集 → IQ-TREE + 種区分解析")

    # --- ステート初期化 ---
    if 'phylo_step' not in st.session_state: st.session_state.phylo_step = 1
    if 'phylo_aligned_df' not in st.session_state: st.session_state.phylo_aligned_df = None

    # --- ツールパス設定 ---
    mafft_def = find_tool_path("mafft") or "mafft"
    trimal_def = find_tool_path("trimal") or "trimal"
    iqtree_def = find_tool_path("iqtree") or find_tool_path("iqtree2") or "iqtree2"
    
    # --- ツール詳細設定 ---
    # --- ツール詳細設定 ---
    with st.expander("🔧 ツール詳細設定 (Tool Settings)", expanded=False):
        t_mafft, t_trimal, t_iqtree = st.tabs(["MAFFT", "trimAl", "IQ-TREE"])
        
        with t_mafft:
            st.markdown("#### MAFFT Settings")
            c1, c2 = st.columns(2)
            mafft_bin = c1.text_input("Path", value=mafft_def, key="m_path_v2")
            mafft_algo = c2.selectbox("Algo", ["--auto", "--linsi", "--fftnsi"], key="m_algo")
            c3, c4 = st.columns(2)
            mafft_op = c3.text_input("Op (Gap Open)", value="1.53", key="m_op")
            mafft_ep = c4.text_input("Ep (Offset)", value="0.0", key="m_ep")
            
        with t_trimal:
            st.markdown("#### trimAl Settings")
            use_trimal = st.checkbox("Use trimAl", value=False, key="t_use")
            c1, c2 = st.columns(2)
            trimal_bin = c1.text_input("Path", value=trimal_def, key="t_path_v2")
            trimal_met = c2.selectbox("Method", ["automated1", "gappyout"], key="t_met")
            
        with t_iqtree:
            st.markdown("#### IQ-TREE Settings")
            iqtree_bin = st.text_input("Path", value=iqtree_def, key="i_path_v2")
            boot = st.number_input("Bootstrap", 1000, step=100, key="i_boot")
            
            # Phase 2: Outgroup Selection
            st.markdown("##### Outgroup")
            outgroup_list = ["(None)"]
            if st.session_state.phylo_aligned_df is not None:
                # Include=Trueのもののみ
                sel_og = st.session_state.phylo_aligned_df[st.session_state.phylo_aligned_df["Include"]==True]
                outgroup_list += sel_og["ID"].tolist()
            
            outgroup_sel = st.selectbox("Select Outgroup", outgroup_list, key="i_outgroup")
            
            st.markdown("##### Model Selection")
            model_list = ["Auto (ModelFinder)", "GTR+G", "HKY+G", "TIM2+I+G", "GTR+I+G", "Custom"]
            model_sel_ui = st.radio("Model", model_list, key="i_model_sel", horizontal=True)
            
            model_str = ""
            if "Auto" in model_sel_ui:
                model_str = ""
            elif model_sel_ui == "Custom":
                model_str = st.text_input("Enter Model String (e.g. TVM+F+R3)", value="", key="i_model_custom")
            else:
                model_str = model_sel_ui

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
                            cmd = [mafft_bin, "--quiet", mafft_algo, "--op", mafft_op, "--ep", mafft_ep, inp]
                            # stderrをPIPEにリダイレクトして /dev/stderr への書き込みエラーを回避
                            with open(out_aln, "w") as f: 
                                res = run_command(cmd, stdout=f, stderr=subprocess.PIPE)
                                if res.returncode != 0:
                                     st.error(f"MAFFT Error: {res.stderr}")
                        
                        final_aln = out_aln
                        if use_trimal:
                            trim = os.path.join(td, "trim.fa")
                            cmd_t = [trimal_bin, "-in", out_aln, "-out", trim, "-" + trimal_met]
                            run_command(cmd_t)
                            final_aln = trim

                        # check file size
                        if os.path.getsize(final_aln) == 0:
                             st.error("アラインメント結果が空です。MAFFTが正しく実行されなかった可能性があります。")
                             if 'res' in locals() and res.stderr:
                                 with st.expander("MAFFT Error Log"):
                                     st.code(res.stderr)
                        else:
                            recs = list(SeqIO.parse(final_aln, "fasta"))
                            if not recs:
                                st.error("アラインメント結果のパースに失敗しました。")
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
                            if outgroup_sel and outgroup_sel != "(None)":
                                cmd.extend(["-o", outgroup_sel])
                            
                            with st.spinner("Running IQ-TREE..."):
                                res = run_command(cmd)
                                if os.path.exists(os.path.join(td,"out.treefile")):
                                    with open(os.path.join(td,"out.treefile")) as f: st.session_state.ptree = f.read()
                                if os.path.exists(os.path.join(td,"out.iqtree")):
                                    with open(os.path.join(td,"out.iqtree")) as f: st.session_state.preport = f.read()
                                st.session_state.plog = res.stdout
                                
                                # Phase 2: Generate Log
                                params = {
                                    "Tool": "IQ-TREE",
                                    "Bootstrap": boot,
                                    "Model": model_str if model_str else "Auto",
                                    "Outgroup": outgroup_sel,
                                    "MAFFT Algo": mafft_algo,
                                    "trimAl": "Yes" if use_trimal else "No"
                                }
                                tools = {"IQ-TREE": iqtree_bin, "MAFFT": mafft_bin} 
                                st.session_state.method_log = generate_methods_log(tools, params)
                                
                                st.session_state.phylo_step = 3
                                st.rerun()

                # ASAP
                with c_asap:
                    st.markdown("### 🧬 種区分解析 (ASAP-like)")
                    st.info("閾値をスキャンして最適な種区分を探します。")
                    
                    # Phase 2: Distance Model Selection
                    dist_model = st.selectbox("Distance Model", ["p-dist", "k2p"], index=1, help="K2P is recommended for barcoding")

                    if st.button("Run ASAP Scan", use_container_width=True):
                        sel = st.session_state.phylo_aligned_df[st.session_state.phylo_aligned_df["Include"]==True]
                        with tempfile.TemporaryDirectory() as td:
                            aln = os.path.join(td, "aln_asap.fa")
                            SeqIO.write([SeqRecord(Seq(r["Sequence"]), id=r["ID"], description="") for i,r in sel.iterrows()], aln, "fasta")
                            
                            # モジュール呼び出し
                            df_scan, dist_mat, Z, ids = run_asap_scan(aln, model=dist_model)
                            
                            if df_scan is not None:
                                st.session_state.asap_scan = df_scan
                                st.session_state.asap_dist = dist_mat
                                st.session_state.asap_Z = Z
                                st.session_state.asap_ids = ids
                                st.session_state.phylo_step = 3
                                st.rerun()
                            else:
                                st.error(ids) # idsがエラーメッセージになる仕様


        # === Step 3: 結果 ===
        elif st.session_state.phylo_step == 3:
            st.subheader("3. 解析結果")
            t1, t2 = st.tabs(["IQ-TREE", "ASAP"])
            
            with t1:
                if 'ptree' in st.session_state:
                    st.success("Finished!")
                    c1, c2, c3 = st.columns(3)
                    c1.download_button("📥 Treefile", st.session_state.ptree, "phylo.treefile")
                    c2.download_button("📄 Report", st.session_state.preport, "report.iqtree")
                    
                    if 'method_log' in st.session_state:
                         c3.download_button("📝 Methods Log", st.session_state.method_log, "methods_log.txt")
                    
                    with st.expander("Log"): st.code(st.session_state.get('plog'))
                else:
                    st.info("IQ-TREE results not available.")

            with t2:
                if 'asap_scan' in st.session_state:
                    st.success("ASAP Scan Finished!")
                    
                    c_scan, c_det = st.columns([1, 1.5])
                    with c_scan:
                        st.subheader("1. Select Threshold")
                        st.caption("スキャン結果から閾値を選択してください。")
                        
                        # スキャン結果テーブル
                        st.dataframe(st.session_state.asap_scan, use_container_width=True, hide_index=True)
                        
                        # 閾値選択スライダー（スキャン範囲に基づき設定）
                        min_t = st.session_state.asap_scan["Threshold (Distance)"].min()
                        max_t = st.session_state.asap_scan["Threshold (Distance)"].max()
                        sel_t = st.slider("Select Threshold", float(min_t), float(max_t), 0.02, 0.005)
                    
                    with c_det:
                        st.subheader("2. Partition Result")
                        # 選択された閾値でパーティション取得
                        Z = st.session_state.asap_Z
                        ids = st.session_state.asap_ids
                        res_df = get_partition_by_threshold(Z, ids, sel_t)
                        
                        st.write(f"**Threshold: {sel_t}** -> **{len(res_df['Cluster'].unique())} Species**")
                        st.dataframe(res_df, use_container_width=True, hide_index=True)
                        
                        csv = res_df.to_csv(index=False).encode('utf-8')
                        st.download_button("📥 Download CSV", csv, f"species_t{sel_t}.csv")

                    st.divider()
                    if 'asap_dist' in st.session_state:
                        st.subheader("Distance Matrix")
                        fig, ax = plt.subplots(figsize=(8, 6))
                        sns.heatmap(st.session_state.asap_dist, ax=ax, cmap="viridis")
                        st.pyplot(fig)
                else:
                    st.info("ASAP results not available.")
            
            if st.button("最初からやり直す", key="p_rst"): st.session_state.phylo_step=1; st.rerun()
