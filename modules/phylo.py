import streamlit as st
import pandas as pd
import os
import tempfile
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from io import StringIO
from .common import find_tool_robust, run_command, generate_alignment_html_from_df

def app_phylo():
    st.header("🌳 PhyloPipeline Pro")
    st.info("MAFFT → (trimAl) → 編集 → IQ-TREE")

    mafft_def = find_tool_robust(["mafft"])
    trimal_def = find_tool_robust(["trimal"])
    iqtree_def = find_tool_robust(["iqtree2", "iqtree", "iqtree3"]) 
    
    with st.expander("🔧 ツール詳細設定", expanded=False):
        tab1, tab2, tab3 = st.tabs(["MAFFT", "trimAl", "IQ-TREE"])
        with tab1:
            c1, c2 = st.columns(2)
            mafft_binary = c1.text_input("MAFFT Path", value=mafft_def if mafft_def else "")
            mafft_algo = c2.selectbox("Algorithm", ["--auto", "--linsi", "--fftnsi"], index=0)
            c3, c4 = st.columns(2)
            mafft_op = c3.text_input("Op", value="1.53")
            mafft_ep = c4.text_input("Ep", value="0.0")
        with tab2:
            trimal_binary = st.text_input("trimAl Path", value=trimal_def if trimal_def else "")
            use_trimal = st.checkbox("trimAlを使用する", value=False)
            trimal_method = st.selectbox("Method", ["automated1", "gappyout", "strictplus"], index=0)
        with tab3:
            iqtree_binary = st.text_input("IQ-TREE Path", value=iqtree_def if iqtree_def else "")
            c5, c6 = st.columns(2)
            bootstrap = c5.number_input("Bootstrap", value=1000, step=100)
            model_sel = c6.text_input("Model", placeholder="e.g. GTR+G")

    if not mafft_binary or not iqtree_binary:
        st.error("⚠️ ツールが見つかりません。")
        return

    # ステート管理 (キー名を固有のものに変更)
    if 'phylo_step' not in st.session_state: st.session_state.phylo_step = 1
    if 'phylo_aligned_df' not in st.session_state: st.session_state.phylo_aligned_df = None

    uploaded_file = st.file_uploader("FASTAファイルをアップロード", type=["fasta", "fas", "fa"], key="phylo_uploader")

    if uploaded_file:
        if st.session_state.get('phylo_current_file') != uploaded_file.name:
            st.session_state.phylo_step = 1
            st.session_state.phylo_current_file = uploaded_file.name
            st.session_state.phylo_aligned_df = None
            stringio = uploaded_file.getvalue().decode("utf-8")
            raw_seqs = list(SeqIO.parse(StringIO(stringio), "fasta"))
            data = [{"Include": True, "ID": s.id, "Sequence": str(s.seq)} for s in raw_seqs]
            st.session_state.phylo_initial_df = pd.DataFrame(data)

        # Step 1
        if st.session_state.phylo_step == 1:
            st.subheader("1. アラインメント")
            input_df = st.data_editor(st.session_state.phylo_initial_df, key="phylo_ed1", hide_index=True)
            
            if st.button("🚀 解析を実行", key="phylo_btn_run"):
                sel = input_df[input_df["Include"]==True]
                if len(sel) < 2: st.error("最低2配列必要です"); st.stop()
                
                with tempfile.TemporaryDirectory() as td:
                    inp = os.path.join(td, "in.fa")
                    out_aln = os.path.join(td, "out.aln")
                    SeqIO.write([SeqRecord(Seq(r["Sequence"]), id=r["ID"], description="") for i,r in sel.iterrows()], inp, "fasta")
                    
                    # MAFFT
                    cmd_mafft = [mafft_binary, mafft_algo, "--op", mafft_op, "--ep", mafft_ep, inp]
                    with st.spinner("Running MAFFT..."):
                        with open(out_aln, "w") as f:
                            run_command(cmd_mafft, stdout=f)

                    # trimAl
                    final_aln = out_aln
                    if use_trimal and trimal_binary:
                        out_trim = os.path.join(td, "out_trim.aln")
                        cmd_trimal = [trimal_binary, "-in", out_aln, "-out", out_trim, "-" + trimal_method]
                        run_command(cmd_trimal)
                        final_aln = out_trim

                    recs = list(SeqIO.parse(final_aln, "fasta"))
                    st.session_state.phylo_aligned_df = pd.DataFrame([{"Include":True,"ID":s.id,"Sequence":str(s.seq)} for s in recs])
                    st.session_state.phylo_step = 2
                    st.rerun()

        # Step 2
        elif st.session_state.phylo_step == 2:
            st.subheader("2. 確認・編集")
            edf = st.data_editor(st.session_state.phylo_aligned_df, key="phylo_ed2", hide_index=True)
            st.markdown(generate_alignment_html_from_df(edf), unsafe_allow_html=True)
            
            c1, c2 = st.columns(2)
            with c1:
                if st.button("🔄 再アラインメント", key="phylo_btn_realign"):
                    sel = edf[edf["Include"]==True]
                    new_data = [{"Include":True,"ID":r["ID"],"Sequence":r["Sequence"].replace("-","")} for i,r in sel.iterrows()]
                    st.session_state.phylo_initial_df = pd.DataFrame(new_data)
                    st.session_state.phylo_step = 1
                    st.rerun()
            with c2:
                if st.button("🌳 IQ-TREEを実行", key="phylo_btn_iqtree"):
                    sel = edf[edf["Include"]==True]
                    with tempfile.TemporaryDirectory() as td:
                        aln = os.path.join(td, "aln.fa")
                        SeqIO.write([SeqRecord(Seq(r["Sequence"]), id=r["ID"], description="") for i,r in sel.iterrows()], aln, "fasta")
                        
                        cmd = [iqtree_binary, "-s", aln, "-bb", str(bootstrap), "-pre", os.path.join(td,"out"), "-nt", "AUTO"]
                        if model_sel: cmd.extend(["-m", model_sel])
                        
                        with st.spinner("Running IQ-TREE..."):
                            res = run_command(cmd)
                            if os.path.exists(os.path.join(td,"out.treefile")):
                                with open(os.path.join(td,"out.treefile"), "r") as f: st.session_state['phylo_tree'] = f.read()
                            if os.path.exists(os.path.join(td,"out.iqtree")):
                                with open(os.path.join(td,"out.iqtree"), "r") as f: st.session_state['phylo_report'] = f.read()
                            st.session_state['phylo_log'] = res.stdout
                            st.session_state.phylo_step = 3
                            st.rerun()

        # Step 3
        elif st.session_state.phylo_step == 3:
            st.subheader("3. 完了")
            if 'phylo_tree' in st.session_state: st.download_button("📥 系統樹 (.treefile)", st.session_state['phylo_tree'], "phylo.treefile")
            if 'phylo_report' in st.session_state: st.download_button("📄 レポート (.iqtree)", st.session_state['phylo_report'], "report.iqtree")
            with st.expander("Log"): st.code(st.session_state.get('phylo_log'))
            if st.button("最初から", key="phylo_btn_reset"): st.session_state.phylo_step=1; st.rerun()