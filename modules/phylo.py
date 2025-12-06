import streamlit as st
import pandas as pd
import os
import tempfile
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from io import StringIO
# 共通関数をインポート
from .common import find_tool_path, generate_alignment_html_from_df, run_command

def app_phylo():
    st.header("🌳 PhyloPipeline Pro")
    st.info("MAFFT → (trimAl) → 編集 → IQ-TREE")

    # クラウド環境ではパスが見つからなくてもコマンド名そのままで動くことが多い
    mafft_def = find_tool_path("mafft") or "mafft"
    trimal_def = find_tool_path("trimal") or "trimal"
    iqtree_def = find_tool_path("iqtree") or find_tool_path("iqtree2") or "iqtree2"
    
    with st.expander("🔧 ツール詳細設定", expanded=False):
        t1, t2, t3 = st.tabs(["MAFFT", "trimAl", "IQ-TREE"])
        with t1:
            mafft_bin = st.text_input("MAFFT Path", value=mafft_def)
            mafft_algo = st.selectbox("Algorithm", ["--auto", "--linsi", "--fftnsi"])
        with t2:
            trimal_bin = st.text_input("trimAl Path", value=trimal_def)
            use_trimal = st.checkbox("trimAlを使用", value=False)
            trimal_method = st.selectbox("Method", ["automated1", "gappyout", "strictplus"])
        with t3:
            iqtree_bin = st.text_input("IQ-TREE Path", value=iqtree_def)
            bootstrap = st.number_input("Bootstrap", 1000, step=100)
            model_sel = st.text_input("Model", placeholder="e.g. GTR+G")

    if 'phylo_step' not in st.session_state: st.session_state.phylo_step = 1
    
    uploaded = st.file_uploader("FASTAアップロード", type=["fasta", "fas", "fa"], key="phylo_up")
    
    if uploaded:
        if st.session_state.get('last_phylo_file') != uploaded.name:
            st.session_state.phylo_step = 1
            st.session_state.last_phylo_file = uploaded.name
            seqs = list(SeqIO.parse(StringIO(uploaded.getvalue().decode("utf-8")), "fasta"))
            st.session_state.phylo_init_df = pd.DataFrame([{"Include":True, "ID":s.id, "Sequence":str(s.seq)} for s in seqs])

        # Step 1
        if st.session_state.phylo_step == 1:
            st.subheader("1. アラインメント")
            edf = st.data_editor(st.session_state.phylo_init_df, hide_index=True)
            if st.button("🚀 解析実行", key="p_run"):
                sel = edf[edf["Include"]==True]
                with tempfile.TemporaryDirectory() as td:
                    inp = os.path.join(td, "in.fa")
                    out = os.path.join(td, "out.fa")
                    SeqIO.write([SeqRecord(Seq(r["Sequence"]), id=r["ID"], description="") for i,r in sel.iterrows()], inp, "fasta")
                    
                    with st.spinner("Running MAFFT..."):
                        with open(out, "w") as f:
                            run_command([mafft_bin, mafft_algo, inp], stdout=f)
                    
                    final = out
                    if use_trimal:
                        trim = os.path.join(td, "trim.fa")
                        run_command([trimal_bin, "-in", out, "-out", trim, "-" + trimal_method])
                        final = trim
                    
                    recs = list(SeqIO.parse(final, "fasta"))
                    st.session_state.phylo_aln_df = pd.DataFrame([{"Include":True, "ID":s.id, "Sequence":str(s.seq)} for s in recs])
                    st.session_state.phylo_step = 2
                    st.rerun()

        # Step 2
        elif st.session_state.phylo_step == 2:
            st.subheader("2. 確認・編集")
            edf = st.data_editor(st.session_state.phylo_aln_df, hide_index=True)
            st.markdown(generate_alignment_html_from_df(edf), unsafe_allow_html=True)
            c1, c2 = st.columns(2)
            if c1.button("🔄 再アラインメント", key="p_re"):
                sel = edf[edf["Include"]==True]
                st.session_state.phylo_init_df = pd.DataFrame([{"Include":True, "ID":r["ID"], "Sequence":r["Sequence"].replace("-","")} for i,r in sel.iterrows()])
                st.session_state.phylo_step = 1
                st.rerun()
            if c2.button("🌳 IQ-TREE実行", key="p_iq"):
                sel = edf[edf["Include"]==True]
                with tempfile.TemporaryDirectory() as td:
                    aln = os.path.join(td, "aln.fa")
                    SeqIO.write([SeqRecord(Seq(r["Sequence"]), id=r["ID"], description="") for i,r in sel.iterrows()], aln, "fasta")
                    cmd = [iqtree_bin, "-s", aln, "-bb", str(bootstrap), "-pre", os.path.join(td,"out"), "-nt", "AUTO"]
                    if model_sel: cmd.extend(["-m", model_sel])
                    with st.spinner("Running IQ-TREE..."):
                        res = run_command(cmd)
                        if os.path.exists(os.path.join(td,"out.treefile")):
                            with open(os.path.join(td,"out.treefile")) as f: st.session_state.ptree = f.read()
                        st.session_state.plog = res.stdout
                        st.session_state.phylo_step = 3
                        st.rerun()

        # Step 3
        elif st.session_state.phylo_step == 3:
            st.subheader("3. 完了")
            if 'ptree' in st.session_state: st.download_button("📥 Download Tree", st.session_state.ptree, "tree.treefile")
            with st.expander("Log"): st.code(st.session_state.get('plog'))
            if st.button("最初から", key="p_rst"): st.session_state.phylo_step=1; st.rerun()
