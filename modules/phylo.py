import streamlit as st
import pandas as pd
import os
import tempfile
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from io import StringIO
from modules.common import find_tool_path, generate_alignment_html_from_df, run_command

def app_phylo():
    st.header("🌳 PhyloPipeline Pro")
    
    mafft = find_tool_path("mafft") or "mafft"
    trimal = find_tool_path("trimal") or "trimal"
    iqtree = find_tool_path("iqtree") or find_tool_path("iqtree2") or "iqtree2"
    
    with st.expander("🔧 Tools"):
        t1, t2, t3 = st.tabs(["MAFFT", "trimAl", "IQ-TREE"])
        with t1:
            m_bin = st.text_input("MAFFT Path", mafft)
            m_algo = st.selectbox("Algo", ["--auto", "--linsi"])
        with t2:
            t_bin = st.text_input("trimAl Path", trimal)
            use_t = st.checkbox("Use trimAl")
            t_met = st.selectbox("Method", ["automated1", "gappyout"])
        with t3:
            i_bin = st.text_input("IQ-TREE Path", iqtree)
            boot = st.number_input("Bootstrap", 1000, step=100)
            model = st.text_input("Model", placeholder="GTR+G")

    if 'p_step' not in st.session_state: st.session_state.p_step = 1
    up = st.file_uploader("FASTA", type=["fasta","fas","fa"], key="phylo_up")
    
    if up:
        if st.session_state.get('p_file') != up.name:
            st.session_state.p_step = 1
            st.session_state.p_file = up.name
            seqs = list(SeqIO.parse(StringIO(up.getvalue().decode("utf-8")), "fasta"))
            st.session_state.p_df = pd.DataFrame([{"Include":True, "ID":s.id, "Sequence":str(s.seq)} for s in seqs])

        if st.session_state.p_step == 1:
            st.subheader("1. Alignment")
            edf = st.data_editor(st.session_state.p_df, key="p_ed1")
            if st.button("Run"):
                sel = edf[edf["Include"]==True]
                with tempfile.TemporaryDirectory() as td:
                    inp = os.path.join(td, "in.fa")
                    out = os.path.join(td, "out.fa")
                    SeqIO.write([SeqRecord(Seq(r["Sequence"]), id=r["ID"], description="") for i,r in sel.iterrows()], inp, "fasta")
                    with open(out, "w") as f: run_command([m_bin, m_algo, inp], stdout=f)
                    final = out
                    if use_t:
                        trim = os.path.join(td, "trim.fa")
                        run_command([t_bin, "-in", out, "-out", trim, "-"+t_met])
                        final = trim
                    recs = list(SeqIO.parse(final, "fasta"))
                    st.session_state.p_adf = pd.DataFrame([{"Include":True, "ID":s.id, "Sequence":str(s.seq)} for s in recs])
                    st.session_state.p_step = 2
                    st.rerun()

        elif st.session_state.p_step == 2:
            st.subheader("2. Check")
            edf = st.data_editor(st.session_state.p_adf, key="p_ed2")
            st.markdown(generate_alignment_html_from_df(edf), unsafe_allow_html=True)
            if st.button("Run IQ-TREE"):
                sel = edf[edf["Include"]==True]
                with tempfile.TemporaryDirectory() as td:
                    aln = os.path.join(td, "aln.fa")
                    SeqIO.write([SeqRecord(Seq(r["Sequence"]), id=r["ID"], description="") for i,r in sel.iterrows()], aln, "fasta")
                    cmd = [i_bin, "-s", aln, "-bb", str(boot), "-pre", os.path.join(td,"out"), "-nt", "AUTO"]
                    if model: cmd.extend(["-m", model])
                    res = run_command(cmd)
                    if os.path.exists(os.path.join(td,"out.treefile")):
                        with open(os.path.join(td,"out.treefile")) as f: st.session_state.tree = f.read()
                    st.session_state.plog = res.stdout
                    st.session_state.p_step = 3
                    st.rerun()

        elif st.session_state.p_step == 3:
            st.subheader("3. Done")
            if 'tree' in st.session_state: st.download_button("Tree", st.session_state.tree, "tree.treefile")
            with st.expander("Log"): st.code(st.session_state.plog)
            if st.button("Reset"): st.session_state.p_step=1; st.rerun()
