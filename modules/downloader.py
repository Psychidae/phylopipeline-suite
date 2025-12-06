import streamlit as st
import pandas as pd
import time
from Bio import Entrez, SeqIO
from io import StringIO

def app_downloader():
    st.header("🧬 GenBank Sequence Downloader")
    st.info("CSVファイルから配列を一括取得")

    c1, c2 = st.columns([1, 2])
    with c1:
        st.subheader("設定")
        email = st.text_input("Email", placeholder="email@example.com")
        gene = st.text_input("Gene", value="Internal Transcribed Spacer")
        db = st.selectbox("DB", ["nucleotide", "protein"])
        ret = st.number_input("Max Count", 1, 100, 1)

    with c2:
        up = st.file_uploader("List (CSV/TXT)", type=["csv", "txt"])
        if up and email:
            if st.button("Download"):
                try:
                    df = pd.read_csv(up, header=None)
                    sps = df[0].tolist()
                    Entrez.email = email
                    recs = []
                    log = ""
                    bar = st.progress(0)
                    for i, sp in enumerate(sps):
                        bar.progress((i+1)/len(sps))
                        try:
                            h = Entrez.esearch(db=db, term=f'"{sp}"[Organism] AND {gene}[All Fields]', retmax=ret)
                            ids = Entrez.read(h)["IdList"]
                            if ids:
                                h2 = Entrez.efetch(db=db, id=ids, rettype="fasta", retmode="text")
                                for r in SeqIO.parse(StringIO(h2.read()), "fasta"):
                                    r.description = f"{sp} | {r.id}"
                                    r.id = f"{sp.replace(' ', '_')}_{r.id}"
                                    recs.append(r)
                                log += f"✅ {sp}: OK\n"
                            else: log += f"❌ {sp}: None\n"
                            time.sleep(0.5)
                        except Exception as e: log += f"⚠️ {sp}: {e}\n"
                    st.text_area("Log", log)
                    if recs:
                        out = StringIO()
                        SeqIO.write(recs, out, "fasta")
                        st.download_button("Download FASTA", out.getvalue(), "seqs.fasta")
                except Exception as e: st.error(e)
        elif up: st.warning("Enter Email")
