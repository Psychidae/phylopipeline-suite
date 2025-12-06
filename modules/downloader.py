import streamlit as st
import pandas as pd
import time
from Bio import Entrez, SeqIO
from io import StringIO

def app_downloader():
    st.header("🧬 GenBank Sequence Downloader")
    st.info("CSVファイル（種名リスト）から配列を一括取得します。")

    col1, col2 = st.columns([1, 2])
    with col1:
        st.subheader("設定")
        email = st.text_input("Email (必須)", placeholder="email@example.com")
        target_gene = st.text_input("ターゲット遺伝子", value="Internal Transcribed Spacer")
        db_select = st.selectbox("データベース", ["nucleotide", "protein"], index=0)
        max_ret = st.number_input("取得数", 1, 100, 1)

    with col2:
        uploaded_file = st.file_uploader("種名リスト (CSV/TXT)", type=["csv", "txt"])
        
        if uploaded_file and email:
            if st.button("🚀 ダウンロード開始"):
                try:
                    df = pd.read_csv(uploaded_file, header=None)
                    species_list = df[0].tolist()
                    Entrez.email = email
                    fasta_records = []
                    log_text = ""
                    prog = st.progress(0)
                    
                    for i, sp in enumerate(species_list):
                        prog.progress((i + 1) / len(species_list))
                        term = f'"{sp}"[Organism] AND {target_gene}[All Fields]' 
                        try:
                            handle = Entrez.esearch(db=db_select, term=term, retmax=max_ret)
                            record = Entrez.read(handle)
                            if record["IdList"]:
                                handle = Entrez.efetch(db=db_select, id=record["IdList"], rettype="fasta", retmode="text")
                                for r in SeqIO.parse(StringIO(handle.read()), "fasta"):
                                    r.description = f"{sp} | {r.id}"
                                    r.id = f"{sp.replace(' ', '_')}_{r.id}"
                                    fasta_records.append(r)
                                log_text += f"✅ {sp}: OK\n"
                            else:
                                log_text += f"❌ {sp}: None\n"
                            time.sleep(0.5)
                        except Exception as e:
                            log_text += f"⚠️ {sp}: Error {e}\n"
                    
                    st.text_area("Log", log_text)
                    if fasta_records:
                        out = StringIO()
                        SeqIO.write(fasta_records, out, "fasta")
                        st.download_button("📥 Download FASTA", out.getvalue(), "genbank_seqs.fasta")
                except Exception as e:
                    st.error(f"Error: {e}")
