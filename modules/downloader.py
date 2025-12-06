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
        email = st.text_input("Email (必須)", placeholder="your_email@example.com", help="NCBI利用規約により必須")
        target_gene = st.text_input("ターゲット遺伝子", value="Internal Transcribed Spacer", help="例: COI, 16S, rbcL")
        db_select = st.selectbox("データベース", ["nucleotide", "protein"], index=0)
        max_ret = st.number_input("1種あたりの最大取得数", min_value=1, value=1)

    with col2:
        st.subheader("実行")
        uploaded_file = st.file_uploader("種名リスト (CSV/TXT, ヘッダーなし)", type=["csv", "txt"])
        
        if uploaded_file and email:
            try:
                df = pd.read_csv(uploaded_file, header=None)
                species_list = df[0].tolist()
                st.write(f"✅ {len(species_list)} 種のリストを読み込みました。")
                
                if st.button("🚀 ダウンロード開始", key="btn_download"):
                    Entrez.email = email
                    fasta_records = []
                    log_text = ""
                    prog_bar = st.progress(0)
                    
                    for i, sp in enumerate(species_list):
                        prog_bar.progress((i + 1) / len(species_list))
                        term = f'"{sp}"[Organism] AND {target_gene}[All Fields]' 
                        
                        try:
                            handle = Entrez.esearch(db=db_select, term=term, retmax=max_ret)
                            record = Entrez.read(handle)
                            
                            if record["IdList"]:
                                handle = Entrez.efetch(db=db_select, id=record["IdList"], rettype="fasta", retmode="text")
                                seq_data = handle.read()
                                count = 0
                                for r in SeqIO.parse(StringIO(seq_data), "fasta"):
                                    clean_sp = sp.replace(" ", "_")
                                    r.description = f"{sp} | {r.id}"
                                    r.id = f"{clean_sp}_{r.id}"
                                    fasta_records.append(r)
                                    count += 1
                                log_text += f"✅ {sp}: {count}件取得\n"
                            else:
                                log_text += f"❌ {sp}: なし\n"
                            time.sleep(0.5)
                        except Exception as e:
                            log_text += f"⚠️ {sp}: エラー {e}\n"
                    
                    st.text_area("処理ログ", log_text, height=150)
                    
                    if fasta_records:
                        out = StringIO()
                        SeqIO.write(fasta_records, out, "fasta")
                        st.download_button("📥 FASTAをダウンロード", out.getvalue(), f"downloaded_seqs.fasta")
                    else:
                        st.warning("配列が取得できませんでした。")
            except Exception as e:
                st.error(f"ファイル読み込みエラー: {e}")
        elif uploaded_file:
            st.warning("メールアドレスを入力してください。")