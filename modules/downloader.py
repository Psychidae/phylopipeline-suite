import streamlit as st
import pandas as pd
import time
from Bio import Entrez, SeqIO
from io import StringIO, BytesIO
import xml.etree.ElementTree as ET

def app_downloader():
    st.header("🧬 GenBank Sequence Downloader")
    st.info("配列データ(FASTA)と、詳細なメタデータ(Excel)を一括取得します。")

    col1, col2 = st.columns([1, 1.5])
    with col1:
        st.subheader("設定")
        email = st.text_input("Email (必須)", placeholder="your_email@example.com", help="NCBIの利用規約により必須です")
        target_gene = st.text_input("ターゲット遺伝子", value="COI", help="例: COI, 16S, NADH dehydrogenase subunit 1")
        max_ret = st.number_input("1種あたりの最大取得数", 1, 100, 1)
        
        st.markdown("---")
        st.caption("絞り込みオプション (任意)")
        filter_author = st.text_input("登録者/著者名", placeholder="Smith J", help="この著者が関わった配列のみ取得")
        filter_journal = st.text_input("論文/雑誌名", placeholder="Nature", help="この雑誌/論文に含まれる配列のみ取得")

    with col2:
        st.subheader("リストアップロード")
        
        # テンプレートダウンロード
        template_csv = "Homo sapiens\nMus musculus\nDrosophila melanogaster\n"
        st.download_button(
            "📥 テンプレートCSVをダウンロード",
            template_csv,
            "species_list_template.csv",
            "text/csv",
            help="種名を1行に1つ記載したCSV/TXTファイルを作成してください"
        )
        
        uploaded_file = st.file_uploader("種名リスト (CSV/TXT, ヘッダーなし)", type=["csv", "txt"])
        
        if uploaded_file and email:
            if st.button("🚀 ダウンロード開始", type="primary"):
                try:
                    df_input = pd.read_csv(uploaded_file, header=None)
                    species_list = df_input[0].tolist()
                    Entrez.email = email
                    
                    fasta_records = []
                    metadata_list = []
                    
                    log_text = ""
                    prog_bar = st.progress(0)
                    status_text = st.empty()
                    
                    for i, sp in enumerate(species_list):
                        prog_bar.progress((i + 1) / len(species_list))
                        status_text.text(f"Searching: {sp}...")
                        
                        # 検索クエリ構築
                        term = f'"{sp}"[Organism] AND {target_gene}[All Fields]'
                        if filter_author:
                            term += f' AND {filter_author}[Author]'
                        if filter_journal:
                            term += f' AND {filter_journal}[Journal]'
                            
                        try:
                            # 1. ID検索
                            handle = Entrez.esearch(db="nucleotide", term=term, retmax=max_ret)
                            record = Entrez.read(handle)
                            id_list = record["IdList"]
                            
                            if not id_list:
                                log_text += f"❌ {sp}: なし (条件に一致するデータがありません)\n"
                                continue
                            
                            # 2. メタデータ取得 (GB形式XML)
                            # FASTAとメタデータを別々に取るのは効率が悪いので、GBファイルを取得してパースする
                            handle_gb = Entrez.efetch(db="nucleotide", id=id_list, rettype="gb", retmode="xml")
                            gb_records = Entrez.parse(handle_gb)
                            
                            count = 0
                            for rec in gb_records:
                                # 基本情報
                                accession = rec.get("GBSeq_primary-accession", "")
                                definition = rec.get("GBSeq_definition", "")
                                length = rec.get("GBSeq_length", "")
                                sequence = rec.get("GBSeq_sequence", "").upper()
                                create_date = rec.get("GBSeq_create-date", "")
                                update_date = rec.get("GBSeq_update-date", "")
                                
                                # 論文情報 (最初のReference)
                                refs = rec.get("GBSeq_references", [])
                                journal = ""
                                authors = ""
                                title = ""
                                if refs:
                                    first_ref = refs[0]
                                    journal = first_ref.get("GBReference_journal", "")
                                    title = first_ref.get("GBReference_title", "")
                                    auth_list = first_ref.get("GBReference_authors", [])
                                    authors = ", ".join(auth_list) if auth_list else ""

                                # 特徴テーブル (Sourceから採取地情報を探す)
                                country = ""
                                lat_lon = ""
                                collection_date = ""
                                collector = ""
                                isolation_source = ""
                                geo_loc_name = ""
                                
                                features = rec.get("GBSeq_feature-table", [])
                                for feat in features:
                                    if feat["GBFeature_key"] == "source":
                                        quals = feat.get("GBFeature_quals", [])
                                        for q in quals:
                                            if q["GBQualifier_name"] == "country": country = q["GBQualifier_value"]
                                            if q["GBQualifier_name"] == "lat_lon": lat_lon = q["GBQualifier_value"]
                                            if q["GBQualifier_name"] == "collection_date": collection_date = q["GBQualifier_value"]
                                            if q["GBQualifier_name"] == "collected_by": collector = q["GBQualifier_value"]
                                            if q["GBQualifier_name"] == "isolation_source": isolation_source = q["GBQualifier_value"]
                                            if q["GBQualifier_name"] == "geo_loc_name": geo_loc_name = q["GBQualifier_value"]

                                # FASTA用レコード作成
                                clean_sp = sp.replace(" ", "_")
                                seq_id = f"{clean_sp}_{accession}"
                                fasta_records.append(f">{seq_id} {definition}\n{sequence}\n")
                                
                                # メタデータリスト追加
                                metadata_list.append({
                                    "Sequence_ID": seq_id,
                                    "Species": sp,
                                    "Accession": accession,
                                    "Definition": definition,
                                    "Length": length,
                                    "Country": country,
                                    "Geo_Loc_Name": geo_loc_name,
                                    "Isolation_Source": isolation_source,
                                    "Lat_Lon": lat_lon,
                                    "Collection_Date": collection_date,
                                    "Collector": collector,
                                    "Authors": authors,
                                    "Journal": journal,
                                    "Title": title,
                                    "Create_Date": create_date,
                                    "Update_Date": update_date
                                })
                                count += 1
                                
                            log_text += f"✅ {sp}: {count}件取得\n"
                            time.sleep(0.5) # API制限回避
                            
                        except Exception as e:
                            log_text += f"⚠️ {sp}: Error {e}\n"
                    
                    status_text.success("完了！")
                    st.text_area("ログ", log_text, height=100)
                    
                    if fasta_records:
                        # FASTAダウンロード
                        fasta_str = "".join(fasta_records)
                        st.download_button("📥 FASTAをダウンロード", fasta_str, "sequences.fasta")
                        
                        # Excelダウンロード
                        if metadata_list:
                            df_meta = pd.DataFrame(metadata_list)
                            excel_buffer = BytesIO()
                            with pd.ExcelWriter(excel_buffer, engine='openpyxl') as writer:
                                df_meta.to_excel(writer, index=False, sheet_name='Metadata')
                            
                            st.download_button(
                                "📥 メタデータ(Excel)をダウンロード",
                                excel_buffer.getvalue(),
                                "metadata.xlsx",
                                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                            )
                    else:
                        st.warning("データが見つかりませんでした。絞り込み条件が厳しすぎる可能性があります。")
                        
                except Exception as e:
                    st.error(f"エラーが発生しました: {e}")
        elif uploaded_file:
            st.warning("メールアドレスを入力してください。")
