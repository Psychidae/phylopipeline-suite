import streamlit as st
import pandas as pd
from Bio import SeqIO
from io import StringIO
from modules.phylo_editor import open_alignment_editor

# Re-use the editor logic but embedded in the main page (or launch dialog immediately)
# Since PhyloPipeline uses dialog for editor, we can use a button to launch it
# OR we can extract the plotting logic from dialog to a component?
# The user wants "Editor Functionality" + "Save".
# open_alignment_editor is a @st.dialog.
# We can just call it? Yes.

def app_viewer():
    st.header("🧬 Alignment Viewer")
    st.info("FASTAファイルをアップロードして、アラインメントの確認・編集・保存を行います。複数ファイルの管理が可能です。")
    
    # --- 1. Session Manager ---
    if "file_manager" not in st.session_state:
        st.session_state.file_manager = {} # {filename: df}
    
    # Upload
    uploaded_files = st.file_uploader("Upload FASTA (Multiple allowed)", type=["fasta", "fas", "aln"], accept_multiple_files=True, key="viewer_uploader")
    
    if uploaded_files:
        for f in uploaded_files:
            if f.name not in st.session_state.file_manager:
                try:
                    stringio = StringIO(f.getvalue().decode("utf-8"))
                    recs = list(SeqIO.parse(stringio, "fasta"))
                    if recs:
                        df = pd.DataFrame([{"Include":True, "ID":s.id, "Sequence":str(s.seq)} for s in recs])
                        st.session_state.file_manager[f.name] = df
                        # 新規ファイルを自動選択するためにキーをセットしたいが、selectboxのindex制御が必要
                        # 簡易的に、最後にアップロードされたものを記憶する変数を設定可能だが、今回は手動選択に任せる
                except Exception as e:
                    st.error(f"Error reading {f.name}: {e}")

    # Cross-tool Import (PhyloPipeline)
    if "phylo_aligned_df" in st.session_state and st.session_state.phylo_aligned_df is not None:
        if st.sidebar.button("Import from PhyloPipeline"):
            name = "Phylo_Result.aln"
            st.session_state.file_manager[name] = st.session_state.phylo_aligned_df.copy()
            st.success(f"Imported as {name}")

    # --- 2. File Selection & Management ---
    if not st.session_state.file_manager:
        st.warning("No files loaded.")
        return

    file_list = list(st.session_state.file_manager.keys())
    
    # 選択状態の管理
    if "viewer_selected_file" not in st.session_state:
        st.session_state.viewer_selected_file = file_list[0]
    
    # 削除操作などで選択中のファイルが消えた場合のケア
    if st.session_state.viewer_selected_file not in file_list:
        st.session_state.viewer_selected_file = file_list[0]

    col_sel, col_del = st.columns([3, 1])
    with col_sel:
        selected_file = st.selectbox("Select File", file_list, key="file_selector", index=file_list.index(st.session_state.viewer_selected_file))
        st.session_state.viewer_selected_file = selected_file
    
    with col_del:
        if st.button("🗑️ Delete"):
            del st.session_state.file_manager[selected_file]
            st.rerun()

    # --- 3. Sync Logic ---
    # 選択されたファイルのDFを `viewer_df` にロード（エディタ用）
    # エディタで編集されると `viewer_df` が更新されるので、それを `file_manager` に書き戻す
    
    # ロード: 以前のループと違うファイルが選択されたらロード
    # しかし `viewer_df` と `manager` の同期をどう保つか？
    # 単純化: 毎回 `viewer_df` を `file_manager` からコピーし、エディタ保存時は `rerun` されるので
    # 保存後にここに来たときは `viewer_df` が新しい。
    # したがって、「ファイル切り替え時」のみロードし、それ以外は `viewer_df` を正とする？
    
    # State管理
    if "current_viewing_file_name" not in st.session_state:
        st.session_state.current_viewing_file_name = None
    
    # ファイル切り替え検知
    if st.session_state.current_viewing_file_name != selected_file:
        st.session_state.viewer_df = st.session_state.file_manager[selected_file].copy()
        st.session_state.current_viewing_file_name = selected_file
    else:
        # 同じファイルを見ている場合、viewer_dfが更新されている可能性があるのでManagerに書き戻す
        # (EditorでSave -> Rerun -> ここに来る -> viewer_dfは更新済み)
        if "viewer_df" in st.session_state:
            st.session_state.file_manager[selected_file] = st.session_state.viewer_df.copy()

    # --- 4. Main Viewer UI ---
    if "viewer_df" in st.session_state:
        df = st.session_state.viewer_df
        st.caption(f"Editing: {selected_file} ({len(df)} seqs)")
        
        if st.button("Open Alignment Editor", type="primary"):
            open_alignment_editor(df, target_key="viewer_df")
        
        st.divider()
        st.subheader("Data Preview")
        st.dataframe(df, height=200)
        
        # Download
        st.subheader("Download Selected")
        final_df = df[df["Include"] == True]
        fasta_str = ""
        for _, row in final_df.iterrows():
            fasta_str += f">{row['ID']}\n{row['Sequence']}\n"
        
        st.download_button("Download FASTA", fasta_str, f"edited_{selected_file}")

def open_alignment_editor_wrapper(df):
    open_alignment_editor(df)
