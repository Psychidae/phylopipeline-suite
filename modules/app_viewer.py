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
    st.info("FASTAファイルをアップロードして、アラインメントの確認・編集・保存を行います。")
    
    # 1. File Upload
    # Persistence: Use session_state to store file info
    if "viewer_uploaded_file" not in st.session_state:
        st.session_state.viewer_uploaded_file = None
        
    uploaded_file = st.file_uploader("Upload FASTA", type=["fasta", "fas", "aln"], key="viewer_uploader")
    
    if uploaded_file:
        # Check if new file
        if st.session_state.viewer_uploaded_file != uploaded_file:
            stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
            recs = list(SeqIO.parse(stringio, "fasta"))
            df = pd.DataFrame([{"Include":True, "ID":s.id, "Sequence":str(s.seq)} for s in recs])
            st.session_state.viewer_df = df
            st.session_state.viewer_uploaded_file = uploaded_file
    
    # Check if we have data from other tools (PhyloPipeline)
    if "phylo_aligned_df" in st.session_state and st.sidebar.checkbox("Use data from PhyloPipeline", value=False):
         st.session_state.viewer_df = st.session_state.phylo_aligned_df.copy()

    # 2. Display / Edit Button
    if "viewer_df" in st.session_state:
        st.success(f"Loaded {len(st.session_state.viewer_df)} sequences.")
        
        if st.button("Open Alignment Editor", type="primary"):
            open_alignment_editor(st.session_state.viewer_df, target_key="viewer_df")
        
        st.divider()
        st.subheader("Data Preview")
        st.dataframe(st.session_state.viewer_df)
        
        # Download
        st.subheader("Download Selected")
        final_df = st.session_state.viewer_df[st.session_state.viewer_df["Include"] == True]
        fasta_str = ""
        for _, row in final_df.iterrows():
            fasta_str += f">{row['ID']}\n{row['Sequence']}\n"
        
        st.download_button("Download FASTA", fasta_str, "edited_alignment.fasta")

def open_alignment_editor_wrapper(df):
    open_alignment_editor(df)
