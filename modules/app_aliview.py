import streamlit as st
import pandas as pd
from Bio import AlignIO
from io import StringIO
from modules.aliview_model import load_alignment_data
from modules.aliview_renderer import render_overview_image

def app_aliview():
    st.header("🔄 AliView Integration Bridge")
    st.info("Desktop版 **AliView** で編集を行うための連携モードです。")

    # --- 1. Select Target File ---
    session_files = list(st.session_state.get("file_manager", {}).keys())
    
    if not session_files:
        st.warning("No files loaded in the session. Please go to 'Alignment Viewer' or 'GenBank Downloader' to load data first.")
        # Optional: Allow quick upload here too?
        uploaded_new = st.file_uploader("Or upload a new file to start:", type=["fasta", "fas", "txt"])
        if uploaded_new:
            # Add to session logic would be replicated... maybe just show it.
            # For bridge simplicity, assume usage *within* a workflow.
            pass
        return

    selected_filename = st.selectbox("1. Select File to Edit:", session_files)
    
    # Get current data (DataFrame)
    df_current = st.session_state.file_manager[selected_filename]
    
    # Convert DataFrame to FASTA string (for download/preview)
    # Assumes df has 'id' and 'sequence'
    fasta_str_current = ""
    try:
        if 'id' in df_current.columns and 'sequence' in df_current.columns:
            recs = []
            for _, row in df_current.iterrows():
                recs.append(f">{row['id']}\n{row['sequence']}")
            fasta_str_current = "\n".join(recs)
    except Exception as e:
        st.error(f"Error reading file data: {e}")
        return

    # --- 2. Preview & Download ---
    col1, col2 = st.columns([2, 1])
    
    with col1:
        st.subheader("Current State")
        # Generate Preview Image (Virtual Alignment Object needed for renderer)
        # aliview_renderer needs Bio.MultipleSeqAlignment
        # create one from string
        try:
            alignment_current = AlignIO.read(StringIO(fasta_str_current), "fasta")
            img_current = render_overview_image(alignment_current, width=600, height=80, color_scheme='Default')
            st.image(img_current, use_container_width=True, caption=f"Before Edit ({len(alignment_current)} seqs, {alignment_current.get_alignment_length()} bp)")
        except Exception as e:
            st.warning("Could not render preview image.")

    with col2:
        st.subheader("Action")
        st.download_button(
            label="⬇️ Download FASTA",
            data=fasta_str_current,
            file_name=f"for_aliview_{selected_filename}",
            mime="text/plain",
            type="primary",
            help="Download this file, open it in AliView, edit, and save."
        )
        st.markdown("""
        **How to edit:**
        1. Download file.
        2. Open in **AliView**.
        3. Edit (delete seqs, trim).
        4. **File > Save as Fasta**.
        """)

    st.markdown("---")

    # --- 3. Upload & Update ---
    st.subheader("2. Import Edited File")
    
    uploaded_edited = st.file_uploader("Upload the file saved from AliView:", type=["fasta", "fas", "txt"], key="aliview_import")
    
    if uploaded_edited:
        # Load and Preview
        alignment_edited = load_alignment_data(uploaded_edited, uploaded_edited.name)
        
        if alignment_edited:
            st.subheader("New State Preview")
            img_edited = render_overview_image(alignment_edited, width=600, height=80, color_scheme='Default')
            st.image(img_edited, use_container_width=True, caption=f"After Edit ({len(alignment_edited)} seqs, {alignment_edited.get_alignment_length()} bp)")
            
            # Update Button
            if st.button("🔄 Update Session with Edited Data", type="primary"):
                # Convert Alignment back to DataFrame for Session Manager
                new_rows = []
                for record in alignment_edited:
                    new_rows.append({
                        "id": record.id,
                        "sequence": str(record.seq),
                        "length": len(record.seq)
                    })
                df_new = pd.DataFrame(new_rows)
                
                # UPDATE Session State
                # We overwrite the *original* key with new data? Or create new?
                # User flow suggests updating the current file is mostly what corresponds to "Save".
                st.session_state.file_manager[selected_filename] = df_new
                
                # Also update active viewer_df if it matches
                st.session_state.viewer_df = df_new.copy()
                
                st.success(f"Updated '{selected_filename}' successfully! You can now proceed to PhyloPipeline.")
                st.balloons()

