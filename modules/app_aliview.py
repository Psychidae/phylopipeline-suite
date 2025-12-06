import streamlit as st
import pandas as pd
from Bio import AlignIO
from io import StringIO
from modules.aliview_model import load_alignment_data, prepare_aggrid_data
from modules.aliview_renderer import render_overview_image

# Try importing st_aggrid
try:
    from st_aggrid import AgGrid, GridOptionsBuilder, JsCode
    HAS_AGGRID = True
except ImportError:
    HAS_AGGRID = False

def app_aliview():
    st.header("🧬 AliView (Web Edition)")
    
    if not HAS_AGGRID:
        st.error("`streamlit-aggrid` library is missing.")
        return

    # --- 1. Data Loading ---
    # Check if we have files in session from the main viewer
    session_files = list(st.session_state.get("file_manager", {}).keys())
    
    col1, col2 = st.columns([1, 2])
    with col1:
        if session_files:
            selected_file = st.selectbox("Select File from Session", ["(Upload New)"] + session_files)
        else:
            selected_file = "(Upload New)"
            
        if selected_file == "(Upload New)":
            uploaded_file = st.file_uploader("Upload FASTA", type=["fasta", "fas", "txt"])
            alignment = None
            if uploaded_file:
                alignment = load_alignment_data(uploaded_file, uploaded_file.name)
        else:
            # Load from session manager (it stores DataFrames, not Alignments... needs conversion or re-parsing)
            # Actually, reusing the raw file content or dataframe is tricky if we need Alignment object.
            # Simplified: Just ask for upload for this "Prototype" phase to ensure stability.
            st.info("Currently, please upload FASTA directly for this experimental mode.")
            uploaded_file = st.file_uploader("Upload FASTA", type=["fasta", "fas", "txt"], key="aliview_upload")
            alignment = None
            if uploaded_file:
                alignment = load_alignment_data(uploaded_file, uploaded_file.name)

    if not alignment:
        st.info("Please upload an alignment file to start.")
        return

    # --- 2. Overview (Minimap) ---
    st.subheader("Overview")
    # Generate image
    with st.spinner("Rendering overview..."):
        # Width: Use container width? Hard to get exact pixel width in Streamlit.
        # Fixed width for now or full length?
        # Let's try full length but resized to fits typical screen: 1000px?
        overview_img = render_overview_image(alignment, width=1000, height=100, color_scheme='Default')
        st.image(overview_img, use_container_width=True, caption="Alignment Overview (Condensed)")

    # --- 3. Navigation (Slider) ---
    seq_len = alignment.get_alignment_length()
    window_size = st.slider("Window Size (bases)", 20, 100, 50)
    start_pos = st.slider("Position", 0, max(0, seq_len - window_size), 0)
    end_pos = start_pos + window_size
    
    st.caption(f"Showing bases {start_pos+1} - {end_pos} (Total: {seq_len})")

    # --- 4. Detail View (AgGrid) ---
    st.subheader("Detail View")
    
    df_window = prepare_aggrid_data(alignment, start=start_pos, end=end_pos)
    
    # Configure Grid
    gb = GridOptionsBuilder.from_dataframe(df_window)
    gb.configure_column("ID", pinned="left", width=150)
    
    # JsCode for coloring
    cell_style_jscode = JsCode("""
    function(params) {
        if (params.value == 'A') { return {'backgroundColor': '#FF4444', 'color': 'black', 'textAlign': 'center'}; }
        if (params.value == 'T') { return {'backgroundColor': '#44FF44', 'color': 'black', 'textAlign': 'center'}; }
        if (params.value == 'G') { return {'backgroundColor': '#FFD700', 'color': 'black', 'textAlign': 'center'}; }
        if (params.value == 'C') { return {'backgroundColor': '#4444FF', 'color': 'white', 'textAlign': 'center'}; }
        if (params.value == '-') { return {'backgroundColor': '#E0E0E0', 'color': 'black', 'textAlign': 'center'}; }
        return {'textAlign': 'center'};
    }
    """)
    
    # Apply style to all base columns (1, 2, 3...)
    # Columns are strings in df
    base_cols = [c for c in df_window.columns if c != "ID"]
    for c in base_cols:
        gb.configure_column(c, width=30, cellStyle=cell_style_jscode, sortable=False)
    
    gb.configure_grid_options(domLayout='autoHeight')
    gridOptions = gb.build()
    
    with st.spinner("Rendering grid..."):
        AgGrid(
            df_window,
            gridOptions=gridOptions,
            allow_unsafe_jscode=True,
            height=500,
            theme="alpine"
        )
