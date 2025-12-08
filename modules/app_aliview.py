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

    # --- 3. Controls (Zoom & Scroll) ---
    st.markdown("### Controls")
    c1, c2 = st.columns([1, 1])
    
    # Zoom Level (1: Far view, 10: Close up)
    zoom_level = c1.slider("🔍 Zoom Level", min_value=1, max_value=10, value=5, help="Change character size and window width")
    
    # Calculate grid parameters based on zoom
    # Zoom 10: 20 bases, 40px width, 18px font
    # Zoom 1: 200 bases, 10px width, 0px font (color only)
    
    # Linear interpolation or specific steps
    base_window_map = {
        10: 20, 9: 30, 8: 40, 7: 50, 6: 60, 
        5: 80, 4: 100, 3: 120, 2: 150, 1: 200
    }
    width_map = {
        10: 40, 9: 35, 8: 30, 7: 25, 6: 20, 
        5: 18, 4: 15, 3: 12, 2: 10, 1: 8
    }
    font_map = {
        10: 18, 9: 16, 8: 14, 7: 12, 6: 12, 
        5: 11, 4: 10, 3: 0, 2: 0, 1: 0
    }
    
    window_size = base_window_map[zoom_level]
    col_width = width_map[zoom_level]
    font_size = font_map[zoom_level]
    
    seq_len = alignment.get_alignment_length()
    max_start = max(0, seq_len - window_size)
    
    # Scroll Slider
    start_pos = st.slider("↔️ Scroll Position", 0, max_start, 0)
    end_pos = start_pos + window_size
    
    st.caption(f"Showing bases **{start_pos+1}** to **{end_pos}** (Total: {seq_len})")

    # --- 4. Detail View (AgGrid) ---
    st.subheader("Detail View")
    
    df_window = prepare_aggrid_data(_alignment=alignment, start=start_pos, end=end_pos)
    
    # Configure Grid
    gb = GridOptionsBuilder.from_dataframe(df_window)
    # ID Column Pinned and WIDER
    gb.configure_column("ID", pinned="left", width=200, resizable=True, lockPosition=True)
    
    # JsCode for coloring
    # Inject font-size dynamically
    cell_style_jscode = JsCode(f"""
    function(params) {{
        var fontSize = '{font_size}px';
        var color = 'black';
        var bg = 'white';
        
        if (params.value == 'A') {{ bg = '#FF4444'; }}
        else if (params.value == 'T') {{ bg = '#44FF44'; }}
        else if (params.value == 'G') {{ bg = '#FFD700'; }}
        else if (params.value == 'C') {{ bg = '#4444FF'; color = 'white'; }}
        else if (params.value == '-') {{ bg = '#E0E0E0'; }}
        
        if ({font_size} == 0) {{ color = 'rgba(0,0,0,0)'; }} // Hide text
        
        return {{'backgroundColor': bg, 'color': color, 'textAlign': 'center', 'fontSize': fontSize, 'fontWeight': 'bold'}};
    }}
    """)
    
    # Apply style to all base columns (1, 2, 3...)
    # Columns are strings in df
    base_cols = [c for c in df_window.columns if c != "ID"]
    for c in base_cols:
        gb.configure_column(c, width=col_width, cellStyle=cell_style_jscode, sortable=False, suppressMenu=True)
    
    # Global Options
    gb.configure_grid_options(
        domLayout='autoHeight', 
        headerHeight=30,
        rowHeight=25,
        suppressMovableColumns=True
    )
    
    gridOptions = gb.build()
    
    # Use key to force redraw when zoom changes to avoid column width glitches
    grid_key = f"grid_{zoom_level}_{start_pos}"
    
    with st.spinner("Rendering grid..."):
        AgGrid(
            df_window,
            gridOptions=gridOptions,
            allow_unsafe_jscode=True,
            height=500,
            theme="alpine",
            key=grid_key
        )
