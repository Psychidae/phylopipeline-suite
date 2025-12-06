import streamlit as st
import pandas as pd
import numpy as np

# Try importing st_aggrid, show warning if missing
try:
    from st_aggrid import AgGrid, GridOptionsBuilder, JsCode
    HAS_AGGRID = True
    IMPORT_ERR = ""
except ImportError as e:
    HAS_AGGRID = False
    IMPORT_ERR = str(e)

def app_aliview():
    st.header("🧬 AliView Prototype (AgGrid)")
    
    if not HAS_AGGRID:
        st.error(f"`streamlit-aggrid` library is missing. Please run `pip install streamlit-aggrid`.\nError detail: {IMPORT_ERR}")
        return

    st.info("Performance test using AgGrid + Client-side JS Rendering.")

    # 1. Create Dummy Data
    # Short lines x 20 bases
    data = [
        {"ID": "Seq1", "Seq": "ATGCATGCATGCATGCATGC"},
        {"ID": "Seq2", "Seq": "ATGCATGCAT-CATGCATGC"},
        {"ID": "Seq3", "Seq": "ATGCATGCATGCATGCATGG"},
        {"ID": "Seq4", "Seq": "AAACATGCATGCATGCATGC"},
        {"ID": "Seq5", "Seq": "ATGCATGCAT-CATGCATGC"},
    ]
    
    # AgGrid handles columns. For alignment, we usually split sequence into columns? 
    # OR we treat the whole sequence as one string and use HTML/JS to render it?
    # AgGrid is a grid. AliView is a grid of BASES.
    # So we should split "Seq" into "Pos1", "Pos2"... or just use a custom cell renderer for the whole string?
    # The user prompt said: "Column Width column width ... 25px".
    # This implies 1 Base = 1 Column.
    
    # Transform to Base Columns
    rows = []
    for d in data:
        row = {"ID": d["ID"]}
        for i, char in enumerate(d["Seq"]):
            row[f"{i+1}"] = char
        rows.append(row)
    
    df = pd.DataFrame(rows)

    # 2. Configure Grid
    gb = GridOptionsBuilder.from_dataframe(df)
    
    # ID Column Pinned
    gb.configure_column("ID", pinned="left", width=120)
    
    # Base Columns Configuration
    # Define JS for styling
    # A=Red, T=Green, G=Yellow, C=Blue, Gap=Grey
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
    
    # Apply to all base columns (keys that are digits)
    cols = [c for c in df.columns if c != "ID"]
    for c in cols:
        gb.configure_column(c, width=25, cellStyle=cell_style_jscode)
        
    gb.configure_grid_options(domLayout='normal') # Scrollable
    gridOptions = gb.build()
    
    # 3. Render
    st.write("### Mockup View")
    AgGrid(
        df,
        gridOptions=gridOptions,
        allow_unsafe_jscode=True, # Need this for JS colors
        height=400,
        theme="alpine" # clean theme
    )

