import streamlit as st
import streamlit.components.v1 as components
from io import StringIO

def app_aliview():
    st.header("🧬 AliView (BioJS Mode)")
    st.info("AliViewのような操作感を実現するため、BioJS MSA Viewerを採用しました。マウスでスクロール・ズームが可能です。")

    # --- 1. Data Selection ---
    session_files = list(st.session_state.get("file_manager", {}).keys())
    
    col1, col2 = st.columns([1, 2])
    fasta_str = ""
    
    with col1:
        if session_files:
            selected_file = st.selectbox("Select File", ["(Upload New)"] + session_files)
        else:
            selected_file = "(Upload New)"
            
        if selected_file == "(Upload New)":
            uploaded_file = st.file_uploader("Upload FASTA", type=["fasta", "fas", "txt"], key="biojs_upload")
            if uploaded_file:
                fasta_str = uploaded_file.getvalue().decode("utf-8")
        else:
            # Retrieve from session (DataFrame -> FASTA String conversion needed?)
            # The file manager stores DataFrames. We need to convert back or store raw str.
            # Ideally, we should check if we can reconstruct it.
            # For robustness in this prototype phase, let's ask for upload if we can't easily get str.
            # BUT, we can convert simple dataframe back to fasta.
            try:
                df = st.session_state.file_manager[selected_file]
                # Assuming standard columns 'id', 'seq' or similar from our parser
                # Let's check columns.
                # If dataframe format is known (from app_viewer.py), it has 'id' and 'sequence'.
                if 'id' in df.columns and 'sequence' in df.columns:
                    recs = []
                    for _, row in df.iterrows():
                        recs.append(f">{row['id']}\n{row['sequence']}")
                    fasta_str = "\n".join(recs)
                else:
                    st.warning("Selected file format is not compatible. Please upload raw FASTA.")
            except Exception as e:
                st.error(f"Error loading file: {e}")

    if not fasta_str:
        st.info("👈 Please upload or select a file to view.")
        return

    # --- 2. Render BioJS MSA ---
    # We inject the FASTA string into the JS.
    # Escape newlines for JS string safety.
    clean_fasta = fasta_str.replace("\n", "\\n").replace("'", "\\'")
    
    # BioJS MSA Template
    # Using msa@1.0.3 (Stable)
    html_code = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <script src="https://cdn.jsdelivr.net/npm/msa@1.0.3/dist/msa.min.js"></script>
    </head>
    <body>
        <div id="msa_menu"></div>
        <div id="msa_div"></div>
        
        <script>
            var sequences = [
                {{name: "seq1", seq: "ACGT..."}} // Placeholder logic? No, we use msa.io.fasta
            ];
            
            var fastaData = '{clean_fasta}';
            
            var opts = {{
                el: document.getElementById("msa_div"),
                vis: {{
                    conserv: false,
                    overviewbox: true,
                    seqlogo: true
                }},
                conf: {{
                    dropImport: true
                }},
                zoomer: {{
                    menuFontsize: "12px",
                    autoResize: true
                }}
            }};
            
            var m = new msa.msa(opts);
            
            // Import FASTA
            var seqs = msa.io.fasta.parse(fastaData);
            m.seqs.reset(seqs);
            m.render();
            
            // Adjust height
            // m.el.style.height = "500px"; 
        </script>
    </body>
    </html>
    """
    
    # Refined Template with full width and proper importing
    html_full = f"""
    <script src="https://cdn.jsdelivr.net/npm/msa@1.0.3/dist/msa.min.js"></script>
    <div id="msa_app"></div>
    <script>
        var fasta = '{clean_fasta}';
        var opts = {{
            el: document.getElementById('msa_app'),
            vis: {{
                conserv: true,
                overviewbox: true,
                seqlogo: true,
                labelId: true
            }},
            zoomer: {{
                boxRectHeight: 1,
                boxRectWidth: 1,
                alignmentHeight: 600,
                labelNameLength: 200
            }},
            colorscheme: {{
                scheme: "nucleotide"
            }}
        }};
        var m = new msa.msa(opts);
        var seqs = msa.io.fasta.parse(fasta);
        m.seqs.reset(seqs);
        m.render();
    </script>
    """

    components.html(html_full, height=700, scrolling=True)
