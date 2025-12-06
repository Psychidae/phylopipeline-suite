import pandas as pd
from Bio import AlignIO
from io import StringIO
import streamlit as st

@st.cache_data(show_spinner=False)
def load_alignment_data(file_content, filename):
    """
    Parse uploaded FASTA file and return a MultipleSeqAlignment object.
    Cached to prevent re-parsing on every interaction.
    """
    try:
        stringio = StringIO(file_content.getvalue().decode("utf-8"))
        alignment = AlignIO.read(stringio, "fasta")
        return alignment
    except Exception as e:
        print(f"Error loading alignment: {e}")
        return None

@st.cache_data(show_spinner=False)
def prepare_aggrid_data(alignment, start=0, end=None):
    """
    Convert a slice of alignment to a DataFrame suitable for AgGrid.
    Optimized: Returns a list of dicts directly to avoid Pandas overhead if possible,
    but st-aggrid expects a DataFrame.
    """
    if not alignment:
        return pd.DataFrame()

    num_seqs = len(alignment)
    seq_len = alignment.get_alignment_length()

    if end is None:
        end = min(start + 100, seq_len) # Default view window
    
    # Ensure bounds
    start = max(0, start)
    end = min(seq_len, end)

    # Effective slice
    sliced_alignment = alignment[:, start:end]
    
    rows = []
    for record in sliced_alignment:
        row_data = {"ID": record.id}
        # Expand sequence into columns: "1", "2", "3"... (relative to start)
        # 1-based index for columns
        seq_str = str(record.seq)
        for i, char in enumerate(seq_str):
            col_name = str(start + i + 1)
            row_data[col_name] = char
        rows.append(row_data)
        
    df = pd.DataFrame(rows)
    return df
