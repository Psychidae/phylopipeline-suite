import numpy as np
from PIL import Image
import streamlit as st

# Color Schemes (RGB)
COLORS = {
    'Default': {
        'A': (255, 68, 68),  # Red
        'T': (68, 255, 68),  # Green
        'G': (255, 215, 0),  # Gold/Yellow
        'C': (68, 68, 255),  # Blue
        '-': (224, 224, 224), # Grey
        'N': (128, 128, 128), # Dark Grey
        '?': (255, 255, 255)  # White
    },
    'Nucleotide': { # Similar to AliView
        'A': (230, 10, 10),
        'C': (10, 10, 230),
        'G': (230, 230, 0),
        'T': (10, 230, 10),
        '-': (240, 240, 240),
        'N': (128, 128, 128)
    }
}

# Pre-compute integer map for faster lookup
def _create_color_lut(scheme_name='Default'):
    scheme = COLORS.get(scheme_name, COLORS['Default'])
    # Create an array mapping ASCII values (0-127) to RGB
    lut = np.full((128, 3), 255, dtype=np.uint8) # Default white
    
    for char, rgb in scheme.items():
        # Handle both upper and lower case
        lut[ord(char.upper())] = rgb
        lut[ord(char.lower())] = rgb
    
    return lut

def render_overview_image(alignment, width=None, height=150, color_scheme='Default'):
    """
    Generate a PIL Image for the alignment overview (minimap).
    
    Args:
        alignment: Bio.Align.MultipleSeqAlignment object
        width: Output width (optional). If None, uses alignment length (1px/base).
        height: Output height (fixed for UI).
        color_scheme: Name of color scheme.
        
    Returns:
        PIL.Image
    """
    if not alignment:
        return Image.new('RGB', (100, 20), color='white')

    # 1. Convert Alignment to Numpy Array (Characters)
    # Using np.array([list(rec) for rec in alignment]) is standard but can be slow for huge alignments
    # Optimization: Read directly into bytearray if possible, but BioPython objects are distinct.
    # We'll stick to the list comprehension for now as it's reasonably fast for <10MB.
    
    # Extract sequences as strings
    seqs = [str(record.seq) for record in alignment]
    n_seq = len(seqs)
    seq_len = len(seqs[0]) if n_seq > 0 else 0
    
    if n_seq == 0 or seq_len == 0:
        return Image.new('RGB', (100, 20), color='white')

    # Convert to 2D array of ASCII codes (literals)
    # Viewing as uint8 (ASCII)
    # Note: Pad sequences if lengths differ? (BioPython alignment assumes equal length)
    char_array = np.array([list(s) for s in seqs], dtype='U1')
    # Convert unicode chars to ascii integers for LUT
    # .view(np.int32) is for UCS4, we want simple byte mapping if possible.
    # Safer: use ord vectorization or just loop setup.
    # Fast approach:
    # ord_map = {c: ord(c) for c in 'ATGCNRYSWKMBDHV.-?'}
    # But numpy handles View casting ok if we encode to 'S1' (bytes)
    
    byte_array = np.array(seqs, dtype='S'+str(seq_len)) # Array of bytestrings
    # This creates 1D array of strings. We need 2D bytes.
    # Actually, fast conversion:
    # stack rows
    
    # Robust method:
    # 1. Create mapping array
    lut = _create_color_lut(color_scheme)
    
    # 2. Iterate and fill buffer (or use frombuffer implementation)
    # Since alignment constructs are list of SeqRecord, let's just build the buffer.
    # Join all rows into one giant bytes object
    full_bytes = b"".join([s.encode('ascii', errors='replace') for s in seqs])
    
    # Convert to numpy uint8 array
    int_view = np.frombuffer(full_bytes, dtype=np.uint8)
    
    # 3. Apply LUT (Vectorized Color Mapping)
    # int_view contains ASCII codes. lut is indexed by ASCII code.
    # This is extremely fast.
    try:
        rgb_data = lut[int_view]
    except IndexError:
        # Fallback for non-ascii characters > 127
        # Clip to 127 or handle error
        int_view = np.clip(int_view, 0, 127)
        rgb_data = lut[int_view]
    
    # Reshape to (rows, cols, 3)
    rgb_image = rgb_data.reshape((n_seq, seq_len, 3))
    
    # 4. Create Image
    img = Image.fromarray(rgb_image, mode='RGB')
    
    # 5. Resize if requested (e.g. to fit UI width)
    # Nearest neighbor is best for preserving sharp pixels, but Bilinear looks smoother for overview.
    # AliView uses distinct pixels.
    # If width is specified, resize.
    if width and width != seq_len:
        img = img.resize((width, height), resample=Image.NEAREST)
    else:
        # If no width specified, keep native width but scale height?
        # Usually we want to fit the container.
        # Let's just return the native image for logic flexibility.
        pass

    return img
