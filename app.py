import os
import sys
import streamlit as st
from modules.constants import TOOLS_DIR

# --- 1. Environment Setup ---
# Add tools directory to PATH so subprocess calls can find binaries (e.g. iqtree2)
if os.path.exists(TOOLS_DIR):
    os.environ["PATH"] += os.pathsep + os.path.abspath(TOOLS_DIR)

# --- 2. Tool Installation (Linux/Cloud) ---
# Check and install iqtree2 if necessary
try:
    from modules.installer import install_tools_linux
    success, msg = install_tools_linux()
    if not success:
        st.error(f"Setup Error: {msg}")
    # Optional: Log success to console or debug log, but keep UI clean
    # print(msg) 
except Exception as e:
    st.error(f"Critical Installer Error: {e}")

# --- 3. Run Application ---
try:
    from modules.ui_main import main
    main()
except ImportError as e:
    st.error(f"Failed to load application modules: {e}")
    st.info("Please ensure you are running this script from the project root: `streamlit run app.py`")
except Exception as e:
    st.error(f"Application Runtime Error: {e}")
