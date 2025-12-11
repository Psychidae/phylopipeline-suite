import streamlit as st
import sys
import os

# Add project root to sys.path
root_dir = os.path.join(os.path.dirname(__file__), '..')
sys.path.append(root_dir)

# --- Install Dependencies (Linux/Cloud) ---
installer_log = []
try:
    from modules.installer import install_tools
    # Capture print output from installer (simple mock)
    # Ideally installer should return logs, but for now we trust it works or fails
    install_tools()
    installer_msg = "Installer ran successfully."
except Exception as e:
    installer_msg = f"Installer failed: {e}"

# Add tools directory to PATH (Must be AFTER installer creates it)
tools_dir = os.path.join(root_dir, 'tools')
if os.path.exists(tools_dir):
    os.environ["PATH"] += os.pathsep + os.path.abspath(tools_dir)

from modules.downloader import app_downloader
from modules.phylo import app_phylo
from modules.waveform_ui import app_waveform_main

# --- アプリ設定 ---
st.set_page_config(page_title="PhyloPipeline Suite", layout="wide", page_icon="🧬")

# --- DEBUG INFO ---
with st.expander("🛠 System Diagnostics (Debug)", expanded=True):
    import platform
    import shutil
    st.write(f"**OS:** {platform.system()} {platform.release()} ({platform.machine()})")
    st.write(f"**Installer Status:** {installer_msg}")
    
    tools_dir_debug = os.path.join(os.path.dirname(__file__), "..", "tools")
    st.write(f"**Tools Dir:** `{os.path.abspath(tools_dir_debug)}` (Exists: {os.path.exists(tools_dir_debug)})")
    
    if os.path.exists(tools_dir_debug):
        st.write(f"**Tools Content:** {os.listdir(tools_dir_debug)}")
        iq_path = os.path.join(tools_dir_debug, "iqtree2")
        st.write(f"**iqtree2 exists:** {os.path.exists(iq_path)}")
        if os.path.exists(iq_path):
             st.write(f"**iqtree2 permissions:** {oct(os.stat(iq_path).st_mode)[-3:]}")
    
    st.write(f"**PATH:** `{os.environ.get('PATH')}`")
    
    # Try finding tool
    from modules.common import find_tool_path
    found_path = find_tool_path("iqtree2")
    st.write(f"**find_tool_path('iqtree2'):** `{found_path}`")
    
# ------------------

# --- サイドバーメニュー ---
st.sidebar.title("🧬 PhyloPipeline")
st.sidebar.caption("Integrated Analysis Suite")

app_mode = st.sidebar.radio(
    "Select Mode",
    ["1. Waveform Validator (波形解析)", 
     "2. GenBank Downloader (配列取得)", 
     "3. PhyloPipeline (系統解析)"]
)

st.sidebar.markdown("---")

# --- モード切替 ---
if "Waveform" in app_mode:
    app_waveform_main()
elif "Downloader" in app_mode:
    app_downloader()
elif "Phylo" in app_mode:
    app_phylo()
