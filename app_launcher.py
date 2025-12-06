import streamlit as st
from modules.downloader import app_downloader
from modules.phylo import app_phylo
# 波形解析は既存のものをそのまま利用
from modules.waveform_ui import app_waveform_main
from modules.app_viewer import app_viewer
from modules.app_aliview import app_aliview # Added import for app_aliview
# Assuming app_settings and app_help exist or will be added
# from modules.app_settings import app_settings
# from modules.app_help import app_help

# --- メイン設定 ---
st.set_page_config(page_title="PhyloPipeline Suite Ultimate", layout="wide", page_icon="🧬")

# --- サイドバーメニュー ---
st.sidebar.title("🧬 PhyloPipeline")
st.sidebar.caption("Integrated Analysis Suite")

# Replaced the old app_mode radio and if/elif block with the new structure
tool_mode = st.sidebar.radio(
    "Select Tool:",
    ("PhyloPipeline Pro", "Alignment Viewer", "AliView Prototype", "Settings", "Help / Walkthrough")
)

st.sidebar.markdown("---")

# --- モード切替 ---
if tool_mode == "PhyloPipeline Pro":
    app_phylo()
elif tool_mode == "Alignment Viewer":
    app_viewer()
