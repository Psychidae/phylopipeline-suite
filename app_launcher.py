import streamlit as st
from modules.downloader import app_downloader
from modules.phylo import app_phylo
# 波形解析は既存のものをそのまま利用
from modules.waveform_ui import app_waveform_main
from modules.app_viewer import app_viewer

# --- メイン設定 ---
st.set_page_config(page_title="PhyloPipeline Suite Ultimate", layout="wide", page_icon="🧬")

# --- サイドバーメニュー ---
st.sidebar.title("🧬 PhyloPipeline")
st.sidebar.caption("Integrated Analysis Suite")

app_mode = st.sidebar.radio(
    "Select Mode",
    ["1. Waveform Validator (波形解析)",
     "2. GenBank Downloader (配列取得)",
     "3. PhyloPipeline (系統解析)",
     "4. Alignment Viewer (アライメント表示)"]
)

st.sidebar.markdown("---")

# --- モード切替 ---
if "Waveform Validator" in app_mode:
    app_waveform_main()
elif "GenBank Downloader" in app_mode:
    app_downloader()
elif "Phylo" in app_mode:
    app_phylo()
elif "Alignment Viewer" in app_mode:
    app_viewer()
