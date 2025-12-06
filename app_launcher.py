# app_launcher.py
import streamlit as st
from modules.downloader import app_downloader
from modules.phylo import app_phylo
from modules.waveform_ui import app_waveform_main

# --- メイン設定 ---
st.set_page_config(page_title="PhyloPipeline Suite Ultimate", layout="wide", page_icon="🧬")

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

# --- モード切替実行 ---
if "Waveform" in app_mode:
    app_waveform_main()
elif "Downloader" in app_mode:
    app_downloader()
elif "Phylo" in app_mode:
    app_phylo()
