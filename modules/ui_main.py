import streamlit as st
import os

from modules.downloader import app_downloader
from modules.phylo import app_phylo
from modules.waveform_ui import app_waveform_main


def main():
    # set_page_config は app.py で呼ばれているためここでは呼ばない

    # --- サイドバーメニュー ---
    st.sidebar.title("🧬 PhyloPipeline")
    st.sidebar.caption("Integrated Analysis Suite")

    app_mode = st.sidebar.radio(
        "Select Mode",
        [
            "1. Waveform Validator (波形解析)",
            "2. GenBank Downloader (配列取得)",
            "3. PhyloPipeline (系統解析)",
        ]
    )

    st.sidebar.markdown("---")

    # --- モード切替 ---
    if "Waveform" in app_mode:
        app_waveform_main()
    elif "Downloader" in app_mode:
        app_downloader()
    elif "Phylo" in app_mode:
        app_phylo()
