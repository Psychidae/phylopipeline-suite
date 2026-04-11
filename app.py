import os
import sys
import streamlit as st
from modules.constants import TOOLS_DIR

# --- 1. ページ設定（必ず最初のStreamlitコマンドとして呼ぶ）---
st.set_page_config(page_title="PhyloPipeline Suite", layout="wide", page_icon="🧬")

# --- 2. ツールディレクトリを PATH に追加 ---
if os.path.exists(TOOLS_DIR):
    os.environ["PATH"] += os.pathsep + os.path.abspath(TOOLS_DIR)

# --- 3. ツールインストール（Linux / クラウド環境のみ） ---
try:
    from modules.installer import install_tools_linux
    success, msg = install_tools_linux()
    if not success:
        st.warning(f"Tool Setup Warning: {msg}")
except Exception as e:
    st.warning(f"Installer could not run: {e}")

# --- 4. アプリケーション起動 ---
try:
    from modules.ui_main import main
    main()
except ImportError as e:
    st.error(f"Failed to load application modules: {e}")
    st.info("プロジェクトルートから実行してください: `streamlit run app.py`")
except Exception as e:
    st.error(f"Application Runtime Error: {e}")
    st.exception(e)
