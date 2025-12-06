import streamlit as st
import streamlit as st

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
     "4. Alignment Viewer (アライメント表示)",
     "5. AliView Prototype (実験機能)"]
)

st.sidebar.markdown("---")

# --- モード切替 ---
if "Waveform Validator" in app_mode:
    try:
        from modules.waveform_ui import app_waveform_main
        app_waveform_main()
    except Exception as e:
        st.error(f"Failed to load Waveform Validator: {e}")

elif "GenBank Downloader" in app_mode:
    try:
        from modules.downloader import app_downloader
        app_downloader()
    except Exception as e:
        st.error(f"Failed to load Downloader: {e}")

elif "PhyloPipeline" in app_mode:
    try:
        from modules.phylo import app_phylo
        app_phylo()
    except Exception as e:
        st.error(f"Failed to load PhyloPipeline: {e}")

elif "Alignment Viewer" in app_mode:
    try:
        from modules.app_viewer import app_viewer
        app_viewer()
    except Exception as e:
        st.error(f"Failed to load Viewer: {e}")

elif "AliView Prototype" in app_mode:
    try:
        from modules.app_aliview import app_aliview
        app_aliview()
    except Exception as e:
        st.error(f"Failed to load AliView Prototype: {e}")
