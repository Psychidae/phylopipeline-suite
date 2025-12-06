import streamlit as st
from modules.common import generate_alignment_html_from_df

@st.dialog("配列エディタ (Alignment Editor)", width="large")
def open_alignment_editor(initial_df):
    """
    アラインメントデータの確認・編集を行うダイアログ
    """
    st.write("配列の選択（Include）や内容の確認ができます。")
    
    # UIステート管理（ダイアログ内での操作用）
    if 'editor_show_dots' not in st.session_state:
        st.session_state.editor_show_dots = False
        
    # 表示オプション
    show_dots = st.toggle("同一塩基をドット(.)で表示", value=st.session_state.editor_show_dots, key="toggle_dots")
    st.session_state.editor_show_dots = show_dots
    
    # データエディタ
    edited_df = st.data_editor(
        initial_df, 
        hide_index=True, 
        use_container_width=True,
        height=400,
        num_rows="dynamic",
        column_config={
            "Include": st.column_config.CheckboxColumn("有効", width="small"),
            "ID": st.column_config.TextColumn("ID", width="medium"),
            "Sequence": st.column_config.TextColumn("配列", width="large"),
        }
    )
    
    # プレビュー
    st.markdown("###### プレビュー (先頭20配列)")
    st.markdown(
        generate_alignment_html_from_df(edited_df, max_seqs=20, show_dots=show_dots), 
        unsafe_allow_html=True
    )
    
    # アクションボタン
    col1, col2 = st.columns([1, 3])
    with col1:
        if st.button("変更を保存して閉じる", type="primary"):
            st.session_state.phylo_aligned_df = edited_df
            st.rerun()
            
    with col2:
        st.caption("※「保存」を押すとメイン画面に反映されます。")
