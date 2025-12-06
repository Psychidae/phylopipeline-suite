import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from modules.common import generate_alignment_html_from_df

def alignment_to_matrix(df):
    """アラインメントDFを数値マトリックスに変換"""
    # マッピング: A=1, T=2, G=3, C=4, N=0, -=0, Others=0
    # 色分けのため: -/Nは薄い色(0)、塩基は濃い色(1以上)とする
    # より詳細に: -:0, A:1, T:2, G:3, C:4, N:5
    seqs = df["Sequence"].tolist()
    ids = df["ID"].tolist()
    if not seqs: return None, None
    
    # 最大長に合わせる（通常はアラインメント済みなので同じはずだが念のため）
    max_len = max(len(s) for s in seqs)
    
    # 文字列を数値に変換するための辞書
    base_map = {
        '-': 0, 'N': 0, '?': 0, 
        'A': 1, 'T': 2, 'G': 3, 'C': 4,
        'a': 1, 't': 2, 'g': 3, 'c': 4
    }
    
    matrix = np.zeros((len(seqs), max_len), dtype=np.int8)
    
    for i, seq in enumerate(seqs):
        # numpy配列化して置換するほうが高速だが、可読性重視でループ処理
        # (BioPython等を使えば早いが、依存減らすため手動)
        for j, char in enumerate(seq):
            if j < max_len:
                matrix[i, j] = base_map.get(char, 5) # 5=Other
                
    return matrix, ids

def plot_alignment_heatmap(matrix, ids, start_pos=None, end_pos=None, show_text=False, show_consensus=False, diff_mode=False):
    """Plotlyでヒートマップを描画"""
    if matrix is None: return None
    
    # 逆マッピング
    int_to_char = {0: '-', 1: 'A', 2: 'T', 3: 'G', 4: 'C', 5: 'N'}
    
    # コンセンサス計算 (単純な多数決)
    # axis=0 (列) ごとに最頻値を求める
    # 0(Gap)も含めて多数決をとるか、Gap無視するか？
    # ここではGapも含めて単純多数決とする（視覚的なコンセンサスラインなので）
    consensus_row = []
    if matrix.size > 0:
        for col in range(matrix.shape[1]):
            counts = np.bincount(matrix[:, col], minlength=6)
            mode_val = counts.argmax()
            consensus_row.append(mode_val)
    consensus_row = np.array(consensus_row, dtype=np.int8)

    plot_matrix = matrix.copy()
    plot_ids = list(ids)

    # コンセンサス行追加
    if show_consensus:
        plot_matrix = np.vstack([consensus_row, plot_matrix])
        plot_ids = ["Consensus"] + plot_ids

    # テキスト行列作成 (表示用)
    text_matrix = None
    if show_text or diff_mode:
        text_matrix = np.empty(plot_matrix.shape, dtype=object)
        
        # 1. 基本文字埋め
        for val, char in int_to_char.items():
            text_matrix[plot_matrix == val] = char
            
        # 2. Diff Mode: コンセンサスと一致する場所を '*' に置換
        # Comparison with consensus_row
        if diff_mode:
            # 行ごとに比較
            for i in range(plot_matrix.shape[0]):
                # コンセンサス行自体は置換しない
                if show_consensus and i == 0: continue
                
                # 一致箇所マスク
                match_mask = (plot_matrix[i] == consensus_row)
                # Gap(-) 同士の一致も '*' にするか？ -> 通常はGapは空白っぽく見せたいが、
                # ユーザー要望は「変異箇所を明瞭に」なので、一致は全部アスタリスク
                # ただしGapは薄い色なので '*' があっても目立たないかも?
                # ここでは一致＝'*'とする
                text_matrix[i][match_mask] = '.' # '*'だとごちゃつくのでドット、または要望通りアスタリスク
                # ユーザー要望: "同じところはアスタリスク"
                text_matrix[i][match_mask] = '*'

    # トリム範囲
    shapes = []
    if start_pos is not None and end_pos is not None:
        max_col = plot_matrix.shape[1]
        y_len = len(plot_ids)
        # 左側の除外エリア
        if start_pos > 0:
            shapes.append(dict(type="rect", x0=0.5, y0=-0.5, x1=start_pos+0.5, y1=y_len-0.5, 
                               fillcolor="rgba(100,100,100,0.5)", line=dict(width=0)))
        # 右側の除外エリア
        if end_pos < max_col:
            shapes.append(dict(type="rect", x0=end_pos+0.5, y0=-0.5, x1=max_col+0.5, y1=y_len-0.5, 
                               fillcolor="rgba(100,100,100,0.5)", line=dict(width=0)))

    # カラースケール
    colors = [
        [0.0, '#eeeeee'], [0.16, '#eeeeee'], # 0 (Gap)
        [0.16, '#ff9999'], [0.33, '#ff9999'], # 1 (A)
        [0.33, '#99ff99'], [0.50, '#99ff99'], # 2 (T)
        [0.50, '#ffff99'], [0.66, '#ffff99'], # 3 (G)
        [0.66, '#9999ff'], [0.83, '#9999ff'], # 4 (C)
        [0.83, '#cccccc'], [1.0, '#cccccc']   # 5 (Other)
    ]
    
    # Heatmap引数
    heatmap_kwargs = dict(
        z=plot_matrix,
        x=list(range(1, plot_matrix.shape[1] + 1)),
        y=plot_ids,
        colorscale=colors,
        showscale=False,
        ygap=1,
        xgap=0,
    )
    
    if show_text:
        heatmap_kwargs["text"] = text_matrix
        heatmap_kwargs["texttemplate"] = "%{text}"
        # フォントサイズ調整（文字が見えるように）
        # 自動調整されるが、小さすぎると見えない
    else:
        # Hover info用にテキストを渡しておく（Diffが表示されると便利）
        if diff_mode:
             heatmap_kwargs["hovertext"] = text_matrix

    fig = go.Figure(data=go.Heatmap(**heatmap_kwargs))
    
    fig.update_layout(
        title="Alignment Overview",
        xaxis_title="Position",
        yaxis_title="Sequence ID",
        height=400 + (len(plot_ids) * 15), # テキストが入るので少し高く
        margin=dict(l=20, r=20, t=40, b=20),
        shapes=shapes,
        yaxis=dict(autorange="reversed")
    )
    
    return fig

@st.dialog("Advanced Alignment Editor", width="large")
def open_alignment_editor(initial_df):
    """
    高機能アラインメントエディタ
    - 俯瞰ヒートマップズーム
    - トリムスライダー
    - フィルタリング
    """
    st.info("全体を俯瞰して、不要な配列の削除や、両端のトリミングを行えます。")
    
    # セッションステート初期化（一時編集用）
    if 'editor_df' not in st.session_state:
        st.session_state.editor_df = initial_df.copy()
        
    df = st.session_state.editor_df
    
    # --- 1. ヒートマップ可視化 (Placeholder) ---
    st.subheader("1. Visualization & Trimming")
    
    # 表示オプション
    c_opt1, c_opt2, c_opt3 = st.columns(3)
    show_consensus = c_opt1.checkbox("Show Consensus Row", value=True)
    show_text = c_opt2.checkbox("Show ATGC Text", value=True)
    diff_mode = c_opt3.checkbox("Highlight Diff (*)", value=False, help="Consensusと同じ塩基をアスタリスクで表示します")

    plot_placeholder = st.empty() # 後で描画することで、テーブルの変更を即座に反映
    
    # --- 2. フィルタリング & 編集 ---
    st.subheader("2. Filtering & Editing")
    c1, c2 = st.columns([2, 1])
    
    with c1:
        # データエディタ（削除・ID編集用）
        st.caption("IDの編集や、チェックを外して削除が可能です。")
        edited_df = st.data_editor(
            df, 
            hide_index=True, 
            use_container_width=True,
            height=300,
            column_config={
                "Include": st.column_config.CheckboxColumn("Keep", width="small"),
                "ID": st.column_config.TextColumn("ID", disabled=False), # 編集可能に変更
                "Sequence": st.column_config.TextColumn("Sequence", disabled=True),
            }
        )
        # 編集結果をStateに同期（次回の基準にするため）
        # ただし無限ループに注意。editor_dfを更新するとdata_editorのinputが変わる。
        # data_editorはinputが変わるとリセットされることがあるが、experimental_data_editorの挙動による。
        # ここでは単純に `edited_df` を使ってプロットを描くことで表示上の同期を図る。
    
    with c2:
        st.markdown("##### Filter Actions")
        # ギャップ率でフィルタ
        max_gap = st.slider("Max Gap Ratio", 0.0, 1.0, 1.0, 0.05, help="Remove sequences with more than X% gaps")
        if st.button("Apply Gap Filter"):
            def calc_gap(s): return s.count("-") / len(s) if len(s) > 0 else 1.0
            # edited_dfをベースにフィルタ
            edited_df["Include"] = edited_df["Sequence"].apply(lambda s: calc_gap(s) <= max_gap)
            st.session_state.editor_df = edited_df
            st.rerun()

    # --- 遅延描画: ヒートマップ ---
    # テーブルでの変更（ID変更やInclude変更）を反映したデータをプロットに渡す
    matrix_df = edited_df[edited_df["Include"] == True] if not edited_df.empty else edited_df
    matrix, ids = alignment_to_matrix(matrix_df)
    max_len = matrix.shape[1] if matrix is not None else 0
    
    with plot_placeholder.container():
        # トリム用スライダー（この場所で定義すると上部に表示される）
        trim_range = st.slider("Trim Range (Keep Region)", 1, max_len, (1, max_len), key="trim_slider")
        start_pos, end_pos = trim_range
        
        if matrix is not None:
            fig = plot_alignment_heatmap(matrix, ids, start_pos, end_pos, show_text, show_consensus, diff_mode)
            # ズーム有効化 (scrollZoom)
            st.plotly_chart(fig, use_container_width=True, config={'scrollZoom': True, 'displayModeBar': True})

    # --- 3. 確定ボタン ---
    st.divider()
    col_l, col_r = st.columns([1, 1])
    with col_l:
        if st.button("Save & Close", type="primary"):
            # 1. Include=Falseを削除
            final_df = edited_df[edited_df["Include"] == True].copy()
            # 2. トリム実行 (0-indexed adjustment)
            # スライダーは1-based, python sliceは0-based
            s_idx = start_pos - 1
            e_idx = end_pos
            final_df["Sequence"] = final_df["Sequence"].apply(lambda s: s[s_idx:e_idx])
            
            # メインステートに保存
            st.session_state.phylo_aligned_df = final_df
            del st.session_state.editor_df # 掃除
            st.rerun()
            
    with col_r:
        if st.button("Cancel / Reset"):
            st.session_state.editor_df = initial_df.copy()
            st.rerun()
