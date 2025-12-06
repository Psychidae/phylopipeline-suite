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

from plotly.subplots import make_subplots

def get_colorscale(scheme="Default"):
    # 0:-/N, 1:A, 2:T, 3:G, 4:C, 5:Other
    # 戻り値: (colorscale_list, bar_color_map)
    if scheme == "AliView (Classic)":
        # A=Red, C=Green, G=Yellow, T=Blue (Typical simple)
        # Note: Clustal uses A=Red, G=Orange, C=Blue, T=Green
        # Let's use a high contrast set
        c_map = {
            0: '#FFFFFF', # Gap
            1: '#FF4444', # A Red
            2: '#44FF44', # T Green 
            3: '#FFD700', # G Yellow/Gold
            4: '#4444FF', # C Blue
            5: '#888888'  # Other
        }
    elif scheme == "Clustal":
        c_map = {
            0: '#FFFFFF', 1: '#80a0f0', 2: '#f01505', 3: '#f08040', 4: '#00ff00', 5: '#808080'
            # Note: Clustal is context dependent but we use static approximation
            # A: Red, G: Orange, C: Blue, T: Green
            # A=1(Red), T=2(Green), G=3(Orange), C=4(Blue) -> My mapping is A=1,T=2,G=3,C=4
        }
        c_map = {0: '#FFFFFF', 1: '#E41A1C', 2: '#4DAF4A', 3: '#FF7F00', 4: '#377EB8', 5: '#999999'}
    else: # Default (Pastel)
        c_map = {
            0: '#eeeeee', 1: '#ff9999', 2: '#99ff99', 3: '#ffff99', 4: '#9999ff', 5: '#cccccc'
        }
        
    # Plotly colorscale format
    colors = []
    # 6 discrete values: 0..5
    # Range 0.0-1.0 mapped to 0-5.
    # 0: 0.0 - 0.166
    # 1: 0.166 - 0.333
    # ...
    step = 1.0/6.0
    for i in range(6):
        c = c_map.get(i, '#ffffff')
        colors.append([i*step, c])
        colors.append([(i+1)*step, c])
        
    return colors, c_map

def plot_alignment_heatmap(matrix, ids, start_pos=None, end_pos=None, show_text=False, show_consensus_row=False, diff_mode=False, color_scheme="Default", dense_view=True):
    """Plotlyでヒートマップとコンセンサスバーを描画"""
    if matrix is None: return None
    
    # 逆マッピング
    int_to_char = {0: '-', 1: 'A', 2: 'T', 3: 'G', 4: 'C', 5: 'N'}
    
    # コンセンサス計算
    seq_len = matrix.shape[1]
    consensus_indices = []
    consensus_scores = []
    
    if matrix.size > 0:
        for col in range(seq_len):
            col_data = matrix[:, col]
            counts = np.bincount(col_data, minlength=6)
            total = len(col_data)
            mode_val = counts.argmax()
            score = counts[mode_val] / total if total > 0 else 0
            
            consensus_indices.append(mode_val)
            consensus_scores.append(score)
            
    consensus_row = np.array(consensus_indices, dtype=np.int8)

    # プロット用データ
    plot_matrix = matrix.copy()
    plot_ids = list(ids)

    # Heatmap上のコンセンサス行（オプション）
    if show_consensus_row:
        plot_matrix = np.vstack([consensus_row, plot_matrix])
        plot_ids = ["Consensus"] + plot_ids

    # テキスト行列作成
    text_matrix = None
    if show_text or diff_mode:
        text_matrix = np.empty(plot_matrix.shape, dtype=object)
        for val, char in int_to_char.items():
            text_matrix[plot_matrix == val] = char
            
        if diff_mode:
            for i in range(plot_matrix.shape[0]):
                if show_consensus_row and i == 0: continue
                # コンセンサス配列と比較
                match_mask = (plot_matrix[i] == consensus_row)
                text_matrix[i][match_mask] = '.' # AliView like dot

    # --- Plotly Subplots ---
    fig = make_subplots(
        rows=2, cols=1,
        shared_xaxes=True,
        vertical_spacing=0.01,
        row_heights=[0.15, 0.85]
    )

    # Colors
    colorscale, color_map = get_colorscale(color_scheme)

    # 1. Consensus Bar Chart (Top)
    x_axis = list(range(1, seq_len + 1))
    
    bar_colors = [color_map.get(idx, '#cccccc') for idx in consensus_indices]

    fig.add_trace(go.Bar(
        x=x_axis,
        y=consensus_scores,
        marker_color=bar_colors,
        showlegend=False,
        name="Consensus %",
        hoverinfo="y+text",
        text=[int_to_char[i] for i in consensus_indices]
    ), row=1, col=1)

    # 2. Heatmap (Bottom)
    heatmap_kwargs = dict(
        z=plot_matrix,
        x=x_axis,
        y=plot_ids,
        colorscale=colorscale,
        showscale=False,
        ygap=0 if dense_view else 1,
        xgap=0 if dense_view else 1,
    )
    
    if show_text:
        heatmap_kwargs["text"] = text_matrix
        heatmap_kwargs["texttemplate"] = "%{text}"
        # フォント設定 for readability
        # plotlyは自動調整するが、明示的に大きくしたい場合はlayout.font？
        # ここではtextfontを指定
        heatmap_kwargs["textfont"] = dict(family="mplus, monospace", color="black") 
    else:
        if diff_mode: heatmap_kwargs["hovertext"] = text_matrix

    fig.add_trace(go.Heatmap(**heatmap_kwargs), row=2, col=1)

    # トリム範囲 (Heatmap上のみ)
    if start_pos is not None and end_pos is not None:
        y_len = len(plot_ids)
        if start_pos > 0:
            fig.add_shape(type="rect", x0=0.5, y0=-0.5, x1=start_pos+0.5, y1=y_len-0.5,
                          fillcolor="rgba(50,50,50,0.6)", line=dict(width=0),
                          row=2, col=1)
        if end_pos < seq_len:
            fig.add_shape(type="rect", x0=end_pos+0.5, y0=-0.5, x1=seq_len+0.5, y1=y_len-0.5,
                          fillcolor="rgba(50,50,50,0.6)", line=dict(width=0),
                          row=2, col=1)

    # Layout
    cell_height = 20 # AliViewっぽく少し高さを確保
    total_height = 400 + (len(plot_ids) * cell_height * 0.5) 
    # PlotlyのHeatmapは高さ固定だとセルが伸び縮みする。
    # 正方形に近づけるには、heightをデータ数に合わせて動的に変える必要がある。
    # しかしStreamlit枠内での限界もある。
    
    fig.update_layout(
        title="Alignment Editor",
        height=max(500, len(plot_ids) * 15 + 150), # 最低500, 配列数に応じて縦に伸ばす
        margin=dict(l=20, r=20, t=40, b=20),
        dragmode='pan',
        yaxis2=dict(autorange="reversed")
    )
    
    fig.update_xaxes(range=[0.5, seq_len+0.5])
    
    return fig

@st.dialog("Advanced Alignment Editor", width="large")
def open_alignment_editor(initial_df, target_key="phylo_aligned_df"):
    """
    高機能アラインメントエディタ
    target_key: 保存時に更新するst.session_stateのキー
    """
    st.caption("全体を俯瞰して、不要な配列の削除や、両端のトリミングを行えます。")
    
    if 'editor_df' not in st.session_state:
        st.session_state.editor_df = initial_df.copy()
        
    df = st.session_state.editor_df
    
    # --- 1. ヒートマップ可視化 ---
    st.subheader("1. Visualization")
    
    # UI Control Row
    r1c1, r1c2, r1c3, r1c4 = st.columns(4)
    show_consensus = r1c1.checkbox("Consensus Row", value=True)
    show_text = r1c2.checkbox("Show Text", value=True)
    diff_mode = r1c3.checkbox("Diff Mode", value=False)
    dense_view = r1c4.checkbox("Dense View", value=True, help="Remove gaps between cells")
    
    r2c1, r2c2 = st.columns([1, 1])
    color_scheme = r2c1.selectbox("Color Scheme", ["Default", "AliView (Classic)", "Clustal"])
    
    plot_placeholder = st.empty()
    
    # --- 2. フィルタリング & 編集 ---
    st.subheader("2. Filtering")
    c1, c2 = st.columns([2, 1])
    
    with c1:
        edited_df = st.data_editor(
            df, 
            hide_index=True, 
            use_container_width=True,
            height=200,
            column_config={
                "Include": st.column_config.CheckboxColumn("Keep", width="small"),
                "ID": st.column_config.TextColumn("ID"), 
                "Sequence": st.column_config.TextColumn("Sequence", disabled=True),
            }
        )
    
    with c2:
        max_gap = st.slider("Max Gap Ratio", 0.0, 1.0, 1.0, 0.05)
        if st.button("Apply Gap Filter"):
            def calc_gap(s): return s.count("-") / len(s) if len(s) > 0 else 1.0
            edited_df["Include"] = edited_df["Sequence"].apply(lambda s: calc_gap(s) <= max_gap)
            st.session_state.editor_df = edited_df
            st.rerun()

    # --- 遅延描画 ---
    matrix_df = edited_df[edited_df["Include"] == True] if not edited_df.empty else edited_df
    matrix, ids = alignment_to_matrix(matrix_df)
    max_len = matrix.shape[1] if matrix is not None else 0
    
    with plot_placeholder.container():
        trim_range = st.slider("Trim Range (Keep Region)", 1, max_len, (1, max_len), key="trim_slider")
        start_pos, end_pos = trim_range
        
        if matrix is not None:
            fig = plot_alignment_heatmap(matrix, ids, start_pos, end_pos, show_text, show_consensus, diff_mode, color_scheme, dense_view)
            st.plotly_chart(fig, use_container_width=True, config={'scrollZoom': True, 'displayModeBar': True})

    # --- 3. 確定ボタン ---
    st.divider()
    col_l, col_r = st.columns([1, 1])
    with col_l:
        if st.button("Save & Close", type="primary"):
            final_df = edited_df[edited_df["Include"] == True].copy()
            s_idx = start_pos - 1
            e_idx = end_pos
            final_df["Sequence"] = final_df["Sequence"].apply(lambda s: s[s_idx:e_idx])
            
            st.session_state[target_key] = final_df
            del st.session_state.editor_df
            st.rerun()
            
    with col_r:
        if st.button("Cancel / Reset"):
            st.session_state.editor_df = initial_df.copy()
            st.rerun()
