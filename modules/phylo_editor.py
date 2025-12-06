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
    """
    カラースキーム定義
    戻り値: (colorscale_list, bar_color_map, text_color)
    """
    # 0:-/N, 1:A, 2:T, 3:G, 4:C, 5:Other
    text_color = "black"
    
    if scheme == "AliView (Classic)":
        # High contrast: A=Red, C=Green, G=Yellow, T=Blue
        c_map = {
            0: '#FFFFFF', # Gap (White)
            1: '#FF4444', # A Red
            2: '#44FF44', # T Green 
            3: '#FFD700', # G Yellow
            4: '#4444FF', # C Blue
            5: '#AAAAAA'  # Other
        }
    elif scheme == "Clustal":
        # Clustal approximation
        c_map = {0: '#FFFFFF', 1: '#E41A1C', 2: '#4DAF4A', 3: '#FF7F00', 4: '#377EB8', 5: '#999999'}
        text_color = "black"
    elif scheme == "Nucleotide (Dark)":
        # Dark mode style
        c_map = {0: '#222222', 1: '#FF5555', 2: '#55FF55', 3: '#FFFF55', 4: '#5555FF', 5: '#666666'}
        text_color = "white"
    else: # Default (Pastel)
        c_map = {
            0: '#eeeeee', 1: '#ff9999', 2: '#99ff99', 3: '#ffff99', 4: '#9999ff', 5: '#cccccc'
        }
    
    # Plotly Colorscale construction
    colors = []
    step = 1.0/6.0
    for i in range(6):
        c = c_map.get(i, '#ffffff')
        colors.append([i*step, c])
        colors.append([(i+1)*step, c])
        
    return colors, c_map, text_color

def plot_alignment_heatmap(matrix, ids, start_pos=None, end_pos=None, show_text=False, show_consensus_row=False, diff_mode=False, color_scheme="Default", dense_view=True):
    """
    Plotlyでアラインメントを描画 (3部構成: Bar / Consensus / Alignment)
    """
    if matrix is None: return None
    
    # 1. データの準備
    int_to_char = {0: '-', 1: 'A', 2: 'T', 3: 'G', 4: 'C', 5: 'N'}
    seq_len = matrix.shape[1]
    n_seq = matrix.shape[0]

    # コンセンサス計算
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

    # 2. サブプロット構成
    # Row 1: Bar Chart (0.15)
    # Row 2: Consensus Seq (0.05 or fixed height)
    # Row 3: Main Alignment (0.80)
    fig = make_subplots(
        rows=3, cols=1,
        shared_xaxes=True,
        vertical_spacing=0.01, # 隙間を詰める
        row_heights=[0.15, 0.05, 0.80],
        subplot_titles=("", "", "")
    )

    colorscale, c_map, text_color = get_colorscale(color_scheme)
    x_axis = list(range(1, seq_len + 1))

    # --- Row 1: Consensus Strength (Bar) ---
    bar_colors = [c_map.get(idx, '#cccccc') for idx in consensus_indices]
    fig.add_trace(go.Bar(
        x=x_axis, y=consensus_scores,
        marker_color=bar_colors,
        showlegend=False,
        name="Consensus %",
        hoverinfo="y+text",
        text=[int_to_char[i] for i in consensus_indices]
    ), row=1, col=1)

    # --- Row 2: Consensus Sequence (Fixed Header) ---
    # コンセンサス配列のみのヒートマップ
    # 常にテキストを表示して見やすくする
    cons_text = [int_to_char[v] for v in consensus_indices]
    fig.add_trace(go.Heatmap(
        z=consensus_row.reshape(1, -1),
        x=x_axis,
        y=["Consensus"],
        colorscale=colorscale,
        showscale=False,
        text=[cons_text],
        texttemplate="%{text}",
        textfont=dict(family="monospace", color=text_color, size=14),
        xgap=0 if dense_view else 1,
        ygap=0,
        zsmooth=False, # Pixelated rendering
        hoverinfo="x+text"
    ), row=2, col=1)

    # --- Row 3: Main Alignment ---
    # コンセンサス行を含めない（別枠にしたので）
    # Diff Modeの計算用テキスト
    text_matrix = None
    if show_text or diff_mode:
        text_matrix = np.empty(matrix.shape, dtype=object)
        for val, char in int_to_char.items():
            text_matrix[matrix == val] = char
            
        if diff_mode:
            # コンセンサスと比較して一致を '.' に
            for i in range(n_seq):
                match_mask = (matrix[i] == consensus_row)
                text_matrix[i][match_mask] = '.'

    heatmap_kwargs = dict(
        z=matrix,
        x=x_axis,
        y=ids,
        colorscale=colorscale,
        showscale=False,
        xgap=0 if dense_view else 1,
        ygap=0 if dense_view else 1,
        zsmooth=False, # 重要: ズームアウト時の消失（エイリアシング）を防ぐ
    )

    # テキスト表示制御
    # 大規模データで全セルテキスト表示は重い＆見えないので、閾値を設けるかユーザー設定に従う
    # ここではshow_textに従うが、フォントサイズを調整
    if show_text:
        heatmap_kwargs["text"] = text_matrix
        heatmap_kwargs["texttemplate"] = "%{text}"
        heatmap_kwargs["textfont"] = dict(family="monospace", color=text_color)
    else:
        # Hoverには常に情報を出す
        if diff_mode:
            heatmap_kwargs["hovertext"] = text_matrix
        else:
             # 通常モードでもHoverで塩基が見えると嬉しい
             # しかしtext_matrixを作ると重いので、diff_mode or show_textの時だけに生成している
             pass

    fig.add_trace(go.Heatmap(**heatmap_kwargs), row=3, col=1)

    # --- Shapes (Trimming Mask) ---
    # Row 3 (Main) にのみ適用
    if start_pos is not None and end_pos is not None:
        if start_pos > 0:
            fig.add_shape(type="rect", x0=0.5, y0=-0.5, x1=start_pos+0.5, y1=len(ids)-0.5,
                          fillcolor="rgba(0,0,0,0.5)", line=dict(width=0), 
                          row=3, col=1)
        if end_pos < seq_len:
            fig.add_shape(type="rect", x0=end_pos+0.5, y0=-0.5, x1=seq_len+0.5, y1=len(ids)-0.5,
                          fillcolor="rgba(0,0,0,0.5)", line=dict(width=0),
                          row=3, col=1)

    # --- Layout ---
    # 高さ調整: Main Alignmentの行数依存
    # Row1(Bar) + Row2(Consensus) は固定高
    # Row2は1行分なので非常に小さい。
    min_height = 500
    dynamic_height = 150 + (len(ids) * 20) # Header + list
    
    fig.update_layout(
        title="Alignment Editor",
        height=max(min_height, dynamic_height),
        margin=dict(l=20, r=20, t=40, b=20),
        dragmode='pan',
        # Axes settings
        yaxis=dict(domain=[0.85, 1.0]),     # Row 1
        yaxis2=dict(domain=[0.80, 0.85], fixedrange=True), # Row 2 (Fixed Y)
        yaxis3=dict(domain=[0.0, 0.79], autorange="reversed"), # Row 3 (Scrollable)
    )
    
    # 軸連動
    fig.update_xaxes(range=[0.5, seq_len+0.5])
    
    # LayoutのGrid調整がmake_subplotsで上書きされる場合があるので、必要ならupdate_layoutで調整
    # make_subplotsのrow_heightsは相対比率なので、厳密なピクセル固定は難しい。
    # しかしdomain指定で上書き可能。
    
    # Row 2 (Consensus) のラベルを見やすく
    fig.update_yaxes(showticklabels=False, row=1, col=1) # Bar chart Y axis hidden
    fig.update_yaxes(showticklabels=True, row=2, col=1)  # "Consensus" label
    
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
    color_scheme = r2c1.selectbox("Color Scheme", ["Default", "AliView (Classic)", "Clustal", "Nucleotide (Dark)"])
    
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
