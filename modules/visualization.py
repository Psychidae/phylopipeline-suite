import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
from modules.bio_logic import find_orfs, calculate_translation_map

def create_main_figure(visible_results, alignment_ref, consensus_data, mismatches, target_pos, trace_win, 
                       show_orf_map, show_translation, table_id, quality_thresh):
    
    ref_len = len(consensus_data)
    ref_track = alignment_ref["display"]
    
    # --- 座標変換ヘルパー ---
    ref_peaks = np.array(ref_track["peak_locations"])
    if len(ref_peaks) > 1:
        avg_spacing = np.mean(np.diff(ref_peaks))
    else:
        avg_spacing = 10.0

    def get_ref_trace_x(base_idx):
        if 0 <= base_idx < len(ref_peaks):
            return ref_peaks[base_idx]
        elif base_idx < 0:
            return ref_peaks[0] + base_idx * avg_spacing
        else:
            return ref_peaks[-1] + (base_idx - len(ref_peaks) + 1) * avg_spacing

    # --- 表示範囲 ---
    center_x = get_ref_trace_x(target_pos)
    view_min = center_x - trace_win
    view_max = center_x + trace_win
    
    # --- グラフ構成 ---
    has_orf = show_orf_map
    consensus_row = 1 if not has_orf else 2
    waveform_start_row = consensus_row + 1
    
    num_visible = len(visible_results)
    total_rows = num_visible + (2 if has_orf else 1)
    if num_visible == 0: total_rows += 1

    row_heights = []
    if has_orf: row_heights.append(0.15)
    row_heights.append(0.20 if show_translation else 0.10)
    
    remain_h = 0.85 - (0.15 if has_orf else 0) - (0.20 if show_translation else 0.10)
    if num_visible > 0:
        row_heights.extend([remain_h / num_visible] * num_visible)
    else:
        row_heights.append(remain_h)
    
    fig = make_subplots(
        rows=total_rows, cols=1,
        row_heights=row_heights,
        vertical_spacing=0.02,
        shared_xaxes=True
    )

    all_x_coords = [get_ref_trace_x(i) for i in range(ref_len)]
    x_inf_min = get_ref_trace_x(-100)
    x_inf_max = get_ref_trace_x(ref_len + 100)

    # ==========================================
    # 1. ORF Map
    # ==========================================
    if has_orf:
        final_seq_str = "".join([d["base"] for d in consensus_data])
        orfs = find_orfs(final_seq_str, table_id=table_id)
        
        # --- Logic to find a "Clean" single frame ---
        # If one ORF covers a significant portion of the sequence (e.g. > 70%), show only that frame.
        total_len = len(final_seq_str.replace("-", ""))
        best_frame = None
        max_coverage = 0
        
        for orf in orfs:
            coverage = orf['len'] / total_len if total_len > 0 else 0
            if coverage > 0.70: # Threshold for "Clean"
                if coverage > max_coverage:
                    max_coverage = coverage
                    best_frame = orf['frame']
        
        if best_frame is not None:
            # Filter to show only the best frame
            orfs = [o for o in orfs if o['frame'] == best_frame]

        ymap = {1:3, 2:2, 3:1, -1:-1, -2:-2, -3:-3}
        corf = {1:'blue', 2:'blue', 3:'blue', -1:'red', -2:'red', -3:'red'}
        
        for orf in orfs:
            x0 = get_ref_trace_x(orf['start'])
            x1 = get_ref_trace_x(orf['end'])
            y = ymap[orf['frame']]
            fig.add_shape(type="rect", x0=x0, y0=y-0.4, x1=x1, y1=y+0.4, fillcolor=corf[orf['frame']], line_width=0, opacity=0.5, row=1, col=1)
        
        for y in [0.5, 1.5, 2.5, -0.5, -1.5, -2.5]:
            fig.add_hline(y=y, line_color="lightgray", line_width=1, row=1, col=1)
            
        fig.update_yaxes(tickvals=[3,2,1,-1,-2,-3], ticktext=['+1','+2','+3','-1','-2','-3'], range=[-3.5,3.5], fixedrange=True, row=1, col=1)
        fig.add_vline(x=center_x, line_color="orange", line_width=2, row=1, col=1)

    # ==========================================
    # 2. Consensus & Translation
    # ==========================================
    fig.add_shape(type="rect", x0=x_inf_min, x1=x_inf_max, y0=0.6, y1=1.2, fillcolor="rgba(240, 240, 240, 0.5)", line_width=0, layer="below", row=consensus_row, col=1)
    
    if show_translation:
        fig.add_shape(type="rect", x0=x_inf_min, x1=x_inf_max, y0=-0.1, y1=0.6, fillcolor="rgba(255, 250, 240, 0.5)", line_width=0, layer="below", row=consensus_row, col=1)
        fig.add_shape(type="line", x0=x_inf_min, x1=x_inf_max, y0=0.6, y1=0.6, line=dict(color="gray", width=1, dash="dot"), row=consensus_row, col=1)

    cons_bases = [d["base"] for d in consensus_data]
    cons_colors = []
    for i, d in enumerate(consensus_data):
        base = d["base"]
        is_mm = i in mismatches
        if is_mm: cons_colors.append("red") 
        elif base.islower(): cons_colors.append("magenta") 
        else: cons_colors.append("black")

    fig.add_trace(go.Scatter(
        x=all_x_coords, y=[0.9] * ref_len, mode="text", text=cons_bases,
        textfont=dict(size=14, color=cons_colors, family="monospace", weight="bold"),
        hoverinfo="text", name="Consensus DNA"
    ), row=consensus_row, col=1)

    if show_translation:
        for frame in [1, 2, 3]:
            aa_map = calculate_translation_map(consensus_data, table_id, frame)
            if aa_map:
                aa_x = [get_ref_trace_x(a["seq_idx"]) for a in aa_map]
                aa_chars = [a["aa"] for a in aa_map]
                aa_colors = ["red" if a["is_stop"] else "black" for a in aa_map]
                y_pos = 0.5 - ((frame-1) * 0.2)
                fig.add_trace(go.Scatter(
                    x=aa_x, y=[y_pos] * len(aa_x), mode="text", text=aa_chars,
                    textfont=dict(size=11, color=aa_colors, family="monospace"),
                    hoverinfo="text", showlegend=False, name=f"Frame +{frame}"
                ), row=consensus_row, col=1)

    fig.add_trace(go.Scatter(
        x=[center_x], y=[1.1], mode="text", text=["▼"],
        textfont=dict(color="red", size=14),
        showlegend=False, hoverinfo="skip"
    ), row=consensus_row, col=1)

    # 修正: line_width=2
    fig.add_vline(x=center_x, line_color="orange", line_width=2, row=consensus_row, col=1)
    
    y_min_cons = -0.1 if show_translation else 0.5
    fig.update_yaxes(range=[y_min_cons, 1.3], visible=False, fixedrange=True, row=consensus_row, col=1)

    # ==========================================
    # 3. Waveforms
    # ==========================================
    colors = {'A': '#2ca02c', 'T': '#d62728', 'C': '#1f77b4', 'G': '#000000'}
    
    if num_visible == 0:
        fig.add_annotation(text="No waveforms selected", showarrow=False, row=waveform_start_row, col=1)
    else:
        bases_in_view = (view_max - view_min) / avg_spacing
        show_waveform_text = bases_in_view < 300

        for i, track in enumerate(visible_results):
            row_idx = waveform_start_row + i
            d = track["display"]
            offset = track["offset"]
            
            q_start_ref = offset
            q_end_ref = offset + len(d["sequence"])
            map_start = max(0, q_start_ref)
            map_end = min(len(ref_peaks), q_end_ref)
            
            if map_end <= map_start:
                 fig.add_annotation(x=center_x, y=0.5, text="No alignment overlap", showarrow=False, row=row_idx, col=1)
                 continue

            ref_indices = list(range(map_start, map_end))
            query_peaks_mapped = [d["peak_locations"][k - offset] for k in ref_indices]
            fp = [ref_peaks[k] for k in ref_indices]
            xp = query_peaks_mapped
            
            q_trace_start = max(0, min(xp) - 100)
            q_trace_end = min(len(d["trace"]['A']), max(xp) + 100)
            
            query_x_original = np.arange(q_trace_start, q_trace_end)
            query_x_warped = np.interp(query_x_original, xp, fp)
            
            step = max(1, len(query_x_warped) // 8000)
            x_plot = query_x_warped[::step]
            
            track_max_y = 100
            for b in ['A','T','G','C']:
                y_raw = np.array(d["trace"][b])[q_trace_start:q_trace_end]
                if len(y_raw) > 0: 
                    track_max_y = max(track_max_y, np.max(y_raw))
                y_plot = y_raw[::step]
                fig.add_trace(go.Scatter(
                    x=x_plot, y=y_plot, mode='lines', line=dict(color=colors[b], width=1),
                    showlegend=False, hoverinfo='y+name', name=b
                ), row=row_idx, col=1)

            if show_waveform_text:
                text_y_pos = track_max_y * 1.15
                buffer_idx = 50
                k_start_view = int((view_min - ref_peaks[0]) / avg_spacing) if len(ref_peaks) > 0 else 0
                k_end_view = int((view_max - ref_peaks[0]) / avg_spacing) if len(ref_peaks) > 0 else 0
                
                k_start = max(map_start, k_start_view - buffer_idx)
                k_end = min(map_end, k_end_view + buffer_idx)
                
                track_x, track_c, track_col = [], [], []
                for k in range(k_start, k_end):
                    q_idx = k - offset
                    if 0 <= q_idx < len(d["sequence"]):
                        px = ref_peaks[k]
                        char = d["sequence"][q_idx]
                        qual = d["qualities"][q_idx]
                        t_col = colors.get(char.upper(), "black")
                        if k == target_pos: t_col = "red"
                        elif qual < quality_thresh: t_col = "gray"
                        track_x.append(px); track_c.append(char); track_col.append(t_col)

                if track_x:
                    fig.add_trace(go.Scatter(
                        x=track_x, y=[text_y_pos]*len(track_x), mode="text", text=track_c,
                        textfont=dict(size=12, color=track_col, weight="bold"),
                        hoverinfo="text", showlegend=False
                    ), row=row_idx, col=1)

            # 修正: line_width=2
            fig.add_vline(x=center_x, line_color="orange", line_dash="dot", row=row_idx, col=1)
            
            fig.update_yaxes(title_text=track["name"][:8], range=[-track_max_y*0.1, track_max_y * 1.3], showticklabels=False, fixedrange=True, row=row_idx, col=1)

    fig.update_xaxes(range=[view_min, view_max], row=total_rows, col=1) 
    for r in range(1, total_rows+1):
        fig.update_xaxes(showticklabels=False, row=r, col=1)

    fig.update_layout(height=250 + 180 * max(1, num_visible), margin=dict(t=10, b=10, l=10, r=10), dragmode="pan", hovermode="closest", clickmode="event+select", plot_bgcolor="white")
    
    return fig
