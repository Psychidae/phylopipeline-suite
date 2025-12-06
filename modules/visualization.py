import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
# 絶対インポートを使用
from modules.bio_logic import find_orfs, calculate_translation_map

def create_main_figure(visible_results, alignment_ref, consensus_data, mismatches, target_pos, trace_win, 
                       show_orf_map, show_translation, table_id, quality_thresh):
    
    ref_len = len(consensus_data)
    ref_track = alignment_ref["display"]
    ref_peaks = np.array(ref_track["peak_locations"])
    avg_spacing = np.mean(np.diff(ref_peaks)) if len(ref_peaks)>1 else 10.0

    def get_ref_trace_x(base_idx):
        if 0 <= base_idx < len(ref_peaks): return ref_peaks[base_idx]
        elif base_idx < 0: return ref_peaks[0] + base_idx * avg_spacing
        else: return ref_peaks[-1] + (base_idx - len(ref_peaks) + 1) * avg_spacing

    center_x = get_ref_trace_x(target_pos)
    view_min = center_x - trace_win; view_max = center_x + trace_win
    
    has_orf = show_orf_map
    consensus_row = 1 if not has_orf else 2
    waveform_start_row = consensus_row + 1
    num_visible = len(visible_results)
    total_rows = num_visible + (2 if has_orf else 1)
    if num_visible == 0: total_rows += 1

    row_heights = []
    if has_orf: row_heights.append(0.15)
    row_heights.append(0.20 if show_translation else 0.10)
    remain = 1.0 - sum(row_heights)
    if num_visible > 0: row_heights.extend([remain/num_visible]*num_visible)
    else: row_heights.append(remain)

    fig = make_subplots(rows=total_rows, cols=1, row_heights=row_heights, vertical_spacing=0.02, shared_xaxes=True)
    all_x = [get_ref_trace_x(i) for i in range(ref_len)]
    x_inf_min = get_ref_trace_x(-100); x_inf_max = get_ref_trace_x(ref_len+100)

    # 1. ORF
    if has_orf:
        seq_str = "".join([d["base"] for d in consensus_data])
        orfs = find_orfs(seq_str, table_id=table_id)
        ymap = {1:3, 2:2, 3:1, -1:-1, -2:-2, -3:-3}
        corf = {1:'blue', 2:'blue', 3:'blue', -1:'red', -2:'red', -3:'red'}
        for o in orfs:
            fig.add_shape(type="rect", x0=get_ref_trace_x(o['start']), x1=get_ref_trace_x(o['end']), y0=ymap[o['frame']]-0.4, y1=ymap[o['frame']]+0.4, fillcolor=corf[o['frame']], line_width=0, opacity=0.5, row=1, col=1)
        for y in [0.5,1.5,2.5,-0.5,-1.5,-2.5]: fig.add_hline(y=y, line_color="lightgray", row=1, col=1)
        fig.update_yaxes(tickvals=[3,2,1,-1,-2,-3], ticktext=['+1','+2','+3','-1','-2','-3'], range=[-3.5,3.5], fixedrange=True, row=1, col=1)
        fig.add_vline(x=center_x, line_color="orange", width=2, row=1, col=1)

    # 2. Consensus
    fig.add_shape(type="rect", x0=x_inf_min, x1=x_inf_max, y0=0.6, y1=1.2, fillcolor="rgba(240,240,240,0.5)", line_width=0, layer="below", row=consensus_row, col=1)
    if show_translation:
        fig.add_shape(type="rect", x0=x_inf_min, x1=x_inf_max, y0=-0.1, y1=0.6, fillcolor="rgba(255,250,240,0.5)", line_width=0, layer="below", row=consensus_row, col=1)
        fig.add_shape(type="line", x0=x_inf_min, x1=x_inf_max, y0=0.6, y1=0.6, line=dict(color="gray", width=1, dash="dot"), row=consensus_row, col=1)

    cons_bases = [d["base"] for d in consensus_data]
    cons_cols = ["red" if i in mismatches else ("magenta" if d["base"].islower() else "black") for i,d in enumerate(consensus_data)]
    
    fig.add_trace(go.Scatter(x=all_x, y=[0.9]*ref_len, mode="text", text=cons_bases, textfont=dict(size=14, color=cons_cols, family="monospace", weight="bold"), hoverinfo="text", name="DNA"), row=consensus_row, col=1)
    
    if show_translation:
        for f in [1,2,3]:
            amap = calculate_translation_map(consensus_data, table_id, f)
            if amap:
                ax = [get_ref_trace_x(a["seq_idx"]) for a in amap]
                ac = [a["aa"] for a in amap]
                acol = ["red" if a["is_stop"] else "black" for a in amap]
                fig.add_trace(go.Scatter(x=ax, y=[0.5-(f-1)*0.2]*len(ax), mode="text", text=ac, textfont=dict(size=11, color=acol, family="monospace"), hoverinfo="text", showlegend=False, name=f"F{f}"), row=consensus_row, col=1)

    fig.add_trace(go.Scatter(x=[center_x], y=[1.1], mode="text", text=["▼"], textfont=dict(color="red", size=14), showlegend=False, hoverinfo="skip"), row=consensus_row, col=1)
    fig.add_vline(x=center_x, line_color="orange", width=2, row=consensus_row, col=1)
    fig.update_yaxes(range=[-0.1, 1.3] if show_translation else [0.5, 1.3], visible=False, fixedrange=True, row=consensus_row, col=1)

    # 3. Waveforms
    colors = {'A':'#2ca02c', 'T':'#d62728', 'C':'#1f77b4', 'G':'#000000'}
    if num_visible == 0:
        fig.add_annotation(text="No waveforms selected", showarrow=False, row=waveform_start_row, col=1)
    else:
        bases_in_view = (view_max - view_min) / avg_spacing
        show_wf_text = bases_in_view < 300

        for i, track in enumerate(visible_results):
            row = waveform_start_row + i
            d = track["display"]
            offset = track["offset"]
            
            # Warping logic
            q_start_ref = offset; q_end_ref = offset + len(d["sequence"])
            map_s = max(0, q_start_ref); map_e = min(len(ref_peaks), q_end_ref)
            
            if map_e <= map_s:
                fig.add_annotation(x=center_x, y=0.5, text="No overlap", showarrow=False, row=row, col=1)
                continue

            ref_idxs = list(range(map_s, map_e))
            q_peaks = [d["peak_locations"][k-offset] for k in ref_idxs]
            fp = [ref_peaks[k] for k in ref_idxs]
            
            q_tr_s = max(0, min(q_peaks)-100); q_tr_e = min(len(d["trace"]['A']), max(q_peaks)+100)
            q_x_orig = np.arange(q_tr_s, q_tr_e)
            q_x_warp = np.interp(q_x_orig, q_peaks, fp)
            
            step = max(1, len(q_x_warp)//8000)
            x_pl = q_x_warp[::step]
            
            ymax = 100
            for b in ['A','T','G','C']:
                y_raw = np.array(d["trace"][b])[q_tr_s:q_tr_e]
                if len(y_raw)>0: ymax = max(ymax, np.max(y_raw))
                fig.add_trace(go.Scatter(x=x_pl, y=y_raw[::step], mode='lines', line=dict(color=colors[b], width=1), showlegend=False, hoverinfo='y+name', name=b), row=row, col=1)

            if show_wf_text:
                txt_y = ymax * 1.15
                buf = 50
                k_s_v = int((view_min-ref_peaks[0])/avg_spacing) if len(ref_peaks)>0 else 0
                k_e_v = int((view_max-ref_peaks[0])/avg_spacing) if len(ref_peaks)>0 else 0
                k_start = max(map_s, k_s_v-buf); k_end = min(map_e, k_e_v+buf)
                
                tx, tc, tcol = [], [], []
                for k in range(k_start, k_end):
                    q_idx = k - offset
                    if 0 <= q_idx < len(d["sequence"]):
                        tx.append(ref_peaks[k])
                        tc.append(d["sequence"][q_idx])
                        qual = d["qualities"][q_idx]
                        if k == target_pos: tcol.append("red")
                        elif qual < quality_thresh: tcol.append("gray")
                        else: tcol.append(colors.get(tc[-1].upper(), "black"))
                
                if tx:
                    fig.add_trace(go.Scatter(x=tx, y=[txt_y]*len(tx), mode="text", text=tc, textfont=dict(size=12, color=tcol, weight="bold"), hoverinfo="text", showlegend=False), row=row, col=1)

            fig.add_vline(x=center_x, line_color="orange", dash="dot", row=row, col=1)
            fig.update_yaxes(title_text=track["name"][:8], range=[-ymax*0.1, ymax*1.3], showticklabels=False, fixedrange=True, row=row, col=1)

    fig.update_xaxes(range=[view_min, view_max], row=total_rows, col=1)
    for r in range(1, total_rows+1): fig.update_xaxes(showticklabels=False, row=r, col=1)
    fig.update_layout(height=250 + 180 * max(1, num_visible), margin=dict(t=10, b=10, l=10, r=10), dragmode="pan", hovermode="closest", clickmode="event+select", plot_bgcolor="white")
    return fig
