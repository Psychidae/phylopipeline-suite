# modules/waveform_ui.py
import streamlit as st
import pandas as pd
from modules.constants import IUPAC_CODES, GENETIC_CODES
from modules.bio_logic import (
    parse_ab1, trim_sequence, align_pair, get_rev_comp, 
    calculate_initial_consensus, process_contig_group
)
from modules.file_processing import create_grouping_dataframe
from modules.visualization import create_main_figure

def app_waveform_main():
    st.header("📈 Waveform Validator Pro")

    with st.sidebar:
        st.markdown("---")
        st.markdown("### Waveform Settings")
        qt = st.slider("Quality", 0, 60, 20, key="wf_qt")
        at = st.checkbox("Auto Trim", value=True, key="wf_at")
        orf = st.checkbox("ORF Map", value=False, key="wf_orf")
        trans = st.checkbox("Translation", value=True, key="wf_trans")
        cname = st.selectbox("Genetic Code", list(GENETIC_CODES.keys()), key="wf_code")
        tid = GENETIC_CODES[cname]

    if "all_contigs" not in st.session_state: st.session_state.all_contigs = {}
    if "wf_groups" not in st.session_state: st.session_state.wf_groups = None
    if "wf_curr" not in st.session_state: st.session_state.wf_curr = None
    if "wf_pos" not in st.session_state: st.session_state.wf_pos = 0
    if "wf_zoom" not in st.session_state: st.session_state.wf_zoom = 300
    if "wf_last" not in st.session_state: st.session_state.wf_last = []
    if "wf_sel_tr" not in st.session_state: st.session_state.wf_sel_tr = []

    up = st.file_uploader("AB1 Files", type=["ab1"], accept_multiple_files=True, key="wf_up")

    if up:
        if up != st.session_state.wf_last:
            st.session_state.wf_last = up
            st.session_state.all_contigs = {}
            df = create_grouping_dataframe(up)
            df.insert(0, "Include", True)
            st.session_state.wf_groups = df

        if st.session_state.wf_groups is not None:
            with st.expander("📂 Groups", expanded=not bool(st.session_state.all_contigs)):
                edf = st.data_editor(st.session_state.wf_groups, hide_index=True, use_container_width=True)
                if st.button("🚀 Analyze", key="wf_go"):
                    st.session_state.wf_groups = edf
                    groups = edf[edf["Include"]==True].groupby("Group")
                    res = {}
                    prog = st.progress(0); tot = len(groups)
                    for i, (gn, gd) in enumerate(groups):
                        tfs = [f for f in up if f.name in set(gd["Filename"])]
                        if tfs:
                            r = process_contig_group(tfs, at, qt)
                            if r: res[gn] = r
                        prog.progress((i+1)/tot)
                    prog.empty()
                    st.session_state.all_contigs = res
                    if res:
                        st.session_state.wf_curr = list(res.keys())[0]
                        st.session_state.wf_pos = 0
                        st.session_state.wf_sel_tr = [t["name"] for t in res[st.session_state.wf_curr]["results"]]
                    st.rerun()

    if st.session_state.all_contigs:
        st.divider()
        snames = list(st.session_state.all_contigs.keys())
        if st.session_state.wf_curr not in snames: st.session_state.wf_curr = snames[0]
        
        sel = st.selectbox("Sample", snames, index=snames.index(st.session_state.wf_curr), key="wf_s_sel")
        if sel != st.session_state.wf_curr:
            st.session_state.wf_curr = sel
            st.session_state.wf_pos = 0
            st.session_state.wf_sel_tr = [t["name"] for t in st.session_state.all_contigs[sel]["results"]]
            st.rerun()

        data = st.session_state.all_contigs[sel]
        results = data["results"]
        
        with st.expander("👁 View Tracks", expanded=False):
            sel_tr = st.multiselect("Tracks", [r["name"] for r in results], default=st.session_state.wf_sel_tr, key="wf_t_sel")
        
        if sel_tr != st.session_state.wf_sel_tr:
            st.session_state.wf_sel_tr = sel_tr
            vis = [r for r in results if r["name"] in sel_tr]
            if vis:
                from modules.bio_logic import calculate_initial_consensus
                rl = len(results[0]["display"]["sequence"])
                st.session_state.all_contigs[sel]["consensus"] = calculate_initial_consensus(vis, ref_length=rl)
            st.rerun()

        cons_data = data["consensus"]
        rlen = len(cons_data)
        vis_res = [r for r in results if r["name"] in st.session_state.wf_sel_tr]
        align_ref = results[0]
        max_zoom = max(1000, len(align_ref["display"]["trace"]["A"])//2 + 500)

        mms = []
        for i in range(rlen):
            cb = cons_data[i]["base"].upper()
            diff = False
            for t in vis_res:
                loc = i - t["offset"]
                if 0 <= loc < len(t["display"]["sequence"]):
                    b = t["display"]["sequence"][loc].upper()
                    if b not in ['N','-'] and b != cb: diff = True; break
            if diff: mms.append(i)

        st.markdown(f"### {sel}")
        c1, c2, c3 = st.columns([1,2,1])
        tpos = max(0, min(rlen-1, st.session_state.wf_pos))
        cbase = cons_data[tpos]["base"]

        def f_next(c, d):
            if not mms: return c
            if d == "n":
                cand = [m for m in mms if m > c]
                return cand[0] if cand else mms[0]
            cand = [m for m in mms if m < c]
            return cand[-1] if cand else mms[-1]

        def on_jump(): 
            if st.session_state.wf_j is not None: st.session_state.wf_pos = st.session_state.wf_j

        with c1:
            if mms: st.selectbox("Jump", mms, format_func=lambda x:f"Pos {x+1}", key="wf_j", index=None, on_change=on_jump)
            else: st.success("No Mismatches")
            cp, cn = st.columns(2)
            if cp.button("◀"): st.session_state.wf_pos = f_next(tpos, "p"); st.rerun()
            if cn.button("▶"): st.session_state.wf_pos = f_next(tpos, "n"); st.rerun()

        with c2:
            st.write(f"Pos: **{tpos+1}** | **{cbase}**")
            if "wf_auto" not in st.session_state: st.session_state.wf_auto = True
            def on_k():
                v = st.session_state.wf_k.upper()
                if v in IUPAC_CODES:
                    st.session_state.all_contigs[sel]["consensus"][tpos]["base"] = v.lower()
                    st.session_state.wf_k = ""
                    if st.session_state.wf_auto:
                        np = f_next(tpos, "n")
                        st.session_state.wf_pos = np if np != tpos else min(rlen-1, tpos+1)
            c_in, c_ch = st.columns([2,1])
            with c_in: st.text_input("Key", max_chars=1, key="wf_k", on_change=on_k)
            with c_ch: st.checkbox("Auto", key="wf_auto")

        with c3:
            def do_ed(b):
                st.session_state.all_contigs[sel]["consensus"][tpos]["base"] = b.lower()
                if st.session_state.wf_auto:
                    np = f_next(tpos, "n")
                    st.session_state.wf_pos = np if np != tpos else min(rlen-1, tpos+1)
                st.rerun()
            with st.popover("Palette"):
                c = st.columns(5)
                for i,b in enumerate(['A','T','G','C','N']): 
                    if c[i].button(b, key=f"p1{i}"): do_ed(b)
                c = st.columns(6)
                for i,b in enumerate(['R','Y','K','M','S','W']): 
                    if c[i].button(b, key=f"p2{i}"): do_ed(b)

        fig = create_main_figure(vis_res, align_ref, cons_data, mms, tpos, st.session_state.wf_zoom, orf, trans, tid, qt)
        ev = st.plotly_chart(fig, use_container_width=True, on_select="rerun", selection_mode="points", config={'scrollZoom':True})
        if ev and ev["selection"]["points"]:
            cx = int(round(ev["selection"]["points"][0]["x"]))
            rt = align_ref["display"]
            if rt["peak_locations"]:
                dst = [abs(p-cx) for p in rt["peak_locations"]]
                nr = dst.index(min(dst))
                if nr != st.session_state.wf_pos: st.session_state.wf_pos = nr; st.rerun()

        # --- Genetic Analysis & Stop Codon Report ---
        if trans:
            st.divider()
            st.markdown(f"### 🧬 Genetic Analysis (Code: {cname})")
            from modules.bio_logic import calculate_translation_map
            
            # Check all frames
            found_issues = False
            
            # Dictionary to hold stop codons for dropdown: key=Label, value=Position
            stop_codon_options = {}
            
            cols = st.columns(3)
            for i, frame in enumerate([1, 2, 3]):
                with cols[i]:
                    aa_map = calculate_translation_map(cons_data, tid, frame)
                    stops = [a for a in aa_map if a["is_stop"]]
                    
                    if stops:
                        st.error(f"Frame +{frame}: {len(stops)} Stops")
                        found_issues = True
                        for s in stops:
                            # s['seq_idx'] is the middle base of the codon
                            mid_idx = s["seq_idx"]
                            codon_rng = range(mid_idx-1, mid_idx+2)
                            
                            # Check for mismatches in this codon
                            related_mms = [m for m in mms if m in codon_rng]
                            
                            # Determine status
                            status_icon = "🔴"
                            status_msg = "Stop"
                            if related_mms:
                                status_icon = "⚠️"
                                status_msg = f"Mismatch@{related_mms[0]+1}"
                            
                            label = f"Fr+{frame}: Pos {mid_idx+1} ({status_msg})"
                            stop_codon_options[label] = mid_idx
                    else:
                        st.success(f"Frame +{frame}: Clean")
            
            if found_issues:
                st.info("💡 Select a Stop Codon below to inspect the waveform.")
                
                # Callback for jump
                def on_stop_select():
                    selected_label = st.session_state.wf_stop_sel
                    if selected_label:
                        pos = stop_codon_options[selected_label]
                        st.session_state.wf_pos = pos
                        # We don't need explicit rerun here as on_change triggers it? 
                        # Actually Streamlit flow sometimes needs it or just let it rerun.
                
                # Dropdown
                st.selectbox(
                    "Jump to Stop Codon:",
                    options=list(stop_codon_options.keys()),
                    index=None,
                    key="wf_stop_sel",
                    on_change=on_stop_select,
                    placeholder="Select a stop codon..."
                )
                
                # Show Suggestions based on Current Cursor Position
                current_stop_target = None
                current_stop_frame = None
                
                # Scan all frames to see if current pos is inside a stop codon
                for frame in [1, 2, 3]:
                    aa_map = calculate_translation_map(cons_data, tid, frame)
                    stops = [a for a in aa_map if a["is_stop"]]
                    for s in stops:
                        mid_idx = s["seq_idx"]
                        # Check if tpos is within this codon (mid-1 to mid+1)
                        if abs(tpos - mid_idx) <= 1:
                            current_stop_target = mid_idx
                            current_stop_frame = frame
                            break
                    if current_stop_target is not None:
                        break
                
                if current_stop_target is not None:
                    with st.expander(f"💡 Fix Suggestions (Frame +{current_stop_frame})", expanded=True):
                        mid_idx = current_stop_target
                        
                        # Get exact codon from consensus
                        if 0 < mid_idx < len(cons_data)-1:
                            codon_bases = [
                                cons_data[mid_idx-1]["base"],
                                cons_data[mid_idx]["base"],
                                cons_data[mid_idx+1]["base"]
                            ]
                            codon_str = "".join(codon_bases).upper()
                            
                            st.write(f"selected codon: **{codon_str}** (at {mid_idx+1})")
                            
                            if "N" in codon_str or "-" in codon_str:
                                st.warning("Codon contains gaps or N.")
                            else:
                                from modules.bio_logic import suggest_stop_fixes
                                fixes = suggest_stop_fixes(codon_str, tid)
                                
                                if fixes:
                                    st.caption("Single-base changes that resolve this stop:")
                                    cols_fix = st.columns(4)
                                    for i, fix in enumerate(fixes):
                                        c = cols_fix[i % 4]
                                        abs_pos = mid_idx - 1 + fix['pos']
                                        is_mm = abs_pos in mms
                                        
                                        msg = f"{fix['from']}➔**{fix['to']}** ({fix['aa']})"
                                        if is_mm: msg += " ⚠️" 
                                        
                                        if c.button(msg, key=f"fix_{mid_idx}_{i}", help=f"Mutate Pos {abs_pos+1} to {fix['to']}"):
                                            st.session_state.all_contigs[sel]["consensus"][abs_pos]["base"] = fix['to'].lower()
                                            st.session_state.wf_pos = abs_pos
                                            st.session_state.wf_k = ""
                                            st.rerun()
                                else:
                                    st.info("No single-base mutation resolves this stop.")

        st.divider()
        cd, cn = st.columns([1, 2])
        with cd:
            seq = "".join([d["base"] for d in cons_data]).replace("-","")
            st.download_button("Save", f">{sel}\n{seq}", "cons.fasta", type="primary")
        with cn:
            st.session_state.wf_pos = st.slider("P", 0, rlen-1, st.session_state.wf_pos)
            st.session_state.wf_zoom = st.slider("Z", 50, max_zoom, st.session_state.wf_zoom)

    else: st.info("Upload .ab1 files")
