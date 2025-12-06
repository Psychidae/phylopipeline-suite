import streamlit as st
import pandas as pd
# モジュール内参照 (同じmodulesフォルダ内なので . で参照)
from .constants import IUPAC_CODES, GENETIC_CODES
from .bio_logic import (
    parse_ab1, trim_sequence, align_pair, get_rev_comp, 
    calculate_initial_consensus, process_contig_group
)
from .file_processing import create_grouping_dataframe
from .visualization import create_main_figure

def app_waveform_main():
    st.header("📈 Waveform Validator Pro")

    # --- サイドバー (Waveform専用設定) ---
    with st.sidebar:
        st.markdown("---")
        st.markdown("**Waveform Settings**")
        quality_thresh = st.slider("品質スコア閾値", 0, 60, 20)
        auto_trim = st.checkbox("自動トリミング", value=True)
        show_orf_map = st.checkbox("ORFマップを表示", value=False)
        show_translation = st.checkbox("アミノ酸翻訳を表示", value=True)
        code_name = st.selectbox("遺伝暗号表", list(GENETIC_CODES.keys()), index=0)
        table_id = GENETIC_CODES[code_name]

    # --- セッション状態 ---
    if "all_contigs" not in st.session_state: st.session_state.all_contigs = {}
    if "file_groups_df" not in st.session_state: st.session_state.file_groups_df = None
    if "current_sample" not in st.session_state: st.session_state.current_sample = None
    if "wf_pos" not in st.session_state: st.session_state.wf_pos = 0
    if "zoom_level" not in st.session_state: st.session_state.zoom_level = 300 
    if "last_uploaded_files" not in st.session_state: st.session_state.last_uploaded_files = []
    if "selected_tracks" not in st.session_state: st.session_state.selected_tracks = []

    # ファイルアップロード (メイン画面上部に配置)
    uploaded_files = st.file_uploader("AB1ファイル (複数可)", type=["ab1"], accept_multiple_files=True)

    # --- 1. ファイル読み込みとグルーピング ---
    if uploaded_files:
        if uploaded_files != st.session_state.last_uploaded_files:
            st.session_state.last_uploaded_files = uploaded_files
            st.session_state.all_contigs = {}
            df = create_grouping_dataframe(uploaded_files)
            df.insert(0, "Include", True) 
            st.session_state.file_groups_df = df

        if st.session_state.file_groups_df is not None:
            expander_open = not bool(st.session_state.all_contigs)
            
            with st.expander("📂 サンプルグループ設定 (編集可能)", expanded=expander_open):
                edited_df = st.data_editor(
                    st.session_state.file_groups_df,
                    column_config={
                        "Include": st.column_config.CheckboxColumn("有効", width="small"),
                        "Filename": st.column_config.TextColumn("ファイル名", disabled=True, width="medium"),
                        "Group": st.column_config.TextColumn("グループ名 (Sample ID)", width="medium")
                    },
                    hide_index=True,
                    use_container_width=True
                )
                
                valid_rows = edited_df[edited_df["Include"] == True]
                if not valid_rows.empty:
                    group_summary = valid_rows.groupby("Group")["Filename"].count().reset_index()
                    group_summary.columns = ["Group Name", "Count"]
                    st.dataframe(group_summary, hide_index=True, use_container_width=True)
                
                if st.button("🚀 解析開始 (Contig作成)", type="primary"):
                    st.session_state.file_groups_df = edited_df
                    groups = valid_rows.groupby("Group")
                    results_store = {}
                    prog_bar = st.progress(0)
                    total_groups = len(groups)
                    
                    for idx, (group_name, group_data) in enumerate(groups):
                        target_filenames = set(group_data["Filename"])
                        target_files = [f for f in uploaded_files if f.name in target_filenames]
                        if target_files:
                            res = process_contig_group(target_files, auto_trim, quality_thresh)
                            if res: results_store[group_name] = res
                        prog_bar.progress((idx + 1) / total_groups)
                    
                    prog_bar.empty()
                    st.session_state.all_contigs = results_store
                    if results_store:
                        st.session_state.current_sample = list(results_store.keys())[0]
                        st.session_state.wf_pos = 0
                        first_res = results_store[st.session_state.current_sample]["results"]
                        st.session_state.selected_tracks = [r["name"] for r in first_res]
                    st.rerun()

    # --- 2. 結果ビューワー ---
    if st.session_state.all_contigs:
        st.divider()
        sample_names = list(st.session_state.all_contigs.keys())
        
        # サンプルが削除された場合などの安全性チェック
        if st.session_state.current_sample not in sample_names:
            st.session_state.current_sample = sample_names[0] if sample_names else None
            
        if st.session_state.current_sample:
            selected_sample = st.selectbox("表示サンプル", sample_names, 
                                           index=sample_names.index(st.session_state.current_sample))
            
            if selected_sample != st.session_state.current_sample:
                st.session_state.current_sample = selected_sample
                st.session_state.wf_pos = 0
                new_res = st.session_state.all_contigs[selected_sample]["results"]
                st.session_state.selected_tracks = [r["name"] for r in new_res]
                st.rerun()

            contig_data = st.session_state.all_contigs[selected_sample]
            results = contig_data["results"]
            # 最新のコンセンサスを使用（トラック変更反映のため）
            if "consensus" not in contig_data: # 安全策
                 contig_data["consensus"] = calculate_initial_consensus(results)

            # --- トラック選択 ---
            all_track_names = [r["name"] for r in results]
            with st.expander("👁 表示波形選択", expanded=False):
                selected_tracks = st.multiselect(
                    "表示ファイル", all_track_names, 
                    default=st.session_state.selected_tracks
                )
            
            if selected_tracks != st.session_state.selected_tracks:
                st.session_state.selected_tracks = selected_tracks
                # 表示対象のみでコンセンサス再計算
                vis_res = [r for r in results if r["name"] in selected_tracks]
                if vis_res:
                    ref_len_orig = len(results[0]["display"]["sequence"])
                    new_cons = calculate_initial_consensus(vis_res, ref_length=ref_len_orig)
                    st.session_state.all_contigs[selected_sample]["consensus"] = new_cons
                st.rerun()

            consensus_data = contig_data["consensus"]
            ref_len = len(consensus_data)
            visible_results = [r for r in results if r["name"] in st.session_state.selected_tracks]
            alignment_ref = results[0]
            
            ref_trace_len = len(alignment_ref["display"]["trace"]["A"])
            max_zoom_val = max(1000, ref_trace_len // 2 + 500)

            # 不一致検出
            mismatches = []
            for i in range(ref_len):
                cons_base = consensus_data[i]["base"]
                comp_base = cons_base.upper()
                has_diff = False
                for track in visible_results:
                    loc = i - track["offset"]
                    if 0 <= loc < len(track["display"]["sequence"]):
                        b = track["display"]["sequence"][loc].upper()
                        if b not in ['N', '-'] and b != comp_base:
                            has_diff = True; break
                if has_diff: mismatches.append(i)

            # --- 編集エリア ---
            st.markdown(f"### 🧬 {selected_sample}")
            col_edit1, col_edit2, col_edit3 = st.columns([1, 2, 1])
            target_pos = st.session_state.wf_pos
            target_pos = max(0, min(ref_len - 1, target_pos))
            cur_base = consensus_data[target_pos]["base"]
            
            def find_next_mismatch(curr, direc):
                if not mismatches: return curr
                if direc == "next":
                    c = [m for m in mismatches if m > curr]
                    return c[0] if c else mismatches[0]
                else:
                    c = [m for m in mismatches if m < curr]
                    return c[-1] if c else mismatches[-1]
            def on_mm_sel():
                if st.session_state._mm_sel is not None: st.session_state.wf_pos = st.session_state._mm_sel

            with col_edit1:
                st.caption(f"不一致: {len(mismatches)} 箇所")
                if mismatches:
                    st.selectbox("Jump", mismatches, format_func=lambda x: f"Pos {x+1}", key="_mm_sel", index=None, on_change=on_mm_sel, label_visibility="collapsed", placeholder="Select mismatch...")
                c_p, c_n = st.columns(2)
                if c_p.button("◀ Prev"): st.session_state.wf_pos = find_next_mismatch(target_pos, "prev"); st.rerun()
                if c_n.button("Next ▶"): st.session_state.wf_pos = find_next_mismatch(target_pos, "next"); st.rerun()

            with col_edit2:
                st.write(f"Pos: **{target_pos+1}** | Base: **{cur_base}**")
                if "auto_next" not in st.session_state: st.session_state.auto_next = True
                def on_key():
                    v = st.session_state._k_in.upper()
                    if v in IUPAC_CODES:
                        st.session_state.all_contigs[selected_sample]["consensus"][target_pos]["base"] = v.lower()
                        st.session_state._k_in = ""
                        if st.session_state.auto_next:
                            np = find_next_mismatch(target_pos, "next")
                            st.session_state.wf_pos = np if np != target_pos else min(ref_len-1, target_pos+1)
                c_i, c_c = st.columns([2,1])
                with c_i: st.text_input("Input", max_chars=1, key="_k_in", on_change=on_key, label_visibility="collapsed", placeholder="Key input...")
                with c_c: st.checkbox("Auto Jump", key="auto_next")

            with col_edit3:
                def do_edit(b):
                    st.session_state.all_contigs[selected_sample]["consensus"][target_pos]["base"] = b.lower()
                    if st.session_state.auto_next:
                        np = find_next_mismatch(target_pos, "next")
                        st.session_state.wf_pos = np if np != target_pos else min(ref_len-1, target_pos+1)
                    st.rerun()
                with st.popover("Palette"):
                    c = st.columns(5)
                    for i,b in enumerate(['A','T','G','C','N']): 
                        if c[i].button(b): do_edit(b)
                    c2 = st.columns(6)
                    for i,b in enumerate(['R','Y','K','M','S','W']): 
                        if c2[i].button(b): do_edit(b)

            # --- グラフ描画 ---
            fig = create_main_figure(
                visible_results, alignment_ref, 
                consensus_data, mismatches, target_pos, 
                st.session_state.zoom_level, show_orf_map, show_translation, 
                table_id, quality_thresh
            )
            
            event = st.plotly_chart(
                fig, use_container_width=True, 
                on_select="rerun", selection_mode="points", 
                config={'scrollZoom': True, 'displayModeBar': True}
            )
            
            if event and event["selection"]["points"]:
                cx = int(round(event["selection"]["points"][0]["x"]))
                ref_track = results[0]["display"]
                if ref_track["peak_locations"]:
                    dists = [abs(p - cx) for p in ref_track["peak_locations"]]
                    nearest = dists.index(min(dists))
                    if nearest != st.session_state.wf_pos:
                        st.session_state.wf_pos = nearest
                        st.rerun()

            # --- フッター ---
            st.divider()
            c_dl, c_nav = st.columns([1, 2])
            with c_dl:
                final_seq_str = "".join([d["base"] for d in consensus_data]).replace("-", "")
                st.download_button(f"💾 Save {selected_sample}", data=f">{selected_sample}\n{final_seq_str}", file_name=f"{selected_sample}.fasta", use_container_width=True, type="primary")
            with c_nav:
                st.session_state.wf_pos = st.slider("Pos", 0, ref_len-1, st.session_state.wf_pos, label_visibility="collapsed")
                st.session_state.zoom_level = st.slider("Zoom", 50, max_zoom_val, st.session_state.zoom_level, label_visibility="collapsed")
    else:
        st.info("👈 AB1ファイルをアップロードしてください")