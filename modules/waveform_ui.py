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


# ---------------------------------------------------------------------------
# Keyboard-shortcut JS snippet (injected via st.html)
# ---------------------------------------------------------------------------
_KEYBOARD_JS = """
<script>
(function() {
    // Avoid re-registering on Streamlit re-render
    if (window.__wf_kb_registered) return;
    window.__wf_kb_registered = true;

    document.addEventListener('keydown', function(e) {
        // Skip when user is typing in an input / textarea
        const tag = (e.target.tagName || '').toLowerCase();
        if (tag === 'input' || tag === 'textarea' || e.target.isContentEditable) return;

        const keyMap = {
            'ArrowLeft':  'kb_prev',
            'ArrowRight': 'kb_next',
            '=':          'kb_zoomin',
            '-':          'kb_zoomout',
            '+':          'kb_ins',
            'Delete':     'kb_del',
            'Backspace':  'kb_del',
        };

        // Base keys (A/T/G/C/N etc.)
        const baseKeys = 'atgcnrykmsw';

        let targetId = keyMap[e.key];
        if (!targetId && baseKeys.includes(e.key.toLowerCase())) {
            targetId = 'kb_base_' + e.key.toLowerCase();
        }

        if (targetId) {
            // Find the hidden button inside the Streamlit DOM
            const btns = window.parent.document.querySelectorAll('button[kind="secondary"]');
            for (const btn of btns) {
                if (btn.id === targetId || btn.getAttribute('data-testid') === targetId) {
                    btn.click();
                    e.preventDefault();
                    return;
                }
            }
            // Fallback: search by aria-label
            const btn2 = window.parent.document.querySelector('button[aria-label="' + targetId + '"]');
            if (btn2) { btn2.click(); e.preventDefault(); }
        }
    });
})();
</script>
"""


def app_waveform_main():
    st.header("📈 Waveform Validator Pro")

    # ================================================================
    # Sidebar — Settings
    # ================================================================
    with st.sidebar:
        st.markdown("---")
        st.markdown("### ⚙️ 設定")
        qt = st.slider("Quality Threshold", 0, 60, 20, key="wf_qt",
                        help="塩基のクオリティ閾値。低品質領域のトリム基準に使用")
        at = st.checkbox("Auto Trim", value=True, key="wf_at",
                         help="低品質の末端を自動的にトリム")
        orf = st.checkbox("ORF Map", value=False, key="wf_orf",
                          help="Open Reading Frame マップを波形上部に表示")
        trans = st.checkbox("Translation", value=True, key="wf_trans",
                            help="3フレームのアミノ酸翻訳を表示")
        cname = st.selectbox("Genetic Code", list(GENETIC_CODES.keys()), key="wf_code",
                             help="翻訳に使用するコドンテーブル")
        tid = GENETIC_CODES[cname]

        st.markdown("---")
        st.markdown("### ⌨️ ショートカット")
        st.markdown("""
        | キー | 操作 |
        |------|------|
        | `←` `→` | ミスマッチ間を移動 |
        | `A`~`N` | 塩基を直接入力 |
        | `+` | 塩基を挿入 |
        | `Del` | 塩基を削除 |
        | `-` `=` | ズーム調整 |
        """, unsafe_allow_html=True)

    # ================================================================
    # Session State Initialization
    # ================================================================
    defaults = {
        "all_contigs": {}, "wf_groups": None, "wf_curr": None,
        "wf_pos": 0, "wf_zoom": 300, "wf_last": [], "wf_sel_tr": [],
        "wf_stop_opts": {}, "wf_auto": True,
    }
    for k, v in defaults.items():
        if k not in st.session_state:
            st.session_state[k] = v

    # ================================================================
    # 操作ガイド (collapsible)
    # ================================================================
    with st.expander("📖 操作ガイド", expanded=False):
        st.markdown("""
**基本的な使い方**
1. `.ab1` ファイルをアップロードし、グループを確認して **🚀 Analyze** をクリック
2. 波形チャート上でクリックするか、◀ ▶ ボタンでミスマッチ位置を順にチェック
3. 塩基の編集：キーボードで直接入力（A/T/G/C/N）するか、パレットボタンを使用
4. 塩基の挿入/削除：＋ / ✕ ボタン、またはキーボード `+` / `Delete`
5. 完了したら **💾 FASTA保存** でコンセンサス配列をダウンロード

**色の意味**
- 🟠 オレンジ = 現在のカーソル位置
- 🔴 赤 = トラック間のミスマッチ
- 🟣 紫/マゼンタ = 手動で編集された塩基
- 灰色 = 低品質 (Quality Threshold 以下)
        """)

    # ================================================================
    # File Upload
    # ================================================================
    up = st.file_uploader("AB1 ファイルをアップロード", type=["ab1"],
                          accept_multiple_files=True, key="wf_up",
                          help="Sanger シーケンスの .ab1 ファイルを1つ以上選択")

    if up:
        if up != st.session_state.wf_last:
            st.session_state.wf_last = up
            st.session_state.all_contigs = {}
            df = create_grouping_dataframe(up)
            df.insert(0, "Include", True)
            st.session_state.wf_groups = df

        if st.session_state.wf_groups is not None:
            with st.expander("📂 ファイルグループ", expanded=not bool(st.session_state.all_contigs)):
                st.caption("💡 ファイル名から自動グループ化されています。グループ名や含める/除外を編集できます。")
                edf = st.data_editor(st.session_state.wf_groups, hide_index=True, use_container_width=True)
                if st.button("🚀 Analyze", key="wf_go", type="primary"):
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

    # ================================================================
    # Main Editor View (only when data is loaded)
    # ================================================================
    if not st.session_state.all_contigs:
        st.info("📁 `.ab1` ファイルをアップロードして解析を開始してください。")
        return

    st.divider()
    snames = list(st.session_state.all_contigs.keys())
    if st.session_state.wf_curr not in snames:
        st.session_state.wf_curr = snames[0]

    # --- Sample & Track Selection ---
    sel_col, track_col = st.columns([2, 1])
    with sel_col:
        sel = st.selectbox("🧪 サンプル", snames,
                           index=snames.index(st.session_state.wf_curr), key="wf_s_sel",
                           help="解析するサンプルを選択")
    if sel != st.session_state.wf_curr:
        st.session_state.wf_curr = sel
        st.session_state.wf_pos = 0
        st.session_state.wf_sel_tr = [t["name"] for t in st.session_state.all_contigs[sel]["results"]]
        st.rerun()

    data = st.session_state.all_contigs[sel]
    results = data["results"]

    with track_col:
        with st.popover("👁 トラック表示"):
            sel_tr = st.multiselect("表示トラック", [r["name"] for r in results],
                                    default=st.session_state.wf_sel_tr, key="wf_t_sel",
                                    help="波形に表示するトラックを選択")

    if sel_tr != st.session_state.wf_sel_tr:
        st.session_state.wf_sel_tr = sel_tr
        vis = [r for r in results if r["name"] in sel_tr]
        if vis:
            rl = len(results[0]["display"]["sequence"])
            st.session_state.all_contigs[sel]["consensus"] = calculate_initial_consensus(vis, ref_length=rl)
        st.rerun()

    cons_data = data["consensus"]
    rlen = len(cons_data)
    vis_res = [r for r in results if r["name"] in st.session_state.wf_sel_tr]
    align_ref = results[0]
    max_zoom = max(1000, len(align_ref["display"]["trace"]["A"])//2 + 500)

    # Calculate mismatches
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

    tpos = max(0, min(rlen-1, st.session_state.wf_pos))
    cbase = cons_data[tpos]["base"]

    # --- Helper functions ---
    def f_next(c, d):
        if not mms: return c
        if d == "n":
            cand = [m for m in mms if m > c]
            return cand[0] if cand else mms[0]
        cand = [m for m in mms if m < c]
        return cand[-1] if cand else mms[-1]

    def do_edit(b):
        st.session_state.all_contigs[sel]["consensus"][tpos]["base"] = b.lower()
        if st.session_state.wf_auto:
            np_pos = f_next(tpos, "n")
            st.session_state.wf_pos = np_pos if np_pos != tpos else min(rlen-1, tpos+1)
        st.rerun()

    def do_prev():
        st.session_state.wf_pos = f_next(tpos, "p")

    def do_next():
        st.session_state.wf_pos = f_next(tpos, "n")

    def do_insert():
        st.session_state.all_contigs[sel]["consensus"].insert(tpos, {"base": "n"})

    def do_delete():
        cdata = st.session_state.all_contigs[sel]["consensus"]
        if len(cdata) > 1:
            cdata.pop(tpos)
            if st.session_state.wf_pos >= len(cdata):
                st.session_state.wf_pos = len(cdata) - 1

    def zoom_out():
        st.session_state.wf_zoom = min(max_zoom, st.session_state.wf_zoom + 100)

    def zoom_in():
        st.session_state.wf_zoom = max(50, st.session_state.wf_zoom - 100)

    # ================================================================
    # Section 1: Navigation & Position
    # ================================================================
    st.markdown(f"### 🔍 ナビゲーション — **{sel}**")
    st.caption("💡 ◀ ▶ でミスマッチ間を移動。スライダーで任意位置にジャンプ。ズームで表示範囲を調整。")

    nav_row = st.columns([1, 1, 3, 1, 1])

    with nav_row[0]:
        st.button("◀ Prev", key="kb_prev", on_click=do_prev, use_container_width=True,
                  help="前のミスマッチへ (← キー)")
    with nav_row[1]:
        st.button("Next ▶", key="kb_next", on_click=do_next, use_container_width=True,
                  help="次のミスマッチへ (→ キー)")
    with nav_row[2]:
        # Position info
        mm_info = f"  |  ⚠️ ミスマッチ: **{len(mms)}** 箇所" if mms else "  |  ✅ ミスマッチなし"
        st.markdown(f"📍 位置: **{tpos+1}** / {rlen}　塩基: **{cbase.upper()}**{mm_info}")
    with nav_row[3]:
        st.button("🔎−", key="kb_zoomout", on_click=zoom_out, use_container_width=True,
                  help="ズームアウト (- キー)")
    with nav_row[4]:
        st.button("🔎+", key="kb_zoomin", on_click=zoom_in, use_container_width=True,
                  help="ズームイン (= キー)")

    # Position slider + Mismatch jump
    pos_col, jump_col = st.columns([3, 1])
    with pos_col:
        st.slider("位置", 0, rlen-1, key="wf_pos", label_visibility="collapsed")
    with jump_col:
        def on_jump():
            if st.session_state.wf_j is not None:
                st.session_state.wf_pos = st.session_state.wf_j
        if mms:
            st.selectbox("ジャンプ", mms, format_func=lambda x: f"Pos {x+1}",
                         key="wf_j", index=None, on_change=on_jump,
                         label_visibility="collapsed", placeholder="ミスマッチへジャンプ...")
        else:
            st.empty()

    # ================================================================
    # Section 2: Editing Toolbar
    # ================================================================
    st.markdown("### ✏️ 塩基編集")
    st.caption("💡 キーボードで直接塩基を入力するか、ボタンで選択。Autoモードでは編集後に次のミスマッチへ自動移動。")

    edit_row = st.columns([2, 5, 2, 2, 1])

    with edit_row[0]:
        # Key input
        def on_k():
            v = st.session_state.wf_k.upper()
            if v in IUPAC_CODES:
                st.session_state.all_contigs[sel]["consensus"][tpos]["base"] = v.lower()
                st.session_state.wf_k = ""
                if st.session_state.wf_auto:
                    np_pos = f_next(tpos, "n")
                    st.session_state.wf_pos = np_pos if np_pos != tpos else min(rlen-1, tpos+1)
        st.text_input("塩基入力", max_chars=1, key="wf_k", on_change=on_k,
                      placeholder="A/T/G/C...", label_visibility="collapsed",
                      help="塩基を1文字入力してEnter。IUPACコードも使用可")

    with edit_row[1]:
        # Base palette buttons — arranged as a single row
        base_cols = st.columns(11)
        all_bases = ['A','T','G','C','N','R','Y','K','M','S','W']
        for i, b in enumerate(all_bases):
            if base_cols[i].button(b, key=f"kb_base_{b.lower()}",
                                   help=f"{b} を入力", use_container_width=True):
                do_edit(b)

    with edit_row[2]:
        st.button("＋ 挿入", key="kb_ins", on_click=do_insert, use_container_width=True,
                  help="カーソル位置にN塩基を挿入 (+ キー)")

    with edit_row[3]:
        st.button("✕ 削除", key="kb_del", on_click=do_delete, use_container_width=True,
                  help="カーソル位置の塩基を削除 (Delete キー)")

    with edit_row[4]:
        st.checkbox("Auto", key="wf_auto", help="編集後に次のミスマッチへ自動移動",
                    label_visibility="collapsed")
        st.caption("Auto")

    # Recalculate after possible insert/delete
    cons_data = st.session_state.all_contigs[sel]["consensus"]
    rlen = len(cons_data)
    tpos = max(0, min(rlen - 1, st.session_state.wf_pos))
    vis_res = [r for r in results if r["name"] in st.session_state.wf_sel_tr]
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

    # ================================================================
    # Section 3: Waveform Chart
    # ================================================================
    st.markdown("### 📊 波形チャート")
    st.caption("💡 チャートをクリックで位置選択。マウスホイールで一時的ズーム（スライダーのズームのみ保持されます）。")

    fig = create_main_figure(vis_res, align_ref, cons_data, mms, tpos,
                             st.session_state.wf_zoom, orf, trans, tid, qt)

    ev = st.plotly_chart(fig, use_container_width=True, on_select="rerun",
                         selection_mode="points",
                         config={'scrollZoom': True, 'displayModeBar': True})

    # Handle click events
    if ev and ev["selection"]["points"]:
        cx = int(round(ev["selection"]["points"][0]["x"]))
        rt = align_ref["display"]
        if rt["peak_locations"]:
            dst = [abs(p-cx) for p in rt["peak_locations"]]
            nr = dst.index(min(dst))
            if nr != st.session_state.wf_pos:
                st.session_state.wf_pos = nr
                st.rerun()

    # ================================================================
    # Section 4: Inline Sequence Editor
    # ================================================================
    with st.expander("📝 Sequence Editor", expanded=True):
        st.caption("💡 カーソル周辺の配列を表示。🟠カーソル　🔴ミスマッチ　🟣編集済み — ホバーで位置番号を確認")

        seq_win_half = 50
        seq_start = max(0, tpos - seq_win_half)
        seq_end = min(rlen, tpos + seq_win_half + 1)

        html_parts = []
        html_parts.append(
            '<div style="font-family:\'Courier New\',monospace;font-size:13px;'
            'line-height:1.8;word-break:break-all;background:#1a1a2e;color:#d4d4d4;'
            'padding:14px 16px;border-radius:10px;overflow-x:auto;'
            'border:1px solid #333;">'
        )

        # Position ruler
        ruler_parts = []
        for idx in range(seq_start, seq_end):
            if (idx + 1) % 10 == 0:
                num_str = str(idx + 1)
                ruler_parts.append(
                    f'<span style="color:#666;font-size:9px;letter-spacing:-1px;">{num_str}</span>'
                )
            elif (idx + 1) % 5 == 0:
                ruler_parts.append('<span style="color:#444;font-size:9px;">·</span>')
            else:
                ruler_parts.append('<span style="font-size:9px;">&nbsp;</span>')
        html_parts.append('<div style="margin-bottom:2px;">' + ''.join(ruler_parts) + '</div>')

        # Base characters
        base_colors_map = {
            'a': '#4ec94e', 't': '#e85555', 'c': '#5599e8', 'g': '#cccccc',
            'A': '#4ec94e', 'T': '#e85555', 'C': '#5599e8', 'G': '#cccccc',
            'N': '#777777', 'n': '#777777', '-': '#555555',
        }
        mms_set = set(mms)  # faster lookup

        for idx in range(seq_start, seq_end):
            base = cons_data[idx]["base"]
            display_base = base.upper()
            is_cursor = (idx == tpos)
            is_mm = (idx in mms_set)
            is_edited = base.islower()

            color = base_colors_map.get(base, '#d4d4d4')
            bg = 'transparent'
            fw = 'normal'
            extra = ''

            if is_cursor:
                bg = '#e67e22'
                color = '#fff'
                fw = 'bold'
                extra = 'text-decoration:underline;'
            elif is_mm:
                bg = 'rgba(231,76,60,0.35)'
                color = '#ff7b7b'
                fw = 'bold'
            elif is_edited:
                bg = 'rgba(180,60,200,0.25)'
                color = '#d78bff'

            style = (f'background:{bg};color:{color};font-weight:{fw};'
                     f'padding:1px 2px;border-radius:2px;{extra}')
            html_parts.append(
                f'<span style="{style}" title="Pos {idx+1}: {display_base}">{display_base}</span>'
            )

        html_parts.append('</div>')
        st.html(''.join(html_parts))
        st.caption(f"表示範囲: {seq_start+1}–{seq_end} / {rlen} bp")

    # ================================================================
    # Section 5: Genetic Analysis & Stop Codon Report
    # ================================================================
    if trans:
        st.divider()
        st.markdown(f"### 🧬 遺伝子解析 (コドンテーブル: {cname})")
        st.caption("💡 3フレームの終止コドンを検出。ミスマッチが含まれるコドンには ⚠️ マーク。クリックで位置にジャンプ。")
        from modules.bio_logic import calculate_translation_map

        found_issues = False
        stop_codon_options = {}

        cols = st.columns(3)
        for i, frame in enumerate([1, 2, 3]):
            with cols[i]:
                aa_map = calculate_translation_map(cons_data, tid, frame)
                stops = [a for a in aa_map if a["is_stop"]]

                if stops:
                    st.error(f"Frame +{frame}: {len(stops)} 個の終止コドン")
                    found_issues = True
                    for s in stops:
                        mid_idx = s["seq_idx"]
                        codon_rng = range(mid_idx-1, mid_idx+2)
                        related_mms = [m for m in mms if m in codon_rng]
                        status_msg = f"Mismatch@{related_mms[0]+1}" if related_mms else "Stop"
                        label = f"Fr+{frame}: Pos {mid_idx+1} ({status_msg})"
                        stop_codon_options[label] = mid_idx
                else:
                    st.success(f"Frame +{frame}: Clean ✓")

        # Store in session_state for safe callback access
        st.session_state.wf_stop_opts = stop_codon_options

        if found_issues:
            # Callback — reads from session_state to avoid KeyError
            def on_stop_select():
                selected_label = st.session_state.wf_stop_sel
                if selected_label:
                    opts = st.session_state.get("wf_stop_opts", {})
                    if selected_label in opts:
                        st.session_state.wf_pos = opts[selected_label]

            st.selectbox(
                "終止コドンへジャンプ:",
                options=list(stop_codon_options.keys()),
                index=None,
                key="wf_stop_sel",
                on_change=on_stop_select,
                placeholder="終止コドンを選択..."
            )

            # Fix suggestions for current position
            current_stop_target = None
            current_stop_frame = None
            for frame in [1, 2, 3]:
                aa_map = calculate_translation_map(cons_data, tid, frame)
                stops = [a for a in aa_map if a["is_stop"]]
                for s in stops:
                    mid_idx = s["seq_idx"]
                    if abs(tpos - mid_idx) <= 1:
                        current_stop_target = mid_idx
                        current_stop_frame = frame
                        break
                if current_stop_target is not None:
                    break

            if current_stop_target is not None:
                with st.expander(f"💡 修正候補 (Frame +{current_stop_frame})", expanded=True):
                    mid_idx = current_stop_target
                    if 0 < mid_idx < len(cons_data)-1:
                        codon_bases = [
                            cons_data[mid_idx-1]["base"],
                            cons_data[mid_idx]["base"],
                            cons_data[mid_idx+1]["base"],
                        ]
                        codon_str = "".join(codon_bases).upper()
                        st.write(f"対象コドン: **{codon_str}** (位置 {mid_idx+1})")

                        if "N" in codon_str or "-" in codon_str:
                            st.warning("コドンにギャップまたはNが含まれています。")
                        else:
                            from modules.bio_logic import suggest_stop_fixes
                            fixes = suggest_stop_fixes(codon_str, tid)
                            if fixes:
                                st.caption("1塩基の変更で終止コドンを解消できる候補:")
                                cols_fix = st.columns(min(4, len(fixes)))
                                for i, fix in enumerate(fixes):
                                    c = cols_fix[i % len(cols_fix)]
                                    abs_pos = mid_idx - 1 + fix['pos']
                                    is_mm = abs_pos in mms_set
                                    msg = f"{fix['from']}→**{fix['to']}** ({fix['aa']})"
                                    if is_mm:
                                        msg += " ⚠️"
                                    if c.button(msg, key=f"fix_{mid_idx}_{i}",
                                                help=f"位置 {abs_pos+1} を {fix['to']} に変更"):
                                        st.session_state.all_contigs[sel]["consensus"][abs_pos]["base"] = fix['to'].lower()
                                        st.session_state.wf_pos = abs_pos
                                        st.session_state.wf_k = ""
                                        st.rerun()
                            else:
                                st.info("1塩基の変異では解消できない終止コドンです。")

    # ================================================================
    # Section 6: Export
    # ================================================================
    st.divider()
    st.markdown("### 💾 保存 & エクスポート")

    exp_col1, exp_col2 = st.columns([1, 3])
    with exp_col1:
        seq = "".join([d["base"] for d in cons_data]).replace("-", "")
        st.download_button("📥 FASTA保存", f">{sel}\n{seq}", f"{sel}_consensus.fasta",
                           type="primary", use_container_width=True,
                           help="編集済みコンセンサス配列をFASTA形式でダウンロード")
    with exp_col2:
        st.caption(f"配列長: {len(seq)} bp（ギャップ除く）| 編集済み塩基: {sum(1 for d in cons_data if d['base'].islower())} 箇所 | ミスマッチ: {len(mms)} 箇所")

    # ================================================================
    # Inject Keyboard Shortcuts JS
    # ================================================================
    st.html(_KEYBOARD_JS)
