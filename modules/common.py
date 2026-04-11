import shutil
import os
import sys
import subprocess
import streamlit as st

from modules.constants import TOOLS_DIR

def find_tool_path(tool_name):
    """ツールの実行パスを探す（Windows/Mac/Linux対応）"""
    target_names = [tool_name]
    if tool_name == "iqtree2": target_names.append("iqtree")
    if tool_name == "mafft": target_names.append("mafft.bat")

    # 1. Check Standard System Paths (shutil.which)
    for name in target_names:
        path = shutil.which(name)
        if path: return path

    # 2. Check Local TOOLS_DIR (Centralized Logic)
    for name in target_names:
        candidates = [name, name + ".exe", name + ".bat", name + ".cmd"]
        for cand in candidates:
             tool_path = os.path.join(TOOLS_DIR, cand)
             if os.path.exists(tool_path):
                 return tool_path

    # 3. Python Bin (Fallback)
    python_dir = os.path.dirname(sys.executable)
    for name in target_names:
        candidate = os.path.join(python_dir, name)
        if os.path.exists(candidate): return candidate
        if os.path.exists(candidate + ".exe"): return candidate + ".exe"



    # 4. Check Common System Paths (Homebrew/Local)
    extra_paths = ["/opt/homebrew/bin", "/usr/local/bin"]
    for path_dir in extra_paths:
        for name in target_names:
            candidate = os.path.join(path_dir, name)
            if os.path.exists(candidate): return candidate
            if os.path.exists(candidate + ".exe"): return candidate + ".exe"

    return None

def run_command(cmd, **kwargs):
    """外部コマンド実行"""
    
    # stdout/stderrが指定されていない場合のみキャプチャする
    if 'stdout' not in kwargs and 'stderr' not in kwargs:
         kwargs['capture_output'] = True

    if cmd[0] is None:
        raise FileNotFoundError("ツールが見つかりません。")

    try:
        return subprocess.run(
            cmd, text=True, encoding='utf-8', errors='replace', **kwargs
        )
    except FileNotFoundError:
         import os
         raise RuntimeError(f"コマンドが見つかりません: {cmd[0]}\nPATH: {os.environ.get('PATH')}")
    except Exception as e:
        raise RuntimeError(f"実行エラー: {cmd[0]}\n{e}")

def generate_alignment_html_from_df(df, max_seqs=50, display_width=80, show_dots=False, reference_seq=None):
    """
    アラインメント結果のHTML生成（ドット表示モード対応）
    show_dots=True の場合、1行目（リファレンス）と同じ塩基は '.' で表示する
    """
    if df.empty:
        return "<p>No sequences.</p>"

    colors = {'A': '#ffc7ce', 'C': '#c7e5ff', 'G': '#ffebc7', 'T': '#d4ffc7', '-': '#f0f0f0', 'N': '#e0e0e0'}

    target_df = df[df["Include"] == True].head(max_seqs) if "Include" in df.columns else df.head(max_seqs)

    # リファレンス配列の取得（指定がなければ最初の表示行の配列）
    if show_dots and reference_seq is None and not target_df.empty:
        reference_seq = str(target_df.iloc[0]["Sequence"]).upper()

    html = (
        '<div style="font-family: Consolas, monospace; line-height: 1.2; '
        'overflow-x: auto; white-space: nowrap; background-color: #fff; '
        'padding: 10px; border: 1px solid #ddd; border-radius: 5px;">'
    )

    # enumerate で行カウンタを管理（DataFrameインデックスに依存しない）
    for row_counter, (_, row) in enumerate(target_df.iterrows()):
        seq_id_raw = str(row["ID"])
        seq_id = seq_id_raw[:20] + "..." if len(seq_id_raw) > 23 else seq_id_raw
        seq_str = str(row["Sequence"]).upper()

        # ドット変換ロジック（row_counter > 0 = リファレンス以外の行）
        if show_dots and reference_seq and row_counter > 0:
            display_seq = ""
            for i, char in enumerate(seq_str):
                if i < len(reference_seq) and char == reference_seq[i] and char not in ('-', 'N'):
                    display_seq += '.'
                else:
                    display_seq += char
        else:
            display_seq = seq_str

        # 表示幅でカット
        display_seq = display_seq[:display_width]

        row_html = (
            f'<div style="margin: 2px;">'
            f'<span style="display:inline-block;width:150px;font-size:12px;font-weight:bold;">'
            f'{seq_id}</span>'
        )

        for i, char in enumerate(display_seq):
            original_char = seq_str[i] if i < len(seq_str) else '-'
            bg = colors.get(original_char, '#fff')

            if char == '.':
                bg = '#fff'
                char_style = 'color: #999; font-weight: bold;'
            else:
                char_style = 'color: #000;'

            row_html += (
                f'<span style="background:{bg};{char_style}'
                f'display:inline-block;width:10px;text-align:center;font-size:12px;">'
                f'{char}</span>'
            )

        html += row_html + '</div>'

    html += '</div>'
    return html
