import shutil
import os
import subprocess
import sys
import streamlit as st

def find_tool_path(tool_name):
    """
    ツールの実行パスを探す（Windows/Mac/Linux/Streamlit Cloud完全対応）
    """
    target_names = [tool_name]
    if tool_name == "iqtree2": target_names.append("iqtree")
    if tool_name == "mafft": target_names.append("mafft.bat")

    # 1. 現在のPython環境のbinフォルダ (Streamlit Cloudではここにあるはず)
    current_python_bin = os.path.dirname(sys.executable)
    
    # 2. 探索パスリスト
    search_paths = [
        current_python_bin, # 最優先
        f"/home/appuser/.conda/bin",
        f"/home/adminuser/.conda/bin",
        f"/opt/conda/bin",
        f"/usr/local/bin",
        f"/usr/bin",
        os.path.expanduser(f"~/.conda/bin"),
        # Local tools folder
        os.path.join(os.getcwd(), "tools"),
        os.path.join(os.getcwd(), "tools", f"{tool_name}-win"),
        os.path.join(os.getcwd(), "tools", f"{tool_name}-win", "bin"),
        os.path.join(os.getcwd(), "tools", "trimal"),
        os.path.join(os.getcwd(), "tools", "iqtree2", "bin"),
    ]

    for base_path in search_paths:
        if os.path.exists(base_path):
            for name in target_names:
                full_path = os.path.join(base_path, name)
                if os.path.exists(full_path):
                    return full_path
                # Windows exe check
                if os.path.exists(full_path + ".exe"):
                    return full_path + ".exe"

    # 3. システムパス (shutil.which)
    for name in target_names:
        path = shutil.which(name)
        if path: return path
    
    return None

def run_command(cmd, **kwargs):
    """
    外部コマンド実行
    """
    executable = cmd[0]
    if not executable:
        raise FileNotFoundError("ツールが見つかりません。設定画面のパスを確認してください。")

    # 競合回避: ファイル出力(stdout)指定がある場合はcapture_outputしない
    if 'stdout' not in kwargs and 'stderr' not in kwargs:
        kwargs['capture_output'] = True
    elif 'capture_output' in kwargs:
        del kwargs['capture_output']

    try:
        return subprocess.run(
            cmd, 
            text=True, 
            encoding='utf-8', 
            errors='replace', 
            **kwargs
        )
    except FileNotFoundError:
         raise RuntimeError(f"コマンドが見つかりません: {executable}\n環境設定(environment.yml)またはパス設定を確認してください。")
    except Exception as e:
        raise RuntimeError(f"実行エラー: {executable}\n{e}")

def generate_alignment_html_from_df(df, max_seqs=50, display_width=80, show_dots=False, reference_seq=None):
    """アラインメント結果のHTML生成 (ドット表示対応)"""
    if df is None or df.empty: return "<p>No sequences.</p>"
    
    colors = {'A':'#ffc7ce','C':'#c7e5ff','G':'#ffebc7','T':'#d4ffc7','-':'#f0f0f0','N':'#e0e0e0'}
    
    html = '<div style="font-family: Consolas, monospace; line-height: 1.2; overflow-x: auto; white-space: nowrap; background-color: #fff; padding: 10px; border: 1px solid #ddd; border-radius: 5px;">'
    
    if "Include" in df.columns:
        target_df = df[df["Include"] == True].head(max_seqs)
    else:
        target_df = df.head(max_seqs)

    # リファレンス配列の取得
    if show_dots and reference_seq is None and not target_df.empty:
        reference_seq = str(target_df.iloc[0].get("Sequence", "")).upper()

    for index, row in target_df.iterrows():
        seq_id = str(row.get("ID", ""))[:20]
        if len(str(row.get("ID", ""))) > 23: seq_id += "..."
        
        seq_str = str(row.get("Sequence", "")).upper()
        
        # ドット変換
        display_seq = ""
        if show_dots and reference_seq and index > 0:
            for i, char in enumerate(seq_str):
                if i < len(reference_seq) and char == reference_seq[i] and char not in ['-', 'N']:
                    display_seq += '.'
                else:
                    display_seq += char
        else:
            display_seq = seq_str
            
        display_seq = display_seq[:display_width]
        
        row_html = f'<div style="margin: 2px;"><span style="display:inline-block;width:150px;font-size:12px;font-weight:bold;">{seq_id}</span>'
        
        for i, char in enumerate(display_seq):
            orig = seq_str[i] if i < len(seq_str) else '-'
            bg = colors.get(orig, '#fff')
            style = 'color: #999; font-weight: bold;' if char == '.' else 'color: #000;'
            if char == '.': bg = '#fff'
            row_html += f'<span style="background:{bg};{style}display:inline-block;width:10px;text-align:center;font-size:12px;">{char}</span>'
        html += row_html + '</div>'
    
    html += '</div>'
    return html
