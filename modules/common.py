import shutil
import os
import subprocess
import streamlit as st

def find_tool_path(tool_name):
    """
    ツールの実行パスを探す（Windows/Mac/Linux/Streamlit Cloud完全対応）
    """
    target_names = [tool_name]
    if tool_name == "iqtree2": target_names.append("iqtree")
    if tool_name == "mafft": target_names.append("mafft.bat")

    # 1. システムパス (shutil.which)
    for name in target_names:
        path = shutil.which(name)
        if path: return path

    # 2. Streamlit Cloud / Conda 環境の絶対パス探索 (【重要】ここを追加)
    # 環境変数PATHに入っていなくても、ここにあれば使う
    cloud_paths = [
        f"/home/appuser/.conda/bin/{tool_name}",
        f"/home/adminuser/.conda/bin/{tool_name}",
        f"/opt/conda/bin/{tool_name}",
        f"/usr/local/bin/{tool_name}",
        f"/usr/bin/{tool_name}",
        os.path.expanduser(f"~/.conda/bin/{tool_name}")
    ]
    # iqtree2の場合は iqtree も探すなど柔軟に
    expanded_paths = []
    for p in cloud_paths:
        expanded_paths.append(p)
        if "iqtree2" in p: expanded_paths.append(p.replace("iqtree2", "iqtree"))

    for p in expanded_paths:
        if os.path.exists(p): return p

    # 3. toolsフォルダ探索 (ローカル用)
    base_dir = os.getcwd()
    search_dirs = [os.path.join(base_dir, "tools"), base_dir]
    for search_dir in search_dirs:
        if os.path.exists(search_dir):
            for root, dirs, files in os.walk(search_dir):
                for name in target_names:
                    candidates = [name, name + ".exe", name + ".bat", name + ".cmd"]
                    for f in files:
                        if f.lower() in candidates:
                            return os.path.join(root, f)
    
    return None # 見つからなければ None

def run_command(cmd, **kwargs):
    """
    外部コマンド実行（シンプル版）
    """
    # 競合する引数を削除
    kwargs.pop('stdout', None)
    kwargs.pop('stderr', None)
    kwargs['capture_output'] = True

    # コマンドの先頭が None (見つかっていない) の場合のガード
    if cmd[0] is None:
        raise FileNotFoundError(f"ツールが見つかりません。設定を確認してください。")

    try:
        return subprocess.run(
            cmd, 
            text=True, 
            encoding='utf-8', 
            errors='replace', 
            **kwargs
        )
    except FileNotFoundError:
         raise RuntimeError(f"コマンドが見つかりません: {cmd[0]}")
    except Exception as e:
        raise RuntimeError(f"実行エラー: {cmd[0]}\n{e}")

def generate_alignment_html_from_df(df, max_seqs=50, display_width=80):
    """HTML生成"""
    if df is None or df.empty: return "<p>No sequences.</p>"
    colors = {'A':'#ffc7ce','C':'#c7e5ff','G':'#ffebc7','T':'#d4ffc7','-':'#f0f0f0'}
    html = '<div style="font-family: monospace; overflow-x: auto; white-space: nowrap;">'
    
    if "Include" in df.columns:
        target_df = df[df["Include"] == True].head(max_seqs)
    else:
        target_df = df.head(max_seqs)
        
    for index, row in target_df.iterrows():
        seq_id = str(row.get("ID", ""))[0:15]
        seq = str(row.get("Sequence", "")).upper()[0:display_width]
        row_html = f'<div style="margin: 2px;"><span style="display:inline-block;width:120px;">{seq_id}</span>'
        for char in seq:
            bg = colors.get(char, '#fff')
            row_html += f'<span style="background:{bg};display:inline-block;width:10px;text-align:center;">{char}</span>'
        html += row_html + '</div>'
    html += '</div>'
    return html
