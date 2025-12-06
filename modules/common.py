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

    # 1. 現在実行中のPythonと同じ場所(bin)を探す (これがCloudで最も確実)
    # Streamlit Cloudでは /home/adminuser/.conda/bin/python 等で動いている
    current_bin_dir = os.path.dirname(sys.executable)
    
    # 2. 探索パスリスト
    search_paths = [
        current_bin_dir, # 最優先
        f"/home/appuser/.conda/bin",
        f"/home/adminuser/.conda/bin",
        f"/opt/conda/bin",
        f"/usr/local/bin",
        f"/usr/bin",
        # ローカル環境用
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
    
    return None # 見つからなければ None

def run_command(cmd, **kwargs):
    """
    外部コマンド実行
    """
    executable = cmd[0]
    
    # パスが見つかっていない場合のガード
    if not executable:
        # パスがNoneなら、コマンド名そのものを試す（パスが通っていることに賭ける）
        # cmd[0] が Noneだと subprocessで落ちるため、本来のツール名に戻す必要があるが、
        # ここではエラーとして通知する
        raise FileNotFoundError("ツールのパスが特定できませんでした。")

    # 競合回避: ファイル出力(stdout)指定がある場合はcapture_outputしない
    if 'stdout' in kwargs or 'stderr' in kwargs:
        if 'capture_output' in kwargs:
            del kwargs['capture_output']
    else:
        # 指定がない場合はキャプチャする
        kwargs['capture_output'] = True

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
        seq_id = str(row.get("ID", ""))[:15]
        seq = str(row.get("Sequence", "")).upper()[0:display_width]
        row_html = f'<div style="margin: 2px;"><span style="display:inline-block;width:120px;">{seq_id}</span>'
        for char in seq:
            bg = colors.get(char, '#fff')
            row_html += f'<span style="background:{bg};display:inline-block;width:10px;text-align:center;">{char}</span>'
        html += row_html + '</div>'
    html += '</div>'
    return html
