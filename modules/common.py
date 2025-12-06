import shutil
import os
import subprocess
import sys

def find_tool_path(tool_name):
    """ツールの実行パスを探す（環境変数PATHへの追加も行う）"""
    # 探索パス候補
    search_paths = [
        os.path.dirname(sys.executable),
        "/home/appuser/.conda/bin",
        "/home/adminuser/.conda/bin",
        "/opt/conda/bin",
        "/usr/local/bin"
    ]
    
    # 環境変数PATHにこれらがなければ追加しておく（これでshutil.whichで見つかるようになる）
    current_path = os.environ.get("PATH", "")
    for p in search_paths:
        if p not in current_path and os.path.exists(p):
            os.environ["PATH"] = f"{p}:{os.environ['PATH']}"

    # バリエーション検索
    targets = [tool_name]
    if tool_name == "mafft": targets.append("mafft.bat")
    if tool_name == "iqtree2": targets.append("iqtree")

    for t in targets:
        path = shutil.which(t)
        if path: return path
        
    return None

def run_command(cmd, **kwargs):
    """
    外部コマンド実行（MAFFTエラー対策版）
    """
    executable = cmd[0]
    
    # パス解決を試みる
    if not shutil.which(executable):
        resolved = find_tool_path(executable)
        if resolved:
            cmd[0] = resolved
        else:
            # 見つからなくても、パスが通っていることを期待してそのまま実行させる
            pass

    # 【重要】Streamlit CloudでのMAFFT Permission denied対策
    # stderrを subprocess.DEVNULL に捨てることで、書き込みエラーを回避する
    if 'stderr' not in kwargs:
        kwargs['stderr'] = subprocess.DEVNULL

    # stdoutが指定されていない場合はキャプチャする
    if 'stdout' not in kwargs:
        kwargs['capture_output'] = False # stdout/stderrを個別に指定するためFalse
        kwargs['stdout'] = subprocess.PIPE
        # stderrは上でDEVNULLに設定済み

    try:
        return subprocess.run(
            cmd, 
            text=True, 
            encoding='utf-8', 
            errors='replace', 
            **kwargs
        )
    except Exception as e:
        raise RuntimeError(f"実行エラー: {cmd[0]}\n{e}")

# (以下 generate_alignment_html_from_df は変更なし)
def generate_alignment_html_from_df(df, max_seqs=50, display_width=80, show_dots=False, reference_seq=None):
    if df is None or df.empty: return "<p>No sequences.</p>"
    colors = {'A':'#ffc7ce','C':'#c7e5ff','G':'#ffebc7','T':'#d4ffc7','-':'#f0f0f0','N':'#e0e0e0'}
    html = '<div style="font-family: Consolas, monospace; line-height: 1.2; overflow-x: auto; white-space: nowrap; background-color: #fff; padding: 10px; border: 1px solid #ddd; border-radius: 5px;">'
    
    target_df = df[df["Include"] == True].head(max_seqs) if "Include" in df.columns else df.head(max_seqs)
    if show_dots and reference_seq is None and not target_df.empty:
        reference_seq = str(target_df.iloc[0].get("Sequence", "")).upper()
        
    for index, row in target_df.iterrows():
        seq_id = str(row.get("ID", ""))[:20]
        if len(str(row.get("ID", ""))) > 23: seq_id += "..."
        seq_str = str(row.get("Sequence", "")).upper()[:display_width]
        
        display_seq = ""
        if show_dots and reference_seq and index > 0:
            for i, char in enumerate(seq_str):
                if i < len(reference_seq) and char == reference_seq[i] and char not in ['-', 'N']: display_seq += '.'
                else: display_seq += char
        else: display_seq = seq_str
            
        row_html = f'<div style="margin: 2px;"><span style="display:inline-block;width:150px;font-size:12px;font-weight:bold;">{seq_id}</span>'
        for char in display_seq:
            bg = colors.get(seq_str[display_seq.index(char)] if char == '.' else char, '#fff')
            if char == '.': bg = '#fff'; style = 'color: #999; font-weight: bold;'
            else: style = 'color: #000;'
            row_html += f'<span style="background:{bg};{style}display:inline-block;width:10px;text-align:center;font-size:12px;">{char}</span>'
        html += row_html + '</div>'
    html += '</div>'
    return html
