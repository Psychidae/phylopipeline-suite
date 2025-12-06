import shutil
import os
import subprocess

def find_tool_path(tool_name):
    """
    指定されたツール名(例: 'iqtree2')を
    システムパスまたは 'tools' フォルダ内から再帰的に探索して返す
    """
    target_names = [tool_name]
    if tool_name == "iqtree2": target_names.append("iqtree")
    if tool_name == "mafft": target_names.append("mafft.bat")

    # 1. システムパス
    for name in target_names:
        path = shutil.which(name)
        if path: return path

    # 2. toolsフォルダ探索
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
    return None

def run_command(cmd, **kwargs):
    """外部コマンド実行（Windowsの文字化け対策込み）"""
    return subprocess.run(
        cmd, 
        capture_output=True, 
        text=True, 
        encoding='utf-8', 
        errors='replace', 
        **kwargs
    )

def generate_alignment_html_from_df(df, max_seqs=50, display_width=80):
    """系統解析用: アラインメント結果の簡易HTML生成"""
    if df.empty: return "<p>No sequences.</p>"
    colors = {'A':'#ffc7ce','C':'#c7e5ff','G':'#ffebc7','T':'#d4ffc7','-':'#f0f0f0'}
    html = '<div style="font-family: monospace; overflow-x: auto; white-space: nowrap;">'
    
    if "Include" in df.columns:
        target_df = df[df["Include"] == True].head(max_seqs)
    else:
        target_df = df.head(max_seqs)
        
    for index, row in target_df.iterrows():
        seq_id = str(row["ID"])[:15]
        seq = str(row["Sequence"]).upper()[:display_width]
        row_html = f'<div style="margin: 2px;"><span style="display:inline-block;width:120px;">{seq_id}</span>'
        for char in seq:
            bg = colors.get(char, '#fff')
            row_html += f'<span style="background:{bg};display:inline-block;width:10px;text-align:center;">{char}</span>'
        html += row_html + '</div>'
    html += '</div>'
    return html
