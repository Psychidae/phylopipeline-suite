import shutil
import os
import subprocess
import sys

def find_tool_path(tool_name):
    """ツールの実行パスを探す（Windows/Mac/Linux対応）"""
    # ターゲット名のバリエーション
    target_names = [tool_name]
    if tool_name == "iqtree2":
        target_names.append("iqtree")
    if tool_name == "mafft":
        target_names.append("mafft.bat") # Windows用

    # 1. システムパス (PATH) から探す
    for name in target_names:
        path = shutil.which(name)
        if path: return path

    # 2. カレントディレクトリの 'tools' フォルダ内をくまなく探す
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
    """
    外部コマンド実行（Windowsの文字化け対策込み）
    Streamlit Cloudでの権限エラー回避のため、stderrを明示的にPIPEする
    """
    # コマンドの存在チェック（簡易）
    executable = cmd[0]
    if shutil.which(executable) is None and not os.path.exists(executable):
        # パスが通っていない場合でも動くことがあるため、警告しつつ続行
        pass

    # 設定の競合を防ぐロジック
    # stdout/stderrが指定されていない場合は capture_output=True にしたいが、
    # Windowsや特定環境での挙動安定のため、個別にPIPEを指定する方式に変更
    
    if 'stdout' not in kwargs:
        kwargs['stdout'] = subprocess.PIPE
    
    # 【重要】MAFFTの /dev/stderr Permission denied 対策
    # stderrを明示的にキャプチャすることで、ファイルへの直接書き込みを防ぐ
    if 'stderr' not in kwargs:
        kwargs['stderr'] = subprocess.PIPE

    # capture_output引数は使わない（stdout/stderrと競合するため）
    if 'capture_output' in kwargs:
        del kwargs['capture_output']

    try:
        return subprocess.run(
            cmd, 
            text=True, 
            encoding='utf-8', 
            errors='replace', 
            **kwargs
        )
    except Exception as e:
        raise Exception(f"コマンド実行エラー: {e}")

def generate_alignment_html_from_df(df, max_seqs=50, display_width=80):
    """アラインメント結果の簡易HTML生成"""
    # データフレームが空、または必要な列がない場合のガード
    if df is None or df.empty or "Sequence" not in df.columns: 
        return "<p style='color:gray;'>No sequences to display. (Alignment might have failed)</p>"
    
    colors = {
        'A': '#ffc7ce', 'C': '#c7e5ff', 'G': '#ffebc7', 'T': '#d4ffc7',
        'U': '#d4ffc7', '-': '#f0f0f0', 'N': '#e0e0e0'
    }

    html = '<div style="font-family: Consolas, monospace; line-height: 1.2; overflow-x: auto; white-space: nowrap; background-color: #fff; padding: 10px; border: 1px solid #ddd; border-radius: 5px;">'
    
    # Includeカラムがある場合はフィルタリング、なければ全件
    if "Include" in df.columns:
        target_df = df[df["Include"] == True].head(max_seqs)
    else:
        target_df = df.head(max_seqs)
        
    for index, row in target_df.iterrows():
        # ID列があるか確認
        seq_id = str(row.get("ID", f"Seq_{index}"))
        if len(seq_id) > 20: seq_id = seq_id[:17] + "..."
        
        seq_str = str(row["Sequence"]).upper()[:display_width]
        
        row_html = f'<div style="margin-bottom: 2px;"><span style="display: inline-block; width: 160px; font-weight: bold; color: #333; font-size: 12px;">{seq_id}</span>'
        for char in seq_str:
            bg = colors.get(char, '#ffffff')
            row_html += f'<span style="background-color: {bg}; color: #000; display: inline-block; width: 12px; text-align: center; margin: 0; font-size: 12px;">{char}</span>'
        html += row_html + '</div>'

    if len(df) > max_seqs:
        html += f'<div style="color: #888; margin-top: 10px; font-size: 12px;">... 他 {len(df) - max_seqs} 配列は省略 ...</div>'
    
    html += '</div>'
    return html
