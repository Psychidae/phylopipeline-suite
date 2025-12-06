import shutil
import os
import subprocess
import sys

def find_tool_robust(target_names):
    """
    指定されたツール名リスト(例: ['iqtree2', 'iqtree'])のいずれかを
    システムパスまたは 'tools' フォルダ内から再帰的に探索して返す
    """
    # 1. システムパス (PATH) から探す
    # クラウド環境(Streamlit Cloud)ではここで見つかり、即座に終了するため負荷はない
    for name in target_names:
        path = shutil.which(name)
        if path: return path

    # 2. カレントディレクトリの 'tools' フォルダ内をくまなく探す
    # ローカル環境専用のロジック
    base_dir = os.getcwd()
    search_dirs = [os.path.join(base_dir, "tools"), base_dir]
    
    for search_dir in search_dirs:
        if os.path.exists(search_dir):
            for root, dirs, files in os.walk(search_dir):
                for name in target_names:
                    # Windows用に拡張子ありパターンも作成
                    candidates = [name, name + ".exe", name + ".bat", name + ".cmd"]
                    for f in files:
                        if f.lower() in candidates:
                            return os.path.join(root, f)
    return None

def run_command(cmd, **kwargs):
    """
    外部コマンド実行（Windowsの文字化け対策込み）
    """
    return subprocess.run(
        cmd, 
        capture_output=True, 
        text=True, 
        encoding='utf-8', 
        errors='replace', 
        **kwargs
    )

def generate_alignment_html_from_df(df, max_seqs=50, display_width=80):
    """
    データフレームからアラインメント結果の簡易HTMLを生成する関数
    (PhyloPipeline系統解析機能で使用)
    """
    if df.empty: return "<p>No sequences to display.</p>"
    
    colors = {
        'A': '#ffc7ce', 'C': '#c7e5ff', 'G': '#ffebc7', 'T': '#d4ffc7',
        'U': '#d4ffc7', '-': '#f0f0f0', 'N': '#e0e0e0'
    }

    html = '<div style="font-family: Consolas, monospace; line-height: 1.2; overflow-x: auto; white-space: nowrap; background-color: #fff; padding: 10px; border: 1px solid #ddd; border-radius: 5px;">'
    
    # 表示対象のみ抽出 (IncludeがTrueのもの)
    if "Include" in df.columns:
        target_df = df[df["Include"] == True].head(max_seqs)
    else:
        target_df = df.head(max_seqs)
    
    for index, row in target_df.iterrows():
        seq_id = str(row["ID"])
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
