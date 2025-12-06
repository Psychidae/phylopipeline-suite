import shutil
import os
import subprocess
import sys

def find_tool_path(tool_name):
    """ツールの実行パスを探す（Windows/Mac/Linux対応）"""
    # 1. システムパス (Cloud環境やHomebrew/Condaはこれで見つかる)
    path = shutil.which(tool_name)
    if path: return path
    
    # 2. 一般的なパス (Cloud/Local)
    common_paths = [
        f"/opt/homebrew/bin/{tool_name}",
        f"/usr/local/bin/{tool_name}",
        f"/usr/bin/{tool_name}",
        # Streamlit Cloud / Conda paths
        f"/home/appuser/.conda/bin/{tool_name}",
        f"/home/adminuser/.conda/bin/{tool_name}",
        os.path.expanduser(f"~/.conda/bin/{tool_name}"),
        os.path.expanduser(f"~/anaconda3/bin/{tool_name}"),
        os.path.expanduser(f"~/miniconda3/bin/{tool_name}")
    ]
    
    if os.name == 'nt': # Windows
         common_paths.extend([
            f"C:\\Program Files\\{tool_name}\\{tool_name}.exe",
            f"C:\\BioTools\\{tool_name}\\{tool_name}.exe",
            # toolsフォルダ内の検索
            os.path.join(os.getcwd(), "tools", f"{tool_name}-win", f"{tool_name}.bat"), 
            os.path.join(os.getcwd(), "tools", f"{tool_name}-win", "bin", f"{tool_name}.exe"),
            os.path.join(os.getcwd(), "tools", "trimal", "trimal.exe"),
            os.path.join(os.getcwd(), "tools", "iqtree2", "bin", "iqtree2.exe"),
        ])
    
    for p in common_paths:
        if os.path.exists(p): return p
    return None

def run_command(cmd, **kwargs):
    """
    外部コマンド実行（Windowsの文字化け対策込み）
    stdout/stderr引数がある場合はcapture_outputを無効化して競合を防ぐ
    """
    # ファイル出力指定がない場合のみ、結果を変数にキャプチャする
    if 'stdout' not in kwargs and 'stderr' not in kwargs:
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
        # 実行時に見つからなかった場合
        raise FileNotFoundError(
            f"ツールが見つかりません: '{cmd[0]}'\n"
            "Cloud環境の場合は environment.yml のインストール待ちか、設定を確認してください。"
        )
    except Exception as e:
        raise Exception(f"コマンド実行エラー: {e}")

def generate_alignment_html_from_df(df, max_seqs=50, display_width=80):
    """系統解析用: アラインメント結果の簡易HTML生成"""
    if df.empty: return "<p>No sequences to display.</p>"
    colors = {'A':'#ffc7ce','C':'#c7e5ff','G':'#ffebc7','T':'#d4ffc7','U':'#d4ffc7','-':'#f0f0f0','N':'#e0e0e0'}
    html = '<div style="font-family: Consolas, monospace; line-height: 1.2; overflow-x: auto; white-space: nowrap; background-color: #fff; padding: 10px; border: 1px solid #ddd; border-radius: 5px;">'
    
    # Includeカラムがある場合はフィルタリング
    if "Include" in df.columns:
        target_df = df[df["Include"] == True].head(max_seqs)
    else:
        target_df = df.head(max_seqs)
        
    for index, row in target_df.iterrows():
        seq_id = str(row["ID"])[:17] + "..." if len(str(row["ID"])) > 20 else str(row["ID"])
        seq_str = str(row["Sequence"]).upper()[:display_width]
        row_html = f'<div style="margin-bottom: 2px;"><span style="display: inline-block; width: 160px; font-weight: bold; color: #333; font-size: 12px;">{seq_id}</span>'
        for char in seq_str:
            bg = colors.get(char, '#ffffff')
            row_html += f'<span style="background-color: {bg}; color: #000; display: inline-block; width: 12px; text-align: center; margin: 0; font-size: 12px;">{char}</span>'
        html += row_html + '</div>'
    html += '</div>'
    return html
