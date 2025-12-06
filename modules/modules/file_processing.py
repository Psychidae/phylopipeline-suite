import os
import re
import pandas as pd

def guess_group_name(filename):
    """
    ファイル名からサンプル名（グループ）を推測する。
    拡張子やプライマー名、方向タグを取り除いた共通部分を抽出する。
    """
    # 1. 拡張子除去
    name = os.path.splitext(filename)[0]
    
    # 2. 【強力な判定】特定の頻出キーワードがあれば、その手前までを採用
    # HCO/LCO 以外にも、COI, ITS, 16S, M13 などの一般的なマーカー名に対応
    # 例: "Sample01_HCO-HCO" -> "Sample01"
    #     "MyData_ITS1_F" -> "MyData"
    known_tags = r'(HCO|LCO|COI|ITS|16S|18S|28S|M13|VF|VR|T7|SP6)'
    match = re.search(f'(.+?)[_\.-]+{known_tags}', name, re.IGNORECASE)
    if match:
        return match.group(1).rstrip('_- .')

    # 3. 【汎用的な判定】末尾のタグ除去
    # 以下のステップで末尾から順に不要な情報を削っていく
    
    clean_name = name
    
    # (A) 方向タグの削除 (F, R, Fwd, Rev, etc.)
    clean_name = re.sub(r'[_\.-]?(fwd|rev|forward|reverse|f|r)$', '', clean_name, flags=re.IGNORECASE)
    
    # (B) 末尾の数字連番の削除 (Sample_001 の _001 を消すかどうかはケースバイケースだが、ペアリングのために一旦消してみる)
    # ここでは「(1)」のようなコピー表記は消すが、IDとしての数字は残したい場合が多い。
    clean_name = re.sub(r'\s*\([0-9]+\)$', '', clean_name) # "Sample (1)" -> "Sample"
    
    # (C) 汎用プライマー名除去ロジック
    # 「区切り文字 + 英字で始まる英数字」が末尾にあれば、それはプライマー名やIDタグとみなして削除
    # 例: "Sample_01_PrimerA" -> "Sample_01"
    #     "Sample_01_001" -> "Sample_01_001" (数字始まりはサンプルIDの可能性が高いので消さない)
    # 繰り返し適用して、複数のタグ（_T7_F など）に対応
    while True:
        prev_name = clean_name
        # _ または - または . の後に、英字で始まる文字列が続く場合
        clean_name = re.sub(r'[_\.-][a-zA-Z][a-zA-Z0-9\+\-\.]*$', '', clean_name)
        
        # 変化がなければ終了
        if clean_name == prev_name:
            break
            
    # 末尾の記号削除
    clean_name = clean_name.rstrip('_- .')
    
    # 4. 安全策 (名前が空になってしまった場合や、短すぎる場合への対処)
    if not clean_name:
        return name[:10]
        
    return clean_name

def create_grouping_dataframe(uploaded_files):
    """アップロードされたファイルリストからグルーピング用のDataFrameを作成"""
    data = []
    for f in uploaded_files:
        g_name = guess_group_name(f.name)
        data.append({
            "Filename": str(f.name),
            "Group": str(g_name)
        })
    
    # DataFrame作成 (カラム名と型を明示)
    df = pd.DataFrame(data, columns=["Filename", "Group"])
    df["Filename"] = df["Filename"].astype(str)
    df["Group"] = df["Group"].astype(str)
    
    return df