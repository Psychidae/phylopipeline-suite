import os
import re
import pandas as pd

def guess_group_name(filename):
    """ファイル名からサンプル名（グループ）を推測する"""
    name = os.path.splitext(filename)[0]
    
    # 特定キーワードがあればその手前まで
    known_tags = r'(HCO|LCO|COI|ITS|16S|18S|28S|M13|VF|VR|T7|SP6)'
    match = re.search(f'(.+?)[_\.-]+{known_tags}', name, re.IGNORECASE)
    if match:
        return match.group(1).rstrip('_- .')

    # 汎用タグ除去
    clean_name = name
    clean_name = re.sub(r'[_\.-]?(fwd|rev|forward|reverse|f|r)$', '', clean_name, flags=re.IGNORECASE)
    clean_name = re.sub(r'\s*\([0-9]+\)$', '', clean_name)
    
    while True:
        prev = clean_name
        clean_name = re.sub(r'[_\.-][a-zA-Z][a-zA-Z0-9\+\-\.]*$', '', clean_name)
        if clean_name == prev: break
            
    clean_name = clean_name.rstrip('_- .')
    if not clean_name: return name[:10]
    return clean_name

def create_grouping_dataframe(uploaded_files):
    """グルーピング用DataFrame作成"""
    data = []
    for f in uploaded_files:
        g_name = guess_group_name(f.name)
        data.append({"Filename": str(f.name), "Group": str(g_name)})
    
    df = pd.DataFrame(data, columns=["Filename", "Group"])
    df["Filename"] = df["Filename"].astype(str)
    df["Group"] = df["Group"].astype(str)
    return df
