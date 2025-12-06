import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from Bio import SeqIO

def run_simple_asap_logic(aligned_fasta_path, threshold=0.02):
    """
    簡易的な距離ベースの種区分解析 (Barcode Gap Analysis) のロジック
    """
    # 配列読み込み
    seqs = list(SeqIO.parse(aligned_fasta_path, "fasta"))
    if len(seqs) < 3: return None, "配列数が少なすぎます（最低3配列必要）"
    
    ids = [s.id for s in seqs]
    
    # 配列をnumpy配列化 (文字コードに変換)
    max_len = max(len(s.seq) for s in seqs)
    # パディングして長さを揃える
    matrix = []
    for s in seqs:
        seq_str = str(s.seq).upper().ljust(max_len, '-')
        matrix.append(np.array(list(seq_str)))
    matrix = np.array(matrix)
    
    # 距離行列計算 (p-distance: 異なる塩基の割合)
    # 高速化のために単純なループで実装（大規模データだと遅くなる可能性あり）
    n = len(seqs)
    dist_matrix = np.zeros((n, n))
    
    for i in range(n):
        for j in range(i+1, n):
            s1 = matrix[i]
            s2 = matrix[j]
            
            # 比較可能なサイト（両方ともギャップやNでない）を抽出
            valid_mask = (s1 != '-') & (s2 != '-') & (s1 != 'N') & (s2 != 'N')
            valid_count = np.sum(valid_mask)
            
            if valid_count == 0:
                dist = 0.0 # 比較不能なら距離0とする（あるいは1.0）
            else:
                diff_count = np.sum((s1[valid_mask] != s2[valid_mask]))
                dist = diff_count / valid_count
                
            dist_matrix[i, j] = dist_matrix[j, i] = dist
    
    # クラスタリング (UPGMA: method='average')
    # 距離行列はcondensed formである必要があるためsquareformを使用
    Z = linkage(squareform(dist_matrix), method='average')
    
    # 閾値でグループ分け
    clusters = fcluster(Z, t=threshold, criterion='distance')
    
    # 結果整形
    result_df = pd.DataFrame({"ID": ids, "Cluster": clusters})
    result_df = result_df.sort_values("Cluster")
    
    return result_df, dist_matrix
