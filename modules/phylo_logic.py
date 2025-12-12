import pandas as pd
import numpy as np
import math
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from Bio import SeqIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo.Consensus import bootstrap_trees, get_support

def run_phylo_bootstrap(msa, method="nj", model="identity", replicates=100):
    """
    Run NJ or UPGMA with bootstrap support.
    
    Args:
        msa: MultipleSeqAlignment object
        method: "nj" or "upgma"
        model: "identity", "kimura80" (k2p), etc. (BioPython models)
        replicates: Number of bootstrap replicates (0 to disable)
        
    Returns:
        tree: Bio.Phylo BaseTree object (with support values if replicates > 0)
    """
    # 1. Calculator Setup
    # BioPython model names: 'identity', 'blastn', 'trans', 'benner6', 'k2p', 'jukes-cantor', 'kimura80'
    # Note: 'kimura80' is alias for 'k2p' in some versions, or vice versa. BioPython usually uses 'identity', 'jukes-cantor', 'kimura80'.
    calculator = DistanceCalculator(model)
    constructor = DistanceTreeConstructor(calculator, method)
    
    # 2. Build Main Tree
    main_tree = constructor.build_tree(msa)
    
    # 3. Bootstrap (if requested)
    if replicates > 0:
        # Generate bootstrap trees
        # bootstrap_trees generates replicates of MSA, then builds trees for each
        boot_trees = bootstrap_trees(msa, replicates, constructor)
        
        # Calculate support
        # get_support maps support values onto the 'target_tree' (main_tree)
        # Returns the tree with confidence values (numbers 0-100 or 0.0-1.0 depending on implementation)
        # BioPython usually calculates percentage (0-100) or fraction. 
        # get_support modifies the tree in-place (naming clades with support) OR returns consensus.
        # We want to map support onto our main topology.
        consensus_tree = get_support(main_tree, boot_trees)
        return consensus_tree
    else:
        return main_tree


def generate_methods_log(tool_versions, params):
    """
    解析手法のログテキストを生成する
    """
    log = []
    log.append("=== PhyloPipeline Analysis Log ===")
    log.append(f"Date: {pd.Timestamp.now()}")
    log.append("")
    log.append("[Tools]")
    for k, v in tool_versions.items():
        log.append(f"{k}: {v}")
    log.append("")
    log.append("[Parameters]")
    for k, v in params.items():
        log.append(f"{k}: {v}")
    log.append("")
    log.append("[Citations]")
    log.append("Please cite the following papers if you use these results:")
    log.append("- MAFFT: Katoh et al. (2013) MSE 30:772-780")
    if params.get("Use trimAl"):
        log.append("- trimAl: Capella-Gutierrez et al. (2009) Bioinformatics 25:1972-1973")
    log.append("- IQ-TREE: Minh et al. (2020) MBE 37:1530-1534")
    log.append("- ModelFinder: Kalyaanamoorthy et al. (2017) Nat Methods 14:587-589")
    log.append("- UFBoot: Hoang et al. (2018) MBE 35:518-522")
    log.append("")
    return "\n".join(log)

def calculate_distance_matrix(aligned_fasta_path, model='p-dist'):
    """
    アラインメント済みFASTAから距離行列を計算する
    model: 'p-dist' (Hamming) or 'k2p' (Kimura 2-Parameter)
    """
    seqs = list(SeqIO.parse(aligned_fasta_path, "fasta"))
    if len(seqs) < 3: return None, None, None, "配列数が少なすぎます（最低3配列必要）"
    
    ids = [s.id for s in seqs]
    
    # 配列をnumpy配列化
    max_len = max(len(s.seq) for s in seqs)
    matrix = []
    for s in seqs:
        seq_str = str(s.seq).upper().ljust(max_len, '-')
        matrix.append(np.array(list(seq_str)))
    matrix = np.array(matrix)
    
    # 距離行列計算
    n = len(seqs)
    dist_matrix = np.zeros((n, n))
    
    # 定義: Transition (A<->G, C<->T)
    def is_transition(n1, n2):
        s = {n1, n2}
        return s == {'A', 'G'} or s == {'C', 'T'}

    for i in range(n):
        for j in range(i+1, n):
            s1 = matrix[i]
            s2 = matrix[j]
            valid_mask = (s1 != '-') & (s2 != '-') & (s1 != 'N') & (s2 != 'N')
            valid_count = np.sum(valid_mask)
            
            if valid_count == 0:
                dist = 0.0
            else:
                mismatches = (s1[valid_mask] != s2[valid_mask])
                diff_count = np.sum(mismatches)
                p_dist = diff_count / valid_count
                
                if model == 'k2p':
                    # K2P: d = -0.5 ln(1 - 2P - Q) - 0.25 ln(1 - 2Q)
                    # P = Transitions / valid, Q = Transversions / valid
                    # ベクトル化は複雑なので簡易ループでP, Q数える（速度改善の余地あり）
                    pairs = np.stack((s1[valid_mask], s2[valid_mask]), axis=1)
                    ts_count = 0
                    tv_count = 0
                    for row in pairs:
                        if row[0] != row[1]:
                            if is_transition(row[0], row[1]):
                                ts_count += 1
                            else:
                                tv_count += 1
                    
                    P = ts_count / valid_count
                    Q = tv_count / valid_count
                    
                    try:
                        w1 = 1 - 2*P - Q
                        w2 = 1 - 2*Q
                        if w1 <= 0 or w2 <= 0:
                            dist = 10.0 # 定義不能なほど遠い
                        else:
                            dist = -0.5 * math.log(w1) - 0.25 * math.log(w2)
                    except:
                         dist = 1.0
                else:
                    # p-distance
                    dist = p_dist
                
            dist_matrix[i, j] = dist_matrix[j, i] = dist
            
    # UPGMA
    Z = linkage(squareform(dist_matrix), method='average')
    
    return ids, dist_matrix, Z, None

def run_asap_scan(aligned_fasta_path, model='p-dist', start=0.00, end=0.15, step=0.005):
    """
    閾値をスキャンして、それぞれの閾値での種数（クラスター数）を算出する
    """
    ids, dist_matrix, Z, err = calculate_distance_matrix(aligned_fasta_path, model=model)
    if err: return None, None, None, err
    
    scan_results = []
    # 閾値スキャン
    thresholds = np.arange(start, end + step, step)
    
    for t in thresholds:
        t = round(t, 4)
        clusters = fcluster(Z, t=t, criterion='distance')
        num_clusters = len(np.unique(clusters))
        scan_results.append({
            "Threshold (Distance)": t,
            "OTU Count": num_clusters
        })
        
    df_scan = pd.DataFrame(scan_results)
    
    # ASAP Score計算（簡易版: 閾値の幅が広いパーティションを優先）
    # ここでは単純に表示するのみとし、ユーザーに選択させる
    
    return df_scan, dist_matrix, Z, ids

def get_partition_by_threshold(Z, ids, threshold):
    """
    計算済みのZとIDを使って、特定の閾値でのパーティションを取得する
    """
    clusters = fcluster(Z, t=threshold, criterion='distance')
    result_df = pd.DataFrame({"ID": ids, "Cluster": clusters})
    result_df = result_df.sort_values("Cluster")
    return result_df
