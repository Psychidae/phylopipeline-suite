import pandas as pd
import numpy as np
import math
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor, DistanceMatrix
from Bio.Phylo.Consensus import bootstrap_trees, get_support
from Bio import SeqIO, Phylo

class CustomDistanceCalculator:
    """
    Custom Distance Calculator to handle models missing in some BioPython environments (e.g. kimura80).
    Mimics Bio.Phylo.TreeConstruction.DistanceCalculator interface.
    """
    def __init__(self, model='identity'):
        self.model = model

    def get_distance(self, msa):
        """
        Calculate distance matrix for MSA.
        Returns Bio.Phylo.TreeConstruction.DistanceMatrix
        """
        names = [s.id for s in msa]
        n = len(msa)
        
        # Extract sequences as list of strings
        seqs = [str(s.seq).upper() for s in msa]
        
        # Calculate numpy matrix
        dist_mat_np = _calc_dist_matrix_numpy(seqs, self.model)
        
        # Convert to list of lists (lower triangular)
        # BioPython DistanceMatrix expects: [ [], [d21], [d31, d32], ... ]
        matrix_list = []
        for i in range(n):
            row = []
            for j in range(i):
                row.append(dist_mat_np[i, j])
            matrix_list.append(row)
            
        return DistanceMatrix(names, matrix_list)

def _calc_dist_matrix_numpy(seqs, model):
    """
    Core numpy calculation logic
    """
    max_len = max(len(s) for s in seqs)
    # Pad sequences
    matrix = []
    for s in seqs:
        s_pad = s.ljust(max_len, '-')
        matrix.append(np.array(list(s_pad)))
    matrix = np.array(matrix)
    
    n = len(seqs)
    dist_matrix = np.zeros((n, n))
    
    # Check if protein model requested in DNA context, fallback to p-dist?
    # For now, only implement k2p and p-dist explicitly.
    
    # 定義: Transition (A<->G, C<->T)
    def is_transition_vec(c1, c2):
         # Vectorized transition check is hard with chars.
         # Use simplified check or loops.
         # Since this is "Custom", we assume user wants the logic we have.
         pass
         
    # ... (Reusing logic from original calculate_distance_matrix) ...
    # Optimization: Use indices for calculation
    
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
                
                if model == 'kimura80' or model == 'k2p':
                    # K2P Logic
                    pairs = np.stack((s1[valid_mask], s2[valid_mask]), axis=1)
                    
                    # Transition: {A, G} or {C, T}
                    # A=65, C=67, G=71, T=84
                    # A+G = 136, C+T = 151. No overlap.
                    # Or just:
                    # ts = ((p0=='A')&(p1=='G')) | ((p0=='G')&(p1=='A')) | ((p0=='C')&(p1=='T')) | ((p0=='T')&(p1=='C'))
                    
                    p0 = pairs[:, 0]
                    p1 = pairs[:, 1]
                    
                    mask_diff = (p0 != p1)
                    if np.sum(mask_diff) == 0:
                        ts_count = 0
                        tv_count = 0
                    else:
                        pd0 = p0[mask_diff]
                        pd1 = p1[mask_diff]
                        
                        is_ts = (
                            ((pd0 == 'A') & (pd1 == 'G')) | 
                            ((pd0 == 'G') & (pd1 == 'A')) |
                            ((pd0 == 'C') & (pd1 == 'T')) | 
                            ((pd0 == 'T') & (pd1 == 'C'))
                        )
                        ts_count = np.sum(is_ts)
                        tv_count = np.sum(~is_ts)
                    
                    P = ts_count / valid_count
                    Q = tv_count / valid_count
                    
                    try:
                        w1 = 1.0 - 2.0*P - Q
                        w2 = 1.0 - 2.0*Q
                        if w1 <= 0 or w2 <= 0:
                            dist = 10.0
                        else:
                            dist = -0.5 * np.log(w1) - 0.25 * np.log(w2)
                    except:
                        dist = 1.0
                else:
                    # Default to p-distance for everything else
                    dist = p_dist
            
            dist_matrix[i, j] = dist_matrix[j, i] = dist
            
    return dist_matrix

def run_phylo_bootstrap(msa, method="nj", model="identity", replicates=100):
    """
    Run NJ or UPGMA with bootstrap support.
    """
    # Use CustomCalculator for 'kimura80' to avoid platform issues.
    # For standard 'identity', we could use BioPython's, but consistently using ours is safe too.
    # Let's use ours if model is k2p, else try standard? 
    # Actually, CustomCalculator is more robust given the error seen.
    if model in ['kimura80', 'k2p', 'identity']:
         calculator = CustomDistanceCalculator(model)
    else:
         # Fallback to BioPython standard (if available)
         calculator = DistanceCalculator(model)

    constructor = DistanceTreeConstructor(calculator, method)
    
    # 2. Build Main Tree
    main_tree = constructor.build_tree(msa)
    
    # 3. Bootstrap
    if replicates > 0:
        # boostrap_trees returns a generator. get_support needs a list or len_trees.
        boot_trees = list(bootstrap_trees(msa, replicates, constructor))
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
    Compatible with ASAP scan
    """
    seqs = list(SeqIO.parse(aligned_fasta_path, "fasta"))
    if len(seqs) < 3: return None, None, None, "配列数が少なすぎます（最低3配列必要）"
    
    ids = [s.id for s in seqs]
    seq_strs = [str(s.seq).upper() for s in seqs]
    
    # Use shared logic
    # map model names: 'p-dist' -> 'identity' for consistency in core logic?
    # Our core logic supports 'kimura80'/'k2p' and defaults to identity.
    # ASAP calls pass 'p-dist' or 'k2p'.
    core_model = 'kimura80' if model == 'k2p' else 'identity'
    
    dist_matrix = _calc_dist_matrix_numpy(seq_strs, core_model)
    
    # UPGMA for ASAP (scipy linkage)
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
