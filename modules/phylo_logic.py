import pandas as pd
import numpy as np
import math
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor, DistanceMatrix
from Bio.Phylo.Consensus import bootstrap_trees, get_support
from Bio import SeqIO, Phylo
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform


class CustomDistanceCalculator(DistanceCalculator):
    """
    kimura80 など BioPython 環境によっては存在しないモデルを扱うカスタム計算機。
    Bio.Phylo.TreeConstruction.DistanceCalculator インターフェースに準拠。
    """
    def __init__(self, model='identity'):
        super().__init__('identity')
        self.model = model

    def get_distance(self, msa):
        names = [s.id for s in msa]
        n = len(msa)
        seqs = [str(s.seq).upper() for s in msa]
        dist_mat_np = _calc_dist_matrix_numpy(seqs, self.model)

        # BioPython DistanceMatrix: row i には i+1 要素（対角含む）
        matrix_list = []
        for i in range(n):
            row = [dist_mat_np[i, j] for j in range(i + 1)]
            matrix_list.append(row)

        return DistanceMatrix(names, matrix_list)


def _calc_dist_matrix_numpy(seqs, model):
    """numpy による距離行列計算（p-dist / Kimura 2-Parameter 対応）"""
    max_len = max(len(s) for s in seqs)
    matrix = np.array([list(s.ljust(max_len, '-')) for s in seqs])
    n = len(seqs)
    dist_matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(i + 1, n):
            s1, s2 = matrix[i], matrix[j]
            valid_mask = (s1 != '-') & (s2 != '-') & (s1 != 'N') & (s2 != 'N')
            valid_count = np.sum(valid_mask)

            if valid_count == 0:
                dist = 0.0
            else:
                p_dist = np.sum(s1[valid_mask] != s2[valid_mask]) / valid_count

                if model in ('kimura80', 'k2p'):
                    pairs_0 = s1[valid_mask]
                    pairs_1 = s2[valid_mask]
                    mask_diff = pairs_0 != pairs_1

                    if np.sum(mask_diff) == 0:
                        P = Q = 0.0
                    else:
                        pd0, pd1 = pairs_0[mask_diff], pairs_1[mask_diff]
                        is_ts = (
                            ((pd0 == 'A') & (pd1 == 'G')) |
                            ((pd0 == 'G') & (pd1 == 'A')) |
                            ((pd0 == 'C') & (pd1 == 'T')) |
                            ((pd0 == 'T') & (pd1 == 'C'))
                        )
                        P = np.sum(is_ts) / valid_count
                        Q = np.sum(~is_ts) / valid_count

                    w1 = 1.0 - 2.0 * P - Q
                    w2 = 1.0 - 2.0 * Q
                    if w1 <= 0 or w2 <= 0:
                        dist = 10.0  # 飽和距離
                    else:
                        dist = -0.5 * np.log(w1) - 0.25 * np.log(w2)
                else:
                    dist = p_dist

            dist_matrix[i, j] = dist_matrix[j, i] = dist

    return dist_matrix


def run_phylo_bootstrap(msa, method="nj", model="identity", replicates=100):
    """NJ / UPGMA + ブートストラップを実行してツリーを返す"""
    if model in ('kimura80', 'k2p', 'identity'):
        calculator = CustomDistanceCalculator(model)
    else:
        calculator = DistanceCalculator(model)

    constructor = DistanceTreeConstructor(calculator, method)
    main_tree = constructor.build_tree(msa)

    if replicates > 0:
        boot_trees = list(bootstrap_trees(msa, replicates, tree_constructor=constructor))
        consensus_tree = get_support(main_tree, boot_trees)

        # ブートストラップ値をノード名に変換（FigTree 等のビューア向け）
        for clade in consensus_tree.find_clades():
            if clade.confidence is not None:
                clade.name = (
                    str(int(clade.confidence))
                    if isinstance(clade.confidence, float) and clade.confidence.is_integer()
                    else str(clade.confidence)
                )
                clade.confidence = None

        return consensus_tree
    else:
        return main_tree


def generate_methods_log(tool_versions, params):
    """解析手法ログを生成する"""
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
    if params.get("trimAl") == "Yes":
        log.append("- trimAl: Capella-Gutierrez et al. (2009) Bioinformatics 25:1972-1973")
    tool_name = params.get("Tool", "")
    if "IQ-TREE" in tool_name:
        log.append("- IQ-TREE: Minh et al. (2020) MBE 37:1530-1534")
        log.append("- ModelFinder: Kalyaanamoorthy et al. (2017) Nat Methods 14:587-589")
        log.append("- UFBoot: Hoang et al. (2018) MBE 35:518-522")
    elif "NJ" in tool_name or "UPGMA" in tool_name:
        log.append("- BioPython: Cock et al. (2009) Bioinformatics 25:1422-1423")
    log.append("")
    return "\n".join(log)


def calculate_distance_matrix(aligned_fasta_path, model='p-dist'):
    """
    アラインメント済みFASTAから距離行列・linkage行列・IDリストを返す。
    model: 'p-dist' or 'k2p'
    エラー時: (None, None, None, error_message)
    """
    seqs = list(SeqIO.parse(aligned_fasta_path, "fasta"))
    if len(seqs) < 3:
        return None, None, None, "配列数が少なすぎます（最低3配列必要）"

    ids = [s.id for s in seqs]
    seq_strs = [str(s.seq).upper() for s in seqs]

    core_model = 'kimura80' if model == 'k2p' else 'identity'
    dist_matrix = _calc_dist_matrix_numpy(seq_strs, core_model)

    condensed = squareform(dist_matrix)
    Z = linkage(condensed, method='average')

    return ids, dist_matrix, Z, None


def run_asap_scan(aligned_fasta_path, model='p-dist', start=0.00, end=0.15, step=0.005):
    """
    閾値をスキャンして各閾値でのクラスター数（OTU数）を算出する。
    戻り値: (df_scan, dist_matrix, Z, ids)  エラー時は (None, None, None, error_message)
    """
    ids, dist_matrix, Z, err = calculate_distance_matrix(aligned_fasta_path, model=model)
    if err:
        return None, None, None, err

    scan_results = []
    for t in np.arange(start, end + step, step):
        t = round(t, 4)
        clusters = fcluster(Z, t=t, criterion='distance')
        scan_results.append({
            "Threshold (Distance)": t,
            "OTU Count": len(np.unique(clusters))
        })

    df_scan = pd.DataFrame(scan_results)
    return df_scan, dist_matrix, Z, ids


def get_partition_by_threshold(Z, ids, threshold):
    """計算済みの Z と IDs を使い、指定閾値でのパーティションを返す"""
    clusters = fcluster(Z, t=threshold, criterion='distance')
    result_df = pd.DataFrame({"ID": ids, "Cluster": clusters})
    result_df = result_df.sort_values("Cluster").reset_index(drop=True)
    return result_df
