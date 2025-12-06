import io
import numpy as np
import re
from collections import Counter
from Bio import SeqIO, Align
from Bio.Seq import Seq

# --- 基本解析 ---
def parse_ab1(uploaded_file):
    try:
        uploaded_file.seek(0)
        bytes_io = io.BytesIO(uploaded_file.read())
        record = SeqIO.read(bytes_io, "abi")
        raw = record.annotations['abif_raw']
        trace_data = {
            'G': list(raw.get('DATA9', [])),
            'A': list(raw.get('DATA10', [])),
            'T': list(raw.get('DATA11', [])),
            'C': list(raw.get('DATA12', []))
        }
        ploc = list(raw.get('PLOC1', []))
        if not ploc: ploc = list(raw.get('PLOC2', []))
        quals = record.letter_annotations.get("phred_quality", [0]*len(record.seq))
        return {"record": record, "sequence": list(str(record.seq)), "qualities": quals, "trace": trace_data, "peak_locations": ploc, "file_name": uploaded_file.name, "original_length": len(record.seq)}
    except Exception as e:
        return {"error": str(e), "file_name": uploaded_file.name}

def trim_sequence(data, quality_threshold=20, window_size=10):
    quals = data["qualities"]
    seq_len = len(quals)
    if seq_len == 0: return data, 0, 0
    start, end = 0, seq_len
    for i in range(seq_len - window_size):
        if np.mean(quals[i : i + window_size]) >= quality_threshold: start = i; break
    for i in range(seq_len, window_size, -1):
        if np.mean(quals[i - window_size : i]) >= quality_threshold: end = i; break
    if start >= end: start, end = 0, seq_len 
    new_data = data.copy()
    new_data["sequence"] = data["sequence"][start:end]
    new_data["qualities"] = data["qualities"][start:end]
    new_data["peak_locations"] = data["peak_locations"][start:end]
    return new_data, start, (seq_len - end)

def get_rev_comp(data):
    seq_obj = Seq("".join(data["sequence"]))
    rc_seq = list(str(seq_obj.reverse_complement()))
    rc_qual = data["qualities"][::-1]
    trace_len = len(data["trace"]['A'])
    rc_trace = {'G': data["trace"]['C'][::-1], 'C': data["trace"]['G'][::-1], 'A': data["trace"]['T'][::-1], 'T': data["trace"]['A'][::-1]}
    rc_peaks = [trace_len - 1 - p for p in data["peak_locations"][::-1]]
    return {"sequence": rc_seq, "qualities": rc_qual, "trace": rc_trace, "peak_locations": rc_peaks, "file_name": data["file_name"] + " (Rev)", "is_rc": True}

def align_pair(ref_data, query_data):
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 2; aligner.mismatch_score = -3; aligner.open_gap_score = -5; aligner.extend_gap_score = -2
    seq_ref = "".join(ref_data["sequence"])
    seq_query = "".join(query_data["sequence"])
    score_fwd = aligner.score(seq_ref, seq_query)
    seq_query_rc = str(Seq(seq_query).reverse_complement())
    score_rev = aligner.score(seq_ref, seq_query_rc)
    is_rc = score_rev > score_fwd
    final_query_seq = seq_query_rc if is_rc else seq_query
    alignments = aligner.align(seq_ref, final_query_seq)
    if not alignments: return False, 0
    best = alignments[0]
    try:
        t_start = best.coordinates[0][0]
        q_start = best.coordinates[1][0]
        offset = int(t_start - q_start)
    except: offset = 0
    return is_rc, offset

def calculate_initial_consensus(results, ref_length=None):
    if not results and ref_length is None: return []
    ref_len = len(results[0]["display"]["sequence"]) if ref_length is None else ref_length
    consensus = []
    for i in range(ref_len):
        bases = []
        for track in results:
            loc = i - track["offset"]
            if 0 <= loc < len(track["display"]["sequence"]):
                b = track["display"]["sequence"][loc]
                if b not in ['N', '-']: bases.append(b)
        if bases:
            count = Counter(bases)
            cons = count.most_common(1)[0][0]
            consensus.append({"base": cons})
        else:
            consensus.append({"base": 'N'})
    return consensus

def find_orfs(seq_str, table_id=1, min_len_aa=10):
    seq_obj = Seq(seq_str.replace("-", "").upper()) 
    orfs = []
    for frame in range(3):
        try:
            trans = str(seq_obj[frame:].translate(table=table_id))
            for m in re.finditer(r'M[^*]*\*', trans):
                if len(m.group()) - 1 >= min_len_aa:
                    s = frame + m.start()*3; e = frame + m.end()*3 
                    orfs.append({'frame': frame+1, 'start': s, 'end': e, 'len': len(m.group())-1})
        except: pass
    rc_seq = seq_obj.reverse_complement()
    slen = len(seq_obj)
    for frame in range(3):
        try:
            trans = str(rc_seq[frame:].translate(table=table_id))
            for m in re.finditer(r'M[^*]*\*', trans):
                if len(m.group()) - 1 >= min_len_aa:
                    rc_s = frame + m.start()*3; rc_e = frame + m.end()*3
                    orfs.append({'frame': -(frame+1), 'start': slen-rc_e, 'end': slen-rc_s, 'len': len(m.group())-1})
        except: pass
    return orfs

def calculate_translation_map(consensus_data, table_id=1, frame=1):
    map_data = []
    valid_indices = [i for i, d in enumerate(consensus_data) if d["base"] not in ['-']]
    start_offset = frame - 1
    if len(valid_indices) <= start_offset: return []
    valid_seq_indices = valid_indices[start_offset:]
    for i in range(0, len(valid_seq_indices), 3):
        chunk = valid_seq_indices[i : i+3]
        if len(chunk) < 3: break
        codon_str = "".join([consensus_data[idx]["base"] for idx in chunk]).upper()
        try: aa = str(Seq(codon_str).translate(table=table_id))
        except: aa = "?"
        center_seq_idx = chunk[1]
        map_data.append({"aa": aa, "seq_idx": center_seq_idx, "is_stop": aa == "*"})
    return map_data

def get_linearized_trace(data, start_base, end_base):
    # (使用しなくなったため省略可だが、念のため残すなら古い実装)
    return [], {}

def process_contig_group(target_files, auto_trim, quality_thresh):
    raw_list = []
    for f in target_files:
        d = parse_ab1(f)
        if "error" not in d:
            if auto_trim: d, s, e = trim_sequence(d, quality_thresh)
            raw_list.append(d)
    if not raw_list: return None
    raw_list.sort(key=lambda x: len(x["sequence"]), reverse=True)
    ref = raw_list[0]
    res = [{"name": ref["file_name"], "data": ref, "display": ref, "offset": 0, "is_rc": False}]
    for query in raw_list[1:]:
        is_rc, offset = align_pair(ref, query)
        disp = get_rev_comp(query) if is_rc else query
        res.append({"name": query["file_name"], "data": query, "display": disp, "offset": offset, "is_rc": is_rc})
    cons = calculate_initial_consensus(res)
    return {"results": res, "consensus": cons}