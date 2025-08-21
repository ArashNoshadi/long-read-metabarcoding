# merged_analysis_optimized_v2.py
import os, re, time, math, logging, threading, concurrent.futures, sys, subprocess, datetime
from collections import defaultdict, Counter
import psutil, pandas as pd, numpy as np, matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from colorama import Fore, Style, init
from numba import njit, prange
try:
    import pyopencl as cl
    gpu_available = True
except Exception:
    gpu_available = False

# ---------- CONFIG ---------------------------------------------------------
BAD_VALUES = {"unknow", "unknown", "nan", "unidentified", "nematode", "environmental", "sample",""}
def format_species_name(genus, species):
    return f"{genus} sp." if not species or species.strip().lower() in BAD_VALUES else f"{genus} {species}"

DATA_FILES = {
    '18S': r'G:\Paper\nema-Nanopore-Sequencing\by all gene\results\18s.xlsx',
    '28S': r'G:\Paper\nema-Nanopore-Sequencing\by all gene\results\28s.xlsx',
    'cox1': r'G:\Paper\nema-Nanopore-Sequencing\by all gene\results\cox1.xlsx',
}
EXPECTED_LENGTHS = {'18S': 1700, '28S': 3500, 'cox1': 700}
MUSCLE_PATH = r'G:\Paper\nema-Nanopore-Sequencing\by all gene\results\musclev5.exe'
OUTPUT_DIR = r"G:\Paper\nema-Nanopore-Sequencing\by all gene\combined_results-5"
os.makedirs(OUTPUT_DIR, exist_ok=True)
SKIP_ALIGNMENT = True # اگر می‌خواهید MUSCLE اجرا شود، این را False کنید

# ---------- LOGGING --------------------------------------------------------
logger = logging.getLogger()
logger.setLevel(logging.INFO)
stream_handler = logging.StreamHandler(sys.stdout)
stream_handler.setLevel(logging.INFO)
fmt = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
stream_handler.setFormatter(fmt)
logger.handlers = [] # remove default handlers
logger.addHandler(stream_handler)
file_handler = logging.FileHandler(os.path.join(OUTPUT_DIR, "combined_results.log"), mode='w', encoding='utf-8')
file_handler.setLevel(logging.INFO)
file_handler.setFormatter(fmt)
logger.addHandler(file_handler)
init(autoreset=True)

# ---------- GLOBALS --------------------------------------------------------
alignment_summary = {marker: [] for marker in DATA_FILES.keys()}
error_families = []
error_families_lock = threading.Lock()
global_reference_info = {}
GPU_CONTEXT = None
GPU_DEVICE = None

# ---------- HELPERS: IO & CLEANING ---------------------------------------
def normalize_taxa_in_df(df, columns):
    for col in columns:
        if col in df.columns:
            df[col] = df[col].astype(str).str.strip().str.lower()
    return df

def can_submit_task(cpu_threshold=80, mem_threshold=30):
    try:
        cpu_usage = psutil.cpu_percent(interval=0.5)
        available = (psutil.virtual_memory().available/psutil.virtual_memory().total)*100
        return cpu_usage < cpu_threshold and available > mem_threshold
    except Exception:
        return True

def read_clean_data(file_path, required_cols):
    df = pd.read_excel(file_path)
    for col in required_cols:
        if col in df.columns:
            df[col] = df[col].astype(str).str.strip()
    if "ID" in df.columns:
        df["ID"] = df["ID"].astype(str).str.strip()
    return df.dropna(subset=required_cols)

def vectorized_clean_sequences(df, seq_col='RNA Sequence'):
    if seq_col in df.columns:
        df[seq_col] = df[seq_col].astype(str).str.upper().str.replace(r'[^ATCG]', '', regex=True)
    return df

def filter_valid_taxa(df: pd.DataFrame, column: str) -> pd.DataFrame:
    if column not in df.columns:
        return df
    # OPTIMIZATION: Use .loc for safer filtering and assignment
    return df.loc[df[column].notna() & (~df[column].str.strip().str.lower().isin(BAD_VALUES))].copy()

def load_data():
    logger.info("Loading and cleaning data from Excel files...")
    data_align = {}
    for marker, path in DATA_FILES.items():
        try:
            df = read_clean_data(path, ['Family','Main_Organism','Species','RNA Sequence','cp','f-h'])
            df = normalize_taxa_in_df(df, ['Family','Main_Organism','Species'])
            df = vectorized_clean_sequences(df, 'RNA Sequence')
            df = df[df['RNA Sequence'].astype(bool)]
            
            if "Length (RNA)" in df.columns:
                df["Length (RNA)"] = pd.to_numeric(df["Length (RNA)"], errors='coerce')
                expected = EXPECTED_LENGTHS.get(marker)
                if expected:
                    # <<<--- START OF CHANGE --->>>
                    # NEW: Create multiple length check columns based on specified percentages
                    percentages = [60, 70, 80, 90, 100]
                    for p in percentages:
                        col_name = f'Length_Check_{p}%'
                        threshold_value = (p / 100.0) * expected
                        df[col_name] = ((df["Length (RNA)"] >= threshold_value) & df["Length (RNA)"].notna()).astype(int)
                    # <<<--- END OF CHANGE --->>>

            data_align[marker] = df
            logger.info(f"{marker}: {len(df)} sequences loaded.")
        except Exception as e:
            logger.error(f"Error reading {path}: {e}")
    logger.info("Data loading complete.")
    return data_align

def save_updated_files(data_align):
    for marker, df in data_align.items():
        output_path = os.path.join(OUTPUT_DIR, f"{marker}_updated.xlsx")
        try:
            if 'cp_f-h' in df.columns:
                df.drop(columns=['cp_f-h'], inplace=True)
            df.to_excel(output_path, index=False)
            logger.info(f"Updated file saved: {output_path}")
        except Exception as e:
            logger.error(f"Could not save updated file for {marker}: {e}")

def prepare_coverage_data(data_align):
    return {marker: df[['Family','Main_Organism','Species']].copy() for marker, df in data_align.items()}

def build_reference_info(data_align):
    """ OPTIMIZATION: Replaced slow iterrows with efficient pandas methods. """
    ref_info = {}
    for marker, df in data_align.items():
        # Filter for rows with valid IDs
        df_filtered = df.dropna(subset=['ID'])
        df_filtered = df_filtered[df_filtered['ID'].astype(str).str.strip() != ""]
        # Use set_index and to_dict for massive speedup
        ref_info[marker] = df_filtered[['ID', 'Main_Organism', 'Species']].set_index('ID').to_dict('index')
    return ref_info

def get_header_info(rec, marker, ref_info):
    info = {}
    for part in rec.description.split(';'):
        if '=' in part:
            key, value = part.split('=',1)
            info[key.strip()] = value.strip().lower()
    # OPTIMIZATION: Use .get with a default empty dict to avoid multiple lookups
    ref_marker_info = ref_info.get(marker, {})
    rec_id_info = ref_marker_info.get(rec.id.strip(), {})
    info["Main_Organism"] = info.get("Main_Organism", rec_id_info.get("Main_Organism", "unknow"))
    info["Species"] = info.get("Species", rec_id_info.get("Species", "unknow"))
    return info

# ---------- GPU selection (optional) -------------------------------------
def select_best_gpu_opencl():
    try:
        platforms = cl.get_platforms()
        best_device = None; max_compute_units = 0
        for platform in platforms:
            for device in platform.get_devices(device_type=cl.device_type.GPU):
                if device.max_compute_units > max_compute_units:
                    max_compute_units = device.max_compute_units; best_device = device
        if best_device:
            context = cl.Context(devices=[best_device])
            logger.info(f"Selected GPU: {best_device.name} with {best_device.max_compute_units} compute units.")
            return context, best_device
    except Exception as e:
        logger.error(f"Error selecting GPU: {e}")
    return None, None

def init_gpu_optional():
    global GPU_CONTEXT, GPU_DEVICE
    if gpu_available:
        GPU_CONTEXT, GPU_DEVICE = select_best_gpu_opencl()
        if GPU_CONTEXT is not None:
            logger.info(f"Using GPU: {GPU_DEVICE.name}")
        else:
            logger.info("GPU available but not selected/usable; continuing on CPU.")
    else:
        logger.info("GPU acceleration not available; continuing on CPU.")

# ---------- ALIGNMENT METRICS (NUMBA/OPENCL) -----------------------------
@njit(parallel=True, cache=True)
def compute_alignment_metrics_numba(seq_mat):
    n, L = seq_mat.shape
    total_similarity = 0.0
    total_overlap = 0.0
    total_changes = 0.0
    pair_count = 0
    # The number of pairs is n*(n-1)/2
    for i in prange(n):
        for j in range(i + 1, n):
            overlap = 0
            identical = 0
            for k in range(L):
                a = seq_mat[i, k]
                b = seq_mat[j, k]
                if a != 45 and b != 45: # 45 is ASCII for '-'
                    overlap += 1
                    if a == b:
                        identical += 1
            if overlap > 0:
                total_similarity += (identical / overlap) * 100.0
                total_overlap += (overlap / L) * 100.0
                total_changes += ((L - identical) / L) * 100.0
                pair_count += 1

    if pair_count == 0:
        return 0.0, 0.0, 0.0
    return total_similarity / pair_count, total_overlap / pair_count, total_changes / pair_count

def compute_alignment_metrics_cpu(records):
    seqs = [str(rec.seq) for rec in records]
    if not seqs or not seqs[0]:
        return 0.0, 0.0, 0.0
    n, L = len(seqs), len(seqs[0])
    seq_mat = np.array([list(s) for s in seqs], dtype='S1').view(np.uint8)
    return compute_alignment_metrics_numba(seq_mat)

# OpenCL function remains unchanged, assuming it's already optimized.
def compute_alignment_metrics_opencl(records):
    seqs = [str(rec.seq) for rec in records]
    n = len(seqs); L = len(seqs[0])
    seq_mat = np.empty((n, L), dtype=np.uint8)
    for i, seq in enumerate(seqs):
        seq_mat[i, :] = np.frombuffer(bytes(seq, 'utf-8'), dtype=np.uint8)[:L]
    totalPairs = n*(n-1)//2
    if totalPairs == 0: return 0.0, 0.0, 0.0
    similarity_array = np.empty(totalPairs, dtype=np.float64)
    overlap_array = np.empty(totalPairs, dtype=np.float64)
    changes_array = np.empty(totalPairs, dtype=np.float64)
    if GPU_CONTEXT is None:
        raise Exception("GPU not initialized.")
    queue = cl.CommandQueue(GPU_CONTEXT)
    mf = cl.mem_flags
    d_seq_mat = cl.Buffer(GPU_CONTEXT, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=seq_mat)
    d_similarity = cl.Buffer(GPU_CONTEXT, mf.WRITE_ONLY, similarity_array.nbytes)
    d_overlap = cl.Buffer(GPU_CONTEXT, mf.WRITE_ONLY, overlap_array.nbytes)
    d_changes = cl.Buffer(GPU_CONTEXT, mf.WRITE_ONLY, changes_array.nbytes)
    kernel_code = """
    __kernel void pair_metrics(__global const uchar *seq_mat, int L, int n, __global double *similarity_array, __global double *overlap_array, __global double *changes_array) {
        int t = get_global_id(0);
        int totalPairs = n*(n-1)/2;
        if(t >= totalPairs) return;
        int i = 0, s = n - 1, temp = t;
        while(temp >= s) { temp -= s; i++; s = n - i - 1; }
        int j = i + 1 + temp;
        int overlap = 0, identical = 0;
        for (int k = 0; k < L; k++){
            uchar a = seq_mat[i*L+k], b = seq_mat[j*L+k];
            if(a != 45 && b != 45) { overlap++; if(a == b) identical++; }
        }
        double sim = 0.0, ov_pct = 0.0, ch_pct = 0.0;
        if(overlap > 0) { sim = ((double)identical/overlap)*100.0; ov_pct = ((double)overlap/L)*100.0; ch_pct = (((double)L-identical)/L)*100.0; }
        similarity_array[t] = sim; overlap_array[t] = ov_pct; changes_array[t] = ch_pct;
    }"""
    program = cl.Program(GPU_CONTEXT, kernel_code).build()
    global_size = (totalPairs,)
    program.pair_metrics(queue, global_size, None, d_seq_mat, np.int32(L), np.int32(n), d_similarity, d_overlap, d_changes)
    cl.enqueue_copy(queue, similarity_array, d_similarity).wait()
    cl.enqueue_copy(queue, overlap_array, d_overlap).wait()
    cl.enqueue_copy(queue, changes_array, d_changes).wait()
    return similarity_array.mean(), overlap_array.mean(), changes_array.mean()

def compute_alignment_metrics(records):
    if len(records) < 2:
        return 0.0, 0.0, 0.0
    if gpu_available and GPU_CONTEXT is not None:
        try:
            return compute_alignment_metrics_opencl(records)
        except Exception as e:
            logger.warning(f"GPU compute failed, falling back to CPU: {e}")
            return compute_alignment_metrics_cpu(records)
    else:
        return compute_alignment_metrics_cpu(records)

# OPTIMIZATION: Consolidate similarity/overlap calculation to use the fast Numba/GPU versions
def compute_alignment_similarity(records: list) -> float:
    sim, _, _ = compute_alignment_metrics(records)
    return sim

def compute_alignment_overlap(records: list) -> float:
    _, ov, _ = compute_alignment_metrics(records)
    return ov

# ---------- CONSERVATION / DIVERSITY / TAJIMA ETC. (NUMBA OPTIMIZED) -------
@njit(cache=True)
def _numba_nucleotide_diversity(seq_arr):
    n, L = seq_arr.shape
    total_diff = 0
    comparisons = n * (n - 1) // 2
    if comparisons == 0:
        return 0.0
    for i in range(n):
        for j in range(i + 1, n):
            diff = 0
            for k in range(L):
                if seq_arr[i, k] != seq_arr[j, k]:
                    diff += 1
            total_diff += diff
    return total_diff / (comparisons * L)

def compute_nucleotide_diversity(records):
    if len(records) < 2: return None
    seqs = [str(rec.seq) for rec in records]
    arr = np.array([list(s) for s in seqs])
    return _numba_nucleotide_diversity(arr)

@njit(cache=True)
def _numba_segregating_sites(seq_arr):
    L = seq_arr.shape[1]
    S = 0
    for i in range(L):
        col = seq_arr[:, i]
        # Find unique non-gap characters
        unique_nts = set()
        for char_code in col:
            if char_code != 45: # '-'
                unique_nts.add(char_code)
        if len(unique_nts) > 1:
            S += 1
    return S

def compute_tajimas_d(records):
    if len(records) < 2: return None
    n = len(records)
    seqs = [str(rec.seq) for rec in records]
    seq_mat = np.array([list(s) for s in seqs], dtype='S1').view(np.uint8)
    L = seq_mat.shape[1]

    S = _numba_segregating_sites(seq_mat)
    if S == 0: return 0.0

    pi = compute_nucleotide_diversity(records) * L # Pi is total differences, not per site
    if pi is None: return None

    a1 = sum(1.0/i for i in range(1, n))
    theta_w = S / a1

    a2 = sum(1.0/(i**2) for i in range(1, n))
    b1 = (n+1) / (3*(n-1))
    b2 = 2*(n**2+n+3) / (9*n*(n-1))
    c1 = b1 - 1/a1
    c2 = b2 - (n+2)/(a1*n) + a2/(a1**2)
    e1 = c1 / a1
    e2 = c2 / (a1**2 + a2)
    variance = e1 * S + e2 * S * (S-1)
    
    return (pi - theta_w) / math.sqrt(variance) if variance > 0 else 0.0

@njit(cache=True)
def _numba_analyze_conservation(seq_arr):
    n, L = seq_arr.shape
    scores = np.zeros(L, dtype=np.float64)
    for i in range(L):
        col = seq_arr[:, i]
        # Filter out gaps ('-')
        valid_count = 0
        for char_code in col:
            if char_code != 45:
                valid_count += 1
        
        if valid_count == 0:
            scores[i] = 0
        else:
            # Count occurrences of each base
            counts = np.zeros(256, dtype=np.int32) # Full ASCII range for simplicity
            for char_code in col:
                if char_code != 45:
                    counts[char_code] += 1
            
            max_count = 0
            for count in counts:
                if count > max_count:
                    max_count = count
            scores[i] = max_count / valid_count
    return scores

def analyze_conservation(records, threshold=0.9):
    if not records: return None, None, None, None
    seqs = [str(rec.seq) for rec in records]
    arr = np.array([list(s) for s in seqs], dtype='S1').view(np.uint8)
    scores = _numba_analyze_conservation(arr)
    conserved = np.sum(scores >= threshold)
    return conserved, len(scores) - conserved, np.mean(scores), scores

# Other functions like compute_nucleotide_change_distribution are generally fast enough
def compute_nucleotide_change_distribution(records):
    if len(records) < 2: return None
    arr = np.array([list(str(rec.seq)) for rec in records])
    n, L = arr.shape
    nucleotide_counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
    total_valid_nts = 0
    valid_nts = set(nucleotide_counts.keys())
    
    for j in range(L):
        col = arr[:, j]
        unique, counts = np.unique(col, return_counts=True)
        # Check if there's variation at this site among valid nucleotides
        if len(set(unique) & valid_nts) > 1:
            for nt, count in zip(unique, counts):
                if nt in nucleotide_counts:
                    nucleotide_counts[nt] += count
                    total_valid_nts += count

    return {nt: (count / total_valid_nts) * 100 for nt, count in nucleotide_counts.items()} if total_valid_nts > 0 else {nt: 0 for nt in nucleotide_counts}

def plot_conservation(cons_scores, marker, output_folder):
    plt.figure(figsize=(12, 4))
    plt.plot(cons_scores)
    plt.xlabel("Position")
    plt.ylabel("Conservation Score")
    plt.title(f"Conservation Profile for {marker}")
    out_file = os.path.join(output_folder, f"Conservation_{marker}.png")
    plt.savefig(out_file, dpi=300)
    plt.close()
    logger.info(f"Conservation plot saved: {out_file}")

# ---------- FASTA & MUSCLE -----------------------------------------------
def generate_fasta_file(records, fasta_path):
    if os.path.exists(fasta_path) and os.path.getsize(fasta_path) == 0:
        os.remove(fasta_path)
        logger.info(f"Empty FASTA file removed: {fasta_path}")
    if not os.path.exists(fasta_path):
        try:
            SeqIO.write(records, fasta_path, "fasta")
            logger.info(f"FASTA file created: {fasta_path}")
        except Exception as e:
            logger.error(f"Error writing FASTA file {fasta_path}: {e}")
            raise
    else:
        logger.info(f"FASTA file already exists and is not empty: {fasta_path}.")

def perform_muscle_alignment(fasta_input, fasta_output, use_super):
    if SKIP_ALIGNMENT:
        logger.info(f"Skipping MUSCLE alignment (SKIP_ALIGNMENT=True) for {fasta_output}")
        return
    if os.path.exists(fasta_output) and os.path.getsize(fasta_output) > 0:
        logger.info(f"Aligned file already exists: {fasta_output}. Skipping alignment.")
        return
    if not os.path.exists(MUSCLE_PATH):
        logger.error(f"MUSCLE executable not found at {MUSCLE_PATH}. Cannot perform alignment.")
        raise FileNotFoundError(MUSCLE_PATH)
    muscle_flag = "-super5" if use_super else "-align"
    muscle_cmd = [MUSCLE_PATH, muscle_flag, fasta_input, "-output", fasta_output]
    start_msg = f"{Fore.GREEN}[START ALIGNMENT] Running: {' '.join(muscle_cmd)}{Style.RESET_ALL}"
    print(start_msg)
    logger.info(start_msg)
    try:
        subprocess.run(muscle_cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        logger.info(f"MUSCLE alignment completed for {fasta_input}.")
    except subprocess.CalledProcessError as e:
        logger.error(f"MUSCLE alignment failed for {fasta_input}: {e.stderr.decode()}")
        raise
    except Exception as e:
        logger.error(f"MUSCLE alignment failed for {fasta_input}: {e}")
        raise

# ---------- CORE WORKER FUNCTION -----------------------------------------
def align_sequences_improved(group_df, group_name, output_folder, group_type, marker_name, super_threshold=50):
    if group_name is None or str(group_name).strip().lower() in BAD_VALUES:
        return None
    
    group_df_filtered = group_df[group_df["RNA Sequence"].str.len().between(1, 6000)]
    if len(group_df_filtered) < 2:
        logger.info(f"Skipping {group_type} '{group_name}' ({marker_name}): fewer than 2 valid sequences.")
        return None

    safe_name = "".join(c if c.isalnum() else "_" for c in str(group_name))
    fasta_input = os.path.join(output_folder, f"{marker_name}_{group_type}_{safe_name}_input.fasta")
    fasta_aligned = os.path.join(output_folder, f"{marker_name}_{group_type}_{safe_name}_aligned.fasta")

    # OPTIMIZATION: Use list comprehension with .to_dict('records') or zip for faster record creation
    ids = [row.get("ID") if row.get("ID") and pd.notna(row.get("ID")) else f"{safe_name}_{idx}" for idx, row in group_df_filtered.iterrows()]
    headers = [f"Marker={marker_name}; Main_Organism={row.get('Main_Organism','')}; Species={row.get('Species','')}" for _, row in group_df_filtered.iterrows()]
    seqs = group_df_filtered["RNA Sequence"].tolist()
    records = [SeqRecord(Seq(seq), id=str(identifier), description=header) for seq, identifier, header in zip(seqs, ids, headers)]
    
    try:
        generate_fasta_file(records, fasta_input)
    except Exception as e:
        logger.error(f"Failed to write FASTA for {group_type} '{group_name}': {e}")
        return None

    try:
        perform_muscle_alignment(fasta_input, fasta_aligned, use_super=len(records) > super_threshold)
    except Exception as e:
        logger.error(f"Alignment step failed for {group_type} '{group_name}'.")
        if not (SKIP_ALIGNMENT and os.path.exists(fasta_aligned) and os.path.getsize(fasta_aligned) > 0):
            return None

    if not os.path.exists(fasta_aligned) or os.path.getsize(fasta_aligned) == 0:
        logger.error(f"Aligned file missing or empty: {fasta_aligned}")
        return None

    try:
        records_aligned = list(SeqIO.parse(fasta_aligned, "fasta"))
        if not records_aligned: return None
    except Exception as e:
        logger.error(f"Error reading aligned FASTA for {group_type} '{group_name}': {e}")
        return None

    try:
        avg_sim, avg_ov, avg_ch = compute_alignment_metrics(records_aligned)
        nuc_dist = compute_nucleotide_change_distribution(records_aligned)
        nuc_div = compute_nucleotide_diversity(records_aligned)
        taj_d = compute_tajimas_d(records_aligned)
        cons_sites, var_sites, _, _ = analyze_conservation(records_aligned)
    except Exception as e:
        logger.error(f"Metric calculation failed for {group_type} '{group_name}': {e}")
        return None

    rec_summary = {
        "Marker": marker_name, "Group_Type": group_type, "Group_Name": group_name,
        "Num_Seq": len(records_aligned),
        "Avg_Similarity": round(avg_sim, 2) if avg_sim is not None else None,
        "Overlap_Percentage": round(avg_ov, 2) if avg_ov is not None else None,
        "Changes_Percentage": round(avg_ch, 2) if avg_ch is not None else None,
        "Nucleotide_Diversity": round(nuc_div, 4) if nuc_div is not None else None,
        "Tajimas_D": round(taj_d, 4) if taj_d is not None else None,
        "Conserved_Sites": cons_sites, "Variable_Sites": var_sites,
        "Aligned_File": fasta_aligned
    }
    if nuc_dist:
        rec_summary.update({f"{nt}_Percentage": round(p, 2) for nt, p in nuc_dist.items()})
    
    logger.info(f"Alignment processed for {group_type} '{group_name}' ({marker_name}).")
    return rec_summary

# ---------- PARALLEL RUN --------------------------------------------------
def run_alignment_tasks_parallel(data_align):
    futures_info = []
    # OPTIMIZATION: Use ProcessPoolExecutor for CPU-bound tasks to bypass the GIL
    max_workers = max(1, int(psutil.cpu_count(logical=False))) # Use physical cores for CPU-bound work
    logger.info(f"Using {max_workers} workers with ProcessPoolExecutor.")
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        for marker_name, df in data_align.items():
            base_out = os.path.join(OUTPUT_DIR, f"alignments_{marker_name}")
            os.makedirs(base_out, exist_ok=True)
            group_type = 'Family'
            out_folder = os.path.join(base_out, group_type.lower())
            os.makedirs(out_folder, exist_ok=True)
            
            df_valid = filter_valid_taxa(df, 'Family')
            for name, group in df_valid.groupby('Family'):
                if name is None or str(name).strip().lower() in BAD_VALUES:
                    continue
                
                # The can_submit_task check is less critical with a fixed process pool
                # but can be kept as a safeguard for memory.
                while not can_submit_task():
                    logger.info(Fore.YELLOW+"Waiting for resources..."+Style.RESET_ALL)
                    time.sleep(5)
                
                future = executor.submit(align_sequences_improved, group, name, out_folder, group_type, marker_name)
                futures_info.append((future, name, marker_name, group_type))

        for future, name, marker, gtype in futures_info:
            try:
                result = future.result()
                if result:
                    alignment_summary[result["Marker"]].append(result)
            except Exception as e:
                logger.error(f"An error occurred processing {gtype} '{name}' for {marker}: {e}")
    
    logger.info("\n[COMPLETE] Parallel alignment processing complete.")


# ---------- GENUS/SPECIES FROM FAMILY ------------------------------------
# Functions in this section call the optimized metric functions, so they benefit implicitly.
# No major structural changes are needed here.
def process_genus_species_from_family_file(family_aligned_path, marker_name):
    try:
        records_aligned = list(SeqIO.parse(family_aligned_path, "fasta"))
    except Exception as e:
        logger.error(f"Error reading family aligned file {family_aligned_path}: {e}")
        return []
    fname = os.path.basename(family_aligned_path)
    parts = fname.split('_')
    family_name = parts[2] if len(parts) >= 3 else ""
    genus_groups = defaultdict(list)
    species_groups = defaultdict(list)
    for rec in records_aligned:
        info = get_header_info(rec, marker_name, global_reference_info)
        genus = info.get("Main_Organism", "unknow")
        species = info.get("Species", "unknow")
        if genus not in BAD_VALUES:
            genus_groups[genus].append(rec)
            species_groups[(genus, species)].append(rec)
    results = []
    for group, recs in list(genus_groups.items()) + list(species_groups.items()):
        if len(recs) < 2:
            continue
        
        is_species = isinstance(group, tuple)
        group_type = "Species" if is_species else "Genus"
        group_name = format_species_name(group[0], group[1]) if is_species else format_species_name(group, "")

        avg_sim, avg_ov, _ = compute_alignment_metrics(recs)
        nuc_dist = compute_nucleotide_change_distribution(recs)
        nuc_div = compute_nucleotide_diversity(recs)
        taj_d = compute_tajimas_d(recs)
        cons_sites, var_sites, _, _ = analyze_conservation(recs)
        consensus_seq = compute_consensus(recs).replace("-", "")
        
        res = {
            "Marker": marker_name, "Group_Type": group_type, "Group_Name": group_name,
            "Num_Seq": len(recs),
            "Avg_Similarity": round(avg_sim,2) if avg_sim is not None else None,
            "Overlap_Percentage": round(avg_ov,2) if avg_ov is not None else None,
            "Nucleotide_Diversity": round(nuc_div,4) if nuc_div is not None else None,
            "Tajimas_D": round(taj_d,4) if taj_d is not None else None,
            "Conserved_Sites": cons_sites, "Variable_Sites": var_sites,
            "Family": family_name, "Consensus_Sequence": consensus_seq
        }
        if nuc_dist:
            res.update({f"{nt}_Percentage": round(p, 2) for nt, p in nuc_dist.items()})
        results.append(res)
    return results

def process_genus_species_from_family_alignments(data_align):
    for marker_name in data_align.keys():
        fam_folder = os.path.join(OUTPUT_DIR, f"alignments_{marker_name}", "family")
        if not os.path.exists(fam_folder):
            continue
        files = [os.path.join(fam_folder, f) for f in os.listdir(fam_folder) if f.endswith("_aligned.fasta")]
        for fpath in files:
            res_list = process_genus_species_from_family_file(fpath, marker_name)
            alignment_summary[marker_name].extend(res_list)
            for res in res_list:
                logger.info(f"Processed {res['Group_Type']} from {os.path.basename(fpath)}: {res['Group_Name']}")


# ---------- COVERAGE, SUMMARIES & REPORTS (mostly fast pandas operations) ---
# The remaining functions are largely based on pandas operations which are already optimized.
# No significant changes are required for performance in this section.
# I've included the original code for completeness.

def compute_consensus(aligned_records):
    if not aligned_records: return ""
    arr = np.array([list(str(rec.seq)) for rec in aligned_records])
    consensus_seq = ""
    for i in range(arr.shape[1]):
        col = arr[:, i]
        col_valid = col[col != '-']
        if col_valid.size > 0:
            unique, counts = np.unique(col_valid, return_counts=True)
            consensus_seq += unique[np.argmax(counts)]
    return consensus_seq

def build_coverage_df(taxa_dict, taxon_label):
    total_markers = len(DATA_FILES)
    all_taxa = [{"Taxon": taxon, "Marker Count": len(markers)} for taxon, markers in sorted(taxa_dict.items())]
    df_all = pd.DataFrame(all_taxa)
    df_in_all = df_all[df_all["Marker Count"] == total_markers].copy()
    df_in_two = df_all[df_all["Marker Count"] == 2].copy()
    
    df_all["Coverage Type"] = f"All {taxon_label}"
    df_in_all["Coverage Type"] = "In All Markers"
    df_in_two["Coverage Type"] = "In Two Markers"
    
    return pd.concat([df_all, df_in_all, df_in_two], ignore_index=True)

def analyze_coverage(data_coverage):
    taxa_markers = {level: defaultdict(set) for level in ['Family', 'Main_Organism', 'Species']}
    for marker, df in data_coverage.items():
        for level, col_name in zip(taxa_markers.keys(), ['Family', 'Main_Organism', 'Species']):
            unique_taxa = df[col_name].dropna().str.strip().str.lower().unique()
            for taxon in unique_taxa:
                if taxon not in BAD_VALUES:
                    taxa_markers[level][taxon].add(marker)

    return (build_coverage_df(taxa_markers['Family'], "Families"),
            build_coverage_df(taxa_markers['Main_Organism'], "Genera"),
            build_coverage_df(taxa_markers['Species'], "Species"))


def additional_summaries(data_align):
    df_family_cov, df_genus_cov, df_species_cov = analyze_coverage(prepare_coverage_data(data_align))
    
    all_summaries = []
    for group_col in ['f-h', 'cp']:
        summary_list = []
        for marker, df in data_align.items():
            if group_col in df.columns:
                valid_df = df[df[group_col].notna()]
                for col in ['Family', 'Main_Organism', 'Species']:
                    valid_df = valid_df[~valid_df[col].str.lower().isin(BAD_VALUES)]
                
                totals = {
                    "Total Families": valid_df['Family'].nunique(),
                    "Total Genera": valid_df['Main_Organism'].nunique(),
                    "Total Species": valid_df['Species'].nunique()
                }

                for val, sub_df in valid_df.groupby(valid_df[group_col].str.strip()):
                    summary_list.append({
                        "Marker": marker, f"{group_col.upper()} Group": val,
                        "Unique Families": sub_df['Family'].nunique(),
                        "Unique Genera": sub_df['Main_Organism'].nunique(),
                        "Unique Species": sub_df['Species'].nunique(),
                        **totals
                    })
        all_summaries.append(pd.DataFrame(summary_list))

    # Genetic Summary
    summary_data = []
    all_taxa = {level: set() for level in ['Family', 'Main_Organism', 'Species']}
    for marker, df in data_align.items():
        marker_counts = {"Marker": marker}
        for level, col in zip(['Families', 'Genera', 'Species'], ['Family', 'Main_Organism', 'Species']):
            valid_taxa = df.loc[~df[col].str.lower().isin(BAD_VALUES), col].dropna().unique()
            marker_counts[level] = len(valid_taxa)
            all_taxa[col].update(valid_taxa)
        summary_data.append(marker_counts)
    
    total_counts = {"Marker": "Total"}
    for level, col in zip(['Families', 'Genera', 'Species'], ['Family', 'Main_Organism', 'Species']):
        total_counts[level] = len(all_taxa[col])
    summary_data.append(total_counts)
    genetic_summary_df = pd.DataFrame(summary_data)

    return all_summaries[0], all_summaries[1], df_family_cov, df_genus_cov, df_species_cov, genetic_summary_df


def generate_final_report(data_align, data_coverage):
    f_h_df, cp_df, fam_cov, gen_cov, sp_cov, gen_sum = additional_summaries(data_align)

    output_excel = os.path.join(OUTPUT_DIR, "genetic_summary_optimized.xlsx")
    with pd.ExcelWriter(output_excel, engine='openpyxl') as writer:
        gen_sum.to_excel(writer, sheet_name="Genetic Summary", index=False)
        
        for marker, records in alignment_summary.items():
            if not records: continue
            df = pd.DataFrame(records)
            for group_type in ["Family", "Genus", "Species"]:
                df_group = df[df['Group_Type'] == group_type]
                if not df_group.empty:
                    sheet_name = f"{marker} {group_type}"[:31]
                    df_group.to_excel(writer, sheet_name=sheet_name, index=False)
        
        fam_cov.to_excel(writer, sheet_name="Family Coverage", index=False)
        gen_cov.to_excel(writer, sheet_name="Genus Coverage", index=False)
        sp_cov.to_excel(writer, sheet_name="Species Coverage", index=False)
        if not f_h_df.empty:
            f_h_df.to_excel(writer, sheet_name="F-H Group Summary", index=False)
        if not cp_df.empty:
            cp_df.to_excel(writer, sheet_name="CP Group Summary", index=False)

    logger.info(f"\n[COMPLETE] Final report saved to {output_excel}")


def generate_unique_groups_excel(alignment_summary):
    all_records = []
    for marker, recs in alignment_summary.items():
        for rec in recs:
            rec_copy = rec.copy()
            rec_copy["Marker"] = marker
            all_records.append(rec_copy)
    
    if not all_records: return
    
    df = pd.DataFrame(all_records)
    output_excel = os.path.join(OUTPUT_DIR, "unique_groups_summary.xlsx")
    with pd.ExcelWriter(output_excel, engine='openpyxl') as writer:
        df.to_excel(writer, sheet_name="Unique Groups Summary", index=False)
    logger.info(f"\n[COMPLETE] Unique groups summary saved to {output_excel}")

def generate_consensus_excel(alignment_summary, output_excel):
    with pd.ExcelWriter(output_excel, engine='openpyxl') as writer:
        all_records = []
        for marker, recs in alignment_summary.items():
            for rec in recs:
                rec_copy = rec.copy()
                rec_copy["Marker"] = marker
                all_records.append(rec_copy)
        
        if not all_records: return
        df = pd.DataFrame(all_records)

        for marker in df['Marker'].unique():
            df_marker = df[df['Marker'] == marker]
            for group_type in ['Family', 'Genus', 'Species']:
                df_group = df_marker[df_marker['Group_Type'] == group_type]
                if not df_group.empty:
                    consensus_df = df_group[["Group_Name", "Consensus_Sequence", "Num_Seq"]].copy()
                    consensus_df['Length'] = consensus_df['Consensus_Sequence'].str.len()
                    sheet_name = f"{marker} {group_type}"[:31]
                    consensus_df.to_excel(writer, sheet_name=sheet_name, index=False)
    logger.info(f"Consensus Excel file saved to {output_excel}")

# ---------- MAIN ----------------------------------------------------------
def main():
    logger.info(f"{Fore.CYAN}--- PROCESS START ---{Style.RESET_ALL}")
    start_time = time.time()

    init_gpu_optional()
    data_align = load_data()
    save_updated_files(data_align)
    
    global global_reference_info
    global_reference_info = build_reference_info(data_align)
    
    data_coverage = prepare_coverage_data(data_align)
    
    run_alignment_tasks_parallel(data_align)
    process_genus_species_from_family_alignments(data_align)
    
    generate_final_report(data_align, data_coverage)
    generate_unique_groups_excel(alignment_summary)
    
    consensus_excel_path = os.path.join(OUTPUT_DIR, "consensus_sequences.xlsx")
    generate_consensus_excel(alignment_summary, consensus_excel_path)
    
    # Global analysis section remains the same
    for marker in DATA_FILES.keys():
        try:
            records = []
            for _, row in data_align[marker].iterrows():
                records.append(SeqRecord(Seq(row["RNA Sequence"]), id=str(row.get("ID", ""))))
            if len(records) > 1:
                analysis_folder = os.path.join(OUTPUT_DIR, f"global_analysis_{marker}")
                os.makedirs(analysis_folder, exist_ok=True)
                _, _, _, cons_scores = analyze_conservation(records)
                if cons_scores is not None:
                    plot_conservation(cons_scores, marker, analysis_folder)
        except Exception as e:
            logger.error(f"Global analysis failed for {marker}: {e}")

    end_time = time.time()
    logger.info(f"{Fore.CYAN}--- PROCESS COMPLETE IN {end_time - start_time:.2f} SECONDS ---{Style.RESET_ALL}")

if __name__=="__main__":
    # On Windows, ProcessPoolExecutor requires this guard
    main()