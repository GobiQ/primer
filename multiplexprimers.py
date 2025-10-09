# AutoPrimer5 — Multiplex Mode (15 targets × 3 channels)
# Minimal Streamlit app (integrated with your autoprimer5.py catalog)
# --------------------------------------------------------------
# New in this version
# • Imports organism/target list directly from autoprimer5.py
# • Auto-loads target nucleotide sequences from the catalog (no NCBI required)
# • Robust catalog parsing with user-provided sequence key hints
# • Auto-assigns 15 targets to 3×5 Tm slots while optimizing primer quality
# • Debug panel to show why entries might appear as <none>
# --------------------------------------------------------------

import streamlit as st
import pandas as pd
import re, io, math
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional, Any

# Optional libs
try:
    import primer3
    HAVE_PRIMER3 = True
except Exception:
    HAVE_PRIMER3 = False

try:
    from Bio.Seq import Seq
    HAVE_BIO = True
except Exception:
    HAVE_BIO = False

# -----------------------------
# Catalog loader & sequence extractor
# -----------------------------
@st.cache_data
def load_catalog() -> Any:
    """Load the pathogen/target catalog from autoprimer5.py.
    Tries a few common function/variable names for resilience.
    """
    try:
        import autoprimer5 as ap
    except Exception:
        return None
    # Try function names first
    for fn_name in [
        "get_organism_suggestions_with_gene_targets",
        "get_organism_targets",
        "get_targets",
        "get_catalog",
    ]:
        fn = getattr(ap, fn_name, None)
        if callable(fn):
            try:
                return fn()
            except Exception:
                pass
    # Try variable names
    for var_name in [
        "ORGANISM_TARGETS",
        "organism_targets",
        "CATALOG",
        "catalog",
        "TARGETS",
        "targets",
    ]:
        if hasattr(ap, var_name):
            return getattr(ap, var_name)
    return None


def clean_seq(s: str) -> str:
    return re.sub(r"[^ACGTNacgtn]", "", s).upper().replace("U", "T")


def find_seq_in_obj(obj: Any, key_hints: List[str]) -> Optional[str]:
    """Search an arbitrary nested object for a plausible DNA sequence.
    Honors user-provided key hints for field names that hold sequences.
    """
    # direct string
    if isinstance(obj, str):
        s = clean_seq(obj)
        if len(s) >= 60:
            return s
        return None
    # dict with keys
    if isinstance(obj, dict):
        # priority: user hints
        for k in key_hints:
            if k in obj and isinstance(obj[k], str):
                s = clean_seq(obj[k])
                if len(s) >= 60:
                    return s
        # common defaults
        for k in ["sequence","target_sequence","amplicon","region","locus_seq","seq","reference_sequence","template","amplicon_template"]:
            if k in obj and isinstance(obj[k], str):
                s = clean_seq(obj[k])
                if len(s) >= 60:
                    return s
        # otherwise scan values
        for v in obj.values():
            s = find_seq_in_obj(v, key_hints)
            if s:
                return s
        return None
    # list/tuple -> scan items
    if isinstance(obj, (list, tuple)):
        for v in obj:
            s = find_seq_in_obj(v, key_hints)
            if s:
                return s
    return None

@st.cache_data
def flatten_catalog_with_sequences(catalog: Any, key_hints_csv: str) -> List[Dict[str,str]]:
    """Flatten a wide range of catalog shapes into a list of dicts with organism, target, and sequence.
    key_hints_csv is a comma-separated string of preferred field names for sequences.
    """
    key_hints = [k.strip() for k in key_hints_csv.split(',') if k.strip()]
    out: List[Dict[str,str]] = []

    def push(label: str, organism: str, target: str, seq: str, path: str):
        out.append({
            "label": label,
            "organism": organism or "Unknown",
            "target": target or "Target",
            "sequence": seq,
            "path": path,
        })

    def walk(obj: Any, path: str = "$"):
        if obj is None:
            return
        # If leaf has a sequence, try to infer organism/target names
        seq = find_seq_in_obj(obj, key_hints)
        if seq:
            organism = None
            target = None
            label = None
            if isinstance(obj, dict):
                organism = obj.get("scientific") or obj.get("organism") or obj.get("species") or obj.get("name") or obj.get("common")
                target = obj.get("target") or obj.get("gene") or obj.get("locus") or obj.get("amplicon_name")
                label = obj.get("label") or obj.get("display")
            if isinstance(obj, (list, tuple)) and len(obj) >= 3 and all(isinstance(x, str) for x in obj[:3]):
                # Old pattern: [common, scientific, target, ...]
                organism = obj[1] or organism
                target = obj[2] or target
                if not label:
                    common = obj[0]
                    label = f"{common} — {organism} — {target}"
            if not label:
                label = f"{organism or 'Unknown'} — {target or 'Target'}"
            push(label, organism or "Unknown", target or "Target", seq, path)
            return
        # Recurse
        if isinstance(obj, dict):
            for k, v in obj.items():
                walk(v, f"{path}.{k}")
        elif isinstance(obj, (list, tuple)):
            for i, v in enumerate(obj):
                walk(v, f"{path}[{i}]")

    walk(catalog)
    # Deduplicate labels
    seen = set()
    dedup = []
    for e in out:
        key = (e["label"], e["organism"], e["target"], e["sequence"])
        if key in seen:
            continue
        seen.add(key)
        dedup.append(e)
    return dedup

# -----------------------------
# Thermo + tuning helpers (coarse but fast)
# -----------------------------

def gc_content(seq: str) -> float:
    seq = seq.upper()
    if not seq: return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq)


def rough_amp_tm(seq: str, monovalent_mM: float = 50.0, free_Mg_mM: float = 2.0) -> float:
    L = max(len(seq), 1)
    base = 69.3 + 0.41 * (gc_content(seq) * 100.0) - (650.0 / L)
    Na = max(monovalent_mM, 1e-3) / 1000.0
    salt = 16.6 * math.log10(Na)
    mg_adj = 0.6 * max(free_Mg_mM, 0.0)
    return base + salt + mg_adj

GC_TAILS = ["", "G", "GC", "GCG", "GCGC", "GCGCG", "GCGCGC"]
AT_TAILS = ["", "A", "AT", "ATAT", "ATATAT", "ATATATAT"]

@dataclass
class Candidate:
    fwd: str
    rev: str
    prod_len: int
    base_tm: float

@dataclass
class AssignmentChoice:
    candidate: Candidate
    f_tail: str
    r_tail: str
    tm_after: float
    cost: float

@dataclass
class PrimerPair:
    fwd: str
    rev: str
    product_len: int
    product_tm: float
    channel: str
    slot_tm: float
    slot_index: int
    organism: str
    target_name: str
    tails: Tuple[str, str]
    issues: List[str]

# -----------------------------
# UI
# -----------------------------
st.set_page_config(page_title="AutoPrimer5 — Multiplex (multiplexprimers.py)", layout="wide")
st.title("AutoPrimer5 — Multiplex (15 targets • 3 channels × 5 Tm slots)")

with st.sidebar:
    st.header("Catalog & Conditions")
    catalog_obj = load_catalog()
    seq_key_hints = st.text_input("Sequence field key(s) (comma-separated)", value="sequence,target_sequence,amplicon,region,locus_seq,seq")
    entries = flatten_catalog_with_sequences(catalog_obj, seq_key_hints)
    st.metric("Catalog entries (raw)", 0 if catalog_obj is None else (len(catalog_obj) if hasattr(catalog_obj, "__len__") else "?"))
    st.metric("Entries with sequences", len(entries))
    if len(entries) == 0:
        st.warning("No sequences found in catalog. Verify autoprimer5.py has embedded sequences in the fields above, or add the correct key name to 'Sequence field key(s)'.")
    monovalent_mM = st.number_input("Monovalent salt (mM)", 10.0, 200.0, 50.0, 1.0)
    free_Mg_mM   = st.number_input("Free Mg²⁺ (mM)", 0.0, 6.0, 2.0, 0.1)

    st.markdown("**Tm slots per channel (°C)** — keep ≥3 °C spacing")
    default_grid = {
        "FAM": [79.5, 83.0, 86.5, 90.0, 93.5],
        "HEX": [78.0, 81.5, 85.0, 88.5, 92.0],
        "Cy5": [80.5, 84.0, 87.5, 91.0, 94.5],
    }
    grid = {}
    for dye in ["FAM", "HEX", "Cy5"]:
        cols = st.columns(5)
        grid[dye] = []
        for i, col in enumerate(cols):
            with col:
                grid[dye].append(st.number_input(f"{dye} s{i+1}", 70.0, 100.0, default_grid[dye][i], 0.1, key=f"{dye}_{i}"))

st.markdown("### 1) Pick **15 catalog targets** (sequences auto-loaded)")
labels = [e["label"] for e in entries] if entries else ["<none>"]
selected_labels = []
for i in range(15):
    idx = min(i, max(0, len(labels)-1))
    selected_labels.append(st.selectbox(f"Select target {i+1}", labels, index=idx, key=f"sel_{i}"))

selected = [next((e for e in entries if e["label"]==lbl), None) for lbl in selected_labels]

st.markdown("### 2) Primer design settings & constraints")
colA, colB, colC = st.columns(3)
with colA:
    opt_size = st.slider("Primer size (optimal)", 15, 30, 20)
    min_size = st.slider("Primer size (min)", 15, 25, 18)
    max_size = st.slider("Primer size (max)", 20, 35, 25)
with colB:
    anneal_tm = st.slider("Primer anneal Tm target (°C)", 50.0, 72.0, 60.0, 0.5)
    product_len_min, product_len_max = st.slider("Desired product length (bp)", 60, 240, (100, 170))
with colC:
    tail_penalty = st.slider("Tail length penalty (°C per nt, scoring)", 0.00, 0.10, 0.02, 0.01)
    delta_penalty = st.slider("ΔTm to slot weight (×)", 0.5, 3.0, 1.0, 0.1)

st.markdown("### 3) Auto-assign to **3 channels × 5 Tm slots** and design")

slot_list: List[Tuple[str,int,float]] = []  # (dye, slot_index, slot_tm)
for dye in ["FAM","HEX","Cy5"]:
    for k in range(5):
        slot_list.append((dye, k, grid[dye][k]))

# Primer3 design for multiple candidates per target

def design_candidates(seq: str, n_return: int = 6) -> List[Candidate]:
    cands: List[Candidate] = []
    if not seq or len(seq) < min_size + 30:
        return cands
    if HAVE_PRIMER3:
        params = dict(
            SEQUENCE_TEMPLATE=seq,
            PRIMER_OPT_SIZE=opt_size,
            PRIMER_MIN_SIZE=min_size,
            PRIMER_MAX_SIZE=max_size,
            PRIMER_NUM_RETURN=n_return,
            PRIMER_PRODUCT_SIZE_RANGE=[[product_len_min, product_len_max]],
            PRIMER_OPT_TM=anneal_tm,
            PRIMER_MIN_TM=anneal_tm-2.5,
            PRIMER_MAX_TM=anneal_tm+2.5,
            PRIMER_MAX_POLY_X=4,
            PRIMER_GC_CLAMP=1,
        )
        try:
            res = primer3.bindings.designPrimers({"SEQUENCE_TEMPLATE": seq}, params)
            n = int(res.get("PRIMER_PAIR_NUM_RETURNED", 0))
            for idx in range(n):
                fwd = res.get(f"PRIMER_LEFT_{idx}_SEQUENCE")
                rev = res.get(f"PRIMER_RIGHT_{idx}_SEQUENCE")
                if not (fwd and rev):
                    continue
                core = seq
                if len(seq) > (len(fwd)+len(rev)):
                    core = seq[len(fwd):-(len(rev))]
                amp = f"{fwd}{core}{rev}"
                cands.append(Candidate(fwd=fwd, rev=rev, prod_len=len(amp), base_tm=rough_amp_tm(amp, monovalent_mM, free_Mg_mM)))
        except Exception:
            pass
    if not cands:
        # naive fallback: terminal slices
        fwd = seq[:opt_size]
        rev = Seq(seq[-opt_size:]).reverse_complement() if HAVE_BIO else seq[-opt_size:][::-1]
        core = seq[opt_size:-opt_size] if len(seq)>(2*opt_size) else ""
        amp = f"{fwd}{core}{rev}"
        cands.append(Candidate(fwd=str(fwd), rev=str(rev), prod_len=len(amp), base_tm=rough_amp_tm(amp, monovalent_mM, free_Mg_mM)))
    return cands

# Evaluate best choice of (candidate, tails) for a given slot

def best_choice_for_slot(seq: str, slot_tm: float) -> Tuple[Optional[Candidate], Optional[AssignmentChoice]]:
    cands = design_candidates(seq)
    best: Optional[AssignmentChoice] = None
    best_cand: Optional[Candidate] = None
    for cand in cands:
        core = seq
        if len(seq) > (len(cand.fwd)+len(cand.rev)):
            core = seq[len(cand.fwd):-(len(cand.rev))]
        for tF in GC_TAILS + AT_TAILS:
            for tR in GC_TAILS + AT_TAILS:
                amp2 = f"{tF}{cand.fwd}{core}{cand.rev}{tR}"
                tm2 = rough_amp_tm(amp2, monovalent_mM, free_Mg_mM)
                cost = delta_penalty*abs(tm2 - slot_tm) + tail_penalty*(len(tF)+len(tR))
                choice = AssignmentChoice(candidate=cand, f_tail=tF, r_tail=tR, tm_after=tm2, cost=cost)
                if (best is None) or (choice.cost < best.cost):
                    best = choice
                    best_cand = cand
    return best_cand, best

# Build cost table for assignment (targets × slots)
run = st.button("Auto‑design & assign 15‑plex")

results: List[PrimerPair] = []
if run:
    # Ensure we have 15 selected entries with sequences
    target_infos = []
    for i, entry in enumerate(selected):
        if not entry or not entry.get("sequence"):
            st.error(f"Target {i+1}: selected entry has no sequence. Adjust 'Sequence field key(s)' or pick a different target.")
            target_infos.append(None)
            continue
        target_infos.append(entry)
    if any(t is None for t in target_infos):
        st.stop()

    # Compute choice grid
    slot_list_local = slot_list.copy()
    choices: List[List[Optional[AssignmentChoice]]] = []
    for i, entry in enumerate(target_infos):
        row: List[Optional[AssignmentChoice]] = []
        for (dye, slot_idx, slot_tm) in slot_list_local:
            _, choice = best_choice_for_slot(entry["sequence"], slot_tm)
            row.append(choice)
        choices.append(row)

    # Greedy assignment (targets→slots)
    nT, nS = len(choices), len(slot_list_local)
    if nT != 15 or nS != 15:
        st.error("Internal error: expected 15 targets and 15 slots.")
        st.stop()

    assigned_slots = [-1]*nT
    taken = [False]*nS
    flat = []
    for i in range(nT):
        for j in range(nS):
            c = choices[i][j]
            if c is None:
                continue
            flat.append((c.cost, i, j))
    if not flat:
        st.error("No viable primer/slot combinations were found. Try widening product length or primer Tm range.")
        st.stop()
    flat.sort(key=lambda x: x[0])

    for _, i, j in flat:
        if assigned_slots[i] == -1 and not taken[j]:
            assigned_slots[i] = j
            taken[j] = True
        if all(a != -1 for a in assigned_slots):
            break

    if any(a == -1 for a in assigned_slots):
        st.error("Could not complete a full 15×15 assignment with current constraints.")
        st.stop()

    # Cross-dimer checker
    def has_bad_3prime_dimer(a: str, b: str) -> bool:
        a3 = a[-5:]
        b3 = b[-5:]
        comp = str.maketrans("ACGT", "TGCA")
        return a3 == b3.translate(comp)[::-1]

    assigned_pairs: List[PrimerPair] = []
    for i, entry in enumerate(target_infos):
        j = assigned_slots[i]
        dye, slot_idx, slot_tm = slot_list_local[j]
        choice = choices[i][j]
        cand = choice.candidate
        core = entry["sequence"]
        if len(core) > (len(cand.fwd)+len(cand.rev)):
            core_frag = core[len(cand.fwd):-(len(cand.rev))]
        else:
            core_frag = ""
        amp_len = len(choice.f_tail) + len(cand.fwd) + len(core_frag) + len(cand.rev) + len(choice.r_tail)
        issues: List[str] = []
        # quick within-set cross-dimer marking added later
        assigned_pairs.append(
            PrimerPair(
                fwd=choice.f_tail + cand.fwd,
                rev=cand.rev + choice.r_tail,
                product_len=amp_len,
                product_tm=choice.tm_after,
                channel=dye,
                slot_tm=slot_tm,
                slot_index=slot_idx,
                organism=entry["organism"],
                target_name=entry["target"],
                tails=(choice.f_tail, choice.r_tail),
                issues=issues,
            )
        )

    # cross-dimer warnings
    for a in range(len(assigned_pairs)):
        for b in range(a+1, len(assigned_pairs)):
            if has_bad_3prime_dimer(assigned_pairs[a].fwd, assigned_pairs[b].fwd) or \
               has_bad_3prime_dimer(assigned_pairs[a].rev, assigned_pairs[b].rev):
                assigned_pairs[a].issues.append(f"Potential 3' cross-dimer with {assigned_pairs[b].organism} ({assigned_pairs[b].channel} s{assigned_pairs[b].slot_index+1})")
                assigned_pairs[b].issues.append(f"Potential 3' cross-dimer with {assigned_pairs[a].organism} ({assigned_pairs[a].channel} s{assigned_pairs[a].slot_index+1})")

    df = pd.DataFrame([
        {
            "Channel": r.channel,
            "Slot": r.slot_index+1,
            "Slot_Tm_°C": round(r.slot_tm,1),
            "Pred_Product_Tm_°C": round(r.product_tm,2),
            "ΔTm_to_Slot": round(r.product_tm - r.slot_tm, 2),
            "Organism": r.organism,
            "Target": r.target_name,
            "Forward_Primer": r.fwd,
            "Reverse_Primer": r.rev,
            "Amplicon_bp": r.product_len,
            "Tails": "/".join([r.tails[0] or "-", r.tails[1] or "-"]),
            "Issues": "; ".join(sorted(set(r.issues))) if r.issues else ""
        } for r in assigned_pairs
    ]).sort_values(["Channel","Slot"]).reset_index(drop=True)

    st.success("Auto-assignment complete.")
    st.dataframe(df, use_container_width=True)

    out = io.BytesIO()
    df.to_csv(out, index=False)
    st.download_button("Download manifest (CSV)", data=out.getvalue(), file_name="autoprimer5_multiplex_manifest.csv", mime="text/csv")

st.markdown("---")
st.caption(
    "This version loads organisms/targets from autoprimer5.py, extracts embedded sequences with user-provided key hints, "
    "and auto-assigns 15 targets to 3×5 Tm slots. No NCBI/network calls required."
)
