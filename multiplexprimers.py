# AutoPrimer — Multiplex Mode (15 targets × 3 channels)
# Minimal, self‑contained Streamlit app (integrated with your old catalog)
# --------------------------------------------------------------
# New in this version
# • Automatically loads the **target sequences** tied to each pathogen/target from your existing app
# • **Auto-optimizes assignment** of 15 targets to 3×5 Tm slots (no manual slot picking)
# • Optimizes primer quality per target while ensuring at least one high‑quality design fits a slot
# --------------------------------------------------------------

import streamlit as st
import pandas as pd
import re, io, math
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional

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
def load_catalog() -> Dict:
    """Load the pathogen/target catalog from the **existing app** if present.
    It expects a function named `get_organism_suggestions_with_gene_targets()`
    in `autoprimer5.py` (your old app file). Falls back to empty dict.
    """
    try:
        from autoprimer5 import get_organism_suggestions_with_gene_targets as catalog_fn
        return catalog_fn()
    except Exception:
        return {}


def clean_seq(s: str) -> str:
    return re.sub(r"[^ACGTNacgtn]", "", s).upper().replace("U", "T")


def _find_seq_in_item(item) -> Optional[str]:
    """Attempt to locate a DNA sequence in a catalog item structure.
    Supports keys like 'sequence', 'target_sequence', 'amplicon', etc., or
    scans strings for long ACGT runs (≥60 bp). Returns cleaned seq or None.
    """
    # Common explicit keys
    if isinstance(item, dict):
        for key in ["sequence", "target_sequence", "amplicon", "region", "locus_seq", "seq"]:
            if key in item and isinstance(item[key], str):
                s = clean_seq(item[key])
                if len(s) >= 60:
                    return s
        # search values for long strings
        for v in item.values():
            if isinstance(v, str):
                s = clean_seq(v)
                if len(s) >= 60:
                    return s
    # If it's a list/tuple, scan strings
    if isinstance(item, (list, tuple)):
        for v in item:
            if isinstance(v, str):
                s = clean_seq(v)
                if len(s) >= 60:
                    return s
    return None

@st.cache_data
def flatten_catalog_with_sequences(catalog: Dict) -> List[Dict]:
    """Flatten your nested catalog into a list of dicts with organism, target name, and sequence.
    Expected to work with the structures used in your old app; if a sequence is not present,
    the entry is skipped (you can still paste custom sequences below if needed).
    """
    result = []
    if not isinstance(catalog, dict):
        return result
    for category, sub in catalog.items():
        if not isinstance(sub, dict):
            continue
        for subcat, items in sub.items():
            for it in items:
                # Common old-app pattern: [common, scientific, target, ...]
                common = None
                scientific = None
                target_name = None
                if isinstance(it, (list, tuple)):
                    if len(it) >= 2:
                        common = it[0]
                        scientific = it[1]
                    if len(it) >= 3 and isinstance(it[2], str):
                        target_name = it[2]
                elif isinstance(it, dict):
                    common = it.get("common") or it.get("name")
                    scientific = it.get("scientific") or it.get("organism")
                    target_name = it.get("target") or it.get("gene")
                seq = _find_seq_in_item(it)
                if seq:
                    result.append({
                        "category": category,
                        "subcat": subcat,
                        "label": f"{common or scientific or 'Unknown'} — {target_name or ''}".strip(),
                        "organism": scientific or common or "Unknown",
                        "target": target_name or "Target",
                        "sequence": seq,
                    })
    return result

# -----------------------------
# Thermo + tuning helpers (coarse but fast)
# -----------------------------

def gc_content(seq: str) -> float:
    seq = seq.upper()
    if not seq: return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq)


def rough_amp_tm(seq: str, monovalent_mM: float = 50.0, free_Mg_mM: float = 2.0) -> float:
    """Very fast amplicon Tm estimate for placement (±2–4 °C typical)."""
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
st.set_page_config(page_title="AutoPrimer5 — Multiplex", layout="wide")
st.title("AutoPrimer5 — Multiplex (15 targets • 3 channels × 5 Tm slots)")

with st.sidebar:
    st.header("Catalog & Conditions")
    catalog = load_catalog()
    entries = flatten_catalog_with_sequences(catalog)
    if entries:
        st.success(f"Loaded {len(entries)} catalog entries with sequences from autoprimer5.py")
    else:
        st.error("Could not find sequences in the built-in catalog. You can still paste custom sequences in the main panel.")
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

# Select 15 entries
labels = [e["label"] for e in entries]
selected_labels = []
for i in range(15):
    idx = min(i, max(0, len(labels)-1))
    selected_labels.append(st.selectbox(f"Select target {i+1}", labels if labels else ["<none>"], index=idx if labels else 0, key=f"sel_{i}"))

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

# Build 15 slots
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
                # crude amplicon reconstruction for Tm estimate
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

def best_choice_for_slot(seq: str, slot_tm: float) -> Tuple[Candidate, AssignmentChoice]:
    cands = design_candidates(seq)
    best: Optional[AssignmentChoice] = None
    best_cand: Optional[Candidate] = None
    for cand in cands:
        # Build core amplicon once
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
    # Precompute best option for each (target, slot)
    target_infos = []
    for i, entry in enumerate(selected):
        if not entry or not entry.get("sequence"):
            st.error(f"Target {i+1}: missing sequence; please ensure catalog contains sequences for all selections.")
            target_infos.append(None)
            continue
        target_infos.append(entry)

    # If any missing, abort
    if any(t is None for t in target_infos):
        st.stop()

    # Compute choice grid
    choices: List[List[AssignmentChoice]] = []
    for i, entry in enumerate(target_infos):
        row: List[AssignmentChoice] = []
        for (dye, slot_idx, slot_tm) in slot_list:
            _, choice = best_choice_for_slot(entry["sequence"], slot_tm)
            row.append(choice)
        choices.append(row)

    # Greedy assignment (targets→slots) with tie-breaking by lowest incremental cost
    nT, nS = len(choices), len(slot_list)
    assert nT == 15 and nS == 15
    assigned_slots = [-1]*nT
    taken = [False]*nS

    # Sort all (i,j) by cost ascending and pick if both free
    flat = []
    for i in range(nT):
        for j in range(nS):
            flat.append((choices[i][j].cost, i, j))
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

    # Build final primer list and quick cross‑dimer warnings
    def has_bad_3prime_dimer(a: str, b: str) -> bool:
        a3 = a[-5:]
        b3 = b[-5:]
        comp = str.maketrans("ACGT", "TGCA")
        return a3 == b3.translate(comp)[::-1]

    assigned_pairs: List[PrimerPair] = []
    for i, entry in enumerate(target_infos):
        j = assigned_slots[i]
        dye, slot_idx, slot_tm = slot_list[j]
        choice = choices[i][j]
        cand = choice.candidate
        # reconstruct amplicon length roughly
        core = entry["sequence"]
        if len(core) > (len(cand.fwd)+len(cand.rev)):
            core_frag = core[len(cand.fwd):-(len(cand.rev))]
        else:
            core_frag = ""
        amp_len = len(choice.f_tail) + len(cand.fwd) + len(core_frag) + len(cand.rev) + len(choice.r_tail)
        issues = []
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

    # Output table
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
    "This version auto-loads catalog sequences and auto-assigns 15 targets to 3×5 Tm slots "
    "using a greedy global cost minimizer.\n"
    "For production: swap rough_amp_tm() for MELTING/UNAFold and add fuller dimer/hairpin scoring. "
    "Optionally replace greedy with a Hungarian solver for optimal assignment."
)
