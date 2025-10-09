# autoprimer5 — Multiplex Mode (15 targets × 3 channels)
# Minimal, self‑contained Streamlit app (integrated with your old catalog)
# --------------------------------------------------------------
# Purpose
# • Let user input/select **15 target organisms** at once
# • Keep your existing pathogen/target catalog that lived inside autoprimer5.py
# • Design 15 primer pairs arranged as **3 dye channels × 5 distinct melt Tm slots**
# • Tune **product Tm** with short non‑templated tails (does not change primer anneal Tm)
# • Export a manifest (CSV) with channels, slots, primers, predicted product Tm, and QC flags
# --------------------------------------------------------------

import streamlit as st
import pandas as pd
import re, io, math
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional

# Optional libs (app works in a degraded mode without them)
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
# Catalog loader: import from your existing app file if available
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

# -----------------------------
# Thermo + tuning helpers (coarse but fast)
# -----------------------------

def clean_seq(s: str) -> str:
    return re.sub(r"[^ACGTNacgtn]", "", s).upper().replace("U", "T")


def gc_content(seq: str) -> float:
    seq = seq.upper()
    if not seq: return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq)


def rough_amp_tm(seq: str, monovalent_mM: float = 50.0, free_Mg_mM: float = 2.0) -> float:
    """Very fast amplicon Tm estimate for placement (±2–4 °C typical).
    Based on GC%, length, log([Na+]), and a small Mg2+ boost.
    Replace with MELTING/UNAFold in backends for high accuracy.
    """
    L = max(len(seq), 1)
    base = 69.3 + 0.41 * (gc_content(seq) * 100.0) - (650.0 / L)
    Na = max(monovalent_mM, 1e-3) / 1000.0
    salt = 16.6 * math.log10(Na)
    mg_adj = 0.6 * max(free_Mg_mM, 0.0)
    return base + salt + mg_adj

GC_TAILS = ["G", "GC", "GCG", "GCGC", "GCGCG", "GCGCGC"]
AT_TAILS = ["A", "AT", "ATAT", "ATATAT", "ATATATAT"]

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
    tails: Tuple[str, str]  # (fwd_tail, rev_tail)
    issues: List[str]

# -----------------------------
# UI
# -----------------------------
st.set_page_config(page_title="AutoPrimer5 — Multiplex", layout="wide")
st.title("AutoPrimer5 — Multiplex (15 targets • 3 channels × 5 Tm slots)")

with st.sidebar:
    st.header("Catalog & Conditions")
    catalog = load_catalog()
    if catalog:
        st.success("Loaded pathogen/target catalog from autoprimer5.py")
    else:
        st.info("Could not auto-load catalog from autoprimer5.py. Proceed with manual entries or paste sequences.")
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

st.markdown("### 1) Select **15 target organisms** and provide one sequence per target")
left, right = st.columns([1.2, 1.0])

# Helper: flatten catalog into a selectable list
cat_options = []  # (label, scientific)
if isinstance(catalog, dict):
    for category, sub in catalog.items():
        if not isinstance(sub, dict):
            continue
        for subcat, items in sub.items():
            for it in items:
                if len(it) >= 2:
                    common = it[0]
                    scientific = it[1]
                    label = f"{common} — {scientific} ({category} / {subcat})"
                    cat_options.append((label, scientific))

with left:
    n_targets = 15
    organisms: List[str] = []
    seqs: List[str] = []
    for i in range(n_targets):
        st.markdown(f"**Target {i+1}**")
        if cat_options:
            pick = st.selectbox(
                f"Pick organism {i+1}",
                [lbl for (lbl, _) in cat_options],
                index=min(i, len(cat_options)-1),
                key=f"pick_{i}"
            )
            sci = next(s for (lbl, s) in cat_options if lbl == pick)
            org = sci
        else:
            org = st.text_input(f"Organism {i+1} (scientific/common)", key=f"org_{i}")
        seq = st.text_area(f"Target sequence {i+1} (amplicon region or locus segment)", height=80, key=f"seq_{i}")
        organisms.append(org.strip())
        seqs.append(clean_seq(seq))
        st.divider()

with right:
    st.markdown("**Primer design settings (kept compact)**")
    opt_size = st.slider("Primer size (optimal)", 15, 30, 20)
    min_size = st.slider("Primer size (min)", 15, 25, 18)
    max_size = st.slider("Primer size (max)", 20, 35, 25)
    anneal_tm = st.slider("Primer anneal Tm target (°C)", 50.0, 72.0, 60.0, 0.5)
    product_len_min, product_len_max = st.slider("Desired product length (bp)", 60, 240, (100, 170))
    allow_tails = st.checkbox("Enable automatic Tm tuning with non‑templated tails", True)

st.markdown("### 2) Assign **3 dye channels × 5 Tm slots**")
assign_cols = st.columns(3)
slot_map: Dict[int, Tuple[str, int, float]] = {}
slot_labels = []
slot_i = 0
for j, dye in enumerate(["FAM", "HEX", "Cy5"]):
    with assign_cols[j]:
        st.subheader(dye)
        for k in range(5):
            slot_map[slot_i] = (dye, k, grid[dye][k])
            slot_labels.append(f"{dye} • s{k+1} @ {grid[dye][k]:.1f} °C")
            slot_i += 1

assigned = []
for i in range(15):
    assigned.append(st.selectbox(f"Slot for Target {i+1}", slot_labels, index=i % 15, key=f"slot_{i}"))

# -----------------------------
# Core design routine
# -----------------------------

def design_one(seq: str, size_min: int, size_opt: int, size_max: int) -> Optional[Tuple[str, str, int]]:
    if not seq or len(seq) < size_min + 30:
        return None
    if HAVE_PRIMER3:
        params = dict(
            SEQUENCE_TEMPLATE=seq,
            PRIMER_OPT_SIZE=size_opt,
            PRIMER_MIN_SIZE=size_min,
            PRIMER_MAX_SIZE=size_max,
            PRIMER_NUM_RETURN=3,
            PRIMER_PRODUCT_SIZE_RANGE=[[product_len_min, product_len_max]],
            PRIMER_OPT_TM=anneal_tm,
            PRIMER_MIN_TM=anneal_tm-2.5,
            PRIMER_MAX_TM=anneal_tm+2.5,
            PRIMER_MAX_POLY_X=4,
            PRIMER_GC_CLAMP=1,
        )
        try:
            res = primer3.bindings.designPrimers({"SEQUENCE_TEMPLATE": seq}, params)
            for idx in range(int(res.get("PRIMER_PAIR_NUM_RETURNED", 0))):
                fwd = res[f"PRIMER_LEFT_{idx}_SEQUENCE"]
                rev = res[f"PRIMER_RIGHT_{idx}_SEQUENCE"]
                start = res[f"PRIMER_LEFT_{idx}"][0]
                end   = res[f"PRIMER_RIGHT_{idx}"][0]
                prod_len = (end + len(rev)) - start
                if product_len_min <= prod_len <= product_len_max:
                    return fwd, rev, prod_len
        except Exception:
            pass
    # Fallback (naïve): slice ends
    fwd = seq[:size_opt]
    rev = Seq(seq[-size_opt:]).reverse_complement() if HAVE_BIO else seq[-size_opt:][::-1]
    prod_len = len(seq)
    return fwd, str(rev), prod_len


def tune_product_tm(fwd: str, rev: str, template: str, slot_tm: float) -> Tuple[str, str, int, float, Tuple[str,str], List[str]]:
    issues = []
    # Build amplicon sequence (simple model):
    core = template
    if len(template) > (len(fwd)+len(rev)):
        core = template[len(fwd):-(len(rev))]
    amp = f"{fwd}{core}{rev}"
    base_tm = rough_amp_tm(amp, monovalent_mM, free_Mg_mM)
    best = (fwd, rev, len(amp), base_tm, ("",""), [*issues])
    if not allow_tails:
        return best
    # Try tails to move toward slot
    candidates = [("",""), *[(t,"") for t in GC_TAILS+AT_TAILS], *[("",t) for t in GC_TAILS+AT_TAILS]]
    def score_fn(tm2, tF, tR):
        return abs(tm2 - slot_tm) + 0.02*(len(tF)+len(tR))
    best_score = score_fn(base_tm, "", "")
    for tF, tR in candidates:
        amp2 = f"{tF}{fwd}{core}{rev}{tR}"
        tm2 = rough_amp_tm(amp2, monovalent_mM, free_Mg_mM)
        sc = score_fn(tm2, tF, tR)
        if sc < best_score:
            best_score = sc
            best = (tF+fwd, rev+tR, len(amp2), tm2, (tF,tR), [*issues])
    return best

# Dimer screen (coarse)

def has_bad_3prime_dimer(a: str, b: str) -> bool:
    # very simple 3' complementarity check (last 5 nt)
    a3 = a[-5:]
    b3 = b[-5:]
    comp = str.maketrans("ACGT", "TGCA")
    return a3 == b3.translate(comp)[::-1]

# -----------------------------
# Run design
# -----------------------------
run = st.button("Design 15‑plex now")

results: List[PrimerPair] = []
if run:
    # slot reverse map
    label_to_slot = {f"{dye} • s{idx+1} @ {grid[dye][idx]:.1f} °C": (dye, idx, grid[dye][idx]) for dye in ["FAM","HEX","Cy5"] for idx in range(5)}
    for i in range(15):
        org = organisms[i]
        seq = seqs[i]
        dye, slot_idx, slot_tm = label_to_slot[assigned[i]]
        if not seq:
            st.warning(f"Target {i+1}: missing sequence — skipped")
            continue
        designed = design_one(seq, min_size, opt_size, max_size)
        if not designed:
            st.warning(f"Target {i+1}: design failed — skipped")
            continue
        fwd, rev, prod_len = designed
        fwd2, rev2, prod_len2, prod_tm2, tails, issues = tune_product_tm(fwd, rev, seq, slot_tm)
        issues = list(issues)
        # Cross‑dimer quick check against prior primers
        for prev in results:
            if has_bad_3prime_dimer(fwd2, prev.fwd) or has_bad_3prime_dimer(rev2, prev.rev):
                issues.append(f"Potential 3' cross‑dimer with {prev.organism} ({prev.channel} s{prev.slot_index+1})")
        results.append(PrimerPair(fwd2, rev2, prod_len2, prod_tm2, dye, slot_tm, slot_idx, org or f"Target {i+1}", f"Target{i+1}", tails, issues))

if results:
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
            "Issues": "; ".join(r.issues) if r.issues else ""
        } for r in results
    ]).sort_values(["Channel","Slot"]).reset_index(drop=True)

    st.success("Design complete.")
    st.dataframe(df, use_container_width=True)

    out = io.BytesIO()
    df.to_csv(out, index=False)
    st.download_button("Download manifest (CSV)", data=out.getvalue(), file_name="autoprimer5_multiplex_manifest.csv", mime="text/csv")

st.markdown("---")
st.caption("Tip: for production, replace rough_amp_tm() with MELTING/UNAFold/uMelt backend and add full primer3 dimer/hairpin screens.\nCatalog is auto‑loaded from your original autoprimer5.py if available, preserving your pathogen/target set.")
