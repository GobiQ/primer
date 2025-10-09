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
# Catalog & sequence sources (LOCAL + optional NCBI fetch)
# -----------------------------

# Local, built‑in sample entries (you can keep or replace)
def get_organism_suggestions_with_gene_targets():
    """Get the 14 specific organisms from the image as default selections"""
    return {
        "🦠 Viruses and Viroids": [
            ("Hop Latent Viroid", "Hop Latent Viroid", {
                "Essential genes": ["RNA polymerase", "ribozyme", "hammerhead ribozyme", "viroid-specific sequences"],
                "Pathogenicity genes": ["replication factor", "movement protein", "host interaction factors"],
                "Detection targets": ["viroid RNA", "replication intermediates", "host response genes"]
            }),
            ("Lettuce Chlorosis Virus", "Lettuce Chlorosis Virus", {
                "Essential genes": ["RNA polymerase", "helicase", "protease", "replicase"],
                "Structural genes": ["coat protein", "movement protein", "nucleocapsid protein"],
                "Pathogenicity genes": ["silencing suppressor", "host range determinant", "symptom determinant"]
            }),
            ("Alfalfa Mosaic Virus", "Alfalfa Mosaic Virus", {
                "Essential genes": ["RNA polymerase", "helicase", "protease", "replicase"],
                "Structural genes": ["coat protein", "movement protein", "nucleocapsid protein"],
                "Pathogenicity genes": ["silencing suppressor", "host range determinant", "symptom determinant"]
            }),
            ("Beet Curly Top Virus", "Beet Curly Top Virus", {
                "Essential genes": ["replication protein", "capsid protein", "movement protein"],
                "Pathogenicity genes": ["C4 protein", "V2 protein", "βC1 protein", "host interaction factors"],
                "Detection targets": ["viral DNA", "replication intermediates", "host response genes"]
            }),
            ("Cannabis Cryptic Virus", "Cannabis Cryptic Virus", {
                "Essential genes": ["RNA polymerase", "helicase", "protease", "replicase"],
                "Structural genes": ["coat protein", "movement protein", "nucleocapsid protein"],
                "Pathogenicity genes": ["silencing suppressor", "host range determinant", "symptom determinant"]
            }),
            ("Tomato Ring Spot Virus", "Tomato Ring Spot Virus", {
                "Essential genes": ["RNA polymerase", "helicase", "protease", "replicase"],
                "Structural genes": ["coat protein", "movement protein", "nucleocapsid protein"],
                "Pathogenicity genes": ["silencing suppressor", "host range determinant", "symptom determinant"]
            }),
            ("Tobacco Mosaic Virus", "Tobacco Mosaic Virus", {
                "Essential genes": ["RNA polymerase", "helicase", "protease", "replicase"],
                "Structural genes": ["coat protein", "movement protein", "nucleocapsid protein"],
                "Pathogenicity genes": ["silencing suppressor", "host range determinant", "symptom determinant"]
            }),
            ("Arabis Mosaic Virus", "Arabis Mosaic Virus", {
                "Essential genes": ["RNA polymerase", "helicase", "protease", "replicase"],
                "Structural genes": ["coat protein", "movement protein", "nucleocapsid protein"],
                "Pathogenicity genes": ["silencing suppressor", "host range determinant", "symptom determinant"]
            }),
            ("Tomato Mosaic Virus", "Tomato Mosaic Virus", {
                "Essential genes": ["RNA polymerase", "helicase", "protease", "replicase"],
                "Structural genes": ["coat protein", "movement protein", "nucleocapsid protein"],
                "Pathogenicity genes": ["silencing suppressor", "host range determinant", "symptom determinant"]
            })
        ],
        
        "🍄 Fungus and Pseudofungi": [
            ("Gray mold", "Botrytis cinerea", {
                "Essential genes": ["ACT1 (actin)", "TUB2 (tubulin)", "EF1A (elongation factor)", "RPB2 (RNA polymerase)", "HSP70 (heat shock protein)"],
                "Pathogenicity genes": ["BCG1 (α-galactosidase)", "BMP1 (metalloprotease)", "BCP1 (cerato-platanin)", "BOA1 (botrydial)", "BCR1 (ABC transporter)"],
                "Cell wall degrading": ["BcPG1-6 (polygalacturonases)", "BcPME1 (pectin methylesterase)", "BcXYL1 (xylanase)", "BcCEL1 (cellulase)", "BcCUT1 (cutinase)"],
                "Secondary metabolites": ["BOT1-5 (botrydial cluster)", "DHN1 (1,8-dihydroxynaphthalene)", "PKS1-13 (polyketide synthases)", "NPS1-6 (nonribosomal peptide synthetases)"],
                "Resistance mechanisms": ["ABC1-50 (ABC transporters)", "MFS1-20 (major facilitator superfamily)", "CYP1-20 (cytochrome P450s)", "GST1-10 (glutathione S-transferases)"]
            }),
            ("Pythium root rot", "Pythium myriotylum", {
                "Essential genes": ["ACT1 (actin)", "TUB1 (tubulin)", "EF1A (elongation factor)", "RPB1 (RNA polymerase)", "COX1 (cytochrome oxidase)"],
                "Pathogenicity genes": ["PmNEP1 (necrosis-inducing protein)", "PmPG1 (polygalacturonase)", "PmPME1 (pectin methylesterase)", "PmCUT1 (cutinase)", "PmCP1 (cysteine protease)"],
                "Cell wall degrading": ["PmPG1-5 (polygalacturonases)", "PmPME1 (pectin methylesterase)", "PmXYL1 (xylanase)", "PmCEL1 (cellulase)", "PmCUT1 (cutinase)"],
                "Secondary metabolites": ["PmPKS1-5 (polyketide synthases)", "PmNPS1-3 (nonribosomal peptide synthetases)", "PmTER1-2 (terpene synthases)"],
                "Resistance mechanisms": ["PmCYP1-10 (cytochrome P450s)", "PmGST1-5 (glutathione S-transferases)", "PmABC1-20 (ABC transporters)"]
            }),
            ("Fusarium wilt", "Fusarium oxysporum", {
                "Essential genes": ["ACT1 (actin)", "TUB2 (tubulin)", "EF1A (elongation factor)", "RPB2 (RNA polymerase)", "LSU (large subunit rRNA)"],
                "Pathogenicity genes": ["SIX1-14 (secreted in xylem)", "FTF1 (transcription factor)", "FMK1 (MAPK)", "SGE1 (cutinase)", "PEL1 (pectate lyase)"],
                "Secondary metabolite genes": ["FUM1 (fumonisin biosynthesis)", "TRI5 (trichothecene biosynthesis)", "PKS4 (polyketide synthase)", "BIK1 (bikaverin)", "FUS1 (fusarin)"],
                "Cell wall genes": ["CHI1 (chitinase)", "GEL1 (β-1,3-glucanase)", "CHS1 (chitin synthase)", "FKS1 (β-1,3-glucan synthase)", "PMI1 (mannose-6-phosphate isomerase)"],
                "Resistance targets": ["CYP51 (sterol 14α-demethylase)", "SDH (succinate dehydrogenase)", "QoI (cytochrome bc1)", "MBC (β-tubulin)", "DMI (sterol biosynthesis)"]
            }),
            ("Fusarium root rot", "Fusarium solani", {
                "Essential genes": ["ACT1", "TUB2", "EF1A", "RPB2", "HSP70"],
                "Pathogenicity factors": ["FSOL1-10 (F. solani specific)", "CUT1-5 (cutinases)", "PEL1-3 (pectate lyases)", "XYL1-2 (xylanases)", "CEL1-2 (cellulases)"],
                "Secondary metabolites": ["FUM1-3 (fumonisin)", "TRI1-16 (trichothecene)", "ZEA1-2 (zearalenone)", "FUS1-5 (fusarin)", "BIK1-3 (bikaverin)"],
                "Resistance mechanisms": ["CYP51A1-B1", "SDH1-4", "ABC1-20", "MFS1-15", "GST1-10"]
            }),
            ("Powdery mildew", "Golovinomyces ambrosiae", {
                "Essential genes": ["ACT1", "TUB2", "EF1A", "RPB2", "ITS1-2"],
                "Pathogenicity genes": ["GAH1 (haustorium formation)", "GAC1 (conidiophore development)", "GAS1 (spore germination)", "GAP1 (penetration)", "GAA1 (appressorium formation)"],
                "Effectors": ["CSEP1-100 (candidate secreted effector proteins)", "GAE1-50 (G. ambrosiae effectors)", "AVR1-10 (avirulence candidates)", "HAU1-20 (haustorial expressed)"],
                "Sterol biosynthesis": ["CYP51A1", "CYP51B1", "ERG1 (squalene epoxidase)", "ERG7 (lanosterol synthase)", "ERG11 (sterol 14α-demethylase)"]
            }),
            ("Fusarium ear rot", "Fusarium proliferatum", {
                "Essential genes": ["ACT1", "TUB2", "EF1A", "RPB2", "ITS1-2"],
                "Fumonisin biosynthesis": ["FUM1-21 (fumonisin cluster)", "FUM8 (polyketide synthase)", "FUM6 (aminotransferase)", "FUM3 (C-5 oxygenase)", "FUM19 (transporter)"],
                "Pathogenicity": ["FPR1-10 (F. proliferatum specific)", "CUT1-3", "PEL1-2", "XYL1", "CEL1"],
                "Host interaction": ["HOST1-5 (host-specific)", "ADH1-3 (adhesion)", "INV1-2 (invasion)", "COL1-3 (colonization)"]
            })
        ]
    }

LOCAL_CATALOG = get_organism_suggestions_with_gene_targets()

# No hardcoded sequences - rely entirely on NCBI fetching like autoprimer5.py
LOCAL_SEQUENCES = {}

@st.cache_data
def load_catalog() -> Any:
    """Return the embedded local catalog by default (no imports required)."""
    return LOCAL_CATALOG


def clean_seq(s: str) -> str:
    return re.sub(r"[^ACGTNacgtn]", "", s).upper().replace("U", "T")

# ---- NCBI fetching (optional, like the original app) ----
try:
    from Bio import Entrez
    HAVE_ENTREZ = True
except Exception:
    HAVE_ENTREZ = False

@st.cache_data(show_spinner=False)
def ncbi_search_and_fetch(organism: str, gene: str, email: str, api_key: Optional[str], db: str = "nucleotide", retmax: int = 1) -> Optional[str]:
    """Search NCBI for a nucleotide record for organism+gene and return FASTA sequence (first hit)."""
    import urllib.parse, urllib.request
    import json

    def eutils_get(url, params):
        q = urllib.parse.urlencode(params)
        with urllib.request.urlopen(f"{url}?{q}") as resp:
            return resp.read().decode()

    # Use Biopython if available
    if HAVE_ENTREZ:
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
        handle = Entrez.esearch(db=db, term=f"{organism}[Organism] AND {gene}[Gene]", retmax=retmax)
        rec = Entrez.read(handle)
        ids = rec.get("IdList", [])
        if not ids:
            return None
        fetch = Entrez.efetch(db=db, id=ids[0], rettype="fasta", retmode="text")
        fasta = fetch.read()
        seq = "".join([line.strip() for line in fasta.splitlines() if not line.startswith(">")])
        return clean_seq(seq) if seq else None

    # Fallback to raw E-utilities (no Biopython)
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    search_xml = eutils_get(f"{base}/esearch.fcgi", {
        "db": db, "term": f"{organism}[Organism] AND {gene}[Gene]", "retmax": str(retmax),
        "retmode": "json", **({"api_key": api_key} if api_key else {})
    })
    data = json.loads(search_xml)
    ids = data.get("esearchresult", {}).get("idlist", [])
    if not ids:
        return None
    fasta = eutils_get(f"{base}/efetch.fcgi", {
        "db": db, "id": ids[0], "rettype": "fasta", "retmode": "text", **({"api_key": api_key} if api_key else {})
    })
    seq = "".join([line.strip() for line in fasta.splitlines() if not line.startswith(">")])
    return clean_seq(seq) if seq else None

# -----------------------------
# Sequence extractor (prefers local, then NCBI if enabled)
# -----------------------------

def find_seq_in_obj(obj: Any, key_hints: List[str]) -> Optional[str]:
    # direct string
    if isinstance(obj, str):
        s = clean_seq(obj)
        if len(s) >= 60:
            return s
        return None
    if isinstance(obj, dict):
        for k in key_hints:
            if k in obj and isinstance(obj[k], str):
                s = clean_seq(obj[k])
                if len(s) >= 60:
                    return s
        for v in obj.values():
            s = find_seq_in_obj(v, key_hints)
            if s:
                return s
        return None
    if isinstance(obj, (list, tuple)):
        for v in obj:
            s = find_seq_in_obj(v, key_hints)
            if s:
                return s
    return None

@st.cache_data
def flatten_catalog_with_sequences(catalog: Any, key_hints_csv: str, email: str = "", api_key: Optional[str] = None, enable_ncbi: bool = False) -> List[Dict[str,str]]:
    """Flatten catalog into entries with (organism, target, sequence).
    Order of preference: LOCAL_SEQUENCES -> embedded sequence fields -> NCBI fetch (if enabled).
    """
    key_hints = [k.strip() for k in key_hints_csv.split(',') if k.strip()]
    out: List[Dict[str,str]] = []
    missing: List[Tuple[str,str]] = []

    def push_entry(org: str, tgt: str, seq: Optional[str], path: str, label_hint: Optional[str] = None):
        label = label_hint or f"{org} — {tgt}"
        if seq:
            out.append({"label": label, "organism": org, "target": tgt, "sequence": seq, "path": path})
        else:
            missing.append((org, tgt))

    def walk(obj: Any, path: str = "$", category: Optional[str] = None, subcat: Optional[str] = None):
        if obj is None:
            return
        
        # Handle the new catalog structure: (common_name, scientific_name, gene_targets)
        if isinstance(obj, (list, tuple)) and len(obj) == 3:
            common_name, scientific_name, gene_targets = obj
            if isinstance(gene_targets, dict):
                # Create entries for each gene category
                for gene_category, genes in gene_targets.items():
                    for gene in genes:
                        # Extract gene name (remove description in parentheses)
                        gene_name = gene.split('(')[0].strip()
                        label = f"{scientific_name} — {gene_name}"
                        push_entry(scientific_name, gene_name, None, f"{path}.{gene_category}.{gene_name}", label)
                return
        
        # If this node has an explicit sequence, infer labels
        seq_here = find_seq_in_obj(obj, key_hints)
        if seq_here:
            organism = None
            target = None
            if isinstance(obj, dict):
                organism = obj.get("scientific") or obj.get("organism") or obj.get("species") or obj.get("name") or obj.get("common") or category
                target = obj.get("target") or obj.get("gene") or obj.get("marker") or obj.get("locus") or subcat or "Target"
            elif isinstance(obj, (list, tuple)) and len(obj) >= 3 and all(isinstance(x, str) for x in obj[:3]):
                organism = obj[1]
                target = obj[2]
            push_entry(organism or "Unknown", target or "Target", seq_here, path)
            return
        # Otherwise recurse; when we find a leaf dict with organism + target but no seq, try LOCAL then NCBI
        if isinstance(obj, dict):
            # Detect a leaf with organism/target names
            if any(k in obj for k in ("scientific","organism","species","name","common")) and any(k in obj for k in ("target","gene","marker","locus")):
                org = obj.get("scientific") or obj.get("organism") or obj.get("species") or obj.get("name") or obj.get("common") or category or "Unknown"
                tgt = obj.get("target") or obj.get("gene") or obj.get("marker") or obj.get("locus") or subcat or "Target"
                # 1) local override
                seq = LOCAL_SEQUENCES.get((org, tgt)) or None
                if seq:
                    seq = clean_seq(seq)
                # 2) embedded fields
                if not seq:
                    seq = find_seq_in_obj(obj, key_hints)
                # 3) NCBI fetch
                if not seq and enable_ncbi and email:
                    seq = ncbi_search_and_fetch(org, tgt, email=email, api_key=api_key)
                push_entry(org, tgt, seq, path)
                return
            for k, v in obj.items():
                if category is None:
                    walk(v, f"{path}.{k}", category=k, subcat=None)
                elif subcat is None and isinstance(v, (list, dict, tuple)):
                    walk(v, f"{path}.{k}", category=category, subcat=k)
                else:
                    walk(v, f"{path}.{k}", category=category, subcat=subcat)
        elif isinstance(obj, (list, tuple)):
            for i, v in enumerate(obj):
                walk(v, f"{path}[{i}]", category=category, subcat=subcat)

    walk(catalog)
    # Deduplicate entries
    seen = set()
    dedup = []
    for e in out:
        key = (e["organism"], e["target"], e["sequence"])
        if key in seen:
            continue
        seen.add(key)
        dedup.append(e)
    # Attach missing info for UI
    dedup.append({"_missing": missing})
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
# NCBI Connector
# -----------------------------
import time
from functools import wraps

def retry_with_backoff(max_retries=3, base_delay=0.34):
    """Decorator for retrying functions with exponential backoff"""
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            for attempt in range(max_retries):
                try:
                    return func(*args, **kwargs)
                except Exception as e:
                    if attempt == max_retries - 1:
                        raise e
                    delay = base_delay * (2 ** attempt)
                    time.sleep(delay)
            return None
        return wrapper
    return decorator

class ResilientNCBIConnector:
    """Resilient NCBI connector with exponential backoff and retry logic"""
    
    def __init__(self, email: str, api_key: Optional[str] = None):
        if HAVE_ENTREZ:
            Entrez.email = email
            Entrez.tool = "autoprimertm"
            if api_key:
                Entrez.api_key = api_key
        self.base_delay = 0.34 if not api_key else 0.1
        self.max_retries = 3
        self.timeout = 30
        self.throttle_sec = 0.34
    
    def search_sequences(self, query: str, database: str = "nucleotide", max_results: int = 100) -> List[str]:
        """Search with timeout and retry logic"""
        if not HAVE_ENTREZ:
            return []
        
        time.sleep(self.throttle_sec)
        
        try:
            handle = Entrez.esearch(
                db=database, 
                term=query, 
                retmax=max_results,
                timeout=self.timeout
            )
            search_results = Entrez.read(handle)
            handle.close()
            return search_results["IdList"]
        except Exception as e:
            st.error(f"NCBI search error: {e}")
            return []
    
    @retry_with_backoff(max_retries=3, base_delay=0.34)
    def fetch_sequence(self, seq_id: str, database: str = "nucleotide") -> Optional[str]:
        """Fetch with timeout and retry logic"""
        if not HAVE_ENTREZ:
            return None
            
        try:
            handle = Entrez.efetch(
                db=database, 
                id=seq_id, 
                rettype="fasta", 
                retmode="text",
                timeout=self.timeout
            )
            record = SeqIO.read(handle, "fasta")
            handle.close()
            return str(record.seq)
        except Exception as e:
            st.error(f"NCBI fetch error: {e}")
            return None
    
    def fetch_organism_sequences(self, organism_name: str, max_sequences: int = 10) -> List[Dict]:
        """Fetch multiple sequences for an organism with metadata using autoprimer5.py methodology"""
        if not HAVE_ENTREZ or not organism_name:
            return []
            
        try:
            # Get gene targets for this organism
            suggestions = get_organism_suggestions_with_gene_targets()
            organism_targets = None
            
            # Find matching organism and extract gene targets
            for category, organisms in suggestions.items():
                for item in organisms:
                    if len(item) == 3:
                        common_name, scientific_name, gene_targets = item
                        if (scientific_name.lower() == organism_name.lower() or 
                            organism_name.lower() in scientific_name.lower() or
                            scientific_name.lower() in organism_name.lower()):
                            organism_targets = {
                                'organism': scientific_name,
                                'common_name': common_name,
                                'gene_targets': gene_targets
                            }
                            break
                if organism_targets:
                    break
            
            sequences = []
            
            if organism_targets:
                # Fetch sequences for specific gene targets using enhanced search
                gene_targets = organism_targets['gene_targets']
                for gene_category, genes in gene_targets.items():
                    for gene in genes[:2]:  # Limit to 2 genes per category
                        gene_name = gene.split('(')[0].strip()
                        
                        # Enhanced search queries like autoprimer5.py
                        search_queries = [
                            f'"{organism_name}"[organism] AND {gene_name}[gene]',
                            f'"{organism_name}"[organism] AND {gene_name}[title]',
                            f'"{organism_name}"[organism] AND {gene_name}[All Fields]',
                            f'"{organism_name}"[organism] AND {gene_name}[mRNA]',
                            f'"{organism_name}"[organism] AND {gene_name}[CDS]'
                        ]
                        
                        seq_ids = []
                        for query in search_queries:
                            try:
                                results = self.search_sequences(query, database="nucleotide", max_results=2)
                                if results:
                                    seq_ids = results
                                    break
                            except Exception:
                                continue
                        
                        # If no nucleotide hits, try protein-to-nucleotide linking
                        if not seq_ids:
                            seq_ids = self._find_nuccore_via_protein(organism_name, [gene_name], max_results=2)
                        
                        for seq_id in seq_ids:
                            sequence = self.fetch_sequence(seq_id)
                            if sequence and len(sequence) > 100:
                                # Clean sequence
                                clean_seq = re.sub(r'[^ATGCatgc]', '', sequence.upper())
                                if len(clean_seq) >= 100:
                                    sequences.append({
                                        "id": seq_id,
                                        "organism": organism_name,
                                        "target": gene_name,
                                        "sequence": clean_seq[:10000],  # Limit to 10kb like autoprimer5.py
                                        "label": f"{organism_name} — {gene_name}"
                                    })
                                    if len(sequences) >= max_sequences:
                                        return sequences
            else:
                # Fallback: general organism search
                search_query = f'"{organism_name}"[organism]'
                seq_ids = self.search_sequences(
                    search_query, 
                    database="nucleotide", 
                    max_results=max_sequences
                )
                
                for seq_id in seq_ids[:max_sequences]:
                    sequence = self.fetch_sequence(seq_id)
                    if sequence:
                        clean_seq = re.sub(r'[^ATGCatgc]', '', sequence.upper())
                        if len(clean_seq) >= 100:
                            sequences.append({
                                "id": seq_id,
                                "organism": organism_name,
                                "target": f"sequence_{seq_id}",
                                "sequence": clean_seq[:10000],
                                "label": f"{organism_name} — sequence_{seq_id}"
                            })
            
            return sequences
            
        except Exception as e:
            st.error(f"Error fetching sequences: {e}")
            return []
    
    def _find_nuccore_via_protein(self, organism_name: str, gene_variations: List[str], max_results: int = 5) -> List[str]:
        """Find nucleotide sequences via protein database linking (from autoprimer5.py)"""
        try:
            ids = []
            for gene_var in gene_variations:
                q = f'"{organism_name}"[Organism] AND ({gene_var}[All Fields])'
                handle = Entrez.esearch(db="protein", term=q, retmax=max_results)
                rec = Entrez.read(handle)
                handle.close()
                if not rec.get("IdList"):
                    continue
                for pid in rec["IdList"]:
                    link = Entrez.elink(dbfrom="protein", db="nuccore", id=pid)
                    link_rec = Entrez.read(link)
                    link.close()
                    for lset in link_rec:
                        for linksetdb in lset.get("LinkSetDb", []):
                            if linksetdb.get("DbTo") == "nuccore":
                                ids.extend([x["Id"] for x in linksetdb.get("Link", [])])
            return ids[:max_results]
        except Exception:
            return []

# -----------------------------
# UI
# -----------------------------
st.set_page_config(page_title="AutoPrimer5 — Multiplex (multiplexprimers.py)", layout="wide")
st.title("AutoPrimer5 — Multiplex (15 targets • 3 channels × 5 Tm slots)")

with st.sidebar:
    st.header("NCBI Configuration")
    ncbi_email = st.text_input("NCBI Email", value="your.email@example.com", help="Required for NCBI API access")
    ncbi_api_key = st.text_input("NCBI API Key (optional)", type="password", help="Optional API key for higher rate limits")
    
    st.header("Organism & Target Selection")
    
    # Organism selection buttons
    st.subheader("Quick Select Organisms")
    suggestions = get_organism_suggestions_with_gene_targets()
    
    # Create organism buttons organized by category
    for category, organisms in suggestions.items():
        with st.expander(f"{category}", expanded=True):  # Expanded by default to show all organisms
            cols = st.columns(3)
            for i, organism_item in enumerate(organisms):
                if len(organism_item) == 3:
                    common_name, scientific_name, gene_targets = organism_item
                    col_idx = i % 3
                    with cols[col_idx]:
                        if st.button(
                            f"{common_name}", 
                            key=f"org_btn_{scientific_name.replace(' ', '_')}", 
                            help=f"Select {scientific_name}",
                            use_container_width=True
                        ):
                            st.session_state['selected_organism_name'] = scientific_name
                            st.rerun()
    
    # Organism name input
    default_organism = st.session_state.get('selected_organism_name', '')
    organism_name = st.text_input("Organism Name", value=default_organism, placeholder="e.g., Hop Latent Viroid, Botrytis cinerea, Fusarium oxysporum", help="Enter the scientific name of the organism")
    
    # Clear the selected organism after setting it
    if 'selected_organism_name' in st.session_state:
        del st.session_state['selected_organism_name']
    
    # Add button to fetch sequences from NCBI
    if organism_name and ncbi_email and ncbi_email != "your.email@example.com":
        if st.button("Fetch Sequences from NCBI", type="primary"):
            with st.spinner(f"Fetching sequences for {organism_name}..."):
                try:
                    ncbi = ResilientNCBIConnector(ncbi_email, ncbi_api_key if ncbi_api_key else None)
                    ncbi_sequences = ncbi.fetch_organism_sequences(organism_name, max_sequences=15)
                    if ncbi_sequences:
                        st.success(f"Found {len(ncbi_sequences)} sequences for {organism_name}")
                        # Store in session state for use in the main app
                        st.session_state['ncbi_sequences'] = ncbi_sequences
                    else:
                        st.warning(f"No sequences found for {organism_name}")
                except Exception as e:
                    st.error(f"Error fetching sequences: {e}")
    
    st.header("Catalog & Conditions")
    seq_key_hints = st.text_input("Sequence field key(s) (comma-separated)", value="sequence,target_sequence,amplicon,region,locus_seq,seq")
    catalog_obj = load_catalog()
    
    # Use NCBI sequences if available, otherwise use local catalog
    if 'ncbi_sequences' in st.session_state and st.session_state['ncbi_sequences']:
        entries = st.session_state['ncbi_sequences']
    else:
        entries = flatten_catalog_with_sequences(catalog_obj, seq_key_hints, ncbi_email, ncbi_api_key, enable_ncbi=bool(organism_name))
    st.metric("Catalog entries (raw)", 0 if catalog_obj is None else (len(catalog_obj) if hasattr(catalog_obj, "__len__") else "?"))
    st.metric("Entries with sequences", len(entries))
    if len(entries) == 0:
        if organism_name:
            st.warning("No sequences found. Try clicking 'Fetch Sequences from NCBI' to get sequences for your organism.")
        else:
            st.warning("No sequences found in catalog. Enter an organism name and fetch sequences from NCBI, or verify the local catalog has sequences.")
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
labels = [e.get("label", "Unknown") for e in entries] if entries else ["<none>"]
selected_labels = []
for i in range(15):
    idx = min(i, max(0, len(labels)-1))
    selected_labels.append(st.selectbox(f"Select target {i+1}", labels, index=idx, key=f"sel_{i}"))

selected = [next((e for e in entries if e.get("label")==lbl), None) for lbl in selected_labels]

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
