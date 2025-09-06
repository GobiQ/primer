def get_related_organisms(target_organism):
    """Get related organisms for specificity testing"""
    organism_lower = target_organism.lower()
    
    # Comprehensive organism relationships including viruses
    related_map = {
        # Fungi
        'fusarium': ['Trichoderma harzianum', 'Aspergillus niger', 'Penicillium chrysogenum', 'Alternaria alternata'],
        'botrytis': ['Sclerotinia sclerotiorum', 'Alternaria alternata', 'Fusarium oxysporum', 'Rhizoctonia solani'],
        'pythium': ['Phytophthora infestans', 'Fusarium oxysporum', 'Rhizoctonia solani'],
        'alternaria': ['Stemphylium vesicarium', 'Botrytis cinerea', 'Fusarium oxysporum'],
        'rhizoctonia': ['Fusarium oxysporum', 'Pythium ultimum', 'Sclerotinia sclerotiorum'],
        'aspergillus': ['Penicillium chrysogenum', 'Trichoderma harzianum', 'Fusarium oxysporum'],
        'penicillium': ['Aspergillus niger', 'Trichoderma harzianum', 'Fusarium oxysporum'],
        
        # Arthropods
        'tetranychus': ['Panonychus ulmi', 'Oligonychus ilicis', 'Aculops lycopersici'],
        'bemisia': ['Trialeurodes vaporariorum', 'Aphis gossypii', 'Myzus persicae'],
        'thrips': ['Frankliniella occidentalis', 'Thrips palmi', 'Scirtothrips dorsalis'],
        
        # Viruses - KEY ADDITION
        'coronavirus': ['Influenza A virus', 'Rhinovirus', 'Adenovirus', 'Respiratory syncytial virus'],
        'sars': ['MERS-CoV', 'Human coronavirus 229E', 'Human coronavirus OC43', 'Influenza A virus'],
        'influenza': ['Coronavirus', 'Rhinovirus', 'Parainfluenza virus', 'Respiratory syncytial virus'],
        'tobacco mosaic': ['Tomato mosaic virus', 'Cucumber mosaic virus', 'Potato virus X'],
        'tomato': ['Tobacco mosaic virus', 'Cucumber mosaic virus', 'Potato virus Y'],
        
        # Bacteria
        'ralstonia': ['Pseudomonas syringae', 'Xanthomonas campestris', 'Erwinia amylovora'],
        'erwinia': ['Pseudomonas syringae', 'Xanthomonas campestris', 'Ralstonia solanacearum'],
        'agrobacterium': ['Rhizobium leguminosarum', 'Sinorhizobium meliloti', 'Pseudomonas syringae']
    }
    
    # Find matching genus
    for genus, related in related_map.items():
        if genus in organism_lower:
            return related
    
    # Default related organisms for general testing
    return ['Aspergillus niger', 'Penicillium chrysogenum', 'Trichoderma harzianum', 'Fusarium oxysporum']

def check_session_state_validity():
    """Check if session state has valid data"""
    has_primers = bool(st.session_state.get('primers_designed'))
    has_sequence = bool(st.session_state.get('current_sequence'))
    has_seq_info = bool(st.session_state.get('sequence_info'))
    
    return {
        'has_primers': has_primers,
        'has_sequence': has_sequence,
        'has_seq_info': has_seq_info,
        'primer_count': len(st.session_state.get('primers_designed', [])),
        'sequence_length': len(st.session_state.get('current_sequence', ''))
    }

def get_organism_suggestions():
    """Get agricultural pest and pathogen suggestions organized by category"""
    return {
        "🍄 Fungal Pathogens": {
            "Fusarium species": [
                ("Fusarium wilt", "Fusarium oxysporum"),
                ("Fusarium head blight", "Fusarium graminearum"),
                ("Fusarium crown rot", "Fusarium culmorum"),
                ("Fusarium root rot", "Fusarium solani"),
                ("Fusarium ear rot", "Fusarium proliferatum")
            ],
            "Other fungi": [
                ("Gray mold", "Botrytis cinerea"),
                ("White mold", "Sclerotinia sclerotiorum"),
                ("Late blight", "Phytophthora infestans"),
                ("Powdery mildew", "Erysiphe necator"),
                ("Rust disease", "Puccinia graminis"),
                ("Smut disease", "Ustilago maydis")
            ]
        },
        "🐛 Insect Pests": {
            "Mites": [
                ("Two-spotted spider mite", "Tetranychus urticae"),
                ("European red mite", "Panonychus ulmi"),
                ("Broad mite", "Polyphagotarsonemus latus"),
                ("Russet mite", "Aculops lycopersici")
            ],
            "Sucking insects": [
                ("Silverleaf whitefly", "Bemisia tabaci"),
                ("Greenhouse whitefly", "Trialeurodes vaporariorum"),
                ("Green peach aphid", "Myzus persicae"),
                ("Cotton aphid", "Aphis gossypii"),
                ("Asian citrus psyllid", "Diaphorina citri")
            ],
            "Thrips": [
                ("Western flower thrips", "Frankliniella occidentalis"),
                ("Onion thrips", "Thrips tabaci")
            ]
        },
        "🦠 Oomycetes": {
            "Water molds": [
                ("Pythium root rot", "Pythium"),
                ("Pythium damping-off", "Pythium myriotylum")
            ]
        },
        "🦠 Bacterial Pathogens": {
            "Common bacteria": [
                ("Bacterial wilt", "Ralstonia solanacearum"),
                ("Fire blight", "Erwinia amylovora"),
                ("Crown gall", "Agrobacterium tumefaciens")
            ]
        },
        "🦠 Plant Viruses": {
            "Tobamoviruses": [
                ("Tobacco mosaic virus", "Tobacco mosaic virus"),
                ("Tomato mosaic virus", "Tomato mosaic virus")
            ],
            "Other viruses": [
                ("Beet curly top virus", "Beet curly top virus"),
                ("Arabis mosaic virus", "Arabis mosaic virus"),
                ("Alfalfa mosaic virus", "Alfalfa mosaic virus"),
                ("Cannabis cryptic virus", "Cannabis cryptic virus"),
                ("Lettuce chlorosis virus", "Lettuce chlorosis virus")
            ]
        },
        "🧬 Viroids": {
            "RNA pathogens": [
                ("Hop latent viroid", "Hop latent viroid")
            ]
        }
    }

#!/usr/bin/env python3
"""
Streamlit Web Application for Automated Primer Design - COMPLETE FIXED VERSION
============================================================================

Complete version with T7 dsRNA functionality and all bug fixes.

Key features:
1. Fixed session state management
2. T7 promoter addition for dsRNA production
3. Comprehensive analysis and export options
4. Agricultural pest/pathogen focus

Installation:
pip install streamlit biopython primer3-py requests pandas openpyxl plotly

Run with:
streamlit run autoprimer.py

Author: Automated Primer Design System
"""

import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import requests
import json
import time
import io
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils.MeltingTemp import Tm_NN
from Bio.Align import PairwiseAligner
from Bio.Blast import NCBIWWW, NCBIXML
import primer3
import re
from pathlib import Path
import base64
import numpy as np
from collections import defaultdict
import warnings
warnings.filterwarnings("ignore")

# Configure Streamlit page
st.set_page_config(
    page_title="Automated Primer Designer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

@dataclass
class PrimerPair:
    """Class to store primer pair information"""
    forward_seq: str
    reverse_seq: str
    forward_tm: float
    reverse_tm: float
    product_size: int
    gc_content_f: float
    gc_content_r: float
    forward_start: int
    reverse_start: int
    penalty: float = 0.0

class NCBIConnector:
    """Handles all NCBI database connections and queries"""
    
    def __init__(self, email: str, api_key: Optional[str] = None):
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
        self.rate_limit_delay = 0.1 if api_key else 0.34
    
    def search_sequences(self, query: str, database: str = "nucleotide", 
                        max_results: int = 100) -> List[str]:
        try:
            time.sleep(self.rate_limit_delay)
            handle = Entrez.esearch(db=database, term=query, retmax=max_results)
            search_results = Entrez.read(handle)
            handle.close()
            return search_results["IdList"]
        except Exception as e:
            st.error(f"Error searching NCBI: {e}")
            return []
    
    def fetch_sequence(self, seq_id: str, database: str = "nucleotide") -> Optional[str]:
        try:
            time.sleep(self.rate_limit_delay)
            handle = Entrez.efetch(db=database, id=seq_id, rettype="fasta", retmode="text")
            record = SeqIO.read(handle, "fasta")
            handle.close()
            return str(record.seq)
        except Exception as e:
            st.error(f"Error fetching sequence {seq_id}: {e}")
            return None
    
    def fetch_sequence_info(self, seq_id: str, database: str = "nucleotide") -> Dict:
        try:
            time.sleep(self.rate_limit_delay)
            handle = Entrez.efetch(db=database, id=seq_id, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()
            
            return {
                "id": record.id,
                "description": record.description,
                "length": len(record.seq),
                "organism": record.annotations.get("organism", "Unknown"),
                "sequence": str(record.seq)
            }
        except Exception as e:
            return {}

class PrimerDesigner:
    """Main primer design class with T7 dsRNA functionality"""
    
    def __init__(self):
        self.default_params = {
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_MIN_SIZE': 18,
            'PRIMER_MAX_SIZE': 25,
            'PRIMER_OPT_TM': 60.0,
            'PRIMER_MIN_TM': 57.0,
            'PRIMER_MAX_TM': 63.0,
            'PRIMER_MIN_GC': 40.0,
            'PRIMER_MAX_GC': 60.0,
            'PRIMER_MAX_POLY_X': 4,
            'PRIMER_SALT_MONOVALENT': 50.0,
            'PRIMER_DNA_CONC': 50.0,
            'PRIMER_MAX_NS_ACCEPTED': 0,
            'PRIMER_MAX_SELF_ANY': 12,
            'PRIMER_MAX_SELF_END': 8,
            'PRIMER_PAIR_MAX_COMPL_ANY': 12,
            'PRIMER_PAIR_MAX_COMPL_END': 8,
            'PRIMER_PRODUCT_SIZE_RANGE': [[75, 300], [300, 600], [600, 1000]],
            'PRIMER_NUM_RETURN': 20
        }
        
        # T7 promoter sequence for dsRNA production
        self.t7_promoter = "TAATACGACTCACTATAGGG"
    
    def calculate_gc_content(self, sequence: str) -> float:
        gc_count = sequence.upper().count('G') + sequence.upper().count('C')
        return (gc_count / len(sequence)) * 100 if sequence else 0
    
    def reverse_complement(self, sequence: str) -> str:
        """Generate reverse complement of DNA sequence"""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        return "".join(complement.get(base, base) for base in reversed(sequence.upper()))
    
    def design_primers(self, sequence: str, target_region: Optional[Tuple[int, int]] = None,
                      custom_params: Optional[Dict] = None, add_t7_promoter: bool = False) -> List[PrimerPair]:
        """Design primers with optional T7 promoter for dsRNA production"""
        params = self.default_params.copy()
        if custom_params:
            params.update(custom_params)
        
        seq_args = {
            'SEQUENCE_ID': 'target',
            'SEQUENCE_TEMPLATE': sequence.upper()
        }
        
        if target_region:
            start, end = target_region
            seq_args['SEQUENCE_TARGET'] = [start, end - start]
        
        try:
            primer_results = primer3.bindings.designPrimers(seq_args, params)
            primers = []
            
            num_pairs = primer_results.get('PRIMER_PAIR_NUM_RETURNED', 0)
            
            for i in range(num_pairs):
                forward_seq = primer_results[f'PRIMER_LEFT_{i}_SEQUENCE']
                reverse_seq = primer_results[f'PRIMER_RIGHT_{i}_SEQUENCE']
                forward_start = primer_results[f'PRIMER_LEFT_{i}'][0]
                reverse_start = primer_results[f'PRIMER_RIGHT_{i}'][0]
                product_size = primer_results[f'PRIMER_PAIR_{i}_PRODUCT_SIZE']
                penalty = primer_results[f'PRIMER_PAIR_{i}_PENALTY']
                
                # Add T7 promoter if requested (for dsRNA production)
                if add_t7_promoter:
                    # Add T7 to both forward and reverse primers for bidirectional transcription
                    forward_seq_t7 = self.t7_promoter + forward_seq
                    reverse_seq_t7 = self.t7_promoter + reverse_seq
                    
                    # Calculate Tm for core primer sequence (without T7)
                    forward_tm = Tm_NN(forward_seq)
                    reverse_tm = Tm_NN(reverse_seq)
                    
                    # Store both original and T7-modified sequences
                    primer_pair = PrimerPair(
                        forward_seq=forward_seq_t7,  # T7 + primer
                        reverse_seq=reverse_seq_t7,  # T7 + primer  
                        forward_tm=forward_tm,       # Tm of core primer
                        reverse_tm=reverse_tm,       # Tm of core primer
                        product_size=product_size,
                        gc_content_f=self.calculate_gc_content(forward_seq),  # GC of core primer
                        gc_content_r=self.calculate_gc_content(reverse_seq),  # GC of core primer
                        forward_start=forward_start,
                        reverse_start=reverse_start,
                        penalty=penalty
                    )
                    
                    # Store additional T7 information
                    primer_pair.core_forward_seq = forward_seq
                    primer_pair.core_reverse_seq = reverse_seq
                    primer_pair.has_t7_promoter = True
                    primer_pair.t7_promoter_seq = self.t7_promoter
                    
                else:
                    # Standard primers without T7
                    forward_tm = Tm_NN(forward_seq)
                    reverse_tm = Tm_NN(reverse_seq)
                    
                    primer_pair = PrimerPair(
                        forward_seq=forward_seq,
                        reverse_seq=reverse_seq,
                        forward_tm=forward_tm,
                        reverse_tm=reverse_tm,
                        product_size=product_size,
                        gc_content_f=self.calculate_gc_content(forward_seq),
                        gc_content_r=self.calculate_gc_content(reverse_seq),
                        forward_start=forward_start,
                        reverse_start=reverse_start,
                        penalty=penalty
                    )
                    
                    primer_pair.has_t7_promoter = False
                
                primers.append(primer_pair)
            
            return primers
            
        except Exception as e:
            st.error(f"Error in primer design: {e}")
            return []
    
    def calculate_dsrna_properties(self, primer_pair: PrimerPair, sequence: str) -> Dict:
        """Calculate properties relevant for dsRNA production"""
        if not hasattr(primer_pair, 'has_t7_promoter') or not primer_pair.has_t7_promoter:
            return {}
        
        try:
            # Extract the target region that will be transcribed
            target_start = primer_pair.forward_start
            target_end = primer_pair.reverse_start
            target_sequence = sequence[target_start:target_end + 1]
            
            # Calculate dsRNA properties
            dsrna_length = len(target_sequence)
            gc_content = self.calculate_gc_content(target_sequence)
            
            # Check for optimal dsRNA characteristics
            optimal_length = 100 <= dsrna_length <= 500  # Optimal for RNAi
            moderate_gc = 40 <= gc_content <= 60  # Avoid extreme GC content
            
            # Calculate T7 transcription efficiency factors
            # T7 prefers certain nucleotides at +1 position (G is best)
            transcription_start = target_sequence[0] if target_sequence else 'N'
            t7_efficiency = "High" if transcription_start == 'G' else "Moderate" if transcription_start in ['A', 'C'] else "Low"
            
            return {
                'dsrna_length': dsrna_length,
                'dsrna_gc_content': gc_content,
                'target_sequence': target_sequence[:100] + '...' if len(target_sequence) > 100 else target_sequence,
                'optimal_length': optimal_length,
                'moderate_gc': moderate_gc,
                'transcription_efficiency': t7_efficiency,
                'transcription_start': transcription_start,
                'estimated_yield': "High" if optimal_length and moderate_gc else "Moderate" if optimal_length or moderate_gc else "Low"
            }
        
        except Exception as e:
            return {'error': str(e)}

class ConservationAnalyzer:
    """Analyzes sequence conservation and specificity across multiple sequences"""
    
    def __init__(self, ncbi_connector):
        self.ncbi = ncbi_connector
    
    def analyze_multiple_sequences(self, sequences, min_conservation=0.8, window_size=200, step_size=50):
        """Find conserved regions across multiple sequences"""
        if len(sequences) < 2:
            return []
        
        # Find minimum sequence length
        min_length = min(len(seq) for seq in sequences)
        if min_length < window_size:
            return []
        
        conserved_regions = []
        
        # Sliding window analysis
        for start in range(0, min_length - window_size + 1, step_size):
            end = start + window_size
            
            # Extract windows from all sequences
            windows = [seq[start:end] for seq in sequences]
            
            # Calculate conservation score
            conservation_score = self._calculate_conservation_score(windows)
            
            if conservation_score >= min_conservation:
                # Get consensus sequence
                consensus = self._get_consensus_sequence(windows)
                
                conserved_regions.append({
                    'start': start,
                    'end': end,
                    'length': window_size,
                    'conservation_score': conservation_score,
                    'consensus_sequence': consensus,
                    'sequence_count': len(sequences)
                })
        
        # Merge overlapping regions
        merged_regions = self._merge_overlapping_regions(conserved_regions)
        
        return merged_regions
    
    def _calculate_conservation_score(self, windows):
        """Calculate conservation score for a set of sequence windows"""
        if not windows:
            return 0.0
        
        total_positions = len(windows[0])
        conserved_positions = 0
        
        for pos in range(total_positions):
            # Get nucleotides at this position across all sequences
            nucleotides = [window[pos].upper() for window in windows if pos < len(window)]
            
            if nucleotides:
                # Most common nucleotide
                from collections import Counter
                most_common = Counter(nucleotides).most_common(1)[0]
                frequency = most_common[1] / len(nucleotides)
                
                # Consider position conserved if >80% have same nucleotide
                if frequency >= 0.8:
                    conserved_positions += 1
        
        return conserved_positions / total_positions if total_positions > 0 else 0.0
    
    def _get_consensus_sequence(self, windows):
        """Generate consensus sequence from multiple windows"""
        if not windows:
            return ""
        
        consensus = []
        length = len(windows[0])
        
        for pos in range(length):
            nucleotides = [window[pos].upper() for window in windows if pos < len(window)]
            
            if nucleotides:
                from collections import Counter
                most_common = Counter(nucleotides).most_common(1)[0][0]
                consensus.append(most_common)
            else:
                consensus.append('N')
        
        return ''.join(consensus)
    
    def _merge_overlapping_regions(self, regions, min_gap=30):
        """Merge overlapping conserved regions"""
        if not regions:
            return []
        
        # Sort by start position
        sorted_regions = sorted(regions, key=lambda x: x['start'])
        merged = []
        
        current = sorted_regions[0].copy()
        
        for next_region in sorted_regions[1:]:
            # If regions overlap or are close
            if next_region['start'] - current['end'] <= min_gap:
                # Merge regions
                current['end'] = next_region['end']
                current['length'] = current['end'] - current['start']
                current['conservation_score'] = min(
                    current['conservation_score'], 
                    next_region['conservation_score']
                )
            else:
                merged.append(current)
                current = next_region.copy()
        
        merged.append(current)
        return merged
    
    def test_specificity(self, target_sequence, comparison_organisms, max_similarity=0.7):
        """Test sequence specificity against other organisms"""
        specificity_results = {}
        
        for organism in comparison_organisms:
            try:
                # Search for sequences from comparison organism
                query = f'"{organism}"[organism]'
                seq_ids = self.ncbi.search_sequences(query, max_results=5)
                
                if seq_ids:
                    similarities = []
                    
                    for seq_id in seq_ids[:3]:  # Test against top 3 sequences
                        comparison_seq = self.ncbi.fetch_sequence(seq_id)
                        if comparison_seq:
                            similarity = self._calculate_sequence_similarity(
                                target_sequence, comparison_seq
                            )
                            similarities.append(similarity)
                    
                    if similarities:
                        max_similarity_found = max(similarities)
                        is_specific = max_similarity_found < max_similarity
                        
                        specificity_results[organism] = {
                            'max_similarity': max_similarity_found,
                            'is_specific': is_specific,
                            'sequences_tested': len(similarities)
                        }
                
            except Exception as e:
                specificity_results[organism] = {
                    'error': str(e),
                    'max_similarity': 0.0,
                    'is_specific': True,
                    'sequences_tested': 0
                }
        
        return specificity_results
    
    def _calculate_sequence_similarity(self, seq1, seq2):
        """Calculate best local similarity between two sequences"""
        if not seq1 or not seq2:
            return 0.0
        
        # Simple sliding window approach for local similarity
        best_similarity = 0.0
        window_size = min(len(seq1), 100)  # Use first 100 bp of target
        
        if len(seq2) < window_size:
            return 0.0
        
        target_window = seq1[:window_size].upper()
        
        for i in range(len(seq2) - window_size + 1):
            comparison_window = seq2[i:i + window_size].upper()
            
            # Calculate identity
            matches = sum(1 for a, b in zip(target_window, comparison_window) if a == b)
            similarity = matches / window_size
            
            best_similarity = max(best_similarity, similarity)
        
        return best_similarity

class SequenceManager:
    """Manages sequence selection and analysis workflow"""
    
    def __init__(self, ncbi_connector):
        self.ncbi = ncbi_connector
    
    def fetch_organism_sequences(self, organism_name, max_sequences=10):
        """Fetch multiple sequences for an organism with metadata"""
        try:
            search_query = f'"{organism_name}"[organism]'
            seq_ids = self.ncbi.search_sequences(
                search_query, 
                database="nucleotide", 
                max_results=max_sequences
            )
            
            sequences = []
            for seq_id in seq_ids:
                seq_info = self.ncbi.fetch_sequence_info(seq_id)
                if seq_info and seq_info.get('sequence'):
                    # Clean sequence
                    clean_seq = re.sub(r'[^ATGCatgc]', '', seq_info['sequence'].upper())
                    
                    if len(clean_seq) >= 100:  # Minimum length filter
                        sequences.append({
                            'id': seq_id,
                            'description': seq_info.get('description', ''),
                            'organism': seq_info.get('organism', ''),
                            'length': len(clean_seq),
                            'sequence': clean_seq[:10000],  # Limit to 10kb for performance
                            'full_length': len(clean_seq)
                        })
            
            return sequences
            
        except Exception as e:
            st.error(f"Error fetching sequences: {e}")
            return []

# Streamlit App Functions
def init_session_state():
    """Initialize session state variables"""
    session_vars = {
        'primers_designed': [],
        'current_sequence': "",
        'sequence_info': {},
        'search_results': None,
        'database_used': None,
        'comprehensive_analysis_results': None,
        't7_results': None,
        't7_dsrna_enabled': False,
        't7_settings': {},
        # ADD THESE NEW VARIABLES:
        'found_sequences': [],
        'target_organism': '',
        'conserved_regions': [],
        'conservation_sequences': [],
        'specificity_results': {},
        'analysis_metadata': {}
    }
    
    for var, default_value in session_vars.items():
        if var not in st.session_state:
            st.session_state[var] = default_value

def debug_session_state():
    """Debug function to show session state"""
    with st.expander("🔍 Debug: Session State"):
        st.write("**Session State Variables:**")
        for key, value in st.session_state.items():
            if key.startswith(('primers', 'sequence', 'current', 'search', 'database', 'comprehensive', 't7')):
                if isinstance(value, list):
                    st.write(f"- {key}: {len(value)} items")
                elif isinstance(value, str):
                    st.write(f"- {key}: {len(value)} characters")
                elif isinstance(value, dict):
                    st.write(f"- {key}: {len(value)} keys")
                else:
                    st.write(f"- {key}: {type(value)} - {value}")

def create_primer_visualization(primers: List[PrimerPair]):
    """Create interactive visualizations for primer pairs"""
    if not primers:
        return None
    
    try:
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=('Melting Temperatures', 'GC Content Distribution', 
                           'Product Sizes', 'Penalty Scores'),
            specs=[[{"secondary_y": False}, {"secondary_y": False}],
                   [{"secondary_y": False}, {"secondary_y": False}]]
        )
        
        primer_nums = list(range(1, len(primers) + 1))
        forward_tms = [p.forward_tm for p in primers]
        reverse_tms = [p.reverse_tm for p in primers]
        forward_gcs = [p.gc_content_f for p in primers]
        reverse_gcs = [p.gc_content_r for p in primers]
        product_sizes = [p.product_size for p in primers]
        penalties = [p.penalty for p in primers]
        
        fig.add_trace(
            go.Scatter(x=primer_nums, y=forward_tms, name='Forward Tm', 
                      line=dict(color='blue'), mode='lines+markers'),
            row=1, col=1
        )
        fig.add_trace(
            go.Scatter(x=primer_nums, y=reverse_tms, name='Reverse Tm', 
                      line=dict(color='red'), mode='lines+markers'),
            row=1, col=1
        )
        
        fig.add_trace(
            go.Bar(x=primer_nums, y=forward_gcs, name='Forward GC%', 
                   marker_color='lightblue', opacity=0.7),
            row=1, col=2
        )
        fig.add_trace(
            go.Bar(x=primer_nums, y=reverse_gcs, name='Reverse GC%', 
                   marker_color='lightcoral', opacity=0.7),
            row=1, col=2
        )
        
        fig.add_trace(
            go.Scatter(x=primer_nums, y=product_sizes, name='Product Size', 
                      line=dict(color='green'), mode='lines+markers'),
            row=2, col=1
        )
        
        fig.add_trace(
            go.Bar(x=primer_nums, y=penalties, name='Penalty Score', 
                   marker_color='orange'),
            row=2, col=2
        )
        
        fig.update_layout(height=600, showlegend=True, title_text="Primer Pair Analysis")
        fig.update_xaxes(title_text="Primer Pair Number")
        fig.update_yaxes(title_text="Temperature (°C)", row=1, col=1)
        fig.update_yaxes(title_text="GC Content (%)", row=1, col=2)
        fig.update_yaxes(title_text="Product Size (bp)", row=2, col=1)
        fig.update_yaxes(title_text="Penalty Score", row=2, col=2)
        
        return fig
    except Exception as e:
        st.error(f"Error creating visualization: {e}")
        return None

def create_sequence_diagram(sequence: str, primers: List[PrimerPair], selected_primer: int = 0):
    """Create a sequence diagram showing primer binding sites"""
    if not primers or selected_primer >= len(primers) or not sequence:
        return None
    
    try:
        primer = primers[selected_primer]
        seq_len = len(sequence)
        
        fig = go.Figure()
        
        fig.add_shape(
            type="rect",
            x0=0, y0=0.4, x1=seq_len, y1=0.6,
            fillcolor="lightgray",
            line=dict(color="black", width=1),
        )
        
        # Determine primer length for visualization
        if hasattr(primer, 'has_t7_promoter') and primer.has_t7_promoter:
            forward_len = len(primer.core_forward_seq)
            reverse_len = len(primer.core_reverse_seq)
        else:
            forward_len = len(primer.forward_seq)
            reverse_len = len(primer.reverse_seq)
        
        fig.add_shape(
            type="rect",
            x0=primer.forward_start, y0=0.6, 
            x1=primer.forward_start + forward_len, y1=0.8,
            fillcolor="blue",
            line=dict(color="darkblue", width=2),
        )
        
        fig.add_shape(
            type="rect",
            x0=primer.reverse_start - reverse_len + 1, y0=0.2,
            x1=primer.reverse_start + 1, y1=0.4,
            fillcolor="red",
            line=dict(color="darkred", width=2),
        )
        
        fig.add_annotation(
            x=primer.forward_start + forward_len/2, y=0.7,
            text=f"Forward<br>Tm: {primer.forward_tm:.1f}°C",
            showarrow=True, arrowhead=2, arrowcolor="blue"
        )
        
        fig.add_annotation(
            x=primer.reverse_start - reverse_len/2, y=0.3,
            text=f"Reverse<br>Tm: {primer.reverse_tm:.1f}°C",
            showarrow=True, arrowhead=2, arrowcolor="red"
        )
        
        fig.update_layout(
            title=f"Primer Binding Sites - Pair {selected_primer + 1}",
            xaxis_title="Sequence Position (bp)",
            yaxis=dict(range=[0, 1], showticklabels=False),
            height=300,
            showlegend=False
        )
        
        return fig
    except Exception as e:
        st.error(f"Error creating sequence diagram: {e}")
        return None

def export_to_excel(primers: List[PrimerPair]) -> bytes:
    """Export primer results to Excel format"""
    try:
        data = []
        for i, primer in enumerate(primers):
            if hasattr(primer, 'has_t7_promoter') and primer.has_t7_promoter:
                data.append({
                    'Primer_Pair': i + 1,
                    'Forward_T7_Sequence': primer.forward_seq,
                    'Reverse_T7_Sequence': primer.reverse_seq,
                    'Forward_Core_Sequence': primer.core_forward_seq,
                    'Reverse_Core_Sequence': primer.core_reverse_seq,
                    'Core_Forward_Tm': round(primer.forward_tm, 2),
                    'Core_Reverse_Tm': round(primer.reverse_tm, 2),
                    'dsRNA_Size': primer.product_size,
                    'Core_Forward_GC%': round(primer.gc_content_f, 2),
                    'Core_Reverse_GC%': round(primer.gc_content_r, 2),
                    'Forward_Start': primer.forward_start,
                    'Reverse_Start': primer.reverse_start,
                    'Penalty_Score': round(primer.penalty, 4),
                    'T7_Promoter': primer.t7_promoter_seq,
                    'Primer_Type': 'T7_dsRNA'
                })
            else:
                data.append({
                    'Primer_Pair': i + 1,
                    'Forward_Sequence': primer.forward_seq,
                    'Reverse_Sequence': primer.reverse_seq,
                    'Forward_Tm': round(primer.forward_tm, 2),
                    'Reverse_Tm': round(primer.reverse_tm, 2),
                    'Product_Size': primer.product_size,
                    'Forward_GC%': round(primer.gc_content_f, 2),
                    'Reverse_GC%': round(primer.gc_content_r, 2),
                    'Forward_Start': primer.forward_start,
                    'Reverse_Start': primer.reverse_start,
                    'Penalty_Score': round(primer.penalty, 4),
                    'Primer_Type': 'Standard'
                })
        
        df = pd.DataFrame(data)
        output = io.BytesIO()
        with pd.ExcelWriter(output, engine='openpyxl') as writer:
            df.to_excel(writer, sheet_name='Primer_Results', index=False)
        
        return output.getvalue()
    except Exception as e:
        st.error(f"Error exporting to Excel: {e}")
        return b""


def check_session_state_validity():
    """Check if session state has valid data"""
    has_primers = bool(st.session_state.get('primers_designed'))
    has_sequence = bool(st.session_state.get('current_sequence'))
    has_seq_info = bool(st.session_state.get('sequence_info'))
    
    return {
        'has_primers': has_primers,
        'has_sequence': has_sequence, 
        'has_seq_info': has_seq_info,
        'primer_count': len(st.session_state.get('primers_designed', [])),
        'sequence_length': len(st.session_state.get('current_sequence', ''))
    }

def main():
    """Main Streamlit application"""
    
    init_session_state()
    
    st.title("🧬 Automated Primer Design Tool")
    st.markdown("### Design PCR primers with NCBI database integration and T7 dsRNA functionality")
    
    debug_session_state()
    
    # Sidebar configuration
    st.sidebar.header("⚙️ Configuration")
    
    # NCBI Configuration
    st.sidebar.subheader("NCBI Settings")
    st.sidebar.info("📧 **Email Required**: NCBI requires a valid email address for database access.")
    
    email = st.sidebar.text_input("Email (required for NCBI)", 
                                 placeholder="your.email@example.com",
                                 help="Required by NCBI for database access.")
    
    if st.sidebar.button("🚀 Use test email (demo@example.com)", type="secondary"):
        st.session_state.demo_email = "demo@example.com"
        st.rerun()
    
    if 'demo_email' in st.session_state:
        email = st.session_state.demo_email
        del st.session_state.demo_email
    
    api_key = st.sidebar.text_input("NCBI API Key (optional)", 
                                   type="password",
                                   help="For increased rate limits (not required)")
    
    # Primer Design Parameters
    st.sidebar.subheader("Primer Parameters")
    
    with st.sidebar.expander("Basic Parameters", expanded=True):
        opt_size = st.slider("Optimal primer size", 15, 30, 20)
        min_size = st.slider("Minimum primer size", 15, 25, 18)
        max_size = st.slider("Maximum primer size", 20, 35, 25)
        
        opt_tm = st.slider("Optimal Tm (°C)", 50.0, 70.0, 60.0, 0.5)
        min_tm = st.slider("Minimum Tm (°C)", 45.0, 65.0, 57.0, 0.5)
        max_tm = st.slider("Maximum Tm (°C)", 55.0, 75.0, 63.0, 0.5)
    
    with st.sidebar.expander("Advanced Parameters"):
        min_gc = st.slider("Minimum GC content (%)", 20.0, 50.0, 40.0, 1.0)
        max_gc = st.slider("Maximum GC content (%)", 50.0, 80.0, 60.0, 1.0)
        max_poly_x = st.slider("Max poly-X runs", 3, 6, 4)
        salt_conc = st.slider("Salt concentration (mM)", 10.0, 100.0, 50.0, 1.0)
        num_primers = st.number_input("Number of primer pairs to design", 1, 50, 20, 
                                     help="Maximum number of primer pairs to return")
        
        st.write("Product size ranges:")
        min_product = st.number_input("Minimum product size", 50, 500, 75)
        max_product = st.number_input("Maximum product size", 200, 2000, 1000)
    
    # T7 dsRNA Production Settings
    st.sidebar.subheader("🧬 dsRNA Production")
    enable_t7_dsrna = st.sidebar.checkbox(
        "Add T7 promoters for dsRNA production", 
        value=False,
        help="Adds T7 promoter sequences to both primers for bidirectional transcription and dsRNA synthesis"
    )
    
    if enable_t7_dsrna:
        with st.sidebar.expander("dsRNA Parameters", expanded=True):
            st.info("**T7 dsRNA Production**\nAdds T7 promoter (TAATACGACTCACTATAGGG) to both forward and reverse primers. This enables:\n\n• Bidirectional transcription\n• Double-stranded RNA synthesis\n• RNAi applications\n• Pest management research")
            
            optimal_dsrna_length = st.checkbox(
                "Optimize for dsRNA length (100-500 bp)",
                value=True,
                help="Prioritize primer pairs that produce optimal dsRNA lengths for RNAi"
            )
            
            check_transcription_efficiency = st.checkbox(
                "Check T7 transcription efficiency",
                value=True,
                help="Analyze factors affecting T7 polymerase transcription efficiency"
            )
            
            if optimal_dsrna_length:
                min_product = max(min_product, 100)
                max_product = min(max_product, 500)
                st.write(f"**Adjusted for dsRNA:** {min_product}-{max_product} bp")
    
    custom_params = {
        'PRIMER_OPT_SIZE': opt_size,
        'PRIMER_MIN_SIZE': min_size,
        'PRIMER_MAX_SIZE': max_size,
        'PRIMER_OPT_TM': opt_tm,
        'PRIMER_MIN_TM': min_tm,
        'PRIMER_MAX_TM': max_tm,
        'PRIMER_MIN_GC': min_gc,
        'PRIMER_MAX_GC': max_gc,
        'PRIMER_MAX_POLY_X': max_poly_x,
        'PRIMER_SALT_MONOVALENT': salt_conc,
        'PRIMER_PRODUCT_SIZE_RANGE': [[min_product, max_product]],
        'PRIMER_NUM_RETURN': num_primers
    }
    
    # Main content area
    tab1, tab2, tab3, tab4 = st.tabs([
        "📝 Input", 
        "🔬 Results", 
        "📊 Analysis", 
        "💾 Export"
    ])
    
    state_check = check_session_state_validity()
    
    with tab1:
        st.header("Sequence Input")
        
        if state_check['has_primers']:
            st.success(f"✅ Current session: {state_check['primer_count']} primers designed for {state_check['sequence_length']:,} bp sequence")
        
        input_method = st.radio(
            "Choose input method:",
            ["Organism Name", "GenBank ID", "NCBI Search", "Direct Sequence", "Upload File"]
        )
        
        if input_method == "Organism Name":
            st.subheader("Search by Organism")
            
            if not email:
                st.error("❌ **Email Required**: Please enter an email address in the sidebar first to search for organisms.")
            else:
                st.info("💡 **Tip:** Enter the scientific name (e.g., 'Fusarium oxysporum') for best results.")
            
            col1, col2 = st.columns([2, 1])
            with col1:
                organism_name = st.text_input("Enter organism name:", 
                                            placeholder="e.g., Fusarium oxysporum, Coronavirus, Tetranychus urticae")
            
            with col2:
                max_sequences = st.number_input("Max sequences to search", min_value=5, max_value=50, value=10)
            
            # Analysis mode selection
            st.subheader("Analysis Mode")
            analysis_mode = st.radio(
                "Choose analysis approach:",
                ["Standard (Single Sequence)", "Advanced (Conservation Analysis)"],
                help="Standard: Use first sequence found. Advanced: Analyze conservation across multiple sequences."
            )
            
            # Advanced options (only show for Advanced mode)
            if analysis_mode == "Advanced (Conservation Analysis)":
                with st.expander("Advanced Conservation Analysis Options", expanded=True):
                    st.info("🧬 **Conservation Analysis Mode**: Find conserved regions across multiple sequences for robust primer design")
                    
                    col1, col2 = st.columns(2)
                    with col1:
                        conservation_threshold = st.slider(
                            "Conservation threshold (%)",
                            min_value=60, max_value=95, value=80,
                            help="Minimum conservation required across sequences"
                        ) / 100
                        
                        window_size = st.slider(
                            "Analysis window size (bp)",
                            min_value=100, max_value=500, value=200,
                            help="Size of regions to analyze for conservation"
                        )
                    
                    with col2:
                        enable_specificity_testing = st.checkbox(
                            "Enable specificity testing",
                            value=True,
                            help="Test primer specificity against related organisms"
                        )
                        
                        if enable_specificity_testing:
                            specificity_threshold = st.slider(
                                "Specificity threshold (%)",
                                min_value=50, max_value=90, value=70,
                                help="Maximum allowable similarity to other organisms"
                            ) / 100
                    
                    # Custom comparison organisms
                    comparison_organisms = st.text_input(
                        "Custom comparison organisms (comma-separated, optional):",
                        placeholder="e.g., Aspergillus niger, Penicillium chrysogenum",
                        help="Leave empty for automatic selection based on target organism"
                    )
            
            # Search button
            if st.button("🔍 Search Organism Sequences", type="primary", use_container_width=True):
                if not email:
                    st.error("❌ **Email Required**: Please enter an email address in the sidebar.")
                elif not organism_name:
                    st.error("❌ **Organism Name Required**: Please enter an organism name.")
                else:
                    if analysis_mode == "Standard (Single Sequence)":
                        # EXISTING STANDARD METHOD (keep as-is)
                        with st.spinner(f"Searching for {organism_name} genomes..."):
                            try:
                                ncbi = NCBIConnector(email, api_key)
                                designer = PrimerDesigner()
                                
                                search_query = f'"{organism_name}"[organism]'
                                st.write(f"Searching with query: `{search_query}`")
                                
                                seq_ids = ncbi.search_sequences(search_query, database="nucleotide", max_results=max_sequences)
                                
                                if seq_ids:
                                    st.success(f"Found {len(seq_ids)} sequences!")
                                    
                                    seq_id = seq_ids[0]
                                    st.info(f"Using sequence {seq_id} for primer design...")
                                    
                                    sequence = ncbi.fetch_sequence(seq_id)
                                    seq_info = ncbi.fetch_sequence_info(seq_id)
                                    
                                    if sequence:
                                        clean_sequence = re.sub(r'[^ATGCatgc]', '', sequence.upper())
                                        
                                        if len(clean_sequence) < 50:
                                            st.error("Sequence too short for primer design")
                                        else:
                                            if len(clean_sequence) > 100000:
                                                st.warning(f"Large sequence ({len(clean_sequence):,} bp). Using first 100kb.")
                                                clean_sequence = clean_sequence[:100000]
                                            
                                            st.session_state.current_sequence = clean_sequence
                                            st.session_state.sequence_info = seq_info or {
                                                "id": seq_id,
                                                "description": f"Sequence {seq_id}",
                                                "length": len(clean_sequence),
                                                "organism": organism_name
                                            }
                                            
                                            st.write("Designing primers...")
                                            primers = designer.design_primers(
                                                clean_sequence, 
                                                custom_params=custom_params,
                                                add_t7_promoter=enable_t7_dsrna
                                            )
                                            st.session_state.primers_designed = primers
                                            
                                            if enable_t7_dsrna:
                                                st.session_state.t7_dsrna_enabled = True
                                                st.session_state.t7_settings = {
                                                    'optimal_length': optimal_dsrna_length,
                                                    'check_efficiency': check_transcription_efficiency
                                                }
                                            else:
                                                st.session_state.t7_dsrna_enabled = False
                                            
                                            if primers:
                                                st.success(f"✅ Successfully designed {len(primers)} primer pairs!")
                                                
                                                preview_data = []
                                                for i, primer in enumerate(primers[:5]):
                                                    preview_data.append({
                                                        'Pair': i + 1,
                                                        'Forward': primer.forward_seq[:30] + '...' if len(primer.forward_seq) > 30 else primer.forward_seq,
                                                        'Reverse': primer.reverse_seq[:30] + '...' if len(primer.reverse_seq) > 30 else primer.reverse_seq,
                                                        'Product Size': f"{primer.product_size} bp"
                                                    })
                                                
                                                st.dataframe(pd.DataFrame(preview_data), use_container_width=True)
                                                st.info("📊 Go to other tabs to view detailed analysis!")
                                            else:
                                                st.warning("No suitable primers found. Try adjusting parameters.")
                                    else:
                                        st.error("Failed to fetch sequence")
                                else:
                                    st.warning(f"No sequences found for {organism_name}")
                                    
                            except Exception as e:
                                st.error(f"Error: {e}")
                    
                    else:  # Advanced Conservation Analysis
                        # NEW ADVANCED METHOD
                        with st.spinner(f"Performing advanced conservation analysis for {organism_name}..."):
                            try:
                                # Initialize managers
                                sequence_manager = SequenceManager(NCBIConnector(email, api_key))
                                analyzer = ConservationAnalyzer(NCBIConnector(email, api_key))
                                
                                # Step 1: Fetch sequences
                                st.write("🔍 **Step 1: Fetching sequences...**")
                                sequences = sequence_manager.fetch_organism_sequences(
                                    organism_name, max_sequences
                                )
                                
                                if not sequences:
                                    st.warning(f"No sequences found for {organism_name}")
                                    return
                                
                                st.success(f"Found {len(sequences)} sequences!")
                                
                                # Display sequence summary
                                sequence_data = []
                                for i, seq in enumerate(sequences):
                                    sequence_data.append({
                                        'ID': seq['id'],
                                        'Description': seq['description'][:80] + '...' if len(seq['description']) > 80 else seq['description'],
                                        'Length': f"{seq['length']:,} bp"
                                    })
                                
                                sequence_df = pd.DataFrame(sequence_data)
                                st.dataframe(sequence_df, use_container_width=True)
                                
                                # Step 2: Conservation analysis
                                st.write("🧬 **Step 2: Analyzing conservation...**")
                                seq_list = [seq['sequence'] for seq in sequences]
                                
                                conserved_regions = analyzer.analyze_multiple_sequences(
                                    seq_list,
                                    min_conservation=conservation_threshold,
                                    window_size=window_size
                                )
                                
                                if not conserved_regions:
                                    st.warning("No conserved regions found with current parameters. Try lowering the conservation threshold.")
                                    return
                                
                                st.success(f"Found {len(conserved_regions)} conserved regions!")
                                
                                # Display conserved regions
                                conservation_data = []
                                for i, region in enumerate(conserved_regions):
                                    conservation_data.append({
                                        'Region': i + 1,
                                        'Position': f"{region['start']}-{region['end']}",
                                        'Length': f"{region['length']} bp",
                                        'Conservation': f"{region['conservation_score']:.1%}",
                                        'Sequences': region['sequence_count']
                                    })
                                
                                conservation_df = pd.DataFrame(conservation_data)
                                st.dataframe(conservation_df, use_container_width=True)
                                
                                # Step 3: Specificity testing
                                specificity_results = {}
                                if enable_specificity_testing:
                                    st.write("🎯 **Step 3: Testing specificity...**")
                                    
                                    # Select best conserved region for testing
                                    best_region = max(conserved_regions, key=lambda x: x['conservation_score'])
                                    consensus_seq = best_region['consensus_sequence']
                                    
                                    # Determine comparison organisms
                                    if comparison_organisms.strip():
                                        comp_orgs = [org.strip() for org in comparison_organisms.split(',')]
                                    else:
                                        comp_orgs = get_related_organisms(organism_name)
                                    
                                    st.write(f"Testing against: {', '.join(comp_orgs)}")
                                    
                                    specificity_results = analyzer.test_specificity(
                                        consensus_seq,
                                        comp_orgs,
                                        max_similarity=specificity_threshold
                                    )
                                    
                                    # Display specificity results
                                    specificity_data = []
                                    for org, result in specificity_results.items():
                                        if 'error' not in result:
                                            specificity_data.append({
                                                'Organism': org,
                                                'Max Similarity': f"{result['max_similarity']:.1%}",
                                                'Specific': '✅' if result['is_specific'] else '❌',
                                                'Sequences Tested': result['sequences_tested']
                                            })
                                    
                                    if specificity_data:
                                        specificity_df = pd.DataFrame(specificity_data)
                                        st.dataframe(specificity_df, use_container_width=True)
                                
                                # Step 4: Design primers
                                st.write("🧬 **Step 4: Designing primers from best conserved region...**")
                                
                                # Use the most conserved region
                                best_region = max(conserved_regions, key=lambda x: x['conservation_score'])
                                consensus_seq = best_region['consensus_sequence']
                                
                                # Store sequence and design primers
                                st.session_state.current_sequence = consensus_seq
                                st.session_state.sequence_info = {
                                    "id": f"consensus_region",
                                    "description": f"Conserved region ({best_region['conservation_score']:.1%} conserved) from {len(sequences)} {organism_name} sequences",
                                    "length": len(consensus_seq),
                                    "organism": organism_name,
                                    "conservation_score": best_region['conservation_score'],
                                    "region_info": best_region
                                }
                                
                                # Store conservation analysis results
                                st.session_state.conserved_regions = conserved_regions
                                st.session_state.conservation_sequences = sequences
                                st.session_state.specificity_results = specificity_results
                                
                                # Store analysis metadata
                                st.session_state.analysis_metadata = {
                                    'type': 'conservation_based',
                                    'sequences_analyzed': len(sequences),
                                    'conservation_threshold': conservation_threshold,
                                    'specificity_tested': enable_specificity_testing,
                                    'specificity_threshold': specificity_threshold if enable_specificity_testing else None
                                }
                                
                                # Design primers
                                designer = PrimerDesigner()
                                primers = designer.design_primers(
                                    consensus_seq,
                                    custom_params=custom_params,
                                    add_t7_promoter=enable_t7_dsrna
                                )
                                
                                st.session_state.primers_designed = primers
                                
                                if enable_t7_dsrna:
                                    st.session_state.t7_dsrna_enabled = True
                                
                                if primers:
                                    st.success(f"✅ Designed {len(primers)} conservation-based primers!")
                                    st.info("📊 Go to Results tab to view detailed analysis!")
                                    
                                    # Show preview
                                    preview_data = []
                                    for i, primer in enumerate(primers[:3]):
                                        preview_data.append({
                                            'Pair': i + 1,
                                            'Forward': primer.forward_seq[:30] + '...' if len(primer.forward_seq) > 30 else primer.forward_seq,
                                            'Reverse': primer.reverse_seq[:30] + '...' if len(primer.reverse_seq) > 30 else primer.reverse_seq,
                                            'Product Size': f"{primer.product_size} bp"
                                        })
                                    
                                    st.dataframe(pd.DataFrame(preview_data), use_container_width=True)
                                else:
                                    st.warning("No suitable primers found for conserved region")
                                
                            except Exception as e:
                                st.error(f"Advanced analysis error: {e}")
                                import traceback
                                st.code(traceback.format_exc())
            
            # Agricultural Pests & Pathogens section with improved layout
            st.markdown("---")
            st.markdown("### 🎯 Quick Select: Agricultural Pests & Pathogens")
            st.markdown("*Click any button below to automatically search for that organism*")
            
            suggestions = get_organism_suggestions()
            
            # Show summary of available organisms
            total_organisms = sum(len(orgs) for subcats in suggestions.values() for orgs in subcats.values())
            st.info(f"📊 **{total_organisms} organisms** available across {len(suggestions)} categories")
            
            # Create expandable sections for better organization
            for category, subcategories in suggestions.items():
                with st.expander(f"{category} ({sum(len(orgs) for orgs in subcategories.values())} organisms)", expanded=False):
                    for subcategory, organisms in subcategories.items():
                        st.markdown(f"**{subcategory}**")
                        
                        # Use a more compact grid layout
                        num_cols = min(len(organisms), 4)  # Max 4 columns for better readability
                        cols = st.columns(num_cols)
                        
                        for i, (common_name, latin_name) in enumerate(organisms):
                            with cols[i % num_cols]:
                                if st.button(
                                    common_name, 
                                    key=f"suggest_{category}_{subcategory}_{i}", 
                                    help=f"Search for {latin_name}",
                                    use_container_width=True
                                ):
                                    # Set the organism name and trigger search
                                    st.session_state.selected_organism = latin_name
                                    st.rerun()
                        
                        # Add small spacing between subcategories
                        if subcategory != list(subcategories.keys())[-1]:  # Not the last subcategory
                            st.markdown("")
            
            # Handle organism selection from buttons
            if 'selected_organism' in st.session_state:
                organism_name = st.session_state.selected_organism
                del st.session_state.selected_organism
                
                if not email:
                    st.error("❌ **Email Required**: Please enter an email address in the sidebar.")
                else:
                    with st.spinner(f"Searching for {organism_name} genomes..."):
                        try:
                            ncbi = NCBIConnector(email, api_key)
                            designer = PrimerDesigner()
                            
                            search_query = f'"{organism_name}"[organism]'
                            st.write(f"Searching with query: `{search_query}`")
                            
                            seq_ids = ncbi.search_sequences(search_query, database="nucleotide", max_results=max_sequences)
                            
                            if seq_ids:
                                st.success(f"Found {len(seq_ids)} sequences!")
                                
                                seq_id = seq_ids[0]
                                st.info(f"Using sequence {seq_id} for primer design...")
                                
                                sequence = ncbi.fetch_sequence(seq_id)
                                seq_info = ncbi.fetch_sequence_info(seq_id)
                                
                                if sequence:
                                    clean_sequence = re.sub(r'[^ATGCatgc]', '', sequence.upper())
                                    
                                    if len(clean_sequence) < 50:
                                        st.error("Sequence too short for primer design")
                                    else:
                                        if len(clean_sequence) > 100000:
                                            st.warning(f"Large sequence ({len(clean_sequence):,} bp). Using first 100kb.")
                                            clean_sequence = clean_sequence[:100000]
                                        
                                        st.session_state.current_sequence = clean_sequence
                                        st.session_state.sequence_info = seq_info or {
                                            "id": seq_id,
                                            "description": f"Sequence {seq_id}",
                                            "length": len(clean_sequence),
                                            "organism": organism_name
                                        }
                                        
                                        st.write("Designing primers...")
                                        primers = designer.design_primers(
                                            clean_sequence, 
                                            custom_params=custom_params,
                                            add_t7_promoter=enable_t7_dsrna
                                        )
                                        st.session_state.primers_designed = primers
                                        
                                        if enable_t7_dsrna:
                                            st.session_state.t7_dsrna_enabled = True
                                            st.session_state.t7_settings = {
                                                'optimal_length': optimal_dsrna_length,
                                                'check_efficiency': check_transcription_efficiency
                                            }
                                        else:
                                            st.session_state.t7_dsrna_enabled = False
                                        
                                        if primers:
                                            st.success(f"✅ Successfully designed {len(primers)} primer pairs!")
                                            
                                            preview_data = []
                                            for i, primer in enumerate(primers[:5]):
                                                preview_data.append({
                                                    'Pair': i + 1,
                                                    'Forward': primer.forward_seq[:30] + '...' if len(primer.forward_seq) > 30 else primer.forward_seq,
                                                    'Reverse': primer.reverse_seq[:30] + '...' if len(primer.reverse_seq) > 30 else primer.reverse_seq,
                                                    'Product Size': f"{primer.product_size} bp"
                                                })
                                            
                                            st.dataframe(pd.DataFrame(preview_data), use_container_width=True)
                                            st.info("📊 Go to other tabs to view detailed analysis!")
                                        else:
                                            st.warning("No suitable primers found. Try adjusting parameters.")
                                else:
                                    st.error("Failed to fetch sequence")
                            else:
                                st.warning(f"No sequences found for {organism_name}")
                                
                        except Exception as e:
                            st.error(f"Error: {e}")
        
        elif input_method == "Direct Sequence":
            sequence_input = st.text_area("Enter DNA sequence:", 
                                        placeholder="ATGCGATCGATCG...",
                                        height=150)
            
            if st.button("Design Primers", type="primary"):
                if not sequence_input:
                    st.error("Please provide a DNA sequence")
                else:
                    with st.spinner("Designing primers..."):
                        try:
                            designer = PrimerDesigner()
                            
                            clean_seq = re.sub(r'[^ATGCatgc]', '', sequence_input.upper())
                            
                            st.session_state.current_sequence = clean_seq
                            st.session_state.sequence_info = {
                                "length": len(clean_seq),
                                "description": "User-provided sequence",
                                "organism": "User input",
                                "id": "user_sequence"
                            }
                            
                            primers = designer.design_primers(
                                clean_seq, 
                                custom_params=custom_params,
                                add_t7_promoter=enable_t7_dsrna
                            )
                            st.session_state.primers_designed = primers
                            
                            if enable_t7_dsrna:
                                st.session_state.t7_dsrna_enabled = True
                                st.session_state.t7_settings = {
                                    'optimal_length': optimal_dsrna_length,
                                    'check_efficiency': check_transcription_efficiency
                                }
                            else:
                                st.session_state.t7_dsrna_enabled = False
                            
                            if primers:
                                st.success(f"Successfully designed {len(primers)} primer pairs!")
                                st.info("📊 Go to other tabs to view detailed analysis!")
                            else:
                                st.warning("No suitable primers found with current parameters")
                        except Exception as e:
                            st.error(f"Error: {e}")
    
    with tab2:
        st.header("Primer Design Results")
        
        if not state_check['has_primers']:
            st.info("No primers designed yet. Please use the Input tab to design primers.")
            return
        
        primers = st.session_state.primers_designed
        t7_enabled = st.session_state.get('t7_dsrna_enabled', False)
        
        if t7_enabled:
            st.info("🧬 **T7 dsRNA Mode Active** - Primers include T7 promoter sequences for double-stranded RNA production")
        
        # Check if this is conservation-based analysis
        analysis_metadata = st.session_state.get('analysis_metadata', {})
        is_conservation_based = analysis_metadata.get('type') == 'conservation_based'
        
        if is_conservation_based:
            st.success("🧬 **Conservation-Based Design** - Primers designed from conserved regions across multiple sequences")
            
            col1, col2, col3, col4 = st.columns(4)
            with col1:
                st.metric("Sequences Analyzed", analysis_metadata.get('sequences_analyzed', 'N/A'))
            with col2:
                conservation_thresh = analysis_metadata.get('conservation_threshold', 0)
                st.metric("Conservation Threshold", f"{conservation_thresh:.0%}")
            with col3:
                specificity_tested = analysis_metadata.get('specificity_tested', False)
                st.metric("Specificity Tested", "Yes" if specificity_tested else "No")
            with col4:
                if specificity_tested:
                    spec_thresh = analysis_metadata.get('specificity_threshold', 0)
                    st.metric("Specificity Threshold", f"{spec_thresh:.0%}")
        
        # Sequence information
        if st.session_state.sequence_info:
            st.subheader("Sequence Information")
            info = st.session_state.sequence_info
            
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Length", f"{info.get('length', 'N/A'):,} bp")
            with col2:
                st.metric("Organism", info.get('organism', 'N/A'))
            with col3:
                st.metric("ID", info.get('id', 'N/A'))
            
            if 'description' in info:
                st.write(f"**Description:** {info['description']}")
            
            # Show conservation information if available
            if 'conservation_score' in info:
                st.metric("Conservation Score", f"{info['conservation_score']:.1%}")
        
        # Show conserved regions if available
        if hasattr(st.session_state, 'conserved_regions') and st.session_state.conserved_regions:
            st.subheader("Conserved Regions Analysis")
            
            conserved_regions = st.session_state.conserved_regions
            conservation_data = []
            
            for i, region in enumerate(conserved_regions):
                conservation_data.append({
                    'Region': i + 1,
                    'Position': f"{region['start']}-{region['end']}",
                    'Length': f"{region['length']} bp",
                    'Conservation': f"{region['conservation_score']:.1%}",
                    'Sequences': region['sequence_count']
                })
            
            conservation_df = pd.DataFrame(conservation_data)
            st.dataframe(conservation_df, use_container_width=True)
        
        # Show specificity results if available
        if hasattr(st.session_state, 'specificity_results') and st.session_state.specificity_results:
            st.subheader("Specificity Testing Results")
            
            specificity_results = st.session_state.specificity_results
            specificity_data = []
            
            for organism, result in specificity_results.items():
                if 'error' not in result:
                    specificity_data.append({
                        'Organism': organism,
                        'Max Similarity': f"{result['max_similarity']:.1%}",
                        'Specific': '✅ Yes' if result['is_specific'] else '❌ No',
                        'Sequences Tested': result['sequences_tested']
                    })
            
            if specificity_data:
                specificity_df = pd.DataFrame(specificity_data)
                st.dataframe(specificity_df, use_container_width=True)
                
                # Summary
                total_orgs = len(specificity_data)
                specific_orgs = sum(1 for row in specificity_data if row['Specific'] == '✅ Yes')
                specificity_percentage = (specific_orgs / total_orgs) * 100 if total_orgs > 0 else 0
                
                if specificity_percentage >= 80:
                    st.success(f"🎯 Excellent specificity: {specific_orgs}/{total_orgs} organisms ({specificity_percentage:.0f}%)")
                elif specificity_percentage >= 60:
                    st.info(f"🎯 Good specificity: {specific_orgs}/{total_orgs} organisms ({specificity_percentage:.0f}%)")
                else:
                    st.warning(f"⚠️ Moderate specificity: {specific_orgs}/{total_orgs} organisms ({specificity_percentage:.0f}%)")
        
        # Primer results table
        st.subheader("Primer Pairs")
        
        data = []
        for i, primer in enumerate(primers):
            if hasattr(primer, 'has_t7_promoter') and primer.has_t7_promoter:
                row = {
                    'Pair': i + 1,
                    'Forward (with T7)': primer.forward_seq,
                    'Reverse (with T7)': primer.reverse_seq,
                    'Core Forward': primer.core_forward_seq,
                    'Core Reverse': primer.core_reverse_seq,
                    'Core Tm': f"{primer.forward_tm:.1f}°C / {primer.reverse_tm:.1f}°C",
                    'dsRNA Size': f"{primer.product_size} bp",
                    'Core GC%': f"{primer.gc_content_f:.1f}% / {primer.gc_content_r:.1f}%",
                    'Penalty': f"{primer.penalty:.3f}"
                }
            else:
                row = {
                    'Pair': i + 1,
                    'Forward Sequence': primer.forward_seq,
                    'Reverse Sequence': primer.reverse_seq,
                    'Forward Tm': f"{primer.forward_tm:.1f}°C",
                    'Reverse Tm': f"{primer.reverse_tm:.1f}°C",
                    'Product Size': f"{primer.product_size} bp",
                    'Forward GC%': f"{primer.gc_content_f:.1f}%",
                    'Reverse GC%': f"{primer.gc_content_r:.1f}%",
                    'Penalty': f"{primer.penalty:.3f}"
                }
            
            data.append(row)
        
        df = pd.DataFrame(data)
        st.dataframe(df, use_container_width=True)
        
        # Detailed view for selected primer
        st.subheader("Detailed View")
        selected_primer = st.selectbox("Select primer pair for details:", 
                                     range(len(primers)), 
                                     format_func=lambda x: f"Pair {x+1}")
        
        if selected_primer < len(primers):
            primer = primers[selected_primer]
            
            if hasattr(primer, 'has_t7_promoter') and primer.has_t7_promoter:
                st.info("🧬 **T7 dsRNA Primer Pair** - Includes T7 promoter for double-stranded RNA synthesis")
                
                col1, col2 = st.columns(2)
                
                with col1:
                    st.write("**Forward Primer (with T7)**")
                    st.code(primer.forward_seq, language="text")
                    st.write("**Core Forward Primer**")
                    st.code(primer.core_forward_seq, language="text")
                    st.write(f"- Core Position: {primer.forward_start}")
                    st.write(f"- Core Length: {len(primer.core_forward_seq)} bp")
                    st.write(f"- Core Tm: {primer.forward_tm:.2f}°C")
                    st.write(f"- Core GC Content: {primer.gc_content_f:.1f}%")
                    st.write(f"- Full Length (with T7): {len(primer.forward_seq)} bp")
                
                with col2:
                    st.write("**Reverse Primer (with T7)**")
                    st.code(primer.reverse_seq, language="text")
                    st.write("**Core Reverse Primer**")
                    st.code(primer.core_reverse_seq, language="text")
                    st.write(f"- Core Position: {primer.reverse_start}")
                    st.write(f"- Core Length: {len(primer.core_reverse_seq)} bp")
                    st.write(f"- Core Tm: {primer.reverse_tm:.2f}°C")
                    st.write(f"- Core GC Content: {primer.gc_content_r:.1f}%")
                    st.write(f"- Full Length (with T7): {len(primer.reverse_seq)} bp")
                
                st.write(f"**dsRNA Product Size:** {primer.product_size} bp")
                st.write(f"**T7 Promoter:** {primer.t7_promoter_seq}")
                
                # dsRNA analysis
                if st.session_state.current_sequence:
                    designer = PrimerDesigner()
                    dsrna_props = designer.calculate_dsrna_properties(primer, st.session_state.current_sequence)
                    
                    if dsrna_props:
                        st.subheader("dsRNA Analysis")
                        
                        col1, col2, col3 = st.columns(3)
                        with col1:
                            st.metric("dsRNA Length", f"{dsrna_props.get('dsrna_length', 'N/A')} bp")
                        with col2:
                            st.metric("dsRNA GC Content", f"{dsrna_props.get('dsrna_gc_content', 'N/A'):.1f}%")
                        with col3:
                            st.metric("T7 Efficiency", dsrna_props.get('transcription_efficiency', 'N/A'))
                        
                        st.write("**dsRNA Quality Indicators:**")
                        if dsrna_props.get('optimal_length'):
                            st.success("✅ Optimal length for RNAi (100-500 bp)")
                        else:
                            st.warning("⚠️ Length outside optimal range for RNAi")
                        
                        if dsrna_props.get('moderate_gc'):
                            st.success("✅ Moderate GC content (40-60%)")
                        else:
                            st.warning("⚠️ Extreme GC content may affect efficiency")
                        
                        st.write(f"**Transcription Start:** {dsrna_props.get('transcription_start', 'N/A')} (G is optimal for T7)")
                        st.write(f"**Estimated Yield:** {dsrna_props.get('estimated_yield', 'N/A')}")
                        
                        if 'target_sequence' in dsrna_props:
                            st.write("**Target Sequence (first 100 bp):**")
                            st.code(dsrna_props['target_sequence'], language="text")
            else:
                col1, col2 = st.columns(2)
                
                with col1:
                    st.write("**Forward Primer**")
                    st.code(primer.forward_seq, language="text")
                    st.write(f"- Position: {primer.forward_start}")
                    st.write(f"- Length: {len(primer.forward_seq)} bp")
                    st.write(f"- Tm: {primer.forward_tm:.2f}°C")
                    st.write(f"- GC Content: {primer.gc_content_f:.1f}%")
                
                with col2:
                    st.write("**Reverse Primer**")
                    st.code(primer.reverse_seq, language="text")
                    st.write(f"- Position: {primer.reverse_start}")
                    st.write(f"- Length: {len(primer.reverse_seq)} bp")
                    st.write(f"- Tm: {primer.reverse_tm:.2f}°C")
                    st.write(f"- GC Content: {primer.gc_content_r:.1f}%")
                
                st.write(f"**Product Size:** {primer.product_size} bp")
            
            st.write(f"**Penalty Score:** {primer.penalty:.4f}")
        
        # Conservation and Specificity Summary
        if is_conservation_based:
            with st.expander("Analysis Summary", expanded=False):
                st.write("**Conservation-Based Primer Design Summary:**")
                st.write(f"- Analyzed {analysis_metadata.get('sequences_analyzed', 'N/A')} sequences from {st.session_state.get('target_organism', 'target organism')}")
                st.write(f"- Used {analysis_metadata.get('conservation_threshold', 0):.0%} conservation threshold")
                
                if analysis_metadata.get('specificity_tested'):
                    st.write(f"- Tested specificity against related organisms")
                    st.write(f"- Used {analysis_metadata.get('specificity_threshold', 0):.0%} specificity threshold")
                    
                    if hasattr(st.session_state, 'specificity_results'):
                        total_tested = len(st.session_state.specificity_results)
                        specific_count = sum(1 for result in st.session_state.specificity_results.values() 
                                           if result.get('is_specific', False))
                        st.write(f"- Specificity results: {specific_count}/{total_tested} organisms passed")
                
                st.write(f"- Final result: {len(primers)} primer pairs designed")
        
        # dsRNA Production Protocol
        if t7_enabled and primers:
            st.subheader("dsRNA Production Protocol")
            
            with st.expander("Step-by-Step dsRNA Synthesis Protocol", expanded=False):
                st.markdown("""
                **Materials Required:**
                - T7 RNA Polymerase
                - NTP mix (ATP, CTP, GTP, UTP)
                - T7 transcription buffer
                - RNase-free water
                - DNase I
                - Phenol-chloroform (optional)
                - Ethanol precipitation reagents
                
                **Protocol:**
                
                **Step 1: PCR Amplification**
                1. Use the T7-tagged primers to amplify your target region
                2. Confirm PCR product size by gel electrophoresis
                3. Purify PCR product using standard methods
                
                **Step 2: In Vitro Transcription**
                1. Set up T7 transcription reaction:
                   - 1-2 μg PCR template
                   - 2 mM each NTP
                   - 1× T7 transcription buffer
                   - 20-40 units T7 RNA Polymerase
                   - RNase-free water to 20 μL
                2. Incubate at 37°C for 2-4 hours
                
                **Step 3: DNase Treatment**
                1. Add 2 units DNase I
                2. Incubate at 37°C for 15 minutes
                
                **Step 4: dsRNA Formation**
                1. Heat to 95°C for 5 minutes
                2. Cool slowly to room temperature (30-60 minutes)
                3. This allows sense and antisense strands to anneal
                
                **Step 5: Purification**
                1. Phenol-chloroform extraction (optional)
                2. Ethanol precipitation
                3. Resuspend in RNase-free water
                
                **Step 6: Quality Control**
                1. Check dsRNA by gel electrophoresis
                2. Quantify using spectrophotometer
                3. Store at -80°C
                """)
            
            # Tips for conservation-based dsRNA
            if is_conservation_based:
                with st.expander("Conservation-Based dsRNA Tips", expanded=False):
                    st.markdown("""
                    **Advantages of Conservation-Based dsRNA:**
                    - Higher likelihood of success across different strains/populations
                    - Reduced chance of resistance development
                    - More robust RNAi response
                    
                    **Additional Considerations:**
                    - Test dsRNA against multiple target populations when possible
                    - Monitor for potential cross-reactivity with beneficial organisms
                    - Consider seasonal/geographical variations in target populations
                    
                    **Quality Control for Conservation-Based Primers:**
                    - Verify PCR amplification across multiple target samples
                    - Test dsRNA efficacy on representative target populations
                    - Monitor for potential off-target effects
                    """)
    
    with tab3:
        st.header("Primer Analysis")
        
        if not state_check['has_primers']:
            st.info("No primers designed yet. Please use the Input tab to design primers.")
            return
        
        primers = st.session_state.primers_designed
        t7_enabled = st.session_state.get('t7_dsrna_enabled', False)
        
        if t7_enabled:
            st.info("🧬 **T7 dsRNA Analysis Mode** - Showing analysis for dsRNA production primers")
        
        fig = create_primer_visualization(primers)
        if fig:
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.warning("Could not create primer visualization")
        
        if state_check['has_sequence']:
            st.subheader("Primer Binding Sites")
            selected_for_diagram = st.selectbox(
                "Select primer pair for binding site visualization:", 
                range(len(primers)), 
                format_func=lambda x: f"Pair {x+1}",
                key="diagram_select"
            )
            
            seq_fig = create_sequence_diagram(
                st.session_state.current_sequence, 
                primers, 
                selected_for_diagram
            )
            if seq_fig:
                st.plotly_chart(seq_fig, use_container_width=True)
        
        # dsRNA-specific analysis
        if t7_enabled and primers:
            st.subheader("dsRNA Production Analysis")
            
            designer = PrimerDesigner()
            dsrna_analysis = []
            
            for i, primer in enumerate(primers):
                if hasattr(primer, 'has_t7_promoter') and primer.has_t7_promoter:
                    props = designer.calculate_dsrna_properties(primer, st.session_state.current_sequence)
                    if props:
                        dsrna_analysis.append({
                            'Pair': i + 1,
                            'dsRNA Length': f"{props.get('dsrna_length', 'N/A')} bp",
                            'GC Content': f"{props.get('dsrna_gc_content', 'N/A'):.1f}%",
                            'T7 Efficiency': props.get('transcription_efficiency', 'N/A'),
                            'Optimal Length': '✅' if props.get('optimal_length') else '❌',
                            'Moderate GC': '✅' if props.get('moderate_gc') else '❌',
                            'Estimated Yield': props.get('estimated_yield', 'N/A')
                        })
            
            if dsrna_analysis:
                dsrna_df = pd.DataFrame(dsrna_analysis)
                st.dataframe(dsrna_df, use_container_width=True)
                
                st.subheader("dsRNA Quality Summary")
                total_pairs = len(dsrna_analysis)
                optimal_length_count = sum(1 for row in dsrna_analysis if row['Optimal Length'] == '✅')
                moderate_gc_count = sum(1 for row in dsrna_analysis if row['Moderate GC'] == '✅')
                
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.metric("Optimal Length", f"{optimal_length_count}/{total_pairs}")
                with col2:
                    st.metric("Moderate GC", f"{moderate_gc_count}/{total_pairs}")
                with col3:
                    quality_score = (optimal_length_count + moderate_gc_count) / (2 * total_pairs) * 100
                    st.metric("Quality Score", f"{quality_score:.0f}%")
        
        st.subheader("Statistics")
        if primers:
            col1, col2, col3, col4 = st.columns(4)
            
            with col1:
                avg_tm = sum(p.forward_tm + p.reverse_tm for p in primers) / (2 * len(primers))
                if t7_enabled:
                    st.metric("Average Core Tm", f"{avg_tm:.1f}°C")
                else:
                    st.metric("Average Tm", f"{avg_tm:.1f}°C")
            
            with col2:
                avg_gc = sum(p.gc_content_f + p.gc_content_r for p in primers) / (2 * len(primers))
                if t7_enabled:
                    st.metric("Average Core GC%", f"{avg_gc:.1f}%")
                else:
                    st.metric("Average GC%", f"{avg_gc:.1f}%")
            
            with col3:
                avg_product = sum(p.product_size for p in primers) / len(primers)
                if t7_enabled:
                    st.metric("Average dsRNA Size", f"{avg_product:.0f} bp")
                else:
                    st.metric("Average Product Size", f"{avg_product:.0f} bp")
            
            with col4:
                best_penalty = min(p.penalty for p in primers)
                st.metric("Best Penalty Score", f"{best_penalty:.3f}")
    
    with tab4:
        st.header("Export Results")
        
        if not state_check['has_primers']:
            st.info("No primers to export. Please design primers first.")
            return
        
        primers = st.session_state.primers_designed
        t7_enabled = st.session_state.get('t7_dsrna_enabled', False)
        
        if t7_enabled:
            st.info("🧬 **T7 dsRNA Export Mode** - Export includes both full T7 primers and core sequences")
        
        st.subheader("Download Options")
        
        col1, col2 = st.columns(2)
        
        with col1:
            if st.button("📊 Download as Excel", type="primary"):
                excel_data = export_to_excel(primers)
                if excel_data:
                    filename = "t7_dsrna_primers.xlsx" if t7_enabled else "primer_results.xlsx"
                    st.download_button(
                        label="Click to Download Excel File",
                        data=excel_data,
                        file_name=filename,
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                    )
        
        with col2:
            try:
                data = []
                for i, primer in enumerate(primers):
                    if hasattr(primer, 'has_t7_promoter') and primer.has_t7_promoter:
                        row = {
                            'Primer_Pair': i + 1,
                            'Forward_T7_Sequence': primer.forward_seq,
                            'Reverse_T7_Sequence': primer.reverse_seq,
                            'Forward_Core_Sequence': primer.core_forward_seq,
                            'Reverse_Core_Sequence': primer.core_reverse_seq,
                            'Core_Forward_Tm': round(primer.forward_tm, 2),
                            'Core_Reverse_Tm': round(primer.reverse_tm, 2),
                            'dsRNA_Size': primer.product_size,
                            'Core_Forward_GC_Percent': round(primer.gc_content_f, 2),
                            'Core_Reverse_GC_Percent': round(primer.gc_content_r, 2),
                            'Forward_Start': primer.forward_start,
                            'Reverse_Start': primer.reverse_start,
                            'Penalty_Score': round(primer.penalty, 4),
                            'T7_Promoter': primer.t7_promoter_seq,
                            'Primer_Type': 'T7_dsRNA'
                        }
                    else:
                        row = {
                            'Primer_Pair': i + 1,
                            'Forward_Sequence': primer.forward_seq,
                            'Reverse_Sequence': primer.reverse_seq,
                            'Forward_Tm': round(primer.forward_tm, 2),
                            'Reverse_Tm': round(primer.reverse_tm, 2),
                            'Product_Size': primer.product_size,
                            'Forward_GC_Percent': round(primer.gc_content_f, 2),
                            'Reverse_GC_Percent': round(primer.gc_content_r, 2),
                            'Forward_Start': primer.forward_start,
                            'Reverse_Start': primer.reverse_start,
                            'Penalty_Score': round(primer.penalty, 4),
                            'Primer_Type': 'Standard'
                        }
                    
                    data.append(row)
                
                df = pd.DataFrame(data)
                csv = df.to_csv(index=False)
                
                filename = "t7_dsrna_primers.csv" if t7_enabled else "primer_results.csv"
                st.download_button(
                    label="📄 Download as CSV",
                    data=csv,
                    file_name=filename,
                    mime="text/csv"
                )
            except Exception as e:
                st.error(f"Error creating CSV: {e}")
        
        st.subheader("Export Preview")
        try:
            if 'df' in locals():
                st.dataframe(df, use_container_width=True)
        except:
            st.warning("Could not create preview")
        
        st.subheader("Primer Ordering Format")
        
        if t7_enabled:
            st.write("**T7 dsRNA Primer Ordering Format** - Optimized for synthesis companies:")
            
            try:
                ordering_data = []
                for i, primer in enumerate(primers):
                    if hasattr(primer, 'has_t7_promoter') and primer.has_t7_promoter:
                        ordering_data.extend([
                            {
                                'Name': f"T7_Forward_Primer_{i+1}",
                                'Sequence': primer.forward_seq,
                                'Length': len(primer.forward_seq),
                                'Core_Tm': round(primer.forward_tm, 1),
                                'Notes': f"T7 promoter + {len(primer.core_forward_seq)}bp core"
                            },
                            {
                                'Name': f"T7_Reverse_Primer_{i+1}",
                                'Sequence': primer.reverse_seq,
                                'Length': len(primer.reverse_seq),
                                'Core_Tm': round(primer.reverse_tm, 1),
                                'Notes': f"T7 promoter + {len(primer.core_reverse_seq)}bp core"
                            }
                        ])
                
                ordering_df = pd.DataFrame(ordering_data)
                st.dataframe(ordering_df, use_container_width=True)
                
                ordering_csv = ordering_df.to_csv(index=False)
                st.download_button(
                    label="📋 Download T7 Ordering Format",
                    data=ordering_csv,
                    file_name="t7_dsrna_primer_ordering.csv",
                    mime="text/csv"
                )
                
            except Exception as e:
                st.error(f"Error creating T7 ordering format: {e}")
        
        else:
            st.write("Format suitable for ordering from synthesis companies:")
            
            try:
                ordering_data = []
                for i, primer in enumerate(primers):
                    ordering_data.extend([
                        {
                            'Name': f"Forward_Primer_{i+1}",
                            'Sequence': primer.forward_seq,
                            'Length': len(primer.forward_seq),
                            'Tm': round(primer.forward_tm, 1)
                        },
                        {
                            'Name': f"Reverse_Primer_{i+1}",
                            'Sequence': primer.reverse_seq,
                            'Length': len(primer.reverse_seq),
                            'Tm': round(primer.reverse_tm, 1)
                        }
                    ])
                
                ordering_df = pd.DataFrame(ordering_data)
                st.dataframe(ordering_df, use_container_width=True)
                
                ordering_csv = ordering_df.to_csv(index=False)
                st.download_button(
                    label="📋 Download Ordering Format",
                    data=ordering_csv,
                    file_name="primer_ordering.csv",
                    mime="text/csv"
                )
            except Exception as e:
                st.error(f"Error creating ordering format: {e}")
        
        # Protocol export for T7 dsRNA
        if t7_enabled:
            st.subheader("Protocol Export")
            
            protocol_text = f"""
T7 dsRNA Production Protocol
============================

Generated for {len(primers)} primer pairs designed for dsRNA synthesis.

MATERIALS REQUIRED:
- T7 RNA Polymerase (40 U/μL)
- 10× T7 Transcription Buffer
- NTP Mix (25 mM each ATP, CTP, GTP, UTP)
- DNase I (2 U/μL)
- RNase-free water
- 0.5 M EDTA
- Phenol:chloroform:isoamyl alcohol (25:24:1)
- 3 M sodium acetate (pH 5.2)
- 100% ethanol
- 70% ethanol

PROTOCOL:

1. PCR AMPLIFICATION:
   - Use T7-tagged primers to amplify target regions
   - Verify products by gel electrophoresis
   - Purify PCR products using standard methods

2. T7 TRANSCRIPTION SETUP (20 μL reaction):
   - 1-2 μg purified PCR template
   - 2 μL 10× T7 Transcription Buffer
   - 2 μL NTP Mix (2 mM final each)
   - 1 μL T7 RNA Polymerase (40 units)
   - RNase-free water to 20 μL

3. INCUBATION:
   - 37°C for 2-4 hours
   - Optional: Add additional 20 units T7 polymerase after 2h

4. DNASE TREATMENT:
   - Add 1 μL DNase I (2 units)
   - Incubate 37°C for 15 minutes
   - Add 1 μL 0.5 M EDTA to stop reaction

5. dsRNA ANNEALING:
   - Heat to 95°C for 5 minutes
   - Cool slowly to room temperature (30-60 minutes)
   - This allows complementary strands to anneal

6. PURIFICATION:
   - Add equal volume phenol:chloroform
   - Vortex and centrifuge 15,000g for 5 minutes
   - Transfer aqueous phase to new tube
   - Add 1/10 volume 3 M sodium acetate
   - Add 2.5 volumes cold 100% ethanol
   - Precipitate at -20°C for 30 minutes
   - Centrifuge 15,000g for 15 minutes at 4°C
   - Wash pellet with 70% ethanol
   - Air dry and resuspend in RNase-free water

7. QUALITY CONTROL:
   - Analyze by agarose gel electrophoresis
   - Quantify using spectrophotometer (A260/A280 ratio ~2.0)
   - Store at -80°C in small aliquots

EXPECTED YIELDS:
- 10-50 μg dsRNA per 20 μL reaction
- Higher yields with optimal templates (G at +1 position)

TROUBLESHOOTING:
- Low yield: Check template quality, extend incubation
- Degradation: Ensure RNase-free conditions
- Poor annealing: Optimize cooling rate

APPLICATION NOTES:
- For RNAi in insects: 100-500 ng dsRNA per organism
- For plant applications: 1-10 μg/mL in infiltration buffer
- Store dsRNA at -80°C in single-use aliquots
- Always use RNase-free conditions throughout protocol
"""
            
            st.download_button(
                label="📄 Download Complete Protocol",
                data=protocol_text,
                file_name="t7_dsrna_protocol.txt",
                mime="text/plain"
            )
    
    # Footer
    st.markdown("---")
    st.markdown(
        """
        <div style='text-align: center'>
            <p>🧬 Automated Primer Design Tool | Built with Streamlit</p>
            <p><small>Powered by Primer3, Biopython, and NCBI databases</small></p>
        </div>
        """, 
        unsafe_allow_html=True
    )

if __name__ == "__main__":
    main()
