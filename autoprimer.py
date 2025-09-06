#!/usr/bin/env python3
"""
Streamlit Web Application for Automated Primer Design
===================================================

A user-friendly web interface for automated primer design with NCBI integration.

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

# Import the primer design classes (assuming they're in the same file or imported)
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
            # Don't show error for every failed sequence, just log it
            return {}
    
    def fetch_genome_assembly_info(self, assembly_id: str) -> Dict:
        """Fetch genome assembly information from NCBI"""
        try:
            time.sleep(self.rate_limit_delay)
            # Use esummary to get assembly information
            handle = Entrez.esummary(db="assembly", id=assembly_id)
            summary = Entrez.read(handle)
            handle.close()
            
            if summary and len(summary) > 0:
                assembly = summary[0]
                return {
                    "id": assembly_id,
                    "description": assembly.get("AssemblyName", "Unknown Assembly"),
                    "length": "N/A",  # Assembly length is complex to get
                    "organism": assembly.get("SpeciesName", "Unknown"),
                    "assembly_accession": assembly.get("AssemblyAccession", assembly_id)
                }
            return {}
        except Exception as e:
            # Don't show error for every failed assembly, just log it
            return {}

class PrimerDesigner:
    """Main primer design class with advanced parameters"""
    
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
            'PRIMER_PRODUCT_SIZE_RANGE': [[75, 300], [300, 600], [600, 1000]]
        }
    
    def calculate_gc_content(self, sequence: str) -> float:
        gc_count = sequence.upper().count('G') + sequence.upper().count('C')
        return (gc_count / len(sequence)) * 100 if sequence else 0
    
    def design_primers(self, sequence: str, target_region: Optional[Tuple[int, int]] = None,
                      custom_params: Optional[Dict] = None) -> List[PrimerPair]:
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
                primers.append(primer_pair)
            
            return primers
            
        except Exception as e:
            st.error(f"Error in primer design: {e}")
            return []

class ComprehensiveTaxonomicAnalyzer:
    """Advanced taxonomic analysis for genus-level primer design"""
    
    def __init__(self, ncbi_connector):
        self.ncbi = ncbi_connector
        self.min_sequences_per_subspecies = 50
        self.max_sequences_per_subspecies = 100
        self.conservation_threshold = 0.90  # 90% conservation across subspecies
        self.specificity_threshold = 0.75   # <75% similarity to other taxa
        
    def discover_subspecies(self, genus_name, max_subspecies=20):
        """Comprehensively discover all subspecies/species within a genus"""
        try:
            st.info(f"🔍 **Step 1: Discovering subspecies within genus {genus_name}**")
            
            # Multiple search strategies to capture all variants
            search_queries = [
                f'"{genus_name}"[organism]',
                f'{genus_name}[organism] AND (subspecies OR subsp OR strain OR isolate)',
                f'{genus_name}[organism] AND (f. sp. OR forma OR variety OR var)',
                f'{genus_name}[organism] AND complete[title]'
            ]
            
            all_species_data = {}
            
            for query in search_queries:
                try:
                    st.write(f"Searching: {query}")
                    seq_ids = self.ncbi.search_sequences(
                        query, 
                        database="nucleotide", 
                        max_results=500  # Cast wider net
                    )
                    
                    st.write(f"Found {len(seq_ids)} sequences for this query")
                    
                    # Sample sequences to extract organism names
                    sample_size = min(100, len(seq_ids))
                    for seq_id in seq_ids[:sample_size]:
                        try:
                            seq_info = self.ncbi.fetch_sequence_info(seq_id)
                            if seq_info and seq_info.get('organism'):
                                organism = seq_info['organism']
                                
                                # Clean and standardize organism name
                                clean_organism = self.standardize_organism_name(organism, genus_name)
                                
                                if clean_organism:
                                    if clean_organism not in all_species_data:
                                        all_species_data[clean_organism] = {
                                            'sequence_ids': [],
                                            'sample_info': seq_info
                                        }
                                    all_species_data[clean_organism]['sequence_ids'].append(seq_id)
                        
                        except Exception as e:
                            continue  # Skip problematic sequences
                            
                except Exception as e:
                    st.warning(f"Search query failed: {query} - {e}")
                    continue
            
            # Filter and rank subspecies by available sequence count
            filtered_subspecies = {}
            for organism, data in all_species_data.items():
                seq_count = len(data['sequence_ids'])
                if seq_count >= 5:  # Minimum threshold for inclusion
                    filtered_subspecies[organism] = {
                        'sequence_count': seq_count,
                        'sample_info': data['sample_info'],
                        'sequence_ids': data['sequence_ids'][:self.max_sequences_per_subspecies]
                    }
            
            # Sort by sequence availability
            sorted_subspecies = dict(sorted(
                filtered_subspecies.items(), 
                key=lambda x: x[1]['sequence_count'], 
                reverse=True
            ))
            
            return dict(list(sorted_subspecies.items())[:max_subspecies])
            
        except Exception as e:
            st.error(f"Error discovering subspecies: {e}")
            return {}
    
    def standardize_organism_name(self, organism_name, genus_name):
        """Standardize organism names and filter for the target genus"""
        if not organism_name or genus_name.lower() not in organism_name.lower():
            return None
        
        # Remove common suffixes and standardize
        organism = organism_name.strip()
        
        # Remove strain identifiers but keep subspecies info
        import re
        
        # Patterns to preserve subspecies information
        subspecies_patterns = [
            r'f\.\s*sp\.\s*\w+',  # forma specialis
            r'subsp\.\s*\w+',     # subspecies
            r'var\.\s*\w+',       # variety
            r'pv\.\s*\w+'         # pathovar
        ]
        
        clean_name = organism
        
        # Extract subspecies information
        subspecies_info = ""
        for pattern in subspecies_patterns:
            match = re.search(pattern, organism, re.IGNORECASE)
            if match:
                subspecies_info = " " + match.group()
                break
        
        # Get base species name (first two words)
        words = organism.split()
        if len(words) >= 2 and words[0].lower() == genus_name.lower():
            base_name = f"{words[0]} {words[1]}"
            return base_name + subspecies_info
        
        return None
    
    def fetch_representative_sequences(self, subspecies_data, target_count=75, progress_callback=None):
        """Fetch representative sequences for each subspecies with geographic diversity"""
        st.info(f"🧬 **Step 2: Collecting sequences from {len(subspecies_data)} subspecies**")
        results = {}
        
        total_subspecies = len(subspecies_data)
        for current, (subspecies, data) in enumerate(subspecies_data.items(), 1):
            if progress_callback:
                progress_callback(current, total_subspecies)
            
            st.write(f"Fetching sequences for {subspecies}...")
            
            sequence_ids = data['sequence_ids']
            target_fetch = min(target_count, len(sequence_ids))
            
            # Implement sampling strategy for diversity
            sampled_ids = self.diverse_sampling(sequence_ids, target_fetch)
            
            sequences = []
            fetch_count = 0
            
            for seq_id in sampled_ids:
                if fetch_count >= target_fetch:
                    break
                    
                try:
                    sequence = self.ncbi.fetch_sequence(seq_id)
                    if sequence and len(sequence) > 200:  # Minimum sequence length
                        # Clean sequence
                        clean_seq = re.sub(r'[^ATGCatgc]', '', sequence.upper())
                        if len(clean_seq) > 200:
                            sequences.append({
                                'id': seq_id,
                                'sequence': clean_seq[:5000],  # Limit to 5kb for performance
                                'length': len(clean_seq)
                            })
                            fetch_count += 1
                
                except Exception as e:
                    continue
            
            if sequences:
                results[subspecies] = {
                    'sequences': sequences,
                    'count': len(sequences),
                    'metadata': data
                }
                st.write(f"  Retrieved {len(sequences)} sequences")
            else:
                st.warning(f"  No valid sequences retrieved for {subspecies}")
        
        return results
    
    def diverse_sampling(self, sequence_ids, target_count):
        """Sample sequence IDs to maximize diversity"""
        if len(sequence_ids) <= target_count:
            return sequence_ids
        
        # Simple strategy: take evenly spaced samples
        step = len(sequence_ids) // target_count
        sampled = []
        
        for i in range(0, len(sequence_ids), step):
            sampled.append(sequence_ids[i])
            if len(sampled) >= target_count:
                break
        
        return sampled
    
    def perform_multiple_sequence_alignment(self, subspecies_sequences, progress_callback=None):
        """Align sequences across all subspecies to find conserved regions"""
        st.info("🧮 **Step 3: Analyzing sequence conservation across subspecies**")
        
        # Combine all sequences from all subspecies
        all_sequences = []
        sequence_metadata = []
        
        for subspecies, data in subspecies_sequences.items():
            for seq_data in data['sequences']:
                all_sequences.append(seq_data['sequence'])
                sequence_metadata.append({
                    'subspecies': subspecies,
                    'seq_id': seq_data['id'],
                    'length': seq_data['length']
                })
        
        if len(all_sequences) < 10:
            st.warning("Insufficient sequences for robust alignment analysis")
            return {}
        
        st.write(f"Analyzing {len(all_sequences)} sequences across {len(subspecies_sequences)} subspecies")
        
        # Find conserved regions using sliding window approach
        conserved_regions = self.find_conserved_regions_advanced(
            all_sequences, 
            sequence_metadata, 
            subspecies_sequences,
            progress_callback=progress_callback
        )
        
        return conserved_regions
    
    def find_conserved_regions_advanced(self, sequences, metadata, subspecies_data, progress_callback=None):
        """Advanced conserved region analysis with subspecies representation"""
        
        # Find the minimum length across all sequences
        min_length = min(len(seq) for seq in sequences)
        window_size = 150
        step_size = 25
        
        conserved_regions = []
        positions = range(0, min_length - window_size, step_size)
        
        st.write(f"Analyzing {len(positions)} windows of size {window_size} bp")
        
        # Use provided progress callback or create default progress bar
        if progress_callback:
            for i, start_pos in enumerate(positions):
                progress_callback(i + 1, len(positions))
                
                end_pos = start_pos + window_size
                window_sequences = [seq[start_pos:end_pos] for seq in sequences]
                
                # Calculate conservation metrics
                conservation_metrics = self.calculate_conservation_metrics(
                    window_sequences, 
                    metadata, 
                    subspecies_data
                )
                
                # Check if this region meets conservation criteria
                if (conservation_metrics['overall_conservation'] >= self.conservation_threshold and
                    conservation_metrics['subspecies_representation'] >= 0.8):  # 80% of subspecies represented
                    
                    conserved_regions.append({
                        'start': start_pos,
                        'end': end_pos,
                        'sequence': window_sequences[0],  # Representative sequence
                        'conservation_score': conservation_metrics['overall_conservation'],
                        'subspecies_coverage': conservation_metrics['subspecies_representation'],
                        'sequence_count': len(window_sequences),
                        'gc_content': self.calculate_gc_content(window_sequences[0]),
                        'complexity_score': conservation_metrics['complexity_score']
                    })
        else:
            # Fallback to original progress bar implementation
            progress_bar = st.progress(0)
            
            for i, start_pos in enumerate(positions):
                progress_bar.progress((i + 1) / len(positions))
                
                end_pos = start_pos + window_size
                window_sequences = [seq[start_pos:end_pos] for seq in sequences]
                
                # Calculate conservation metrics
                conservation_metrics = self.calculate_conservation_metrics(
                    window_sequences, 
                    metadata, 
                    subspecies_data
                )
                
                # Check if this region meets conservation criteria
                if (conservation_metrics['overall_conservation'] >= self.conservation_threshold and
                    conservation_metrics['subspecies_representation'] >= 0.8):  # 80% of subspecies represented
                    
                    conserved_regions.append({
                        'start': start_pos,
                        'end': end_pos,
                        'sequence': window_sequences[0],  # Representative sequence
                        'conservation_score': conservation_metrics['overall_conservation'],
                        'subspecies_coverage': conservation_metrics['subspecies_representation'],
                        'sequence_count': len(window_sequences),
                        'gc_content': self.calculate_gc_content(window_sequences[0]),
                        'complexity_score': conservation_metrics['complexity_score']
                    })
            
            progress_bar.progress(1.0)
        
        # Merge overlapping regions
        merged_regions = self.merge_overlapping_regions_advanced(conserved_regions)
        
        st.write(f"Found {len(merged_regions)} conserved regions after merging")
        
        return merged_regions
    
    def calculate_conservation_metrics(self, window_sequences, metadata, subspecies_data):
        """Calculate comprehensive conservation metrics"""
        
        # Overall sequence conservation
        conservation_scores = []
        for i in range(len(window_sequences)):
            for j in range(i + 1, len(window_sequences)):
                similarity = self.calculate_sequence_similarity(
                    window_sequences[i], 
                    window_sequences[j]
                )
                conservation_scores.append(similarity)
        
        overall_conservation = np.mean(conservation_scores) if conservation_scores else 0
        
        # Subspecies representation analysis
        represented_subspecies = set()
        for i, seq_meta in enumerate(metadata):
            if i < len(window_sequences):
                represented_subspecies.add(seq_meta['subspecies'])
        
        subspecies_representation = len(represented_subspecies) / len(subspecies_data)
        
        # Sequence complexity (avoid low-complexity regions)
        representative_seq = window_sequences[0]
        complexity_score = self.calculate_sequence_complexity(representative_seq)
        
        return {
            'overall_conservation': overall_conservation,
            'subspecies_representation': subspecies_representation,
            'complexity_score': complexity_score,
            'represented_subspecies': list(represented_subspecies)
        }
    
    def calculate_sequence_similarity(self, seq1, seq2):
        """Calculate pairwise sequence similarity"""
        if len(seq1) != len(seq2):
            return 0.0
        
        matches = sum(1 for a, b in zip(seq1.upper(), seq2.upper()) if a == b)
        return matches / len(seq1)
    
    def calculate_sequence_complexity(self, sequence):
        """Calculate sequence complexity to avoid low-complexity regions"""
        if not sequence:
            return 0
        
        # Count unique dinucleotides
        dinucleotides = set()
        for i in range(len(sequence) - 1):
            dinucleotides.add(sequence[i:i+2])
        
        # Complexity score based on dinucleotide diversity
        max_possible = min(16, len(sequence) - 1)  # 16 possible dinucleotides
        complexity = len(dinucleotides) / max_possible if max_possible > 0 else 0
        
        return complexity
    
    def calculate_gc_content(self, sequence):
        """Calculate GC content of sequence"""
        if not sequence:
            return 0
        gc_count = sequence.upper().count('G') + sequence.upper().count('C')
        return (gc_count / len(sequence)) * 100
    
    def merge_overlapping_regions_advanced(self, regions, min_gap=50):
        """Merge overlapping conserved regions with quality scoring"""
        if not regions:
            return []
        
        # Sort by position
        sorted_regions = sorted(regions, key=lambda x: x['start'])
        merged = []
        
        current_region = sorted_regions[0].copy()
        
        for next_region in sorted_regions[1:]:
            # If regions overlap or are close
            gap = next_region['start'] - current_region['end']
            
            if gap <= min_gap:
                # Merge regions - take best metrics
                current_region['end'] = next_region['end']
                current_region['conservation_score'] = max(
                    current_region['conservation_score'],
                    next_region['conservation_score']
                )
                current_region['subspecies_coverage'] = max(
                    current_region['subspecies_coverage'],
                    next_region['subspecies_coverage']
                )
                current_region['sequence_count'] = max(
                    current_region['sequence_count'],
                    next_region['sequence_count']
                )
            else:
                # Add current region if it meets quality thresholds
                if (current_region['end'] - current_region['start'] >= 100 and
                    current_region['complexity_score'] >= 0.3):
                    merged.append(current_region)
                current_region = next_region.copy()
        
        # Add final region
        if (current_region['end'] - current_region['start'] >= 100 and
            current_region['complexity_score'] >= 0.3):
            merged.append(current_region)
        
        return merged
    
    def compare_against_other_taxa(self, conserved_regions, target_genus, 
                                  test_genera=True, test_species=True, 
                                  comparison_genera=None, comparison_species=None,
                                  progress_callback=None):
        """Compare conserved regions against other genera and/or species for specificity"""
        
        st.info("🎯 **Step 4: Testing specificity against related organisms**")
        
        specificity_results = {
            'genus_results': {},
            'species_results': {},
            'combined_specific_regions': []
        }
        
        # Test against genera if requested
        if test_genera:
            st.write("**Step 4a: Genus-Level Specificity Testing**")
            if not comparison_genera:
                comparison_genera = self.get_related_genera(target_genus)
            
            st.write(f"Testing against genera: {', '.join(comparison_genera)}")
            genus_specific_regions = self._test_specificity_against_taxa(
                conserved_regions, comparison_genera, "genus", progress_callback
            )
            specificity_results['genus_results'] = {
                'tested_against': comparison_genera,
                'specific_regions': genus_specific_regions
            }
            
            st.write(f"Genus-specific regions found: {len(genus_specific_regions)}")
        
        # Test against species if requested
        if test_species:
            st.write("**Step 4b: Species-Level Specificity Testing**")
            if not comparison_species:
                comparison_species = self.get_related_species(target_genus)
            
            st.write(f"Testing against species: {', '.join(comparison_species)}")
            species_specific_regions = self._test_specificity_against_taxa(
                conserved_regions, comparison_species, "species", progress_callback
            )
            specificity_results['species_results'] = {
                'tested_against': comparison_species,
                'specific_regions': species_specific_regions
            }
            
            st.write(f"Species-specific regions found: {len(species_specific_regions)}")
        
        # Combine results based on what was tested
        if test_genera and test_species:
            # Regions must pass both genus and species specificity
            genus_positions = {(r['start'], r['end']) for r in genus_specific_regions}
            species_positions = {(r['start'], r['end']) for r in species_specific_regions}
            common_positions = genus_positions.intersection(species_positions)
            
            combined_regions = []
            for region in genus_specific_regions:
                if (region['start'], region['end']) in common_positions:
                    # Find corresponding species result
                    species_region = next(
                        (r for r in species_specific_regions 
                         if r['start'] == region['start'] and r['end'] == region['end']), 
                        None
                    )
                    if species_region:
                        # Combine specificity scores
                        combined_region = region.copy()
                        combined_region['genus_specificity'] = region['specificity_score']
                        combined_region['species_specificity'] = species_region['specificity_score']
                        combined_region['combined_specificity'] = min(
                            region['specificity_score'], 
                            species_region['specificity_score']
                        )
                        combined_region['specificity_score'] = combined_region['combined_specificity']
                        combined_regions.append(combined_region)
            
            specificity_results['combined_specific_regions'] = combined_regions
            st.write(f"Final regions passing both tests: {len(combined_regions)}")
            
        elif test_genera:
            specificity_results['combined_specific_regions'] = genus_specific_regions
        elif test_species:
            specificity_results['combined_specific_regions'] = species_specific_regions
        else:
            specificity_results['combined_specific_regions'] = conserved_regions
        
        return specificity_results
    
    def _test_specificity_against_taxa(self, conserved_regions, comparison_taxa, test_type, progress_callback=None):
        """Helper method to test specificity against a list of taxa"""
        specific_regions = []
        
        total_regions = len(conserved_regions)
        
        # Initialize progress bar outside the loop if no callback provided
        if not progress_callback:
            progress_bar = st.progress(0)
        
        for idx, region in enumerate(conserved_regions):
            if progress_callback:
                progress_callback(idx + 1, total_regions)
            else:
                progress_bar.progress((idx + 1) / total_regions)
            
            st.write(f"Testing {test_type} specificity for region {region['start']}-{region['end']}")
            
            # Test against each comparison taxon
            specificity_scores = []
            
            for taxon in comparison_taxa:
                try:
                    specificity = self.test_region_specificity(
                        region['sequence'], 
                        taxon
                    )
                    specificity_scores.append(specificity)
                    
                except Exception as e:
                    st.warning(f"Could not test against {taxon}: {e}")
                    continue
            
            if specificity_scores:
                avg_specificity = 1.0 - np.mean(specificity_scores)  # Convert to specificity
                
                if avg_specificity >= self.specificity_threshold:
                    region_copy = region.copy()
                    region_copy['specificity_score'] = avg_specificity
                    region_copy['tested_against'] = comparison_taxa
                    region_copy['test_type'] = test_type
                    specific_regions.append(region_copy)
                    st.write(f"  ✓ Region passed {test_type} specificity test (score: {avg_specificity:.3f})")
                else:
                    st.write(f"  ✗ Region failed {test_type} specificity test (score: {avg_specificity:.3f})")
        
        # Complete progress bar only if we created one
        if not progress_callback:
            progress_bar.progress(1.0)
        
        return specific_regions
    
    def get_related_genera(self, target_genus):
        """Get related genera for specificity testing - expanded for agricultural pathogens"""
        related_genera_map = {
            'Fusarium': ['Trichoderma', 'Aspergillus', 'Penicillium', 'Verticillium', 'Rhizoctonia', 'Pythium', 'Alternaria', 'Botrytis'],
            'Botrytis': ['Sclerotinia', 'Monilinia', 'Alternaria', 'Fusarium', 'Rhizoctonia'],
            'Pythium': ['Phytophthora', 'Peronospora', 'Saprolegnia', 'Fusarium', 'Rhizoctonia'],
            'Alternaria': ['Stemphylium', 'Ulocladium', 'Botrytis', 'Fusarium', 'Cladosporium'],
            'Rhizoctonia': ['Thanatephorus', 'Ceratobasidium', 'Sclerotium', 'Fusarium', 'Pythium'],
            'Aspergillus': ['Penicillium', 'Trichoderma', 'Fusarium', 'Alternaria'],
            'Trichoderma': ['Aspergillus', 'Penicillium', 'Fusarium'],
            'Phytophthora': ['Pythium', 'Peronospora', 'Fusarium', 'Alternaria']
        }
        
        return related_genera_map.get(target_genus, ['Aspergillus', 'Penicillium', 'Trichoderma', 'Fusarium', 'Alternaria'])
    
    def get_related_species(self, target_genus):
        """Get related species for specificity testing"""
        related_species_map = {
            'Fusarium': [
                'Alternaria alternata', 'Botrytis cinerea', 'Rhizoctonia solani',
                'Pythium ultimum', 'Sclerotinia sclerotiorum', 'Verticillium dahliae',
                'Colletotrichum gloeosporioides', 'Phytophthora infestans'
            ],
            'Botrytis': [
                'Sclerotinia sclerotiorum', 'Alternaria alternata', 'Fusarium oxysporum',
                'Rhizoctonia solani', 'Colletotrichum gloeosporioides'
            ],
            'Pythium': [
                'Phytophthora infestans', 'Fusarium oxysporum', 'Rhizoctonia solani',
                'Alternaria alternata', 'Verticillium dahliae'
            ],
            'Alternaria': [
                'Fusarium oxysporum', 'Botrytis cinerea', 'Stemphylium vesicarium',
                'Cladosporium cladosporioides', 'Rhizoctonia solani'
            ],
            'Rhizoctonia': [
                'Fusarium oxysporum', 'Pythium ultimum', 'Sclerotinia sclerotiorum',
                'Alternaria alternata', 'Botrytis cinerea'
            ],
            'Phytophthora': [
                'Pythium ultimum', 'Fusarium oxysporum', 'Alternaria alternata',
                'Rhizoctonia solani', 'Verticillium dahliae'
            ]
        }
        
        return related_species_map.get(target_genus, [
            'Fusarium oxysporum', 'Alternaria alternata', 'Botrytis cinerea',
            'Rhizoctonia solani', 'Pythium ultimum'
        ])
    
    def test_region_specificity(self, target_sequence, comparison_genus, sample_size=20):
        """Test a region's specificity against another genus"""
        try:
            # Search for sequences from comparison genus
            query = f'"{comparison_genus}"[organism]'
            comparison_ids = self.ncbi.search_sequences(
                query, 
                database="nucleotide", 
                max_results=sample_size
            )
            
            if not comparison_ids:
                return 0.5  # Neutral score if no sequences found
            
            similarity_scores = []
            
            # Test against a sample of comparison sequences
            for seq_id in comparison_ids[:sample_size]:
                try:
                    comparison_seq = self.ncbi.fetch_sequence(seq_id)
                    if comparison_seq and len(comparison_seq) > len(target_sequence):
                        # Find best local alignment (simplified)
                        best_similarity = self.find_best_local_similarity(
                            target_sequence, 
                            comparison_seq
                        )
                        similarity_scores.append(best_similarity)
                        
                except Exception:
                    continue
            
            return np.mean(similarity_scores) if similarity_scores else 0.5
            
        except Exception as e:
            return 0.5  # Default neutral score
    
    def find_best_local_similarity(self, query_seq, target_seq):
        """Find best local similarity between sequences"""
        query_len = len(query_seq)
        best_similarity = 0
        
        # Slide query across target sequence
        for i in range(len(target_seq) - query_len + 1):
            window = target_seq[i:i + query_len]
            similarity = self.calculate_sequence_similarity(query_seq, window)
            best_similarity = max(best_similarity, similarity)
        
        return best_similarity

class T7PrimerDesigner:
    """Handles T7 promoter integration for expression primers"""
    
    def __init__(self):
        self.t7_promoter = "TAATACGACTCACTATAGGG"
        self.t7_terminator = "GCTAGTTATTGCTCAGCGG"
        self.kozak_sequences = {
            "optimal": "GCCACC",
            "strong": "GCCGCC", 
            "moderate": "ACCACC"
        }
        self.start_codons = ["ATG", "GTG", "TTG"]
        self.stop_codons = ["TAA", "TAG", "TGA"]
    
    def find_orfs(self, sequence, min_length=150):
        """Find open reading frames in the sequence"""
        orfs = []
        
        for frame in range(3):
            for start_codon in self.start_codons:
                start_pos = frame
                while start_pos < len(sequence) - 3:
                    codon = sequence[start_pos:start_pos + 3]
                    if codon == start_codon:
                        # Look for stop codon
                        for stop_pos in range(start_pos + 3, len(sequence) - 2, 3):
                            stop_codon = sequence[stop_pos:stop_pos + 3]
                            if stop_codon in self.stop_codons:
                                orf_length = stop_pos - start_pos + 3
                                if orf_length >= min_length:
                                    orfs.append({
                                        'start': start_pos,
                                        'end': stop_pos + 3,
                                        'length': orf_length,
                                        'frame': frame,
                                        'sequence': sequence[start_pos:stop_pos + 3]
                                    })
                                break
                    start_pos = sequence.find(start_codon, start_pos + 1)
                    if start_pos == -1:
                        break
        
        return sorted(orfs, key=lambda x: x['length'], reverse=True)
    
    def design_t7_expression_primers(self, sequence, target_orf=None, 
                                   kozak_type="optimal", add_his_tag=False,
                                   restriction_sites=None):
        """Design primers for T7 expression with promoter integration"""
        
        if target_orf is None:
            # Find the longest ORF
            orfs = self.find_orfs(sequence)
            if not orfs:
                return None, "No suitable ORFs found"
            target_orf = orfs[0]
        
        # Design forward primer with T7 promoter
        kozak = self.kozak_sequences.get(kozak_type, self.kozak_sequences["optimal"])
        
        # Start primer components
        forward_components = [self.t7_promoter]
        
        # Add restriction site if specified
        if restriction_sites and 'forward' in restriction_sites:
            forward_components.append(restriction_sites['forward'])
        
        # Add Kozak sequence
        forward_components.append(kozak)
        
        # Add start of ORF (usually includes ATG)
        orf_start = max(0, target_orf['start'] - 20)  # Include some upstream context
        orf_binding = sequence[orf_start:target_orf['start'] + 20]
        forward_components.append(orf_binding)
        
        forward_primer = "".join(forward_components)
        
        # Design reverse primer
        reverse_components = []
        
        # Add restriction site if specified
        if restriction_sites and 'reverse' in restriction_sites:
            reverse_components.append(restriction_sites['reverse'])
        
        # Add His tag if requested
        if add_his_tag:
            his_tag_reverse = "GTGATGGTGATGGTG"  # Reverse complement of His6 tag
            reverse_components.append(his_tag_reverse)
        
        # Add end of ORF
        orf_end_start = max(target_orf['end'] - 20, target_orf['start'])
        orf_end = min(target_orf['end'] + 20, len(sequence))
        orf_binding_rev = self.reverse_complement(sequence[orf_end_start:orf_end])
        reverse_components.append(orf_binding_rev)
        
        reverse_primer = "".join(reverse_components)
        
        return {
            'forward_primer': forward_primer,
            'reverse_primer': reverse_primer,
            'target_orf': target_orf,
            'expression_features': {
                't7_promoter': True,
                'kozak_sequence': kozak_type,
                'his_tag': add_his_tag,
                'restriction_sites': restriction_sites or {}
            }
        }, None
    
    def reverse_complement(self, sequence):
        """Get reverse complement of DNA sequence"""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        return "".join(complement.get(base, base) for base in reversed(sequence.upper()))
    
    def validate_expression_construct(self, forward_primer, reverse_primer, template):
        """Validate that the expression construct will work"""
        issues = []
        
        # Check for T7 promoter
        if self.t7_promoter not in forward_primer:
            issues.append("T7 promoter not found in forward primer")
        
        # Check for start codon after Kozak
        kozak_pos = -1
        for kozak in self.kozak_sequences.values():
            if kozak in forward_primer:
                kozak_pos = forward_primer.find(kozak)
                break
        
        if kozak_pos >= 0:
            remaining_seq = forward_primer[kozak_pos + 6:]  # After Kozak
            if not any(start in remaining_seq[:10] for start in self.start_codons):
                issues.append("No start codon found near Kozak sequence")
        
        # Check primer lengths
        if len(forward_primer) > 80:
            issues.append("Forward primer very long (>80 bp) - may have synthesis issues")
        if len(reverse_primer) > 80:
            issues.append("Reverse primer very long (>80 bp) - may have synthesis issues")
        
        return {
            'valid': len(issues) == 0,
            'issues': issues,
            'warnings': []
        }

# Streamlit App Functions
def init_session_state():
    """Initialize session state variables"""
    if 'primers_designed' not in st.session_state:
        st.session_state.primers_designed = []
    if 'current_sequence' not in st.session_state:
        st.session_state.current_sequence = ""
    if 'sequence_info' not in st.session_state:
        st.session_state.sequence_info = {}
    if 'search_results' not in st.session_state:
        st.session_state.search_results = None
    if 'database_used' not in st.session_state:
        st.session_state.database_used = None

def get_organism_suggestions():
    """Get agricultural pest and pathogen suggestions for the search"""
    return {
        "Arthropods": {
            "🕷️ Mites": ["Aculops lycopersici", "Tetranychus urticae", "Polyphagotarsonemus latus"],
            "🐛 Other, walking": ["Pemphigus betae", "Thrips tabaci", "Liriomyza trifolii", "Planococcus citri"],
            "🦟 Other, flying": ["Bemisia tabaci", "Bradysia impatiens", "Diaphorina citri"]
        },
        "Fungi": {
            "🍄 Fusarium species": ["Fusarium oxysporum", "Fusarium solani", "Fusarium proliferatum", "Fusarium"],
            "🍄 Botrytis species": ["Botrytis cinerea", "Botrytis"],
            "🍄 Powdery mildew": ["Golovinomyces ambrosiae", "Golovinomyces"]
        },
        "Pseudofungi/Oomycetes": {
            "💧 Water molds": ["Pythium myriotylum", "Pythium"]
        },
        "Viruses": {
            "🦠 Plant viruses": ["Beet curly top virus", "Alfalfa mosaic virus", "Arabis mosaic virus", 
                               "Lettuce chlorosis virus", "Cannabis cryptic virus", "Tomato ringspot virus",
                               "Tomato mosaic virus", "Tobacco mosaic virus"]
        },
        "Viroids": {
            "🧬 RNA pathogens": ["Hop latent viroid"]
        }
    }

def create_primer_visualization(primers: List[PrimerPair]):
    """Create interactive visualizations for primer pairs"""
    if not primers:
        return None
    
    # Create subplot with secondary y-axis
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=('Melting Temperatures', 'GC Content Distribution', 
                       'Product Sizes', 'Penalty Scores'),
        specs=[[{"secondary_y": False}, {"secondary_y": False}],
               [{"secondary_y": False}, {"secondary_y": False}]]
    )
    
    # Data preparation
    primer_nums = list(range(1, len(primers) + 1))
    forward_tms = [p.forward_tm for p in primers]
    reverse_tms = [p.reverse_tm for p in primers]
    forward_gcs = [p.gc_content_f for p in primers]
    reverse_gcs = [p.gc_content_r for p in primers]
    product_sizes = [p.product_size for p in primers]
    penalties = [p.penalty for p in primers]
    
    # Melting temperatures
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
    
    # GC Content
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
    
    # Product sizes
    fig.add_trace(
        go.Scatter(x=primer_nums, y=product_sizes, name='Product Size', 
                  line=dict(color='green'), mode='lines+markers'),
        row=2, col=1
    )
    
    # Penalty scores
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

def create_sequence_diagram(sequence: str, primers: List[PrimerPair], selected_primer: int = 0):
    """Create a sequence diagram showing primer binding sites"""
    if not primers or selected_primer >= len(primers):
        return None
    
    primer = primers[selected_primer]
    seq_len = len(sequence)
    
    # Create figure
    fig = go.Figure()
    
    # Add sequence as background
    fig.add_shape(
        type="rect",
        x0=0, y0=0.4, x1=seq_len, y1=0.6,
        fillcolor="lightgray",
        line=dict(color="black", width=1),
    )
    
    # Add forward primer
    fig.add_shape(
        type="rect",
        x0=primer.forward_start, y0=0.6, 
        x1=primer.forward_start + len(primer.forward_seq), y1=0.8,
        fillcolor="blue",
        line=dict(color="darkblue", width=2),
    )
    
    # Add reverse primer
    fig.add_shape(
        type="rect",
        x0=primer.reverse_start - len(primer.reverse_seq) + 1, y0=0.2,
        x1=primer.reverse_start + 1, y1=0.4,
        fillcolor="red",
        line=dict(color="darkred", width=2),
    )
    
    # Add annotations
    fig.add_annotation(
        x=primer.forward_start + len(primer.forward_seq)/2, y=0.7,
        text=f"Forward<br>Tm: {primer.forward_tm:.1f}°C",
        showarrow=True, arrowhead=2, arrowcolor="blue"
    )
    
    fig.add_annotation(
        x=primer.reverse_start - len(primer.reverse_seq)/2, y=0.3,
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

def export_to_excel(primers: List[PrimerPair]) -> bytes:
    """Export primer results to Excel format"""
    data = []
    for i, primer in enumerate(primers):
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
            'Penalty_Score': round(primer.penalty, 4)
        })
    
    df = pd.DataFrame(data)
    output = io.BytesIO()
    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        df.to_excel(writer, sheet_name='Primer_Results', index=False)
    
    return output.getvalue()

def create_species_specific_visualization(sequence, conserved_regions, primers=None):
    """Create visualization showing species-specific regions and primer locations"""
    if not conserved_regions:
        return None
    
    fig = go.Figure()
    
    seq_len = len(sequence)
    
    # Add sequence background
    fig.add_shape(
        type="rect",
        x0=0, y0=0.3, x1=seq_len, y1=0.7,
        fillcolor="lightgray",
        line=dict(color="black", width=1),
    )
    
    # Add conserved regions
    for i, region in enumerate(conserved_regions):
        fig.add_shape(
            type="rect",
            x0=region['start'], y0=0.25, 
            x1=region['end'], y1=0.75,
            fillcolor="green",
            opacity=0.6,
            line=dict(color="darkgreen", width=2),
        )
        
        # Add region annotation
        fig.add_annotation(
            x=(region['start'] + region['end']) / 2, y=0.8,
            text=f"Region {i+1}<br>Cons: {region['conservation_score']:.1%}<br>Spec: {region['specificity_score']:.1%}",
            showarrow=True, arrowhead=2, arrowcolor="green"
        )
    
    # Add primers if provided
    if primers:
        for i, primer in enumerate(primers[:3]):  # Show first 3 primers
            # Forward primer
            fig.add_shape(
                type="rect",
                x0=primer.forward_start, y0=0.75, 
                x1=primer.forward_start + len(primer.forward_seq), y1=0.9,
                fillcolor="blue",
                opacity=0.8,
                line=dict(color="darkblue", width=1),
            )
            
            # Reverse primer
            fig.add_shape(
                type="rect",
                x0=primer.reverse_start - len(primer.reverse_seq) + 1, y0=0.1,
                x1=primer.reverse_start + 1, y1=0.25,
                fillcolor="red",
                opacity=0.8,
                line=dict(color="darkred", width=1),
            )
    
    fig.update_layout(
        title="Species-Specific Regions and Primer Locations",
        xaxis_title="Sequence Position (bp)",
        yaxis=dict(range=[0, 1], showticklabels=False),
        height=400,
        showlegend=False
    )
    
    return fig

def create_t7_construct_diagram(sequence, t7_results):
    """Create diagram showing T7 expression construct"""
    if not t7_results:
        return None
    
    fig = go.Figure()
    
    target_orf = t7_results['target_orf']
    features = t7_results['expression_features']
    
    # Sequence background
    seq_len = len(sequence)
    fig.add_shape(
        type="rect",
        x0=0, y0=0.4, x1=seq_len, y1=0.6,
        fillcolor="lightgray",
        line=dict(color="black", width=1),
    )
    
    # Target ORF
    fig.add_shape(
        type="rect",
        x0=target_orf['start'], y0=0.35, 
        x1=target_orf['end'], y1=0.65,
        fillcolor="green",
        line=dict(color="darkgreen", width=2),
    )
    
    # T7 promoter (conceptual position)
    fig.add_shape(
        type="rect",
        x0=max(0, target_orf['start'] - 30), y0=0.7, 
        x1=target_orf['start'], y1=0.8,
        fillcolor="blue",
        line=dict(color="darkblue", width=2),
    )
    
    # Add annotations
    fig.add_annotation(
        x=target_orf['start'] - 15, y=0.75,
        text="T7 Promoter",
        showarrow=True, arrowhead=2, arrowcolor="blue"
    )
    
    fig.add_annotation(
        x=(target_orf['start'] + target_orf['end']) / 2, y=0.5,
        text=f"Target ORF<br>{target_orf['length']} bp",
        showarrow=False
    )
    
    if features.get('his_tag'):
        fig.add_annotation(
            x=target_orf['end'] + 10, y=0.3,
            text="His6 Tag",
            showarrow=True, arrowhead=2, arrowcolor="orange"
        )
    
    fig.update_layout(
        title="T7 Expression Construct Design",
        xaxis_title="Sequence Position (bp)",
        yaxis=dict(range=[0, 1], showticklabels=False),
        height=300,
        showlegend=False
    )
    
    return fig

def create_conservation_heatmap(conserved_regions, subspecies_data):
    """Create heatmap showing conservation across subspecies"""
    if not conserved_regions:
        return None
    
    import plotly.express as px
    
    # Prepare data for heatmap
    heatmap_data = []
    for i, region in enumerate(conserved_regions[:20]):  # Limit to top 20 regions
        row = {
            'Region': f"Region {i+1}\n({region['start']}-{region['end']})",
            'Conservation Score': region['conservation_score'],
            'Subspecies Coverage': region['subspecies_coverage'],
            'GC Content': region['gc_content'] / 100,  # Normalize to 0-1
            'Complexity': region['complexity_score']
        }
        heatmap_data.append(row)
    
    df = pd.DataFrame(heatmap_data)
    
    if df.empty:
        return None
    
    fig = px.imshow(
        df.set_index('Region').T,
        aspect="auto",
        color_continuous_scale="RdYlBu_r",
        title="Conservation Metrics Across Top Regions"
    )
    
    fig.update_layout(
        height=400,
        xaxis_title="Genomic Regions",
        yaxis_title="Conservation Metrics"
    )
    
    return fig

def create_specificity_comparison_chart(specificity_results):
    """Create chart comparing genus vs species specificity"""
    if not specificity_results:
        return None
    
    genus_results = specificity_results.get('genus_results', {})
    species_results = specificity_results.get('species_results', {})
    
    if not genus_results and not species_results:
        return None
    
    fig = go.Figure()
    
    # Create comparison data
    regions = specificity_results['combined_specific_regions']
    
    if not regions:
        return None
    
    x_labels = [f"Region {i+1}" for i in range(len(regions))]
    
    if 'genus_specificity' in regions[0]:
        genus_scores = [r['genus_specificity'] for r in regions]
        fig.add_trace(go.Bar(
            name='Genus Specificity',
            x=x_labels,
            y=genus_scores,
            marker_color='lightblue'
        ))
    
    if 'species_specificity' in regions[0]:
        species_scores = [r['species_specificity'] for r in regions]
        fig.add_trace(go.Bar(
            name='Species Specificity', 
            x=x_labels,
            y=species_scores,
            marker_color='lightcoral'
        ))
    
    fig.update_layout(
        title='Specificity Comparison: Genus vs Species Level',
        xaxis_title='Genomic Regions',
        yaxis_title='Specificity Score',
        barmode='group',
        height=400
    )
    
    return fig

def create_analysis_workflow_diagram():
    """Create a diagram showing the analysis workflow"""
    fig = go.Figure()
    
    # Workflow steps
    steps = [
        "1. Subspecies\nDiscovery",
        "2. Sequence\nCollection", 
        "3. Conservation\nAnalysis",
        "4a. Genus\nSpecificity",
        "4b. Species\nSpecificity",
        "5. Primer\nDesign"
    ]
    
    colors = ['lightblue', 'lightgreen', 'lightyellow', 'lightcoral', 'lightpink', 'lightgray']
    
    for i, (step, color) in enumerate(zip(steps, colors)):
        fig.add_shape(
            type="rect",
            x0=i, y0=0, x1=i+0.8, y1=1,
            fillcolor=color,
            line=dict(color="black", width=2)
        )
        
        fig.add_annotation(
            x=i+0.4, y=0.5,
            text=step,
            showarrow=False,
            font=dict(size=10)
        )
        
        # Add arrows between steps
        if i < len(steps) - 1:
            fig.add_annotation(
                x=i+0.9, y=0.5,
                text="→",
                showarrow=False,
                font=dict(size=20)
            )
    
    fig.update_layout(
        title="Comprehensive Taxonomic Analysis Workflow",
        xaxis=dict(range=[-0.2, len(steps)-0.2], showticklabels=False),
        yaxis=dict(range=[-0.2, 1.2], showticklabels=False),
        height=150,
        showlegend=False
    )
    
    return fig

def main():
    """Main Streamlit application"""
    
    # Initialize session state
    init_session_state()
    
    # Initialize email variable
    email = ""
    
    # Header
    st.title("🧬 Automated Primer Design Tool")
    st.markdown("### Design PCR primers with NCBI database integration")
    
    # Sidebar configuration
    st.sidebar.header("⚙️ Configuration")
    
    # NCBI Configuration
    st.sidebar.subheader("NCBI Settings")
    
    # Make email requirement more prominent
    st.sidebar.info("📧 **Email Required**: NCBI requires a valid email address for database access. This is mandatory for all searches.")
    
    email = st.sidebar.text_input("Email (required for NCBI)", 
                                 placeholder="your.email@example.com",
                                 help="Required by NCBI for database access. Use any valid email address.")
    
    # Add a default email option for testing
    if st.sidebar.button("🚀 Use test email (demo@example.com)", type="secondary"):
        st.session_state.demo_email = "demo@example.com"
        st.rerun()
    
    # Use demo email if set
    if 'demo_email' in st.session_state:
        email = st.session_state.demo_email
        del st.session_state.demo_email
    
    api_key = st.sidebar.text_input("NCBI API Key (optional)", 
                                   type="password",
                                   help="For increased rate limits (not required)")
    
    # Quick start guide (moved after email is defined)
    if not email:
        st.warning("🚀 **Quick Start**: Enter an email address in the sidebar to begin searching for organisms and designing primers. NCBI requires this for database access.")
    
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
        
        # Product size ranges
        st.write("Product size ranges:")
        min_product = st.number_input("Minimum product size", 50, 500, 75)
        max_product = st.number_input("Maximum product size", 200, 2000, 1000)
    
    # Comprehensive Taxonomic Analysis
    st.sidebar.subheader("Comprehensive Analysis")
    enable_comprehensive = st.sidebar.checkbox(
        "Enable comprehensive genus analysis", 
        value=False,
        help="Discover all subspecies and design genus-specific primers"
    )

    if enable_comprehensive:
        with st.sidebar.expander("Comprehensive Analysis Parameters", expanded=True):
            max_subspecies = st.slider(
                "Maximum subspecies to analyze:", 
                5, 50, 15,
                help="Number of subspecies/species to include in analysis"
            )
            
            sequences_per_subspecies = st.slider(
                "Target sequences per subspecies:", 
                20, 150, 75,
                help="Number of sequences to fetch for each subspecies"
            )
            
            conservation_threshold = st.slider(
                "Conservation threshold (%):", 
                85, 98, 90,
                help="Minimum conservation across subspecies"
            ) / 100
            
            specificity_threshold = st.slider(
                "Specificity threshold (%):", 
                70, 90, 75,
                help="Minimum specificity against other taxa"
            ) / 100
            
            st.subheader("Specificity Testing Options")
            test_genus_specificity = st.checkbox(
                "Test genus-level specificity",
                value=True,
                help="Test against related genera (e.g., Fusarium vs Aspergillus)"
            )
            
            test_species_specificity = st.checkbox(
                "Test species-level specificity", 
                value=True,
                help="Test against related species (e.g., Fusarium vs Alternaria alternata)"
            )
            
            if not test_genus_specificity and not test_species_specificity:
                st.warning("⚠️ At least one specificity test should be enabled")
                test_genus_specificity = True  # Force at least one to be true
            
            custom_comparison_genera = st.text_input(
                "Custom comparison genera (comma-separated):",
                placeholder="Trichoderma, Aspergillus, Penicillium",
                help="Leave empty for automatic selection based on target genus"
            )
            
            custom_comparison_species = st.text_input(
                "Custom comparison species (comma-separated):",
                placeholder="Alternaria alternata, Botrytis cinerea",
                help="Leave empty for automatic selection based on target genus"
            )
    
    # T7 Expression System
    st.sidebar.subheader("T7 Expression System")
    enable_t7_design = st.sidebar.checkbox(
        "Design for T7 expression", 
        value=False,
        help="Add T7 promoter and expression elements to primers"
    )

    if enable_t7_design:
        with st.sidebar.expander("T7 Expression Parameters", expanded=True):
            kozak_type = st.selectbox(
                "Kozak sequence strength:",
                ["optimal", "strong", "moderate"],
                help="Affects translation efficiency"
            )
            
            add_his_tag = st.checkbox(
                "Add His6 tag",
                value=True,
                help="Add C-terminal His6 tag for purification"
            )
            
            add_restriction_sites = st.checkbox(
                "Add restriction sites",
                value=False,
                help="Add sites for cloning into vectors"
            )
            
            if add_restriction_sites:
                forward_site = st.text_input(
                    "Forward restriction site:",
                    value="GAATTC",  # EcoRI
                    help="Restriction site for 5' end"
                )
                reverse_site = st.text_input(
                    "Reverse restriction site:",
                    value="AAGCTT",  # HindIII
                    help="Restriction site for 3' end"
                )
            
            min_orf_length = st.slider(
                "Minimum ORF length (bp):",
                100, 1000, 150, 50,
                help="Minimum size for open reading frames"
            )
            
            auto_select_orf = st.checkbox(
                "Auto-select longest ORF",
                value=True,
                help="Automatically target the longest open reading frame"
            )
    
    # Create custom parameters
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
        'PRIMER_PRODUCT_SIZE_RANGE': [[min_product, max_product]]
    }
    
    # Main content area
    tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs([
        "📝 Input", 
        "🔬 Results", 
        "📊 Analysis", 
        "🧬 Conservation", 
        "🎯 Specificity", 
        "💾 Export"
    ])
    
    with tab1:
        st.header("Sequence Input")
        
        input_method = st.radio(
            "Choose input method:",
            ["Organism Name", "GenBank ID", "NCBI Search", "Direct Sequence", "Upload File"]
        )
        
        if input_method == "Organism Name":
            st.subheader("Search by Organism")
            
            if not email:
                st.error("❌ **Email Required**: Please enter an email address in the sidebar first to search for organisms.")
            else:
                st.info("💡 **Tip:** Enter the scientific name (e.g., 'Fusarium oxysporum') for best results. You can also search for specific genes within an organism using the advanced options below.")
            
            # Display stored search results if they exist
            if st.session_state.search_results:
                st.subheader("Previous Search Results")
                genome_df = pd.DataFrame(st.session_state.search_results)
                st.dataframe(genome_df, use_container_width=True)
                
                # Let user select a sequence from stored results
                database_used = st.session_state.database_used
                if database_used == "genome":
                    select_text = "Select a genome assembly to design primers for:"
                    button_text = "Design Primers for Selected Assembly"
                else:
                    select_text = "Select a nucleotide sequence to design primers for:"
                    button_text = "Design Primers for Selected Sequence"
                
                selected_genome = st.selectbox(
                    select_text,
                    range(len(st.session_state.search_results)),
                    format_func=lambda x: f"{st.session_state.search_results[x]['ID']} - {st.session_state.search_results[x]['Description'][:100]}..."
                )
                
                if st.button(button_text, type="primary"):
                    with st.spinner("Fetching sequence and designing primers..."):
                        sequence_id = st.session_state.search_results[selected_genome]['ID']
                        
                        try:
                            # Initialize NCBI connector for stored results
                            ncbi = NCBIConnector(email, api_key)
                            designer = PrimerDesigner()
                            
                            if database_used == "genome":
                                # For genome assemblies, search for nucleotide sequences
                                assembly_query = f'"{sequence_id}"[Assembly]'
                                nucleotide_ids = ncbi.search_sequences(assembly_query, database="nucleotide", max_results=10)
                                
                                if nucleotide_ids:
                                    st.info(f"Found {len(nucleotide_ids)} sequences from this assembly. Using the first one for primer design.")
                                    sequence_id = nucleotide_ids[0]
                                else:
                                    st.warning("No nucleotide sequences found for this assembly.")
                                    return
                            
                            # Fetch the sequence
                            st.write(f"🔍 Fetching sequence {sequence_id} from {database_used} database...")
                            
                            if database_used == "nucleotide":
                                # For nucleotide sequences, fetch directly
                                st.write("📡 Using nucleotide database for direct fetch...")
                                sequence = ncbi.fetch_sequence(sequence_id, database="nucleotide")
                                seq_info = ncbi.fetch_sequence_info(sequence_id, database="nucleotide")
                            else:
                                # For genome assemblies, use the nucleotide ID we found
                                st.write("📡 Using default database for genome assembly...")
                                sequence = ncbi.fetch_sequence(sequence_id)
                                seq_info = ncbi.fetch_sequence_info(sequence_id)
                            
                            if sequence:
                                st.write(f"✅ Successfully fetched sequence: {len(sequence)} bp")
                                st.write(f"📝 First 100 characters: {sequence[:100]}...")
                            else:
                                st.write("❌ Failed to fetch sequence")
                                st.write("🔧 Trying alternative fetch method...")
                                
                                # Try alternative fetch method
                                try:
                                    import requests
                                    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id={sequence_id}&rettype=fasta&retmode=text"
                                    response = requests.get(url)
                                    if response.status_code == 200:
                                        fasta_text = response.text
                                        lines = fasta_text.split('\n')
                                        sequence = ''.join(line for line in lines[1:] if not line.startswith('>'))
                                        st.write(f"✅ Alternative method succeeded: {len(sequence)} bp")
                                    else:
                                        st.write(f"❌ Alternative method failed: HTTP {response.status_code}")
                                except Exception as alt_e:
                                    st.write(f"❌ Alternative method error: {alt_e}")
                            
                            if sequence:
                                # Clean the sequence (remove non-DNA characters)
                                clean_sequence = re.sub(r'[^ATGCatgc]', '', sequence.upper())
                                
                                if len(clean_sequence) < 50:
                                    st.error("Sequence is too short for primer design (minimum 50 bp required)")
                                    return
                                
                                # For very large sequences, limit the length
                                if len(clean_sequence) > 1000000:  # 1MB limit
                                    st.warning(f"Sequence is very large ({len(clean_sequence):,} bp). Using first 1MB for primer design.")
                                    clean_sequence = clean_sequence[:1000000]
                                
                                # If we couldn't get detailed info, create basic info
                                if not seq_info:
                                    seq_info = {
                                        "id": sequence_id,
                                        "description": f"Sequence {sequence_id}",
                                        "length": len(clean_sequence),
                                        "organism": organism_name
                                    }
                                
                                st.session_state.sequence_info = seq_info
                                st.session_state.current_sequence = clean_sequence
                                
                                # Note: Species-specific analysis is now handled in comprehensive analysis
                                # Standard primer design will proceed unless comprehensive analysis is enabled
                                target_regions = None
                                
                                st.write(f"🧬 Designing primers for {len(clean_sequence)} bp sequence...")
                                
                                # Add progress bar for primer design
                                design_progress_bar = st.progress(0)
                                design_status_text = st.empty()
                                
                                # Design primers
                                try:
                                    primers = []
                                    t7_results = None
                                    
                                    if enable_t7_design:
                                        # T7 expression primer design
                                        design_status_text.text("🧬 Designing T7 expression primers...")
                                        design_progress_bar.progress(0.1)
                                        
                                        t7_designer = T7PrimerDesigner()
                                        
                                        # Find ORFs first
                                        design_status_text.text("🔍 Finding open reading frames...")
                                        design_progress_bar.progress(0.2)
                                        orfs = t7_designer.find_orfs(clean_sequence, min_length=min_orf_length)
                                        
                                        if orfs:
                                            st.write(f"Found {len(orfs)} potential ORFs")
                                            
                                            # Show ORF selection if not auto-selecting
                                            target_orf = None
                                            if auto_select_orf:
                                                target_orf = orfs[0]  # Longest ORF
                                                st.info(f"Auto-selected ORF: {target_orf['start']}-{target_orf['end']} ({target_orf['length']} bp)")
                                            else:
                                                orf_options = [f"ORF {i+1}: {orf['start']}-{orf['end']} ({orf['length']} bp)" 
                                                              for i, orf in enumerate(orfs[:10])]  # Show top 10
                                                selected_orf_idx = st.selectbox("Select ORF for expression:", range(len(orf_options)), 
                                                                               format_func=lambda x: orf_options[x])
                                                target_orf = orfs[selected_orf_idx]
                                            
                                            # Prepare restriction sites
                                            restriction_sites = None
                                            if add_restriction_sites:
                                                restriction_sites = {
                                                    'forward': forward_site,
                                                    'reverse': reverse_site
                                                }
                                            
                                            # Design T7 primers
                                            design_status_text.text("🧬 Designing T7 expression primers...")
                                            design_progress_bar.progress(0.4)
                                            t7_results, error = t7_designer.design_t7_expression_primers(
                                                clean_sequence,
                                                target_orf=target_orf,
                                                kozak_type=kozak_type,
                                                add_his_tag=add_his_tag,
                                                restriction_sites=restriction_sites
                                            )
                                            
                                            if t7_results:
                                                st.success("T7 expression primers designed!")
                                                
                                                # Validate the construct
                                                design_status_text.text("🔍 Validating expression construct...")
                                                design_progress_bar.progress(0.6)
                                                validation = t7_designer.validate_expression_construct(
                                                    t7_results['forward_primer'],
                                                    t7_results['reverse_primer'],
                                                    clean_sequence
                                                )
                                                
                                                if validation['valid']:
                                                    st.success("Expression construct validation passed!")
                                                else:
                                                    st.warning("Validation issues found:")
                                                    for issue in validation['issues']:
                                                        st.write(f"- {issue}")
                                                
                                                # Store T7 results
                                                st.session_state.t7_results = t7_results
                                                
                                                # Convert to PrimerPair format for compatibility
                                                from Bio.SeqUtils.MeltingTemp import Tm_NN
                                                
                                                forward_tm = Tm_NN(t7_results['forward_primer'][-20:])  # Use last 20bp for Tm calc
                                                reverse_tm = Tm_NN(t7_results['reverse_primer'][-20:])
                                                
                                                t7_primer = PrimerPair(
                                                    forward_seq=t7_results['forward_primer'],
                                                    reverse_seq=t7_results['reverse_primer'],
                                                    forward_tm=forward_tm,
                                                    reverse_tm=reverse_tm,
                                                    product_size=target_orf['length'] + len(t7_results['forward_primer']) + len(t7_results['reverse_primer']),
                                                    gc_content_f=designer.calculate_gc_content(t7_results['forward_primer']),
                                                    gc_content_r=designer.calculate_gc_content(t7_results['reverse_primer']),
                                                    forward_start=target_orf['start'],
                                                    reverse_start=target_orf['end'],
                                                    penalty=0.0
                                                )
                                                
                                                primers.append(t7_primer)
                                                
                                            else:
                                                st.error(f"T7 primer design failed: {error}")
                                        else:
                                            st.warning("No suitable ORFs found for T7 expression")
                                    
                                    if target_regions:
                                        # Design primers for each target region
                                        for i, (start, end) in enumerate(target_regions):
                                            st.write(f"Designing primers for region {i+1}: {start}-{end}")
                                            region_primers = designer.design_primers(
                                                clean_sequence, 
                                                target_region=(start, end),
                                                custom_params=custom_params
                                            )
                                            if region_primers:
                                                primers.extend(region_primers[:3])  # Take top 3 from each region
                                    elif not enable_t7_design:  # Only do standard design if T7 is not enabled
                                        # Standard primer design
                                        design_status_text.text("🧬 Designing standard primers...")
                                        design_progress_bar.progress(0.8)
                                        primers = designer.design_primers(clean_sequence, custom_params=custom_params)
                                    
                                    design_status_text.text("✅ Primer design completed!")
                                    design_progress_bar.progress(1.0)
                                    st.write(f"🔬 Primer design completed. Found {len(primers) if primers else 0} primer pairs.")
                                    
                                    # SET ALL REQUIRED SESSION STATE
                                    st.session_state.primers_designed = primers
                                    st.session_state.current_sequence = clean_sequence  # ADD THIS LINE
                                    
                                    # DEBUG OUTPUT
                                    st.success(f"✅ **Session state set successfully:**")
                                    st.write(f"- Primers designed: {len(primers) if primers else 0}")
                                    st.write(f"- Sequence length: {len(st.session_state.current_sequence)}")
                                    st.write(f"- Sequence info: {st.session_state.sequence_info}")
                                    
                                    if primers:
                                        st.success(f"✅ Successfully designed {len(primers)} primer pairs!")
                                        
                                        # Show all primers
                                        st.subheader("🎯 All Primer Pairs")
                                        preview_data = []
                                        for i, primer in enumerate(primers):  # Show all pairs
                                            preview_data.append({
                                                'Pair': i + 1,
                                                'Forward': primer.forward_seq,
                                                'Reverse': primer.reverse_seq,
                                                'Tm': f"{primer.forward_tm:.1f}°C / {primer.reverse_tm:.1f}°C",
                                                'Size': f"{primer.product_size} bp",
                                                'GC%': f"{primer.gc_content_f:.1f}% / {primer.gc_content_r:.1f}%",
                                                'Penalty': f"{primer.penalty:.3f}"
                                            })
                                        
                                        preview_df = pd.DataFrame(preview_data)
                                        st.dataframe(preview_df, use_container_width=True)
                                        
                                        st.info("📊 Go to the 'Results' tab for detailed primer analysis and export options!")
                                        
                                        # Clear session state after successful search
                                        if 'organism_name' in st.session_state:
                                            del st.session_state.organism_name
                                    else:
                                        st.warning("⚠️ No suitable primers found with current parameters. Try adjusting the primer parameters in the sidebar.")
                                except Exception as primer_e:
                                    st.error(f"❌ Error during primer design: {primer_e}")
                                    st.write("🔧 This might be due to sequence quality or primer parameters.")
                            else:
                                st.error("Failed to fetch sequence. The sequence might be too large or unavailable.")
                                
                        except Exception as e:
                            st.error(f"Error fetching sequence: {e}")
                
                # Add a button to clear search results
                if st.button("🔄 New Search", type="secondary"):
                    st.session_state.search_results = None
                    st.session_state.database_used = None
                    st.rerun()
                
                return  # Skip the search form if we have stored results
            
            col1, col2 = st.columns([2, 1])
            with col1:
                # Initialize organism name from session state if available
                if 'organism_name' in st.session_state:
                    default_organism = st.session_state.organism_name
                else:
                    default_organism = ''
                
                organism_name = st.text_input("Enter organism name:", 
                                            value=default_organism,
                                            key="organism_input",
                                            placeholder="e.g., Fusarium oxysporum, Tetranychus urticae, Beet curly top virus")
                
                # Update session state with current input value
                if organism_name:
                    st.session_state.organism_name = organism_name
                
                # Show organism suggestions
                suggestions = get_organism_suggestions()
                st.write("**Agricultural Pests & Pathogens:**")
                
                for category, subcategories in suggestions.items():
                    st.write(f"**{category}**")
                    
                    for subcategory, organisms in subcategories.items():
                        st.write(f"*{subcategory}*")
                        
                        # Create columns for organisms in this subcategory
                        cols = st.columns(min(len(organisms), 4))
                        for i, organism in enumerate(organisms):
                            with cols[i % 4]:
                                if st.button(organism, key=f"suggest_{category}_{subcategory}_{i}", help=f"Click to search for {organism}"):
                                    st.session_state.organism_name = organism
                                    st.rerun()
                        st.write("")  # Add spacing between subcategories
            with col2:
                max_genomes = st.number_input("Max genomes to search", min_value=1, max_value=20, value=5)
            
            # Additional search options
            with st.expander("Advanced Search Options"):
                col1, col2 = st.columns(2)
                with col1:
                    genome_type = st.selectbox("Genome type:", 
                                             ["Any", "Complete genome", "Chromosome", "Scaffold", "Contig"])
                with col2:
                    assembly_level = st.selectbox("Assembly level:", 
                                                ["Any", "Complete", "Chromosome", "Scaffold", "Contig"])
                
                # Optional gene search
                gene_name = st.text_input("Optional: Search for specific gene (e.g., BRCA1, COI, 16S):", 
                                        placeholder="Leave empty to search entire genome")
            
            if enable_comprehensive:
                if st.button("Comprehensive Genus Analysis", type="primary"):
                    if not email:
                        st.error("Email required for NCBI access")
                    elif not organism_name:
                        st.error("Please enter a genus name (e.g., 'Botrytis', 'Fusarium')")
                    else:
                        with st.spinner("Performing comprehensive taxonomic analysis..."):
                            try:
                                # Extract genus name
                                genus = organism_name.split()[0]
                                
                                # Initialize comprehensive analyzer
                                analyzer = ComprehensiveTaxonomicAnalyzer(NCBIConnector(email, api_key))
                                analyzer.conservation_threshold = conservation_threshold
                                analyzer.specificity_threshold = specificity_threshold
                                
                                # Step 1: Discover subspecies
                                subspecies_data = analyzer.discover_subspecies(genus, max_subspecies)
                                
                                if not subspecies_data:
                                    st.error("No subspecies found for analysis")
                                    return
                                
                                # Display discovered subspecies
                                subspecies_df = pd.DataFrame([
                                    {
                                        'Subspecies': subspecies,
                                        'Available Sequences': data['sequence_count'],
                                        'Organism': data['sample_info'].get('organism', 'Unknown')
                                    }
                                    for subspecies, data in subspecies_data.items()
                                ])
                                
                                st.dataframe(subspecies_df, use_container_width=True)
                                st.success(f"✅ Discovered {len(subspecies_data)} subspecies for analysis")
                                
                                # Step 2: Fetch representative sequences
                                st.info("🔄 **Step 2: Fetching representative sequences from each subspecies**")
                                progress_bar = st.progress(0)
                                status_text = st.empty()
                                
                                subspecies_sequences = analyzer.fetch_representative_sequences(
                                    subspecies_data, 
                                    sequences_per_subspecies,
                                    progress_callback=lambda current, total: (
                                        progress_bar.progress(current / total),
                                        status_text.text(f"Fetching sequences: {current}/{total} subspecies")
                                    )
                                )
                                
                                if not subspecies_sequences:
                                    st.error("No sequences retrieved for analysis")
                                    return
                                
                                # Display sequence collection summary
                                collection_df = pd.DataFrame([
                                    {
                                        'Subspecies': subspecies,
                                        'Sequences Retrieved': data['count'],
                                        'Average Length': int(np.mean([s['length'] for s in data['sequences']]))
                                    }
                                    for subspecies, data in subspecies_sequences.items()
                                ])
                                
                                st.dataframe(collection_df, use_container_width=True)
                                
                                total_sequences = sum(data['count'] for data in subspecies_sequences.values())
                                st.success(f"✅ Collected {total_sequences} sequences for conservation analysis")
                                
                                # Step 3: Multiple sequence alignment and conservation analysis
                                st.info("🧮 **Step 3: Analyzing sequence conservation across subspecies**")
                                progress_bar = st.progress(0)
                                status_text = st.empty()
                                
                                conserved_regions = analyzer.perform_multiple_sequence_alignment(
                                    subspecies_sequences,
                                    progress_callback=lambda current, total: (
                                        progress_bar.progress(current / total),
                                        status_text.text(f"Analyzing conservation: {current}/{total} windows")
                                    )
                                )
                                
                                if not conserved_regions:
                                    st.warning("No conserved regions found with current thresholds")
                                    st.info("Try lowering the conservation threshold or increasing the number of subspecies")
                                    return
                                
                                # Display conserved regions
                                conservation_df = pd.DataFrame([
                                    {
                                        'Region': i + 1,
                                        'Start': region['start'],
                                        'End': region['end'],
                                        'Length': region['end'] - region['start'],
                                        'Conservation': f"{region['conservation_score']:.1%}",
                                        'Subspecies Coverage': f"{region['subspecies_coverage']:.1%}",
                                        'GC Content': f"{region['gc_content']:.1f}%",
                                        'Complexity': f"{region['complexity_score']:.3f}"
                                    }
                                    for i, region in enumerate(conserved_regions)
                                ])
                                
                                st.dataframe(conservation_df, use_container_width=True)
                                st.success(f"✅ Found {len(conserved_regions)} highly conserved regions")
                                
                                # Step 4: Specificity testing against other taxa
                                st.info("🎯 **Step 4: Testing specificity against related taxa**")
                                progress_bar = st.progress(0)
                                status_text = st.empty()
                                
                                comparison_genera = None
                                if custom_comparison_genera.strip():
                                    comparison_genera = [g.strip() for g in custom_comparison_genera.split(',')]

                                comparison_species = None  
                                if custom_comparison_species.strip():
                                    comparison_species = [s.strip() for s in custom_comparison_species.split(',')]

                                try:
                                    specificity_results = analyzer.compare_against_other_taxa(
                                        conserved_regions, 
                                        genus,
                                        test_genera=test_genus_specificity,
                                        test_species=test_species_specificity,
                                        comparison_genera=comparison_genera,
                                        comparison_species=comparison_species,
                                        progress_callback=lambda current, total: (
                                            progress_bar.progress(current / total),
                                            status_text.text(f"Testing specificity: {current}/{total} regions")
                                        )
                                    )
                                except Exception as e:
                                    st.error(f"❌ Specificity testing failed: {e}")
                                    st.write("This might be due to NCBI API issues or network problems.")
                                    return

                                specific_regions = specificity_results['combined_specific_regions']
                                
                                # Debug information
                                st.write(f"🔍 **Debug Info:**")
                                st.write(f"- Conserved regions found: {len(conserved_regions)}")
                                st.write(f"- Specific regions found: {len(specific_regions)}")
                                st.write(f"- Specificity results keys: {list(specificity_results.keys())}")

                                if not specific_regions:
                                    st.warning("No regions found that pass specificity requirements")
                                    st.info("Try lowering the specificity threshold or adjusting comparison organisms")
                                    return

                                # Display detailed specificity results
                                if test_genus_specificity and specificity_results['genus_results']:
                                    genus_results = specificity_results['genus_results']
                                    st.write(f"**Genus-level tested against:** {', '.join(genus_results['tested_against'])}")
                                    st.write(f"**Regions passing genus test:** {len(genus_results['specific_regions'])}")

                                if test_species_specificity and specificity_results['species_results']:
                                    species_results = specificity_results['species_results']
                                    st.write(f"**Species-level tested against:** {', '.join(species_results['tested_against'])}")
                                    st.write(f"**Regions passing species test:** {len(species_results['specific_regions'])}")

                                # Display final specific regions with enhanced information
                                specificity_df_data = []
                                for i, region in enumerate(specific_regions):
                                    row = {
                                        'Region': i + 1,
                                        'Position': f"{region['start']}-{region['end']}",
                                        'Length': region['end'] - region['start'],
                                        'Conservation': f"{region['conservation_score']:.1%}",
                                        'Subspecies Coverage': f"{region['subspecies_coverage']:.1%}",
                                    }
                                    
                                    # Add specificity columns based on what was tested
                                    if 'genus_specificity' in region:
                                        row['Genus Specificity'] = f"{region['genus_specificity']:.1%}"
                                    if 'species_specificity' in region:
                                        row['Species Specificity'] = f"{region['species_specificity']:.1%}"
                                    if 'combined_specificity' in region:
                                        row['Combined Specificity'] = f"{region['combined_specificity']:.1%}"
                                    else:
                                        row['Overall Specificity'] = f"{region['specificity_score']:.1%}"
                                    
                                    row['Quality Score'] = f"{region['conservation_score'] * region['specificity_score']:.3f}"
                                    specificity_df_data.append(row)

                                specificity_df = pd.DataFrame(specificity_df_data)
                                st.dataframe(specificity_df, use_container_width=True)
                                st.success(f"✅ Final high-quality regions: {len(specific_regions)}")
                                
                                # Step 5: Design primers for best regions
                                st.info("🧬 **Step 5: Designing genus-specific primers**")
                                
                                # Select top regions for primer design
                                top_regions = sorted(
                                    specific_regions, 
                                    key=lambda x: x['conservation_score'] * x['specificity_score'], 
                                    reverse=True
                                )[:5]  # Top 5 regions
                                
                                all_primers = []
                                designer = PrimerDesigner()
                                
                                # Use the first sequence as template (they should be similar in conserved regions)
                                template_sequence = list(subspecies_sequences.values())[0]['sequences'][0]['sequence']
                                
                                st.write(f"Designing primers for top {len(top_regions)} regions...")
                                
                                # Add progress bar for primer design
                                primer_progress_bar = st.progress(0)
                                primer_status_text = st.empty()
                                
                                for i, region in enumerate(top_regions):
                                    primer_progress_bar.progress((i + 1) / len(top_regions))
                                    primer_status_text.text(f"Designing primers for region {i+1}/{len(top_regions)}: {region['start']}-{region['end']}")
                                    
                                    try:
                                        region_primers = designer.design_primers(
                                            template_sequence,
                                            target_region=(region['start'], region['end']),
                                            custom_params=custom_params
                                        )
                                        
                                        if region_primers:
                                            # Add metadata to primers
                                            for primer in region_primers:
                                                primer.region_info = {
                                                    'conservation_score': region['conservation_score'],
                                                    'specificity_score': region['specificity_score'],
                                                    'quality_score': region['conservation_score'] * region['specificity_score'],
                                                    'region_number': i + 1
                                                }
                                            
                                            all_primers.extend(region_primers[:2])  # Top 2 from each region
                                            st.write(f"  ✅ Designed {len(region_primers[:2])} primer pairs")
                                        else:
                                            st.write(f"  ⚠️ No primers found for this region")
                                    
                                    except Exception as e:
                                        st.warning(f"Primer design failed for region {i+1}: {e}")

                                if all_primers:
                                    # Store results in session state for other tabs
                                    st.session_state.primers_designed = all_primers
                                    st.session_state.current_sequence = template_sequence
                                    st.session_state.sequence_info = {
                                        "description": f"Comprehensive analysis: {genus} genus",
                                        "length": len(template_sequence),
                                        "organism": genus,
                                        "id": f"{genus}_comprehensive_analysis"
                                    }
                                    
                                    # Store comprehensive analysis results
                                    st.session_state.comprehensive_analysis_results = {
                                        'subspecies_data': subspecies_data,
                                        'conserved_regions': conserved_regions,
                                        'specific_regions': specific_regions,
                                        'specificity_results': specificity_results,
                                        'analysis_summary': {
                                            'total_subspecies': len(subspecies_data),
                                            'total_sequences': total_sequences,
                                            'conserved_regions_found': len(conserved_regions),
                                            'specific_regions_found': len(specific_regions),
                                            'genus_specificity_tested': test_genus_specificity,
                                            'species_specificity_tested': test_species_specificity,
                                            'primers_designed': len(all_primers)
                                        }
                                    }
                                    
                                    # Debug: Verify session state storage
                                    st.write(f"🔍 **Session State Debug:**")
                                    st.write(f"- Comprehensive analysis results stored: {hasattr(st.session_state, 'comprehensive_analysis_results')}")
                                    st.write(f"- Primers designed: {len(st.session_state.primers_designed)}")
                                    st.write(f"- Analysis summary: {st.session_state.comprehensive_analysis_results['analysis_summary']}")
                                    
                                    st.success(f"🎉 **Comprehensive analysis complete!** Designed {len(all_primers)} high-quality genus-specific primers for {genus}")
                                    st.info("📊 Go to the 'Results', 'Conservation', and 'Specificity' tabs to view detailed analysis and primer information!")
                                    
                                    # Show quick preview of results
                                    preview_data = []
                                    for i, primer in enumerate(all_primers):
                                        preview_data.append({
                                            'Pair': i + 1,
                                            'Region': primer.region_info['region_number'],
                                            'Forward': primer.forward_seq,
                                            'Reverse': primer.reverse_seq,
                                            'Quality': f"{primer.region_info['quality_score']:.3f}"
                                        })
                                    
                                    preview_df = pd.DataFrame(preview_data)
                                    st.dataframe(preview_df, use_container_width=True)
                                    
                                else:
                                    st.error("No primers could be designed from the identified regions")
                                    st.info("Try adjusting primer parameters or specificity thresholds")
                                
                            except Exception as e:
                                st.error(f"Comprehensive analysis failed: {e}")
                                import traceback
                                st.code(traceback.format_exc())
            else:
                if st.button("Search Organism Genomes", type="primary"):
                    if not email:
                        st.error("❌ **Email Required**: Please enter an email address in the sidebar to access NCBI databases. You can use any valid email address or click 'Use test email' for demo purposes.")
                    elif not organism_name:
                        st.error("❌ **Organism Name Required**: Please enter an organism name or click one of the suggested organisms above.")
                    else:
                        # Clear previous search results when starting a new search
                        st.session_state.search_results = None
                        st.session_state.database_used = None
                        with st.spinner(f"Searching for {organism_name} genomes..."):
                            try:
                                ncbi = NCBIConnector(email, api_key)
                                designer = PrimerDesigner()
                                
                                # Build search query
                                query_parts = [f'"{organism_name}"[organism]']
                                
                                if gene_name:
                                    query_parts.append(f'"{gene_name}"[gene]')
                                
                                if genome_type != "Any":
                                    query_parts.append(f'"{genome_type}"[title]')
                                
                                if assembly_level != "Any":
                                    query_parts.append(f'"{assembly_level}"[assembly_level]')
                                
                                search_query = " AND ".join(query_parts)
                                
                                # Search for genome sequences
                                st.write(f"Searching with query: `{search_query}`")
                                
                                # Search in genome database first
                                genome_ids = ncbi.search_sequences(search_query, database="genome", max_results=max_genomes)
                                
                                # If no genomes found, try nucleotide database
                                if not genome_ids:
                                    st.info("No genomes found, searching nucleotide database...")
                                    genome_ids = ncbi.search_sequences(search_query, database="nucleotide", max_results=max_genomes)
                                    database_used = "nucleotide"
                                else:
                                    database_used = "genome"
                                
                                if genome_ids:
                                    if database_used == "genome":
                                        st.success(f"Found {len(genome_ids)} genome assemblies!")
                                    else:
                                        st.success(f"Found {len(genome_ids)} nucleotide sequences!")
                                    
                                    # Store search results in session state
                                    st.session_state.database_used = database_used
                                    
                                    # Display found sequences
                                    if database_used == "genome":
                                        st.subheader("Available Genome Assemblies")
                                    else:
                                        st.subheader("Available Nucleotide Sequences")
                                    genome_info = []
                                    
                                    for i, genome_id in enumerate(genome_ids):
                                        with st.spinner(f"Fetching info {i+1}/{len(genome_ids)}..."):
                                            try:
                                                if database_used == "genome":
                                                    # For genome database, get assembly info
                                                    info = ncbi.fetch_genome_assembly_info(genome_id)
                                                else:
                                                    # For nucleotide database, get sequence info
                                                    info = ncbi.fetch_sequence_info(genome_id, database="nucleotide")
                                                
                                                if info:
                                                    genome_info.append({
                                                        'ID': genome_id,
                                                        'Description': info.get('description', 'N/A'),
                                                        'Length': info.get('length', 'N/A'),
                                                        'Organism': info.get('organism', 'N/A')
                                                    })
                                            except Exception as e:
                                                st.warning(f"Could not fetch info for {genome_id}: {e}")
                                                # Add basic info even if detailed fetch fails
                                                genome_info.append({
                                                    'ID': genome_id,
                                                    'Description': f'Sequence {genome_id}' if database_used == "nucleotide" else f'Genome Assembly {genome_id}',
                                                    'Length': 'N/A',
                                                    'Organism': organism_name
                                                })
                                    
                                    if genome_info:
                                        # Store search results in session state
                                        st.session_state.search_results = genome_info
                                        
                                        genome_df = pd.DataFrame(genome_info)
                                        st.dataframe(genome_df, use_container_width=True)
                                        
                                        # Let user select a sequence
                                        if database_used == "genome":
                                            select_text = "Select a genome assembly to design primers for:"
                                            button_text = "Design Primers for Selected Assembly"
                                        else:
                                            select_text = "Select a nucleotide sequence to design primers for:"
                                            button_text = "Design Primers for Selected Sequence"
                                        
                                        selected_genome = st.selectbox(
                                            select_text,
                                            range(len(genome_info)),
                                            format_func=lambda x: f"{genome_info[x]['ID']} - {genome_info[x]['Description'][:100]}..."
                                        )
                                        
                                        if st.button(button_text, type="primary"):
                                            with st.spinner("Fetching sequence and designing primers..."):
                                                sequence_id = genome_info[selected_genome]['ID']
                                                
                                                try:
                                                    if database_used == "genome":
                                                        # For genome assemblies, search for nucleotide sequences
                                                        assembly_query = f'"{sequence_id}"[Assembly]'
                                                        nucleotide_ids = ncbi.search_sequences(assembly_query, database="nucleotide", max_results=10)
                                                        
                                                        if nucleotide_ids:
                                                            st.info(f"Found {len(nucleotide_ids)} sequences from this assembly. Using the first one for primer design.")
                                                            sequence_id = nucleotide_ids[0]
                                                        else:
                                                            st.warning("No nucleotide sequences found for this assembly.")
                                                            return
                                                    
                                                    # Fetch the sequence
                                                    st.write(f"🔍 Fetching sequence {sequence_id} from {database_used} database...")
                                                    
                                                    if database_used == "nucleotide":
                                                        # For nucleotide sequences, fetch directly
                                                        st.write("📡 Using nucleotide database for direct fetch...")
                                                        sequence = ncbi.fetch_sequence(sequence_id, database="nucleotide")
                                                        seq_info = ncbi.fetch_sequence_info(sequence_id, database="nucleotide")
                                                    else:
                                                        # For genome assemblies, use the nucleotide ID we found
                                                        st.write("📡 Using default database for genome assembly...")
                                                        sequence = ncbi.fetch_sequence(sequence_id)
                                                        seq_info = ncbi.fetch_sequence_info(sequence_id)
                                                    
                                                    if sequence:
                                                        st.write(f"✅ Successfully fetched sequence: {len(sequence)} bp")
                                                        st.write(f"📝 First 100 characters: {sequence[:100]}...")
                                                    else:
                                                        st.write("❌ Failed to fetch sequence")
                                                        st.write("🔧 Trying alternative fetch method...")
                                                    
                                                        # Try alternative fetch method
                                                        try:
                                                            import requests
                                                            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id={sequence_id}&rettype=fasta&retmode=text"
                                                            response = requests.get(url)
                                                            if response.status_code == 200:
                                                                fasta_text = response.text
                                                                lines = fasta_text.split('\n')
                                                                sequence = ''.join(line for line in lines[1:] if not line.startswith('>'))
                                                                st.write(f"✅ Alternative method succeeded: {len(sequence)} bp")
                                                            else:
                                                                st.write(f"❌ Alternative method failed: HTTP {response.status_code}")
                                                        except Exception as alt_e:
                                                            st.write(f"❌ Alternative method error: {alt_e}")
                                                    
                                                    if sequence:
                                                        # Clean the sequence (remove non-DNA characters)
                                                        clean_sequence = re.sub(r'[^ATGCatgc]', '', sequence.upper())
                                                        
                                                        if len(clean_sequence) < 50:
                                                            st.error("Sequence is too short for primer design (minimum 50 bp required)")
                                                            return
                                                        
                                                        # For very large sequences, limit the length
                                                        if len(clean_sequence) > 1000000:  # 1MB limit
                                                            st.warning(f"Sequence is very large ({len(clean_sequence):,} bp). Using first 1MB for primer design.")
                                                            clean_sequence = clean_sequence[:1000000]
                                                        
                                                        # If we couldn't get detailed info, create basic info
                                                        if not seq_info:
                                                            seq_info = {
                                                                "id": sequence_id,
                                                                "description": f"Sequence {sequence_id}",
                                                                "length": len(clean_sequence),
                                                                "organism": organism_name
                                                            }
                                                        
                                                        st.session_state.sequence_info = seq_info
                                                        st.session_state.current_sequence = clean_sequence
                                                        
                                                        st.write(f"🧬 Designing primers for {len(clean_sequence)} bp sequence...")
                                                        st.write(f"⚙️ Using parameters: min_size={custom_params.get('PRIMER_MIN_SIZE')}, max_size={custom_params.get('PRIMER_MAX_SIZE')}, min_tm={custom_params.get('PRIMER_MIN_TM')}")
                                                        
                                                        # Design primers
                                                        try:
                                                            primers = designer.design_primers(clean_sequence, custom_params=custom_params)
                                                            st.write(f"🔬 Primer design completed. Found {len(primers) if primers else 0} primer pairs.")
                                                            
                                                            st.session_state.primers_designed = primers
                                                            
                                                            if primers:
                                                                st.success(f"✅ Successfully designed {len(primers)} primer pairs!")
                                                                # Clear session state after successful search
                                                                if 'organism_name' in st.session_state:
                                                                    del st.session_state.organism_name
                                                            else:
                                                                st.warning("⚠️ No suitable primers found with current parameters. Try adjusting the primer parameters in the sidebar.")
                                                        except Exception as primer_e:
                                                            st.error(f"❌ Error during primer design: {primer_e}")
                                                            st.write("🔧 This might be due to sequence quality or primer parameters.")
                                                    else:
                                                        st.error("Failed to fetch sequence. The sequence might be too large or unavailable.")
                                                        
                                                except Exception as e:
                                                    st.error(f"Error fetching sequence: {e}")
                                        else:
                                            st.warning("No genome information could be retrieved")
                                    else:
                                        st.warning(f"No genomes found for organism: {organism_name}")
                                        st.info("Try adjusting your search terms or check the spelling of the organism name.")
                                        
                            except Exception as e:
                                st.error(f"Error searching for organism: {e}")
        
        elif input_method == "GenBank ID":
            genbank_id = st.text_input("Enter GenBank Accession Number:", 
                                     placeholder="e.g., NM_000314.6")
            
            if st.button("Fetch and Design Primers", type="primary"):
                if not email:
                    st.error("❌ **Email Required**: Please enter an email address in the sidebar to access NCBI databases. You can use any valid email address or click 'Use test email' for demo purposes.")
                elif not genbank_id:
                    st.error("Please provide a GenBank ID")
                else:
                    with st.spinner("Fetching sequence and designing primers..."):
                        try:
                            ncbi = NCBIConnector(email, api_key)
                            designer = PrimerDesigner()
                            
                            # Fetch sequence
                            sequence = ncbi.fetch_sequence(genbank_id)
                            if sequence:
                                seq_info = ncbi.fetch_sequence_info(genbank_id)
                                st.session_state.sequence_info = seq_info
                                st.session_state.current_sequence = sequence
                                
                                # Design primers
                                primers = designer.design_primers(sequence, custom_params=custom_params)
                                st.session_state.primers_designed = primers
                                
                                # DEBUG OUTPUT
                                st.success(f"✅ **Session state set successfully:**")
                                st.write(f"- Primers designed: {len(primers) if primers else 0}")
                                st.write(f"- Sequence length: {len(st.session_state.current_sequence)}")
                                st.write(f"- Sequence info: {st.session_state.sequence_info}")
                                
                                if primers:
                                    st.success(f"Successfully designed {len(primers)} primer pairs!")
                                    st.info("📊 Go to the 'Results', 'Analysis', and 'Export' tabs to view detailed primer information!")
                                else:
                                    st.warning("No suitable primers found with current parameters")
                            else:
                                st.error("Failed to fetch sequence")
                        except Exception as e:
                            st.error(f"Error: {e}")
        
        elif input_method == "NCBI Search":
            search_query = st.text_input("Enter NCBI search query:", 
                                       placeholder="e.g., Homo sapiens[organism] AND BRCA1[gene]")
            max_results = st.slider("Maximum sequences to process", 1, 10, 3)
            
            if st.button("Search and Design Primers", type="primary"):
                if not email:
                    st.error("❌ **Email Required**: Please enter an email address in the sidebar to access NCBI databases. You can use any valid email address or click 'Use test email' for demo purposes.")
                elif not search_query:
                    st.error("Please provide a search query")
                else:
                    with st.spinner("Searching NCBI and designing primers..."):
                        try:
                            ncbi = NCBIConnector(email, api_key)
                            designer = PrimerDesigner()
                            
                            # Search sequences
                            seq_ids = ncbi.search_sequences(search_query, max_results=max_results)
                            
                            if seq_ids:
                                all_primers = []
                                all_info = []
                                
                                for seq_id in seq_ids:
                                    seq_info = ncbi.fetch_sequence_info(seq_id)
                                    if seq_info and seq_info.get('sequence'):
                                        primers = designer.design_primers(
                                            seq_info['sequence'], 
                                            custom_params=custom_params
                                        )
                                        if primers:
                                            all_primers.extend(primers)
                                            all_info.append(seq_info)
                                
                                st.session_state.primers_designed = all_primers
                                st.session_state.sequence_info = {"multiple": all_info}
                                # Note: current_sequence not set for multiple sequences
                                
                                # DEBUG OUTPUT
                                st.success(f"✅ **Session state set successfully:**")
                                st.write(f"- Primers designed: {len(all_primers) if all_primers else 0}")
                                st.write(f"- Sequences processed: {len(all_info)}")
                                st.write(f"- Sequence info: {st.session_state.sequence_info}")
                                
                                if all_primers:
                                    st.success(f"Found {len(all_primers)} primer pairs from {len(all_info)} sequences!")
                                    st.info("📊 Go to the 'Results', 'Analysis', and 'Export' tabs to view detailed primer information!")
                                else:
                                    st.warning("No suitable primers found")
                            else:
                                st.warning("No sequences found for the search query")
                        except Exception as e:
                            st.error(f"Error: {e}")
        
        elif input_method == "Direct Sequence":
            sequence_input = st.text_area("Enter DNA sequence:", 
                                        placeholder="ATGCGATCGATCG...",
                                        height=150)
            
            # Optional target region
            st.subheader("Target Region (Optional)")
            col1, col2 = st.columns(2)
            with col1:
                target_start = st.number_input("Start position", min_value=0, value=0)
            with col2:
                target_end = st.number_input("End position", min_value=0, value=0)
            
            if st.button("Design Primers", type="primary"):
                if not sequence_input:
                    st.error("Please provide a DNA sequence")
                else:
                    with st.spinner("Designing primers..."):
                        try:
                            designer = PrimerDesigner()
                            
                            # Clean sequence
                            clean_seq = re.sub(r'[^ATGCatgc]', '', sequence_input.upper())
                            
                            target_region = None
                            if target_start > 0 and target_end > target_start:
                                target_region = (target_start, target_end)
                            
                            primers = designer.design_primers(
                                clean_seq, 
                                target_region=target_region,
                                custom_params=custom_params
                            )
                            
                            # SET ALL REQUIRED SESSION STATE
                            st.session_state.primers_designed = primers
                            st.session_state.current_sequence = clean_seq
                            st.session_state.sequence_info = {
                                "length": len(clean_seq),
                                "description": "User-provided sequence"
                            }
                            
                            # DEBUG OUTPUT
                            st.success(f"✅ **Session state set successfully:**")
                            st.write(f"- Primers designed: {len(primers) if primers else 0}")
                            st.write(f"- Sequence length: {len(st.session_state.current_sequence)}")
                            st.write(f"- Sequence info: {st.session_state.sequence_info}")
                            st.write(f"- Session state keys: {list(st.session_state.keys())}")
                            
                            if primers:
                                st.success(f"Successfully designed {len(primers)} primer pairs!")
                                st.info("📊 Go to the 'Results', 'Analysis', and 'Export' tabs to view detailed primer information!")
                            else:
                                st.warning("No suitable primers found with current parameters")
                        except Exception as e:
                            st.error(f"Error: {e}")
        
        elif input_method == "Upload File":
            uploaded_file = st.file_uploader("Choose a FASTA file", type=['fasta', 'fa', 'txt'])
            
            if uploaded_file is not None:
                if st.button("Process File", type="primary"):
                    with st.spinner("Processing file..."):
                        try:
                            # Read file content
                            content = uploaded_file.read().decode('utf-8')
                            
                            # Parse FASTA
                            if content.startswith('>'):
                                lines = content.split('\n')
                                sequence = ''.join(line for line in lines[1:] if not line.startswith('>'))
                            else:
                                sequence = content
                            
                            # Clean sequence
                            clean_seq = re.sub(r'[^ATGCatgc]', '', sequence.upper())
                            
                            designer = PrimerDesigner()
                            primers = designer.design_primers(clean_seq, custom_params=custom_params)
                            
                            # SET ALL REQUIRED SESSION STATE
                            st.session_state.primers_designed = primers
                            st.session_state.current_sequence = clean_seq
                            st.session_state.sequence_info = {
                                "length": len(clean_seq),
                                "description": f"Uploaded file: {uploaded_file.name}"
                            }
                            
                            # DEBUG OUTPUT
                            st.success(f"✅ **Session state set successfully:**")
                            st.write(f"- Primers designed: {len(primers) if primers else 0}")
                            st.write(f"- Sequence length: {len(st.session_state.current_sequence)}")
                            st.write(f"- Sequence info: {st.session_state.sequence_info}")
                            
                            if primers:
                                st.success(f"Successfully designed {len(primers)} primer pairs!")
                                st.info("📊 Go to the 'Results', 'Analysis', and 'Export' tabs to view detailed primer information!")
                            else:
                                st.warning("No suitable primers found")
                        except Exception as e:
                            st.error(f"Error processing file: {e}")
    
    with tab2:
        st.header("Primer Design Results")
        
        if st.session_state.primers_designed:
            primers = st.session_state.primers_designed
            
            # Show analysis type
            if hasattr(st.session_state, 'comprehensive_analysis_results'):
                st.info("🧬 These primers were designed using **Comprehensive Taxonomic Analysis** with enhanced specificity testing")
                
                # Show analysis summary
                analysis_summary = st.session_state.comprehensive_analysis_results['analysis_summary']
                col1, col2, col3, col4 = st.columns(4)
                with col1:
                    st.metric("Subspecies Analyzed", analysis_summary['total_subspecies'])
                with col2:
                    st.metric("Total Sequences", analysis_summary['total_sequences'])
                with col3:
                    st.metric("Conserved Regions", analysis_summary['conserved_regions_found'])
                with col4:
                    st.metric("Final Specific Regions", analysis_summary['specific_regions_found'])
            
            # Sequence information
            if st.session_state.sequence_info:
                st.subheader("Sequence Information")
                info = st.session_state.sequence_info
                
                if "multiple" in info:
                    st.write(f"Processed {len(info['multiple'])} sequences")
                else:
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        st.metric("Length", f"{info.get('length', 'N/A')} bp")
                    with col2:
                        st.metric("Organism", info.get('organism', 'N/A'))
                    with col3:
                        st.metric("ID", info.get('id', 'N/A'))
                    
                    if 'description' in info:
                        st.write(f"**Description:** {info['description']}")
            
            # Primer results table
            st.subheader("Primer Pairs")
            
            # Create DataFrame for display
            data = []
            for i, primer in enumerate(primers):
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
                
                # Add quality information if available from comprehensive analysis
                if hasattr(primer, 'region_info'):
                    row['Region'] = primer.region_info['region_number']
                    row['Quality Score'] = f"{primer.region_info['quality_score']:.3f}"
                
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
            
            # T7 Expression Results
            if hasattr(st.session_state, 't7_results') and st.session_state.t7_results:
                try:
                    st.subheader("T7 Expression Design")
                    
                    t7_res = st.session_state.t7_results
                    
                    col1, col2 = st.columns(2)
                    with col1:
                        st.write("**Expression Forward Primer:**")
                        st.code(t7_res['forward_primer'], language="text")
                        st.write(f"Length: {len(t7_res['forward_primer'])} bp")
                        
                        features = t7_res['expression_features']
                        st.write("**Features included:**")
                        if features.get('t7_promoter'):
                            st.write("✅ T7 Promoter")
                        if features.get('kozak_sequence'):
                            st.write(f"✅ Kozak sequence ({features['kozak_sequence']})")
                        if features.get('his_tag'):
                            st.write("✅ His6 tag")
                        if features.get('restriction_sites'):
                            sites = features['restriction_sites']
                            if sites:
                                st.write(f"✅ Restriction sites: {', '.join(sites.values())}")
                    
                    with col2:
                        st.write("**Expression Reverse Primer:**")
                        st.code(t7_res['reverse_primer'], language="text")
                        st.write(f"Length: {len(t7_res['reverse_primer'])} bp")
                        
                        orf = t7_res['target_orf']
                        st.write("**Target ORF:**")
                        st.write(f"Position: {orf['start']}-{orf['end']}")
                        st.write(f"Length: {orf['length']} bp")
                        st.write(f"Reading frame: {orf['frame']}")
                    
                    # T7 construct diagram
                    t7_fig = create_t7_construct_diagram(st.session_state.current_sequence, t7_res)
                    if t7_fig:
                        st.plotly_chart(t7_fig, use_container_width=True)
                        
                except Exception as e:
                    st.warning(f"T7 results display error: {e}")
            
            # Comprehensive Analysis Results
            if hasattr(st.session_state, 'comprehensive_analysis_results') and st.session_state.comprehensive_analysis_results:
                st.subheader("Comprehensive Taxonomic Analysis Results")
                
                analysis_results = st.session_state.comprehensive_analysis_results
                summary = analysis_results['analysis_summary']
                
                # Display analysis summary
                col1, col2, col3, col4 = st.columns(4)
                with col1:
                    st.metric("Subspecies Analyzed", summary['total_subspecies'])
                with col2:
                    st.metric("Total Sequences", summary['total_sequences'])
                with col3:
                    st.metric("Conserved Regions", summary['conserved_regions_found'])
                with col4:
                    st.metric("Specific Regions", summary['specific_regions_found'])
                
                # Display subspecies data
                st.subheader("Subspecies Analysis")
                subspecies_df = pd.DataFrame([
                    {
                        'Subspecies': subspecies,
                        'Available Sequences': data['sequence_count'],
                        'Organism': data['sample_info'].get('organism', 'Unknown')
                    }
                    for subspecies, data in analysis_results['subspecies_data'].items()
                ])
                st.dataframe(subspecies_df, use_container_width=True)
                
                # Display specific regions
                if analysis_results['specific_regions']:
                    st.subheader("Genus-Specific Regions")
                    specific_df = pd.DataFrame([
                        {
                            'Region': i + 1,
                            'Position': f"{region['start']}-{region['end']}",
                            'Length': region['end'] - region['start'],
                            'Conservation': f"{region['conservation_score']:.1%}",
                            'Specificity': f"{region['specificity_score']:.1%}",
                            'Quality Score': region['conservation_score'] * region['specificity_score']
                        }
                        for i, region in enumerate(analysis_results['specific_regions'])
                    ])
                    st.dataframe(specific_df, use_container_width=True)
                
                # Display primer quality information
                if st.session_state.primers_designed:
                    st.subheader("Primer Quality Analysis")
                    primer_quality_data = []
                    for i, primer in enumerate(st.session_state.primers_designed):
                        if hasattr(primer, 'region_info'):
                            primer_quality_data.append({
                                'Primer Pair': i + 1,
                                'Region': primer.region_info['region_number'],
                                'Conservation Score': f"{primer.region_info['conservation_score']:.1%}",
                                'Specificity Score': f"{primer.region_info['specificity_score']:.1%}",
                                'Quality Score': f"{primer.region_info['quality_score']:.3f}",
                                'Forward Tm': f"{primer.forward_tm:.1f}°C",
                                'Reverse Tm': f"{primer.reverse_tm:.1f}°C",
                                'Product Size': f"{primer.product_size} bp"
                            })
                    
                    if primer_quality_data:
                        quality_df = pd.DataFrame(primer_quality_data)
                        st.dataframe(quality_df, use_container_width=True)
                        
                        # Quality score distribution
                        quality_scores = [float(data['Quality Score']) for data in primer_quality_data]
                        avg_quality = np.mean(quality_scores)
                        st.metric("Average Quality Score", f"{avg_quality:.3f}")
                        
                        if avg_quality >= 0.8:
                            st.success("🎯 Excellent primer quality - suitable for diagnostic applications")
                        elif avg_quality >= 0.6:
                            st.info("✅ Good primer quality - suitable for research applications")
                        else:
                            st.warning("⚠️ Moderate primer quality - consider adjusting parameters")
        else:
            st.info("No primers designed yet. Please use the Input tab to design primers.")
    
    with tab3:
        st.header("Primer Analysis")
        
        # Debug session state
        st.write("🔍 Debug - Analysis tab session state:")
        st.write(f"Primers designed: {len(st.session_state.primers_designed) if st.session_state.primers_designed else 0}")
        
        if st.session_state.primers_designed:
            primers = st.session_state.primers_designed
            
            # Create visualizations
            fig = create_primer_visualization(primers)
            if fig:
                st.plotly_chart(fig, use_container_width=True)
            
            # Comprehensive analysis visualization
            if hasattr(st.session_state, 'comprehensive_analysis_results') and st.session_state.comprehensive_analysis_results:
                if 'specific_regions' in st.session_state.comprehensive_analysis_results:
                    st.subheader("Comprehensive Analysis Visualization")
                    
                    species_fig = create_species_specific_visualization(
                        st.session_state.current_sequence,
                        st.session_state.comprehensive_analysis_results.get('specific_regions', []),
                        primers
                    )
                    
                    if species_fig:
                        st.plotly_chart(species_fig, use_container_width=True)
            
            # Sequence diagram
            if st.session_state.current_sequence:
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
            
            # Statistics
            st.subheader("Statistics")
            col1, col2, col3, col4 = st.columns(4)
            
            with col1:
                avg_tm = sum(p.forward_tm + p.reverse_tm for p in primers) / (2 * len(primers))
                st.metric("Average Tm", f"{avg_tm:.1f}°C")
            
            with col2:
                avg_gc = sum(p.gc_content_f + p.gc_content_r for p in primers) / (2 * len(primers))
                st.metric("Average GC%", f"{avg_gc:.1f}%")
            
            with col3:
                avg_product = sum(p.product_size for p in primers) / len(primers)
                st.metric("Average Product Size", f"{avg_product:.0f} bp")
            
            with col4:
                best_penalty = min(p.penalty for p in primers)
                st.metric("Best Penalty Score", f"{best_penalty:.3f}")
        else:
            st.info("No primers designed yet. Please use the Input tab to design primers.")
    
    with tab4:  # Conservation tab
        st.header("Conservation Analysis")
        
        if hasattr(st.session_state, 'comprehensive_analysis_results') and st.session_state.comprehensive_analysis_results:
            analysis_results = st.session_state.comprehensive_analysis_results
            
            # Workflow diagram
            workflow_fig = create_analysis_workflow_diagram()
            if workflow_fig:
                st.plotly_chart(workflow_fig, use_container_width=True)
            
            # Conservation overview
            st.subheader("Conservation Overview")
            conserved_regions = analysis_results.get('conserved_regions', [])
            
            if conserved_regions:
                # Conservation heatmap
                heatmap_fig = create_conservation_heatmap(
                    conserved_regions, 
                    analysis_results['subspecies_data']
                )
                if heatmap_fig:
                    st.plotly_chart(heatmap_fig, use_container_width=True)
                
                # Conservation statistics
                st.subheader("Conservation Statistics")
                col1, col2, col3, col4 = st.columns(4)
                
                conservation_scores = [r['conservation_score'] for r in conserved_regions]
                subspecies_coverage = [r['subspecies_coverage'] for r in conserved_regions]
                
                with col1:
                    st.metric("Average Conservation", f"{np.mean(conservation_scores):.1%}")
                with col2:
                    st.metric("Max Conservation", f"{np.max(conservation_scores):.1%}")
                with col3:
                    st.metric("Average Subspecies Coverage", f"{np.mean(subspecies_coverage):.1%}")
                with col4:
                    st.metric("Regions Found", len(conserved_regions))
                
                # Detailed conservation table
                st.subheader("Detailed Conservation Results")
                conservation_df = pd.DataFrame([
                    {
                        'Region': i + 1,
                        'Position': f"{region['start']}-{region['end']}",
                        'Length': region['end'] - region['start'],
                        'Conservation Score': f"{region['conservation_score']:.1%}",
                        'Subspecies Coverage': f"{region['subspecies_coverage']:.1%}",
                        'GC Content': f"{region['gc_content']:.1f}%",
                        'Complexity Score': f"{region['complexity_score']:.3f}",
                        'Sequence Count': region['sequence_count']
                    }
                    for i, region in enumerate(conserved_regions)
                ])
                
                st.dataframe(conservation_df, use_container_width=True)
        else:
            st.info("No conservation analysis results available. Run comprehensive analysis first.")

    with tab5:  # Specificity tab
        st.header("Specificity Analysis")
        
        if hasattr(st.session_state, 'comprehensive_analysis_results') and st.session_state.comprehensive_analysis_results:
            analysis_results = st.session_state.comprehensive_analysis_results
            specificity_results = analysis_results.get('specificity_results', {})
            
            if specificity_results:
                # Specificity comparison chart
                specificity_fig = create_specificity_comparison_chart(specificity_results)
                if specificity_fig:
                    st.plotly_chart(specificity_fig, use_container_width=True)
                
                # Specificity overview
                st.subheader("Specificity Test Results")
                
                col1, col2 = st.columns(2)
                
                with col1:
                    if 'genus_results' in specificity_results and specificity_results['genus_results']:
                        st.subheader("Genus-Level Testing")
                        genus_data = specificity_results['genus_results']
                        st.write(f"**Tested Against:** {', '.join(genus_data['tested_against'])}")
                        st.write(f"**Regions Passing:** {len(genus_data['specific_regions'])}")
                        
                        if genus_data['specific_regions']:
                            avg_genus_spec = np.mean([r['specificity_score'] for r in genus_data['specific_regions']])
                            st.metric("Average Genus Specificity", f"{avg_genus_spec:.1%}")
                
                with col2:
                    if 'species_results' in specificity_results and specificity_results['species_results']:
                        st.subheader("Species-Level Testing")
                        species_data = specificity_results['species_results']
                        st.write(f"**Tested Against:** {', '.join(species_data['tested_against'])}")
                        st.write(f"**Regions Passing:** {len(species_data['specific_regions'])}")
                        
                        if species_data['specific_regions']:
                            avg_species_spec = np.mean([r['specificity_score'] for r in species_data['specific_regions']])
                            st.metric("Average Species Specificity", f"{avg_species_spec:.1%}")
                
                # Final regions summary
                final_regions = specificity_results['combined_specific_regions']
                if final_regions:
                    st.subheader("Final Specific Regions")
                    st.write(f"Regions passing all specificity tests: **{len(final_regions)}**")
                    
                    # Create detailed specificity table
                    specificity_table_data = []
                    for i, region in enumerate(final_regions):
                        row = {
                            'Region': i + 1,
                            'Position': f"{region['start']}-{region['end']}",
                            'Conservation': f"{region['conservation_score']:.1%}",
                        }
                        
                        if 'genus_specificity' in region:
                            row['Genus Specificity'] = f"{region['genus_specificity']:.1%}"
                        if 'species_specificity' in region:
                            row['Species Specificity'] = f"{region['species_specificity']:.1%}"
                        if 'combined_specificity' in region:
                            row['Combined Specificity'] = f"{region['combined_specificity']:.1%}"
                        else:
                            row['Overall Specificity'] = f"{region['specificity_score']:.1%}"
                        
                        row['Quality Score'] = f"{region['conservation_score'] * region['specificity_score']:.3f}"
                        specificity_table_data.append(row)
                    
                    specificity_df = pd.DataFrame(specificity_table_data)
                    st.dataframe(specificity_df, use_container_width=True)
        else:
            st.info("No specificity analysis results available. Run comprehensive analysis first.")

    with tab6:  # Export tab (moved from tab4)
        st.header("Export Results")
        
        if st.session_state.primers_designed:
            primers = st.session_state.primers_designed
            
            st.subheader("Download Options")
            
            col1, col2 = st.columns(2)
            
            with col1:
                # Excel export
                if st.button("📊 Download as Excel", type="primary"):
                    excel_data = export_to_excel(primers)
                    st.download_button(
                        label="Click to Download Excel File",
                        data=excel_data,
                        file_name="primer_results.xlsx",
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                    )
            
            with col2:
                # CSV export
                data = []
                for i, primer in enumerate(primers):
                    data.append({
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
                        'Penalty_Score': round(primer.penalty, 4)
                    })
                
                df = pd.DataFrame(data)
                csv = df.to_csv(index=False)
                
                st.download_button(
                    label="📄 Download as CSV",
                    data=csv,
                    file_name="primer_results.csv",
                    mime="text/csv"
                )
            
            # Preview export data
            st.subheader("Export Preview")
            df_preview = pd.DataFrame(data)
            st.dataframe(df_preview, use_container_width=True)
            
            # Custom export options
            st.subheader("Custom Export")
            
            with st.expander("Select columns to export"):
                all_columns = list(df_preview.columns)
                selected_columns = st.multiselect(
                    "Choose columns:",
                    all_columns,
                    default=all_columns
                )
                
                if selected_columns:
                    custom_df = df_preview[selected_columns]
                    custom_csv = custom_df.to_csv(index=False)
                    
                    st.download_button(
                        label="Download Custom CSV",
                        data=custom_csv,
                        file_name="custom_primer_results.csv",
                        mime="text/csv"
                    )
            
            # Primer ordering format
            st.subheader("Primer Ordering Format")
            st.write("Format suitable for ordering from synthesis companies:")
            
            ordering_data = []
            for i, primer in enumerate(primers):
                ordering_data.append({
                    'Name': f"Forward_Primer_{i+1}",
                    'Sequence': primer.forward_seq,
                    'Length': len(primer.forward_seq),
                    'Tm': round(primer.forward_tm, 1)
                })
                ordering_data.append({
                    'Name': f"Reverse_Primer_{i+1}",
                    'Sequence': primer.reverse_seq,
                    'Length': len(primer.reverse_seq),
                    'Tm': round(primer.reverse_tm, 1)
                })
            
            ordering_df = pd.DataFrame(ordering_data)
            st.dataframe(ordering_df, use_container_width=True)
            
            ordering_csv = ordering_df.to_csv(index=False)
            st.download_button(
                label="📋 Download Ordering Format",
                data=ordering_csv,
                file_name="primer_ordering.csv",
                mime="text/csv"
            )
            
        else:
            st.info("No primers to export. Please design primers first.")
    
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
