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
import primer3
import re
from pathlib import Path
import base64

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
            st.error(f"Error fetching sequence info for {seq_id}: {e}")
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

# Streamlit App Functions
def init_session_state():
    """Initialize session state variables"""
    if 'primers_designed' not in st.session_state:
        st.session_state.primers_designed = []
    if 'current_sequence' not in st.session_state:
        st.session_state.current_sequence = ""
    if 'sequence_info' not in st.session_state:
        st.session_state.sequence_info = {}

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
    tab1, tab2, tab3, tab4 = st.tabs(["📝 Input", "🔬 Results", "📊 Analysis", "💾 Export"])
    
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
            
            col1, col2 = st.columns([2, 1])
            with col1:
                # Initialize organism name from session state if available
                default_organism = st.session_state.get('organism_name', '')
                
                organism_name = st.text_input("Enter organism name:", 
                                            value=default_organism,
                                            placeholder="e.g., Fusarium oxysporum, Tetranychus urticae, Beet curly top virus")
                
                # Clear session state after using it
                if 'organism_name' in st.session_state:
                    del st.session_state.organism_name
                
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
            
            if st.button("Search Organism Genomes", type="primary"):
                if not email:
                    st.error("❌ **Email Required**: Please enter an email address in the sidebar to access NCBI databases. You can use any valid email address or click 'Use test email' for demo purposes.")
                elif not organism_name:
                    st.error("Please provide an organism name")
                else:
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
                            
                            if genome_ids:
                                st.success(f"Found {len(genome_ids)} genome entries!")
                                
                                # Display found genomes
                                st.subheader("Available Genomes")
                                genome_info = []
                                
                                for i, genome_id in enumerate(genome_ids):
                                    with st.spinner(f"Fetching genome {i+1}/{len(genome_ids)}..."):
                                        info = ncbi.fetch_sequence_info(genome_id, database="genome")
                                        if info:
                                            genome_info.append({
                                                'ID': genome_id,
                                                'Description': info.get('description', 'N/A'),
                                                'Length': info.get('length', 'N/A'),
                                                'Organism': info.get('organism', 'N/A')
                                            })
                                
                                if genome_info:
                                    genome_df = pd.DataFrame(genome_info)
                                    st.dataframe(genome_df, use_container_width=True)
                                    
                                    # Let user select a genome
                                    selected_genome = st.selectbox(
                                        "Select a genome to design primers for:",
                                        range(len(genome_info)),
                                        format_func=lambda x: f"{genome_info[x]['ID']} - {genome_info[x]['Description'][:100]}..."
                                    )
                                    
                                    if st.button("Design Primers for Selected Genome", type="primary"):
                                        with st.spinner("Fetching genome sequence and designing primers..."):
                                            # Get the full genome sequence
                                            genome_id = genome_info[selected_genome]['ID']
                                            sequence = ncbi.fetch_sequence(genome_id, database="genome")
                                            
                                            if sequence:
                                                # For very large genomes, we might want to limit the sequence length
                                                if len(sequence) > 1000000:  # 1MB limit
                                                    st.warning(f"Genome is very large ({len(sequence):,} bp). Using first 1MB for primer design.")
                                                    sequence = sequence[:1000000]
                                                
                                                seq_info = ncbi.fetch_sequence_info(genome_id, database="genome")
                                                st.session_state.sequence_info = seq_info
                                                st.session_state.current_sequence = sequence
                                                
                                                # Design primers
                                                primers = designer.design_primers(sequence, custom_params=custom_params)
                                                st.session_state.primers_designed = primers
                                                
                                                if primers:
                                                    st.success(f"Successfully designed {len(primers)} primer pairs!")
                                                else:
                                                    st.warning("No suitable primers found with current parameters")
                                            else:
                                                st.error("Failed to fetch genome sequence")
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
                                
                                if primers:
                                    st.success(f"Successfully designed {len(primers)} primer pairs!")
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
                                
                                if all_primers:
                                    st.success(f"Found {len(all_primers)} primer pairs from {len(all_info)} sequences!")
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
                            
                            st.session_state.primers_designed = primers
                            st.session_state.current_sequence = clean_seq
                            st.session_state.sequence_info = {
                                "length": len(clean_seq),
                                "description": "User-provided sequence"
                            }
                            
                            if primers:
                                st.success(f"Successfully designed {len(primers)} primer pairs!")
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
                            
                            st.session_state.primers_designed = primers
                            st.session_state.current_sequence = clean_seq
                            st.session_state.sequence_info = {
                                "length": len(clean_seq),
                                "description": f"Uploaded file: {uploaded_file.name}"
                            }
                            
                            if primers:
                                st.success(f"Successfully designed {len(primers)} primer pairs!")
                            else:
                                st.warning("No suitable primers found")
                        except Exception as e:
                            st.error(f"Error processing file: {e}")
    
    with tab2:
        st.header("Primer Design Results")
        
        if st.session_state.primers_designed:
            primers = st.session_state.primers_designed
            
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
                data.append({
                    'Pair': i + 1,
                    'Forward Sequence': primer.forward_seq,
                    'Reverse Sequence': primer.reverse_seq,
                    'Forward Tm': f"{primer.forward_tm:.1f}°C",
                    'Reverse Tm': f"{primer.reverse_tm:.1f}°C",
                    'Product Size': f"{primer.product_size} bp",
                    'Forward GC%': f"{primer.gc_content_f:.1f}%",
                    'Reverse GC%': f"{primer.gc_content_r:.1f}%",
                    'Penalty': f"{primer.penalty:.3f}"
                })
            
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
        else:
            st.info("No primers designed yet. Please use the Input tab to design primers.")
    
    with tab3:
        st.header("Primer Analysis")
        
        if st.session_state.primers_designed:
            primers = st.session_state.primers_designed
            
            # Create visualizations
            fig = create_primer_visualization(primers)
            if fig:
                st.plotly_chart(fig, use_container_width=True)
            
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
    
    with tab4:
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
