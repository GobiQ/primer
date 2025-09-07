def get_related_organisms(target_organism, max_organisms=100):
    """Get comprehensive list of related organisms for specificity testing (100+ organisms)"""
    organism_lower = target_organism.lower()
    
    # Comprehensive organism database with 100+ organisms across all major taxonomic groups
    comprehensive_organisms = {
        # Fungi (40+ organisms)
        'fungi': [
            'Fusarium oxysporum', 'Fusarium solani', 'Fusarium graminearum', 'Fusarium verticillioides',
            'Aspergillus niger', 'Aspergillus flavus', 'Aspergillus fumigatus', 'Aspergillus terreus',
            'Penicillium chrysogenum', 'Penicillium expansum', 'Penicillium digitatum', 'Penicillium italicum',
            'Botrytis cinerea', 'Sclerotinia sclerotiorum', 'Alternaria alternata', 'Alternaria solani',
            'Rhizoctonia solani', 'Pythium ultimum', 'Phytophthora infestans', 'Phytophthora capsici',
            'Trichoderma harzianum', 'Trichoderma viride', 'Trichoderma reesei', 'Trichoderma atroviride',
            'Erysiphe necator', 'Puccinia graminis', 'Ustilago maydis', 'Magnaporthe oryzae',
            'Colletotrichum gloeosporioides', 'Verticillium dahliae', 'Verticillium albo-atrum',
            'Cercospora beticola', 'Septoria tritici', 'Mycosphaerella graminicola', 'Cladosporium fulvum',
            'Monilinia fructicola', 'Monilinia laxa', 'Monilinia fructigena', 'Sclerotium rolfsii',
            'Armillaria mellea', 'Ganoderma lucidum', 'Pleurotus ostreatus', 'Agaricus bisporus'
        ],
        
        # Bacteria (30+ organisms)
        'bacteria': [
            'Pseudomonas syringae', 'Pseudomonas aeruginosa', 'Pseudomonas fluorescens', 'Pseudomonas putida',
            'Xanthomonas campestris', 'Xanthomonas oryzae', 'Xanthomonas axonopodis', 'Xanthomonas citri',
            'Erwinia amylovora', 'Erwinia carotovora', 'Erwinia chrysanthemi', 'Erwinia pyrifoliae',
            'Ralstonia solanacearum', 'Agrobacterium tumefaciens', 'Rhizobium leguminosarum', 'Sinorhizobium meliloti',
            'Bacillus subtilis', 'Bacillus thuringiensis', 'Bacillus cereus', 'Bacillus anthracis',
            'Staphylococcus aureus', 'Streptococcus pneumoniae', 'Escherichia coli', 'Salmonella enterica',
            'Listeria monocytogenes', 'Clostridium botulinum', 'Mycobacterium tuberculosis', 'Corynebacterium diphtheriae',
            'Streptomyces coelicolor', 'Actinomyces israelii', 'Propionibacterium acnes', 'Lactobacillus acidophilus'
        ],
        
        # Viruses (25+ organisms)
        'viruses': [
            'Tobacco mosaic virus', 'Tomato mosaic virus', 'Cucumber mosaic virus', 'Potato virus X',
            'Potato virus Y', 'Beet curly top virus', 'Arabis mosaic virus', 'Alfalfa mosaic virus',
            'Lettuce chlorosis virus', 'Cannabis cryptic virus', 'Hop latent viroid', 'Potato spindle tuber viroid',
            'Influenza A virus', 'Coronavirus', 'Rhinovirus', 'Adenovirus', 'Respiratory syncytial virus',
            'MERS-CoV', 'Human coronavirus 229E', 'Human coronavirus OC43', 'Parainfluenza virus',
            'Measles virus', 'Mumps virus', 'Rubella virus', 'Varicella-zoster virus', 'Herpes simplex virus'
        ],
        
        # Arthropods (25+ organisms)
        'arthropods': [
            'Tetranychus urticae', 'Panonychus ulmi', 'Oligonychus ilicis', 'Aculops lycopersici',
            'Bemisia tabaci', 'Trialeurodes vaporariorum', 'Aphis gossypii', 'Myzus persicae',
            'Diaphorina citri', 'Frankliniella occidentalis', 'Thrips tabaci', 'Thrips palmi',
            'Scirtothrips dorsalis', 'Polyphagotarsonemus latus', 'Tarsonemus pallidus',
            'Drosophila melanogaster', 'Musca domestica', 'Lucilia sericata', 'Calliphora vicina',
            'Tribolium castaneum', 'Sitophilus oryzae', 'Plodia interpunctella', 'Ephestia kuehniella',
            'Helicoverpa armigera', 'Spodoptera frugiperda', 'Plutella xylostella', 'Pieris rapae'
        ],
        
        # Plants (20+ organisms)
        'plants': [
            'Arabidopsis thaliana', 'Oryza sativa', 'Zea mays', 'Triticum aestivum',
            'Hordeum vulgare', 'Solanum lycopersicum', 'Solanum tuberosum', 'Capsicum annuum',
            'Nicotiana tabacum', 'Gossypium hirsutum', 'Glycine max', 'Phaseolus vulgaris',
            'Vitis vinifera', 'Malus domestica', 'Prunus persica', 'Citrus sinensis',
            'Brassica oleracea', 'Daucus carota', 'Lactuca sativa', 'Spinacia oleracea',
            'Helianthus annuus', 'Beta vulgaris', 'Saccharum officinarum', 'Manihot esculenta'
        ],
        
        # Nematodes (10+ organisms)
        'nematodes': [
            'Caenorhabditis elegans', 'Meloidogyne incognita', 'Meloidogyne javanica', 'Meloidogyne hapla',
            'Heterodera glycines', 'Globodera rostochiensis', 'Pratylenchus penetrans', 'Radopholus similis',
            'Ditylenchus dipsaci', 'Aphelenchoides fragariae', 'Bursaphelenchus xylophilus', 'Steinernema carpocapsae'
        ]
    }
    
    # Try to find specific matches first
    specific_matches = []
    for category, organisms in comprehensive_organisms.items():
        for organism in organisms:
            organism_lower_check = organism.lower()
            # Check for genus or species matches
            if any(part in organism_lower_check for part in organism_lower.split()):
                specific_matches.append(organism)
    
    # If we found specific matches, return them plus related organisms
    if specific_matches:
        result = specific_matches[:]
        # Add organisms from the same category
        for category, organisms in comprehensive_organisms.items():
            if any(match.lower() in [org.lower() for org in organisms] for match in specific_matches):
                result.extend([org for org in organisms if org not in result])
                break
        return result[:max_organisms]
    
    # If no specific matches, return a diverse set from all categories
    diverse_set = []
    for category, organisms in comprehensive_organisms.items():
        diverse_set.extend(organisms[:max_organisms//len(comprehensive_organisms)])
    
    # Add more organisms to reach the target number
    remaining_needed = max_organisms - len(diverse_set)
    if remaining_needed > 0:
        all_organisms = []
        for organisms in comprehensive_organisms.values():
            all_organisms.extend(organisms)
        
        # Add random selection from remaining organisms
        import random
        remaining_organisms = [org for org in all_organisms if org not in diverse_set]
        if remaining_organisms:
            diverse_set.extend(random.sample(remaining_organisms, min(remaining_needed, len(remaining_organisms))))
    
    return diverse_set[:max_organisms]

def check_session_state_validity():
    """Check if session state has valid data"""
    # Ensure all required keys exist with safe defaults
    primers = st.session_state.get('primers_designed', [])
    sequence = st.session_state.get('current_sequence', '')
    seq_info = st.session_state.get('sequence_info', {})
    
    has_primers = bool(primers)
    has_sequence = bool(sequence)
    has_seq_info = bool(seq_info)
    
    return {
        'has_primers': has_primers,
        'has_sequence': has_sequence,
        'has_seq_info': has_seq_info,
        'primer_count': len(primers),
        'sequence_length': len(sequence)
    }

def get_gene_priority(category):
    """Get priority level for gene category"""
    category_lower = category.lower()
    
    # High priority categories
    if any(keyword in category_lower for keyword in ['essential', 'housekeeping', 'actin', 'tubulin', 'rna polymerase', 'ribosomal']):
        return "High"
    
    # Medium priority categories
    elif any(keyword in category_lower for keyword in ['pathogenicity', 'virulence', 'effector', 'resistance', 'detoxification']):
        return "Medium"
    
    # Low priority categories
    else:
        return "Low"

def get_gene_use_recommendation(category):
    """Get usage recommendation for gene category"""
    category_lower = category.lower()
    
    if 'essential' in category_lower or 'housekeeping' in category_lower:
        return "Recommended for reliable detection and quantification. These genes are highly conserved and expressed constitutively."
    
    elif 'pathogenicity' in category_lower or 'virulence' in category_lower:
        return "Useful for detecting pathogenic strains and studying disease mechanisms. May not be present in all isolates."
    
    elif 'resistance' in category_lower:
        return "Important for monitoring resistance development and screening for resistant strains. Critical for management decisions."
    
    elif 'detoxification' in category_lower:
        return "Useful for studying pesticide resistance and detoxification mechanisms. Important for resistance management."
    
    elif 'secondary metabolite' in category_lower:
        return "Useful for detecting toxin-producing strains and studying secondary metabolism. May be strain-specific."
    
    elif 'cell wall' in category_lower:
        return "Good for structural studies and cell wall targeting. Generally well-conserved across strains."
    
    elif 'development' in category_lower:
        return "Useful for studying life cycle stages and development processes. May be stage-specific."
    
    else:
        return "General purpose gene category. Consider specific research objectives when selecting targets."

def get_organism_suggestions_with_gene_targets():
    """Get agricultural pest and pathogen suggestions organized by category with comprehensive gene targets"""
    return {
        "🍄 Fungal Pathogens": {
            "Fusarium species": [
                ("Fusarium wilt", "Fusarium oxysporum", {
                    "Essential genes": ["ACT1 (actin)", "TUB2 (tubulin)", "EF1A (elongation factor)", "RPB2 (RNA polymerase)", "LSU (large subunit rRNA)"],
                    "Pathogenicity genes": ["SIX1-14 (secreted in xylem)", "FTF1 (transcription factor)", "FMK1 (MAPK)", "SGE1 (cutinase)", "PEL1 (pectate lyase)"],
                    "Secondary metabolite genes": ["FUM1 (fumonisin biosynthesis)", "TRI5 (trichothecene biosynthesis)", "PKS4 (polyketide synthase)", "BIK1 (bikaverin)", "FUS1 (fusarin)"],
                    "Cell wall genes": ["CHI1 (chitinase)", "GEL1 (β-1,3-glucanase)", "CHS1 (chitin synthase)", "FKS1 (β-1,3-glucan synthase)", "PMI1 (mannose-6-phosphate isomerase)"],
                    "Resistance targets": ["CYP51 (sterol 14α-demethylase)", "SDH (succinate dehydrogenase)", "QoI (cytochrome bc1)", "MBC (β-tubulin)", "DMI (sterol biosynthesis)"]
                }),
                ("Fusarium head blight", "Fusarium graminearum", {
                    "Essential genes": ["ACT1 (actin)", "TUB1 (α-tubulin)", "TUB2 (β-tubulin)", "EF1A (elongation factor)", "RPB2 (RNA polymerase II)"],
                    "Trichothecene pathway": ["TRI1-TRI16 (trichothecene cluster)", "TRI101 (trichothecene resistance)", "TRI201 (regulatory)", "TRI301 (transport)", "TRI401 (efflux pump)"],
                    "DON biosynthesis": ["TRI4 (trichodiene oxygenase)", "TRI11 (isotrichodermin C-15 hydroxylase)", "TRI3 (15-O-acetyltransferase)", "TRI7 (C-4 hydroxylase)", "TRI13 (isotrichodermol dehydrogenase)"],
                    "Pathogenicity factors": ["MGV1 (major facilitator)", "TRI12 (efflux pump)", "FgCHS1 (chitin synthase)", "FgPKS4 (polyketide synthase)", "FgLAE1 (global regulator)"],
                    "Stress response": ["HOG1 (osmoregulation)", "FgOS2 (osmotic stress)", "FgSLT2 (cell wall integrity)", "ATF1 (oxidative stress)", "FgAP1 (AP-1 transcription factor)"]
                }),
                ("Fusarium crown rot", "Fusarium culmorum", {
                    "Housekeeping genes": ["ACT1", "TUB2", "EF1A", "RPB2", "LSU rRNA"],
                    "Pathogenicity": ["CUT1 (cutinase)", "PEL2 (pectinase)", "XYL1 (xylanase)", "CEL1 (cellulase)", "AMY1 (amylase)"],
                    "Secondary metabolites": ["ZEA1-2 (zearalenone)", "TRI cluster", "AUR1 (aurofusarin)", "CUL1 (culmorin)", "BUT1 (butenolide)"],
                    "Survival genes": ["HSP70 (heat shock)", "SOD1 (superoxide dismutase)", "CAT1 (catalase)", "GPX1 (glutathione peroxidase)", "TRX1 (thioredoxin)"]
                }),
                ("Fusarium root rot", "Fusarium solani", {
                    "Essential genes": ["ACT1", "TUB2", "EF1A", "RPB2", "HSP70"],
                    "Pathogenicity factors": ["FSOL1-10 (F. solani specific)", "CUT1-5 (cutinases)", "PEL1-3 (pectate lyases)", "XYL1-2 (xylanases)", "CEL1-2 (cellulases)"],
                    "Secondary metabolites": ["FUM1-3 (fumonisin)", "TRI1-16 (trichothecene)", "ZEA1-2 (zearalenone)", "FUS1-5 (fusarin)", "BIK1-3 (bikaverin)"],
                    "Resistance mechanisms": ["CYP51A1-B1", "SDH1-4", "ABC1-20", "MFS1-15", "GST1-10"]
                }),
                ("Fusarium ear rot", "Fusarium proliferatum", {
                    "Essential genes": ["ACT1", "TUB2", "EF1A", "RPB2", "ITS1-2"],
                    "Fumonisin biosynthesis": ["FUM1-21 (fumonisin cluster)", "FUM8 (polyketide synthase)", "FUM6 (aminotransferase)", "FUM3 (C-5 oxygenase)", "FUM19 (transporter)"],
                    "Pathogenicity": ["FPR1-10 (F. proliferatum specific)", "CUT1-3", "PEL1-2", "XYL1", "CEL1"],
                    "Host interaction": ["HOST1-5 (host-specific)", "ADH1-3 (adhesion)", "INV1-2 (invasion)", "COL1-3 (colonization)"]
                })
            ],
            "Other fungi": [
                ("Gray mold", "Botrytis cinerea", {
                    "Essential genes": ["ACT1", "TUB2", "EF1A", "RPB2", "HSP70"],
                    "Pathogenicity factors": ["BCG1 (α-galactosidase)", "BMP1 (metalloprotease)", "BCP1 (cerato-platanin)", "BOA1 (botrydial)", "BCR1 (ABC transporter)"],
                    "Cell wall degrading": ["BcPG1-6 (polygalacturonases)", "BcPME1 (pectin methylesterase)", "BcXYL1 (xylanase)", "BcCEL1 (cellulase)", "BcCUT1 (cutinase)"],
                    "Secondary metabolites": ["BOT1-5 (botrydial cluster)", "DHN1 (1,8-dihydroxynaphthalene)", "PKS1-13 (polyketide synthases)", "NPS1-6 (nonribosomal peptide synthetases)"],
                    "Resistance mechanisms": ["ABC1-50 (ABC transporters)", "MFS1-20 (major facilitator superfamily)", "CYP1-20 (cytochrome P450s)", "GST1-10 (glutathione S-transferases)"]
                }),
                ("White mold", "Sclerotinia sclerotiorum", {
                    "Essential genes": ["ACT1", "TUB2", "EF1A", "RPB2", "CAL1"],
                    "Pathogenicity factors": ["SSPG1 (polygalacturonase)", "SsPME1 (pectin methylesterase)", "SsCUT1 (cutinase)", "SsNEP1 (necrosis-inducing protein)", "SsCP1 (cysteine protease)"],
                    "Oxalic acid pathway": ["OAH1 (oxaloacetate hydrolase)", "OXA1 (oxalate oxidase)", "OXD1 (oxalate decarboxylase)", "OMT1 (oxalate metabolism)", "OXT1 (oxalate transporter)"],
                    "Sclerotia formation": ["SMK1 (sclerotial development)", "SLT1 (sclerotial maturation)", "SCL1 (sclerotial pigmentation)", "MBF1 (melanin biosynthesis)", "LAC1 (laccase)"]
                }),
                ("Late blight", "Phytophthora infestans", {
                    "Essential genes": ["ACT1", "TUB1", "EF1A", "RPB1", "COX1 (cytochrome c oxidase)"],
                    "Effector genes": ["AVR1-11 (avirulence genes)", "PexRD1-36 (Phytophthora expressed)", "Pi02860-Pi17316 (candidate effectors)", "RXLR1-100 (effector motif)", "CRN1-50 (crinkler effectors)"],
                    "RXLR effectors": ["PiAVR3a", "PiAVR4", "PiAVRblb1", "PiAVRblb2", "PiAVRvnt1", "Pi02860", "Pi03192", "Pi04314", "Pi07569", "Pi09585"],
                    "Cell wall degrading": ["PiCBEL1-15 (cellulose-binding elicitor lectin)", "PiGH12-1 (endoglucanase)", "PiPL1-20 (pectate lyases)", "PiXEG1 (xyloglucanase)", "PiPG1-5 (polygalacturonases)"],
                    "Pathogenicity factors": ["PiINF1 (infestin)", "PiNPP1 (necrosis-inducing protein)", "PiEPI1-10 (epidermis-specific)", "PiHAM34 (haustorial membrane)", "PiMCF1 (mycelium-cyst transition)"]
                }),
                ("Powdery mildew", "Erysiphe necator", {
                    "Essential genes": ["ACT1", "TUB2", "EF1A", "RPB2", "ITS1-2"],
                    "Pathogenicity": ["ENH1 (haustorium formation)", "ENC1 (conidiophore development)", "ENS1 (spore germination)", "ENP1 (penetration)", "ENA1 (appressorium formation)"],
                    "Effectors": ["CSEP1-100 (candidate secreted effector proteins)", "ENE1-50 (E. necator effectors)", "AVR1-10 (avirulence candidates)", "HAU1-20 (haustorial expressed)"],
                    "Sterol biosynthesis": ["CYP51A1", "CYP51B1", "ERG1 (squalene epoxidase)", "ERG7 (lanosterol synthase)", "ERG11 (sterol 14α-demethylase)"]
                }),
                ("Rust disease", "Puccinia graminis", {
                    "Essential genes": ["ACT1", "TUB1", "EF1A", "RPB1", "COI"],
                    "Pathogenicity factors": ["PG1-50 (P. graminis specific)", "AVR1-20 (avirulence genes)", "EFF1-30 (effector proteins)", "SEC1-10 (secreted proteins)", "HST1-5 (host-specific toxins)"],
                    "Life cycle genes": ["TEL1-5 (teliospore formation)", "BAS1-3 (basidiospore)", "PYC1-2 (pycniospore)", "AEC1-3 (aeciospore)", "URE1-5 (urediniospore)"],
                    "Host interaction": ["HOST1-10 (host recognition)", "PEN1-5 (penetration)", "COL1-8 (colonization)", "SPR1-5 (sporulation)"]
                }),
                ("Smut disease", "Ustilago maydis", {
                    "Essential genes": ["ACT1", "TUB1", "EF1A", "RPB1", "HSP70"],
                    "Pathogenicity factors": ["UM1-100 (U. maydis specific)", "EFF1-50 (effector proteins)", "SEC1-20 (secreted proteins)", "CWP1-10 (cell wall proteins)", "ENZ1-15 (enzymes)"],
                    "Mating and development": ["MAT1-2 (mating type)", "B1-6 (b locus)", "A1-4 (a locus)", "CLP1 (clamp connection)", "DIC1-3 (dikaryon formation)"],
                    "Host specificity": ["HOST1-8 (host recognition)", "MAI1-5 (maize interaction)", "TIS1-3 (tissue specificity)", "SIZ1-2 (size control)"]
                })
            ]
        },
        
        "🐛 Insect Pests": {
            "Mites": [
                ("Two-spotted spider mite", "Tetranychus urticae", {
                    "Essential genes": ["ACT1 (actin)", "TUB1 (tubulin)", "EF1A (elongation factor)", "RPL32 (ribosomal protein L32)", "RPS3 (ribosomal protein S3)"],
                    "Detoxification genes": ["CYP1-100 (cytochrome P450s)", "GST1-30 (glutathione S-transferases)", "EST1-20 (esterases)", "UGT1-15 (UDP-glucuronosyltransferases)", "ABC1-50 (ABC transporters)"],
                    "Acaricide resistance": ["AChE (acetylcholinesterase)", "VGSC (voltage-gated sodium channel)", "RDL (GABA receptor)", "nAChR (nicotinic acetylcholine receptor)", "GluCl (glutamate-gated chloride channel)"],
                    "Development genes": ["JH1-3 (juvenile hormone)", "ECR (ecdysone receptor)", "USP (ultraspiracle)", "E74 (ecdysone response)", "BR-C (broad complex)"],
                    "Reproduction genes": ["VIT1-6 (vitellogenin)", "VTG1-3 (vitellogenin)", "CHR (chorion)", "EGG1-5 (egg development)", "EMB1-10 (embryogenesis)"]
                }),
                ("European red mite", "Panonychus ulmi", {
                    "Essential genes": ["ACT1", "TUB1", "EF1A", "RPL32", "COI (cytochrome oxidase I)"],
                    "Resistance genes": ["CYP1-50", "GST1-20", "EST1-15", "P450-1-25", "MFO1-10 (mixed function oxidases)"],
                    "Cold tolerance": ["HSP70 (heat shock protein)", "AFP1-3 (antifreeze proteins)", "TRE1 (trehalose)", "GLY1-2 (glycerol)", "CRY1-2 (cryoprotectants)"],
                    "Diapause genes": ["DIA1-5 (diapause-associated)", "CLK (clock)", "PER (period)", "TIM (timeless)", "CYC (cycle)"]
                }),
                ("Broad mite", "Polyphagotarsonemus latus", {
                    "Essential genes": ["ACT1", "TUB1", "EF1A", "RPL32", "COI"],
                    "Host range genes": ["HOST1-10 (host specificity)", "DET1-5 (detoxification)", "ADH1-3 (adhesion)", "PEN1-2 (penetration)", "COL1-3 (colonization)"],
                    "Resistance mechanisms": ["CYP1-30", "GST1-15", "EST1-10", "ABC1-20", "MFS1-10"],
                    "Development": ["JH1-2", "ECR", "USP", "E74", "BR-C"]
                }),
                ("Russet mite", "Aculops lycopersici", {
                    "Essential genes": ["ACT1", "TUB1", "EF1A", "RPL32", "COI"],
                    "Tomato adaptation": ["TOM1-5 (tomato-specific)", "LYC1-3 (lycopersicon)", "SOL1-2 (solanum)", "ADH1-2 (adhesion)", "PEN1 (penetration)"],
                    "Resistance genes": ["CYP1-25", "GST1-12", "EST1-8", "ABC1-15", "MFS1-8"],
                    "Feeding behavior": ["FED1-5 (feeding)", "SAL1-3 (salivary)", "GUT1-3 (gut-specific)", "DIG1-2 (digestion)"]
                })
            ],
            "Sucking insects": [
                ("Silverleaf whitefly", "Bemisia tabaci", {
                    "Essential genes": ["ACT1", "TUB1", "EF1A", "COI", "16S rRNA"],
                    "Insecticide resistance": ["CYP6CM1", "CYP4C64", "CYP4G61", "CYP4G70", "GST1-15", "EST1-10", "ABC1-30", "nAChR (nicotinic receptor)", "VGSC (sodium channel)"],
                    "Biotype markers": ["mtCOI (mitochondrial COI)", "ITS1", "16S rRNA", "28S rRNA", "RAPD markers"],
                    "Virus transmission": ["GroEL (endosymbiont)", "HSP70", "cyclophilin", "importin-α", "karyopherin", "VP1-4 (viral proteins)"],
                    "Endosymbiont genes": ["Portiera (P-endosymbiont)", "Hamiltonella", "Arsenophonus", "Cardinium", "Wolbachia", "Rickettsia"]
                }),
                ("Greenhouse whitefly", "Trialeurodes vaporariorum", {
                    "Essential genes": ["ACT1", "TUB1", "EF1A", "COI", "COII"],
                    "Development markers": ["JH (juvenile hormone)", "ECR (ecdysone receptor)", "CHI1-3 (chitinase)", "CHS1-2 (chitin synthase)", "TRE1-2 (trehalase)"],
                    "Host plant interaction": ["SUC1-3 (sucrase)", "APH1-2 (aphid-like stylet)", "SAL1-3 (salivary)", "GUT1-5 (gut-specific)", "PHE1-3 (phenoloxidase)"]
                }),
                ("Green peach aphid", "Myzus persicae", {
                    "Essential genes": ["ACT1", "TUB1", "EF1A", "COI", "COII"],
                    "Insecticide resistance": ["MACE (modified acetylcholinesterase)", "kdr (knockdown resistance)", "RDL (GABA receptor)", "CYP6CY3", "CYP4CJ1", "GST1-10", "EST1-8", "ABC1-20"],
                    "Morph determination": ["APH1 (apterous)", "WIN1 (wingless)", "VG1 (vestigial)", "DSX1 (doublesex)", "FRU1 (fruitless)"],
                    "Virus transmission": ["PLRV (potato leafroll virus)", "PVY (potato virus Y)", "CMV (cucumber mosaic virus)", "receptor proteins", "helper factors"],
                    "Endosymbiont": ["Buchnera aphidicola", "trpA-E (tryptophan synthesis)", "aroA-Q (aromatic amino acids)", "ilv (isoleucine-valine)", "leu (leucine)", "phe (phenylalanine)"]
                }),
                ("Cotton aphid", "Aphis gossypii", {
                    "Essential genes": ["ACT1", "TUB1", "EF1A", "COI", "16S rRNA"],
                    "Host specialization": ["HSP1-10 (host selection)", "DET1-5 (detoxification)", "FED1-3 (feeding behavior)", "GOT1-5 (gossypol tolerance)", "TAN1-3 (tannin tolerance)"],
                    "Polyphenism": ["WIN1-5 (wing development)", "ALT1-3 (alternate morph)", "ENV1-5 (environmental response)", "DEN1-3 (density-dependent)", "PHO1-3 (photoperiod)"]
                }),
                ("Asian citrus psyllid", "Diaphorina citri", {
                    "Essential genes": ["ACT1", "TUB1", "EF1A", "COI", "16S rRNA"],
                    "CLas transmission": ["CLas1-10 (Candidatus Liberibacter asiaticus)", "HSP1-5 (heat shock)", "SEC1-5 (secretion)", "TRP1-3 (transporter)", "BIND1-3 (binding)"],
                    "Citrus adaptation": ["CIT1-5 (citrus-specific)", "RUT1-3 (rutaceae)", "FLAV1-2 (flavonoid)", "TERP1-2 (terpene)", "OIL1-2 (essential oil)"],
                    "Development": ["JH1-2", "ECR", "USP", "E74", "BR-C"]
                })
            ],
            "Thrips": [
                ("Western flower thrips", "Frankliniella occidentalis", {
                    "Essential genes": ["ACT1", "TUB1", "EF1A", "COI", "18S rRNA"],
                    "Insecticide resistance": ["CYP1-30", "GST1-15", "EST1-10", "AChE", "VGSC", "RDL", "GluCl", "nAChR"],
                    "Virus transmission": ["TSWV (tomato spotted wilt virus)", "INSV (impatiens necrotic spot virus)", "receptor1-3", "helper1-2", "vector1-5"],
                    "Host preference": ["ORN1-10 (olfactory)", "GRN1-5 (gustatory)", "CHE1-8 (chemoreception)", "HOST1-5 (host selection)", "FED1-3 (feeding)"],
                    "Development": ["JHE (juvenile hormone esterase)", "JH1-3", "ECR", "USP", "E74", "BR-C", "CHI1-2", "CHS1-2"]
                }),
                ("Onion thrips", "Thrips tabaci", {
                    "Essential genes": ["ACT1", "TUB1", "EF1A", "COI", "COII"],
                    "Arrhenotoky": ["SEX1-5 (sex determination)", "HAP1-3 (haplodiploidy)", "CSD (complementary sex determiner)", "FEM1-2 (feminizer)", "DSX"],
                    "Onion adaptation": ["ALL1-5 (alliin tolerance)", "SUL1-3 (sulfur metabolism)", "ONI1-5 (onion-specific)", "LAC1-2 (lachrymatory factor)", "THI1-3 (thiosulfinate)"]
                })
            ]
        },

        "🦠 Bacterial Pathogens": {
            "Erwinia species": [
                ("Fire blight", "Erwinia amylovora", {
                    "Essential genes": ["16S rRNA", "23S rRNA", "gyrA", "gyrB", "rpoB", "rpoD"],
                    "Pathogenicity factors": ["hrpA-W (type III secretion)", "dspA/E (disease-specific)", "amsA-K (amylovoran synthesis)", "rcsA-C (regulation)", "avrRpt2EA"],
                    "Type III effectors": ["eop1-4 (effector of pathogenicity)", "hopA1EA", "hopC1EA", "avrE1", "dspA/E", "eop2", "eop3", "eop4"],
                    "Exopolysaccharide": ["amsA-K (amylovoran)", "galE", "galF", "ugpA-E", "wza-c", "gmd", "fcl", "rmlA-D"],
                    "Virulence regulation": ["rcsA-C", "kdpD/E", "ompR/envZ", "phoP/Q", "rpoS", "fur", "crp", "ihfA/B"]
                })
            ],
            "Pseudomonas species": [
                ("Bacterial speck", "Pseudomonas syringae", {
                    "Essential genes": ["16S rRNA", "gyrA", "gyrB", "rpoB", "rpoD"],
                    "Type III secretion": ["hrpA-U", "hrcC", "hrcJ", "hrcN", "hrcQ-U", "hrpZ", "hrpA", "hrpG", "hrpL"],
                    "Effector proteins": ["avrE1", "hopA1", "hopM1", "hopZ1-3", "avrPto", "avrPtoB", "avrRpm1", "avrRpt2", "hopF2", "hopG1"],
                    "Toxin production": ["coronatine biosynthesis", "cmaA-T", "cfa1-9", "syringomycin", "syrA-E", "syringopeptin", "sypA-C"],
                    "Ice nucleation": ["inaA-Z", "ice nucleation proteins", "ina genes", "frost injury", "epiphytic survival"]
                })
            ],
            "Xanthomonas species": [
                ("Bacterial blight", "Xanthomonas campestris", {
                    "Essential genes": ["16S rRNA", "gyrA", "gyrB", "rpoB", "atpD"],
                    "Xanthan production": ["gumB-M", "xanA-B", "xanthan gum", "EPS production", "biofilm formation"],
                    "Type III system": ["hrpA-F", "hrcC", "hrcJ", "hrcN", "hrcQ-U", "hpaA-F", "hpaP"],
                    "TAL effectors": ["avrBs3 family", "PthA1-4", "TAL1-20", "transcription activator-like", "DNA-binding repeats"],
                    "Pathogenicity factors": ["rpfA-G (regulation)", "clp (cellulase)", "prt (protease)", "man (mannanase)", "xps (secretion)"]
                })
            ],
            "Ralstonia species": [
                ("Bacterial wilt", "Ralstonia solanacearum", {
                    "Essential genes": ["16S rRNA", "gyrA", "gyrB", "rpoB", "egl"],
                    "Type III secretion": ["hrpA-Y", "hrcC", "hrcJ", "hrcN", "hrcQ-U", "popA-F", "ripA-Z"],
                    "Rip effectors": ["ripA1-Z9", "popP1-2", "popF1-4", "avrA", "ripG1-7", "ripP1-2", "ripE1", "ripF1"],
                    "EPS production": ["epsA-R", "exopolysaccharide", "biofilm", "wilt induction", "vascular plugging"],
                    "Phylotypes": ["phylotype I-IV", "sequevars", "biovars 1-5", "race 1-5", "geographic distribution"]
                })
            ],
            "Common bacteria": [
                ("Crown gall", "Agrobacterium tumefaciens", {
                    "Essential genes": ["16S rRNA", "gyrA", "gyrB", "rpoB", "atpD"],
                    "Ti plasmid genes": ["virA-G (virulence)", "tms1-2 (tumor morphology)", "tmr (tumor morphology)", "tss (tumor size)", "ipt (isopentenyl transferase)"],
                    "T-DNA transfer": ["T-DNA border sequences", "overdrive sequences", "virD1-4 (T-DNA processing)", "virE1-2 (T-DNA protection)", "virB1-11 (T4SS)"],
                    "Opine catabolism": ["occ (octopine)", "nop (nopaline)", "agr (agropine)", "man (mannopine)", "suc (succinamopine)"],
                    "Plant interaction": ["chvA-B (chromosomal virulence)", "exoC (exopolysaccharide)", "cel (cellulose)", "att (attachment)", "biofilm formation"]
                })
            ]
        },

        "🦠 Plant Viruses": {
            "Tobamoviruses": [
                ("Tobacco mosaic virus", "Tobacco mosaic virus", {
                    "Structural genes": ["MP (movement protein)", "CP (coat protein)", "Rep (replicase)", "helicase domain", "polymerase domain"],
                    "Movement": ["MP (30K movement protein)", "cell-to-cell movement", "plasmodesmata gating", "viral transport", "cytoskeleton interaction"],
                    "Replication": ["RdRp (RNA-dependent RNA polymerase)", "126K protein", "183K protein", "methyltransferase", "helicase"],
                    "Host interaction": ["elicitor recognition", "hypersensitive response", "systemic acquired resistance", "pathogenesis-related proteins"],
                    "Resistance genes": ["N gene", "Tm-1", "Tm-2", "Tm-22", "resistance-breaking strains"]
                }),
                ("Tomato mosaic virus", "Tomato mosaic virus", {
                    "Essential genes": ["CP (coat protein)", "MP (movement protein)", "Rep (replicase)", "126K protein", "183K protein"],
                    "Pathogenicity": ["virulence determinants", "host range", "symptom expression", "systemic movement"],
                    "Resistance breaking": ["resistance-breaking strains", "Tm-2 resistance", "point mutations", "strain classification"]
                })
            ],
            "Other viruses": [
                ("Beet curly top virus", "Beet curly top virus", {
                    "Essential genes": ["C1 (replication initiator)", "C2 (transcription activator)", "C3 (replication enhancer)", "C4 (pathogenicity)", "V1 (coat protein)", "V2 (movement protein)"],
                    "Geminivirus features": ["circular ssDNA", "bipartite genome", "rolling circle replication", "transcription activation", "silencing suppression"],
                    "Host range": ["beet", "tomato", "pepper", "bean", "spinach", "squash", "cucumber"],
                    "Vector transmission": ["beet leafhopper", "Circulifer tenellus", "persistent transmission", "propagative transmission"],
                    "Pathogenicity factors": ["C4 (symptom determinant)", "C2 (transcription activator)", "V2 (movement protein)", "silencing suppression", "host range determination"]
                }),
                ("Arabis mosaic virus", "Arabis mosaic virus", {
                    "Essential genes": ["RNA1 (replication)", "RNA2 (movement)", "CP (coat protein)", "MP (movement protein)", "Rep (replicase)"],
                    "Nepovirus features": ["bipartite RNA genome", "polyprotein processing", "3' poly(A) tail", "5' VPg", "icosahedral particles"],
                    "Host range": ["arabis", "strawberry", "raspberry", "grapevine", "tobacco", "cucumber"],
                    "Vector transmission": ["nematode vectors", "Xiphinema diversicaudatum", "X. coxi", "soil transmission", "seed transmission"],
                    "Pathogenicity": ["systemic infection", "mosaic symptoms", "stunting", "yield reduction", "latent infection"]
                }),
                ("Alfalfa mosaic virus", "Alfalfa mosaic virus", {
                    "Essential genes": ["RNA1 (P1, P2)", "RNA2 (P3)", "RNA3 (MP, CP)", "P1 (replicase)", "P2 (helicase)", "P3 (polymerase)", "MP (movement protein)", "CP (coat protein)"],
                    "Alfamovirus features": ["tripartite RNA genome", "bacilliform particles", "coat protein requirement", "genome activation", "replication"],
                    "Host range": ["alfalfa", "tobacco", "tomato", "pepper", "bean", "cucumber", "lettuce"],
                    "Vector transmission": ["aphid transmission", "non-persistent", "stylet-borne", "Myzus persicae", "Aphis gossypii"],
                    "Pathogenicity": ["mosaic symptoms", "yellowing", "stunting", "systemic infection", "yield reduction"]
                }),
                ("Cannabis cryptic virus", "Cannabis cryptic virus", {
                    "Essential genes": ["RNA1 (RdRp)", "RNA2 (CP)", "RdRp (RNA-dependent RNA polymerase)", "CP (coat protein)", "MP (movement protein)"],
                    "Cryptic virus features": ["persistent infection", "latent infection", "no symptoms", "vertical transmission", "seed transmission"],
                    "Host specificity": ["Cannabis sativa", "hemp", "marijuana", "endophytic", "systemic infection"],
                    "Transmission": ["seed transmission", "pollen transmission", "no vector", "vertical transmission", "graft transmission"],
                    "Molecular features": ["dsRNA genome", "icosahedral particles", "persistent infection", "no cell-to-cell movement", "replication in cytoplasm"]
                }),
                ("Lettuce chlorosis virus", "Lettuce chlorosis virus", {
                    "Essential genes": ["RNA1 (replication)", "RNA2 (movement)", "P1 (replicase)", "P2 (helicase)", "P3 (polymerase)", "MP (movement protein)", "CP (coat protein)"],
                    "Crinivirus features": ["bipartite RNA genome", "whitefly transmission", "phloem-limited", "long flexuous particles", "genome activation"],
                    "Host range": ["lettuce", "tomato", "pepper", "cucumber", "melon", "squash", "bean"],
                    "Vector transmission": ["whitefly transmission", "Bemisia tabaci", "Trialeurodes vaporariorum", "semi-persistent", "circulative"],
                    "Pathogenicity": ["chlorosis", "yellowing", "stunting", "phloem necrosis", "yield reduction"]
                }),
                ("Potato virus Y", "Potato virus Y", {
                    "Essential genes": ["P1", "HC-Pro", "P3", "6K1", "CI", "6K2", "VPg", "NIa-Pro", "NIb", "CP"],
                    "Silencing suppression": ["HC-Pro (helper component proteinase)", "P1", "RNA silencing suppression", "siRNA binding", "RISC complex"],
                    "Strain differentiation": ["PVY-O (ordinary)", "PVY-N (necrotic)", "PVY-C", "PVY-Z", "PVYNTN", "PVYN-Wi"],
                    "Vector transmission": ["aphid transmission", "stylet-borne", "non-persistent", "helper component", "virion retention"],
                    "Recombination": ["recombinant strains", "PVYNTN", "breakpoints", "fitness advantage", "emergence"]
                }),
                ("Cucumber mosaic virus", "Cucumber mosaic virus", {
                    "Essential genes": ["1a", "2a", "2b", "3a", "MP", "CP"],
                    "Satellite RNA": ["satRNA", "symptom modulation", "D-satRNA", "Y-satRNA", "WL1-satRNA"],
                    "Subgroups": ["subgroup I", "subgroup II", "subgroup IA", "subgroup IB", "phylogenetic classification"],
                    "Host range": ["wide host range", "over 1000 species", "monocots", "dicots", "woody plants"],
                    "Symptom determinants": ["2b protein", "satellite RNA", "strain-specific", "host-dependent", "environmental factors"]
                })
            ]
        },

        "🦠 Oomycetes": {
            "Water molds": [
                ("Pythium root rot", "Pythium ultimum", {
                    "Essential genes": ["ACT1", "TUB1", "EF1A", "RPB1", "COX1"],
                    "Pathogenicity factors": ["PyCBEL1-10 (cellulose-binding elicitor lectin)", "PyPL1-8 (pectate lyases)", "PyCUT1-3 (cutinases)", "PyPRO1-5 (proteases)", "PyLIP1-3 (lipases)"],
                    "Zoospore motility": ["FLA1-20 (flagellar proteins)", "DYN1-5 (dynein)", "KIN1-8 (kinesin)", "TUB1-3 (tubulins)", "MOT1-10 (motility)"],
                    "Oospore formation": ["OOS1-15 (oospore genes)", "SEX1-5 (sexual reproduction)", "MAT1-3 (mating type)", "GER1-5 (germination)", "DOR1-3 (dormancy)"],
                    "Cell wall synthesis": ["CEL1-5 (cellulose synthase)", "CHI1-3 (chitinase)", "GEL1-3 (β-1,3-glucanase)", "CHS1 (chitin synthase)", "CAL1 (callose synthase)"]
                }),
                ("Pythium damping-off", "Pythium myriotylum", {
                    "Essential genes": ["ACT1", "TUB1", "EF1A", "RPB1", "COX1"],
                    "Pathogenicity factors": ["PmCBEL1-8 (cellulose-binding elicitor lectin)", "PmPL1-6 (pectate lyases)", "PmCUT1-2 (cutinases)", "PmPRO1-4 (proteases)", "PmLIP1-2 (lipases)"],
                    "Damping-off factors": ["DAMP1-5 (damping-off specific)", "ROOT1-3 (root infection)", "SEED1-2 (seed rot)", "STEM1-2 (stem rot)", "CROWN1 (crown rot)"],
                    "Host range": ["HOST1-10 (host specificity)", "CROP1-5 (crop-specific)", "SOL1-3 (solanaceous)", "CUC1-2 (cucurbit)", "LEG1-2 (legume)"],
                    "Environmental adaptation": ["TEMP1-3 (temperature)", "MOIST1-2 (moisture)", "pH1-2 (pH tolerance)", "SALT1 (salt tolerance)", "ANA1 (anaerobic)"]
                })
            ]
        },

        "🧬 Viroids": {
            "RNA pathogens": [
                ("Hop latent viroid", "Hop latent viroid", {
                    "Structural domains": ["P (pathogenicity)", "C (central)", "V (variable)", "T1 (terminal left)", "T2 (terminal right)"],
                    "Secondary structure": ["rod-like structure", "base-pairing", "hairpin loops", "bulges", "thermodynamic stability"],
                    "Variants": ["HLVd variants", "sequence diversity", "geographic strains", "hop cultivar adaptation"],
                    "Pathogenicity": ["latent infection", "hop stunt", "yield reduction", "brewing quality", "bitter compound"],
                    "Detection": ["RT-PCR", "northern blot", "dot blot", "in situ hybridization", "high-throughput sequencing"]
                })
            ]
        }
    }

def get_organism_suggestions():
    """Get agricultural pest and pathogen suggestions organized by category with comprehensive gene targets"""
    # Import the enhanced function
    return get_organism_suggestions_with_gene_targets()

def search_organism_with_gene_targets(organism_name, email, api_key=None, max_sequences=10):
    """Enhanced organism search that includes gene target information"""
    
    # Get gene targets for the organism
    suggestions = get_organism_suggestions_with_gene_targets()
    organism_targets = None
    
    # Normalize the input organism name
    organism_name_lower = organism_name.lower().strip()
    
    # Find matching organism and extract gene targets
    for category, subcategories in suggestions.items():
        for subcategory, organisms in subcategories.items():
            for item in organisms:
                if len(item) == 3:  # New format with gene targets
                    common_name, scientific_name, gene_targets = item
                else:  # Old format without gene targets
                    common_name, scientific_name = item
                    gene_targets = {"Essential genes": ["16S rRNA", "18S rRNA", "ACT1", "TUB1", "EF1A"]}
                
                # Improved matching logic
                scientific_lower = scientific_name.lower().strip()
                common_lower = common_name.lower().strip()
                
                # Check for exact matches or partial matches
                if (organism_name_lower == scientific_lower or 
                    organism_name_lower == common_lower or
                    organism_name_lower in scientific_lower or 
                    scientific_lower in organism_name_lower or
                    organism_name_lower in common_lower or
                    common_lower in organism_name_lower):
                    
                    organism_targets = {
                        'organism': scientific_name,
                        'common_name': common_name,
                        'category': category,
                        'subcategory': subcategory,
                        'gene_targets': gene_targets
                    }
                    break
            if organism_targets:
                break
        if organism_targets:
            break
    
    return organism_targets

def display_gene_targets_interface(organism_targets):
    """Display gene targets interface in Streamlit with improved state management"""
    
    if organism_targets:
        st.success(f"📋 **Gene Targets Available for {organism_targets['common_name']} ({organism_targets['organism']})**")
        
        with st.expander("🧬 View Available Gene Targets", expanded=True):
            st.write(f"**Category:** {organism_targets['category']} → {organism_targets['subcategory']}")
            
            # Display gene targets by category
            for gene_category, genes in organism_targets['gene_targets'].items():
                st.write(f"**{gene_category}:**")
                for gene in genes:
                    st.write(f"  • {gene}")
                st.write("")
            
            # Gene target selection
            st.subheader("🎯 Select Gene Targets for Primer Design")
            
            # Priority-based selection
            essential_cats = [cat for cat in organism_targets['gene_targets'].keys() if 'essential' in cat.lower() or 'housekeeping' in cat.lower()]
            pathogen_cats = [cat for cat in organism_targets['gene_targets'].keys() if any(x in cat.lower() for x in ['pathogen', 'virulence', 'effector'])]
            resistance_cats = [cat for cat in organism_targets['gene_targets'].keys() if 'resistance' in cat.lower()]
            
            # Initialize default categories only once per organism
            organism_key = f"gene_categories_{organism_targets['organism']}"
            if organism_key not in st.session_state:
                default_categories = []
                if essential_cats:
                    default_categories.extend(essential_cats[:2])  # First 2 essential categories
                if pathogen_cats:
                    default_categories.extend(pathogen_cats[:1])   # First pathogenicity category
                if resistance_cats:
                    default_categories.extend(resistance_cats[:1]) # First resistance category
                
                # Fallback if no defaults found
                if not default_categories:
                    default_categories = [list(organism_targets['gene_targets'].keys())[0]]
                
                st.session_state[organism_key] = default_categories
            
            # Use organism-specific stored selection
            selected_categories = st.multiselect(
                "Choose gene categories to target:",
                list(organism_targets['gene_targets'].keys()),
                default=st.session_state.get(organism_key, []),
                help="Select which gene categories to focus on for primer design. Essential genes are recommended for reliable detection.",
                key=f"gene_category_selector_{organism_targets['organism']}"
            )
            
            # Update organism-specific session state only if selection changed
            if selected_categories != st.session_state.get(organism_key, []):
                st.session_state[organism_key] = selected_categories
            
            if selected_categories:
                # Display selection summary
                total_genes = sum(len(organism_targets['gene_targets'][cat]) for cat in selected_categories)
                
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.metric("Selected Categories", len(selected_categories))
                with col2:
                    st.metric("Total Gene Targets", total_genes)
                with col3:
                    priority_score = sum(1 if get_gene_priority(cat) == "High" else 0.5 if get_gene_priority(cat) == "Medium" else 0.1 for cat in selected_categories)
                    st.metric("Priority Score", f"{priority_score:.1f}")
                
                # Show selected targets
                with st.expander("📋 Selected Gene Targets", expanded=False):
                    selected_genes = []
                    for category in selected_categories:
                        st.write(f"**{category}:** ({get_gene_priority(category)} Priority)")
                        genes = organism_targets['gene_targets'][category]
                        for gene in genes:
                            st.write(f"  • {gene}")
                            selected_genes.append(f"{category}: {gene}")
                        st.write(f"*{get_gene_use_recommendation(category)}*")
                        st.write("")
                
                # Store in session state - always update when selections change
                st.session_state.selected_gene_targets = {
                    'organism_info': organism_targets,
                    'selected_categories': selected_categories,
                    'selected_genes': selected_genes,
                    'total_targets': total_genes
                }
                
                return True
            else:
                # Clear session state when no categories are selected
                if 'selected_gene_targets' in st.session_state:
                    del st.session_state.selected_gene_targets
    
    return False

def display_results_with_gene_context():
    """Display primer results with gene target context - workflow aware"""
    
    # Check design mode
    design_mode = st.session_state.get('sequence_info', {}).get('design_mode', 'unknown')
    
    if design_mode == 'gene_targeted' and 'selected_gene_targets' in st.session_state:
        target_info = st.session_state.get('selected_gene_targets', {})
        
        st.subheader("🎯 Gene-Targeted Design Context")
        
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Target Organism", target_info['organism_info']['common_name'])
        with col2:
            st.metric("Gene Categories", len(target_info['selected_categories']))
        with col3:
            st.metric("Total Gene Targets", target_info['total_targets'])
        with col4:
            high_priority = sum(1 for cat in target_info['selected_categories'] if get_gene_priority(cat) == "High")
            st.metric("High Priority Categories", high_priority)
        
        # Target-specific recommendations
        with st.expander("📋 Gene Target Recommendations", expanded=False):
            recommendations = []
            for category in target_info['selected_categories']:
                rec = get_gene_use_recommendation(category)
                priority = get_gene_priority(category)
                gene_count = len(target_info['organism_info']['gene_targets'][category])
                recommendations.append({
                    'Category': category,
                    'Priority': priority,
                    'Gene Count': gene_count,
                    'Recommendation': rec[:100] + "..." if len(rec) > 100 else rec
                })
            
            rec_df = pd.DataFrame(recommendations)
            st.dataframe(rec_df, use_container_width=True)
    
    elif design_mode == 'conservation_based':
        st.subheader("🧬 Conservation-Based Design Context")
        analysis_metadata = st.session_state.get('analysis_metadata', {})
        
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
    
    elif design_mode == 'standard':
        st.subheader("⚡ Standard Design Context")
        st.info("Primers designed from single sequence using standard approach. Consider Gene-Targeted or Conservation-Based design for enhanced specificity.")

def export_with_gene_targets(primers, format_type="excel"):
    """Export primers with gene target information"""
    
    # Get organism and gene target info
    organism_info = st.session_state.get('sequence_info', {})
    gene_targets = st.session_state.get('selected_gene_targets', {})
    
    if format_type == "excel":
        return export_to_excel_with_targets(primers, gene_targets.get('organism_info'))
    else:
        # CSV export with gene target context
        data = []
        for i, primer in enumerate(primers):
            row_data = {
                'Primer_Pair': i + 1,
                'Target_Organism': gene_targets.get('organism_info', {}).get('organism', 'Unknown'),
                'Common_Name': gene_targets.get('organism_info', {}).get('common_name', 'Unknown'),
                'Organism_Category': gene_targets.get('organism_info', {}).get('category', 'Unknown'),
                'Selected_Gene_Categories': '; '.join(gene_targets.get('selected_categories', [])),
                'Total_Gene_Targets': gene_targets.get('total_targets', 0)
            }
            
            # Add primer-specific data
            if hasattr(primer, 'has_t7_promoter') and primer.has_t7_promoter:
                row_data.update({
                    'Forward_T7_Sequence': primer.forward_seq,
                    'Reverse_T7_Sequence': primer.reverse_seq,
                    'Forward_Core_Sequence': primer.core_forward_seq,
                    'Reverse_Core_Sequence': primer.core_reverse_seq,
                    'Core_Forward_Tm': round(primer.forward_tm, 2),
                    'Core_Reverse_Tm': round(primer.reverse_tm, 2),
                    'dsRNA_Size': primer.product_size,
                    'Primer_Type': 'T7_dsRNA'
                })
            else:
                row_data.update({
                    'Forward_Sequence': primer.forward_seq,
                    'Reverse_Sequence': primer.reverse_seq,
                    'Forward_Tm': round(primer.forward_tm, 2),
                    'Reverse_Tm': round(primer.reverse_tm, 2),
                    'Product_Size': primer.product_size,
                    'Primer_Type': 'Standard'
                })
            
            data.append(row_data)
        
        df = pd.DataFrame(data)
        return df.to_csv(index=False)

def get_gene_priority(category):
    """Get priority level for gene category"""
    high_priority = ["Essential genes", "Pathogenicity genes", "Resistance targets", "Insecticide resistance", "Acaricide resistance"]
    medium_priority = ["Secondary metabolite genes", "Detoxification genes", "Development genes", "Effector genes", "Type III secretion"]
    
    if category in high_priority:
        return "High"
    elif category in medium_priority:
        return "Medium"
    else:
        return "Low"

def get_gene_use_recommendation(category):
    """Get recommendation for gene category use"""
    recommendations = {
        "Essential genes": "Universal targets - always effective for species identification and basic primer design",
        "Pathogenicity genes": "Ideal for pathogen detection and virulence studies",
        "Resistance targets": "Critical for resistance monitoring and management strategies", 
        "Secondary metabolite genes": "Excellent for toxin detection and food safety applications",
        "Detoxification genes": "Essential for resistance monitoring in pest management",
        "Acaricide resistance": "Monitor for resistance development in mite populations",
        "Insecticide resistance": "Track resistance evolution in insect pest populations",
        "Development genes": "Target for growth regulation and lifecycle disruption",
        "Reproduction genes": "Fertility control and population management applications",
        "Effector genes": "Highly specific pathogenicity targets",
        "Cell wall degrading": "Virulence factors for disease mechanism studies",
        "Type III secretion": "Pathogenicity system targets",
        "Virulence regulation": "Regulatory control of disease development"
    }
    return recommendations.get(category, "General research and diagnostic applications")

def export_to_excel_with_targets(primers, organism_info):
    """Export primers to Excel with gene target information"""
    try:
        # Create main primer data
        primer_data = []
        for i, primer in enumerate(primers):
            row_data = {
                'Primer_Pair': i + 1,
                'Target_Organism': organism_info.get('organism', 'Unknown') if organism_info else 'Unknown',
                'Common_Name': organism_info.get('common_name', 'Unknown') if organism_info else 'Unknown',
                'Organism_Category': organism_info.get('category', 'Unknown') if organism_info else 'Unknown',
                'Organism_Subcategory': organism_info.get('subcategory', 'Unknown') if organism_info else 'Unknown'
            }
            
            # Add primer-specific data
            if hasattr(primer, 'has_t7_promoter') and primer.has_t7_promoter:
                row_data.update({
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
                row_data.update({
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
            
            primer_data.append(row_data)
        
        # Create gene target data if available
        gene_target_data = []
        if organism_info and 'gene_targets' in organism_info:
            targets = organism_info['gene_targets']
            for category, genes in targets.items():
                for gene in genes:
                    gene_target_data.append({
                        'Gene_Category': category,
                        'Gene_Target': gene,
                        'Target_Organism': organism_info.get('organism', 'Unknown'),
                        'Priority': get_gene_priority(category),
                        'Recommended_Use': get_gene_use_recommendation(category)
                    })
        
        # Create Excel file
        output = io.BytesIO()
        with pd.ExcelWriter(output, engine='openpyxl') as writer:
            # Primer results sheet
            primer_df = pd.DataFrame(primer_data)
            primer_df.to_excel(writer, sheet_name='Primer_Results', index=False)
            
            # Gene targets sheet
            if gene_target_data:
                targets_df = pd.DataFrame(gene_target_data)
                targets_df.to_excel(writer, sheet_name='Gene_Targets', index=False)
        
        return output.getvalue()
    except Exception as e:
        st.error(f"Error exporting to Excel: {e}")
        return b""

def create_comprehensive_gene_target_export():
    """Create comprehensive gene target database export"""
    try:
        suggestions = get_organism_suggestions_with_gene_targets()
        data = []
        
        for category, subcategories in suggestions.items():
            for subcategory, organisms in subcategories.items():
                for item in organisms:
                    if len(item) == 3:
                        common_name, scientific_name, gene_targets = item
                    else:
                        common_name, scientific_name = item
                        gene_targets = {"Essential genes": ["16S rRNA", "18S rRNA", "ACT1", "TUB1", "EF1A"]}
                    
                    for gene_category, genes in gene_targets.items():
                        for gene in genes:
                            data.append({
                                'Organism_Category': category,
                                'Organism_Subcategory': subcategory,
                                'Common_Name': common_name,
                                'Scientific_Name': scientific_name,
                                'Gene_Category': gene_category,
                                'Gene_Name': gene,
                                'Priority': get_gene_priority(gene_category),
                                'Recommendation': get_gene_use_recommendation(gene_category),
                                'Target_Type': 'Pathogen' if any(x in category for x in ['🍄', '🦠']) else 'Pest'
                            })
        
        return pd.DataFrame(data)
    except Exception as e:
        st.error(f"Error creating gene target export: {e}")
        return pd.DataFrame()

def generate_gene_target_statistics():
    """Generate comprehensive statistics about available gene targets"""
    try:
        df = create_comprehensive_gene_target_export()
        
        if df.empty:
            return {}, {}, {}
        
        stats = {
            'Total Organisms': df['Scientific_Name'].nunique(),
            'Total Gene Targets': len(df),
            'Organism Categories': df['Organism_Category'].nunique(),
            'Gene Categories': df['Gene_Category'].nunique(),
            'High Priority Targets': len(df[df['Priority'] == 'High']),
            'Pathogen Targets': len(df[df['Target_Type'] == 'Pathogen']),
            'Pest Targets': len(df[df['Target_Type'] == 'Pest'])
        }
        
        category_counts = df['Organism_Category'].value_counts().to_dict()
        gene_category_counts = df['Gene_Category'].value_counts().to_dict()
        
        return stats, category_counts, gene_category_counts
    except Exception as e:
        st.error(f"Error generating statistics: {e}")
        return {}, {}, {}

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
    gene_target: str = "Standard Design"  # Specific gene target for this primer pair

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
                      custom_params: Optional[Dict] = None, add_t7_promoter: bool = False, 
                      gene_target: str = "Standard Design") -> List[PrimerPair]:
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
                        penalty=penalty,
                        gene_target=gene_target
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
                        penalty=penalty,
                        gene_target=gene_target
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
    
    def fetch_complete_genomes(self, organism_name, max_sequences=10):
        """Fetch complete genomes/sequences for conservation analysis"""
        try:
            # Search for complete genomes and complete sequences only
            complete_genome_query = f'"{organism_name}"[organism] AND (complete genome[title] OR complete sequence[title] OR whole genome[title] OR chromosome[title] OR genome[title])'
            
            seq_ids = self.ncbi.search_sequences(
                complete_genome_query, 
                database="nucleotide", 
                max_results=max_sequences * 2  # Get more results to filter
            )
            
            sequences = []
            for seq_id in seq_ids:
                seq_info = self.ncbi.fetch_sequence_info(seq_id)
                if seq_info and seq_info.get('sequence'):
                    description = seq_info.get('description', '').lower()
                    
                    # Filter for complete genomes/sequences only
                    is_complete = any(keyword in description for keyword in [
                        'complete genome', 'complete sequence', 'whole genome', 
                        'chromosome', 'genome sequence', 'complete'
                    ])
                    
                    # Also check sequence length - complete genomes are typically much longer
                    clean_seq = re.sub(r'[^ATGCatgc]', '', seq_info['sequence'].upper())
                    is_long_enough = len(clean_seq) >= 100000  # At least 100kb for complete genomes
                    
                    if is_complete and is_long_enough:
                        sequences.append({
                            'id': seq_id,
                            'description': seq_info.get('description', ''),
                            'organism': seq_info.get('organism', ''),
                            'length': len(clean_seq),
                            'sequence': clean_seq[:50000],  # Use larger sample for conservation analysis
                            'full_length': len(clean_seq),
                            'type': 'complete_genome'
                        })
                        
                        if len(sequences) >= max_sequences:
                            break
            
            # If we didn't find enough complete genomes, try a broader search but still filter
            if len(sequences) < max_sequences:
                st.warning(f"Found only {len(sequences)} complete genomes. Searching for additional complete sequences...")
                
                # Broader search for complete sequences
                broader_query = f'"{organism_name}"[organism] AND complete[title]'
                additional_ids = self.ncbi.search_sequences(
                    broader_query, 
                    database="nucleotide", 
                    max_results=max_sequences
                )
                
                for seq_id in additional_ids:
                    if len(sequences) >= max_sequences:
                        break
                        
                    seq_info = self.ncbi.fetch_sequence_info(seq_id)
                    if seq_info and seq_info.get('sequence'):
                        description = seq_info.get('description', '').lower()
                        clean_seq = re.sub(r'[^ATGCatgc]', '', seq_info['sequence'].upper())
                        
                        # More lenient filtering for complete sequences
                        is_complete = 'complete' in description and len(clean_seq) >= 10000  # At least 10kb
                        
                        if is_complete and not any(seq['id'] == seq_id for seq in sequences):
                            sequences.append({
                                'id': seq_id,
                                'description': seq_info.get('description', ''),
                                'organism': seq_info.get('organism', ''),
                                'length': len(clean_seq),
                                'sequence': clean_seq[:50000],
                                'full_length': len(clean_seq),
                                'type': 'complete_sequence'
                            })
            
            return sequences
            
        except Exception as e:
            st.error(f"Error fetching complete genomes: {e}")
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
        'analysis_metadata': {},
        # GENE TARGET VARIABLES:
        'selected_gene_targets': {},
        'gene_target_stats': {},
        # ADDITIONAL PERSISTENCE VARIABLES:
        'session_initialized': True,
        'last_activity': None
    }
    
    for var, default_value in session_vars.items():
        if var not in st.session_state:
            st.session_state[var] = default_value
    
    # Ensure critical variables are never None
    if st.session_state.get('primers_designed') is None:
        st.session_state.primers_designed = []
    if st.session_state.get('sequence_info') is None:
        st.session_state.sequence_info = {}
    if st.session_state.get('specificity_results') is None:
        st.session_state.specificity_results = {}

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
            # Handle existing primer pairs that don't have gene_target attribute
            gene_target = getattr(primer, 'gene_target', 'Standard Design')
            
            if hasattr(primer, 'has_t7_promoter') and primer.has_t7_promoter:
                data.append({
                    'Primer_Pair': i + 1,
                    'Gene_Target': gene_target,
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
                    'Gene_Target': gene_target,
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



def search_organism_with_gene_targets(organism_name, email, api_key=None):
    """Search for organism and return gene targets if available"""
    try:
        # Get gene targets for the organism
        suggestions = get_organism_suggestions_with_gene_targets()
        organism_targets = None
        
        # Normalize the input organism name
        organism_name_lower = organism_name.lower().strip()
        
        # Debug: Print what we're searching for
        print(f"DEBUG: Searching for organism: '{organism_name}' (normalized: '{organism_name_lower}')")
        
        # Find matching organism and extract gene targets
        for category, subcategories in suggestions.items():
            for subcategory, organisms in subcategories.items():
                for item in organisms:
                    if len(item) == 3:  # New format with gene targets
                        common_name, scientific_name, gene_targets = item
                    else:  # Old format without gene targets
                        common_name, scientific_name = item
                        gene_targets = {"Essential genes": ["16S rRNA", "18S rRNA", "ACT1", "TUB1", "EF1A"]}
                    
                    # Improved matching logic - more strict
                    scientific_lower = scientific_name.lower().strip()
                    common_lower = common_name.lower().strip()
                    
                    # Check for exact matches first, then more flexible partial matches
                    if (organism_name_lower == scientific_lower or 
                        organism_name_lower == common_lower or
                        # Allow partial matches for genus names (e.g., "fusarium" matches "Fusarium oxysporum")
                        (organism_name_lower in scientific_lower and len(organism_name_lower) > 2) or
                        (organism_name_lower in common_lower and len(organism_name_lower) > 2) or
                        # Allow reverse partial matches (e.g., "oxysporum" matches "Fusarium oxysporum")
                        (scientific_lower in organism_name_lower and len(scientific_lower) > 2) or
                        (common_lower in organism_name_lower and len(common_lower) > 2)):
                        
                        print(f"DEBUG: Found match - Scientific: '{scientific_name}', Common: '{common_name}'")
                        
                        organism_targets = {
                            'organism': scientific_name,
                            'common_name': common_name,
                            'category': category,
                            'subcategory': subcategory,
                            'gene_targets': gene_targets
                        }
                        break
                if organism_targets:
                    break
            if organism_targets:
                break
        
        if not organism_targets:
            print(f"DEBUG: No match found for '{organism_name}'")
            # Try to find similar organisms for suggestions
            similar_organisms = []
            for category, subcategories in suggestions.items():
                for subcategory, organisms in subcategories.items():
                    for item in organisms:
                        if len(item) == 3:
                            common_name, scientific_name, gene_targets = item
                            # Check if any word from the input matches any word from the organism names
                            input_words = set(organism_name_lower.split())
                            scientific_words = set(scientific_name.lower().split())
                            common_words = set(common_name.lower().split())
                            
                            if input_words.intersection(scientific_words) or input_words.intersection(common_words):
                                similar_organisms.append(f"{common_name} ({scientific_name})")
            
            if similar_organisms:
                print(f"DEBUG: Found similar organisms: {similar_organisms[:5]}")
        
        return organism_targets
        
    except Exception as e:
        print(f"Error in search_organism_with_gene_targets: {e}")
        return None

def perform_gene_targeted_design(organism_name, email, api_key, max_sequences, custom_params, enable_t7_dsrna, optimal_dsrna_length, check_transcription_efficiency):
    """Perform gene-targeted primer design workflow with better error handling"""
    
    # Debug: Print what organism we're designing primers for
    print(f"DEBUG: perform_gene_targeted_design called with organism: '{organism_name}'")
    
    with st.spinner(f"Designing gene-targeted primers for {organism_name}..."):
        try:
            # Check if gene targets have been selected
            if 'selected_gene_targets' not in st.session_state:
                st.error("Please select gene targets first.")
                return
            
            gene_targets = st.session_state.selected_gene_targets
            selected_genes = gene_targets.get('selected_genes', [])
            
            if not selected_genes:
                st.error("No gene targets selected.")
                return
            
            # Initialize NCBI connector with timeout handling
            try:
                ncbi = NCBIConnector(email, api_key)
                designer = PrimerDesigner()
            except Exception as e:
                st.error(f"Failed to initialize NCBI connection: {e}")
                return
            
            st.write("🔍 **Step 1: Searching for gene-specific sequences...**")
            
            # Try to find sequences for selected genes
            gene_sequences = []
            for gene_info in selected_genes[:5]:  # Limit to first 5 genes
                try:
                    # Parse gene info: "Category: Gene name"
                    if ': ' in gene_info:
                        category, gene_name = gene_info.split(': ', 1)
                        # Clean gene name for search
                        clean_gene = gene_name.split('(')[0].strip()
                        
                        st.write(f"Searching for {clean_gene} in {organism_name}...")
                        
                        # Search for gene-specific sequences
                        search_query = f'"{organism_name}"[organism] AND "{clean_gene}"'
                        seq_ids = ncbi.search_sequences(search_query, database="nucleotide", max_results=2)
                        
                        if seq_ids:
                            for seq_id in seq_ids[:1]:  # Take first sequence per gene
                                sequence = ncbi.fetch_sequence(seq_id)
                                if sequence and len(sequence) > 100:
                                    gene_sequences.append({
                                        'gene': clean_gene,
                                        'category': category,
                                        'sequence': sequence,
                                        'id': seq_id
                                    })
                                    st.write(f"✅ Found {clean_gene} sequence: {len(sequence)} bp")
                                    break
                        else:
                            st.write(f"⚠️ No specific sequence found for {clean_gene}")
                except Exception as e:
                    st.write(f"⚠️ Error searching for {gene_info}: {e}")
                    continue
            
            if not gene_sequences:
                st.warning("No gene-specific sequences found. Falling back to general organism search...")
                # Fallback to general organism search
                try:
                    search_query = f'"{organism_name}"[organism]'
                    seq_ids = ncbi.search_sequences(search_query, database="nucleotide", max_results=min(max_sequences, 5))
                    
                    if not seq_ids:
                        st.error(f"No sequences found for {organism_name}")
                        return
                    
                    # Use the first sequence
                    sequence = ncbi.fetch_sequence(seq_ids[0])
                    if not sequence:
                        st.error("Failed to fetch sequence")
                        return
                    
                    # Clean and prepare sequence
                    clean_sequence = re.sub(r'[^ATGCatgc]', '', sequence.upper())
                    if len(clean_sequence) < 100:
                        st.error("Sequence too short for primer design")
                        return
                    
                    # Limit sequence size for performance
                    if len(clean_sequence) > 10000:
                        clean_sequence = clean_sequence[:10000]
                        st.info("Using first 10kb of sequence for primer design")
                    
                    st.success(f"Found general sequence: {len(clean_sequence)} bp")
                    sequence_source = "General organism sequence (fallback)"
                    
                except Exception as e:
                    st.error(f"Error fetching sequences: {e}")
                    return
            else:
                st.success(f"Found {len(gene_sequences)} gene-specific sequences!")
                
                # Display found genes
                gene_data = []
                for gene_seq in gene_sequences:
                    gene_data.append({
                        'Gene': gene_seq['gene'],
                        'Category': gene_seq['category'],
                        'Sequence ID': gene_seq['id'],
                        'Length': f"{len(gene_seq['sequence']):,} bp"
                    })
                
                gene_df = pd.DataFrame(gene_data)
                st.dataframe(gene_df, use_container_width=True)
                
                # Use the longest gene sequence for primer design
                best_gene = max(gene_sequences, key=lambda x: len(x['sequence']))
                clean_sequence = re.sub(r'[^ATGCatgc]', '', best_gene['sequence'].upper())
                if len(clean_sequence) > 10000:
                    clean_sequence = clean_sequence[:10000]
                sequence_source = f"Gene-specific sequence: {best_gene['gene']}"
                selected_gene_info = best_gene
            
            st.write("🧬 **Step 2: Designing primers...**")
            
            try:
                # Determine the gene target for primer design
                if 'selected_gene_info' in locals():
                    gene_target_name = f"{selected_gene_info['gene']} ({selected_gene_info['category']})"
                else:
                    gene_target_name = "Gene-Targeted Design (Fallback)"
                
                st.write(f"Designing primers for: {gene_target_name}")
                
                # Design primers
                primers = designer.design_primers(
                    clean_sequence, 
                    custom_params=custom_params,
                    add_t7_promoter=enable_t7_dsrna,
                    gene_target=gene_target_name
                )
                
                if not primers:
                    st.warning("No suitable primers found. Try adjusting parameters.")
                    return
                
                # Store results with gene target information
                st.session_state.current_sequence = clean_sequence
                
                # Determine sequence ID and description based on whether we found gene-specific sequences
                if 'selected_gene_info' in locals():
                    seq_id = selected_gene_info['id']
                    description = f"Gene-targeted design for {selected_gene_info['gene']} in {organism_name}"
                    gene_target_info = {
                        'gene_name': selected_gene_info['gene'],
                        'gene_category': selected_gene_info['category'],
                        'sequence_source': sequence_source
                    }
                else:
                    seq_id = seq_ids[0] if 'seq_ids' in locals() else "unknown"
                    description = f"Gene-targeted design for {organism_name} (fallback)"
                    gene_target_info = {
                        'gene_name': 'General organism sequence',
                        'gene_category': 'Fallback',
                        'sequence_source': sequence_source
                    }
                
                st.session_state.sequence_info = {
                    "id": seq_id,
                    "description": description,
                    "length": len(clean_sequence),
                    "organism": organism_name,
                    "design_mode": "gene_targeted",
                    "gene_target_info": gene_target_info
                }
                st.session_state.primers_designed = primers
                
                if enable_t7_dsrna:
                    st.session_state.t7_dsrna_enabled = True
                    st.session_state.t7_settings = {
                        'optimal_length': optimal_dsrna_length,
                        'check_efficiency': check_transcription_efficiency
                    }
                
                st.success(f"✅ Successfully designed {len(primers)} primer pairs for {gene_target_name}!")
                st.info("📊 Go to Results tab to view detailed analysis!")
                
                # Show preview with gene target information
                preview_data = []
                for i, primer in enumerate(primers[:3]):
                    preview_data.append({
                        'Pair': i + 1,
                        'Gene Target': gene_target_name,
                        'Forward': primer.forward_seq[:30] + '...' if len(primer.forward_seq) > 30 else primer.forward_seq,
                        'Reverse': primer.reverse_seq[:30] + '...' if len(primer.reverse_seq) > 30 else primer.reverse_seq,
                        'Product Size': f"{primer.product_size} bp"
                    })
                
                st.dataframe(pd.DataFrame(preview_data), use_container_width=True)
                
            except Exception as e:
                st.error(f"Error in primer design: {e}")
                return
                
        except Exception as e:
            st.error(f"Gene-targeted design error: {e}")

def perform_conservation_based_design(organism_name, email, api_key, max_sequences, conservation_threshold, window_size, enable_specificity_testing, specificity_threshold, comparison_organisms, custom_params, enable_t7_dsrna):
    """Perform conservation-based primer design workflow"""
    with st.spinner(f"Performing conservation analysis for {organism_name}..."):
        try:
            # Clear any existing gene targets to avoid confusion
            if 'selected_gene_targets' in st.session_state:
                del st.session_state['selected_gene_targets']
            
            # Initialize managers
            sequence_manager = SequenceManager(NCBIConnector(email, api_key))
            analyzer = ConservationAnalyzer(NCBIConnector(email, api_key))
            
            # Step 1: Fetch complete genomes/sequences for conservation analysis
            st.write("🔍 **Step 1: Fetching complete genomes/sequences...**")
            st.info("🧬 **Conservation Analysis**: Only complete genomes and complete sequences are used for accurate conservation analysis. Partial sequences and specific genes are excluded.")
            sequences = sequence_manager.fetch_complete_genomes(organism_name, max_sequences)
            
            if not sequences:
                st.warning(f"No sequences found for {organism_name}")
                return
            
            st.success(f"Found {len(sequences)} sequences!")
            
            # Display sequence summary
            sequence_data = []
            for i, seq in enumerate(sequences):
                sequence_type = seq.get('type', 'unknown')
                type_emoji = "🧬" if sequence_type == 'complete_genome' else "📄" if sequence_type == 'complete_sequence' else "❓"
                
                sequence_data.append({
                    'ID': seq['id'],
                    'Type': f"{type_emoji} {sequence_type.replace('_', ' ').title()}",
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
                    comp_orgs = get_related_organisms(organism_name, max_organisms=100)
                
                st.write(f"Testing against {len(comp_orgs)} organisms for comprehensive specificity analysis...")
                
                # Add progress bar for specificity testing
                progress_bar = st.progress(0)
                st.write("Running specificity analysis...")
                
                specificity_results = analyzer.test_specificity(
                    consensus_seq,
                    comp_orgs,
                    max_similarity=specificity_threshold
                )
                
                progress_bar.progress(1.0)
                
                # Display specificity results
                specificity_data = []
                for org, result in specificity_results.items():
                    if 'error' not in result and isinstance(result, dict):
                        # Safely access result keys with defaults
                        max_similarity = result.get('max_similarity', 0.0)
                        is_specific = result.get('is_specific', False)
                        sequences_tested = result.get('sequences_tested', 0)
                        
                        specificity_data.append({
                            'Organism': org,
                            'Max Similarity': f"{max_similarity:.1%}",
                            'Specific': '✅' if is_specific else '❌',
                            'Sequences Tested': sequences_tested
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
                "region_info": best_region,
                "design_mode": "conservation_based"
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
                add_t7_promoter=enable_t7_dsrna,
                gene_target="Conservation-Based Design"
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
            st.error(f"Conservation analysis error: {e}")

def perform_standard_design(organism_name, email, api_key, max_sequences, custom_params, enable_t7_dsrna, optimal_dsrna_length, check_transcription_efficiency):
    """Perform standard primer design workflow"""
    with st.spinner(f"Searching for {organism_name} genomes..."):
        try:
            # Clear any existing gene targets and conservation data
            for key in ['selected_gene_targets', 'conserved_regions', 'specificity_results', 'analysis_metadata']:
                if key in st.session_state:
                    del st.session_state[key]
            
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
                            "organism": organism_name,
                            "design_mode": "standard"
                        }
                        
                        st.write("Designing primers...")
                        primers = designer.design_primers(
                            clean_sequence, 
                            custom_params=custom_params,
                            add_t7_promoter=enable_t7_dsrna,
                            gene_target="Standard Design"
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
                            
                            # Step 2: Specificity testing for standard primers
                            st.write("🎯 **Step 2: Testing primer specificity...**")
                            
                            # Get related organisms for specificity testing (100+ organisms)
                            related_organisms = get_related_organisms(organism_name, max_organisms=100)
                            st.write(f"Testing against {len(related_organisms)} related organisms for comprehensive specificity analysis...")
                            
                            # Test specificity for the best primers
                            specificity_results = {}
                            analyzer = ConservationAnalyzer(NCBIConnector(email, api_key))
                            
                            # Add progress bar for specificity testing
                            progress_bar = st.progress(0)
                            total_tests = 3 * len(related_organisms)  # Testing top 3 primers
                            current_test = 0
                            
                            # Test the top 3 primers (best penalty scores)
                            best_primers = sorted(primers, key=lambda p: p.penalty)[:3]
                            
                            for i, primer in enumerate(best_primers):
                                # Create a test sequence from the primer binding region
                                test_sequence = clean_sequence[primer.forward_start:primer.reverse_start + 1]
                                
                                if len(test_sequence) > 50:  # Ensure we have enough sequence for testing
                                    try:
                                        primer_specificity = analyzer.test_specificity(
                                            test_sequence,
                                            related_organisms,
                                            max_similarity=0.7  # 70% similarity threshold
                                        )
                                        specificity_results[f"Primer Pair {i+1}"] = primer_specificity
                                        
                                        # Update progress bar
                                        current_test += len(related_organisms)
                                        progress_bar.progress(min(current_test / total_tests, 1.0))
                                        
                                    except Exception as e:
                                        st.warning(f"Could not test specificity for Primer Pair {i+1}: {e}")
                                        current_test += len(related_organisms)
                                        progress_bar.progress(min(current_test / total_tests, 1.0))
                            
                            # Store specificity results
                            st.session_state.specificity_results = specificity_results
                            
                            # Display specificity results
                            if specificity_results:
                                st.subheader("Specificity Testing Results")
                                
                                specificity_data = []
                                for primer_name, results in specificity_results.items():
                                    for organism, result in results.items():
                                        if 'error' not in result:
                                            # Safely access result keys with defaults
                                            max_similarity = result.get('max_similarity', 0.0)
                                            is_specific = result.get('is_specific', False)
                                            sequences_tested = result.get('sequences_tested', 0)
                                            
                                            specificity_data.append({
                                                'Primer Pair': primer_name,
                                                'Test Organism': organism,
                                                'Max Similarity': f"{max_similarity:.1%}",
                                                'Specific': '✅ Yes' if is_specific else '❌ No',
                                                'Sequences Tested': sequences_tested
                                            })
                                
                                if specificity_data:
                                    specificity_df = pd.DataFrame(specificity_data)
                                    st.dataframe(specificity_df, use_container_width=True)
                                    
                                    # Summary statistics
                                    total_tests = len(specificity_data)
                                    specific_tests = sum(1 for row in specificity_data if row['Specific'] == '✅ Yes')
                                    specificity_percentage = (specific_tests / total_tests) * 100 if total_tests > 0 else 0
                                    
                                    if specificity_percentage >= 80:
                                        st.success(f"🎯 Excellent specificity: {specific_tests}/{total_tests} tests passed ({specificity_percentage:.0f}%)")
                                    elif specificity_percentage >= 60:
                                        st.info(f"🎯 Good specificity: {specific_tests}/{total_tests} tests passed ({specificity_percentage:.0f}%)")
                                    else:
                                        st.warning(f"⚠️ Moderate specificity: {specific_tests}/{total_tests} tests passed ({specificity_percentage:.0f}%)")
                            
                            preview_data = []
                            for i, primer in enumerate(primers[:5]):
                                preview_data.append({
                                    'Pair': i + 1,
                                    'Forward': primer.forward_seq[:30] + '...' if len(primer.forward_seq) > 30 else primer.forward_seq,
                                    'Reverse': primer.reverse_seq[:30] + '...' if len(primer.reverse_seq) > 30 else primer.reverse_seq,
                                    'Product Size': f"{primer.product_size} bp"
                                })
                            
                            st.dataframe(pd.DataFrame(preview_data), use_container_width=True)
                            st.info("📊 Go to other tabs to view detailed analysis with specificity results!")
                        else:
                            st.warning("No suitable primers found. Try adjusting parameters.")
                else:
                    st.error("Failed to fetch sequence")
            else:
                st.warning(f"No sequences found for {organism_name}")
                
        except Exception as e:
            st.error(f"Error: {e}")

def main():
    """Main Streamlit application"""
    
    # Initialize session state first
    init_session_state()
    
    # Ensure ALL critical session state variables exist with safe defaults
    required_vars = {
        't7_dsrna_enabled': False,
        'primers_designed': [],
        'sequence_info': {},
        'specificity_results': {},
        'selected_gene_targets': {},
        'analysis_metadata': {},
        'found_sequences': [],
        'target_organism': '',
        'conserved_regions': [],
        'conservation_sequences': [],
        'gene_target_stats': {},
        'search_results': None,
        'database_used': None,
        'comprehensive_analysis_results': None,
        't7_results': None,
        't7_settings': {},
        'current_sequence': '',
        'session_initialized': True,
        'last_activity': None,
        # ORGANISM SELECTION VARIABLES:
        'organism_name_input': '',
        'trigger_gene_search': False,
        'current_organism_targets': None,
        'selected_organism_name': '',
        'gene_targets_loaded': False,
        'stored_organism_name': '',
        'processing_gene_design': False,
        'selected_gene_categories': []
    }
    
    for var, default_value in required_vars.items():
        if var not in st.session_state:
            st.session_state[var] = default_value
    
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
    
    # Initialize T7 dsRNA parameters with default values
    optimal_dsrna_length = False
    check_transcription_efficiency = False
    
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
                # Check if an organism was selected from buttons or preserve current value
                default_organism = st.session_state.get('selected_organism_name', '')
                if not default_organism and st.session_state.get('stored_organism_name'):
                    # If no new selection but we have a stored organism, use that
                    default_organism = st.session_state.stored_organism_name
                
                organism_name = st.text_input("Enter organism name:", 
                                            value=default_organism,
                                            placeholder="e.g., Fusarium oxysporum, Coronavirus, Tetranychus urticae",
                                            key="organism_input")
                
                # Debug: Show current organism name state
                st.write(f"🔍 Debug: organism_name = '{organism_name}'")
                st.write(f"🔍 Debug: selected_organism_name = '{st.session_state.get('selected_organism_name', 'None')}'")
                st.write(f"🔍 Debug: stored_organism_name = '{st.session_state.get('stored_organism_name', 'None')}'")
                
                # Store the current organism name in session state for persistence
                if organism_name:
                    st.session_state.stored_organism_name = organism_name
                
                # Clear the selected organism after setting it, but preserve the input
                if 'selected_organism_name' in st.session_state:
                    del st.session_state['selected_organism_name']
            
            with col2:
                max_sequences = st.number_input("Max sequences to search", min_value=5, max_value=50, value=10)
            
            # ==========================================
            # NEW: WORKFLOW SELECTION
            # ==========================================
            st.subheader("🔬 Primer Design Strategy")
            
            workflow_choice = st.radio(
                "Choose your primer design approach:",
                [
                    "🎯 Gene-Targeted Design (Recommended for specific genes)",
                    "🧬 Conservation-Based Design (Recommended for robust primers)",
                    "⚡ Standard Design (Single sequence, fastest)"
                ],
                help="Select based on your research goals: specific gene detection vs. broad applicability vs. speed"
            )
            
            # Initialize workflow variables
            gene_targets_workflow = "🎯 Gene-Targeted Design" in workflow_choice
            conservation_workflow = "🧬 Conservation-Based Design" in workflow_choice
            standard_workflow = "⚡ Standard Design" in workflow_choice
            
            # ==========================================
            # WORKFLOW 1: GENE-TARGETED DESIGN (FIXED)
            # ==========================================
            if gene_targets_workflow:
                st.info("🎯 **Gene-Targeted Design Mode**\nDesign primers for specific genes with known biological functions.")
                
                # Debug: Show available organisms
                with st.expander("🔍 Debug: Available Organisms", expanded=False):
                    suggestions = get_organism_suggestions_with_gene_targets()
                    st.write("**Available organisms with gene targets:**")
                    for category, subcategories in suggestions.items():
                        st.write(f"**{category}:**")
                        for subcategory, organisms in subcategories.items():
                            st.write(f"  {subcategory}:")
                            for item in organisms[:3]:  # Show first 3 organisms per subcategory
                                if len(item) == 3:
                                    common_name, scientific_name, gene_targets = item
                                    st.write(f"    - {common_name} ({scientific_name})")
                                else:
                                    st.write(f"    - {item}")
                
                if organism_name and email:
                    # Simple gene target loading
                    if st.button("🔍 Load Gene Targets", type="secondary"):
                        st.write(f"🔍 Debug: Button clicked for organism: '{organism_name}'")
                        st.write(f"🔍 Debug: organism_name type: {type(organism_name)}")
                        st.write(f"🔍 Debug: organism_name length: {len(organism_name) if organism_name else 0}")
                        with st.spinner(f"Loading gene targets for {organism_name}..."):
                            organism_targets = search_organism_with_gene_targets(organism_name, email, api_key)
                            
                            st.write(f"🔍 Debug: search_organism_with_gene_targets returned: {organism_targets}")
                            
                            if organism_targets:
                                st.session_state.current_organism_targets = organism_targets
                                st.success(f"Found gene targets for {organism_targets['common_name']}")
                            else:
                                st.warning("No specific gene targets found for this organism.")
                                st.write("🔍 Debug: No organism targets found. Check the console for debug output.")
                                
                                # Show some suggestions
                                st.info("💡 **Try these examples:**")
                                st.write("- Fusarium oxysporum")
                                st.write("- Fusarium graminearum") 
                                st.write("- Aspergillus niger")
                                st.write("- Botrytis cinerea")
                                st.write("- Or click on one of the organism buttons above")
                
                    # Display gene targets if available
                    if 'current_organism_targets' in st.session_state and st.session_state.current_organism_targets:
                        organism_targets = st.session_state.current_organism_targets
                        
                        # Show available gene categories
                        st.subheader("Available Gene Categories")
                        gene_categories = list(organism_targets['gene_targets'].keys())
                        
                        selected_categories = st.multiselect(
                            "Select gene categories:",
                            gene_categories,
                            default=gene_categories[:2] if len(gene_categories) >= 2 else gene_categories,
                            key="gene_category_select"
                        )
                        
                        if selected_categories:
                            # Store selection in session state
                            st.session_state.selected_gene_categories = selected_categories
                            
                            # Show selected genes
                            st.write("**Selected Gene Targets:**")
                            for category in selected_categories:
                                st.write(f"**{category}:**")
                                for gene in organism_targets['gene_targets'][category]:
                                    st.write(f"  • {gene}")
                            
                            # Design button
                            if st.button("🎯 Design Gene-Targeted Primers", type="primary", key="gene_design_final"):
                                # Set up gene target info for design
                                st.session_state.selected_gene_targets = {
                                    'organism_info': organism_targets,
                                    'selected_categories': selected_categories,
                                    'selected_genes': [f"{cat}: {gene}" for cat in selected_categories 
                                                     for gene in organism_targets['gene_targets'][cat]]
                                }
                                
                                # Debug: Print what organism name we're passing
                                print(f"DEBUG: Calling perform_gene_targeted_design with organism_name: '{organism_name}'")
                                
                                # Call the design function
                                perform_gene_targeted_design(
                                    organism_name, email, api_key, max_sequences, 
                                    custom_params, enable_t7_dsrna, 
                                    optimal_dsrna_length, check_transcription_efficiency
                                )
                        else:
                            st.info("Select at least one gene category to proceed.")
                else:
                    if not email:
                        st.error("Email required for gene target search")
                    if not organism_name:
                        st.info("Enter an organism name above to search for gene targets.")
            
            # ==========================================
            # WORKFLOW 2: CONSERVATION-BASED DESIGN  
            # ==========================================
            elif conservation_workflow:
                st.info("🧬 **Conservation-Based Design Mode**\nAnalyze multiple complete genomes/sequences to find conserved regions for robust, broad-spectrum primers. Ideal for species identification and population studies.")
                st.warning("⚠️ **Important**: Conservation analysis requires complete genomes or complete sequences only. Partial sequences and specific genes are excluded to ensure accurate conservation patterns.")
                
                with st.expander("Conservation Analysis Parameters", expanded=True):
                    col1, col2 = st.columns(2)
                    with col1:
                        conservation_threshold = st.slider(
                            "Conservation threshold (%)",
                            min_value=0, max_value=95, value=80,
                            help="Minimum conservation required across sequences (0% = no conservation required)"
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
                            help="Test primer specificity against 100+ related organisms across diverse taxonomic groups for comprehensive analysis"
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
                
                # Search button for conservation-based design
                if st.button("🧬 Design Conservation-Based Primers", type="primary", use_container_width=True):
                    if not email or not organism_name:
                        st.error("❌ Please provide email and organism name.")
                    else:
                        perform_conservation_based_design(
                            organism_name, email, api_key, max_sequences, 
                            conservation_threshold, window_size, enable_specificity_testing, 
                            specificity_threshold if enable_specificity_testing else None,
                            comparison_organisms, custom_params, enable_t7_dsrna
                        )
            
            # ==========================================
            # WORKFLOW 3: STANDARD DESIGN
            # ==========================================
            else:  # standard_workflow
                st.info("⚡ **Standard Design Mode**\nQuick primer design from the first available sequence. Ideal for rapid prototyping and basic applications.")
                
                # Search button for standard design
                if st.button("⚡ Design Standard Primers", type="primary", use_container_width=True):
                    if not email or not organism_name:
                        st.error("❌ Please provide email and organism name.")
                    else:
                        perform_standard_design(organism_name, email, api_key, max_sequences, custom_params, enable_t7_dsrna, optimal_dsrna_length, check_transcription_efficiency)
            
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
                        
                        for i, organism_item in enumerate(organisms):
                            # Handle both old format (common_name, latin_name) and new format (common_name, latin_name, gene_targets)
                            if len(organism_item) == 3:
                                common_name, latin_name, gene_targets = organism_item
                            else:
                                common_name, latin_name = organism_item
                            
                            with cols[i % num_cols]:
                                # Create a unique key and use callback to set organism name directly
                                button_key = f"suggest_{category}_{subcategory}_{i}_{latin_name.replace(' ', '_')}"
                                
                                if st.button(
                                    common_name, 
                                    key=button_key, 
                                    help=f"Search for {latin_name}",
                                    use_container_width=True
                                ):
                                    # Set the organism name to appear in the text input
                                    st.session_state.selected_organism_name = latin_name
                                    st.rerun()
                        
                        # Add small spacing between subcategories
                        if subcategory != list(subcategories.keys())[-1]:  # Not the last subcategory
                            st.markdown("")
            
        
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
                                add_t7_promoter=enable_t7_dsrna,
                                gene_target="User Input Sequence"
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
                                
                                # Note about specificity testing for direct sequence input
                                st.info("ℹ️ **Note**: Specificity testing is not available for direct sequence input since no organism information is provided. For specificity testing, please use the 'Organism Name' input method.")
                                
                                st.info("📊 Go to other tabs to view detailed analysis!")
                            else:
                                st.warning("No suitable primers found with current parameters")
                        except Exception as e:
                            st.error(f"Error: {e}")
    
    with tab2:
        st.header("Primer Design Results")
        
        # Debug information
        with st.expander("🔍 Debug Information", expanded=True):
            st.write("**Session State Check:**")
            st.write(f"- Has primers: {state_check['has_primers']}")
            st.write(f"- Has sequence: {state_check['has_sequence']}")
            st.write(f"- Has sequence info: {state_check['has_seq_info']}")
            st.write(f"- Primer count: {state_check['primer_count']}")
            st.write(f"- Sequence length: {state_check['sequence_length']}")
            
            if 'primers_designed' in st.session_state:
                st.write(f"- Primers in session state: {len(st.session_state['primers_designed'])}")
                if st.session_state['primers_designed']:
                    st.write(f"- First primer: {st.session_state['primers_designed'][0]}")
            else:
                st.write("- No primers in session state")
            
            st.write("**All session state keys:**")
            st.write(list(st.session_state.keys()))
        
        if not state_check['has_primers']:
            st.info("No primers designed yet. Please use the Input tab to design primers.")
            st.stop()
        
        primers = st.session_state.get('primers_designed', [])
        t7_enabled = st.session_state.get('t7_dsrna_enabled', False)
        
        # Display gene target context if available
        try:
            display_results_with_gene_context()
        except Exception as e:
            st.warning(f"Gene target context display unavailable: {e}")
            # Continue with basic results display even if gene context fails
        
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
                if 'error' not in result and isinstance(result, dict):
                    # Safely access result keys with defaults
                    max_similarity = result.get('max_similarity', 0.0)
                    is_specific = result.get('is_specific', False)
                    sequences_tested = result.get('sequences_tested', 0)
                    
                    specificity_data.append({
                        'Organism': organism,
                        'Max Similarity': f"{max_similarity:.1%}",
                        'Specific': '✅ Yes' if is_specific else '❌ No',
                        'Sequences Tested': sequences_tested
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
            # Handle existing primer pairs that don't have gene_target attribute
            gene_target = getattr(primer, 'gene_target', 'Standard Design')
            
            if hasattr(primer, 'has_t7_promoter') and primer.has_t7_promoter:
                row = {
                    'Pair': i + 1,
                    'Gene Target': gene_target,
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
                    'Gene Target': gene_target,
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
        
        # Debug information
        with st.expander("🔍 Analysis Debug Information", expanded=True):
            st.write("**Session State Check:**")
            st.write(f"- Has primers: {state_check['has_primers']}")
            st.write(f"- Primer count: {state_check['primer_count']}")
            st.write(f"- All session state keys: {list(st.session_state.keys())}")
        
        if not state_check['has_primers']:
            st.info("No primers designed yet. Please use the Input tab to design primers.")
            st.stop()
        
        primers = st.session_state.get('primers_designed', [])
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
        
        # Debug information
        with st.expander("🔍 Export Debug Information", expanded=True):
            st.write("**Session State Check:**")
            st.write(f"- Has primers: {state_check['has_primers']}")
            st.write(f"- Primer count: {state_check['primer_count']}")
            st.write(f"- All session state keys: {list(st.session_state.keys())}")
        
        if not state_check['has_primers']:
            st.info("No primers to export. Please design primers first.")
            st.stop()
        
        primers = st.session_state.get('primers_designed', [])
        t7_enabled = st.session_state.get('t7_dsrna_enabled', False)
        gene_targets_available = 'selected_gene_targets' in st.session_state
        
        if gene_targets_available:
            st.success("🧬 **Gene Target Information Available** - Export will include comprehensive gene target data")
        
        if t7_enabled:
            st.info("🧬 **T7 dsRNA Export Mode** - Export includes both full T7 primers and core sequences")
        
        st.subheader("Download Options")
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            if st.button("📊 Download Enhanced Excel", type="primary"):
                try:
                    excel_data = export_with_gene_targets(primers, "excel")
                    if excel_data:
                        filename = "enhanced_primer_results_with_targets.xlsx"
                        st.download_button(
                            label="Click to Download Enhanced Excel File",
                            data=excel_data,
                            file_name=filename,
                            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                        )
                except Exception as e:
                    st.error(f"Error creating enhanced Excel: {e}")
        
        with col2:
            try:
                csv_data = export_with_gene_targets(primers, "csv")
                filename = "primer_results_with_gene_targets.csv"
                st.download_button(
                    label="📄 Download Enhanced CSV",
                    data=csv_data,
                    file_name=filename,
                    mime="text/csv"
                )
            except Exception as e:
                st.error(f"Error creating enhanced CSV: {e}")
        
        with col3:
            if gene_targets_available:
                # Export comprehensive gene target database
                try:
                    all_targets_df = create_comprehensive_gene_target_export()
                    targets_csv = all_targets_df.to_csv(index=False)
                    st.download_button(
                        label="🗃️ Download Complete Gene Database",
                        data=targets_csv,
                        file_name="comprehensive_gene_targets_database.csv",
                        mime="text/csv"
                    )
                except:
                    st.info("Gene database export not available")
        
        # Gene Target Statistics
        if gene_targets_available:
            st.subheader("📊 Gene Target Statistics")
            
            try:
                stats, category_counts, gene_category_counts = generate_gene_target_statistics()
                
                col1, col2, col3, col4 = st.columns(4)
                with col1:
                    st.metric("Total Organisms", stats['Total Organisms'])
                with col2:
                    st.metric("Total Gene Targets", stats['Total Gene Targets'])
                with col3:
                    st.metric("High Priority Targets", stats['High Priority Targets'])
                with col4:
                    st.metric("Organism Categories", stats['Organism Categories'])
                
                # Category breakdown
                with st.expander("📈 Detailed Statistics", expanded=False):
                    col1, col2 = st.columns(2)
                    
                    with col1:
                        st.write("**Organism Categories:**")
                        for category, count in list(category_counts.items())[:10]:
                            st.write(f"• {category}: {count} organisms")
                    
                    with col2:
                        st.write("**Top Gene Categories:**")
                        for gene_cat, count in list(gene_category_counts.items())[:10]:
                            st.write(f"• {gene_cat}: {count} targets")
            except:
                st.info("Gene target statistics not available")
        
        # Original export functionality (fallback)
        st.subheader("📋 Standard Export Options")
        
        col1, col2 = st.columns(2)
        
        with col1:
            if st.button("📊 Download Standard Excel", type="secondary"):
                excel_data = export_to_excel(primers)
                if excel_data:
                    filename = "t7_dsrna_primers.xlsx" if t7_enabled else "primer_results.xlsx"
                    st.download_button(
                        label="Click to Download Standard Excel File",
                        data=excel_data,
                        file_name=filename,
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                    )
        
        with col2:
            try:
                data = []
                for i, primer in enumerate(primers):
                    # Handle existing primer pairs that don't have gene_target attribute
                    gene_target = getattr(primer, 'gene_target', 'Standard Design')
                    
                    if hasattr(primer, 'has_t7_promoter') and primer.has_t7_promoter:
                        row = {
                            'Primer_Pair': i + 1,
                            'Gene_Target': gene_target,
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
                            'Gene_Target': gene_target,
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
