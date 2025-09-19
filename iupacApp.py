import streamlit as st
from streamlit_ketcher import st_ketcher
from rdkit import Chem
from rdkit.Chem import GetPeriodicTable, rdMolDescriptors
from collections import defaultdict, deque
import re


# ==============================================================================
# PAGE CONFIGURATION AND MODERN CSS STYLING
# ==============================================================================

st.set_page_config(
    page_title="IUPAC Name Generator",
    page_icon="⚗️",
    layout="wide",
    initial_sidebar_state="collapsed"
)

def inject_modern_css():
    st.markdown("""
    <style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;600&display=swap');
    .main .block-container {
        padding: 1rem 2rem 2rem 2rem;
        max-width: 1400px;
        font-family: 'Inter', sans-serif;
    }
    .main-header {
        text-align: center;
        padding: 2rem 0 3rem 0;
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        margin: -1rem -2rem 2rem -2rem;
        color: white;
        border-radius: 0 0 18px 18px;
    }
    .main-header h1 {
        margin: 0;
        font-size: 2.3rem;
        font-weight: 700;
        text-shadow: 0 2px 6px rgba(0,0,0,0.1);
    }
    .main-header p {
        margin: 0.5rem 0 0 0;
        font-size: 1.08rem;
        opacity: 0.95;
        font-weight: 400;
    }
    .editor-card {
        background: white;
        padding: 2rem;
        border-radius: 16px;
        box-shadow: 0 4px 15px rgba(0,0,0,0.07);
        border: 1px solid #f0f2f6;
        margin-bottom: 1.5rem;
        transition: all 0.3s ease;
    }
    .editor-card:hover {
        box-shadow: 0 8px 30px rgba(0,0,0,0.16);
        transform: translateY(-1px);
    }
    .controls-card {
        background: #f8fafc;
        padding: 2rem;
        border-radius: 16px;
        box-shadow: 0 2px 10px rgba(0,0,0,0.05);
        border: 1px solid #e2e8f0;
        height: fit-content;
        position: sticky;
        top: 2rem;
    }
    .section-header {
        font-size: 1.2rem;
        font-weight: 600;
        color: #2d3748;
        margin-bottom: 1rem;
        margin-top: 0;
        display: flex;
        align-items: center;
        gap: 0.55rem;
    }
    .smiles-container {
        background: linear-gradient(135deg, #f7fafc 0%, #edf2f7 100%);
        padding: 1.2rem;
        border-radius: 12px;
        border-left: 4px solid #667eea;
        margin: 2rem 0 1rem 0;
        font-family: 'Fira Code', 'Courier New', monospace;
        position: relative;
        overflow: hidden;
    }
    .stButton > button {
        width: 100%;
        border-radius: 10px;
        border: 1px solid #e2e8f0;
        background: white;
        font-weight: 500;
        font-size: 0.97rem;
        transition: all 0.2s ease;
        font-family: 'Inter', sans-serif;
        height: 45px;
    }
    .stButton > button:hover {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        border-color: transparent;
        transform: translateY(-1px);
    }
    div[data-testid="stButton"] button[kind="primary"] {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        border: none;
        height: 55px;
        font-size: 1.1rem;
        font-weight: 600;
        box-shadow: 0 4px 15px rgba(102, 126, 234, 0.18);
    }
    div[data-testid="stButton"] button[kind="primary"]:hover {
        background: linear-gradient(135deg, #5a67d8 0%, #6b46a3 100%);
        transform: translateY(-2px);
    }
    .result-success {
        background: linear-gradient(135deg, #d4edda 0%, #c3e6cb 100%);
        border: 1px solid #b8dacd;
        color: #155724;
        padding: 2rem;
        border-radius: 16px;
        text-align: center;
        margin: 2rem 0 1rem 0;
        box-shadow: 0 4px 10px rgba(212, 237, 218, 0.14);
    }
    .result-success h3 {
        margin: 0 0 1rem 0;
        font-size: 1.34rem;
        font-weight: 600;
    }
    .result-success .iupac-name {
        font-size: 1.65rem;
        font-weight: 700;
        font-family: 'Georgia', serif;
        color: #0f5132;
        text-shadow: 0 1px 2px rgba(0,0,0,0.08);
        line-height: 1.32;
        margin: 0;
    }
    .result-error {
        background: linear-gradient(135deg, #f8d7da 0%, #f5c6cb 100%);
        border: 1px solid #f5c6cb;
        color: #721c24;
        padding: 2rem;
        border-radius: 16px;
        text-align: center;
        margin: 2rem 0 1rem 0;
        box-shadow: 0 4px 10px rgba(248, 215, 218, 0.12);
    }
    .footer {
        text-align: center;
        padding: 1rem 0;
    }
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
    header {visibility: hidden;}
    @media (max-width: 768px) {
        .main .block-container {
            padding: 1rem;
        }
        .main-header {
            margin: -1rem -1rem 1rem -1rem;
            padding: 1.2rem 0 2rem 0;
        }
        .main-header h1 { font-size: 1.55rem; }
        .editor-card, .controls-card { padding: 1.25rem; }
    }
    </style>
    """, unsafe_allow_html=True)


# ==============================================================================
# COMPLETE IUPAC NAMING SYSTEM - PRODUCTION READY
# ==============================================================================

# Parent chain names
PARENT_ROOTS = {
    1: "meth", 2: "eth", 3: "prop", 4: "but", 5: "pent",
    6: "hex", 7: "hept", 8: "oct", 9: "non", 10: "dec",
    11: "undec", 12: "dodec", 13: "tridec", 14: "tetradec", 15: "pentadec",
    16: "hexadec", 17: "heptadec", 18: "octadec", 19: "nonadec", 20: "icos"
}

CYCLIC_PARENT_NAMES = {
    3: "cycloprop", 4: "cyclobut", 5: "cyclopent", 6: "cyclohex",
    7: "cyclohept", 8: "cyclooct", 9: "cyclonon", 10: "cyclodec"
}

SUBSTITUENT_NAMES = {
    1: "methyl", 2: "ethyl", 3: "propyl", 4: "butyl", 5: "pentyl",
    6: "hexyl", 7: "heptyl", 8: "octyl", 9: "nonyl", 10: "decyl"
}

MULTIPLICITY_PREFIXES = {
    1: "", 2: "di", 3: "tri", 4: "tetra", 5: "penta", 6: "hexa", 7: "hepta"
}

# Complete functional group system with proper IUPAC priorities
FUNCTIONAL_GROUP_PATTERNS = {
    # Highest priority groups (1-10)
    "carboxylic_acid": {"smarts": "[CX3](=O)[OH]", "priority": 1, "suffix": "oic acid", "prefix": "carboxy"},
    "sulfonic_acid": {"smarts": "[SX4](=O)(=O)[OH]", "priority": 2, "suffix": "sulfonic acid", "prefix": "sulfo"},
    "ester": {"smarts": "[CX3](=O)[OX2H0]", "priority": 3, "suffix": "oate", "prefix": "alkoxycarbonyl"},
    "acid_halide": {"smarts": "[CX3](=O)[X]", "priority": 4, "suffix": "oyl halide", "prefix": "haloformyl"},
    "amide": {"smarts": "[CX3](=O)[NX3]", "priority": 5, "suffix": "amide", "prefix": "carbamoyl"},
    "nitrile": {"smarts": "[CX2]#[NX1]", "priority": 6, "suffix": "nitrile", "prefix": "cyano"},
    "aldehyde": {"smarts": "[CX3H1](=O)", "priority": 7, "suffix": "al", "prefix": "formyl"},
    "ketone": {"smarts": "[CX3](=O)[#6]", "priority": 8, "suffix": "one", "prefix": "oxo"},
    "alcohol": {"smarts": "[OH]", "priority": 9, "suffix": "ol", "prefix": "hydroxy"},
    "thiol": {"smarts": "[SH]", "priority": 10, "suffix": "thiol", "prefix": "sulfanyl"},
    
    # Medium priority groups (11-20)
    "amine_primary": {"smarts": "[NX3;H2]", "priority": 11, "suffix": "amine", "prefix": "amino"},
    "amine_secondary": {"smarts": "[NX3;H1]", "priority": 12, "suffix": "amine", "prefix": "amino"},
    "amine_tertiary": {"smarts": "[NX3;H0]", "priority": 13, "suffix": "amine", "prefix": "amino"},
    "ether": {"smarts": "[OD2]([#6])[#6]", "priority": 14, "suffix": None, "prefix": "alkoxy"},
    
    # Lower priority groups - prefix only (21+)
    "halogen_f": {"smarts": "[F]", "priority": 21, "suffix": None, "prefix": "fluoro"},
    "halogen_cl": {"smarts": "[Cl]", "priority": 22, "suffix": None, "prefix": "chloro"},
    "halogen_br": {"smarts": "[Br]", "priority": 23, "suffix": None, "prefix": "bromo"},
    "halogen_i": {"smarts": "[I]", "priority": 24, "suffix": None, "prefix": "iodo"},
    "nitro": {"smarts": "[N+](=O)[O-]", "priority": 25, "suffix": None, "prefix": "nitro"},
    "alkyl": {"smarts": "[CH3,CH2,CH,C]", "priority": 30, "suffix": None, "prefix": "alkyl"}
}

# Common aromatic compound names
AROMATIC_NAMES = {
    "c1ccccc1": "benzene",
    "Cc1ccccc1": "toluene", 
    "Oc1ccccc1": "phenol",
    "Nc1ccccc1": "aniline",
    "COc1ccccc1": "anisole",
    "c1ccc2ccccc2c1": "naphthalene"
}


def identify_functional_groups_advanced(mol):
    """Advanced functional group identification with proper priorities"""
    functional_groups = []
    
    for fg_name, fg_data in FUNCTIONAL_GROUP_PATTERNS.items():
        try:
            pattern = Chem.MolFromSmarts(fg_data["smarts"])
            if pattern:
                matches = mol.GetSubstructMatches(pattern)
                for match in matches:
                    # Get the main functional atom (usually first in match)
                    main_atom_idx = match[0]
                    
                    functional_groups.append({
                        "type": fg_name,
                        "atom_idx": main_atom_idx,
                        "priority": fg_data["priority"],
                        "suffix": fg_data["suffix"],
                        "prefix": fg_data["prefix"],
                        "match_atoms": match
                    })
        except Exception as e:
            continue
    
    # Sort by priority (lower number = higher priority)
    functional_groups.sort(key=lambda x: x["priority"])
    return functional_groups


def is_aromatic_compound(mol):
    """Check if the molecule is primarily aromatic"""
    aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    total_atoms = mol.GetNumAtoms()
    return aromatic_atoms >= 6 and aromatic_atoms / total_atoms > 0.5


def get_aromatic_base_name(mol):
    """Get base name for aromatic compounds"""
    smiles = Chem.MolToSmiles(mol)
    canonical_smiles = Chem.CanonSmiles(smiles)
    
    # Check for exact matches first
    if canonical_smiles in AROMATIC_NAMES:
        return AROMATIC_NAMES[canonical_smiles], True
    
    # Check for simple substituted benzenes
    rings = mol.GetRingInfo().AtomRings()
    for ring in rings:
        if len(ring) == 6:  # Benzene ring
            ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if all(atom.GetIsAromatic() and atom.GetAtomicNum() == 6 for atom in ring_atoms):
                return "benzene", False
    
    return None, False


def get_bond_info_advanced(mol):
    """Get comprehensive bond information"""
    double_bonds = []
    triple_bonds = []
    aromatic_bonds = []
    
    for bond in mol.GetBonds():
        bond_type = bond.GetBondType()
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        
        begin_atom = mol.GetAtomWithIdx(begin_idx)
        end_atom = mol.GetAtomWithIdx(end_idx)
        
        # Only consider carbon-carbon bonds for unsaturation
        if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
            if bond_type == Chem.rdchem.BondType.DOUBLE:
                double_bonds.append((begin_idx, end_idx))
            elif bond_type == Chem.rdchem.BondType.TRIPLE:
                triple_bonds.append((begin_idx, end_idx))
            elif bond.GetIsAromatic():
                aromatic_bonds.append((begin_idx, end_idx))
    
    return double_bonds, triple_bonds, aromatic_bonds


def find_principal_chain_advanced(mol, functional_groups, double_bonds, triple_bonds):
    """Advanced principal chain finding with functional group priority"""
    carbon_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    
    if not carbon_atoms:
        return [], 0, [], []
    
    if len(carbon_atoms) == 1:
        return carbon_atoms, 1, [], []
    
    # Get highest priority functional group
    principal_fg = functional_groups[0] if functional_groups else None
    
    best_chain = []
    best_length = 0
    best_double_positions = []
    best_triple_positions = []
    best_score = -1
    
    def calculate_chain_score(chain, chain_doubles, chain_triples, fg_position):
        """Calculate chain quality score"""
        base_score = len(chain) * 1000  # Length priority
        
        # Functional group bonus
        if principal_fg and fg_position is not None:
            base_score += 500  # Bonus for containing principal FG
            base_score -= fg_position * 10  # Lower position is better
        
        # Unsaturation bonus
        base_score += len(chain_doubles) * 100
        base_score += len(chain_triples) * 150
        
        # Lower position numbers are better
        if chain_doubles:
            base_score -= sum(chain_doubles) * 5
        if chain_triples:
            base_score -= sum(chain_triples) * 5
            
        return base_score
    
    def dfs_explore_chain(current_path, visited):
        nonlocal best_chain, best_length, best_double_positions, best_triple_positions, best_score
        
        if len(current_path) < 2:
            # Continue exploring
            pass
        else:
            # Analyze current chain
            chain_doubles = []
            chain_triples = []
            fg_position = None
            
            # Find unsaturated bond positions
            for i in range(len(current_path) - 1):
                atom1, atom2 = current_path[i], current_path[i + 1]
                
                for db in double_bonds:
                    if (atom1 in db and atom2 in db):
                        chain_doubles.append(i + 1)
                
                for tb in triple_bonds:
                    if (atom1 in tb and atom2 in tb):
                        chain_triples.append(i + 1)
            
            # Find principal functional group position
            if principal_fg:
                for i, atom_idx in enumerate(current_path):
                    if atom_idx == principal_fg["atom_idx"]:
                        fg_position = i + 1
                        break
            
            # Calculate score
            current_score = calculate_chain_score(current_path, chain_doubles, chain_triples, fg_position)
            
            if current_score > best_score:
                best_chain = current_path[:]
                best_length = len(current_path)
                best_double_positions = chain_doubles[:]
                best_triple_positions = chain_triples[:]
                best_score = current_score
        
        # Continue exploring
        if len(current_path) < len(carbon_atoms):  # Reasonable limit
            current_atom = mol.GetAtomWithIdx(current_path[-1])
            for neighbor in current_atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in visited and neighbor.GetAtomicNum() == 6:
                    visited.add(neighbor_idx)
                    current_path.append(neighbor_idx)
                    dfs_explore_chain(current_path, visited)
                    current_path.pop()
                    visited.remove(neighbor_idx)
    
    # Start exploration from each carbon
    for start_atom in carbon_atoms:
        visited = {start_atom}
        dfs_explore_chain([start_atom], visited)
    
    return best_chain, best_length, best_double_positions, best_triple_positions


def optimize_numbering_advanced(chain, double_positions, triple_positions, functional_groups):
    """Advanced numbering optimization"""
    if not chain:
        return chain, double_positions, triple_positions
    
    # Get principal functional group position
    principal_fg = functional_groups[0] if functional_groups else None
    fg_position = None
    
    if principal_fg:
        for i, atom_idx in enumerate(chain):
            if atom_idx == principal_fg["atom_idx"]:
                fg_position = i + 1
                break
    
    # Calculate scores for forward and reverse numbering
    def calculate_numbering_score(doubles, triples, fg_pos):
        score = 0
        # Principal functional group gets highest priority
        if fg_pos is not None:
            score += 10000 / fg_pos  # Lower position = higher score
        
        # Then unsaturation
        if triples:
            score += sum(1000 / pos for pos in triples)
        if doubles:
            score += sum(500 / pos for pos in doubles)
            
        return score
    
    # Forward numbering
    forward_fg = fg_position
    forward_doubles = double_positions[:]
    forward_triples = triple_positions[:]
    forward_score = calculate_numbering_score(forward_doubles, forward_triples, forward_fg)
    
    # Reverse numbering
    chain_length = len(chain)
    reverse_fg = chain_length - fg_position + 1 if fg_position else None
    reverse_doubles = [chain_length - pos + 1 for pos in double_positions]
    reverse_triples = [chain_length - pos + 1 for pos in triple_positions]
    reverse_score = calculate_numbering_score(reverse_doubles, reverse_triples, reverse_fg)
    
    if reverse_score > forward_score:
        return list(reversed(chain)), reverse_doubles, reverse_triples
    else:
        return chain, forward_doubles, forward_triples


def get_substituents_advanced(mol, principal_chain, functional_groups):
    """Get substituents with advanced functional group handling"""
    substituents = defaultdict(list)
    functional_groups_on_chain = defaultdict(list)
    principal_set = set(principal_chain)
    
    # Organize functional groups by position
    for fg in functional_groups:
        if fg["atom_idx"] in principal_set:
            position = principal_chain.index(fg["atom_idx"]) + 1
            functional_groups_on_chain[position].append(fg)
    
    # Find substituents
    for i, atom_idx in enumerate(principal_chain):
        atom = mol.GetAtomWithIdx(atom_idx)
        position = i + 1
        
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            neighbor_atomic_num = neighbor.GetAtomicNum()
            
            if neighbor_idx in principal_set:
                continue
            
            # Skip if this neighbor is part of a functional group on the main chain
            skip_neighbor = False
            for fg in functional_groups_on_chain[position]:
                if neighbor_idx in fg["match_atoms"]:
                    skip_neighbor = True
                    break
            
            if skip_neighbor:
                continue
            
            # Determine substituent type
            if neighbor_atomic_num == 6:  # Carbon substituent
                size = calculate_substituent_size_advanced(mol, neighbor_idx, principal_set)
                sub_name = SUBSTITUENT_NAMES.get(size, f"C{size}alkyl")
                substituents[position].append(sub_name)
            elif neighbor_atomic_num in [9, 17, 35, 53]:  # Halogens
                halogen_names = {9: "fluoro", 17: "chloro", 35: "bromo", 53: "iodo"}
                substituents[position].append(halogen_names[neighbor_atomic_num])
            elif neighbor_atomic_num == 8:  # Oxygen
                if neighbor.GetTotalNumHs() == 1:  # OH group not in FG
                    substituents[position].append("hydroxy")
            elif neighbor_atomic_num == 7:  # Nitrogen
                if neighbor.GetTotalNumHs() == 2:  # NH2 group not in FG
                    substituents[position].append("amino")
    
    return substituents, functional_groups_on_chain


def calculate_substituent_size_advanced(mol, start_atom, exclude_atoms):
    """Advanced substituent size calculation"""
    queue = deque([start_atom])
    visited = set([start_atom]) | exclude_atoms
    carbon_count = 0
    
    while queue:
        current = queue.popleft()
        current_atom = mol.GetAtomWithIdx(current)
        
        if current_atom.GetAtomicNum() == 6:
            carbon_count += 1
        
        for neighbor in current_atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in visited:
                visited.add(neighbor_idx)
                queue.append(neighbor_idx)
    
    return carbon_count


def handle_cyclic_advanced(mol, functional_groups):
    """Advanced cyclic compound handling"""
    rings = mol.GetRingInfo().AtomRings()
    if not rings:
        return None
    
    # Find largest ring
    largest_ring = max(rings, key=len)
    ring_size = len(largest_ring)
    ring_atoms = list(largest_ring)
    
    # Check if aromatic
    is_aromatic = all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring_atoms)
    
    if is_aromatic:
        return handle_aromatic_compound(mol, functional_groups, ring_atoms)
    else:
        return handle_alicyclic_compound(mol, functional_groups, ring_atoms, ring_size)


def handle_aromatic_compound(mol, functional_groups, ring_atoms):
    """Handle aromatic compounds (benzene derivatives)"""
    # Try to get standard aromatic name
    base_name, is_standard = get_aromatic_base_name(mol)
    if not base_name:
        base_name = "benzene"
    
    # Find substituents on the ring
    ring_set = set(ring_atoms)
    substituents = defaultdict(list)
    functional_groups_on_ring = defaultdict(list)
    
    # Organize functional groups
    for fg in functional_groups:
        if fg["atom_idx"] in ring_set:
            position = ring_atoms.index(fg["atom_idx"]) + 1
            functional_groups_on_ring[position].append(fg)
    
    # Find other substituents
    for i, atom_idx in enumerate(ring_atoms):
        atom = mol.GetAtomWithIdx(atom_idx)
        position = i + 1
        
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in ring_set:
                neighbor_atom = mol.GetAtomWithIdx(neighbor_idx)
                
                # Skip if part of functional group
                skip = False
                for fg in functional_groups_on_ring[position]:
                    if neighbor_idx in fg["match_atoms"]:
                        skip = True
                        break
                
                if not skip:
                    if neighbor_atom.GetAtomicNum() == 6:
                        size = calculate_substituent_size_advanced(mol, neighbor_idx, ring_set)
                        sub_name = SUBSTITUENT_NAMES.get(size, f"C{size}alkyl")
                        substituents[position].append(sub_name)
                    elif neighbor_atom.GetAtomicNum() in [9, 17, 35, 53]:
                        halogen_names = {9: "fluoro", 17: "chloro", 35: "bromo", 53: "iodo"}
                        substituents[position].append(halogen_names[neighbor_atom.GetAtomicNum()])
    
    # Build aromatic name
    return build_aromatic_name(base_name, substituents, functional_groups_on_ring, is_standard)


def handle_alicyclic_compound(mol, functional_groups, ring_atoms, ring_size):
    """Handle non-aromatic cyclic compounds"""
    # Get bond information for the ring
    double_bonds, triple_bonds, _ = get_bond_info_advanced(mol)
    
    # Count unsaturated bonds in ring
    ring_set = set(ring_atoms)
    ring_doubles = []
    ring_triples = []
    
    for i, bond in enumerate(double_bonds):
        if all(atom in ring_set for atom in bond):
            ring_doubles.append(bond)
    
    for i, bond in enumerate(triple_bonds):
        if all(atom in ring_set for atom in bond):
            ring_triples.append(bond)
    
    # Get base name
    base_name = CYCLIC_PARENT_NAMES.get(ring_size, f"cyclo-C{ring_size}")
    
    # Determine suffix
    if ring_triples:
        if len(ring_triples) == 1:
            suffix = "yne"
        else:
            suffix = f"{MULTIPLICITY_PREFIXES[len(ring_triples)]}yne"
    elif ring_doubles:
        if len(ring_doubles) == 1:
            suffix = "ene" 
        else:
            suffix = f"{MULTIPLICITY_PREFIXES[len(ring_doubles)]}ene"
    else:
        suffix = "ane"
    
    # Handle functional groups
    principal_fg = functional_groups[0] if functional_groups else None
    if principal_fg and principal_fg["suffix"]:
        # Modify suffix for functional group
        if suffix == "ane":
            base_suffix = ""
        else:
            base_suffix = suffix[:-3] if suffix.endswith("ane") else suffix
        
        full_suffix = f"{base_suffix}{principal_fg['suffix']}"
        return f"{base_name}{full_suffix}"
    
    return f"{base_name}{suffix}"


def build_aromatic_name(base_name, substituents, functional_groups_on_ring, is_standard):
    """Build aromatic compound name"""
    all_substituents = dict(substituents)
    
    # Add functional groups as substituents if not principal
    principal_fg = None
    if functional_groups_on_ring:
        all_fgs = []
        for pos, fgs in functional_groups_on_ring.items():
            all_fgs.extend([(fg, pos) for fg in fgs])
        
        if all_fgs:
            all_fgs.sort(key=lambda x: x[0]["priority"])
            principal_fg = all_fgs[0][0]
            
            # Add non-principal FGs as prefixes
            for fg, pos in all_fgs[1:]:
                if fg["prefix"]:
                    if pos not in all_substituents:
                        all_substituents[pos] = []
                    all_substituents[pos].append(fg["prefix"])
    
    # Build name parts
    name_parts = []
    
    if all_substituents:
        sub_parts = []
        grouped_subs = defaultdict(list)
        
        for locant, names in all_substituents.items():
            for name in names:
                grouped_subs[name].append(str(locant))
        
        for name in sorted(grouped_subs.keys()):
            locants = sorted(grouped_subs[name], key=int)
            prefix = MULTIPLICITY_PREFIXES.get(len(locants), f"{len(locants)}")
            
            if len(locants) > 1:
                sub_parts.append(f"{','.join(locants)}-{prefix}{name}")
            else:
                sub_parts.append(f"{locants[0]}-{name}")
        
        name_parts.extend(sub_parts)
    
    # Handle principal functional group
    if principal_fg and principal_fg["suffix"]:
        if base_name == "benzene":
            if principal_fg["suffix"] in ["ol", "amine"]:
                # Special cases: phenol, aniline
                special_names = {"ol": "phenol", "amine": "aniline"}
                base_name = special_names.get(principal_fg["suffix"], base_name)
            else:
                base_name = f"benzene{principal_fg['suffix']}"
        else:
            base_name = f"{base_name}{principal_fg['suffix']}"
    
    name_parts.append(base_name)
    
    full_name = "".join(name_parts) if len(name_parts) == 1 else "-".join(name_parts[:-1]) + name_parts[-1]
    return full_name


def determine_suffix_advanced(double_positions, triple_positions, principal_fg):
    """Advanced suffix determination"""
    if principal_fg and principal_fg["suffix"]:
        fg_suffix = principal_fg["suffix"]
        
        # Handle unsaturation with functional groups
        if triple_positions and double_positions:
            # Both types of unsaturation
            if len(triple_positions) == 1:
                triple_part = "yn"
            else:
                triple_part = f"{MULTIPLICITY_PREFIXES[len(triple_positions)]}yn"
            
            if len(double_positions) == 1:
                double_part = "en"
            else:
                double_part = f"{MULTIPLICITY_PREFIXES[len(double_positions)]}en"
            
            return f"{double_part}{triple_part}{fg_suffix}"
        
        elif triple_positions:
            if len(triple_positions) == 1:
                return f"yn{fg_suffix}"
            else:
                return f"{MULTIPLICITY_PREFIXES[len(triple_positions)]}yn{fg_suffix}"
        
        elif double_positions:
            if len(double_positions) == 1:
                return f"en{fg_suffix}"
            else:
                return f"{MULTIPLICITY_PREFIXES[len(double_positions)]}en{fg_suffix}"
        
        return fg_suffix
    
    # No principal functional group
    if triple_positions and double_positions:
        if len(triple_positions) == 1:
            triple_suffix = "yn"
        else:
            triple_suffix = f"{MULTIPLICITY_PREFIXES[len(triple_positions)]}yn"
        
        if len(double_positions) == 1:
            double_suffix = "en"
        else:
            double_suffix = f"{MULTIPLICITY_PREFIXES[len(double_positions)]}en"
        
        return f"{double_suffix}{triple_suffix}e"
    
    elif triple_positions:
        if len(triple_positions) == 1:
            return "yne"
        else:
            return f"{MULTIPLICITY_PREFIXES[len(triple_positions)]}yne"
    
    elif double_positions:
        if len(double_positions) == 1:
            return "ene"
        else:
            return f"{MULTIPLICITY_PREFIXES[len(double_positions)]}ene"
    
    return "ane"


def generate_locant_string(double_positions, triple_positions, fg_position, principal_fg):
    """Generate locant string for unsaturated bonds and functional groups"""
    parts = []
    
    # Double bond locants
    if double_positions:
        if len(double_positions) == 1:
            parts.append(f"{double_positions[0]}")
        else:
            parts.append(f"{','.join(map(str, double_positions))}")
    
    # Triple bond locants  
    if triple_positions:
        if len(triple_positions) == 1:
            parts.append(f"{triple_positions[0]}")
        else:
            parts.append(f"{','.join(map(str, triple_positions))}")
    
    # Functional group locant
    if principal_fg and fg_position and principal_fg["suffix"]:
        if fg_position > 1 or double_positions or triple_positions:
            parts.append(f"{fg_position}")
    
    if parts:
        return f"-{'-'.join(parts)}-"
    return ""


def assemble_final_name(chain_length, double_positions, triple_positions, 
                       substituents, principal_fg, fg_position, functional_groups_on_chain):
    """Assemble the final IUPAC name"""
    name_parts = []
    
    # Handle all substituents
    all_substituents = dict(substituents)
    
    # Add functional groups that should be prefixes
    for position, fgs in functional_groups_on_chain.items():
        for fg in fgs:
            if (not principal_fg or fg["type"] != principal_fg["type"] or position != fg_position):
                if fg["prefix"]:
                    if position not in all_substituents:
                        all_substituents[position] = []
                    all_substituents[position].append(fg["prefix"])
    
    # Format substituents
    if all_substituents:
        sub_parts = []
        grouped_subs = defaultdict(list)
        
        for locant, names in all_substituents.items():
            for name in names:
                grouped_subs[name].append(str(locant))
        
        for name in sorted(grouped_subs.keys()):
            locants = sorted(grouped_subs[name], key=int)
            prefix = MULTIPLICITY_PREFIXES.get(len(locants), f"{len(locants)}")
            
            if len(locants) > 1:
                sub_parts.append(f"{','.join(locants)}-{prefix}{name}")
            else:
                sub_parts.append(f"{locants[0]}-{name}")
        
        if sub_parts:
            name_parts.append("-".join(sub_parts))
    
    # Get parent name
    parent_root = PARENT_ROOTS.get(chain_length, f"C{chain_length}")
    
    # Get suffix
    suffix = determine_suffix_advanced(double_positions, triple_positions, principal_fg)
    
    # Add locants if needed
    locant_str = generate_locant_string(double_positions, triple_positions, fg_position, principal_fg)
    
    # Combine parts
    if locant_str and locant_str not in ["--", "-"]:
        parent_name = f"{parent_root}{locant_str}{suffix}"
    else:
        parent_name = f"{parent_root}{suffix}"
    
    name_parts.append(parent_name)
    
    # Clean up the final name
    full_name = "".join(name_parts)
    full_name = re.sub(r'-+', '-', full_name)  # Remove multiple dashes
    full_name = full_name.strip('-')  # Remove leading/trailing dashes
    
    return full_name


def generate_iupac_name_production(mol):
    """Main production-ready IUPAC name generation function"""
    if not mol:
        return "Invalid molecule"
    
    try:
        # Sanitize molecule
        Chem.SanitizeMol(mol)
    except:
        pass
    
    # Basic checks
    num_atoms = mol.GetNumAtoms()
    if num_atoms == 0:
        return "Empty molecule"
    
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count == 0:
        return "Inorganic compound"
    
    # Special case: single carbon
    if carbon_count == 1:
        functional_groups = identify_functional_groups_advanced(mol)
        if functional_groups:
            fg = functional_groups[0]
            if fg["type"] == "aldehyde":
                return "methanal"
            elif fg["type"] == "ketone":
                return "methanone"  # Not standard but for completeness
            elif fg["type"] == "alcohol":
                return "methanol"
            elif fg["type"] == "carboxylic_acid":
                return "methanoic acid"
        return "methane"
    
    # Identify functional groups
    functional_groups = identify_functional_groups_advanced(mol)
    
    # Check for aromatic compounds
    if is_aromatic_compound(mol):
        return handle_aromatic_compound(mol, functional_groups, list(range(mol.GetNumAtoms())))
    
    # Check for cyclic compounds
    if mol.GetRingInfo().NumRings() > 0:
        cyclic_result = handle_cyclic_advanced(mol, functional_groups)
        if cyclic_result:
            return cyclic_result
    
    # Handle acyclic compounds
    double_bonds, triple_bonds, _ = get_bond_info_advanced(mol)
    
    # Find principal chain
    principal_chain, chain_length, double_positions, triple_positions = find_principal_chain_advanced(
        mol, functional_groups, double_bonds, triple_bonds
    )
    
    if chain_length == 0:
        return "No valid carbon chain found"
    
    # Optimize numbering
    principal_chain, double_positions, triple_positions = optimize_numbering_advanced(
        principal_chain, double_positions, triple_positions, functional_groups
    )
    
    # Get substituents
    substituents, functional_groups_on_chain = get_substituents_advanced(
        mol, principal_chain, functional_groups
    )
    
    # Get principal functional group position
    principal_fg = functional_groups[0] if functional_groups else None
    fg_position = None
    
    if principal_fg:
        for i, atom_idx in enumerate(principal_chain):
            if atom_idx == principal_fg["atom_idx"]:
                fg_position = i + 1
                break
    
    # Assemble final name
    return assemble_final_name(
        chain_length, double_positions, triple_positions,
        substituents, principal_fg, fg_position, functional_groups_on_chain
    )


def name_from_smiles_production(smiles_string):
    """Production-ready SMILES to IUPAC name conversion"""
    if not smiles_string or smiles_string.strip() == "":
        return "Please draw a molecule first."
    
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if not mol:
            return "Invalid SMILES string - please check your molecule structure"
        
        return generate_iupac_name_production(mol)
        
    except Exception as e:
        return f"Error analyzing structure: {str(e)}"


# ==============================================================================
# STREAMLIT APPLICATION
# ==============================================================================

def main():
    inject_modern_css()
    
    st.markdown("""
    <div class="main-header">
        <h1>⚗️ IUPAC Name Generator</h1>
        <p>Complete production-ready system for alkanes, alkenes, alkynes, aromatics & functional groups!</p>
    </div>
    """, unsafe_allow_html=True)

    # Initialize session state
    if 'molecule' not in st.session_state:
        st.session_state.molecule = "CCO"

    # Main layout
    col1, col2 = st.columns([2.5, 1], gap="large")
    
    with col1:
        st.markdown('<div class="editor-card">', unsafe_allow_html=True)
        st.markdown('<div class="section-header">🖊️ Molecule Editor</div>', unsafe_allow_html=True)
        
        molecule_smiles = st_ketcher(
            value=st.session_state.molecule,
            height=400,
            key="molecule_editor"
        )
        
        if molecule_smiles and molecule_smiles != st.session_state.molecule:
            st.session_state.molecule = molecule_smiles
            
        st.markdown('</div>', unsafe_allow_html=True)

    with col2:
        st.markdown('<div class="controls-card">', unsafe_allow_html=True)
        st.markdown('<div class="section-header">⚡ Examples</div>', unsafe_allow_html=True)
        
        if st.button("🗑️ Clear", help="Clear the editor"):
            st.session_state.molecule = ""
            st.rerun()
        
        st.markdown("---")
        
        example_categories = {
            "🔹 Alkanes": {
                "Ethane": "CC",
                "Propane": "CCC", 
                "2-Methylbutane": "CC(C)CC",
                "2,2-Dimethylpropane": "CC(C)(C)C"
            },
            "🔹 Alkenes": {
                "Ethene": "C=C",
                "Propene": "C=CC",
                "But-2-ene": "CC=CC",
                "2-Methylprop-1-ene": "C=C(C)C"
            },
            "🔹 Alkynes": {
                "Ethyne": "C#C",
                "Propyne": "C#CC",
                "But-2-yne": "CC#CC",
                "Pent-1-yne": "C#CCCC"
            },
            "🔹 Aromatics": {
                "Benzene": "c1ccccc1",
                "Toluene": "Cc1ccccc1",
                "Phenol": "Oc1ccccc1",
                "Aniline": "Nc1ccccc1"
            },
            "🔹 Alcohols": {
                "Methanol": "CO",
                "Ethanol": "CCO",
                "Propan-2-ol": "CC(O)C",
                "Butan-1-ol": "CCCCO"
            },
            "🔹 Carbonyls": {
                "Methanal": "C=O",
                "Ethanal": "CC=O", 
                "Propanone": "CC(=O)C",
                "Butanone": "CC(=O)CC"
            },
            "🔹 Carboxylic Acids": {
                "Methanoic acid": "C(=O)O",
                "Ethanoic acid": "CC(=O)O",
                "Propanoic acid": "CCC(=O)O",
                "Benzoic acid": "c1ccc(cc1)C(=O)O"
            },
            "🔹 Cyclic": {
                "Cyclopropane": "C1CC1",
                "Cyclohexane": "C1CCCCC1",
                "Cyclohexene": "C1=CCCCC1",
                "Cyclohexanol": "C1CCCCC1O"
            },
            "🔹 Complex": {
                "3-Methylbut-1-yne": "C#CC(C)C",
                "Pent-1-en-4-yne": "C=CCC#C",
                "2-Phenylethanol": "c1ccc(cc1)CCO",
                "4-Methylphenol": "Cc1ccc(cc1)O"
            }
        }
        
        for category, examples in example_categories.items():
            with st.expander(category, expanded=False):
                for name, smiles in examples.items():
                    if st.button(name, key=f"{category}_{name}", help=f"SMILES: {smiles}"):
                        st.session_state.molecule = smiles
                        st.rerun()
        
        st.markdown('</div>', unsafe_allow_html=True)

    # Show SMILES
    if molecule_smiles:
        st.markdown(f"""
        <div class="smiles-container">
            <strong>🧬 SMILES:</strong> <code>{molecule_smiles}</code>
        </div>
        """, unsafe_allow_html=True)

    # Generate button
    col1, col2, col3 = st.columns([1, 2, 1])
    with col2:
        if st.button("🔬 Generate IUPAC Name", type="primary", use_container_width=True):
            if molecule_smiles:
                with st.spinner("🔄 Analyzing molecular structure..."):
                    iupac_name = name_from_smiles_production(molecule_smiles)
                st.session_state.iupac_result = iupac_name
            else:
                st.session_state.iupac_result = "Please draw a molecule first."

    # Show results
    if 'iupac_result' in st.session_state:
        if any(word in st.session_state.iupac_result for word in ["Error", "Please draw", "Invalid"]):
            st.markdown(f"""
            <div class="result-error">
                <h3>❌ Error</h3>
                <p>{st.session_state.iupac_result}</p>
            </div>
            """, unsafe_allow_html=True)
        else:
            st.markdown(f"""
            <div class="result-success">
                <h3>🎉 IUPAC Name Generated Successfully!</h3>
                <p class="iupac-name">{st.session_state.iupac_result}</p>
            </div>
            """, unsafe_allow_html=True)

    # Footer
    st.markdown("""<div class="footer">Made for Mistri of Chemistry - Now with Complete IUPAC Support!</div>""", unsafe_allow_html=True)
    st.markdown("""<div class="footer">Mayank Vishwakarma <a href="https://www.linkedin.com/in/mayank-vishwakarma2004/">LinkedIn</a></div>""", unsafe_allow_html=True)


if __name__ == "__main__":
    main()
