import streamlit as st
from streamlit_ketcher import st_ketcher
from rdkit import Chem
from rdkit.Chem import GetPeriodicTable
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
# CORE IUPAC NAMING FUNCTIONS - IMPROVED FOR ALKENES/ALKYNES
# ==============================================================================

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
    1: "methyl", 2: "ethyl", 3: "propyl", 4: "butyl", 5: "pentyl"
}

PREFIXES = {1: "", 2: "di", 3: "tri", 4: "tetra", 5: "penta", 6: "hexa"}

MULTIPLICITY_PREFIXES = {
    1: "", 2: "di", 3: "tri", 4: "tetra", 5: "penta", 6: "hexa"
}

# Functional group patterns and priorities
FUNCTIONAL_GROUP_SMARTS = {
    "aldehyde": "[CX3H1](=O)",
    "ketone": "[CX3](=O)[#6]",
    "alcohol": "[OH]",
    "carboxylic_acid": "[CX3](=O)[OH]",
    "ester": "[CX3](=O)[OX2H0]",
    "amine_primary": "[NX3;H2]",
    "halogen_f": "[F]",
    "halogen_cl": "[Cl]",
    "halogen_br": "[Br]",
    "halogen_i": "[I]"
}

FUNCTIONAL_GROUP_PRIORITIES = {
    "carboxylic_acid": 1, "ester": 2, "aldehyde": 3, "ketone": 4,
    "alcohol": 5, "amine_primary": 6, "halogen_f": 10,
    "halogen_cl": 10, "halogen_br": 10, "halogen_i": 10,
}

FUNCTIONAL_GROUP_SUFFIXES = {
    "carboxylic_acid": "oic acid", "ester": "oate", "aldehyde": "al",
    "ketone": "one", "alcohol": "ol", "amine_primary": "amine"
}

FUNCTIONAL_GROUP_PREFIXES = {
    "halogen_f": "fluoro", "halogen_cl": "chloro",
    "halogen_br": "bromo", "halogen_i": "iodo", "alcohol": "hydroxy"
}


def identify_functional_groups_smarts(mol):
    """Identify functional groups using SMARTS patterns"""
    functional_groups = []
    for fg_name, smarts_pattern in FUNCTIONAL_GROUP_SMARTS.items():
        try:
            pattern = Chem.MolFromSmarts(smarts_pattern)
            if pattern:
                matches = mol.GetSubstructMatches(pattern)
                for match in matches:
                    functional_groups.append({
                        "type": fg_name,
                        "atom_idx": match[0],
                        "priority": FUNCTIONAL_GROUP_PRIORITIES.get(fg_name, 99),
                        "match_atoms": match
                    })
        except:
            continue
    return functional_groups


def get_bond_info(mol):
    """Extract information about double and triple bonds"""
    double_bonds = []
    triple_bonds = []
    
    for bond in mol.GetBonds():
        bond_type = bond.GetBondType()
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        
        # Only consider carbon-carbon bonds
        begin_atom = mol.GetAtomWithIdx(begin_idx)
        end_atom = mol.GetAtomWithIdx(end_idx)
        
        if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
            if bond_type == Chem.rdchem.BondType.DOUBLE:
                double_bonds.append((begin_idx, end_idx))
            elif bond_type == Chem.rdchem.BondType.TRIPLE:
                triple_bonds.append((begin_idx, end_idx))
    
    return double_bonds, triple_bonds


def is_cyclic_molecule(mol):
    """Check if molecule contains rings"""
    return mol.GetRingInfo().NumRings() > 0


def get_ring_systems(mol):
    """Get information about ring systems"""
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    ring_systems = []
    
    for ring in atom_rings:
        ring_size = len(ring)
        is_aromatic = all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
        ring_systems.append({
            "atoms": list(ring),
            "size": ring_size,
            "is_aromatic": is_aromatic
        })
    
    return ring_systems


def find_longest_chain_with_unsaturation(mol, double_bonds, triple_bonds):
    """Find the longest carbon chain that contains the most unsaturated bonds"""
    carbon_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    
    if not carbon_atoms:
        return [], 0, [], []
    
    if len(carbon_atoms) == 1:
        return carbon_atoms, 1, [], []
    
    # Create a set of all atoms involved in unsaturated bonds
    unsaturated_atoms = set()
    for bond in double_bonds + triple_bonds:
        unsaturated_atoms.update(bond)
    
    best_chain = []
    best_length = 0
    best_double_positions = []
    best_triple_positions = []
    
    # Try all possible chains
    def dfs_chain(current_path, visited):
        nonlocal best_chain, best_length, best_double_positions, best_triple_positions
        
        if len(current_path) > best_length:
            # Calculate unsaturated bond positions in this chain
            chain_doubles = []
            chain_triples = []
            
            for i in range(len(current_path) - 1):
                atom1, atom2 = current_path[i], current_path[i + 1]
                
                # Check for double bonds
                for db in double_bonds:
                    if (atom1 in db and atom2 in db):
                        chain_doubles.append(i + 1)  # Position of first carbon in bond
                
                # Check for triple bonds  
                for tb in triple_bonds:
                    if (atom1 in tb and atom2 in tb):
                        chain_triples.append(i + 1)  # Position of first carbon in bond
            
            # Prioritize chains with more unsaturated bonds
            current_unsaturated_count = len(chain_doubles) + len(chain_triples)
            best_unsaturated_count = len(best_double_positions) + len(best_triple_positions)
            
            if (len(current_path) > best_length or 
                (len(current_path) == best_length and current_unsaturated_count > best_unsaturated_count)):
                best_chain = current_path[:]
                best_length = len(current_path)
                best_double_positions = chain_doubles[:]
                best_triple_positions = chain_triples[:]
        
        # Continue building chain
        current_atom = mol.GetAtomWithIdx(current_path[-1])
        for neighbor in current_atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in visited and neighbor.GetAtomicNum() == 6:
                visited.add(neighbor_idx)
                current_path.append(neighbor_idx)
                dfs_chain(current_path, visited)
                current_path.pop()
                visited.remove(neighbor_idx)
    
    # Start DFS from each carbon atom
    for start_atom in carbon_atoms:
        visited = {start_atom}
        dfs_chain([start_atom], visited)
    
    return best_chain, best_length, best_double_positions, best_triple_positions


def optimize_chain_numbering(chain, double_positions, triple_positions):
    """Optimize chain numbering to give lowest numbers to unsaturated bonds"""
    if not double_positions and not triple_positions:
        return chain, double_positions, triple_positions
    
    # Try numbering from both directions
    forward_doubles = double_positions[:]
    forward_triples = triple_positions[:]
    
    # Reverse numbering
    chain_length = len(chain)
    reverse_doubles = [chain_length - pos for pos in double_positions]
    reverse_triples = [chain_length - pos for pos in triple_positions]
    
    # Compare which gives lower numbers (double bonds have priority)
    forward_sum = sum(forward_doubles) + sum(forward_triples)
    reverse_sum = sum(reverse_doubles) + sum(reverse_triples)
    
    if reverse_sum < forward_sum:
        return list(reversed(chain)), reverse_doubles, reverse_triples
    else:
        return chain, forward_doubles, forward_triples


def get_substituents_and_functional_groups(mol, principal_structure):
    """Get substituents and functional groups attached to the principal chain"""
    substituents = defaultdict(list)
    functional_groups_on_chain = defaultdict(list)
    principal_set = set(principal_structure)
    functional_groups = identify_functional_groups_smarts(mol)
    
    for i, atom_idx in enumerate(principal_structure):
        atom = mol.GetAtomWithIdx(atom_idx)
        position = i + 1
        
        # Find functional groups on this carbon
        carbon_functional_groups = [fg for fg in functional_groups if fg["atom_idx"] == atom_idx]
        for fg in carbon_functional_groups:
            functional_groups_on_chain[position].append(fg)
        
        # Find substituents
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            neighbor_atomic_num = neighbor.GetAtomicNum()
            
            if neighbor_idx in principal_set:
                continue
            
            if neighbor_atomic_num == 6:
                size = calculate_substituent_size(mol, neighbor_idx, principal_set)
                sub_name = SUBSTITUENT_NAMES.get(size, f"C{size}-alkyl")
                substituents[position].append(sub_name)
            elif neighbor_atomic_num in [9, 17, 35, 53]:  # Halogens
                halogen_names = {9: "fluoro", 17: "chloro", 35: "bromo", 53: "iodo"}
                substituents[position].append(halogen_names[neighbor_atomic_num])
            elif neighbor_atomic_num == 8:  # Oxygen
                # Check if it's part of a functional group
                is_part_of_fg = any(neighbor_idx in fg["match_atoms"] for fg in carbon_functional_groups)
                if not is_part_of_fg and neighbor.GetTotalNumHs() == 1:
                    substituents[position].append("hydroxy")
    
    return substituents, functional_groups_on_chain


def calculate_substituent_size(mol, start_atom, exclude_atoms):
    """Calculate the size of a substituent branch"""
    queue = deque([start_atom])
    visited = set([start_atom]) | exclude_atoms
    size = 0
    
    while queue:
        current = queue.popleft()
        current_atom = mol.GetAtomWithIdx(current)
        
        if current_atom.GetAtomicNum() == 6:
            size += 1
        
        for neighbor in current_atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in visited:
                visited.add(neighbor_idx)
                queue.append(neighbor_idx)
    
    return size


def determine_principal_functional_group(functional_groups_on_chain):
    """Determine the highest priority functional group"""
    all_fgs = []
    for position, fgs in functional_groups_on_chain.items():
        for fg in fgs:
            all_fgs.append((fg, position))
    
    if not all_fgs:
        return None, None
    
    # Sort by priority (lower number = higher priority)
    all_fgs.sort(key=lambda x: x[0]["priority"])
    return all_fgs[0][0], all_fgs[0][1]


def generate_iupac_name(mol):
    """Generate complete IUPAC name for the molecule"""
    if not mol:
        return "Invalid molecule"
    
    # Get bond information
    double_bonds, triple_bonds = get_bond_info(mol)
    
    # Check if cyclic
    if is_cyclic_molecule(mol):
        return handle_cyclic_compound(mol, double_bonds, triple_bonds)
    
    # Find principal chain with unsaturation priority
    principal_chain, chain_length, double_positions, triple_positions = find_longest_chain_with_unsaturation(
        mol, double_bonds, triple_bonds
    )
    
    if chain_length == 0:
        return "No valid carbon chain found"
    
    # Optimize numbering for unsaturated bonds
    principal_chain, double_positions, triple_positions = optimize_chain_numbering(
        principal_chain, double_positions, triple_positions
    )
    
    # Get substituents and functional groups
    substituents, functional_groups_on_chain = get_substituents_and_functional_groups(mol, principal_chain)
    
    # Determine principal functional group
    principal_fg, fg_position = determine_principal_functional_group(functional_groups_on_chain)
    
    # Build the name
    return assemble_iupac_name(
        chain_length, double_positions, triple_positions, 
        substituents, principal_fg, fg_position, functional_groups_on_chain
    )


def handle_cyclic_compound(mol, double_bonds, triple_bonds):
    """Handle cyclic compounds (simplified)"""
    ring_systems = get_ring_systems(mol)
    if not ring_systems:
        return "Complex cyclic structure"
    
    largest_ring = max(ring_systems, key=lambda x: x["size"])
    ring_size = largest_ring["size"]
    
    # Get parent name
    parent_root = CYCLIC_PARENT_NAMES.get(ring_size, f"cyclo-C{ring_size}")
    
    # Count unsaturated bonds in the ring
    ring_atoms = set(largest_ring["atoms"])
    double_count = sum(1 for bond in double_bonds if all(atom in ring_atoms for atom in bond))
    triple_count = sum(1 for bond in triple_bonds if all(atom in ring_atoms for atom in bond))
    
    # Determine suffix
    if triple_count > 0:
        if triple_count == 1:
            suffix = "yne"
        else:
            suffix = f"{MULTIPLICITY_PREFIXES[triple_count]}yne"
    elif double_count > 0:
        if double_count == 1:
            suffix = "ene"
        else:
            suffix = f"{MULTIPLICITY_PREFIXES[double_count]}ene"
    else:
        suffix = "ane"
    
    return f"{parent_root}{suffix}"


def assemble_iupac_name(chain_length, double_positions, triple_positions, 
                       substituents, principal_fg, fg_position, functional_groups_on_chain):
    """Assemble the final IUPAC name"""
    name_parts = []
    
    # Handle substituents (including functional groups as prefixes)
    all_substituents = dict(substituents)
    
    # Add functional groups that should be prefixes
    for position, fgs in functional_groups_on_chain.items():
        for fg in fgs:
            if (principal_fg is None or fg["type"] != principal_fg["type"] or position != fg_position):
                if fg["type"] in FUNCTIONAL_GROUP_PREFIXES:
                    prefix_name = FUNCTIONAL_GROUP_PREFIXES[fg["type"]]
                    if position not in all_substituents:
                        all_substituents[position] = []
                    all_substituents[position].append(prefix_name)
    
    # Format substituent names
    if all_substituents:
        substituent_parts = []
        grouped_subs = defaultdict(list)
        
        # Group substituents by name
        for locant, names in all_substituents.items():
            for name in names:
                grouped_subs[name].append(str(locant))
        
        # Create substituent strings
        for name in sorted(grouped_subs.keys()):
            locants = sorted(grouped_subs[name], key=int)
            prefix = MULTIPLICITY_PREFIXES.get(len(locants), f"{len(locants)}")
            
            if len(locants) > 1:
                substituent_parts.append(f"{','.join(locants)}-{prefix}{name}")
            else:
                substituent_parts.append(f"{locants[0]}-{name}")
        
        if substituent_parts:
            name_parts.append("-".join(substituent_parts))
    
    # Get parent chain name
    parent_root = PARENT_ROOTS.get(chain_length, f"C{chain_length}")
    
    # Determine suffix based on unsaturation and functional groups
    suffix = determine_suffix(double_positions, triple_positions, principal_fg)
    
    # Add locants for unsaturated bonds if needed
    locant_part = ""
    if double_positions or triple_positions:
        locant_part = generate_locants(double_positions, triple_positions, principal_fg, fg_position)
    
    # Combine parent name with locants and suffix
    if locant_part:
        parent_name = f"{parent_root}{locant_part}{suffix}"
    else:
        parent_name = f"{parent_root}{suffix}"
    
    name_parts.append(parent_name)
    
    # Join all parts
    full_name = "".join(name_parts)
    full_name = re.sub(r'-+', '-', full_name)  # Remove multiple consecutive dashes
    full_name = full_name.strip('-')  # Remove leading/trailing dashes
    
    return full_name


def determine_suffix(double_positions, triple_positions, principal_fg):
    """Determine the appropriate suffix based on unsaturation and functional groups"""
    # If there's a principal functional group with higher priority than unsaturation
    if principal_fg and principal_fg["type"] in FUNCTIONAL_GROUP_SUFFIXES:
        fg_suffix = FUNCTIONAL_GROUP_SUFFIXES[principal_fg["type"]]
        
        # For alcohols, ketones, etc., we may need to modify based on unsaturation
        if triple_positions:
            if len(triple_positions) == 1:
                base_suffix = "yn"
            else:
                base_suffix = f"{MULTIPLICITY_PREFIXES[len(triple_positions)]}yn"
            
            if double_positions:
                if len(double_positions) == 1:
                    return f"-{base_suffix}en-{fg_suffix}"
                else:
                    return f"-{base_suffix}{MULTIPLICITY_PREFIXES[len(double_positions)]}en-{fg_suffix}"
            else:
                return f"-{base_suffix}-{fg_suffix}"
        
        elif double_positions:
            if len(double_positions) == 1:
                base_suffix = "en"
            else:
                base_suffix = f"{MULTIPLICITY_PREFIXES[len(double_positions)]}en"
            return f"-{base_suffix}-{fg_suffix}"
        
        else:
            return f"-{fg_suffix}"
    
    # No high-priority functional group, use unsaturation-based suffix
    if triple_positions and double_positions:
        # Both double and triple bonds
        if len(triple_positions) == 1:
            triple_suffix = "yn"
        else:
            triple_suffix = f"{MULTIPLICITY_PREFIXES[len(triple_positions)]}yn"
        
        if len(double_positions) == 1:
            double_suffix = "en"
        else:
            double_suffix = f"{MULTIPLICITY_PREFIXES[len(double_positions)]}en"
        
        return f"a-{double_suffix}-{triple_suffix}e"  # e.g., "a-en-yne"
    
    elif triple_positions:
        # Only triple bonds
        if len(triple_positions) == 1:
            return "yne"
        else:
            return f"a{MULTIPLICITY_PREFIXES[len(triple_positions)]}yne"
    
    elif double_positions:
        # Only double bonds
        if len(double_positions) == 1:
            return "ene"
        else:
            return f"a{MULTIPLICITY_PREFIXES[len(double_positions)]}ene"
    
    else:
        # Saturated
        return "ane"


def generate_locants(double_positions, triple_positions, principal_fg, fg_position):
    """Generate locant numbers for unsaturated bonds"""
    locants = []
    
    if double_positions:
        if len(double_positions) == 1:
            locants.append(f"-{double_positions[0]}")
        else:
            locants.append(f"-{','.join(map(str, double_positions))}")
    
    if triple_positions:
        if len(triple_positions) == 1:
            locants.append(f"-{triple_positions[0]}")
        else:
            locants.append(f"-{','.join(map(str, triple_positions))}")
    
    # Add functional group locant if needed
    if principal_fg and fg_position and principal_fg["type"] in FUNCTIONAL_GROUP_SUFFIXES:
        if fg_position > 1 or (double_positions or triple_positions):
            locants.append(f"-{fg_position}")
    
    return "".join(locants)


def name_from_smiles_complete(smiles_string):
    """Main function to generate IUPAC name from SMILES"""
    if not smiles_string or smiles_string.strip() == "":
        return "Please draw a molecule first."
    
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if not mol:
            return "Invalid SMILES string - please check your molecule structure"
        
        # Sanitize molecule
        try:
            Chem.SanitizeMol(mol)
        except:
            pass
        
        # Check for atoms
        num_atoms = mol.GetNumAtoms()
        if num_atoms == 0:
            return "Empty molecule - please draw a structure"
        
        # Check for carbon atoms
        carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        if carbon_count == 0:
            return "Inorganic compound detected"
        
        # Handle special cases
        if carbon_count == 1:
            return "methane"
        
        # Generate IUPAC name
        return generate_iupac_name(mol)
        
    except Exception as e:
        return f"Error analyzing structure: Please try a simpler molecule."


# ==============================================================================
# STREAMLIT APPLICATION
# ==============================================================================

def main():
    inject_modern_css()
    
    st.markdown("""
    <div class="main-header">
        <h1>⚗️ IUPAC Name Generator</h1>
        <p>Production-ready alkane, alkene & alkyne naming system!</p>
    </div>
    """, unsafe_allow_html=True)

    # Initialize session state
    if 'molecule' not in st.session_state:
        st.session_state.molecule = "CC=O"

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
                "Butane": "CCCC",
                "2-Methylpropane": "CC(C)C"
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
            "🔹 Mixed Unsaturated": {
                "But-1-en-3-yne": "C#CC=C",
                "Pent-1-en-4-yne": "C=CCC#C",
                "Hexa-1,5-diene": "C=CCCC=C",
                "But-1,3-diyne": "C#CC#C"
            },
            "🔹 Functional Groups": {
                "Ethanol": "CCO",
                "Ethanal": "CC=O",
                "Propanone": "CC(=O)C",
                "Ethanoic acid": "CC(=O)O"
            },
            "🔹 Cyclic": {
                "Cyclopropane": "C1CC1",
                "Cyclohexene": "C1=CCCCC1",
                "Benzene": "c1ccccc1",
                "Cyclohexanol": "C1CCCCC1O"
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
                with st.spinner("🔄 Analyzing structure..."):
                    iupac_name = name_from_smiles_complete(molecule_smiles)
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
                <h3>🎉 IUPAC Name Generated</h3>
                <p class="iupac-name">{st.session_state.iupac_result}</p>
            </div>
            """, unsafe_allow_html=True)

    # Footer
    st.markdown("""<div class="footer">Made for Mistri of Chemistry</div>""", unsafe_allow_html=True)
    st.markdown("""<div class="footer">Mayank Vishwakarma <a href="https://www.linkedin.com/in/mayank-vishwakarma2004/">LinkedIn</a></div>""", unsafe_allow_html=True)


if __name__ == "__main__":
    main()
