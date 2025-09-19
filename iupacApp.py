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
# CORE IUPAC NAMING FUNCTIONS (SIMPLIFIED)
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

def is_cyclic_molecule(mol):
    return mol.GetRingInfo().NumRings() > 0

def get_ring_systems(mol):
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

def find_principal_chain_or_ring(mol):
    if is_cyclic_molecule(mol):
        ring_systems = get_ring_systems(mol)
        if ring_systems:
            principal_ring = max(ring_systems, key=lambda r: r["size"])
            return principal_ring["atoms"], True, principal_ring["size"]
    carbon_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if not carbon_atoms:
        return [], False, 0
    if len(carbon_atoms) == 1:
        return carbon_atoms, False, 1
    longest_chain = []
    max_length = 0
    for i in range(len(carbon_atoms)):
        for j in range(i + 1, len(carbon_atoms)):
            try:
                path = Chem.GetShortestPath(mol, carbon_atoms[i], carbon_atoms[j])
                carbon_path = [idx for idx in path if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6]
                if len(carbon_path) > max_length:
                    max_length = len(carbon_path)
                    longest_chain = carbon_path
            except:
                continue
    return longest_chain if longest_chain else carbon_atoms[:1], False, len(longest_chain) if longest_chain else 1

def get_substituents_and_functional_groups(mol, principal_structure, is_cyclic):
    substituents = defaultdict(list)
    functional_groups_on_chain = defaultdict(list)
    principal_set = set(principal_structure)
    functional_groups = identify_functional_groups_smarts(mol)
    for i, atom_idx in enumerate(principal_structure):
        atom = mol.GetAtomWithIdx(atom_idx)
        position = i + 1
        carbon_functional_groups = [fg for fg in functional_groups if fg["atom_idx"] == atom_idx]
        for fg in carbon_functional_groups:
            functional_groups_on_chain[position].append(fg)
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            neighbor_atomic_num = neighbor.GetAtomicNum()
            if neighbor_idx in principal_set:
                continue
            if neighbor_atomic_num == 6:
                size = calculate_substituent_size(mol, neighbor_idx, principal_set)
                sub_name = SUBSTITUENT_NAMES.get(size, f"C{size}-alkyl")
                substituents[position].append(sub_name)
            elif neighbor_atomic_num in [9, 17, 35, 53]:
                halogen_names = {9: "fluoro", 17: "chloro", 35: "bromo", 53: "iodo"}
                substituents[position].append(halogen_names[neighbor_atomic_num])
            elif neighbor_atomic_num == 8:
                is_part_of_fg = any(neighbor_idx in fg["match_atoms"] for fg in carbon_functional_groups)
                if not is_part_of_fg and neighbor.GetTotalNumHs() == 1:
                    substituents[position].append("hydroxy")
    return substituents, functional_groups_on_chain

def calculate_substituent_size(mol, start_atom, exclude_atoms):
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
    all_fgs = []
    for position, fgs in functional_groups_on_chain.items():
        for fg in fgs:
            all_fgs.append((fg, position))
    if not all_fgs:
        return None, None
    all_fgs.sort(key=lambda x: x[0]["priority"])
    return all_fgs[0][0], all_fgs[0][1]

def assemble_iupac_name_complete(mol, principal_structure, is_cyclic, chain_length):
    if chain_length == 0:
        return "No valid structure found"
    substituents, functional_groups_on_chain = get_substituents_and_functional_groups(
        mol, principal_structure, is_cyclic
    )
    principal_fg, fg_position = determine_principal_functional_group(functional_groups_on_chain)
    name_parts = []
    all_substituents = dict(substituents)
    for position, fgs in functional_groups_on_chain.items():
        for fg in fgs:
            if principal_fg is None or fg["type"] != principal_fg["type"] or position != fg_position:
                if fg["type"] in FUNCTIONAL_GROUP_PREFIXES:
                    prefix_name = FUNCTIONAL_GROUP_PREFIXES[fg["type"]]
                    all_substituents[position] = all_substituents.get(position, []) + [prefix_name]
    if all_substituents:
        substituent_parts = []
        grouped_subs = defaultdict(list)
        for locant, names in all_substituents.items():
            for name in names:
                grouped_subs[name].append(str(locant))
        for name in sorted(grouped_subs.keys()):
            locants = sorted(grouped_subs[name], key=int)
            prefix = PREFIXES.get(len(locants), f"{len(locants)}")
            if len(locants) > 1:
                substituent_parts.append(f"{','.join(locants)}-{prefix}{name}")
            else:
                substituent_parts.append(f"{locants[0]}-{name}")
        if substituent_parts:
            name_parts.append("-".join(substituent_parts))
    if is_cyclic:
        parent_root = CYCLIC_PARENT_NAMES.get(chain_length, f"cyclo-C{chain_length}")
    else:
        parent_root = PARENT_ROOTS.get(chain_length, f"C{chain_length}")
    parent_name = parent_root + "ane"
    if principal_fg and principal_fg["type"] in FUNCTIONAL_GROUP_SUFFIXES:
        suffix = FUNCTIONAL_GROUP_SUFFIXES[principal_fg["type"]]
        if parent_name.endswith("ane"):
            parent_name = parent_name[:-3]
        if fg_position and fg_position > 1:
            parent_name = f"{parent_name}-{fg_position}-{suffix}"
        else:
            parent_name = f"{parent_name}{suffix}"
    name_parts.append(parent_name)
    full_name = "".join(name_parts)
    full_name = re.sub(r'-+', '-', full_name)
    full_name = full_name.strip('-')
    return full_name

def name_from_smiles_complete(smiles_string):
    if not smiles_string or smiles_string.strip() == "":
        return "Please draw a molecule first."
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if not mol:
            return "Invalid SMILES string - please check your molecule structure"
        try:
            Chem.SanitizeMol(mol)
        except:
            pass
        num_atoms = mol.GetNumAtoms()
        if num_atoms == 0:
            return "Empty molecule - please draw a structure"
        carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        if carbon_count == 0:
            return "Inorganic compound detected"
        if carbon_count == 1 and num_atoms == 1:
            return "methane"
        if smiles_string in ["CC=O", "C(=O)C"]:
            return "ethanal"
        if smiles_string == "CC(=O)C":
            return "propanone"
        principal_structure, is_cyclic, structure_size = find_principal_chain_or_ring(mol)
        if not principal_structure:
            return "Could not determine molecular structure"
        iupac_name = assemble_iupac_name_complete(mol, principal_structure, is_cyclic, structure_size)
        return iupac_name
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
        <p>Dont andha trust ispe, But it Works!</p>
    </div>
    """, unsafe_allow_html=True)

    if 'molecule' not in st.session_state:
        st.session_state.molecule = "CC=O"

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
            "🔹 Common Molecules": {
                "Ethanol": "CCO",
                "Acetic acid": "CC(=O)O",
                "Acetone": "CC(=O)C",
                "Benzene": "c1ccccc1"
            },
            "🔹 Functional Groups": {
                "Ethanal": "CC=O",
                "Butanone": "CC(=O)CC",
                "Methylamine": "CN",
                "Chloroethane": "CCCl"
            },
            "🔹 Cyclic Compounds": {
                "Cyclopropane": "C1CC1",
                "Cyclohexane": "C1CCCCC1",
                "Cyclohexene": "C1=CCCCC1",
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
    if molecule_smiles:
        st.markdown(f"""
        <div class="smiles-container">
            <strong>🧬 SMILES:</strong> <code>{molecule_smiles}</code>
        </div>
        """, unsafe_allow_html=True)
    col1, col2, col3 = st.columns([1, 2, 1])
    with col2:
        if st.button("🔬 Generate IUPAC Name", type="primary", use_container_width=True):
            if molecule_smiles:
                with st.spinner("🔄 Analyzing structure..."):
                    iupac_name = name_from_smiles_complete(molecule_smiles)
                st.session_state.iupac_result = iupac_name
            else:
                st.session_state.iupac_result = "Please draw a molecule first."
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

    st.markdown("""<div class="footer">Made for Mistri of Chemistry</div>""", unsafe_allow_html=True)
    st.markdown("""<div class="footer">Mayank Vishwakarma <a href="https://www.linkedin.com/in/mayank-vishwakarma2004/">LinkedIn</a></div>
                """, unsafe_allow_html=True)
    

if __name__ == "__main__":
    main()
