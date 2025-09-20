
import streamlit as st
import pandas as pd
import requests
import time
import urllib.parse
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Draw
import re
import json
from datetime import datetime, timedelta
import hashlib
from typing import Dict, List, Optional, Tuple
import asyncio
import aiohttp
from concurrent.futures import ThreadPoolExecutor
import io
import base64


# Try to import streamlit_ketcher, provide fallback if not available
try:
    from streamlit_ketcher import st_ketcher
    KETCHER_AVAILABLE = True
except ImportError:
    KETCHER_AVAILABLE = False


# ==============================================================================
# PAGE CONFIGURATION 
# ==============================================================================

st.set_page_config(
    page_title="Advanced IUPAC Name Generator",
    page_icon="⚗️", 
    layout="wide",
    initial_sidebar_state="expanded"
)


# ==============================================================================
# ENHANCED CACHING SYSTEM WITH st.cache_data AND st.cache_resource
# ==============================================================================

@st.cache_data(ttl=3600, max_entries=1000, show_spinner=False)
def cached_iupac_lookup(smiles: str, api_source: str) -> Optional[str]:
    """Cached IUPAC name lookup with 1-hour TTL"""
    if api_source == "nci":
        return smiles_to_iupac_nci_enhanced(smiles)
    elif api_source == "pubchem":
        return smiles_to_iupac_pubchem_enhanced(smiles)
    elif api_source == "chemspider":
        return smiles_to_iupac_chemspider_enhanced(smiles)
    return None

@st.cache_data(ttl=86400, show_spinner=False)  # Cache for 24 hours
def cached_molecular_properties(smiles: str) -> Dict:
    """Cache molecular properties calculation"""
    return calculate_enhanced_molecular_properties(smiles)

@st.cache_resource
def load_example_molecules() -> Dict:
    """Cache example molecules data"""
    return {
        "🔹 Basic Alkanes": {
            "Methane": "C",
            "Ethane": "CC", 
            "Propane": "CCC",
            "Butane": "CCCC",
            "2-Methylpropane": "CC(C)C",
            "Pentane": "CCCCC",
            "Hexane": "CCCCCC"
        },
        "🔹 Alkenes & Alkynes": {
            "Ethene": "C=C",
            "Propene": "C=CC", 
            "But-1-ene": "C=CCC",
            "But-2-ene": "CC=CC",
            "Ethyne": "C#C",
            "Propyne": "C#CC"
        },
        "🔹 Functional Groups": {
            "Methanol": "CO",
            "Ethanol": "CCO",
            "Propan-1-ol": "CCCO", 
            "Propan-2-ol": "CC(O)C",
            "Ethanoic acid": "CC(=O)O",
            "Methanal": "C=O"
        },
        "🔹 Aromatics": {
            "Benzene": "c1ccccc1",
            "Toluene": "Cc1ccccc1",
            "Phenol": "Oc1ccccc1",
            "Benzoic acid": "c1ccc(cc1)C(=O)O"
        },
        "🔹 Complex Examples": {
            "Aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
            "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "Glucose": "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O",
            "Cholesterol": "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C"
        }
    }


# ==============================================================================
# ENHANCED SMILES VALIDATION AND PROCESSING
# ==============================================================================

@st.cache_data(show_spinner=False)
def validate_and_standardize_smiles(smiles: str) -> Tuple[bool, str, Optional[str]]:
    """Enhanced SMILES validation with standardization and error details"""
    try:
        if not smiles or smiles.strip() == "":
            return False, smiles, "Empty SMILES string"

        # Clean the SMILES
        clean_smiles = smiles.strip()

        # Try to parse with RDKit
        mol = Chem.MolFromSmiles(clean_smiles)
        if not mol:
            return False, clean_smiles, "Invalid SMILES syntax - RDKit cannot parse"

        # Try to sanitize
        try:
            Chem.SanitizeMol(mol)
        except Exception as e:
            return False, clean_smiles, f"Sanitization failed: {str(e)}"

        # Check for reasonable molecule size
        if mol.GetNumAtoms() == 0:
            return False, clean_smiles, "No atoms found in molecule"

        if mol.GetNumAtoms() > 200:
            return False, clean_smiles, "Molecule too large (>200 atoms)"

        # Check for carbon atoms (organic chemistry focus)
        carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        if carbon_count == 0:
            return False, clean_smiles, "No carbon atoms - inorganic compound"

        # Standardize the SMILES
        standardized_smiles = Chem.MolToSmiles(mol, canonical=True)

        return True, standardized_smiles, None

    except Exception as e:
        return False, smiles, f"Validation error: {str(e)}"


# ==============================================================================
# ADVANCED API MANAGEMENT WITH RATE LIMITING AND RETRY LOGIC
# ==============================================================================

class APIRateLimiter:
    """Advanced rate limiting with exponential backoff"""

    def __init__(self):
        if 'api_calls' not in st.session_state:
            st.session_state.api_calls = {}

    def can_make_request(self, api_name: str, max_calls_per_minute: int = 30) -> bool:
        """Check if we can make an API request"""
        now = datetime.now()
        if api_name not in st.session_state.api_calls:
            st.session_state.api_calls[api_name] = []

        # Remove calls older than 1 minute
        st.session_state.api_calls[api_name] = [
            call_time for call_time in st.session_state.api_calls[api_name]
            if now - call_time < timedelta(minutes=1)
        ]

        return len(st.session_state.api_calls[api_name]) < max_calls_per_minute

    def record_request(self, api_name: str):
        """Record that we made an API request"""
        if api_name not in st.session_state.api_calls:
            st.session_state.api_calls[api_name] = []
        st.session_state.api_calls[api_name].append(datetime.now())


rate_limiter = APIRateLimiter()


def make_api_request_with_retry(url: str, headers: Dict = None, data: Dict = None, 
                               max_retries: int = 3, timeout: int = 15) -> Optional[requests.Response]:
    """Make API request with exponential backoff retry logic"""
    for attempt in range(max_retries):
        try:
            if data:
                response = requests.post(url, headers=headers, data=data, timeout=timeout)
            else:
                response = requests.get(url, headers=headers, timeout=timeout)

            if response.status_code == 200:
                return response
            elif response.status_code == 429:  # Rate limited
                wait_time = 2 ** attempt  # Exponential backoff
                time.sleep(wait_time)
                continue
            else:
                return None

        except requests.exceptions.RequestException as e:
            if attempt == max_retries - 1:
                return None
            time.sleep(2 ** attempt)

    return None


# ==============================================================================
# ENHANCED API IMPLEMENTATIONS
# ==============================================================================

def encode_smiles_for_url(smiles: str) -> str:
    """Enhanced URL encoding for SMILES"""
    # Use urllib.parse.quote for proper URL encoding
    return urllib.parse.quote(smiles, safe='')


def smiles_to_iupac_nci_enhanced(smiles: str) -> Optional[str]:
    """Enhanced NCI API with better error handling and rate limiting"""
    try:
        if not rate_limiter.can_make_request("nci", 20):
            return None

        encoded_smiles = encode_smiles_for_url(smiles)
        base_url = "https://cactus.nci.nih.gov/chemical/structure"
        url = f"{base_url}/{encoded_smiles}/iupac_name"

        headers = {
            'User-Agent': 'Mozilla/5.0 (Advanced-IUPAC-Generator/2.0)',
            'Accept': 'text/plain, text/html, */*',
            'Accept-Language': 'en-US,en;q=0.9'
        }

        response = make_api_request_with_retry(url, headers=headers)

        if response:
            rate_limiter.record_request("nci")
            iupac_name = response.text.strip()

            # Enhanced validation
            if (iupac_name and 
                len(iupac_name) > 0 and 
                len(iupac_name) < 500 and  # Reasonable length check
                not any(error in iupac_name.lower() for error in 
                       ["error", "not found", "404", "<!doctype", "exception"])):
                return iupac_name

        return None

    except Exception as e:
        return None


def smiles_to_iupac_pubchem_enhanced(smiles: str) -> Optional[str]:
    """Enhanced PubChem API with better handling"""
    try:
        if not rate_limiter.can_make_request("pubchem", 25):
            return None

        # Step 1: Get CID from SMILES
        search_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/cids/JSON"

        data = {"smiles": smiles}
        headers = {
            'Content-Type': 'application/x-www-form-urlencoded',
            'User-Agent': 'Advanced-IUPAC-Generator/2.0'
        }

        response1 = make_api_request_with_retry(search_url, headers=headers, data=data)

        if response1:
            rate_limiter.record_request("pubchem")
            result1 = response1.json()

            if "IdentifierList" in result1 and "CID" in result1["IdentifierList"]:
                cid = result1["IdentifierList"]["CID"][0]

                # Step 2: Get IUPAC name from CID
                name_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IUPACName/JSON"
                response2 = make_api_request_with_retry(name_url, headers=headers)

                if response2:
                    result2 = response2.json()
                    if ("PropertyTable" in result2 and 
                        "Properties" in result2["PropertyTable"] and 
                        len(result2["PropertyTable"]["Properties"]) > 0):

                        properties = result2["PropertyTable"]["Properties"][0]
                        if "IUPACName" in properties:
                            return properties["IUPACName"]

        return None

    except Exception as e:
        return None


def smiles_to_iupac_chemspider_enhanced(smiles: str) -> Optional[str]:
    """Enhanced ChemSpider implementation (placeholder for API key implementation)"""
    # This would require API key setup for full implementation
    return None


def smiles_to_iupac_opsin(smiles: str) -> Optional[str]:
    """OPSIN (Open Parser for Systematic IUPAC nomenclature) API"""
    try:
        if not rate_limiter.can_make_request("opsin", 15):
            return None

        # Convert SMILES to IUPAC using OPSIN web service
        url = "https://opsin.ch.cam.ac.uk/opsin/smilesToIupac"
        data = {"smiles": smiles}

        response = make_api_request_with_retry(url, data=data)

        if response:
            rate_limiter.record_request("opsin")
            result = response.text.strip()
            if result and not result.startswith("Could not"):
                return result

        return None

    except Exception as e:
        return None


# ==============================================================================
# ENHANCED MOLECULAR PROPERTIES CALCULATION
# ==============================================================================

def calculate_enhanced_molecular_properties(smiles: str) -> Dict:
    """Calculate comprehensive molecular properties"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return {}

        properties = {
            # Basic properties
            "molecular_formula": Chem.rdMolDescriptors.CalcMolFormula(mol),
            "molecular_weight": round(Descriptors.MolWt(mol), 2),
            "exact_mass": round(Descriptors.ExactMolWt(mol), 4),

            # Atom counts
            "num_atoms": mol.GetNumAtoms(),
            "num_bonds": mol.GetNumBonds(),
            "num_rings": mol.GetRingInfo().NumRings(),

            # Element counts
            "carbon_atoms": sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6),
            "nitrogen_atoms": sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7),
            "oxygen_atoms": sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8),
            "hetero_atoms": Descriptors.NumHeteroatoms(mol),

            # Chemical properties
            "logp": round(Descriptors.MolLogP(mol), 2),
            "tpsa": round(Descriptors.TPSA(mol), 2),
            "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
            "h_bond_donors": Descriptors.NumHDonors(mol),
            "h_bond_acceptors": Descriptors.NumHAcceptors(mol),

            # Aromaticity
            "aromatic_rings": Descriptors.NumAromaticRings(mol),
            "aromatic_atoms": sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic()),
        }

        # Lipinski Rule of Five
        properties["lipinski_violations"] = sum([
            properties["molecular_weight"] > 500,
            properties["logp"] > 5,
            properties["h_bond_donors"] > 5,
            properties["h_bond_acceptors"] > 10
        ])

        return properties

    except Exception as e:
        return {"error": str(e)}


# ==============================================================================
# ADVANCED IUPAC NAME GENERATION WITH MULTIPLE SOURCES
# ==============================================================================

def get_advanced_iupac_name(smiles: str) -> Dict:
    """Get IUPAC name using multiple sources with detailed results"""

    if not smiles or smiles.strip() == "":
        return {"error": "Please draw a molecule first.", "source": None}

    # Validate and standardize SMILES
    is_valid, clean_smiles, error_msg = validate_and_standardize_smiles(smiles)
    if not is_valid:
        return {"error": f"Invalid SMILES: {error_msg}", "source": None}

    # Try multiple APIs in order of reliability
    api_sources = [
        ("NCI Chemical Identifier Resolver", "nci"),
        ("PubChem Database", "pubchem"),
        ("OPSIN Parser", "opsin")
    ]

    results = {}

    for source_name, source_key in api_sources:
        with st.status(f"🔍 Querying {source_name}...", expanded=False) as status:
            try:
                iupac_name = cached_iupac_lookup(clean_smiles, source_key)

                if iupac_name:
                    status.update(label=f"✅ Found result from {source_name}", state="complete")
                    return {
                        "iupac_name": iupac_name,
                        "source": source_name,
                        "smiles_used": clean_smiles,
                        "confidence": "high"
                    }
                else:
                    status.update(label=f"❌ No result from {source_name}", state="error")

            except Exception as e:
                status.update(label=f"❌ Error from {source_name}: {str(e)}", state="error")

            time.sleep(0.5)  # Brief pause between API calls

    # Fallback to enhanced naming if all APIs fail
    try:
        mol = Chem.MolFromSmiles(clean_smiles)
        fallback_name = get_advanced_fallback_name(mol)
        return {
            "iupac_name": fallback_name,
            "source": "Enhanced Fallback System",
            "smiles_used": clean_smiles,
            "confidence": "estimated"
        }
    except Exception as e:
        return {"error": f"All naming methods failed: {str(e)}", "source": None}


def get_advanced_fallback_name(mol) -> str:
    """Advanced fallback naming with better chemical knowledge"""
    try:
        # Get molecular info
        carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        total_atoms = mol.GetNumAtoms()

        # Enhanced alkane names up to C20
        alkane_names = {
            1: "methane", 2: "ethane", 3: "propane", 4: "butane", 5: "pentane",
            6: "hexane", 7: "heptane", 8: "octane", 9: "nonane", 10: "decane",
            11: "undecane", 12: "dodecane", 13: "tridecane", 14: "tetradecane",
            15: "pentadecane", 16: "hexadecane", 17: "heptadecane", 18: "octadecane",
            19: "nonadecane", 20: "icosane"
        }

        # Analyze molecular structure
        ring_info = mol.GetRingInfo()
        is_cyclic = ring_info.NumRings() > 0
        num_rings = ring_info.NumRings()

        # Bond analysis
        bonds = mol.GetBonds()
        has_double_bond = any(bond.GetBondType() == Chem.rdchem.BondType.DOUBLE for bond in bonds)
        has_triple_bond = any(bond.GetBondType() == Chem.rdchem.BondType.TRIPLE for bond in bonds)

        # Functional group detection
        atoms = mol.GetAtoms()
        has_alcohol = any(atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() >= 1 for atom in atoms)
        has_carbonyl = any(
            atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 0 and
            any(bond.GetBondType() == Chem.rdchem.BondType.DOUBLE for bond in atom.GetBonds())
            for atom in atoms
        )
        has_carboxyl = "C(=O)O" in Chem.MolToSmiles(mol)
        has_amino = any(atom.GetAtomicNum() == 7 for atom in atoms)

        # Aromatic detection
        is_aromatic = any(atom.GetIsAromatic() for atom in atoms)
        aromatic_atoms = sum(1 for atom in atoms if atom.GetIsAromatic())

        # Generate name based on structure
        if is_aromatic and carbon_count == 6 and aromatic_atoms == 6:
            if has_alcohol:
                return "phenol (or substituted phenol)"
            elif has_carboxyl:
                return "benzoic acid (or substituted benzoic acid)"
            elif carbon_count == 7:  # Toluene and derivatives
                return "toluene (or substituted toluene)"
            else:
                return "benzene (or substituted benzene)"

        elif is_cyclic:
            # Cyclic compounds
            ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
            primary_ring_size = ring_sizes[0] if ring_sizes else carbon_count

            if primary_ring_size in [3, 4, 5, 6, 7, 8]:
                cycle_names = {3: "cyclopropane", 4: "cyclobutane", 5: "cyclopentane", 
                              6: "cyclohexane", 7: "cycloheptane", 8: "cyclooctane"}
                base_name = cycle_names.get(primary_ring_size, f"cyclo-C{primary_ring_size}")

                # Modify for unsaturation
                if has_triple_bond:
                    return f"{base_name.replace('ane', 'yne')} derivative"
                elif has_double_bond:
                    return f"{base_name.replace('ane', 'ene')} derivative"
                else:
                    return f"{base_name} derivative"
            else:
                return f"cyclic compound (C{carbon_count})"

        elif carbon_count in alkane_names:
            # Linear molecules
            base_name = alkane_names[carbon_count]

            if has_carboxyl:
                if carbon_count == 1:
                    return "methanoic acid (formic acid)"
                else:
                    return f"{base_name.replace('ane', 'anoic')} acid"
            elif has_carbonyl and not has_carboxyl:
                if carbon_count == 1:
                    return "methanal (formaldehyde)"
                elif carbon_count == 2:
                    return "ethanal (acetaldehyde)"
                else:
                    return f"{base_name.replace('ane', 'anal')} or {base_name.replace('ane', 'anone')}"
            elif has_alcohol:
                if carbon_count == 1:
                    return "methanol"
                else:
                    return f"{base_name.replace('ane', 'anol')}"
            elif has_triple_bond:
                if carbon_count == 2:
                    return "ethyne (acetylene)"
                else:
                    return f"{base_name.replace('ane', 'yne')}"
            elif has_double_bond:
                if carbon_count == 2:
                    return "ethene (ethylene)"
                else:
                    return f"{base_name.replace('ane', 'ene')}"
            else:
                return base_name

        else:
            # Complex molecules
            structure_desc = []
            if is_aromatic:
                structure_desc.append("aromatic")
            if is_cyclic:
                structure_desc.append("cyclic")
            if has_multiple_rings := num_rings > 1:
                structure_desc.append("polycyclic")

            desc = " ".join(structure_desc) if structure_desc else "complex"
            return f"{desc.capitalize()} organic compound (C{carbon_count})"

    except Exception as e:
        return "Complex organic molecule - structure analysis failed"


# ==============================================================================
# ENHANCED UI WITH MODERN STREAMLIT FEATURES
# ==============================================================================

def inject_enhanced_css():
    """Enhanced CSS with modern design"""
    st.markdown("""
    <style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap');

    .main .block-container {
        padding: 1rem 2rem;
        max-width: 1600px;
        font-family: 'Inter', sans-serif;
    }

    .main-header {
        text-align: center;
        padding: 2.5rem 0;
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        margin: -1rem -2rem 2rem -2rem;
        color: white;
        border-radius: 0 0 20px 20px;
        box-shadow: 0 8px 32px rgba(102, 126, 234, 0.15);
    }

    .main-header h1 {
        margin: 0;
        font-size: 2.5rem;
        font-weight: 700;
        text-shadow: 0 2px 4px rgba(0,0,0,0.1);
        letter-spacing: -0.02em;
    }

    .main-header p {
        margin: 0.8rem 0 0 0;
        font-size: 1.1rem;
        opacity: 0.95;
        font-weight: 400;
    }

    .result-success {
        background: linear-gradient(135deg, #d4edda 0%, #c3e6cb 100%);
        border: 2px solid #b8dacd;
        color: #155724;
        padding: 2.5rem;
        border-radius: 16px;
        text-align: center;
        margin: 2rem 0;
        box-shadow: 0 8px 24px rgba(212, 237, 218, 0.2);
    }

    .result-success .iupac-name {
        font-size: 1.8rem;
        font-weight: 700;
        font-family: 'Georgia', serif;
        color: #0f5132;
        line-height: 1.4;
        margin: 1rem 0;
        word-break: break-word;
        background: rgba(255,255,255,0.7);
        padding: 1rem;
        border-radius: 8px;
    }

    .properties-grid {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
        gap: 1rem;
        margin: 1.5rem 0;
    }

    .property-card {
        background: #f8fafc;
        padding: 1rem;
        border-radius: 8px;
        border-left: 4px solid #667eea;
        text-align: center;
    }

    .property-value {
        font-size: 1.2rem;
        font-weight: 600;
        color: #2d3748;
    }

    .property-label {
        font-size: 0.9rem;
        color: #64748b;
        margin-top: 0.25rem;
    }

    .status-badge {
        display: inline-block;
        padding: 0.25rem 0.75rem;
        border-radius: 20px;
        font-size: 0.85rem;
        font-weight: 500;
        margin: 0.25rem;
    }

    .status-high { background: #d4edda; color: #155724; }
    .status-medium { background: #fff3cd; color: #856404; }
    .status-low { background: #f8d7da; color: #721c24; }

    .history-item {
        background: #f8fafc;
        padding: 1rem;
        border-radius: 8px;
        margin: 0.5rem 0;
        border-left: 4px solid #e2e8f0;
    }

    .smiles-code {
        font-family: 'JetBrains Mono', monospace;
        background: #2d3748;
        color: #e2e8f0;
        padding: 0.75rem;
        border-radius: 6px;
        font-size: 0.9rem;
        word-break: break-all;
    }
    </style>
    """, unsafe_allow_html=True)


# ==============================================================================
# HISTORY AND SESSION MANAGEMENT
# ==============================================================================

def initialize_session_state():
    """Initialize session state variables"""
    if 'molecule_history' not in st.session_state:
        st.session_state.molecule_history = []
    if 'current_molecule' not in st.session_state:
        st.session_state.current_molecule = "CCO"
    if 'last_result' not in st.session_state:
        st.session_state.last_result = None


def add_to_history(smiles: str, iupac_name: str, source: str):
    """Add successful result to history"""
    timestamp = datetime.now().strftime("%H:%M:%S")
    entry = {
        "timestamp": timestamp,
        "smiles": smiles,
        "iupac_name": iupac_name,
        "source": source
    }

    # Avoid duplicates
    if not any(h["smiles"] == smiles for h in st.session_state.molecule_history):
        st.session_state.molecule_history.insert(0, entry)  # Add to beginning

    # Keep only last 20 entries
    if len(st.session_state.molecule_history) > 20:
        st.session_state.molecule_history = st.session_state.molecule_history[:20]


# ==============================================================================
# MAIN APPLICATION
# ==============================================================================

def main():
    """Enhanced main application"""
    inject_enhanced_css()
    initialize_session_state()

    # Header
    st.markdown("""
    <div class="main-header">
        <h1>⚗️ Advanced IUPAC Name Generator</h1>
        <p>Professional-grade chemical nomenclature with multiple data sources</p>
    </div>
    """, unsafe_allow_html=True)

    # Sidebar for settings and history
    with st.sidebar:
        st.header("🔧 Settings & History")

        # Settings
        with st.expander("⚙️ API Settings", expanded=False):
            st.info("API rate limits automatically managed")
            col1, col2 = st.columns(2)
            with col1:
                if 'api_calls' in st.session_state:
                    total_calls = sum(len(calls) for calls in st.session_state.api_calls.values())
                    st.metric("API Calls", total_calls)
            with col2:
                cache_info = st.cache_data.clear
                if st.button("Clear Cache", help="Clear cached results"):
                    st.cache_data.clear()
                    st.success("Cache cleared!")

        # History
        with st.expander("📚 Recent History", expanded=True):
            if st.session_state.molecule_history:
                for i, entry in enumerate(st.session_state.molecule_history[:5]):  # Show last 5
                    with st.container():
                        st.markdown(f"**{entry['timestamp']}**")
                        st.markdown(f"`{entry['smiles']}`")
                        st.markdown(f"*{entry['iupac_name'][:50]}{'...' if len(entry['iupac_name']) > 50 else ''}*")
                        if st.button("Load", key=f"load_{i}", help="Load this molecule"):
                            st.session_state.current_molecule = entry['smiles']
                            st.rerun()
                        st.divider()
            else:
                st.info("No history yet")

        # Export functionality
        if st.session_state.molecule_history:
            if st.button("📥 Export History"):
                df = pd.DataFrame(st.session_state.molecule_history)
                csv = df.to_csv(index=False)
                st.download_button(
                    "Download CSV",
                    csv,
                    "iupac_history.csv",
                    "text/csv"
                )

    # Main content
    col1, col2 = st.columns([2.5, 1.5], gap="large")

    # Molecule Editor
    with col1:
        st.subheader("🖊️ Molecule Editor")

        if KETCHER_AVAILABLE:
            try:
                molecule_smiles = st_ketcher(
                    value=st.session_state.current_molecule,
                    height=450,
                    key="advanced_molecule_editor"
                )

                if molecule_smiles and molecule_smiles != st.session_state.current_molecule:
                    st.session_state.current_molecule = molecule_smiles
                    st.session_state.last_result = None

            except Exception as e:
                st.error(f"Ketcher Error: {str(e)}")
                molecule_smiles = st.text_input(
                    "SMILES Input (Ketcher unavailable)",
                    value=st.session_state.current_molecule,
                    placeholder="e.g., CCO for ethanol"
                )
        else:
            st.warning("⚠️ Ketcher editor not available. Install: `pip install streamlit-ketcher`")
            molecule_smiles = st.text_input(
                "Enter SMILES string:",
                value=st.session_state.current_molecule,
                placeholder="e.g., CCO for ethanol, C#C for ethyne"
            )

    # Controls and Examples
    with col2:
        st.subheader("⚡ Quick Actions")

        
        if st.button("🗑️ Clear", use_container_width=True):
            st.session_state.current_molecule = ""
            st.session_state.last_result = None
            st.rerun()


        st.markdown("---")
        st.subheader("📚 Example Molecules")

        examples = load_example_molecules()
        for category, molecules in examples.items():
            with st.expander(category):
                for name, smiles in molecules.items():
                    if st.button(name, key=f"ex_{category}_{name}", use_container_width=True):
                        st.session_state.current_molecule = smiles
                        st.session_state.last_result = None
                        st.rerun()

    # Current molecule display
    current_smiles = st.session_state.current_molecule
    if current_smiles:
        is_valid, clean_smiles, error_msg = validate_and_standardize_smiles(current_smiles)

        col_x, col_y = st.columns([3, 1])
        with col_x:
            if is_valid:
                st.markdown(f'<div class="smiles-code">✅ Valid SMILES: {clean_smiles}</div>', unsafe_allow_html=True)
            else:
                st.markdown(f'<div class="smiles-code">❌ Invalid SMILES: {current_smiles}</div>', unsafe_allow_html=True)
                st.error(f"Validation Error: {error_msg}")

        with col_y:
            if st.button("🔬 Generate IUPAC Name", type="primary", use_container_width=True, disabled=not is_valid):
                if is_valid:
                    result = get_advanced_iupac_name(clean_smiles)
                    st.session_state.last_result = result

                    if "iupac_name" in result:
                        add_to_history(clean_smiles, result["iupac_name"], result["source"])

                    st.rerun()

    # Results display
    if st.session_state.last_result:
        result = st.session_state.last_result

        if "error" in result:
            st.error(f"❌ {result['error']}")
        else:
            # Success display
            confidence_class = f"status-{result['confidence']}" if result['confidence'] in ['high', 'medium', 'low'] else "status-medium"

            st.markdown(f"""
            <div class="result-success">
                <h3>🎉 IUPAC Name Generated Successfully</h3>
                <div class="iupac-name">{result['iupac_name']}</div>
                <div style="margin-top: 1rem;">
                    <span class="status-badge {confidence_class}">
                        Confidence: {result['confidence'].title()}
                    </span>
                    <span class="status-badge status-medium">
                        Source: {result['source']}
                    </span>
                </div>
            </div>
            """, unsafe_allow_html=True)

            # Molecular properties
            if current_smiles:
                try:
                    properties = cached_molecular_properties(result['smiles_used'])

                    if properties and 'error' not in properties:
                        st.subheader("📊 Molecular Properties")

                        # Key properties in columns
                        prop_cols = st.columns(4)
                        key_props = [
                            ("Formula", properties.get("molecular_formula", "N/A")),
                            ("MW (g/mol)", properties.get("molecular_weight", "N/A")),
                            ("LogP", properties.get("logp", "N/A")),
                            ("Rings", properties.get("num_rings", "N/A"))
                        ]

                        for i, (label, value) in enumerate(key_props):
                            with prop_cols[i]:
                                st.metric(label, value)

                        # Detailed properties in expander
                        with st.expander("🔍 Detailed Properties"):
                            detail_cols = st.columns(3)

                            with detail_cols[0]:
                                st.markdown("**Structure**")
                                st.write(f"Atoms: {properties.get('num_atoms', 'N/A')}")
                                st.write(f"Bonds: {properties.get('num_bonds', 'N/A')}")
                                st.write(f"Carbon atoms: {properties.get('carbon_atoms', 'N/A')}")
                                st.write(f"Hetero atoms: {properties.get('hetero_atoms', 'N/A')}")

                            with detail_cols[1]:
                                st.markdown("**Chemical Properties**")
                                st.write(f"TPSA: {properties.get('tpsa', 'N/A')}")
                                st.write(f"Rotatable bonds: {properties.get('rotatable_bonds', 'N/A')}")
                                st.write(f"H-bond donors: {properties.get('h_bond_donors', 'N/A')}")
                                st.write(f"H-bond acceptors: {properties.get('h_bond_acceptors', 'N/A')}")

                            with detail_cols[2]:
                                st.markdown("**Drug-likeness**")
                                violations = properties.get('lipinski_violations', 0)
                                if violations == 0:
                                    st.success(f"Lipinski violations: {violations}")
                                elif violations <= 1:
                                    st.warning(f"Lipinski violations: {violations}")
                                else:
                                    st.error(f"Lipinski violations: {violations}")

                                st.write(f"Aromatic rings: {properties.get('aromatic_rings', 'N/A')}")
                                st.write(f"Exact mass: {properties.get('exact_mass', 'N/A')}")

                except Exception as e:
                    st.error(f"Error calculating properties: {e}")

    # Footer
    st.markdown("---")
    col_foot1, col_foot2, col_foot3 = st.columns(3)

    with col_foot1:
        st.markdown("**🧪 Advanced IUPAC Generator**")
        st.caption("Professional chemical nomenclature tool")

    with col_foot2:
        st.markdown("**📊 Data Sources**")
        st.caption("NCI, PubChem, OPSIN, Enhanced Fallback")

    with col_foot3:
        st.markdown("**👨‍💻 Developer**")
        st.caption("Developed by Mayank Vishwakarma")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        st.error(f"Application Error: {str(e)}")
        st.exception(e)
        if st.button("🔄 Restart Application"):
            st.rerun()
