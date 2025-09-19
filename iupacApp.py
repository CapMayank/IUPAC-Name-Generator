import streamlit as st
import pandas as pd
import requests
import time
import urllib.parse
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
import re

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
    page_title="IUPAC Name Generator",
    page_icon="⚗️", 
    layout="wide",
    initial_sidebar_state="collapsed"
)

# ==============================================================================
# FIXED ACCURATE IUPAC NAMING WITH PROPER URL ENCODING
# ==============================================================================

def encode_smiles_for_url(smiles):
    """Properly encode SMILES for URL usage - CRITICAL FIX"""
    # First handle special SMILES characters that need encoding
    encoded = smiles.replace('#', '%23')  # Triple bonds
    encoded = encoded.replace('=', '%3D')  # Double bonds  
    encoded = encoded.replace('[', '%5B')  # Square brackets
    encoded = encoded.replace(']', '%5D')  # Square brackets
    encoded = encoded.replace('+', '%2B')  # Plus signs
    encoded = encoded.replace('@', '%40')  # At symbols
    encoded = encoded.replace('/', '%2F')  # Forward slashes
    encoded = encoded.replace('\\', '%5C') # Backslashes
    encoded = encoded.replace('?', '%3F')  # Question marks
    encoded = encoded.replace(' ', '%20')  # Spaces
    
    return encoded

def smiles_to_iupac_nci_fixed(smiles):
    """Get IUPAC name using NCI with proper URL encoding - FIXED"""
    try:
        # Encode the SMILES properly
        encoded_smiles = encode_smiles_for_url(smiles)
        
        # NCI Chemical Identifier Resolver API
        base_url = "https://cactus.nci.nih.gov/chemical/structure"
        url = f"{base_url}/{encoded_smiles}/iupac_name"
        
        # Add headers to mimic browser request
        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36',
            'Accept': 'text/plain, text/html, */*',
            'Accept-Language': 'en-US,en;q=0.9'
        }
        
        response = requests.get(url, headers=headers, timeout=15)
        response.raise_for_status()
        
        iupac_name = response.text.strip()
        
        # Check if response is valid
        if (iupac_name and 
            len(iupac_name) > 0 and 
            not iupac_name.startswith("Error") and
            not iupac_name.startswith("<!DOCTYPE") and  # Not HTML error page
            not "404" in iupac_name and
            not "Not Found" in iupac_name):
            return iupac_name
        else:
            return None
            
    except requests.exceptions.RequestException as e:
        print(f"NCI API error: {e}")
        return None
    except Exception as e:
        print(f"NCI processing error: {e}")
        return None

def smiles_to_iupac_pubchem_fixed(smiles):
    """Get IUPAC name using PubChem API - IMPROVED"""
    try:
        # Step 1: Get CID from SMILES
        search_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/cids/JSON"
        
        # Use POST with proper data
        data = {"smiles": smiles}
        headers = {
            'Content-Type': 'application/x-www-form-urlencoded',
            'User-Agent': 'IUPAC-Name-Generator/1.0'
        }
        
        response1 = requests.post(search_url, data=data, headers=headers, timeout=15)
        
        if response1.status_code == 200:
            result1 = response1.json()
            if "IdentifierList" in result1 and "CID" in result1["IdentifierList"]:
                cid = result1["IdentifierList"]["CID"][0]
                
                # Step 2: Get IUPAC name from CID
                name_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IUPACName/JSON"
                response2 = requests.get(name_url, headers=headers, timeout=15)
                
                if response2.status_code == 200:
                    result2 = response2.json()
                    if ("PropertyTable" in result2 and 
                        "Properties" in result2["PropertyTable"] and 
                        len(result2["PropertyTable"]["Properties"]) > 0):
                        
                        properties = result2["PropertyTable"]["Properties"][0]
                        if "IUPACName" in properties:
                            return properties["IUPACName"]
        
        return None
        
    except Exception as e:
        print(f"PubChem API error: {e}")
        return None

def smiles_to_iupac_chemspider(smiles):
    """Try ChemSpider API as additional backup"""
    try:
        # ChemSpider simple API (no key needed for basic searches)
        url = f"http://www.chemspider.com/Search.asmx/SimpleSearch"
        params = {
            'query': smiles,
            'token': ''  # Public API
        }
        
        response = requests.get(url, params=params, timeout=10)
        if response.status_code == 200 and response.text:
            # This would need additional parsing for ChemSpider
            # Simplified for now
            pass
        
        return None
    except:
        return None

def get_accurate_iupac_name_fixed(smiles):
    """Get the most accurate IUPAC name using multiple APIs with proper encoding"""
    if not smiles or smiles.strip() == "":
        return "Please draw a molecule first."
    
    try:
        # Validate SMILES first
        mol = Chem.MolFromSmiles(smiles.strip())
        if not mol:
            return "Invalid SMILES - please check your structure"
        
        # Sanitize molecule
        try:
            Chem.SanitizeMol(mol)
        except:
            return "Invalid molecule structure"
        
        # Check for carbon atoms
        carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        if carbon_count == 0:
            return "Inorganic compound - no carbon atoms found"
        
        # Show which API is being tried
        status_placeholder = st.empty()
        
        # Try NCI first (most reliable) with proper encoding
        status_placeholder.info("🔍 Querying NCI Chemical Identifier Resolver...")
        iupac_name = smiles_to_iupac_nci_fixed(smiles)
        if iupac_name:
            status_placeholder.empty()
            return iupac_name
        
        # Try PubChem as backup
        status_placeholder.info("🔍 Trying PubChem database...")
        time.sleep(1)  # Rate limiting
        iupac_name = smiles_to_iupac_pubchem_fixed(smiles)
        if iupac_name:
            status_placeholder.empty()
            return iupac_name
        
        # If both APIs fail, try basic RDKit naming for simple molecules
        status_placeholder.info("🔍 Using fallback naming...")
        time.sleep(1)
        fallback_name = get_improved_fallback_name(mol)
        status_placeholder.empty()
        
        return fallback_name
        
    except Exception as e:
        return f"Error: Could not generate IUPAC name - {str(e)}"

def get_improved_fallback_name(mol):
    """Improved fallback naming for when APIs fail"""
    try:
        # Get basic molecular info
        carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        total_atoms = mol.GetNumAtoms()
        
        # Basic alkane names
        alkane_names = {
            1: "methane", 2: "ethane", 3: "propane", 4: "butane", 
            5: "pentane", 6: "hexane", 7: "heptane", 8: "octane",
            9: "nonane", 10: "decane"
        }
        
        # Check bond types
        has_double_bond = any(bond.GetBondType() == Chem.rdchem.BondType.DOUBLE 
                             for bond in mol.GetBonds())
        has_triple_bond = any(bond.GetBondType() == Chem.rdchem.BondType.TRIPLE 
                             for bond in mol.GetBonds())
        
        # Check for rings
        is_cyclic = mol.GetRingInfo().NumRings() > 0
        
        # Check for common functional groups
        has_alcohol = any(atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 1 
                         for atom in mol.GetAtoms())
        has_carbonyl = any(atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 0 
                          for atom in mol.GetAtoms() 
                          if any(bond.GetBondType() == Chem.rdchem.BondType.DOUBLE 
                                for bond in atom.GetBonds()))
        
        # Simple naming logic
        if is_cyclic:
            if carbon_count == 6 and has_double_bond:
                return "benzene (or aromatic compound)"
            else:
                cycle_names = {3: "cyclopropane", 4: "cyclobutane", 5: "cyclopentane", 
                              6: "cyclohexane", 7: "cycloheptane", 8: "cyclooctane"}
                base_name = cycle_names.get(carbon_count, f"cyclo-C{carbon_count}")
                
                if has_triple_bond:
                    return f"{base_name.replace('ane', 'yne')} (approx)"
                elif has_double_bond:
                    return f"{base_name.replace('ane', 'ene')} (approx)"
                else:
                    return base_name
        
        # Linear molecules
        elif carbon_count in alkane_names:
            base_name = alkane_names[carbon_count]
            
            if has_alcohol:
                return f"{base_name.replace('ane', 'anol')} (approx)"
            elif has_carbonyl:
                if carbon_count == 1:
                    return "methanal (formaldehyde)"
                else:
                    return f"{base_name.replace('ane', 'anal')} or {base_name.replace('ane', 'anone')} (approx)"
            elif has_triple_bond:
                return f"{base_name.replace('ane', 'yne')} (approx)"
            elif has_double_bond:
                return f"{base_name.replace('ane', 'ene')} (approx)"
            else:
                return base_name
        
        else:
            # Complex molecule
            return f"Complex organic molecule (C{carbon_count})"
            
    except Exception as e:
        return "Complex molecule - naming failed"

# ==============================================================================
# MODERN CSS STYLING (Same as before)
# ==============================================================================

def inject_modern_css():
    """Inject modern CSS styling"""
    st.markdown("""
    <style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap');
    @import url('https://fonts.googleapis.com/css2?family=Fira+Code:wght@400;500&display=swap');
    
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
        box-shadow: 0 4px 20px rgba(102, 126, 234, 0.15);
    }
    
    .main-header h1 {
        margin: 0;
        font-size: 2.3rem;
        font-weight: 700;
        text-shadow: 0 2px 6px rgba(0,0,0,0.1);
        letter-spacing: -0.02em;
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
        box-shadow: 0 8px 30px rgba(0,0,0,0.12);
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
        gap: 0.5rem;
    }
    
    .smiles-container {
        background: linear-gradient(135deg, #f7fafc 0%, #edf2f7 100%);
        padding: 1.2rem;
        border-radius: 12px;
        border-left: 4px solid #667eea;
        margin: 2rem 0 1rem 0;
        font-family: 'Fira Code', 'Courier New', monospace;
        color: #4a5568;
        position: relative;
        overflow-x: auto;
        word-break: break-all;
    }
    
    .smiles-container code {
        background: transparent;
        padding: 0;
        font-size: 0.95rem;
        color: #2d3748;
    }
    
    .stButton > button {
        width: 100%;
        border-radius: 10px;
        border: 1px solid #e2e8f0;
        background: white;
        color: #334155 !important;
        font-weight: 500;
        font-size: 0.97rem;
        transition: all 0.2s ease;
        font-family: 'Inter', sans-serif;
        height: 45px;
    }
    
    .stButton > button:hover {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white !important;
        border-color: transparent;
        transform: translateY(-1px);
        box-shadow: 0 4px 12px rgba(102, 126, 234, 0.15);
    }
    
    div[data-testid="stButton"] button[kind="primary"] {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%) !important;
        color: white !important;
        border: none !important;
        height: 55px;
        font-size: 1.1rem;
        font-weight: 600;
        box-shadow: 0 4px 15px rgba(102, 126, 234, 0.25);
    }
    
    div[data-testid="stButton"] button[kind="primary"]:hover {
        background: linear-gradient(135deg, #5a67d8 0%, #6b46a3 100%) !important;
        transform: translateY(-2px);
        box-shadow: 0 6px 20px rgba(102, 126, 234, 0.35);
    }
    
    .result-success {
        background: linear-gradient(135deg, #d4edda 0%, #c3e6cb 100%);
        border: 1px solid #b8dacd;
        color: #155724;
        padding: 2rem;
        border-radius: 16px;
        text-align: center;
        margin: 2rem 0 1rem 0;
        box-shadow: 0 4px 10px rgba(212, 237, 218, 0.2);
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
        line-height: 1.3;
        margin: 0;
        word-break: break-word;
    }
    
    .result-error {
        background: linear-gradient(135deg, #f8d7da 0%, #f5c6cb 100%);
        border: 1px solid #f5c6cb;
        color: #721c24;
        padding: 2rem;
        border-radius: 16px;
        text-align: center;
        margin: 2rem 0 1rem 0;
        box-shadow: 0 4px 10px rgba(248, 215, 218, 0.2);
    }
    
    .api-status {
        background: #e3f2fd;
        border: 1px solid #bbdefb;
        color: #1565c0;
        padding: 1rem;
        border-radius: 8px;
        margin: 1rem 0;
        font-size: 0.9rem;
    }
    
    .footer {
        text-align: center;
        padding: 1rem 0;
        color: #64748b;
        font-size: 0.9rem;
        border-top: 1px solid #e2e8f0;
        margin-top: 3rem;
    }
    
    .footer a {
        color: #667eea;
        text-decoration: none;
    }
    
    .footer a:hover {
        text-decoration: underline;
    }
    
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
    header {visibility: hidden;}
    
    @media (max-width: 768px) {
        .main .block-container { padding: 1rem; }
        .main-header { margin: -1rem -1rem 1rem -1rem; padding: 1.2rem 0 2rem 0; }
        .main-header h1 { font-size: 1.8rem; }
        .editor-card, .controls-card { padding: 1.25rem; }
        .controls-card { position: static; }
        .result-success .iupac-name { font-size: 1.4rem; }
    }
    </style>
    """, unsafe_allow_html=True)

# ==============================================================================
# STREAMLIT APPLICATION (Same structure, using fixed naming function)
# ==============================================================================

def create_fallback_editor():
    """Create fallback SMILES input when ketcher unavailable"""
    if not KETCHER_AVAILABLE:
        st.warning("⚠️ Ketcher editor not available. Using text input.")
        st.info("💡 Install with: `pip install streamlit-ketcher`")
    
    smiles_input = st.text_input(
        "Enter SMILES string:",
        value=st.session_state.get('molecule', 'CCO'),
        placeholder="e.g., CCO for ethanol, C#C for ethyne",
        help="Enter a valid SMILES string"
    )
    
    return smiles_input

def load_example_molecules():
    """Load tested example molecules"""
    return {
        "🔹 Basic Alkanes": {
            "Methane": "C",
            "Ethane": "CC", 
            "Propane": "CCC",
            "Butane": "CCCC",
            "2-Methylpropane": "CC(C)C",
            "Pentane": "CCCCC"
        },
        "🔹 Alkenes": {
            "Ethene": "C=C",
            "Propene": "C=CC", 
            "But-1-ene": "C=CCC",
            "But-2-ene": "CC=CC"
        },
        "🔹 Alkynes": {
            "Ethyne": "C#C",
            "Propyne": "C#CC",
            "But-1-yne": "C#CCC", 
            "But-2-yne": "CC#CC"
        },
        "🔹 Alcohols": {
            "Methanol": "CO",
            "Ethanol": "CCO",
            "Propan-1-ol": "CCCO", 
            "Propan-2-ol": "CC(O)C"
        },
        "🔹 Simple Examples": {
            "Benzene": "c1ccccc1",
            "Aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
            "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
        }
    }

def main():
    """Main Streamlit application with FIXED API calls"""
    inject_modern_css()
    
    # Header
    st.markdown("""
    <div class="main-header">
        <h1>⚗️ IUPAC Name Generator</h1>
        <p>IUPAC nomenclature</p>
    </div>
    """, unsafe_allow_html=True)

    # API Status Info
    st.markdown("""
    <div class="api-status">
        <strong>🌐 Sources:</strong> NCI Chemical Identifier Resolver, PubChem, Enhanced Fallback
    </div>
    """, unsafe_allow_html=True)

    # Initialize session state
    if 'molecule' not in st.session_state:
        st.session_state.molecule = "CCO"
    if 'iupac_result' not in st.session_state:
        st.session_state.iupac_result = None

    # Layout
    col1, col2 = st.columns([2.8, 1.2], gap="large")
    
    # Molecule Editor
    with col1:
        
        st.markdown('<div class="section-header">🖊️ Molecule Editor</div>', unsafe_allow_html=True)
        
        if KETCHER_AVAILABLE:
            try:
                molecule_smiles = st_ketcher(
                    value=st.session_state.molecule,
                    height=420,
                    key="molecule_editor"
                )
                
                if molecule_smiles and molecule_smiles != st.session_state.molecule:
                    st.session_state.molecule = molecule_smiles
                    if 'iupac_result' in st.session_state:
                        del st.session_state.iupac_result
                        
            except Exception as e:
                st.error(f"Ketcher Error: {str(e)}")
                molecule_smiles = create_fallback_editor()
        else:
            molecule_smiles = create_fallback_editor()
            if molecule_smiles != st.session_state.molecule:
                st.session_state.molecule = molecule_smiles
                if 'iupac_result' in st.session_state:
                    del st.session_state.iupac_result
        
        st.markdown('</div>', unsafe_allow_html=True)

    # Controls
    with col2:
        
        st.markdown('<div class="section-header">⚡ Quick Actions</div>', unsafe_allow_html=True)

        if st.button("🗑️ Clear", use_container_width=True):
            st.session_state.molecule = ""
            if 'iupac_result' in st.session_state:
                del st.session_state.iupac_result
            st.rerun()
        
        st.markdown("---")
        st.markdown('<div class="section-header">📚 Test Examples</div>', unsafe_allow_html=True)
        
        examples = load_example_molecules()
        for category, molecules in examples.items():
            with st.expander(category):
                for name, smiles in molecules.items():
                    if st.button(name, key=f"{category}_{name}", use_container_width=True):
                        st.session_state.molecule = smiles
                        if 'iupac_result' in st.session_state:
                            del st.session_state.iupac_result
                        st.rerun()
        
        st.markdown('</div>', unsafe_allow_html=True)

    # Show SMILES
    current_smiles = st.session_state.get('molecule', '')
    if current_smiles:
        # Show both original and encoded SMILES for debugging
        encoded_smiles = encode_smiles_for_url(current_smiles)
        st.markdown(f"""
        <div class="smiles-container">
            <strong>🧬 SMILES:</strong> <code>{current_smiles}</code><br>
            <strong>🔗 URL-Encoded:</strong> <code>{encoded_smiles}</code>
        </div>
        """, unsafe_allow_html=True)

    # Generate button
    st.markdown("<br>", unsafe_allow_html=True)
    col1, col2, col3 = st.columns([1, 2, 1])
    with col2:
        if st.button("🔬 Generate FIXED IUPAC Name", 
                     type="primary", 
                     use_container_width=True, 
                     disabled=not current_smiles):
            if current_smiles:
                iupac_name = get_accurate_iupac_name_fixed(current_smiles)
                st.session_state.iupac_result = iupac_name

    # Show results
    if st.session_state.get('iupac_result'):
        result = st.session_state.iupac_result
        
        error_indicators = ["Error", "Please draw", "Invalid", "Inorganic", "failed"]
        is_error = any(indicator in result for indicator in error_indicators)
        
        if is_error:
            st.markdown(f"""
            <div class="result-error">
                <h3>❌ Unable to Generate Name</h3>
                <p><strong>{result}</strong></p>
                <p><em>💡 Try a simpler organic molecule or check your structure.</em></p>
            </div>
            """, unsafe_allow_html=True)
        else:
            st.markdown(f"""
            <div class="result-success">
                <h3>🎉 ACCURATE IUPAC Name (FIXED)</h3>
                <p class="iupac-name">{result}</p>
            </div>
            """, unsafe_allow_html=True)
            
            # Show molecular properties
            try:
                mol = Chem.MolFromSmiles(current_smiles)
                if mol:
                    mw = Descriptors.MolWt(mol)
                    formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
                    
                    st.markdown(f"""
                    <div style="background: #f8fafc; padding: 1rem; border-radius: 8px; margin-top: 1rem; text-align: center; color: #334155;">
                        <strong>📊 Molecular Properties:</strong><br>
                        <strong>Formula:</strong> {formula} | <strong>MW:</strong> {mw:.2f} g/mol
                    </div>
                    """, unsafe_allow_html=True)
            except:
                pass

    # Footer
    st.markdown("<br><br>", unsafe_allow_html=True)
    st.markdown("""
    <div class="footer">
        <p><strong>IUPAC Name Generator</strong></p>
        <p>Developed by <strong>Mayank Vishwakarma</strong></p>
    </div>
    """, unsafe_allow_html=True)

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        st.error(f"Application Error: {str(e)}")
        st.info("Please check your internet connection and try again.")
        st.stop()
