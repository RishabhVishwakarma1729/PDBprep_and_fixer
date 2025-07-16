import streamlit as st
from Bio.PDB import PDBParser, PDBIO, Select
from pdbfixer import PDBFixer
from openmm.app import PDBFile
import tempfile
import os

# ========== Chain Filter Class ==========
class ChainSelect(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id

    def accept_residue(self, residue):
        return residue.get_id()[0] == " "  # Exclude HETATM

    def accept_chain(self, chain):
        return chain.get_id() == self.chain_id

# ========== Main Streamlit App ==========
st.title("🔬 Clean PDB with PDBFixer")
st.write("Upload a PDB file, select a chain, remove HETATM residues, and optionally fix issues using PDBFixer.")

uploaded_file = st.file_uploader("Upload PDB file", type=["pdb"])

if uploaded_file:
    # Save uploaded file to temp location
    with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as tmp_input:
        tmp_input.write(uploaded_file.read())
        pdb_path = tmp_input.name

    # Extract available chains
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_path)
    chains = {chain.get_id() for model in structure for chain in model}
    chain_id = st.selectbox("Select chain to retain", sorted(chains))

    run_fixer = st.checkbox("🔧 Run PDBFixer (Add missing atoms, hydrogens, remove ligands, etc.)", value=True)

    if st.button("🧪 Process"):
        # Save cleaned PDB with selected chain
        cleaned_pdb = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb").name
        io = PDBIO()
        io.set_structure(structure)
        io.save(cleaned_pdb, select=ChainSelect(chain_id))

        if run_fixer:
            try:
                # Run PDBFixer on the cleaned structure
                fixer = PDBFixer(filename=cleaned_pdb)
                fixer.findMissingResidues()
                fixer.findNonstandardResidues()
                fixer.replaceNonstandardResidues()
                fixer.removeHeterogens(keepWater=False)
                fixer.findMissingAtoms()
                fixer.addMissingAtoms()
                fixer.addMissingHydrogens(pH=7.4)

                # Save fixed structure
                fixed_pdb = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb").name
                with open(fixed_pdb, "w") as out:
                    PDBFile.writeFile(fixer.topology, fixer.positions, out)

                with open(fixed_pdb, "rb") as f:
                    st.success("✅ PDBFixer applied. Download the fixed PDB:")
                    st.download_button("📥 Download Fixed PDB", f, file_name="fixed_structure.pdb")

            except Exception as e:
                st.error(f"❌ PDBFixer failed: {e}")

        else:
            # Offer cleaned PDB without PDBFixer
            with open(cleaned_pdb, "rb") as f:
                st.success("✅ Cleaned structure without PDBFixer.")
                st.download_button("📥 Download Cleaned PDB", f, file_name="cleaned_chain.pdb")
