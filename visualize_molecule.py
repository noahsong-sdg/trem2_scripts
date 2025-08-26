#!/usr/bin/env python3
"""
Script to visualize molecules from SDF files using RDKit
"""

import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
import matplotlib.pyplot as plt

def visualize_sdf_conformers(sdf_file, max_conformers=5):
    """
    Visualize multiple conformers from an SDF file
    
    Args:
        sdf_file (str): Path to the SDF file
        max_conformers (int): Maximum number of conformers to display
    """
    # Read the SDF file
    suppl = Chem.SDMolSupplier(sdf_file)
    
    if not suppl:
        print(f"Error: Could not read SDF file {sdf_file}")
        return
    
    # Collect molecules and their docking scores
    mols = []
    scores = []
    
    for mol in suppl:
        if mol is not None:
            mols.append(mol)
            # Try to get docking score from SDF properties
            score = mol.GetProp('docking_score') if mol.HasProp('docking_score') else 'N/A'
            scores.append(score)
    
    print(f"Found {len(mols)} conformers in {sdf_file}")
    
    if not mols:
        print("No valid molecules found in the SDF file")
        return
    
    # Limit the number of conformers to display
    mols_to_show = mols[:max_conformers]
    scores_to_show = scores[:max_conformers]
    
    # Create a grid of molecules
    if len(mols_to_show) == 1:
        # Single molecule
        mol = mols_to_show[0]
        print(f"Docking Score: {scores_to_show[0]}")
        
        # Generate 2D coordinates if not present
        if mol.GetNumConformers() == 0:
            AllChem.Compute2DCoords(mol)
        
        # Draw the molecule
        img = Draw.MolToImage(mol, size=(400, 400))
        img.show()
        
    else:
        # Multiple molecules
        print("Docking Scores:")
        for i, score in enumerate(scores_to_show):
            print(f"  Conformer {i+1}: {score}")
        
        # Generate 2D coordinates for all molecules
        for mol in mols_to_show:
            if mol.GetNumConformers() == 0:
                AllChem.Compute2DCoords(mol)
        
        # Create a grid image
        img = Draw.MolsToGridImage(
            mols_to_show,
            molsPerRow=min(3, len(mols_to_show)),
            subImgSize=(300, 300),
            legends=[f"Score: {score}" for score in scores_to_show]
        )
        img.show()

def visualize_3d_conformer(sdf_file, conformer_index=0):
    """
    Visualize a specific 3D conformer
    
    Args:
        sdf_file (str): Path to the SDF file
        conformer_index (int): Index of the conformer to visualize (0-based)
    """
    suppl = Chem.SDMolSupplier(sdf_file)
    
    if not suppl or conformer_index >= len(suppl):
        print(f"Error: Conformer {conformer_index} not found")
        return
    
    mol = suppl[conformer_index]
    if mol is None:
        print(f"Error: Invalid molecule at index {conformer_index}")
        return
    
    # Get docking score
    score = mol.GetProp('docking_score') if mol.HasProp('docking_score') else 'N/A'
    print(f"Conformer {conformer_index + 1}, Docking Score: {score}")
    
    # Generate 2D coordinates for visualization
    mol_2d = Chem.Mol(mol)
    AllChem.Compute2DCoords(mol_2d)
    
    # Draw the molecule
    img = Draw.MolToImage(mol_2d, size=(400, 400))
    img.show()

def main():
    sdf_file = "ZINC000069019136.sdf"
    
    print("ZINC000069019136 Molecule Visualization")
    print("=" * 40)
    
    # Check if file exists
    try:
        with open(sdf_file, 'r') as f:
            pass
    except FileNotFoundError:
        print(f"Error: File {sdf_file} not found")
        return
    
    print("\nOptions:")
    print("1. View first 5 conformers")
    print("2. View specific conformer")
    print("3. View all conformers (may be slow)")
    
    choice = input("\nEnter your choice (1-3): ").strip()
    
    if choice == "1":
        visualize_sdf_conformers(sdf_file, max_conformers=5)
    elif choice == "2":
        conformer_idx = input("Enter conformer index (0-19): ").strip()
        try:
            idx = int(conformer_idx)
            visualize_3d_conformer(sdf_file, idx)
        except ValueError:
            print("Invalid index")
    elif choice == "3":
        visualize_sdf_conformers(sdf_file, max_conformers=20)
    else:
        print("Invalid choice")

if __name__ == "__main__":
    main()
