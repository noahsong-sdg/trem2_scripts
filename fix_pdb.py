from pdbfixer import PDBFixer
from openmm.app import PDBFile

fixer = PDBFixer(filename='../data/receptor/cluster1.pdb')
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens()
PDBFile.writeFile(fixer.topology, fixer.positions, open('../data/receptor/cluster1_fixed.pdb', 'w'))
