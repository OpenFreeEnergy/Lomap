[![CI](https://github.com/OpenFreeEnergy/Lomap/actions/workflows/CI.yaml/badge.svg)](https://github.com/OpenFreeEnergy/Lomap/actions/workflows/CI.yaml)


# Lomap
Alchemical free energy calculations hold increasing promise 
as an aid to drug discovery efforts. However, applications of 
these techniques in discovery projects have been relatively 
rare, partly because of the difficulty of planning and setting up 
calculations. The lead optimization mapper (LOMAP) was 
introduced as an automated algorithm to plan efficient relative 
free energy calculations between potential ligands within 
a substantial of compounds. The original LOMAP code was mainly
based on commercial APIs such as OpenEye and Schrodinger. The aim 
of this project is to develop a new version of LOMAP based on free
avalaible APIs such as RDKit offering the scientific community a 
free tool to plan in advance binding free energy calculations.

## Prerequisites
* RDKit Release > 2021
* NetworkX
* Matplotlib 
* python > 3.8

Authors
-------
* Gaetano Calabro' <gcalabro@uci.edu>
* David Mobley <dmobley@uci.edu>

* Mark Mackey
* Lester Hedges
* Antonia S J S Mey
* Jenke Scheen

## Installation

`conda install -c conda-forge lomap2`

Alternatively, you can install the prerequisites and then install from source
with:

`python setup.py install`

Usage
-----
As a commandline tool LOMAP can be simply used as:
`
lomap test/basic/
`

For a basic example run:
`python examples/example.py`

For generating radial graphs with a hub, run:
`python examples/example_radial.py`

If you would rather use the API directly, try:

```python
import lomap

# Generate the molecule database starting from a directory containing .mol2 files

db_mol = lomap.DBMolecules("python string pointing to a directory with mol2 files", output=True)

    #More graphing options:
    # Use the complete radial graph option. The ligand with the most structural similarity to all of the others will be picked as the 'lead compounds' and used as the central compound.
    db_mol = lomap.DBMolecules("python string pointing to a directory with mol2 files", output=True, radial=True)

    # Use a radial graph with a manually specified hub compound
    db_mol = lomap.DBMolecules("python string pointing to a directory with mol2 files", output=True, radial=True, hub=filename.mol2)

    # Use a radial graph with a manually specified hub compound and fast graphing option
    #the fast graphing option create the initial graph by connecting the hub ligand with the possible surrounding ligands and add surrounding edges based on the similarities accoss surrounding nodes
    db_mol = lomap.DBMolecules("python string pointing to a directory with mol2 files", output=True, radial=True, hub=filename.mol2, fast=True)

# Calculate the similarity matrix betweeen the database molecules. Two molecules are generated
# related to the scrict rule and loose rule 

strict, loose = db_mol.build_matrices()

# Generate the NetworkX graph and output the results
nx_graph = db_mol.build_graph() 


# Calculate the Maximum Common Subgraph (MCS) between 
# the first two molecules in the molecule database 
# ignoring hydrogens and depicting the mapping in a file
    
MC = lomap.MCS.getMapping(db_mol[0].getMolecule(), db_mol[1].getMolecule(), hydrogens=False, fname='mcs.png')


# Alchemical transformation are usually performed between molecules with
# the same charges. However, it is possible to allow this transformation
# manually setting the electrostatic score for the whole set of molecules 
# producing a connected graph. The electrostatic scrore must be in the 
# range [0,1]


db_mol = lomap.DBMolecules("python string pointing to a directory with mol2 files", output=True, ecrscore=0.1)
strict, loose = db_mol.build_matrices()
nx_graph = db_mol.build_graph() 
```

