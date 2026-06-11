"""
Lomap2
======

Alchemical free energy calculations hold increasing promise as an aid to drug
discovery efforts. However, applications of these techniques in discovery
projects have been relatively few, partly because of the difficulty of planning
and setting up calculations. The Lead Optimization Mapper (LOMAP) is an
automated algorithm to plan efficient relative free energy calculations between
potential ligands within a substantial of compounds.

Authors: Gaetano Calabro' <gcalabro@uci.edu>
         David Mobley     <dmobley@uci.edu>


Licence: MIT

URL: https://github.com/OpenFreeEnergy/Lomap
"""

from importlib.metadata import version

__version__ = version("lomap2")

from .dbmol import DBMolecules, Molecule, SMatrix
from .gufe_bindings import (
    LomapAtomMapper,
    default_lomap_score,
    generate_lomap_network,
)
from .mcs import MCS

# Issue #127
del dbmol  # type: ignore[name-defined] # noqa: F821
del mcs  # type: ignore[name-defined] # noqa: F821

from . import _due
