[![CI](https://github.com/OpenFreeEnergy/Lomap/actions/workflows/CI.yaml/badge.svg)](https://github.com/OpenFreeEnergy/Lomap/actions/workflows/CI.yaml)
[![Documentation Status](https://readthedocs.org/projects/lomap/badge/?version=latest)](https://lomap.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16898468.svg)](https://doi.org/10.5281/zenodo.16898468)
# Lomap
The Lead Optimization Mapper (LOMAP) is an automated algorithm for planning
efficient relative free energy calculation networks across a set of ligands,
built on freely available tools such as RDKit. 
With the optional [`gufe`](https://github.com/OpenFreeEnergy/gufe) dependency installed, 
it also integrates with the [Open Free Energy](https://openfree.energy) ecosystem.
The method is described in the original
[LOMAP publication](https://doi.org/10.1007/s10822-013-9678-y).

## Installation

`lomap` is available on conda-forge as the `lomap2` package. See the
[installation documentation](https://lomap.readthedocs.io/en/latest/installation.html)
for install instructions, including the development install and optional
dependencies (`gufe` and `pygraphviz`).

## Quickstart

This example uses LOMAP's optional `gufe` bindings to load two example ligands
bundled with the package and plan a perturbation network between them with the
default mapper and scorer.

```python
# requires the optional `gufe` dependency (see Installation)
import importlib.resources

import lomap
from gufe import SmallMoleculeComponent

# Two example ligands ship with the package under lomap.tests.data
data = importlib.resources.files("lomap.tests.data")
ligands = [
    SmallMoleculeComponent.from_sdf_file(data / name)
    for name in ["lig_41.sdf", "lig_74.sdf"]
]

# Build a LigandNetwork using LOMAP's scoring and network-construction rules
network = lomap.generate_lomap_network(
    ligands=ligands,
    mappers=lomap.LomapAtomMapper(),
    scorer=lomap.default_lomap_score,
)
print(f"{len(network.nodes)} ligands, {len(network.edges)} edges")
```

For proposing individual mappings, customising edge scores, and the full set of
network options, see the [documentation](https://lomap.readthedocs.io).

## Authors

See [AUTHORS.md](https://github.com/OpenFreeEnergy/Lomap/blob/main/AUTHORS.md).

## History

Open Free Energy took over maintenance of LOMAP in 2022. The original
development repository has since been archived, but can be found at
https://github.com/MobleyLab/Lomap.
