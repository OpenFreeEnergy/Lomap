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

### Latest Release

From conda-forge (note the package name is `lomap2`):

```bash
conda install -c conda-forge lomap2
# or, equivalently
mamba install -c conda-forge lomap2
```

## Optional dependencies

`lomap` has optional dependencies that extend its capabilities:

### gufe

The atom mapping, scoring, and network-planning API is provided through optional
[`gufe`](https://gufe.openfree.energy/en/latest/) bindings, which let `lomap`
interoperate seamlessly with the rest of the
[Open Free Energy ecosystem](https://openfree.energy/projects):

```bash
conda install -c conda-forge gufe
```

See the [gufe bindings API
documentation](https://lomap.readthedocs.io/en/latest/api.html#gufe-bindings-api)
for how to use them.

### pygraphviz

The `GraphGen` class can plot the network graph through its `draw()` method,
which requires [pygraphviz](https://pygraphviz.github.io/):

```bash
conda install -c conda-forge pygraphviz
# or
python -m pip install pygraphviz
```

### Development Version
Alternatively, you can install the development version of `lomap` directly from the `main` branch of this repository.

```bash
conda env create -f environment.yaml
conda activate lomap-env
pip install -e .
```

Quickstart
----------

This example loads two ligands from SDF files and plans a perturbation network between
them using LOMAP's default atom mapper and scorer.

```python
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

## Deprecated APIs

The legacy `DBMolecules` API and CLI are deprecated
and will be removed in the next major release; new code should use
`generate_lomap_network` instead. See
[issue #138](https://github.com/OpenFreeEnergy/Lomap/issues/138) and the
[legacy documentation](https://lomap.readthedocs.io/en/latest/legacy.html).

## Authors

See [AUTHORS.md](https://github.com/OpenFreeEnergy/Lomap/blob/main/AUTHORS.md).

## History

Open Free Energy took over maintenance of LOMAP in 2022. The original
development repository has since been archived, but can be found at
https://github.com/MobleyLab/Lomap.
