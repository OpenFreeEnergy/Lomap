import importlib.resources

import pytest
from rdkit import Chem

import lomap

try:
    import gufe

    HAS_GUFE = True
except ImportError:
    HAS_GUFE = False


@pytest.fixture
def basic():
    data_path = importlib.resources.files("lomap.tests.basic")

    fns = [
        data_path / f
        for f in (
            "1,3,7-trimethylnaphthalene.mol2",
            "1-butyl-4-methylbenzene.mol2",
            "2,6-dimethylnaphthalene.mol2",
            "2-methyl-6-propylnaphthalene.mol2",
            "2-methylnaphthalene.mol2",
            "2-naftanol.mol2",
            "methylcyclohexane.mol2",
            "toluene.mol2",
        )
    ]

    mols = [gufe.SmallMoleculeComponent(Chem.MolFromMol2File(str(f), removeHs=False)) for f in fns]

    return mols


@pytest.mark.skipif(HAS_GUFE, reason="requires not having gufe installed")
def test_generate_network_nogufe_failure():
    msg = "gufe is required to use `generate_lomap_network` but is not installed."
    with pytest.raises(ImportError, match=msg):
        _ = lomap.generate_lomap_network(
            ligands=None,
            mappers=None,
            scorer=None,
        )


@pytest.mark.skipif(not HAS_GUFE, reason="requires gufe installed")
def test_generate_network_smoketest(basic):
    with pytest.deprecated_call(match="'molecules' is deprecated, please use 'ligands'"):
        network = lomap.generate_lomap_network(
            molecules=basic,
            mappers=lomap.LomapAtomMapper(),
            scorer=lomap.default_lomap_score,
        )

        assert isinstance(network, gufe.LigandNetwork)


@pytest.mark.skipif(not HAS_GUFE, reason="requires gufe installed")
def test_overdefined_deprecated_generate_network(basic):
    with pytest.raises(ValueError, match="Both 'molecules' and 'ligands' are defined"):
        lomap.generate_lomap_network(
            molecules=basic,
            ligands=basic,
            mappers=lomap.LomapAtomMapper(),
            scorer=lomap.default_lomap_score,
        )
