import importlib.resources
import gufe
import pytest
from rdkit import Chem
from functools import partial
import lomap


@pytest.fixture
def smcs():
    data_path = importlib.resources.files('lomap.tests.basic')
    fns = data_path / 'thrombin_ligands.sdf'
    rdmols = [mol for mol in Chem.SDMolSupplier(str(fns), removeHs=False)]
    mols = [gufe.SmallMoleculeComponent.from_rdkit(mol) for mol in rdmols]
    return mols


def test_unconnected_lomap_network(smcs):
    network = lomap.generate_lomap_network(
        molecules=smcs,
        mappers=lomap.LomapAtomMapper(),
        scorer=partial(lomap.default_lomap_score, charge_changes_score=0.0),
    )

    assert network.is_connected() == False


def test_connected_lomap_network(smcs):
    network = lomap.generate_lomap_network(
        molecules=smcs,
        mappers=lomap.LomapAtomMapper(),
        scorer=partial(lomap.default_lomap_score, charge_changes_score=0.1),
    )

    assert network.is_connected() == True

