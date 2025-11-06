import importlib.resources
import gufe
import pytest
from rdkit import Chem
from functools import partial
from numpy.testing import assert_equal
import lomap
from lomap.gufe_bindings.scorers import ecr_score
from lomap import dbmol


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
    # The default charge_changes_score is now 0.1, so a network of ligands
    # with different net charges should now always be connected
    network = lomap.generate_lomap_network(
        molecules=smcs,
        mappers=lomap.LomapAtomMapper(),
        scorer=partial(lomap.default_lomap_score),
    )
    assert network.is_connected() == True


def test_ecr_consistency(smcs):
    lig_4 = [m for m in smcs if m.name == 'lig_4'][0]
    other_mols = [m for m in smcs if m.name != 'lig_4']

    score_zero = []
    score_nonzero = []
    score_dbmol = []

    mapper = lomap.LomapAtomMapper()

    for m in other_mols:
        mapping = next(mapper.suggest_mappings(lig_4, m))
        score_zero.append(ecr_score(mapping, charge_changes_score=0.0))
        score_nonzero.append(ecr_score(mapping, charge_changes_score=0.1))
        score_dbmol.append(dbmol.ecr(lig_4.to_rdkit(), m.to_rdkit()))

    assert_equal(score_zero, score_dbmol)
    assert not any(score_zero)
    assert all([i == 0.1 for i in score_nonzero])


def test_default_and_explicit_charge_change_score_same(smcs):
    mapper = lomap.LomapAtomMapper()
    mapping = next(mapper.suggest_mappings(smcs[0], smcs[1]))
    result_default = lomap.default_lomap_score(mapping)
    result_explicit = lomap.default_lomap_score(mapping,
                                                charge_changes_score=0.1)

    assert result_default == result_explicit
