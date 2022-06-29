"""
Regression tests for the MCS class

These might start failing if RDKit changes canonical order of atoms.
"""
import pkg_resources
import pytest
from rdkit import Chem

from lomap import mcs


def _rf(fn):
    # get path to file from inside lomap installation
    f = pkg_resources.resource_filename('lomap', 'test/' + fn)
    return f.replace('/lomap/test', '/test')


@pytest.fixture(scope='session')
def toluene():
    return Chem.MolFromMol2File(_rf('basic/toluene.mol2'))


@pytest.fixture(scope='session')
def methylcyclohexane():
    return Chem.MolFromMol2File(_rf('basic/methylcyclohexane.mol2'))


@pytest.fixture(scope='session')
def methylnaphthalene():
    return Chem.MolFromMol2File(_rf('basic/2-methylnaphthalene.mol2'))


@pytest.fixture(scope='session')
def trimethylnaphthalene():
    return Chem.MolFromMol2File(_rf('basic/1,3,7-trimethylnaphthalene.mol2'))


@pytest.fixture(scope='session')
def butylmethylbenzene():
    return Chem.MolFromMol2File(_rf('basic/1-butyl-4-methylbenzene.mol2'))


@pytest.fixture(scope='session')
def dimethylnaphthalene():
    return Chem.MolFromMol2File(_rf('basic/2,6-dimethylnaphthalene.mol2'))


@pytest.fixture(scope='session')
def methylpropylnaphthalene():
    return Chem.MolFromMol2File(_rf('basic/2-methyl-6-propylnaphthalene.mol2'))


@pytest.fixture(scope='session')
def naphthol():
    return Chem.MolFromMol2File(_rf('basic/2-naftanol.mol2'))


@pytest.fixture(scope='session')
def toluene_explicit():
    return Chem.MolFromMol2File(_rf('basic/toluene.mol2'), removeHs=False)


@pytest.fixture(scope='session')
def methylcyclohexane_explicit():
    return Chem.MolFromMol2File(_rf('basic/methylcyclohexane.mol2'),
                                removeHs=False)


@pytest.fixture(scope='session')
def methylnaphthalene_explicit():
    return Chem.MolFromMol2File(_rf('basic/2-methylnaphthalene.mol2'),
                                removeHs=False)


@pytest.fixture(scope='session')
def trimethylnaphthalene_explicit():
    return Chem.MolFromMol2File(_rf('basic/1,3,7-trimethylnaphthalene.mol2'),
                                removeHs=False)


@pytest.fixture(scope='session')
def butylmethylbenzene_explicit():
    return Chem.MolFromMol2File(_rf('basic/1-butyl-4-methylbenzene.mol2'),
                                removeHs=False)


@pytest.fixture(scope='session')
def dimethylnaphthalene_explicit():
    return Chem.MolFromMol2File(_rf('basic/2,6-dimethylnaphthalene.mol2'),
                                removeHs=False)


@pytest.fixture(scope='session')
def methylpropylnaphthalene_explicit():
    return Chem.MolFromMol2File(_rf('basic/2-methyl-6-propylnaphthalene.mol2'),
                                removeHs=False)


@pytest.fixture(scope='session')
def naphthol_explicit():
    return Chem.MolFromMol2File(_rf('basic/2-naftanol.mol2'), removeHs=False)


OTHER_MOLS = [
    'methylcyclohexane', 'methylnaphthalene', 'trimethylnaphthalene',
    'butylmethylbenzene', 'dimethylnaphthalene', 'methylpropylnaphthalene',
    'naphthol',
]


@pytest.fixture(scope='function', params=OTHER_MOLS)
def regression_mappings(toluene, methylcyclohexane, methylnaphthalene,
                        trimethylnaphthalene, butylmethylbenzene,
                        dimethylnaphthalene, methylpropylnaphthalene,
                        naphthol, request):
    regression_mappings = {
        'methylcyclohexane': (methylcyclohexane,
                              [(0, 0), (1, 1), (2, 6), (3, 5), (4, 4), (5, 3),
                               (6, 2)]),
        'methylnaphthalene': (methylnaphthalene,
                              [(0, 0), (1, 1), (2, 2), (3, 3), (4, 8), (5, 9),
                               (6, 10)]),
        'trimethylnaphthalene': (trimethylnaphthalene,
                                 [(0, 0), (1, 1), (2, 10), (3, 9), (4, 8),
                                  (5, 3), (6, 2)]),
        'butylmethylbenzene': (butylmethylbenzene,
                               [(0, 3), (1, 4), (2, 5), (3, 6), (4, 7), (5, 8),
                                (6, 9)]),
        'dimethylnaphthalene': (dimethylnaphthalene,
                                [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4),
                                 (5, 10), (6, 11)]),
        'methylpropylnaphthalene': (methylpropylnaphthalene,
                                    [(0, 2), (1, 3), (2, 8), (3, 7), (4, 6),
                                     (5, 5), (6, 4)]),
        'naphthol': (naphthol,
                     [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4), (5, 9), (6, 10)])
    }
    other, other_mapping = regression_mappings[request.param]
    # a bit paranoid, but return copies
    return Chem.Mol(toluene), Chem.Mol(other), other_mapping


@pytest.fixture(scope='function', params=OTHER_MOLS)
def regression_mappings_explicit(toluene_explicit,
                                 methylcyclohexane_explicit,
                                 methylnaphthalene_explicit,
                                 trimethylnaphthalene_explicit,
                                 butylmethylbenzene_explicit,
                                 dimethylnaphthalene_explicit,
                                 methylpropylnaphthalene_explicit,
                                 naphthol_explicit,
                                 request):
    regression_mappings = {
        'methylcyclohexane': (methylcyclohexane_explicit,
                              [(0, 0),
                               (1, 1),
                               (2, 6),
                               (3, 5),
                               (4, 4),
                               (5, 3),
                               (6, 2),
                               (7, 7),
                               (8, 8),
                               (9, 9),
                               (10, 20),
                               (11, 18),
                               (12, 16),
                               (13, 14),
                               (14, 12)]),
        'methylnaphthalene': (methylnaphthalene_explicit,
                              [(0, 0),
                               (1, 1),
                               (2, 2),
                               (3, 3),
                               (4, 8),
                               (5, 9),
                               (6, 10),
                               (7, 11),
                               (8, 12),
                               (9, 13),
                               (10, 14),
                               (11, 4),
                               (12, 7),
                               (13, 19),
                               (14, 20)]),
        'trimethylnaphthalene': (trimethylnaphthalene_explicit,
                                 [(0, 0),
                                  (1, 1),
                                  (2, 10),
                                  (3, 9),
                                  (4, 8),
                                  (5, 3),
                                  (6, 2),
                                  (7, 13),
                                  (8, 14),
                                  (9, 15),
                                  (10, 20),
                                  (11, 19),
                                  (12, 7),
                                  (13, 4),
                                  (14, 16)]),
        'butylmethylbenzene': (butylmethylbenzene_explicit,
                               [(0, 3),
                                (1, 4),
                                (2, 5),
                                (3, 6),
                                (4, 7),
                                (5, 8),
                                (6, 9),
                                (7, 2),
                                (8, 18),
                                (9, 19),
                                (10, 20),
                                (11, 21),
                                (12, 10),
                                (13, 22),
                                (14, 23)]),
        'dimethylnaphthalene': (dimethylnaphthalene_explicit,
                                [(0, 0),
                                 (1, 1),
                                 (2, 2),
                                 (3, 3),
                                 (4, 4),
                                 (5, 10),
                                 (6, 11),
                                 (7, 12),
                                 (8, 13),
                                 (9, 14),
                                 (10, 15),
                                 (11, 16),
                                 (12, 5),
                                 (13, 9),
                                 (14, 23)]),
        'methylpropylnaphthalene': (methylpropylnaphthalene_explicit,
                                    [(0, 2),
                                     (1, 3),
                                     (2, 8),
                                     (3, 7),
                                     (4, 6),
                                     (5, 5),
                                     (6, 4),
                                     (7, 1),
                                     (8, 19),
                                     (9, 20),
                                     (10, 23),
                                     (11, 22),
                                     (12, 9),
                                     (13, 12),
                                     (14, 21)]),
        'naphthol': (naphthol_explicit,
                     [(0, 0),
                      (1, 1),
                      (2, 2),
                      (3, 3),
                      (4, 4),
                      (5, 9),
                      (6, 10),
                      (9, 11),
                      (10, 12),
                      (11, 13),
                      (12, 5),
                      (13, 8),
                      (14, 18)])
    }

    other, other_mapping = regression_mappings[request.param]

    return Chem.Mol(toluene_explicit), Chem.Mol(other), other_mapping


def test_toluene_to_other(regression_mappings):
    toluene, other, ref_mapping = regression_mappings

    mapper = mcs.MCS(toluene, other)

    mapping = mapper.heavy_atom_mcs_map()

    assert sorted(mapping) == ref_mapping


def test_toluene_to_other_explicit(regression_mappings_explicit):
    toluene, other, ref_mapping = regression_mappings_explicit

    mapper = mcs.MCS(toluene, other)

    mapping = mapper.all_atom_match_list()

    print(mapping)

    mapping = [tuple(map(int, row.split(':'))) for row in mapping.split(',')]

    assert sorted(mapping) == ref_mapping


@pytest.mark.parametrize('element_change', [True, False])
def test_no_element_change(naphthol_explicit, methylnaphthalene_explicit,
                           element_change):
    # map naphthol to methylnaphthalene with and without element changes
    # the -OH to -CH3 shouldn't get mapped when element changes not allowed
    mapper = mcs.MCS(naphthol_explicit,
                     methylnaphthalene_explicit,
                     element_change=element_change)

    expected_heavy = 11 if element_change else 10
    expected_all = expected_heavy + (8 if element_change else 7)

    assert len(mapper.heavy_atom_mcs_map()) == expected_heavy
    assert len(mapper.all_atom_match_list().split(',')) == expected_all
    # paranoid, check the oxygen didn't get mapped
    saw_oxygen = any(naphthol_explicit.GetAtomWithIdx(i).GetAtomicNum() == 8
                     for (i, _) in mapper.heavy_atom_mcs_map())
    assert saw_oxygen == element_change


@pytest.mark.parametrize('element_change', [True, False])
def test_no_element_change_hydrogen_to_heavy(toluene_explicit,
                                             dimethylnaphthalene_explicit,
                                             element_change):
    """
    Checks that hydrogens in toluene are not matches to carbons in
    dimethylnaphthalene
    """
    mapper = mcs.MCS(toluene_explicit, dimethylnaphthalene_explicit,
                     threed=True, element_change=element_change)
    expected_heavy = 7 # 7 + 2 heavy to hydrogen
    expected_all = 15 if element_change else 13

    assert len(mapper.heavy_atom_mcs_map()) == expected_heavy
    assert len(mapper.all_atom_match_list().split(',')) == expected_all
