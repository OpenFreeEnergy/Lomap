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
    f = pkg_resources.resource_filename('lomap', 'tests/' + fn)
    return f


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


@pytest.fixture()
def whole_fragment_pair():
    """files for issue #32

    these molecules previously generated a mapping with multiple fragments
    """
    m1 = Chem.MolFromMolFile(_rf('data/lig_41.sdf'), removeHs=False)
    m2 = Chem.MolFromMolFile(_rf('data/lig_74.sdf'), removeHs=False)

    return m1, m2


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
    expected_heavy = 7  # 7 + 2 heavy to hydrogen
    expected_all = 15 if element_change else 13

    assert len(mapper.heavy_atom_mcs_map()) == expected_heavy
    assert len(mapper.all_atom_match_list().split(',')) == expected_all


@pytest.fixture
def naphthalene_shift():
    block = """\
naphthalene
     RDKit          3D

 18 19  0  0  0  0  0  0  0  0999 V2000
   -0.0251   -0.6892   -0.0803 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1579   -1.4195   -0.1275 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3995   -0.8252   -0.0200 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4400    0.5457    0.1394 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2704    1.2774    0.1872 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0223    0.6814    0.0793 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1469    1.4130    0.1271 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3984    0.8248    0.0200 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4271   -0.5460   -0.1391 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2650   -1.2897   -0.1885 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1134   -2.5024   -0.2538 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.3168   -1.4007   -0.0576 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.4009    1.0420    0.2268 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.2696    2.3501    0.3110 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0667    2.4757    0.2522 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3195    1.4047    0.0580 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4173   -1.0014   -0.2227 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3249   -2.3407   -0.3116 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  7  1  0
  7  8  2  0
  8  9  1  0
  9 10  2  0
  6  1  1  0
 10  1  1  0
  2 11  1  0
  3 12  1  0
  4 13  1  0
  5 14  1  0
  7 15  1  0
  8 16  1  0
  9 17  1  0
 10 18  1  0
M  END
    """

    return Chem.MolFromMolBlock(block, removeHs=False)


@pytest.fixture
def benzoxazole_shift():
    block = """\
benzoxazole
     RDKit          3D

 14 15  0  0  0  0  0  0  0  0999 V2000
    1.1790   -1.1478   -0.0955 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0924   -0.6927   -0.0828 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3456   -1.2599   -0.1876 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5080   -0.5410   -0.1410 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3899    0.8321    0.0211 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1444    1.4227    0.1283 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0254    0.6749    0.0787 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3085    0.9716    0.1531 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.0357   -0.1154    0.0505 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4354   -2.3353   -0.3144 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4840   -1.0008   -0.2247 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3142    1.4080    0.0586 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0701    2.4861    0.2533 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.1347   -0.1915    0.0761 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  3  4  1  0
  4  5  2  0
  5  6  1  0
  6  7  2  0
  7  8  1  0
  8  9  1  0
  9  1  2  0
  7  2  1  0
  3 10  1  0
  4 11  1  0
  5 12  1  0
  6 13  1  0
  9 14  1  0
M  END

"""

    return Chem.MolFromMolBlock(block, removeHs=False)


def mismatch(m1, m2, shift, mapping) -> float:
    """Total geometric mismatch in the mapping"""
    total = 0

    ci = m1.GetConformer()
    cj = m2.GetConformer()

    m1_idx, m2_idx = list(zip(*mapping))

    if shift:
        offset = mcs.substructure_centre(m1, m1_idx) - mcs.substructure_centre(m2, m2_idx)
    else:
        offset = mcs.Point3D(0.0, 0.0, 0.0)

    for i, j in mapping:
        total += (ci.GetAtomPosition(i) - cj.GetAtomPosition(j) - offset).Length()

    return total


@pytest.mark.parametrize('use_shift,refval1,refval2', [
    (True, 0.24450728924894888, 14.749875486694247),  # distance mismatch values with and without CoG shift
    (False, 0.2526494046771085, 0.263024457571414),
    # the fact first case can be 0.24 when the second case is either 0.25 or 0.26 is the crux of issue #26
])
def test_shift_parameter(use_shift, refval1, refval2, naphthalene_shift, benzoxazole_shift):
    mapper = mcs.MCS(naphthalene_shift, benzoxazole_shift, threed=True, shift=use_shift)
    mapping = mapper.heavy_atom_mcs_map()

    assert len(mapper.heavy_atom_mcs_map()) == 6
    assert mismatch(naphthalene_shift, benzoxazole_shift, True, mapping) == pytest.approx(refval1)
    assert mismatch(naphthalene_shift, benzoxazole_shift, False, mapping) == pytest.approx(refval2)


def test_whole_fragment_match_only(whole_fragment_pair):
    # issue #32
    m = mcs.MCS(whole_fragment_pair[0], whole_fragment_pair[1])

    frags = Chem.GetMolFrags(m.mcs_mol)

    assert len(frags) == 1
