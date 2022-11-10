import pkg_resources

from lomap.mcs import MCS
from rdkit import RDLogger, Chem
from lomap.dbmol import DBMolecules
from lomap import dbmol
import multiprocessing
import math
import argparse
import logging

import pytest


def _rf(fn):
    # get path to file from inside lomap installation
    f = pkg_resources.resource_filename('lomap', 'test/' + fn)
    return f.replace('/lomap/test', '/test')


@pytest.mark.parametrize('fn1, fn2, max3d_arg, threed_arg, exp_mcsr', [
    (_rf('transforms/phenyl.sdf'), _rf('transforms/toluyl.sdf'), False, 1000, math.exp(-0.1 * (6 + 7 - 2*6))),
    (_rf('transforms/phenyl.sdf'), _rf('transforms/chlorophenyl.sdf'), False, 1000, math.exp(-0.1 * (6 + 7 - 2*6))),
    (_rf('transforms/toluyl.sdf'), _rf('transforms/chlorophenyl.sdf'), False, 1000, 1),
    (_rf('transforms/toluyl.sdf'), _rf('transforms/chlorophenol.sdf'), False, 1000, math.exp(-0.1 * (7 + 8 - 2*7))),
    (_rf('transforms/phenyl.sdf'), _rf('transforms/napthyl.sdf'), False, 1000, math.exp(-0.1 * (8 + 12 - 2*8))),
   (_rf('transforms/chlorophenyl.sdf'), _rf('transforms/fluorophenyl.sdf'), False, 1000, 1),
   (_rf('transforms/chlorophenyl.sdf'), _rf('transforms/bromophenyl.sdf'), False, 1000, 1),
   (_rf('transforms/chlorophenyl.sdf'), _rf('transforms/iodophenyl.sdf'), False, 1000, 1),

   # Compare with and without 3D pruning
   (_rf('transforms/chlorophenyl.sdf'), _rf('transforms/chlorophenyl2.sdf'), False, 1000, 1),
   (_rf('transforms/chlorophenyl.sdf'), _rf('transforms/chlorophenyl2.sdf'), False, 2, math.exp(-0.1 * (7 + 7 - 2*6))),

   # Compare with and without 3D matching
   (_rf('transforms/chlorotoluyl1.sdf'), _rf('transforms/chlorotoluyl2.sdf'), False, 1000, 1),
   (_rf('transforms/chlorotoluyl1.sdf'), _rf('transforms/chlorotoluyl2.sdf'), True, 1000, 1)
])
@pytest.mark.skip("mcsr issue")
def test_mcsr(fn1, fn2, max3d_arg, threed_arg, exp_mcsr):
    # MolA, molB, 3D?, max3d, mcsr, atomic_number_rule
    logging.basicConfig(format='%(message)s', level=logging.CRITICAL)

    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)

    mola = Chem.MolFromMolFile(fn1, sanitize=False, removeHs=False)
    molb = Chem.MolFromMolFile(fn1, sanitize=False, removeHs=False)
    MC = MCS(mola, molb, time=20, verbose='info', max3d=max3d_arg, threed=threed_arg)
    mcsr = MC.mcsr()

    assert mcsr == pytest.approx(exp_mcsr)

@pytest.mark.parametrize('fn1, fn2, max3d_arg, threed_arg, exp_atnum', [
    (_rf('transforms/phenyl.sdf'), _rf('transforms/toluyl.sdf'), False, 1000, 1),
    (_rf('transforms/phenyl.sdf'), _rf('transforms/chlorophenyl.sdf'), False, 1000, 1),
    (_rf('transforms/toluyl.sdf'), _rf('transforms/chlorophenyl.sdf'), False, 1000, math.exp(-0.1 * 0.5)),
    (_rf('transforms/toluyl.sdf'), _rf('transforms/chlorophenol.sdf'), False, 1000, math.exp(-0.1 * 0.5)),
    (_rf('transforms/phenyl.sdf'), _rf('transforms/napthyl.sdf'), False, 1000, 1),
   (_rf('transforms/chlorophenyl.sdf'), _rf('transforms/fluorophenyl.sdf'), False, 1000, math.exp(-0.1 * 0.5 )),
   (_rf('transforms/chlorophenyl.sdf'), _rf('transforms/bromophenyl.sdf'), False, 1000, math.exp(-0.1 * 0.15)),
   (_rf('transforms/chlorophenyl.sdf'), _rf('transforms/iodophenyl.sdf'), False, 1000, math.exp(-0.1 * 0.35)),

   # Compare with and without 3D pruning
   (_rf('transforms/chlorophenyl.sdf'), _rf('transforms/chlorophenyl2.sdf'), False, 1000, 1),
   (_rf('transforms/chlorophenyl.sdf'), _rf('transforms/chlorophenyl2.sdf'), False, 1000, 1),

   # Compare with and without 3D matching
   (_rf('transforms/chlorotoluyl1.sdf'), _rf('transforms/chlorotoluyl2.sdf'), False, 1000, 1),
   (_rf('transforms/chlorotoluyl1.sdf'), _rf('transforms/chlorotoluyl2.sdf'), True, 1000, math.exp(-0.05 * 2))
])
@pytest.mark.skip('atnum issue')
def test_atnum_rule(fn1, fn2, max3d_arg, threed_arg, exp_atnum):
    # MolA, molB, 3D?, max3d, mcsr, atomic_number_rule
    logging.basicConfig(format='%(message)s', level=logging.CRITICAL)

    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)

    mola = Chem.MolFromMolFile(fn1, sanitize=False, removeHs=False)
    molb = Chem.MolFromMolFile(fn1, sanitize=False, removeHs=False)
    MC = MCS(mola, molb, time=20, verbose='info', max3d=max3d_arg, threed=threed_arg)

    atnum = MC.atomic_number_rule()

    assert atnum == exp_atnum

def test_iter_next():
    logging.basicConfig(format='%(message)s', level=logging.CRITICAL)
    inst = DBMolecules(_rf('basic/'), parallel=1, verbose='off', output=False, time=20, ecrscore=0.0, name='out',
                       display=False, max=6, cutoff=0.4, radial=False, hub=None)
    for i in range(0, inst.nums()):
        inst.next()
    with pytest.raises(StopIteration):
        inst.next()


def test_get_set_add():
    logging.basicConfig(format='%(message)s', level=logging.CRITICAL)
    inst = DBMolecules(_rf('basic/'), parallel=1, verbose='off', output=False, time=20, ecrscore=0.0, name='out',
                       display=False, max=6, cutoff=0.4, radial=False, hub=None)
    with pytest.raises(IndexError):
        inst.__getitem__(inst.nums()+1)
    with pytest.raises(IndexError):
        inst.__setitem__(inst.nums()+1,inst[1])
    with pytest.raises(ValueError):
        inst.__setitem__(0, 'no_mol_obj')
    with pytest.raises(ValueError):
        inst.__add__('no_mol_obj')


# Check serial and parallel mode
def test_serial_parallel():
    logging.basicConfig(format='%(message)s', level=logging.CRITICAL)
    db = DBMolecules(_rf('basic/'))
    s_strict, s_loose = db.build_matrices()
    db.options['parallel'] = multiprocessing.cpu_count()
    p_strict, p_loose = db.build_matrices()

    assert (s_strict == p_strict).all()
    assert (s_loose == p_loose).all()


def test_read_mol2_files():
    logging.basicConfig(format='%(message)s', level=logging.CRITICAL)
    db = DBMolecules(_rf('basic/'))
    db.options['directory'] = _rf('')
    with pytest.raises(IOError):
        db.read_molecule_files()


@pytest.mark.parametrize('fn, expected', [
    (_rf('transforms/phenylfuran.sdf'), 1),
    (_rf('transforms/phenylimidazole.sdf'), math.exp(-0.1*4)),
    (_rf('transforms/phenylisoxazole.sdf'), math.exp(-0.1*4)),
    (_rf('transforms/phenyloxazole.sdf'), math.exp(-0.1*4)),
    (_rf('transforms/phenylpyrazole.sdf'), math.exp(-0.1*4)),
    (_rf('transforms/phenylpyridine1.sdf'), math.exp(-0.1*4)),
    (_rf('transforms/phenylpyridine2.sdf'), math.exp(-0.1*4)),
    (_rf('transforms/phenylpyrimidine.sdf'), math.exp(-0.1*4)),
    (_rf('transforms/phenylpyrrole.sdf'), 1),
    (_rf('transforms/phenylphenyl.sdf'), 1),
])
# Test which heterocycles I can grow (growing off a phenyl)
# Test by Max indicates that growing complex heterocycles tends
# to fail, so only allow growing phenyl, furan and pyrrole
def test_heterocycle_scores(fn, expected):
    logging.basicConfig(format='%(message)s', level=logging.CRITICAL)
    parent=Chem.MolFromMolFile(_rf('transforms/phenyl.sdf'), sanitize=False, removeHs=False)
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)

    comp=Chem.MolFromMolFile(fn, sanitize=False, removeHs=False)
    MC=MCS(parent, comp)
    assert MC.heterocycles_rule(penalty=4) == expected, 'Fail on heterocycle check for ' + fn


# Tests by Max indicate that growing a sulfonamide all in one go is
# dodgy, so disallow it
@pytest.mark.parametrize('fn, expected', [
    (_rf('transforms/cdk2_lig11.sdf'), math.exp(-0.1*4)),
    (_rf('transforms/cdk2_lig1.sdf'), 1),
    (_rf('transforms/cdk2_lig2.sdf'), math.exp(-0.1*4)),
    (_rf('transforms/cdk2_lig13.sdf'), 1),
    (_rf('transforms/cdk2_lig14.sdf'), 1),
    (_rf('transforms/cdk2_lig15.sdf'), 1),
])
def test_sulfonamide_scores(fn, expected):
    logging.basicConfig(format='%(message)s', level=logging.CRITICAL)
    parent = Chem.MolFromMolFile(_rf('transforms/cdk2_lig16.sdf'), sanitize=False, removeHs=False)
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)
    comp = Chem.MolFromMolFile(fn, sanitize=False, removeHs=False)
    MC = MCS(parent, comp)
    assert MC.sulfonamides_rule(penalty=4) == expected, 'Fail on sulfonamide check for ' + fn


# Test to check symmetry equivalence by matching atomic numbers where possible
def test_symmetry_matchheavies():
    logging.basicConfig(format='%(message)s', level=logging.CRITICAL)
    mol1 = Chem.MolFromMolFile(_rf('transforms/chlorophenol.sdf'), sanitize=False, removeHs=False)
    mol2 = Chem.MolFromMolFile(_rf('transforms/chlorophenyl.sdf'), sanitize=False, removeHs=False)
    mol3 = Chem.MolFromMolFile(_rf('transforms/chlorophenyl2.sdf'), sanitize=False, removeHs=False)
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)
    MCS1 = MCS(mol1,mol2)
    MCS2 = MCS(mol2,mol3)
    MCS3 = MCS(mol1,mol3)
    assert MCS1.mcs_mol.GetNumHeavyAtoms() == 9
    assert [int(at.GetProp('to_moli')) for at in MCS1.mcs_mol.GetAtoms()] == [0, 5, 4, 3, 2, 1, 7, 6, 9]
    assert [int(at.GetProp('to_molj')) for at in MCS1.mcs_mol.GetAtoms()] == [0, 5, 4, 3, 2, 1, 7, 6, 8]
    assert [int(at.GetProp('to_moli')) for at in MCS2.mcs_mol.GetAtoms()] == [0, 5, 4, 3, 2, 1, 7, 6, 8]
    assert [int(at.GetProp('to_molj')) for at in MCS2.mcs_mol.GetAtoms()] == [4, 5, 0, 1, 2, 3, 7, 6, 8]
    assert [int(at.GetProp('to_moli')) for at in MCS3.mcs_mol.GetAtoms()] == [4, 5, 0, 1, 2, 3, 7, 6, 9]
    assert [int(at.GetProp('to_molj')) for at in MCS3.mcs_mol.GetAtoms()] == [0, 5, 4, 3, 2, 1, 7, 6, 8]


# Test to check symmetry equivalence by matching 3D coordinates rather than atomic numbers
def test_symmetry_match3d():
    logging.basicConfig(format='%(message)s', level=logging.CRITICAL)
    mol1 = Chem.MolFromMolFile(_rf('transforms/chlorophenol.sdf'), sanitize=False, removeHs=False)
    mol2 = Chem.MolFromMolFile(_rf('transforms/chlorophenyl.sdf'), sanitize=False, removeHs=False)
    mol3 = Chem.MolFromMolFile(_rf('transforms/chlorophenyl2.sdf'), sanitize=False, removeHs=False)
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)
    MCS1 = MCS(mol1, mol2, time=20, verbose='info', max3d=1000, threed=True)
    MCS2 = MCS(mol2, mol3, time=20, verbose='info', max3d=1000, threed=True)
    MCS3 = MCS(mol1, mol3, time=20, verbose='info', max3d=1000, threed=True)
    assert MCS1.mcs_mol.GetNumHeavyAtoms() == 9
    # MCS1 and MCS2 are the same as in the matchheavies case, but MCS3 gives a diffrent answer
    assert [int(at.GetProp('to_moli')) for at in MCS1.mcs_mol.GetAtoms()] == [0, 5, 4, 3, 2, 1, 7, 6, 9]
    assert [int(at.GetProp('to_molj')) for at in MCS1.mcs_mol.GetAtoms()] == [0, 5, 4, 3, 2, 1, 7, 6, 8]
    assert [int(at.GetProp('to_moli')) for at in MCS2.mcs_mol.GetAtoms()] == [0, 5, 4, 3, 2, 1, 7, 6, 8]
    assert [int(at.GetProp('to_molj')) for at in MCS2.mcs_mol.GetAtoms()] == [4, 5, 0, 1, 2, 3, 7, 6, 8]
    assert [int(at.GetProp('to_moli')) for at in MCS3.mcs_mol.GetAtoms()] == [0, 5, 4, 3, 2, 1, 8, 6, 9]
    assert [int(at.GetProp('to_molj')) for at in MCS3.mcs_mol.GetAtoms()] == [0, 5, 4, 3, 2, 1, 7, 6, 8]


# Test to check removing atoms from the MCS when the 3D coords are too far apart
def test_clip_on_3d():
    logging.basicConfig(format='%(message)s', level=logging.CRITICAL)
    mol1 = Chem.MolFromMolFile(_rf('transforms/chlorophenyl.sdf'), sanitize=False, removeHs=False)
    mol2 = Chem.MolFromMolFile(_rf('transforms/chlorophenyl2.sdf'), sanitize=False, removeHs=False)
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)
    MCS1 = MCS(mol1, mol2, time=20, verbose='info', max3d=1000, threed=True)
    MCS2 = MCS(mol1, mol2, time=20, verbose='info', max3d=2, threed=True)
    assert MCS1.mcs_mol.GetNumHeavyAtoms() == 9
    assert MCS2.mcs_mol.GetNumHeavyAtoms() == 8


@pytest.mark.parametrize('fn1, fn2, expected', [
    (_rf('transforms/phenyl.sdf'), _rf('transforms/toluyl3.sdf'), 1),
    (_rf('transforms/toluyl3.sdf'), _rf('transforms/chlorotoluyl1.sdf'), 1),
    (_rf('transforms/toluyl3.sdf'), _rf('transforms/phenylfuran.sdf'), math.exp(-0.1*4)),
    (_rf('transforms/toluyl3.sdf'), _rf('transforms/phenylpyridine1.sdf'), math.exp(-0.1*4)),
    (_rf('transforms/phenyl.sdf'), _rf('transforms/phenylfuran.sdf'), 1),
    (_rf('transforms/phenyl.sdf'), _rf('transforms/phenylpyridine1.sdf'), 1),
    (_rf('transforms/chlorophenol.sdf'), _rf('transforms/phenylfuran.sdf'), 1),
])
# Test disallowing turning a methyl group (or larger) into a ring atom
def test_transmuting_methyl_into_ring_rule(fn1, fn2, expected):
    logging.basicConfig(format='%(message)s', level=logging.CRITICAL)
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)
    parent=Chem.MolFromMolFile(fn1, sanitize=False, removeHs=False)
    comp=Chem.MolFromMolFile(fn2, sanitize=False, removeHs=False)
    MC=MCS(parent,comp)
    assert MC.transmuting_methyl_into_ring_rule(penalty=4) == expected, 'Fail on transmuting-methyl-to-ring check for ' + fn1 + ' ' + fn2


@pytest.mark.parametrize('fn1, fn2, expected', [
    (_rf('transforms/napthyl.sdf'), _rf('transforms/tetrahydronaphthyl.sdf'), math.exp(-0.1 * 4)),
])
# Test penalising hybridization changes
def test_hybridization_rule(fn1, fn2, expected):
    logging.basicConfig(format='%(message)s', level=logging.CRITICAL)

    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)

    parent=Chem.MolFromMolFile(fn1, sanitize=False, removeHs=False)
    comp=Chem.MolFromMolFile(fn2, sanitize=False, removeHs=False)
    MC=MCS(parent,comp)

    assert MC.hybridization_rule(1.0) ==  pytest.approx(expected)


@pytest.mark.parametrize('fn1, fn2, expected', [
    (_rf('transforms/phenyl.sdf'), _rf('transforms/phenylcyclopropyl.sdf'), 1),
    (_rf('transforms/toluyl.sdf'), _rf('transforms/phenylcyclopropyl.sdf'), 1),
    # Disallowed by test_transmuting_methyl_into_ring_rule instead
    (_rf('transforms/phenylcyclopropyl.sdf'), _rf('transforms/phenylcyclobutyl.sdf'), 0.1),
    (_rf('transforms/phenylcyclopropyl.sdf'), _rf('transforms/phenylcyclopentyl.sdf'), 0.1),
    (_rf('transforms/phenylcyclopropyl.sdf'), _rf('transforms/phenylcyclononyl.sdf'), 0.1),
    (_rf('transforms/phenylcyclobutyl.sdf'), _rf('transforms/phenylcyclopentyl.sdf'), 0.1),
    (_rf('transforms/phenylcyclobutyl.sdf'), _rf('transforms/phenylcyclononyl.sdf'), 0.1),
    (_rf('transforms/phenylcyclopentyl.sdf'), _rf('transforms/phenylcyclononyl.sdf'), 1),
])
# Test disallowing turning a ring into an incompatible ring
def test_transmuting_ring_sizes_rule(fn1, fn2, expected):
    logging.basicConfig(format='%(message)s', level=logging.CRITICAL)

    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)
    parent = Chem.MolFromMolFile(fn1, sanitize=False, removeHs=False)
    comp = Chem.MolFromMolFile(fn2, sanitize=False, removeHs=False)
    MC = MCS(parent, comp)
    assert MC.transmuting_ring_sizes_rule() == expected, 'Fail on transmuting-ring-size check for ' + fn1 + ' ' + fn2


@pytest.mark.parametrize('fn1, fn2, expected', [
    (_rf('transforms/phenyl.sdf'), _rf('transforms/toluyl3.sdf'), "0:0,1:1,2:2,3:3,4:4,5:5,6:6,7:7"),
    (_rf('transforms/toluyl2.sdf'), _rf('transforms/chlorotoluyl1.sdf'), "0:0,1:1,2:2,3:3,4:4,5:5,6:6,7:8,8:9"),
    (_rf('transforms/toluyl3.sdf'), _rf('transforms/phenylfuran.sdf'), "0:0,1:1,2:2,3:3,4:4,5:5,6:6,7:7"),
])
# Test getting the mapping string out of the MCS
def test_mapping_string_heavy(fn1, fn2, expected):
    logging.basicConfig(format='%(message)s', level=logging.CRITICAL)

    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)

    parent = Chem.MolFromMolFile(fn1, sanitize=False, removeHs=False)
    comp = Chem.MolFromMolFile(fn2, sanitize=False, removeHs=False)
    MC = MCS(parent, comp)
    assert MC.heavy_atom_match_list() == expected, 'Fail on heavy atom match list for ' + fn1 + ' ' + fn2


@pytest.mark.parametrize('fn1, fn2, expected', [
    (_rf('transforms/phenyl.sdf'), _rf('transforms/toluyl3.sdf'),
     "0:0,1:1,2:2,3:3,4:4,5:5,6:6,7:7,8:9,9:10,10:8,11:11,12:17,13:12,14:13,15:14,16:15,17:16"),
    (_rf('transforms/toluyl2.sdf'), _rf('transforms/chlorotoluyl1.sdf'),
     "0:0,1:1,2:2,3:3,4:4,5:5,6:6,7:8,8:9,9:10,10:7,11:11,12:12,13:13,14:14,15:15,16:16,17:17,18:18,19:19,20:20"),
    (_rf('transforms/toluyl3.sdf'), _rf('transforms/phenylfuran.sdf'),
     "0:0,1:1,2:2,3:3,4:4,5:5,6:6,7:7,9:13,10:14,11:15,12:17,13:18,14:19,15:20,16:21,17:16"),
    (_rf('transforms/toluyl.sdf'), _rf('transforms/phenylmethylamino.sdf'),
     "0:0,1:1,2:2,3:3,4:4,5:5,6:6,7:7,8:8,9:10,10:11,11:12,12:13,13:14,14:15,15:16,16:17,17:18,18:9,19:19,20:20"),
])
# Test getting the mapping string including hydrogens out of the MCS
def test_mapping_string_hydrogen(fn1, fn2, expected):
    logging.basicConfig(format='%(message)s', level=logging.CRITICAL)

    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)

    parent = Chem.MolFromMolFile(fn1, sanitize=False, removeHs=False)
    comp = Chem.MolFromMolFile(fn2, sanitize=False, removeHs=False)
    MC = MCS(parent, comp)
    assert MC.all_atom_match_list() == expected, 'Fail on all-atom match list for ' + fn1 + ' ' + fn2


@pytest.mark.parametrize('fn1, fn2, expected', [
    (_rf('chiral/Chiral1R.sdf'), _rf('chiral/Chiral1S.sdf'), 6),
    (_rf('chiral/Chiral1R.sdf'), _rf('chiral/Chiral2R.sdf'), 7),
    (_rf('chiral/Chiral1S.sdf'), _rf('chiral/Chiral2R.sdf'), 6),
    (_rf('chiral/Chiral3RS.sdf'), _rf('chiral/Chiral3SS.sdf'), 11),
    (_rf('chiral/Chiral3SR.sdf'), _rf('chiral/Chiral3SS.sdf'), 10),
    (_rf('chiral/Chiral3SR.sdf'), _rf('chiral/Chiral3RS.sdf'), 9),
    (_rf('chiral/Chiral4RR.sdf'), _rf('chiral/Chiral4RS.sdf'), 5),
    (_rf('chiral/RingChiralR.sdf'), _rf('chiral/RingChiralS.sdf'), 6),
    (_rf('chiral/SpiroR.sdf'), _rf('chiral/SpiroS.sdf'), 6),
    (_rf('chiral/bace_mk1.sdf'), _rf('chiral/bace_cat_13d.sdf'), 21),  # Bug found in BACE data set
    (_rf('chiral/bace_cat_13d.sdf'), _rf('chiral/bace_mk1.sdf'), 21),  # check both ways round
    (_rf('chiral/bace_mk1.sdf'), _rf('chiral/bace_cat_13d_perm1.sdf'), 21),
    # Check unaffected by atom order
    (_rf('chiral/bace_mk1.sdf'), _rf('chiral/bace_cat_13d_perm2.sdf'), 21),
    (_rf('chiral/bace_mk1.sdf'), _rf('chiral/bace_cat_13d_perm3.sdf'), 21),
    (_rf('chiral/bace_mk1.sdf'), _rf('chiral/bace_cat_13d_perm4.sdf'), 21),
    (_rf('chiral/bace_mk1.sdf'), _rf('chiral/bace_cat_13d_perm5.sdf'), 21),
    (_rf('chiral/bace_cat_13d_inverted.sdf'), _rf('chiral/bace_mk1.sdf'), 13),
    # Check that we do pick up the inverted case OK
    (_rf('chiral/bace_cat_13d.sdf'), _rf('chiral/bace_cat_13d_inverted.sdf'), 13)
    # Check normal vs inverted
])
# Test to check correct handling of chirality
def test_chirality_handling(fn1, fn2, expected):
    logging.basicConfig(format='%(message)s', level=logging.CRITICAL)

    lg = RDLogger.logger()
    lg.setLevel(RDLogger.INFO)

    parent = Chem.MolFromMolFile(fn1, sanitize=False, removeHs=False)
    comp = Chem.MolFromMolFile(fn2, sanitize=False, removeHs=False)
    MC = MCS(parent, comp, time=20, verbose='info', max3d=5, threed=True)
    assert MC.mcs_mol.GetNumHeavyAtoms() == expected, 'Fail on chiral MCS size for ' + fn1 + ' ' + fn2


# Test to check correct trimming of rings when 3D coordinate matching is used
def test_ring_trimming_on_3d_match():
    logging.basicConfig(format='%(message)s', level=logging.CRITICAL)
    parent=Chem.MolFromMolFile(_rf('transforms/phenylcyclopentylmethyl1.sdf'), sanitize=False, removeHs=False)
    comp=Chem.MolFromMolFile(_rf('transforms/phenylcyclopentylmethyl2.sdf'), sanitize=False, removeHs=False)
    MC=MCS(parent, comp, time=20, verbose='info', max3d=2, threed=True)
    assert MC.mcs_mol.GetNumHeavyAtoms() == 9, 'Fail on ring trim on 3D match'


# Test to check handling of the alpha- vs beta-naphthyl bug
def test_rdkit_broken_mcs_fix():
    logging.basicConfig(format='%(message)s', level=logging.CRITICAL)
    parent=Chem.MolFromMolFile(_rf('transforms/napthyl2.sdf'),sanitize=False, removeHs=False)
    comp=Chem.MolFromMolFile(_rf('transforms/napthyl3.sdf'),sanitize=False, removeHs=False)
    MC=MCS(parent, comp, time=20, verbose='info', max3d=0, threed=False)
    assert MC.mcs_mol.GetNumHeavyAtoms() < 25, 'Fail on detecting broken RDkit MCS on fused ring'


# Test to check handling of mapping to prochiral hydrogens
def test_mapping_prochiral_hydrogen():
    logging.basicConfig(format='%(message)s', level=logging.CRITICAL)
    parent=Chem.MolFromMolFile(_rf('chiral/tpbs2_lig1.sdf'), sanitize=False, removeHs=False)
    parent2=Chem.MolFromMolFile(_rf('chiral/tpbs2_lig1a.sdf'), sanitize=False, removeHs=False)   # Atom ordering changed
    comp=Chem.MolFromMolFile(_rf('chiral/tpbs2_lig2.sdf'), sanitize=False, removeHs=False)
    MC=MCS(parent, comp, time=20, verbose='info', max3d=3, threed=True)
    # Check that the correct prochiral hydrogen matches the bridging carbons
    assert "51:12" in MC.all_atom_match_list()
    assert "35:11" in MC.all_atom_match_list()
    MC=MCS(parent2, comp, time=20, verbose='info', max3d=3, threed=True)
    # parent2 is the same mol as parent1, except that atoms 34 and 35 were swapped
    assert "51:12" in MC.all_atom_match_list()
    assert "34:11" in MC.all_atom_match_list()


""" Problem is the executable moves around, so hard to test
def test_insufficient_arguments(self):
    cmd = [executable()]
    error_string = b'error: the following arguments are required: directory'
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    assert(error_string in stderr)
"""


def fields_for_link(mola, molb):
    """ Parse the out_score_with_connection.txt file, find the line for mola to molb, and return its fields. """
    with open('./out_score_with_connection.txt', 'r') as f:
        for line in f.readlines():
            fields = line.replace(",","").split()
            if (fields[2] == mola and fields[3] == molb) or (fields[3] == mola and fields[2] == molb):
                return fields
    pytest.fail("Output not found")


def score_for_link(mola, molb):
    return float(fields_for_link(mola,molb)[4])


@pytest.fixture
def in_tmpdir(tmpdir):
    with tmpdir.as_cwd():
        yield


def test_complete_run(in_tmpdir):
    logging.basicConfig(format='%(message)s', level=logging.CRITICAL)

    #progname=sys.argv[0]
    #sys.argv=[progname,'-o','--output-no-images','--output-no-graph','test/linksfile']
    #dbmol.startup()
    dbmol._startup_inner(
        output=True,
        directory=_rf('linksfile/'),
        output_no_images=True,
        output_no_graph=True,
    )

    # Check scores
    assert score_for_link('phenyl.sdf', 'phenylcyclobutyl.sdf') == pytest.approx(0.67032)
    assert score_for_link('phenyl.sdf', 'phenylfuran.sdf') == pytest.approx(0.60653)
    assert score_for_link('phenyl.sdf', 'toluyl.sdf') == pytest.approx(0.90484)
    assert score_for_link('phenylcyclobutyl.sdf', 'phenylfuran.sdf') == pytest.approx(0.40657)
    assert score_for_link('phenylcyclobutyl.sdf', 'toluyl.sdf') == pytest.approx(0.33287)
    assert score_for_link('phenylfuran.sdf', 'toluyl.sdf') == pytest.approx(0.54881)
    # Check connections
    assert fields_for_link('phenyl.sdf', 'phenylcyclobutyl.sdf')[7] == "Yes"
    assert fields_for_link('phenyl.sdf', 'phenylfuran.sdf')[7] == "No"
    assert fields_for_link('phenyl.sdf', 'toluyl.sdf')[7] == "Yes"
    assert fields_for_link('phenylcyclobutyl.sdf', 'phenylfuran.sdf')[7] == "Yes"
    assert fields_for_link('phenylcyclobutyl.sdf', 'toluyl.sdf')[7] == "No"
    assert fields_for_link('phenylfuran.sdf', 'toluyl.sdf')[7] == "Yes"


def test_complete_run_parallel(in_tmpdir):
    '''Test running in parallel mode with 5 subprocesses.'''
    logging.basicConfig(format='%(message)s', level=logging.CRITICAL)
    #progname=sys.argv[0]
    #sys.argv=[progname,'-o','-p','5','--output-no-images','--output-no-graph','test/linksfile']
    #dbmol.startup()
    dbmol._startup_inner(
        directory=_rf('linksfile/'),
        output=True,
        parallel=5,
        output_no_images=True,
        output_no_graph=True,
    )

    # Check scores
    assert score_for_link('phenyl.sdf', 'phenylcyclobutyl.sdf') == pytest.approx(0.67032)
    assert score_for_link('phenyl.sdf', 'phenylfuran.sdf') == pytest.approx(0.60653)
    assert score_for_link('phenyl.sdf', 'toluyl.sdf') == pytest.approx(0.90484)
    assert score_for_link('phenylcyclobutyl.sdf', 'phenylfuran.sdf') == pytest.approx(0.40657)
    assert score_for_link('phenylcyclobutyl.sdf', 'toluyl.sdf') == pytest.approx(0.33287)
    assert score_for_link('phenylfuran.sdf', 'toluyl.sdf') == pytest.approx(0.54881)
    # Check connections
    assert fields_for_link('phenyl.sdf','phenylcyclobutyl.sdf')[7] == "Yes"
    assert fields_for_link('phenyl.sdf','phenylfuran.sdf')[7] == "No"
    assert fields_for_link('phenyl.sdf','toluyl.sdf')[7] == "Yes"
    assert fields_for_link('phenylcyclobutyl.sdf','phenylfuran.sdf')[7] == "Yes"
    assert fields_for_link('phenylcyclobutyl.sdf','toluyl.sdf')[7] == "No"
    assert fields_for_link('phenylfuran.sdf','toluyl.sdf')[7] == "Yes"


def test_linksfile(in_tmpdir):
    """ Test a linksfile forcing a link from phenyl to phenylfuran. """
    logging.basicConfig(format='%(message)s', level=logging.CRITICAL)

    #progname=sys.argv[0]
    #sys.argv=[progname,'-o','--output-no-images','--output-no-graph','--links-file','test/linksfile/links1.txt','test/linksfile']
    #dbmol.startup()
    dbmol._startup_inner(
        directory=_rf('linksfile/'),
        output=True,
        output_no_images=True,
        output_no_graph=True,
        links_file=_rf('linksfile/links1.txt'),
    )

    # Check scores
    assert score_for_link('phenyl.sdf', 'phenylcyclobutyl.sdf') == pytest.approx(0.67032)
    assert score_for_link('phenyl.sdf', 'phenylfuran.sdf') == pytest.approx(0.60653)
    assert score_for_link('phenyl.sdf', 'toluyl.sdf') == pytest.approx(0.90484)
    assert score_for_link('phenylcyclobutyl.sdf', 'phenylfuran.sdf') == pytest.approx(0.40657)
    assert score_for_link('phenylcyclobutyl.sdf', 'toluyl.sdf') == pytest.approx(0.33287)
    assert score_for_link('phenylfuran.sdf', 'toluyl.sdf') == pytest.approx(0.54881)
    # Check connections
    assert fields_for_link('phenyl.sdf','phenylcyclobutyl.sdf')[7] == "Yes"
    assert fields_for_link('phenyl.sdf','phenylfuran.sdf')[7] == "Yes"
    assert fields_for_link('phenyl.sdf','toluyl.sdf')[7] == "Yes"
    assert fields_for_link('phenylcyclobutyl.sdf','phenylfuran.sdf')[7] == "Yes"
    assert fields_for_link('phenylcyclobutyl.sdf','toluyl.sdf')[7] == "No"
    assert fields_for_link('phenylfuran.sdf','toluyl.sdf')[7] == "Yes"


def test_linksfile_scores(in_tmpdir):
    """ Test a linksfile forcing prespecified scores for some links."""
    logging.basicConfig(format='%(message)s', level=logging.CRITICAL)

    #progname=sys.argv[0]
    #sys.argv=[progname,'-o','--output-no-images','--output-no-graph','--links-file','test/linksfile/links2.txt','test/linksfile']
    #dbmol.startup()
    dbmol._startup_inner(
        directory=_rf('linksfile/'),
        output=True,
        output_no_images=True,
        output_no_graph=True,
        links_file=_rf('linksfile/links2.txt'),
    )

    # Check scores
    assert score_for_link('phenyl.sdf', 'phenylcyclobutyl.sdf') == pytest.approx(0.67032)
    assert score_for_link('phenyl.sdf', 'phenylfuran.sdf') == pytest.approx(0.77777)  # Forced from linksfile
    assert score_for_link('phenyl.sdf', 'toluyl.sdf') == pytest.approx(0.88888)  # Forced from linksfile
    assert score_for_link('phenylcyclobutyl.sdf', 'phenylfuran.sdf') == pytest.approx(0.40657)
    assert score_for_link('phenylcyclobutyl.sdf', 'toluyl.sdf') == pytest.approx(0.33287)
    assert score_for_link('phenylfuran.sdf', 'toluyl.sdf') == pytest.approx(0.54881)
    # Check connections
    assert fields_for_link('phenyl.sdf', 'phenylcyclobutyl.sdf')[7] == "Yes"
    assert fields_for_link('phenyl.sdf', 'phenylfuran.sdf')[7] == "No"
    assert fields_for_link('phenyl.sdf', 'toluyl.sdf')[7] == "Yes"
    assert fields_for_link('phenylcyclobutyl.sdf', 'phenylfuran.sdf')[7] == "Yes"
    assert fields_for_link('phenylcyclobutyl.sdf', 'toluyl.sdf')[7] == "No"
    assert fields_for_link('phenylfuran.sdf', 'toluyl.sdf')[7] == "Yes"


def test_linksfile_scores_force(in_tmpdir):
    """ Test a linksfile forcing prespecified scores and link inclusion for some links."""
    logging.basicConfig(format='%(message)s', level=logging.CRITICAL)

    #progname=sys.argv[0]
    #sys.argv=[progname,'-o','--output-no-images','--output-no-graph','--links-file','test/linksfile/links3.txt','test/linksfile']
    #dbmol.startup()
    dbmol._startup_inner(
        directory=_rf('linksfile/'),
        output=True,
        output_no_images=True,
        output_no_graph=True,
        links_file=_rf('linksfile/links3.txt'),
    )

    # Check scores
    assert score_for_link('phenyl.sdf', 'phenylcyclobutyl.sdf') == pytest.approx(0.1)
    assert score_for_link('phenyl.sdf', 'phenylfuran.sdf') == pytest.approx(0.2)
    assert score_for_link('phenyl.sdf', 'toluyl.sdf') == pytest.approx(0.3)
    assert score_for_link('phenylcyclobutyl.sdf', 'phenylfuran.sdf') == pytest.approx(0.4)
    assert score_for_link('phenylcyclobutyl.sdf', 'toluyl.sdf') == pytest.approx(0.5)
    assert score_for_link('phenylfuran.sdf', 'toluyl.sdf') == pytest.approx(0.6)
    # Check connections
    assert fields_for_link('phenyl.sdf','phenylcyclobutyl.sdf')[7] == "Yes"
    assert fields_for_link('phenyl.sdf','phenylfuran.sdf')[7] == "Yes"
    assert fields_for_link('phenyl.sdf','toluyl.sdf')[7] == "No"
    assert fields_for_link('phenylcyclobutyl.sdf','phenylfuran.sdf')[7] == "Yes"
    assert fields_for_link('phenylcyclobutyl.sdf','toluyl.sdf')[7] == "Yes"
    assert fields_for_link('phenylfuran.sdf','toluyl.sdf')[7] == "Yes"


def test_no_cycle_cover(in_tmpdir):
    logging.basicConfig(format='%(message)s', level=logging.CRITICAL)

    #progname=sys.argv[0]
    #sys.argv=[progname,'-o','-T','--output-no-images','--output-no-graph','test/linksfile']
    #dbmol.startup()
    dbmol._startup_inner(
        directory=_rf('linksfile/'),
        output=True,
        allow_tree=True,
        output_no_images=True,
        output_no_graph=True,
    )
    # Check connections
    assert fields_for_link('phenyl.sdf','phenylcyclobutyl.sdf')[7] == "Yes"
    assert fields_for_link('phenyl.sdf','phenylfuran.sdf')[7] == "Yes"
    assert fields_for_link('phenyl.sdf','toluyl.sdf')[7] == "Yes"
    assert fields_for_link('phenylcyclobutyl.sdf','phenylfuran.sdf')[7] == "No"
    assert fields_for_link('phenylcyclobutyl.sdf','toluyl.sdf')[7] == "No"
    assert fields_for_link('phenylfuran.sdf','toluyl.sdf')[7] == "No"
