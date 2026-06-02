from __future__ import annotations

import math
from collections import defaultdict

from rdkit import Chem

try:
    from gufe import LigandAtomMapping
except ImportError:
    pass

from lomap import dbmol as _dbmol
from lomap import mcs as lomap_mcs
from lomap.utils import requires_package

DEFAULT_ANS_DIFFICULTY = {
    # H to element - not sure this has any effect currently
    1: {9: 0.5, 17: 0.25, 35: 0, 53: -0.5},
    # O to element - methoxy to Cl/Br is easier than expected
    8: {17: 0.85, 35: 0.85},
    # F to element
    9: {17: 0.5, 35: 0.25, 53: 0},
    # Cl to element
    17: {35: 0.85, 53: 0.65},
    # Br to element
    35: {53: 0.85},
}


@requires_package("gufe")
def ecr_score(mapping: LigandAtomMapping, charge_changes_score) -> float:
    """Equal charge rule (ECR) score.

    Returns 1.0 if both molecules have the same formal charge,
    otherwise returns ``charge_changes_score``.

    Parameters
    ----------
    mapping : LigandAtomMapping
      Mapping between the two ligands in the edge.
    charge_changes_score : float
      Score assigned when the two molecules differ in net formal charge.

    Returns
    -------
    score : float
      Value in the range [0, 1.0].
    """
    molA = mapping.componentA.to_rdkit()
    molB = mapping.componentB.to_rdkit()

    # Get formal charges of both molecules
    total_charge_molA = _dbmol.formal_charge(molA)
    total_charge_molB = _dbmol.formal_charge(molB)
    if abs(total_charge_molB - total_charge_molA) < 1e-3:
        return _dbmol.ecr(molA, molB)
    else:
        return charge_changes_score


@requires_package("gufe")
def mcsr_score(mapping: LigandAtomMapping, beta: float = 0.1) -> float:
    """Maximum common substructure rule (MCSR) score.

    This rule is defined as:

    .. math::

        mcsr = exp( - beta * (n1 + n2 - 2 * n_common))

    Where n1 and n2 are the number of heavy atoms in each molecule, and
    n_common is the number of heavy atoms in the MCS. This makes the term
    ``n1 + n2 - 2 * n_common`` the total number of atoms inserted or
    deleted in the transformation.

    The exponential is used to ensure the score ranges between 0 and 1,
    and to strongly favor small structural changes.

    Parameters
    ----------
    mapping : LigandAtomMapping
      Mapping between the two ligands in the edge.
    beta : float, default 0.1
      Scaling factor.

    Returns
    -------
    score : float
      Value in the range [0, 1.0], with 1.0 indicating complete overlap.
    """
    molA = mapping.componentA.to_rdkit()
    molB = mapping.componentB.to_rdkit()
    molA_to_molB = mapping.componentA_to_componentB

    n1 = molA.GetNumHeavyAtoms()
    n2 = molB.GetNumHeavyAtoms()
    # get heavy atom mcs count
    n_common = 0
    for i, j in molA_to_molB.items():
        if (
            molA.GetAtomWithIdx(i).GetAtomicNum() != 1
            and molB.GetAtomWithIdx(j).GetAtomicNum() != 1
        ):
            n_common += 1

    mcsr = math.exp(-beta * (n1 + n2 - 2 * n_common))

    return mcsr


@requires_package("gufe")
def mncar_score(mapping: LigandAtomMapping, ths: int = 4) -> float:
    """Minimum number of common atoms rule (MNCAR) score.

    The two molecules must share at least ``ths`` heavy atoms to be regarded
    as similar. Returns 1.0 if this condition is met, or if either molecule
    is small (fewer than ``ths + 3`` heavy atoms), otherwise returns 0.0.

    Parameters
    ----------
    mapping : LigandAtomMapping
      Mapping between the two ligands in the edge.
    ths : int, default 4
      Minimum number of common heavy atoms required, default 4.

    Returns
    -------
    score : float
      1.0 if the constraint is satisfied, 0.0 otherwise.
    """
    molA = mapping.componentA.to_rdkit()
    molB = mapping.componentB.to_rdkit()
    molA_to_molB = mapping.componentA_to_componentB

    n1 = molA.GetNumHeavyAtoms()
    n2 = molB.GetNumHeavyAtoms()
    n_common = 0
    for i, j in molA_to_molB.items():
        if (
            molA.GetAtomWithIdx(i).GetAtomicNum() != 1
            and molB.GetAtomWithIdx(j).GetAtomicNum() != 1
        ):
            n_common += 1

    ok = (n_common > ths) or (n1 < ths + 3) or (n2 < ths + 3)

    return 1.0 if ok else 0.0


@requires_package("gufe")
def tmcsr_score(mapping: LigandAtomMapping):
    """Trimmed maximum common substructure rule (TMCSR) score.

    .. warning::
       Not yet implemented.
    """
    raise NotImplementedError


@requires_package("gufe")
def atomic_number_score(
    mapping: LigandAtomMapping,
    beta=0.1,
    difficulty: dict[int, dict[int, float]] | None = None
) -> float:
    """A score on the elemental changes happening in the mapping

    For each transmuted atom, a mismatch score is summed, according to the
    difficulty scores (see difficult parameter). The final score is then
    given as:

    .. math::

        score = exp(-beta * mismatch)

    Parameters
    ----------
    mapping : LigandAtomMapping
      Mapping between the two ligands in the edge.
    beta : float, default 0.1
      Scaling factor for this rule, default 0.1
    difficulty : dict[int, dict[int, float] | None, default None
      A dict of dicts, mapping atomic number of one species, to another,
      to a mismatch in the identity of these elements.  1.0 indicates two
      elements are considered interchangeable, 0.0 indicates two elements
      are incompatible, a default of 0.5 is used.
      The scores in openfe.setup.lomap_mapper.DEFAULT_ANS_DIFFICULT are
      used by default

    Returns
    -------
    score : float
    """
    molA = mapping.componentA.to_rdkit()
    molB = mapping.componentB.to_rdkit()
    molA_to_molB = mapping.componentA_to_componentB

    if difficulty is None:
        difficulty = DEFAULT_ANS_DIFFICULTY

    nmismatch = 0
    for i, j in molA_to_molB.items():
        atom_i = molA.GetAtomWithIdx(i)
        atom_j = molB.GetAtomWithIdx(j)

        n_i = atom_i.GetAtomicNum()
        n_j = atom_j.GetAtomicNum()

        if n_i == n_j:
            continue
        elif n_i == 1 or n_j == 1:  # ignore hydrogen switches?
            continue

        try:
            ij = difficulty[n_i][n_j]
        except KeyError:
            ij = -1
        try:
            ji = difficulty[n_j][n_i]
        except KeyError:
            ji = -1
        diff = max(ij, ji)
        if diff == -1:
            diff = 0.5

        nmismatch += 1 - diff

    return math.exp(-beta * nmismatch)


@requires_package("gufe")
def hybridization_score(mapping: LigandAtomMapping, beta=0.15) -> float:
    """Hybridization score — penalizes atom hybridization mismatches in the mapping.

    For each mapped heavy atom pair with differing hybridization states,
    a mismatch is counted. N sp3/sp2 interchanges are permitted. The final
    score is:

    .. math::

        score = exp(-beta * nmismatch)

    Parameters
    ----------
    mapping : LigandAtomMapping
      Mapping between the two ligands in the edge.
    beta : float, default 0.15
      Scaling factor.

    Returns
    -------
    score : float
      Value in the range [0, 1.0], with 1.0 indicating no hybridization
      mismatches.
    """
    mol1 = mapping.componentA.to_rdkit()
    mol2 = mapping.componentB.to_rdkit()
    molA_to_molB = mapping.componentA_to_componentB

    nmismatch = 0
    for i, j in molA_to_molB.items():
        atom_i = mol1.GetAtomWithIdx(i)
        atom_j = mol2.GetAtomWithIdx(j)

        if atom_i.GetAtomicNum() == 1 or atom_j.GetAtomicNum() == 1:
            # skip hydrogen changes
            continue

        hyb_i = lomap_mcs.atom_hybridization(atom_i)
        hyb_j = lomap_mcs.atom_hybridization(atom_j)

        mismatch = hyb_i != hyb_j
        # Allow Nsp3 to match Nsp2, otherwise guanidines etc become painful
        if (
            atom_i.GetAtomicNum() == 7
            and atom_j.GetAtomicNum() == 7
            and hyb_i in [2, 3]
            and hyb_j in [2, 3]
        ):
            mismatch = False

        if mismatch:
            nmismatch += 1

    return math.exp(-beta * nmismatch)


@requires_package("gufe")
def sulfonamides_score(mapping: LigandAtomMapping, beta=0.4) -> float:
    """Sulfonamide score — penalizes mappings that mutate a sulfonamide group in or out.

    Testing has shown that growing a sulfonamide from scratch performs very
    badly. Returns ``math.exp(-beta)`` if a sulfonamide group appears in the
    unmapped remainder of either molecule, otherwise 1.0.

    Parameters
    ----------
    mapping : LigandAtomMapping
      Mapping between the two ligands in the edge.
    beta : float, default 0.4
      Scaling factor controlling the size of the penalty. Smaller values give
      larger penalties.

    Returns
    -------
    score : float
      ``math.exp(-beta)`` if a sulfonamide is mutated in or out, else 1.0.
    """
    molA = mapping.componentA.to_rdkit()
    molB = mapping.componentB.to_rdkit()
    molA_to_molB = mapping.componentA_to_componentB

    def has_sulfonamide(mol):
        return mol.HasSubstructMatch(Chem.MolFromSmarts("S(=O)(=O)N"))

    # create "remainders" of both molA and molB
    remA = Chem.EditableMol(molA)
    # this incremental deletion only works when we go from high to low,
    # as atoms are reindexed as we delete
    for i, j in sorted(molA_to_molB.items(), reverse=True):
        if molA.GetAtomWithIdx(i).GetAtomicNum() != molB.GetAtomWithIdx(j).GetAtomicNum():
            continue
        remA.RemoveAtom(i)
    # loop molB separately, sorted by A indices doesn't necessarily sort
    # the B indices too, so these loops are in different orders
    remB = Chem.EditableMol(molB)
    for i, j in sorted(molA_to_molB.items(), key=lambda x: x[1], reverse=True):
        if molA.GetAtomWithIdx(i).GetAtomicNum() != molB.GetAtomWithIdx(j).GetAtomicNum():
            continue
        remB.RemoveAtom(j)

    if has_sulfonamide(remA.GetMol()) or has_sulfonamide(remB.GetMol()):
        return math.exp(-beta)
    else:
        return 1.0


@requires_package("gufe")
def heterocycles_score(mapping: LigandAtomMapping, beta=0.4) -> float:
    """Heterocycle score — penalizes mappings that form a heterocycle from a hydrogen.

    Returns ``math.exp(-beta)`` if a heterocycle is formed from a hydrogen.
    Testing has shown that growing a pyridine or other heterocycle
    is unlikely to work (better to grow phenyl than mutate).

    Parameters
    ----------
    mapping : LigandAtomMapping
      Mapping between the two ligands in the edge.
    beta : float, default 0.4
      Scaling factor controlling the size of the penalty.

    Returns
    -------
    score : float
      ``math.exp(-beta)`` if a disallowed heterocycle is formed, else 1.0.
    """
    molA = mapping.componentA.to_rdkit()
    molB = mapping.componentB.to_rdkit()
    molA_to_molB = mapping.componentA_to_componentB

    def creates_heterocyle(mol):
        # these patterns are lifted from lomap2 repo
        return (
            mol.HasSubstructMatch(Chem.MolFromSmarts("[n]1[c,n][c,n][c,n][c,n][c,n]1"))
            or mol.HasSubstructMatch(Chem.MolFromSmarts("[o,n,s]1[n][c,n][c,n][c,n]1"))
            or mol.HasSubstructMatch(Chem.MolFromSmarts("[o,n,s]1[c,n][n][c,n][c,n]1"))
        )

    # create "remainders" of both molA and molB
    # create "remainders" of both molA and molB
    remA = Chem.EditableMol(molA)
    # this incremental deletion only works when we go from high to low,
    # as atoms are reindexed as we delete
    for i, j in sorted(molA_to_molB.items(), reverse=True):
        if molA.GetAtomWithIdx(i).GetAtomicNum() != molB.GetAtomWithIdx(j).GetAtomicNum():
            continue
        remA.RemoveAtom(i)
    # loop molB separately, sorted by A indices doesn't necessarily sort
    # the B indices too, so these loops are in different orders
    remB = Chem.EditableMol(molB)
    for i, j in sorted(molA_to_molB.items(), key=lambda x: x[1], reverse=True):
        if molA.GetAtomWithIdx(i).GetAtomicNum() != molB.GetAtomWithIdx(j).GetAtomicNum():
            continue
        remB.RemoveAtom(j)

    if creates_heterocyle(remA.GetMol()) or creates_heterocyle(remB.GetMol()):
        return math.exp(-beta)
    else:
        return 1.0


@requires_package("gufe")
def transmuting_methyl_into_ring_score(mapping: LigandAtomMapping, beta=0.1, penalty=6.0) -> float:
    """Penalises having a non-mapped ring atoms become a non-ring

    This score would for example penalise R-CH3 to R-Ph where R is the same
    mapped atom and both CH3 and Ph are unmapped. Does not penalise R-H to R-Ph.
    If any atoms trigger the rule returns a score of:

    .. math::

        exp(-1 * beta * penalty)

    Parameters
    ----------
    mapping : LigandAtomMapping
      Mapping between the two ligands in the edge.
    beta : float
      Score scaling factor.
    penalty : float
      Score scaling factor.

    Returns
    -------
    score : float
    """
    molA = mapping.componentA.to_rdkit()
    molB = mapping.componentB.to_rdkit()
    # filter out hydrogens from mapping
    # this matches the current implementation in Lomap
    molA_to_molB = {}
    for i, j in mapping.componentA_to_componentB.items():
        if molA.GetAtomWithIdx(i).GetAtomicNum() == 1:
            continue
        if molB.GetAtomWithIdx(j).GetAtomicNum() == 1:
            continue
        molA_to_molB[i] = j

    ringbreak = False
    for i, j in molA_to_molB.items():
        atomA = molA.GetAtomWithIdx(i)

        for bA in atomA.GetBonds():
            otherA = bA.GetOtherAtom(atomA)
            if otherA.GetAtomicNum() == 1:
                continue
            if otherA.GetIdx() in molA_to_molB:
                # if other end of bond in core, ignore
                continue

            # try and find the corresponding atom in molecule B
            atomB = molB.GetAtomWithIdx(j)
            for bB in atomB.GetBonds():
                otherB = bB.GetOtherAtom(atomB)
                if otherB.GetAtomicNum() == 1:
                    continue
                if otherB.GetIdx() in molA_to_molB.values():
                    continue

                if otherA.IsInRing() ^ otherB.IsInRing():
                    ringbreak = True

    if not ringbreak:
        return 1.0
    else:
        return math.exp(-beta * penalty)


@requires_package("gufe")
def transmuting_ring_sizes_score(mapping: LigandAtomMapping) -> float:
    """Ring size score — penalizes mappings that alter a ring size.

    Checks first-degree neighbors of mapped atoms; if a non-mapped neighbor
    is in a ring in both molecules but the ring sizes differ, a value of 0.1
    is returned, otherwise 1.0 is returned.

    Parameters
    ----------
    mapping : LigandAtomMapping
      Mapping between the two ligands in the edge.

    Returns
    -------
    score : float
      0.1 if any ring size change is detected, else 1.0.
    """
    molA = mapping.componentA.to_rdkit()
    molB = mapping.componentB.to_rdkit()
    molA_to_molB = mapping.componentA_to_componentB

    def gen_ringdict(mol):
        # maps atom idx to ring sizes
        ringinfo = mol.GetRingInfo()

        idx_to_ringsizes = defaultdict(list)
        for r in ringinfo.AtomRings():
            for idx in r:
                idx_to_ringsizes[idx].append(len(r))
        return idx_to_ringsizes

    # generate ring size dicts
    ringdictA = gen_ringdict(molA)
    ringdictB = gen_ringdict(molB)

    is_bad = False
    # check first degree neighbours of core atoms to see if their ring
    # sizes are the same
    for i, j in molA_to_molB.items():
        atomA = molA.GetAtomWithIdx(i)

        for bA in atomA.GetBonds():
            otherA = bA.GetOtherAtom(atomA)
            if otherA.GetIdx() in molA_to_molB:
                # if other end of bond in core, ignore
                continue
            # otherA is an atom not in the mapping, but bonded to an
            # atom in the mapping
            if not otherA.IsInRing():
                continue

            # try and find the corresponding atom in molecule B
            atomB = molB.GetAtomWithIdx(j)
            for bB in atomB.GetBonds():
                otherB = bB.GetOtherAtom(atomB)
                if otherB.GetIdx() in molA_to_molB.values():
                    continue
                if not otherB.IsInRing():
                    continue

                # ringdict[idx] will give the list of ringsizes for an atom
                if set(ringdictA[otherA.GetIdx()]) != set(ringdictB[otherB.GetIdx()]):
                    is_bad = True

    return 0.1 if is_bad else 1.0


@requires_package("gufe")
def default_lomap_score(mapping: LigandAtomMapping, charge_changes_score=0.1) -> float:
    """The default score function from Lomap2


    This score is a combination of many rules combined and considers factors such as the
    number of heavy atoms in common, if ring sizes are changed or rings are broken,
    or if other alchemically unwise transformations are attempted.

    Parameters
    ----------
    mapping : LigandAtomMapping
      Mapping between the two ligands in the edge.
    charge_changes_score: float
      The electrostatic score to be assigned for mappings of ligands that
      differ in net charge.
      Default: 0.1 (e.g. allowing net charge changes)

    Returns
    -------
    score : float
       A rating of how good this mapping is, from 0.0 (terrible) to 1.0 (great).
    """
    score = math.prod(
        (
            ecr_score(mapping, charge_changes_score),
            mncar_score(mapping),
            mcsr_score(mapping),
            atomic_number_score(mapping),
            hybridization_score(mapping),
            sulfonamides_score(mapping),
            heterocycles_score(mapping),
            transmuting_methyl_into_ring_score(mapping),
            transmuting_ring_sizes_score(mapping),
        )
    )

    return score
