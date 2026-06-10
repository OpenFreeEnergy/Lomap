"""

LOMAP: Maximum Common Subgraph and scoring calculations
=====

Alchemical free energy calculations hold increasing promise as an aid to drug
discovery efforts. However, applications of these techniques in discovery
projects have been relatively few, partly because of the difficulty of planning
and setting up calculations. The Lead Optimization Mapper (LOMAP) is an
automated algorithm to plan efficient relative free energy calculations between
potential ligands within a substantial set of compounds.

"""

import logging
import math
from typing import Iterator, Literal
import warnings

from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Draw, rdFMCS, rdmolops
from rdkit.Chem.Draw.MolDrawing import DrawingOptions
from rdkit.Geometry.rdGeometry import Point3D

__all__ = ["MCS"]

logger = logging.getLogger(__name__)


def atom_hybridization(a: Chem.Atom) -> int:
    """
    RDKit has an un-useful hybridization definition. Instead, just look at the number
    of multiple bonds from an atom.

    Parameters
    ----------
    a : Chem.Atom
      Atom to get the hybridization for.

    Returns
    -------
    int
      The hybridization of the atom.
    """
    if a.GetIsAromatic():
        return 2

    xs: float = 0
    for b in a.GetBonds():
        if b.GetBondType() == Chem.rdchem.BondType.AROMATIC:
            return 2
        if b.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            xs += 1
        if b.GetBondType() == Chem.rdchem.BondType.TRIPLE:
            xs += 2
        if b.GetBondType() == Chem.rdchem.BondType.ONEANDAHALF:
            xs += 0.5

    # O- is sp2 to avoid problems with carboxylate etc
    if a.GetAtomicNum() == 8 and a.GetFormalCharge() < 0:
        return 2  # sp2

    if xs == 0:
        return 3  # sp3
    if xs > 1.1:
        return 1  # sp
    return 2  # sp2


def substructure_centre(mol: Chem.Mol, mol_sub: tuple[int, ...]) -> Point3D:
    """

    This function takes a molecule and a list of atom indices
    in that molecule and returns an RDKit Point3D representing
    the geometric centre of the atoms in the list

    Parameters
    ----------
    mol: Chem.Mol
      Molecule to get the substructure atoms from.
    mol_sub : tuple[int, ...]
      Atom indices for the substructure.

    Returns
    -------
    Point3D
      Geometric center of the substructure atoms.
    """
    s = Point3D()
    for i in mol_sub:
        s += mol.GetConformer().GetAtomPosition(i)
    return s / len(mol_sub)


class MCS:
    """

    This class is used to compute the Maximum Common Subgraph (MCS) between two
    RDkit molecule objects and to score their similarity by using defined rules

    """

    # Need to declare type to make mypy happy
    mcs_mol: Chem.Mol
    mcs_mol_smarts: str
    beta: float
    moli: Chem.Mol
    molj: Chem.Mol
    _moli_noh: Chem.Mol
    _molj_noh: Chem.Mol
    _map_moli_noh: dict[int, int]
    _map_molj_noh: dict[int, int]
    _map_moli_molj: list[tuple[int, int]]
    options: dict[str, int | str | float | bool]

    def __init__(
        self,
        moli: Chem.Mol,
        molj: Chem.Mol,
        time: int = 20,
        verbose: Literal["debug", "info", "warning", "error", "critical"] = "info",
        max3d: float = 1000.0,
        threed: bool = False,
        element_change: bool = True,
        seed: str = "",
        shift: bool = True,
    ) -> None:
        """
        Initialization function

        Parameters
        ----------
        moli : Chem.Mol
          The first molecule used to perform the MCS calculation.
        molj : Chem.Mol
          The second molecule used to perform the MCS calculation.
        time : int, default 20
          Timeout of MCS algorithm in seconds, passed to RDKit.
        verbose : Literal["debug", "info", "warning", "error", "critical"], default 'info'
          Logging level, defines verbosity of generated logs.
        max3d : float, default 1000.0
          The MCS is trimmed to remove atoms which are further apart than
          this distance (in units of Angstrom). The default, 1000.0, essentially
          tries not to trim.
        threed : bool, default False
          When disambiguating the substructure found back to the original
          molecules, if ``True`` 3D coordinates are used, otherwise the number
          of elemental changes is minimised.
        element_change : bool, default True
          If ``True``, allow elemental changes in mappings.
        seed : string, default ""
          SMARTS string to use as seed for MCS searches.  When used across an
          entire set of ligands, this can speed up calculations considerably.
        shift : bool, default True
          If ``True``, when ``threed`` is also ``True``, translate the
          molecules' coordinates to maximise 3D overlap.

        .. versionchanged:: 2.1.0
           Added element_change kwarg
        .. versionchanged:: 2.2.0
           Added seed option
        .. versionchanged:: 2.3.0
           Added shift option
        """
        self.options = {
            "time": time,
            "verbose": verbose,
            "max3d": max3d,
            "threed": threed,
            "element_change": element_change,
            "seed": seed,
            "shift": shift,
        }

        def best_substruct_match_to_mcs(
            moli: Chem.Mol,
            molj: Chem.Mol,
            by_rmsd: bool,
            use_shift: bool,
        ) -> tuple[tuple[int, ...], tuple[int, ...]]:
            """
            This function looks over all of the substructure matches and returns the one
            with the best 3D correspondence (if ``by_rmsd`` is ``True``), or the fewest number
            of atomic number mismatches (if ``by_rmsd`` is ``False``).

            Parameters
            ----------
            moli, molj : Chem.Mol
              Molecules to get the substructure match for.
            by_rmsd : bool
              If ``True``, prune substructure matches by 3D correspondence.
              If ``False``, prune substructure by the number of atomic mismatches.
            use_shift : bool
              If ``True``, translate the 3D coordinates for the substructure by
              center of geometry distance. Note that this does not do any rotation.

            Returns
            -------
            besti, bestj : tuple[int, ...]
              Best substructure matches for ``moli`` and ``molj``.
            """

            # Sanity checking
            if not moli.HasSubstructMatch(self.mcs_mol):
                raise ValueError("RDkit MCS Subgraph first molecule search failed")

            if not molj.HasSubstructMatch(self.mcs_mol):
                raise ValueError("RDkit MCS Subgraph second molecule search failed")

            moli_sub = moli.GetSubstructMatches(self.mcs_mol, uniquify=False)
            molj_sub = molj.GetSubstructMatches(self.mcs_mol, uniquify=False)
            best_rmsd = float("inf")
            besti: tuple[int, ...] = moli_sub[0]
            bestj: tuple[int, ...] = molj_sub[0]
            for mapi in moli_sub:
                for mapj in molj_sub:
                    # Compute the translation to bring molj's centre over moli
                    if by_rmsd and use_shift:
                        coord_delta = substructure_centre(moli, mapi) - substructure_centre(
                            molj, mapj
                        )
                    else:
                        coord_delta = Point3D(0.0, 0.0, 0.0)
                    rmsd: float = 0
                    for pair in zip(mapi, mapj):
                        if by_rmsd:
                            rmsd += (
                                moli.GetConformer().GetAtomPosition(pair[0])
                                - molj.GetConformer().GetAtomPosition(pair[1])
                                - coord_delta
                            ).LengthSq()
                        elif (
                            moli.GetAtomWithIdx(pair[0]).GetAtomicNum()
                            != molj.GetAtomWithIdx(pair[1]).GetAtomicNum()
                        ):
                            rmsd += 1
                    if rmsd < best_rmsd:
                        besti = mapi
                        bestj = mapj
                        best_rmsd = rmsd

            return besti, bestj

        def trim_mcs_mol(max_deviation: float, use_shift: bool) -> None:
            """

            This function is used to trim the MCS molecule to remove mismatched atoms i.e atoms
            where the topological mapping does not work in 3D coordinates.

            The sets of mapped atoms are optionally translated to bring their geometric centres
            into alignment before trimming

            Parameters
            ----------
            max_deviation : float
              The maximum difference in Angstroms between mapped atoms to allow.
            use_shift : bool
              If ``True``, translate the two molecules by the center of geometry
              of the substructure match to maximise overlap before calculating distances.
            """

            while True:
                mapi, mapj = best_substruct_match_to_mcs(
                    self._moli_noh, self._molj_noh, by_rmsd=True, use_shift=use_shift
                )
                # Compute the translation to bring molj's centre over moli
                if shift:
                    coord_delta = substructure_centre(self._moli_noh, mapi) - substructure_centre(
                        self._molj_noh, mapj
                    )
                else:
                    coord_delta = Point3D(0.0, 0.0, 0.0)
                worstatomidx = -1
                worstdist: float = 0
                atomidx = 0
                for pair in zip(mapi, mapj):
                    dist = (
                        self._moli_noh.GetConformer().GetAtomPosition(pair[0])
                        - self._molj_noh.GetConformer().GetAtomPosition(pair[1])
                        - coord_delta
                    ).Length()
                    if dist > worstdist:
                        worstdist = dist
                        worstatomidx = atomidx
                    atomidx = atomidx + 1

                if worstdist > max_deviation:
                    # Remove the furthest-away atom and try again
                    rwm = Chem.RWMol(self.mcs_mol)
                    rwm.RemoveAtom(worstatomidx)
                    if verbose == "pedantic":
                        logging.info(
                            f"Removing atom {worstatomidx} from MCS based on distance {worstdist}"
                        )
                    self.mcs_mol = Chem.Mol(rwm)
                else:
                    break

        def trim_mcs_fix_broken_rdkit_code() -> None:
            """
            Detect cases where the RDKit has generated an incorrect MCS (for alpha vs beta naphthyl, for
            instance), and delete some atoms to break the ring. The excess atoms will then be
            removed later by delete_broken_ring()

            Can be removed once RDKit is fixed

            Algorithm: find a bond in moli where the atoms are both in the MCS, they are bonded in
            moli, but are not bonded in the MCS.
            """

            to_remove = []
            for ai in self.moli.GetAtoms():
                if ai.HasProp("to_mcs"):  # is ai in the MCS?
                    aimcs = int(ai.GetProp("to_mcs"))
                    for bai in ai.GetNeighbors():
                        if bai.HasProp("to_mcs"):  # Atom bonded to ai is also in the MCS
                            baimcs = int(bai.GetProp("to_mcs"))
                            if aimcs < baimcs:  # only do each bond once!
                                # Check if the corresponding MCS atoms are bonded
                                if not self.mcs_mol.GetBondBetweenAtoms(aimcs, baimcs):
                                    to_remove.append(aimcs)
                                    if verbose == "pedantic":
                                        logging.info(
                                            f"Bond in first mol between atoms {ai.GetIdx()} and {bai.GetIdx()} not matched in MCS"
                                        )

            if to_remove:
                # Delete atoms from the MCS, highest index first
                to_remove.sort(reverse=True)

                if verbose == "pedantic":
                    logging.info(
                        f"Removing {len(to_remove)} atoms from MCS based on detection of broken RDKit ring bond matching"
                    )

                edit_mcs_mol = Chem.EditableMol(self.mcs_mol)
                for i in to_remove:
                    edit_mcs_mol.RemoveAtom(i)

                self.mcs_mol = edit_mcs_mol.GetMol()
                map_mcs_mol()  # Regenerate mappings

        def trim_mcs_chiral_atoms() -> None:
            """
            Remove all atoms in the MCS where there might be a chirality inversion i.e.
            (a) the corresponding atoms in the input molecules are both chiral, and
            (b) the parity of the atom mapping in the input molecules is reversed

            Calls map_mcs_mol as it uses the mappings generated there.

            """

            def permutation_parity(perm: list[int]) -> bool:
                """
                Returns the parity of the provided permutation tuple/array: even parity
                is True, odd parity is False.

                Parameters
                ----------
                perm : list[int]
                  The list of permutations.

                Returns
                -------
                bool
                  ``True`` if parity is even, otherwise ``False``.
                """
                parity = True
                for i in range(len(perm) - 1):
                    for j in range(i + 1, len(perm)):
                        if perm[i] < perm[j]:
                            parity = not parity

                return parity

            def atom_mcs_chiral_parity(a: Chem.Atom) -> Chem.rdchem.ChiralType:
                """
                Take the neighbours of chiral atom a. Get the index of each of these atoms
                in the MCS. Combine the parity of this list with the chirality flag for
                a to determine the "MCS parity".

                Parameters
                ----------
                a : Chem.Atom
                  The chiral atom to inspect.

                Returns
                -------
                Chem.rdchem.ChiralType
                  The MCS chiral parity of the atom.
                """
                nbrs = []
                for aj in a.GetNeighbors():
                    try:
                        nbrs.append(int(aj.GetProp("to_mcs")))
                    except Exception:
                        nbrs.append(1000)  # should not be more than one!

                if not permutation_parity(nbrs):
                    if a.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW:
                        return Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW
                    if a.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW:
                        return Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW
                return a.GetChiralTag()

            def flag_inverted_atoms_in_mcs() -> None:
                """
                Flag all atoms in the MCS where the chirality is inverted between
                moli and molj. Atoms are flagged with CHI_TETRAHEDRAL_CW
                """

                # Generate atommappings as they are useful below
                map_mcs_mol()

                chiral_at_moli = [seq[0] for seq in Chem.FindMolChiralCenters(moli)]
                chiral_at_molj = [seq[0] for seq in Chem.FindMolChiralCenters(molj)]

                invertedatoms = []

                for i in chiral_at_moli:
                    # Is atom i in the MCS?
                    ai = moli.GetAtomWithIdx(i)
                    if ai.HasProp("to_mcs"):
                        for j in chiral_at_molj:
                            # Is atom j in the MCS?
                            aj = molj.GetAtomWithIdx(j)
                            if aj.HasProp("to_mcs"):
                                # Are they the same atom?
                                if ai.GetProp("to_mcs") == aj.GetProp("to_mcs"):
                                    # OK, atoms are both chiral, and match the same MCS atom.
                                    # Take the list of neighbours for ai, and get their indices in
                                    # the MCS. Use the parity of this index list together with the
                                    # chiral parity of ai to work out the "MCS parity". Do the same
                                    # for aj and check if the two are the same.
                                    #
                                    # If not, flag with the CHI_TETRAHEDRAL_CW property.
                                    pi = atom_mcs_chiral_parity(ai)
                                    pj = atom_mcs_chiral_parity(aj)
                                    if pi != pj:
                                        invertedatoms.append(int(aj.GetProp("to_mcs")))

                for i in invertedatoms:
                    mcsat = self.mcs_mol.GetAtomWithIdx(i)
                    mcsat.SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW)
                    if verbose == "pedantic":
                        logging.info(f"Inverted chiral atom detected: {i}")

            # Flag inverted atoms
            flag_inverted_atoms_in_mcs()

            # Trim inverted chiral Atoms. The algorithm is to delete the chiral centre,
            # fragment the molecule, and keep only the two largest fragments. Rinse and
            # repeat until no more flagged chiral centres remain.

            while True:
                atom_idx = -1

                for atom in self.mcs_mol.GetAtoms():
                    # Note that any atom in the MCS which has inverted chirality between the input mols is
                    # flagged with CHI_TETRAHEDRAL_CW
                    if atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW:
                        atom_idx = atom.GetIdx()
                        atom.SetChiralTag(
                            Chem.rdchem.ChiralType.CHI_UNSPECIFIED
                        )  # Remove from consideration for the next loop
                        break

                if atom_idx == -1:  # Not found any more chiral atoms, so done
                    break

                # Move the chiral atom to the end (avoids indexing problems)
                newindexes = list(range(self.mcs_mol.GetNumAtoms()))
                newindexes.remove(atom_idx)
                newindexes.append(atom_idx)
                self.mcs_mol = Chem.RenumberAtoms(self.mcs_mol, newindexes)

                # Now we loop, deleting groups attached to the chiral atom, until the
                # chiral atom has at most two heavy atom connections
                # Note that getAtoms()[-1] returns the first atom not the last if you
                # don't convert it to a list. Grr.
                while list(self.mcs_mol.GetAtoms())[-1].GetDegree() > 2:
                    # Delete the chiral atom in a temporary molecule, and fragment. Since the
                    # chiral atom was the last one, the indexes in the temporary molecule are the
                    # same as in self.mcs_mol
                    edit_mol = Chem.EditableMol(self.mcs_mol)
                    edit_mol.RemoveAtom(self.mcs_mol.GetNumAtoms() - 1)
                    tmp_mol = edit_mol.GetMol()
                    fragments = Chem.rdmolops.GetMolFrags(tmp_mol)

                    # Get index of smallest fragments
                    min_idx = 0
                    lgt_min = 10000

                    for idx in range(0, len(fragments)):
                        lgt = len(fragments[idx])
                        if lgt < lgt_min:
                            lgt_min = lgt
                            min_idx = idx

                    # Get the atoms in this fragment and sort them so we delete the
                    # largest index first
                    min_frag = list(fragments[min_idx])
                    min_frag.sort(reverse=True)

                    if verbose == "pedantic":
                        logging.info(f"Removing {len(min_frag)} atoms to remove chiral inversion")
                    edit_mol = Chem.EditableMol(self.mcs_mol)
                    for idx in min_frag:
                        edit_mol.RemoveAtom(idx)
                    self.mcs_mol = edit_mol.GetMol()

            map_mcs_mol()  # Regenerate mappings after deletion
            # Done!

        def delete_broken_ring() -> None:
            """
            This function checks the MCS to see if there are any
            atoms which are in a ring in the parent molecules, but
            not in a ring in the MCS. This may occur if we have deleted
            some atoms from the MCS in 3D coordinate matching, for
            example.
            """

            to_remove = []
            for at in self.mcs_mol.GetAtoms():
                moli_idx = int(at.GetProp("to_moli"))
                moli_at = self._moli_noh.GetAtomWithIdx(moli_idx)
                molj_idx = int(at.GetProp("to_molj"))
                molj_at = self._molj_noh.GetAtomWithIdx(molj_idx)

                # Testing moli and molj is redundant due to the way that the
                # MCS is calculated, but I'd rather be paranoid here
                if moli_at.IsInRing() and molj_at.IsInRing() and not at.IsInRing():
                    to_remove.append(at.GetIdx())

            if to_remove:
                # Delete atoms from the MCS, highest index first
                to_remove.sort(reverse=True)

                if verbose == "pedantic":
                    logging.info(
                        f"Removing {len(to_remove)} atoms from MCS to clear up partial rings"
                    )

                edit_mcs_mol = Chem.EditableMol(self.mcs_mol)
                for i in to_remove:
                    edit_mcs_mol.RemoveAtom(i)

                self.mcs_mol = edit_mcs_mol.GetMol()

                map_mcs_mol()  # Regenerate mappings after deletion

        def map_mcs_mol() -> None:
            """
            This function is used to define a map between the generated mcs, the
            molecules and vice versa
            """

            # Get self-mapping for the MCS
            mcsi_sub = tuple(range(self.mcs_mol.GetNumAtoms()))

            moli_sub, molj_sub = best_substruct_match_to_mcs(
                self._moli_noh, self._molj_noh, by_rmsd=threed, use_shift=shift
            )

            # mcs to moli
            map_mcs_mol_to_moli_sub = list(zip(mcsi_sub, moli_sub))

            # Clear all properties as we may call this function more than once
            # Here are all the properties we use:
            # `to_moli` and `to_molj`: this gives the atom correspondence
            # between the mcs molecule and the heavy atom only moli and molj
            # molecules.
            # `to_moli_all` and `to_molj_all`: this gives the atom
            # correspondence between the MCS molecules and the input (pre
            # removal of hydrogens) molecules.
            # `to_mcs` this gives the correspondence between the full
            # (including hydrogens) molecule and the mcs molecule
            for a in self.mcs_mol.GetAtoms():
                a.ClearProp("to_moli")
                a.ClearProp("to_moli_all")
                a.ClearProp("to_molj")
                a.ClearProp("to_molj_all")
            for a in self.moli.GetAtoms():
                a.ClearProp("to_mcs")
            for a in self.molj.GetAtoms():
                a.ClearProp("to_mcs")

            # An RDkit atomic property is defined to store the mapping to moli
            for idx in map_mcs_mol_to_moli_sub:
                self.mcs_mol.GetAtomWithIdx(idx[0]).SetProp("to_moli", str(idx[1]))
                self.mcs_mol.GetAtomWithIdx(idx[0]).SetProp(
                    "to_moli_all", str(self._map_moli_noh[idx[1]])
                )
                self.moli.GetAtomWithIdx(self._map_moli_noh[idx[1]]).SetProp("to_mcs", str(idx[0]))

            mcsj_sub = tuple(range(self.mcs_mol.GetNumAtoms()))

            # mcs to molj
            map_mcs_mol_to_molj_sub = list(zip(mcsj_sub, molj_sub))

            # Map between the two molecules
            self._map_moli_molj = list(zip(moli_sub, molj_sub))

            # An RDkit atomic property is defined to store the mapping to molj
            for idx in map_mcs_mol_to_molj_sub:
                self.mcs_mol.GetAtomWithIdx(idx[0]).SetProp("to_molj", str(idx[1]))
                self.mcs_mol.GetAtomWithIdx(idx[0]).SetProp(
                    "to_molj_all", str(self._map_molj_noh[idx[1]])
                )
                self.molj.GetAtomWithIdx(self._map_molj_noh[idx[1]]).SetProp("to_mcs", str(idx[0]))

            # For each mcs atom we save its original index in a specified
            # property. This could be very useful in the code development
            # when deletion or atom insertions are performed
            for at in self.mcs_mol.GetAtoms():
                at.SetProp("org_idx", str(at.GetIdx()))

            return

        def set_ring_counter(mol: Chem.Mol) -> None:
            """
            This function is used to attach to each molecule atom a ring counter
            rc. This parameter is used to assess if a ring has been broken or not
            during the MCS mapping

            Parameters
            ----------
            mol : Chem.Mol
              The RDKit Molecule used to define the atom ring counters.
            """

            # set to zero the atom ring counters
            for at in mol.GetAtoms():
                at.SetProp("rc", "0")

            rginfo = mol.GetRingInfo()

            rings = rginfo.AtomRings()

            rings_set = {idx for ring in rings for idx in ring}

            for idx in rings_set:
                for ring in rings:
                    if idx in ring:
                        val = int(mol.GetAtomWithIdx(idx).GetProp("rc"))
                        val = val + 1
                        mol.GetAtomWithIdx(idx).SetProp("rc", str(val))
            return

        # START of __init__ function
        # Set logging level and format
        logging.basicConfig(format="%(levelname)s:\t%(message)s", level=logging.INFO)

        # Global beta setting for atom penalties
        self.beta = 0.1

        # Local pointers to the passed molecules
        self.moli = moli
        self.molj = molj

        # Sanitize input molecules
        Chem.SanitizeMol(self.moli)
        Chem.SanitizeMol(self.molj)

        # Set chirality flags from 3D coords if working in 3D
        if threed:
            Chem.rdmolops.AssignAtomChiralTagsFromStructure(self.moli, replaceExistingTags=True)
            Chem.rdmolops.AssignAtomChiralTagsFromStructure(self.molj, replaceExistingTags=True)

        if not verbose == "pedantic":
            lg = RDLogger.logger()
            lg.setLevel(RDLogger.CRITICAL)

        # Local pointers to the passed molecules without hydrogens
        # These variables are defined as private
        try:
            self._moli_noh = Chem.RemoveHs(moli)
            self._molj_noh = Chem.RemoveHs(molj)
        except Exception:
            self._moli_noh = Chem.RemoveHs(moli, sanitize=False)
            self._molj_noh = Chem.RemoveHs(molj, sanitize=False)

            Chem.SanitizeMol(self._moli_noh, sanitizeOps=Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
            Chem.SanitizeMol(self._molj_noh, sanitizeOps=Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)

        # Get maps of the atom correspondences between the no hydrogen
        # molecules and the original molecules
        self._map_moli_noh = self._heavy_to_all_pos_remap(self._moli_noh, moli)
        self._map_molj_noh = self._heavy_to_all_pos_remap(self._molj_noh, molj)

        # MCS calculation. In RDKit the MCS is a smart string. Ring atoms are
        # always mapped in ring atoms.
        # Don't add the mcs result as a member variable as it can't be pickled
        if element_change:
            atom_compare = rdFMCS.AtomCompare.CompareAny
        else:
            atom_compare = rdFMCS.AtomCompare.CompareElements

        __mcs = rdFMCS.FindMCS(
            [self._moli_noh, self._molj_noh],
            timeout=time,
            atomCompare=atom_compare,
            bondCompare=rdFMCS.BondCompare.CompareAny,
            matchValences=False,
            ringMatchesRingOnly=True,
            completeRingsOnly=True,
            matchChiralTag=False,
            seedSmarts=seed,
        )

        # Note that we need matchChiralTag=False as we want to match chiral atoms with different
        # parities, we just need to trim the MCS to the largest possible match that doesn't have
        # a mismatched chiral centre in it (eg we want to match FC[C@](C)CO to FC[C@@](C)CO
        # using the MCS FCCCO - this includes the chiral atom but we delete the methyl group
        # that led it to be chiral. The trimming is done in trim_mcs_chiral_atoms()

        # Checking
        if __mcs.canceled:
            logging.warning("Timeout reached to find the MCS between the molecules")

        if __mcs.numAtoms == 0:
            raise ValueError("No MCS was found between the molecules")

        # The found MCS pattern (smart strings) is converted to a RDKit molecule
        self.mcs_mol_smarts = __mcs.smartsString
        mcs_mol = Chem.MolFromSmarts(__mcs.smartsString)
        if mcs_mol is None:
            raise ValueError("Failed to convert MCS SMARTS to molecule")
        self.mcs_mol = mcs_mol

        # There's a symmetry-related bug here: if there was more than one MCS
        # of the same size and score, we'll get only one at random. We then try
        # to choose the mapping that matches 3D coords the best, but one of the
        # not-considered MCSes that we never saw may give a better mapping.
        # We can rescue some of this by converting all partial-query atoms to
        # full query atoms
        testmol = Chem.MolFromSmarts("*")  # Create a "match anything" query atom for us to copy
        for a in self.mcs_mol.GetAtoms():
            if a.DescribeQuery().startswith("AtomOr"):  # Matches more than one element
                a.SetQuery(
                    list(testmol.GetAtoms())[0]
                )  # Set this atom to a copy of the "match anything" atom

        try:  # Try to sanitize the MCS molecule
            Chem.SanitizeMol(self.mcs_mol)
        except Exception:  # if not, try to recover the atom aromaticity which is
            # important for the ring counter
            sanitFail = Chem.SanitizeMol(
                self.mcs_mol,
                sanitizeOps=Chem.SanitizeFlags.SANITIZE_SETAROMATICITY,
                catchErrors=True,
            )
            if sanitFail:  # if not, the MCS is skipped
                raise ValueError("Sanitization Failed...")

        # Trim the MCS to remove atoms with too-large real-space deviations
        if max3d > 0:
            try:
                trim_mcs_mol(max_deviation=max3d, use_shift=shift)
            except Exception as e:
                raise ValueError(str(e))

        # Trim the MCS further to remove chirality mismatches
        trim_mcs_chiral_atoms()

        # Check to see if we've hit the RDKit incorrect-MCS bug
        trim_mcs_fix_broken_rdkit_code()

        # Cleanup any partial rings remaining
        delete_broken_ring()

        # Only single fragment
        mols = Chem.GetMolFrags(self.mcs_mol, asMols=True, sanitizeFrags=False)
        self.mcs_mol = max(mols, key=lambda x: x.GetNumAtoms())

        # Mapping between the found MCS molecule and moli,  molj
        try:
            map_mcs_mol()
        except Exception as e:
            raise ValueError(str(e))

        # Set the ring counters for each molecule
        set_ring_counter(self._moli_noh)
        set_ring_counter(self._molj_noh)
        set_ring_counter(self.mcs_mol)

        # for at in self.mcs_mol.GetAtoms():
        #     print 'at = %d rc = %d' % (at.GetIdx(), int(at.GetProp('rc')))

        if not verbose == "pedantic":
            lg.setLevel(RDLogger.WARNING)

        return

    @staticmethod
    def _heavy_to_all_pos_remap(
        heavy_mol: Chem.Mol, all_mol: Chem.Mol, tolerance: float = 0.5
    ) -> dict[int, int]:
        """
        Convenience method to map a molecule without hydrogens
        (``heavy_mol``) back to its pre RemoveHs version (``all_mol``).

        Parameters
        ----------
        heavy_mol : Chem.Mol
          Heavy atom RDKit Molecule version of ``all_mol``.
        all_mol : Chem.Mol
          RDKit Molecule with all atoms present (i.e. pre-RemoveHs).
        tolerance : float, default 0.5
          Maximum deviation acceptable between corresponding atoms in Angstroms.

        Returns
        -------
        mapping : dict[int, int]
          Dictionary mapping the atom indices of the heavy molecule
          to the all-atom molecule.
        """
        mapping: dict[int, int] = {}
        for at in heavy_mol.GetAtoms():
            best = tolerance
            at_idx = at.GetIdx()
            pos = heavy_mol.GetConformer().GetAtomPosition(at_idx)
            for at2 in all_mol.GetAtoms():
                at2_idx = at2.GetIdx()
                pos2 = all_mol.GetConformer().GetAtomPosition(at2_idx)
                if (pos - pos2).Length() < best:
                    best = (pos - pos2).Length()
                    mapping[at_idx] = at2_idx
        return mapping

    @staticmethod
    def getMapping(
        moli: Chem.Mol,
        molj: Chem.Mol,
        hydrogens: bool = False,
        fname: str | None = None,
        time_out: int = 150,
    ) -> Iterator[tuple[int, int]]:
        """Compute the MCS between two passed molecules

        Parameters
        ----------
        moli : Chem.Mol
          The first molecule used to perform the MCS calculation.
        molj : Chem.Mol
          The second molecule used to perform the MCS calculation.
        hydrogens : bool, default False
          Include or not the hydrogens in the MCS calculation.
        fname : string | None, default None
          The filename used to output a png file depicting the MCS mapping.
        time_out: int, default 150
          The maximum time in seconds which can be used to compute the MCS.

        Returns
        -------
        map_moli_molj: Iterator[tuple[int, int]]
          Iterator of tuples which contains the atom mapping indexes between
          the two molecules. The indexes (i,j) are respectively related to
          the first (moli) and the second (molj) passed molecules.
        """

        # Molecule copies
        moli_c = Chem.Mol(moli)
        molj_c = Chem.Mol(molj)

        if not hydrogens:
            moli_c = Chem.RemoveHs(moli_c)
            molj_c = Chem.RemoveHs(molj_c)

        # MCS calculation. In RDKit the MCS is a smart string. Ring atoms are
        # always mapped in ring atoms.
        mcs = rdFMCS.FindMCS(
            [moli_c, molj_c],
            timeout=time_out,
            atomCompare=rdFMCS.AtomCompare.CompareAny,
            bondCompare=rdFMCS.BondCompare.CompareAny,
            matchValences=False,
            ringMatchesRingOnly=True,
            completeRingsOnly=False,
            matchChiralTag=False,
        )

        # Checking
        if mcs.canceled:
            raise ValueError("Timeout! No MCS found between passed molecules")

        if mcs.numAtoms == 0:
            raise ValueError("No MCS was found between the molecules")

        # The found MCS pattern (smart strings) is converted to a RDKit molecule
        mcs_mol = Chem.MolFromSmarts(mcs.smartsString)

        try:
            Chem.SanitizeMol(mcs_mol)
        except Exception:  # if not try to recover the atom aromaticity which is
            # important for the ring counter
            sanitFail = Chem.SanitizeMol(
                mcs_mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_SETAROMATICITY, catchErrors=True
            )
            if sanitFail:  # if not the MCS is skipped
                raise ValueError("Sanitization Failed...")

        # mcs indexes mapped back to the first molecule moli
        if moli_c.HasSubstructMatch(mcs_mol):
            moli_sub = moli_c.GetSubstructMatch(mcs_mol)
        else:
            raise ValueError("RDkit MCS Subgraph first molecule search failed")
        # mcs indexes mapped back to the second molecule molj
        if molj_c.HasSubstructMatch(mcs_mol):
            molj_sub = molj_c.GetSubstructMatch(mcs_mol)
        else:
            raise ValueError("RDkit MCS Subgraph second molecule search failed")

        if mcs_mol.HasSubstructMatch(mcs_mol):
            mcs_sub = mcs_mol.GetSubstructMatch(mcs_mol)
        else:
            raise ValueError("RDkit MCS Subgraph search failed")

        # Map between the two molecules
        map_moli_to_molj = zip(moli_sub, molj_sub)

        # depict the mapping by using a .png file
        if fname:
            Chem.rdDepictor.Compute2DCoords(moli_c)
            Chem.rdDepictor.Compute2DCoords(molj_c)
            Chem.rdDepictor.Compute2DCoords(mcs_mol)

            DrawingOptions.includeAtomNumbers = True

            moli_fname = "Moli"
            molj_fname = "Molj"
            mcs_fname = "Mcs"

            img = Draw.MolsToGridImage(
                [moli_c, molj_c, mcs_mol],
                molsPerRow=3,
                subImgSize=(400, 400),
                legends=[moli_fname, molj_fname, mcs_fname],
                highlightAtomLists=[moli_sub, molj_sub, mcs_sub],
            )

            img.save(fname)

            DrawingOptions.includeAtomNumbers = False

        return map_moli_to_molj

    def mcsr(self) -> float:
        """
        This rule computes the similarity between the two passed molecules
        used to compute the MCS.

        Returns
        -------
        scr_mcsr : float
          The MCSR rule score.
        """

        # The number of heavy atoms in each molecule
        nha_moli = self.moli.GetNumHeavyAtoms()
        nha_molj = self.molj.GetNumHeavyAtoms()
        # Note that the mcs_mol (a) doesn't contain hydrogens, and (b) does contain
        # wildcard atoms, which don't count as 'heavy'. Use the total atom count instead.
        nha_mcs_mol = self.mcs_mol.GetNumAtoms()

        # score
        scr_mcsr = math.exp(-self.beta * (nha_moli + nha_molj - 2 * nha_mcs_mol))

        logging.info(
            f"MCSR from MCS size {nha_mcs_mol}, molecule sizes {nha_moli}, {nha_molj} is {scr_mcsr}"
        )

        return scr_mcsr

    def mncar(self, ths: int = 4) -> float:
        """
        This rule cuts the similarity score between two molecules if they do
        not share the selected number of atoms.

        Parameters
        ----------
        ths : int
          The minimum number of atoms to share.

        Returns
        -------
        scr_mncar : float
          The MNCAR rule score.
        """

        # This rule has been modified from the rule described in the Lomap paper
        # to match the LOMAP first implementation provided by schrodinger

        nha_mcs_mol = self.mcs_mol.GetNumHeavyAtoms()
        nha_moli = self.moli.GetNumHeavyAtoms()
        nha_molj = self.molj.GetNumHeavyAtoms()

        scr_mncar = float((nha_mcs_mol >= ths) or (nha_moli < ths + 3) or (nha_molj < ths + 3))

        return scr_mncar

    def tmcsr(self, strict_flag: bool = True) -> float:
        """
        TMCSR (Trim) rule.
        This score is no longer implemented and now always returns 1.0.

        Notes
        -----
        * MDM - we don't use this as we don't have the same limitation on
          partial ring deletion as Schrodinger.
        * Removed the chirality check, the MCS is now trimmed to remove
          chirality.
        """
        warnings.warn(
            "The TMCSR is deprecated and will be removed in a future release",
            DeprecationWarning,
            stacklevel=2,
        )
        return 1.0

    def atomic_number_rule(self) -> float:
        """
        This rule checks how many elements have been changed in the MCS
        and returns a score based on the fraction of MCS matches that are the same atomic number.
        When used with beta=0.1 and multiplied by mcsr, this is equivalent to counting
        mismatched atoms at only half weight.

        This has been extended to modify the amount of mismatch according to the
        atoms being mapped.

        Returns
        -------
        an_score : float
          Atomic number rule score.
        """

        # A value of 0.5 is the same behaviour as before, a value of 1 means that the
        # atoms are perfectly equivalent, a value of 0 means that the atoms are perfectly
        # non-equivalent (i.e the penalty should basically remove this atom pair from the
        # MCS). The default for pairs not in this data structure is 0.5.
        #
        # Note that we don't need the symmetry equivalent values: we will use the large of
        # [i][j] and [j][i]
        transform_difficulty = {
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
        nmismatch: float = 0
        for at in self.mcs_mol.GetAtoms():
            moli_idx = int(at.GetProp("to_moli"))
            molj_idx = int(at.GetProp("to_molj"))
            moli_a = self._moli_noh.GetAtomWithIdx(moli_idx)
            molj_a = self._molj_noh.GetAtomWithIdx(molj_idx)

            if moli_a.GetAtomicNum() != molj_a.GetAtomicNum():
                ij: float = -1
                ji: float = -1
                try:
                    ij = transform_difficulty[moli_a.GetAtomicNum()][molj_a.GetAtomicNum()]
                except KeyError:
                    pass
                try:
                    ji = transform_difficulty[molj_a.GetAtomicNum()][moli_a.GetAtomicNum()]
                except KeyError:
                    pass
                diff = max(ij, ji)
                if diff == -1:
                    diff = 0.5  # default for elements not found

                nmismatch += 1 - diff

        an_score = math.exp(-1 * self.beta * nmismatch)
        logging.info(f"atomic number score from {nmismatch} mismatches is {an_score}")
        return an_score

    def hybridization_rule(self, penalty_weight: float = 1.5) -> float:
        """
        This rule checks how many atoms have changed hybridization state.

        Parameters
        ----------
        penalty_weight : float, default 1.5
          How many "atoms" different a hybridization state change has.
          1 means that the atom is effectively removed from the MCS for scoring
          purposes, 0 means that the hybridization changes are free.
          When used with ``beta`` of ``0.1``, and multiplied by ``mcsr``,
          this is equivalent to counting mismatched atoms at a weight of
          ``(1 - penalty_weight)``.

        Returns
        -------
        hyb_score : float
          The hybridization rule score.
        """
        nmismatch: int = 0
        for at in self.mcs_mol.GetAtoms():
            moli_idx = int(at.GetProp("to_moli"))
            molj_idx = int(at.GetProp("to_molj"))
            moli_a = self._moli_noh.GetAtomWithIdx(moli_idx)
            molj_a = self._molj_noh.GetAtomWithIdx(molj_idx)

            hybi = atom_hybridization(moli_a)
            hybj = atom_hybridization(molj_a)
            mismatch = hybi != hybj

            # Allow Nsp3 to match Nsp2, otherwise guanidines etc become painful
            if (
                moli_a.GetAtomicNum() == 7
                and molj_a.GetAtomicNum() == 7
                and (hybi in [2, 3])
                and hybj in [2, 3]
            ):
                mismatch = False

            if mismatch:
                nmismatch += 1
                logging.info(
                    f"Hybridization mismatch {moli_a.GetIdx()} {moli_a.GetSymbol()} {hybi} vs {molj_a.GetIdx()} {molj_a.GetSymbol()} {hybj}"
                )

        hyb_score = math.exp(-1 * self.beta * nmismatch * penalty_weight)
        logging.info(f"hybridization score from {nmismatch} mismatches is {hyb_score}")
        return hyb_score

    def sulfonamides_rule(self, penalty: int = 4) -> float:
        """
        This rule checks to see if we are growing a complete sulfonamide, and
        returns 0 if we are. This means that if this rule is used we effectively disallow
        this transition. Testing has shown that growing -SO2NH2 from scratch performs
        very badly.

        Parameters
        ----------
        penalty : int, default 4
          The number of atom mismatches that failing this rule will lower the score by.

        Returns
        -------
        sulf_score : float
          The sulfonamides rule score.
        """

        def adds_sulfonamide(mol: Chem.Mol) -> bool:
            """
            Check if the transformation would add a new sulfonamide group.

            Parameters
            ----------
            mol : Chem.Mol
              The Molecule to inspect.

            Returns
            -------
            bool
              ``True`` if removing the MCS from the Molecule leaves a sulfonamide.
            """

            if not mol.HasSubstructMatch(self.mcs_mol):
                raise ValueError("RDkit MCS Subgraph molecule search failed in sulfonamide check")

            rwm = rdmolops.DeleteSubstructs(mol, self.mcs_mol)
            return rwm.HasSubstructMatch(Chem.MolFromSmarts("S(=O)(=O)N"))

        fail = 1 if (adds_sulfonamide(self._moli_noh)) else 0
        fail = 1 if (adds_sulfonamide(self._molj_noh)) else fail
        sulf_score = math.exp(-1 * self.beta * fail * penalty)
        logging.info(f"sulfonamide score is {sulf_score}")
        return sulf_score

    def heterocycles_rule(self, penalty: int = 4) -> float:
        """
        This rule checks to see if we are growing a heterocycle from a hydrogen, and
        returns <1 if we are. This means that if this rule is used we penalise
        this transition. Testing has shown that growing a pyridine or other heterocycle
        is unlikely to work (better to grow phenyl then mutate)

        Parameters
        ----------
        penalty : int, default 4
          The number of atom mismatches that failing this rule will lower the score by.

        Returns
        -------
        het_score : float
          The heterocycles rule score.
        """

        def adds_heterocycle(mol: Chem.Mol) -> bool:
            """
            Returns true if the removal of the MCS from the provided molecule
            leaves a heterocycle.

            Parameters
            ----------
            mol : Chem.Mol
              Molecule to inspect.

            Returns
            -------
            bool
              ``True`` if a heterocycle is left after removing the MCS.
            """

            if not mol.HasSubstructMatch(self.mcs_mol):
                raise ValueError("RDkit MCS Subgraph molecule search failed in heterocycle check")

            rwm = rdmolops.DeleteSubstructs(mol, self.mcs_mol)
            # Only picking up N/C containing heterocycles - odd cases like pyran derivatives are not caught
            grow6mheterocycle = rwm.HasSubstructMatch(
                Chem.MolFromSmarts("[n]1[c,n][c,n][c,n][c,n][c,n]1")
            )

            # Note that growing pyrrole, furan or thiophene is allowed
            grow5mheterocycle = rwm.HasSubstructMatch(
                Chem.MolFromSmarts("[o,n,s]1[n][c,n][c,n][c,n]1")
            )
            grow5mheterocycle |= rwm.HasSubstructMatch(
                Chem.MolFromSmarts("[o,n,s]1[c,n][n][c,n][c,n]1")
            )
            return grow6mheterocycle | grow5mheterocycle

        fail = 1 if (adds_heterocycle(self._moli_noh)) else 0
        fail = 1 if (adds_heterocycle(self._molj_noh)) else fail
        het_score = math.exp(-1 * self.beta * fail * penalty)
        logging.info(f"heterocycle score is {het_score}")
        return het_score

    def transmuting_methyl_into_ring_rule(self, penalty: int = 6) -> float:
        """
         Rule to prevent turning a methyl into a ring atom and similar transformations
         (you can grow a ring, but you can't transmute into one)

        Parameters
        ----------
        penalty : int, default 6
          The number of atom mismatches that failing this rule will lower the score by.

        Returns
        -------
        mescore : float
          The methyl to ring transmuting score.
        """
        moli = self._moli_noh
        molj = self._molj_noh

        # Get list of bonds in mol i and j that go from the MCS to a non-MCS atom,
        # arranged in tuples with the index of the MCS atom
        moli_sub = moli.GetSubstructMatch(self.mcs_mol)
        molj_sub = molj.GetSubstructMatch(self.mcs_mol)

        is_bad = False

        for i in range(0, len(moli_sub)):
            edge_bondsi = [
                b.GetBeginAtomIdx()
                for b in moli.GetBonds()
                if (b.GetEndAtomIdx() == moli_sub[i] and b.GetBeginAtomIdx() not in moli_sub)
            ]
            edge_bondsi += [
                b.GetEndAtomIdx()
                for b in moli.GetBonds()
                if (b.GetBeginAtomIdx() == moli_sub[i] and b.GetEndAtomIdx() not in moli_sub)
            ]
            edge_bondsj = [
                b.GetBeginAtomIdx()
                for b in molj.GetBonds()
                if (b.GetEndAtomIdx() == molj_sub[i] and b.GetBeginAtomIdx() not in molj_sub)
            ]
            edge_bondsj += [
                b.GetEndAtomIdx()
                for b in molj.GetBonds()
                if (b.GetBeginAtomIdx() == molj_sub[i] and b.GetEndAtomIdx() not in molj_sub)
            ]

            for edgeAtom_i in edge_bondsi:
                for edgeAtom_j in edge_bondsj:
                    if (
                        moli.GetAtomWithIdx(edgeAtom_i).IsInRing()
                        ^ molj.GetAtomWithIdx(edgeAtom_j).IsInRing()
                    ):
                        is_bad = True

        mescore = math.exp(-1 * self.beta * penalty) if is_bad else 1
        logging.info(f"methyl-to-ring transformation score is {mescore}")
        return mescore

    def transmuting_ring_sizes_rule(self) -> float:
        """
        Rule to prevent turning a ring atom into a ring atom with a different ring size
        (you can grow a ring, but you can't turn a cyclopentyl into a cyclohexyl)

        Hard rule: sets score to near zero if violated

        Returns
        -------
        float
          The ring size transmuting rule score.
        """
        moli = self._moli_noh
        molj = self._molj_noh

        # Get list of bonds in mol i and j that go from the MCS to a non-MCS atom,
        # arranged in tuples with the index of the MCS atom
        moli_sub = moli.GetSubstructMatch(self.mcs_mol)
        molj_sub = molj.GetSubstructMatch(self.mcs_mol)

        is_bad = False

        for i in range(0, len(moli_sub)):
            edge_bondsi = [
                b.GetBeginAtomIdx()
                for b in moli.GetBonds()
                if (b.GetEndAtomIdx() == moli_sub[i] and b.GetBeginAtomIdx() not in moli_sub)
            ]
            edge_bondsi += [
                b.GetEndAtomIdx()
                for b in moli.GetBonds()
                if (b.GetBeginAtomIdx() == moli_sub[i] and b.GetEndAtomIdx() not in moli_sub)
            ]
            edge_bondsj = [
                b.GetBeginAtomIdx()
                for b in molj.GetBonds()
                if (b.GetEndAtomIdx() == molj_sub[i] and b.GetBeginAtomIdx() not in molj_sub)
            ]
            edge_bondsj += [
                b.GetEndAtomIdx()
                for b in molj.GetBonds()
                if (b.GetBeginAtomIdx() == molj_sub[i] and b.GetEndAtomIdx() not in molj_sub)
            ]
            # print("Atom",i,"index",moli_sub[i],"edge atoms on mol 1 are",edge_bondsi);
            # print("Atom",i,"index",molj_sub[i],"edge atoms on mol 2 are",edge_bondsj);

            for edgeAtom_i in edge_bondsi:
                for edgeAtom_j in edge_bondsj:
                    # print("Checking ring for atom",edgeAtom_i,edgeAtom_j,moli.GetAtomWithIdx(edgeAtom_i).IsInRing(),molj.GetAtomWithIdx(edgeAtom_j).IsInRing())
                    if (
                        moli.GetAtomWithIdx(edgeAtom_i).IsInRing()
                        and molj.GetAtomWithIdx(edgeAtom_j).IsInRing()
                    ):
                        for ring_size in range(3, 8):
                            if moli.GetAtomWithIdx(edgeAtom_i).IsInRingSize(
                                ring_size
                            ) ^ molj.GetAtomWithIdx(edgeAtom_j).IsInRingSize(ring_size):
                                logging.info(
                                    f"transforming ring sizes score is 0 based on atom {edgeAtom_i} in moli and {edgeAtom_j} in molj"
                                )
                                is_bad = True
                            if moli.GetAtomWithIdx(edgeAtom_i).IsInRingSize(
                                ring_size
                            ) or molj.GetAtomWithIdx(edgeAtom_j).IsInRingSize(ring_size):
                                break

        return 0.1 if is_bad else 1

    def heavy_atom_mcs_map(self) -> list[tuple[int, int]]:
        """
        Gives a list of tuples mapping atoms from moli to molj
        Heavy atoms only, returned sorted by first index.

        Returns
        -------
        maplist : list[tuple[int, int]]
          List of paired heavy atom indices.
        """
        maplist = []
        for at in self.mcs_mol.GetAtoms():
            moli_idx = int(at.GetProp("to_moli"))
            molj_idx = int(at.GetProp("to_molj"))
            maplist.append((moli_idx, molj_idx))
        maplist.sort()
        return maplist

    def heavy_atom_match_list(self) -> str:
        """
        Gives a string listing the MCS match between the two molecules as
          `atom_m1:atom_m2,atom_m1:atom_m2,...`

        Heavy atoms only.

        Returns
        -------
        str
          Comma separated string of paired heavy atom indices.
        """
        return ",".join([str(i) + ":" + str(j) for (i, j) in self.heavy_atom_mcs_map()])

    def all_atom_match_list(self) -> str:
        """
        Gives a string listing the MCS match between the two molecules as
        `atom_m1:atom_m2,atom_m1:atom_m2,...`

        All atoms including hydrogens. The string is sorted by first index.
        We need to be careful that this function is symmetric, and that hydrogens
        are mapped correctly.

        Returns
        -------
        str
          Comma separated string of paired atom indices, sorted by
          the first index.
        """

        def get_attached_atoms_not_in_mcs(mol: Chem.Mol, i: int) -> list[int]:
            """Get atoms attached to atom i which are not in the MCS"""
            attached = []
            for b in mol.GetBonds():
                if b.GetEndAtomIdx() == i or b.GetBeginAtomIdx() == i:
                    j = b.GetEndAtomIdx()
                    if j == i:
                        j = b.GetBeginAtomIdx()
                    # OK, so j is the atom at the other end of the bond from atom i. Is it in the MCS?
                    inMCS = mol.GetAtomWithIdx(j).HasProp("to_mcs")
                    if not inMCS:
                        attached.append(j)
            return attached

        moli = self.moli
        molj = self.molj

        heavy_maplist = self.heavy_atom_mcs_map()

        maplist = []

        for entry in heavy_maplist:
            new_entry = (self._map_moli_noh[entry[0]], self._map_molj_noh[entry[1]])
            maplist.append(new_entry)

        # OK, this is painful, as the MCS only includes heavies. We now need to match up
        # hydrogens hanging off the MCS

        # Iterate over all atoms in the MCS
        for at in self.mcs_mol.GetAtoms():
            moli_idx = int(at.GetProp("to_moli_all"))
            molj_idx = int(at.GetProp("to_molj_all"))
            attached_i = get_attached_atoms_not_in_mcs(moli, moli_idx)
            attached_j = get_attached_atoms_not_in_mcs(molj, molj_idx)

            # Now, we need to match these up, with the caveat that we *must* not match
            # a heavy to a heavy (as if we were allowed to match these, then they would be
            # in the MCS!
            #
            # In 3D mode, ensure that maps happen to the closest atom in 3D coordinates -
            # this gets mappings to prochiral hydrogens correct (SFT-15791)

            # Match H to H first
            while attached_i and attached_j:
                hidx_i = -1
                hidx_j = -1
                best_dist = 10000
                for ai in attached_i:
                    if moli.GetAtomWithIdx(ai).GetAtomicNum() == 1:
                        for aj in attached_j:
                            if molj.GetAtomWithIdx(aj).GetAtomicNum() == 1:
                                dist = (
                                    moli.GetConformer().GetAtomPosition(ai)
                                    - molj.GetConformer().GetAtomPosition(aj)
                                ).Length()
                                if dist < best_dist or not self.options["threed"]:
                                    hidx_i = ai
                                    hidx_j = aj
                                    best_dist = dist
                if (hidx_i < 0) and self.options["element_change"]:
                    # OK, no hydrogen-hydrogen matches left. Try to match a hydrogen to a non-hydrogen
                    for ai in attached_i:
                        for aj in attached_j:
                            if (
                                moli.GetAtomWithIdx(ai).GetAtomicNum() == 1
                                or molj.GetAtomWithIdx(aj).GetAtomicNum() == 1
                            ):
                                dist = (
                                    moli.GetConformer().GetAtomPosition(ai)
                                    - molj.GetConformer().GetAtomPosition(aj)
                                ).Length()
                                if dist < best_dist or not self.options["threed"]:
                                    hidx_i = ai
                                    hidx_j = aj
                                    best_dist = dist

                if hidx_i >= 0:
                    # Found a mappable pair: add and try again
                    maplist.append((hidx_i, hidx_j))
                    attached_i.remove(hidx_i)
                    attached_j.remove(hidx_j)
                else:
                    break  # No mappable pairs left

        maplist.sort()
        return ",".join([str(i) + ":" + str(j) for (i, j) in maplist])
