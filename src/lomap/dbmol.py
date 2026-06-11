# ******************
# MODULE DOCSTRING
# ******************

"""

LOMAP
=====

Alchemical free energy calculations hold increasing promise as an aid to drug
discovery efforts. However, applications of these techniques in discovery
projects have been relatively few, partly because of the difficulty of planning
and setting up calculations. The Lead Optimization Mapper (LOMAP) is an
automated algorithm to plan efficient relative free energy calculations between
potential ligands within a substantial set of compounds.

"""

# *****************************************************************************
# Lomap2: A toolkit to plan alchemical relative binding affinity calculations
# Copyright 2015 - 2016  UC Irvine and the Authors
#
# Authors: Dr Gaetano Calabro' and Dr David Mobley
#
# *****************************************************************************

from __future__ import annotations

import argparse
import glob
import logging
import math
import multiprocessing
import os
import pickle
import warnings
from collections.abc import Iterator, Sequence
from typing import Any, Literal

import networkx as nx
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFMCS

import lomap
from lomap import graphgen, mcs

__all__ = ["DBMolecules", "SMatrix", "Molecule"]


def formal_charge(mol: Chem.Mol) -> float:
    """
    Compute the total formal charge of a molecule.

    For mol2 files, sums the ``_TriposPartialCharge`` atom property; for
    other file types (e.g. SDF files), sums RDKit obtained ``GetFormalCharge()`` formal charges.

    Parameters
    ----------
    mol : Chem.Mol
      The molecule whose total charge is required.

    Returns
    -------
    float
      Total formal charge of the molecule.
    """
    try:
        # Assume mol2
        total_charge_mol = sum(float(a.GetProp("_TriposPartialCharge")) for a in mol.GetAtoms())
    except KeyError:
        # wasn't mol2, so assume SDF with correct formal charge props for mols
        total_charge_mol = sum(a.GetFormalCharge() for a in mol.GetAtoms())

    return total_charge_mol


def ecr(mol_i: Chem.Mol, mol_j: Chem.Mol) -> float:
    """
    Compute the similarity score between two molecules using the
    EleCtrostatic Rule (ECR).

    Parameters
    ----------
    mol_i : Chem.Mol
      The first molecule used to calculate the ECR score.
    mol_j : Chem.Mol
      The second molecule used to calculate the ECR score.

    Returns
    -------
    float
      1.0 if ``mol_i`` and ``mol_j`` have the same total formal charge,
      0.0 otherwise.
    """
    total_charge_mol_i = formal_charge(mol_i)
    total_charge_mol_j = formal_charge(mol_j)

    if abs(total_charge_mol_j - total_charge_mol_i) < 1e-3:
        scr_ecr = 1.0
    else:
        scr_ecr = 0.0

    return scr_ecr


def _find_common_core(mols: list[Chem.Mol], element_change: bool) -> str:
    """
    Find the common core SMARTS among all input molecules.

    Used to seed pairwise MCS searches and speed up calculations.

    Parameters
    ----------
    mols : list[Chem.Mol]
      Input molecules to search for a common core.
    element_change : bool
      If ``True``, allow element changes when finding the common core.

    Returns
    -------
    str
      SMARTS string for the common core, or an empty string if no common
      core is found within the timeout.
    """
    # strip hydrogens off
    mols2 = [Chem.RemoveHs(m) for m in mols]

    if element_change:
        atom_compare = rdFMCS.AtomCompare.CompareAny
    else:
        atom_compare = rdFMCS.AtomCompare.CompareElements

    res = rdFMCS.FindMCS(
        mols2,
        timeout=60,
        atomCompare=atom_compare,
        bondCompare=rdFMCS.BondCompare.CompareAny,
        matchValences=False,
        ringMatchesRingOnly=True,
        completeRingsOnly=True,
        matchChiralTag=False,
    )

    if res.canceled:  # timeout
        return ""
    else:
        return res.smartsString


class DBMolecules:
    """

    This class is used as a container for all the Molecules

    """

    _list: list[Molecule]

    # Initialization function
    def __init__(
        self,
        directory: str,
        parallel: int = 1,
        verbose: Literal["off", "info", "pedantic"] = "off",
        time: int = 20,
        ecrscore: float = 0.0,
        threed: bool = False,
        max3d: float = 1000.0,
        element_change: bool = True,
        output: bool = False,
        name: str = "out",
        output_no_images: bool = False,
        output_no_graph: bool = False,
        display: bool = False,
        allow_tree: bool = False,
        max: int = 6,
        cutoff: float = 0.4,
        radial: bool = False,
        hub: str | None = None,
        fast: bool = False,
        links_file: str | None = None,
        known_actives_file: str | None = None,
        max_dist_from_actives: int = 2,
        use_common_core: bool = True,
        shift: bool = True,
    ):
        """
        Initialization of the Molecule Database Class.

        Parameters
        ----------
        directory : str
          Path to the directory containing mol2 and/or sdf input files.
        parallel : int, default 1
          Number of processes to use when computing the similarity score matrices.
        verbose : str, default 'off'
          Logging verbosity; one of ``'off'``, ``'info'``, or ``'pedantic'``.
        time : int, default 20
          Maximum time in seconds allowed for each pairwise MCS search.
        ecrscore : float, default 0.0
          Electrostatic score override applied when two molecules have
          different formal charges. A value of 0.0 disables cross-charge
          comparisons.
        threed : bool, default False
          If ``True``, symmetry-equivalent MCSes are filtered to prefer
          the one with the best real-space 3D alignment.
        max3d : float, default 1000.0
          MCS atom pairs further apart than this distance (Angstrom) are
          removed. The default of 1000.0 is effectively no filter.
        element_change : bool, default True
          If ``True``, allow element changes between mapped atoms.
        output : bool, default False
          If ``True``, write output files (graph, images, pickle).
        name : str, default 'out'
          File name prefix used when writing output files.
        output_no_images : bool, default False
          If ``True``, skip generation of 2D image output files.
        output_no_graph : bool, default False
          If ``True``, skip generation of the graph ``.dot`` file.
        display : bool, default False
          If ``True``, display the generated graph with Matplotlib.
        allow_tree : bool, default False
          If ``True``, the final graph does not require a cycle covering
          and may be a tree.
        max : int, default 6
          Maximum diameter of the resulting graph.
        cutoff : float, default 0.4
          Minimum Similarity Score (MSS) used to build the graph.
        radial : bool, default False
          If ``True``, build a radial graph around the hub compound.
        hub : str or None, default None
          Name of the molecule to use as the hub in a radial graph.
        fast : bool, default False
          If ``True``, use fast graph building when creating a radial graph.
          If ``radial`` is not ``True`` and ``hub`` is not set, this
          will be ignored.
        links_file : str or None, default None
          Path to a file listing molecule pairs that should be seeded as
          links in the graph.
        known_actives_file : str or None, default None
          Path to a file listing molecules whose activity is known.
        max_dist_from_actives : int, default 2
          Maximum number of graph edges separating any molecule from a
          known active.
        use_common_core : bool, default True
          If ``True``, search for a common core across all input molecules
          to seed and speed up pairwise MCS calculations.
        shift : bool, default True
          If ``True``, when ``threed`` is also ``True``, translate molecules
          to maximise 3D overlap before evaluating the alignment.
        """
        # Set the Logging
        if verbose == "off":
            logging.basicConfig(format="%(message)s", level=logging.CRITICAL)
        elif verbose == "info":
            logging.basicConfig(format="%(message)s", level=logging.INFO)
        elif verbose == "pedantic":
            logging.basicConfig(format="%(message)s", level=logging.DEBUG)
            # logging.basicConfig(format='%(levelname)s:\t%(message)s', level=logging.DEBUG)

        if not isinstance(output, bool):
            raise TypeError("The output flag is not a bool type")
        elif not isinstance(output_no_images, bool):
            raise TypeError("The output_no_images flag is not a bool type")
        elif not isinstance(output_no_graph, bool):
            raise TypeError("The output_no_graph flag is not a bool type")
        elif not isinstance(display, bool):
            raise TypeError("The display flag is not a bool type")
        elif not isinstance(radial, bool):
            raise TypeError("The radial flag is not a bool type")

        self.options: dict[str, Any] = dict()
        CheckDir._check_directory(directory)
        self.options["directory"] = directory
        CheckPos._check(parallel)
        self.options["parallel"] = parallel
        self.options["verbose"] = verbose

        # MCS settings
        CheckPos._check(time)
        self.options["time"] = time
        CheckEcrscore._check(ecrscore)
        self.options["ecrscore"] = ecrscore
        self.options["threed"] = bool(threed)
        self.options["shift"] = bool(shift)
        self.options["max3d"] = max3d
        self.options["element_change"] = bool(element_change)

        # Output settings
        self.options["output"] = bool(output)
        self.options["name"] = name
        self.options["output_no_images"] = bool(output_no_images)
        self.options["output_no_graph"] = bool(output_no_graph)
        self.options["display"] = bool(display)

        # Graph settings
        self.options["allow_tree"] = bool(allow_tree)
        CheckPos._check(max)
        self.options["max"] = max
        CheckPos._check(max_dist_from_actives)
        self.options["max_dist_from_actives"] = max_dist_from_actives
        CheckCutoff._check(cutoff)
        self.options["cutoff"] = cutoff
        self.options["radial"] = bool(radial)
        self.options["hub"] = str(hub) if hub is not None else hub
        self.options["fast"] = bool(fast)
        self.options["links_file"] = links_file
        self.options["known_actives_file"] = known_actives_file

        # Internal list container used to store the loaded molecule objects
        self._list = self.read_molecule_files()

        if use_common_core:
            self.options["seed"] = _find_common_core(
                [m.getMolecule() for m in self._list], self.options["element_change"]
            )
        else:
            self.options["seed"] = ""

        # Dictionary which holds the mapping between the generated molecule IDs and molecule file names
        self.dic_mapping: dict[int, str] = {}
        self.inv_dic_mapping: dict[str, int] = {}

        # Hold the MCS index map strings for each molecule pair. Indexed by a tuple of molecule IDs (lowest first)
        self.mcs_map_store: dict[tuple[int, int], str] = {}

        # Pre-specified links between molecules - a map of molecule index tuples to score.
        # A value < -1 means "recompute, but force the link to be included"
        # A negative value (>= -1) means "Use the absolute value as the score, but force the link to be included"
        # A positive value means "Use this value as the score, but treat the link as normal in the graph calculation"
        self.prespecified_links: dict[tuple[int, int], float] = {}

        # List of which molecules are "known actives". Note that all pairs of known actives
        # are automatically added as prespecified links with a score of -1 (i.e force score to
        # 1 and force link to be included)
        self.known_actives: list[int] = []

        for mol in self._list:
            self.dic_mapping[mol.getID()] = mol.getName()
            self.inv_dic_mapping[mol.getName()] = mol.getID()

        if self.options["links_file"]:
            self.parse_links_file(self.options["links_file"])

        if self.options["known_actives_file"]:
            self.parse_known_actives_file(self.options["known_actives_file"])

        # Index used to perform index selection by using __iter__ function
        self._ci = 0

        # Symmetric matrices used to store the mcs scoring. The matrices are subclasses of numpy
        self.strict_mtx = SMatrix(shape=(0,))
        self.loose_mtx = SMatrix(shape=(0,))

        # Empty pointer to the networkx graph
        self.Graph: nx.Graph = nx.Graph()

    def __iter__(self) -> Iterator[Molecule]:
        """
        Iterate through the database.

        Returns
        -------
        Iterator[Molecule]
          This database instance as an iterator.
        """
        return self

    def __next__(self) -> Molecule:
        return self.next()

    def next(self) -> Molecule:
        """
        Return the next molecule in the database sequence.

        Returns
        -------
        Molecule
          The next molecule in the database.

        Raises
        ------
        StopIteration
          When all molecules have been yielded.
        """

        if self._ci > len(self._list) - 1:
            self._ci = 0
            raise StopIteration
        else:
            self._ci = self._ci + 1
            return self._list[self._ci - 1]

    def __getitem__(self, index: int) -> Molecule:
        """
        Return the molecule at the given index.

        Parameters
        ----------
        index : int
          Zero-based position in the molecule list.

        Returns
        -------
        Molecule
          The molecule at position ``index``.
        """

        return self._list[index]

    def __setitem__(self, index: int, molecule: Molecule) -> None:
        """
        Replace the molecule at the given index.

        Parameters
        ----------
        index : int
          Zero-based position in the molecule list.
        molecule : Molecule
          The molecule to store at position ``index``.

        Raises
        ------
        ValueError
          If ``molecule`` is not a :class:`Molecule` instance.
        """

        if not isinstance(molecule, Molecule):
            raise ValueError("The passed molecule is not a Molecule object")

        self._list[index] = molecule

    def __add__(self, molecule: Molecule) -> None:
        """
        Append a molecule to the database.

        Parameters
        ----------
        molecule : Molecule
          The molecule to append.

        Raises
        ------
        ValueError
          If ``molecule`` is not a :class:`Molecule` instance.
        """

        if not isinstance(molecule, Molecule):
            raise ValueError("The passed molecule is not a Molecule object")

        self._list.append(molecule)

    def nums(self) -> int:
        """
        Return the total number of molecules in the database.

        Returns
        -------
        int
          Number of molecules currently stored.
        """
        return len(self._list)

    def read_molecule_files(self) -> list[Molecule]:
        """
        Read all mol2 and SDF files from the configured directory.

        Returns
        -------
        list[Molecule]
          All successfully loaded molecules, in sorted filename order.

        Raises
        ------
        OSError
          If the directory contains fewer than two mol2/sdf files.
        """

        # This list is used as container to handle all the molecules read in by using RdKit.
        # All the molecules are instances of  Molecule class
        molid_list = []

        # List of molecule that failed to load in
        mol_error_list_fn: list[str] = []

        logging.info(30 * "-")

        # The .mol2 and .sdf file formats are the only supported so far
        mol_fnames = glob.glob(self.options["directory"] + "/*.mol2")
        mol_fnames += glob.glob(self.options["directory"] + "/*.sdf")

        mol_fnames.sort()

        if len(mol_fnames) < 2:
            raise OSError(
                f"The directory {self.options['directory']} must contain at least two mol2/sdf files"
            )

        print_cnt = 0
        mol_id_cnt = 0

        for fname in mol_fnames:
            # The RDkit molecule object reads in as mol2/sdf file. The molecule is not sanitized and
            # all the hydrogens are kept in place - we are assuming 3D input, correctly charged
            # and prepared in the protein active site
            if fname.endswith(".mol2"):
                rdkit_mol = Chem.MolFromMol2File(fname, sanitize=False, removeHs=False)
            else:
                rdkit_mol = Chem.MolFromMolFile(fname, sanitize=False, removeHs=False)

            # Reading problems
            if rdkit_mol is None:
                logging.warning(f"Error reading the file: {os.path.basename(fname)}")
                mol_error_list_fn.append(os.path.basename(fname))
                continue

            # The Rdkit molecule is stored in a Molecule object
            mol = Molecule(rdkit_mol, mol_id_cnt, os.path.basename(fname))
            mol_id_cnt += 1

            # Cosmetic printing and status
            if print_cnt < 15 or print_cnt == (len(mol_fnames) - 1):
                logging.info(f"ID {mol.getID()}\t{os.path.basename(fname)}")

            if print_cnt == 15:
                logging.info(f"ID {mol.getID()}\t{os.path.basename(fname)}")
                logging.info(3 * "\t.\t.\n")

            print_cnt += 1

            molid_list.append(mol)

        logging.info(30 * "-")

        logging.info(
            f"Finish reading input files. {len(molid_list)} structures in total....skipped {len(mol_error_list_fn)}\n"
        )

        if mol_error_list_fn:
            logging.warning("Skipped molecules:")
            logging.warning(30 * "-")
            for fn in mol_error_list_fn:
                logging.warning(str(fn))
            print(30 * "-")

        return molid_list

    def parse_links_file(self, links_file: str) -> None:
        """
        Parse a links file and register the specified molecule pairs.

        Each line may have the form ``mol1 mol2``, ``mol1 mol2 score``, or
        ``mol1 mol2 score force``.  The ``force`` keyword causes the link to be
        included in the final graph regardless of its score.

        Parameters
        ----------
        links_file : str
          Path to the links file.

        Raises
        ------
        OSError
          If the file has a syntax error or references an unknown molecule file.
        """
        try:
            with open(links_file) as lf:
                for line in lf:
                    mols = line.split()
                    if len(mols) < 2 or len(mols) > 4:
                        raise OSError("Syntax error in links file parsing line:" + line)
                    indexa = self.inv_dic_mapping[mols[0]]
                    indexb = self.inv_dic_mapping[mols[1]]
                    score = -2.0
                    if len(mols) > 2:
                        score = float(mols[2])
                    if len(mols) > 3:
                        if mols[3] != "force":
                            raise OSError(
                                "Syntax error parsing fourth argument in links file on line:" + line
                            )
                        score = -score
                    self.prespecified_links[(indexa, indexb)] = score
                    self.prespecified_links[(indexb, indexa)] = score
                    print(
                        "Added prespecified link for mols",
                        mols,
                        "->",
                        (indexa, indexb),
                        "score",
                        score,
                    )
        except KeyError as e:
            raise OSError(
                'Filename within the links file "' + links_file + '" not found: ' + str(e)
            ) from None

    def parse_known_actives_file(self, actives_file: str) -> None:
        """
        Parse a known-actives file and mark the listed molecules as active.

        All pairs of known actives are automatically added as prespecified
        links with a score of 1.0 forced into the graph.

        Parameters
        ----------
        actives_file : str
          Path to the known-actives file (one molecule filename per line).

        Raises
        ------
        OSError
          If the file references an unknown molecule file name.
        """
        try:
            with open(actives_file) as lf:
                for line in lf:
                    mols = line.split()
                    indexa = self.inv_dic_mapping[mols[0]]
                    self.known_actives.append(indexa)
                    self._list[indexa].setActive(True)
                    logging.info(f"Added known activity for mol {mols[0]} -> {indexa}")
        except KeyError as e:
            raise OSError(
                'Filename within the actives file "' + actives_file + '" not found: ' + str(e)
            ) from None
        # Add all combinations of these to the set of prespecified links
        for t in [(x, y) for x in self.known_actives for y in self.known_actives]:
            logging.info(f"Added prespecified link for {t}")
            self.prespecified_links[t] = -1

    def set_MCSmap(self, i: int, j: int, MCmap: str) -> None:
        """
        Store the MCS atom-index map string for a molecule pair.

        Parameters
        ----------
        i : int
          Index of the first molecule.
        j : int
          Index of the second molecule.
        MCmap : str
          Serialised MCS atom-index map between the two molecules.

        Notes
        -----
        The pair is stored with the lower index first so that
        ``set_MCSmap(i, j, ...)`` and ``set_MCSmap(j, i, ...)`` address the
        same entry.
        """
        if i < j:
            idx = (i, j)
        else:
            idx = (j, i)
        self.mcs_map_store[idx] = MCmap

    def get_MCSmap(self, i: int, j: int) -> str | None:
        """
        Retrieve the MCS atom-index map string for a molecule pair.

        Parameters
        ----------
        i : int
          Index of the first molecule.
        j : int
          Index of the second molecule.

        Returns
        -------
        str or None
          The stored MCS atom-index map string, or ``None`` if no map has
          been stored for this pair.
        """
        if i < j:
            idx = (i, j)
        else:
            idx = (j, i)
        if idx in self.mcs_map_store:
            return self.mcs_map_store[idx]
        return None

    def compute_mtx(
        self,
        a: int,
        b: int,
        strict_mtx: Any,
        loose_mtx: Any,
        true_strict_mtx: Any,
        MCS_map: dict[tuple[int, int], str] | Any,
    ) -> None:
        """
        Compute a chunk of the similarity score matrices.

        The chunk spans linear indices ``a`` to ``b`` (inclusive) of the
        upper-triangle of the symmetric score matrix stored as a flat array.

        Parameters
        ----------
        a : int
          Start index of the chunk to compute.
        b : int
          End index (inclusive) of the chunk to compute.
        strict_mtx : SMatrix or multiprocessing.Array
          Flat array for strict similarity scores, shared across processes.
          Can be ``multiprocessing.Array``.
        loose_mtx : SMatrix | Any
          Flat array for loose similarity scores, shared across processes.
          Can be ``multiprocessing.Array``.
        true_strict_mtx : SMatrix | Any
          Flat array that stores the strict score before any forced-link
          override is applied. Can be ``multiprocecssing.Array``.
        MCS_map : dict[tuple[int, int], str] | Any
          Mapping from molecule-index pairs to their MCS atom-index map
          strings, shared across processes. Can be multiprocessing dict.
        """

        # name = multiprocessing.current_process().name
        # print(name)
        # print('a = %d, b = %d' % (a,b))
        # print('\n')

        # Total number of loaded molecules
        n = self.nums()

        # Looping over all the elements of the selected matrix chunk
        for k in range(a, b + 1):
            # The linear index k is converted into the row and column indexes of
            # a hypothetical bidimensional symmetric matrix
            i = int(n - 2 - math.floor(math.sqrt(-8 * k + 4 * n * (n - 1) - 7) / 2.0 - 0.5))
            j = int(k + i + 1 - n * (n - 1) / 2 + (n - i) * ((n - i) - 1) / 2)
            # print 'k = %d , i = %d , j = %d' % (k,i,j)

            # The Rdkit molecules moli and molj are extracted from the molecule database
            moli = self[i].getMolecule()
            molj = self[j].getMolecule()

            logging.info(f"Processing molecules: {self[i].getName()}-{self[j].getName()}")

            # The Electrostatic score rule is calculated
            ecr_score = ecr(moli, molj)

            # If the prespecified links map has this link, and the value is >-1, then
            # we don't need to compute the score
            if (i, j) in self.prespecified_links and self.prespecified_links[(i, j)] >= -1:
                strict_scr = abs(self.prespecified_links[(i, j)])
                loose_scr = strict_scr
                logging.info(
                    f"MCS molecules: {self[i].getName()} - {self[j].getName()} "
                    f"final score {strict_scr} set in links file"
                )
            else:
                # The MCS is computed only if the passed molecules have the same charges
                if ecr_score or self.options["ecrscore"]:
                    if ecr_score == 0.0 and self.options["ecrscore"]:
                        logging.critical(
                            "WARNING: Mutation between different charge molecules is enabled"
                        )
                        ecr_score = self.options["ecrscore"]

                    try:
                        if self.options["verbose"] == "pedantic":
                            logging.info(50 * "-")
                            logging.info(
                                f"MCS molecules: {self[i].getName()} - {self[j].getName()}"
                            )

                        # Maximum Common Subgraph (MCS) calculation
                        MC = mcs.MCS(
                            moli,
                            molj,
                            time=self.options["time"],
                            verbose=self.options["verbose"],
                            threed=self.options["threed"],
                            max3d=self.options["max3d"],
                            element_change=self.options["element_change"],
                            seed=self.options["seed"],
                            shift=self.options["shift"],
                        )
                        ml = MC.all_atom_match_list()
                        MCS_map[(i, j)] = ml

                    except Exception as e:
                        logging.warning(
                            f"Skipping MCS molecules (exception): {self[i].getName()} - {self[j].getName()}\t\n\n"
                            f"{e}"
                        )
                        logging.info(50 * "-")
                        continue
                else:
                    continue

                # The scoring between the two molecules is performed by using different rules.
                # The total score will be the product of all the single rules
                tmp_scr = (
                    ecr_score
                    * MC.mncar()
                    * MC.mcsr()
                    * MC.atomic_number_rule()
                    * MC.hybridization_rule()
                )
                tmp_scr *= (
                    MC.sulfonamides_rule()
                    * MC.heterocycles_rule()
                    * MC.transmuting_methyl_into_ring_rule()
                )
                tmp_scr *= MC.transmuting_ring_sizes_rule()
                # Note - no longer using tmcsr rule!
                strict_scr = tmp_scr * 1  # MC.tmcsr(strict_flag=True)
                loose_scr = tmp_scr * 1  # MC.tmcsr(strict_flag=False)

            strict_mtx[k] = strict_scr
            loose_mtx[k] = loose_scr
            true_strict_mtx[k] = strict_scr

            # process prespecified links now and overwrite the existing info
            if (i, j) in self.prespecified_links and self.prespecified_links[(i, j)] < 0:
                print(f"Molecule pair {i} {j} forced to be included in the graph - score set to 1")
                strict_mtx[k] = 1.0
                loose_mtx[k] = 1.0
                # Note that true_strict_mtx holds the original strict_scr value
                continue

        return

    def build_matrices(self) -> tuple[SMatrix, SMatrix]:
        """
        Compute the pairwise similarity score matrices for all loaded molecules.

        Work will be distributed between ``self.parallel`` processes.

        Returns
        -------
        tuple[SMatrix, SMatrix]
          A pair ``(strict_mtx, loose_mtx)`` of symmetric score matrices stored
          as flat :class:`SMatrix` arrays.
        """

        logging.info("\nMatrix scoring in progress....\n")

        # The similarity score matrices are defined instances of the class SMatrix
        # which implements a basic class for symmetric matrices
        self.strict_mtx = SMatrix(shape=(self.nums(),))
        self.loose_mtx = SMatrix(shape=(self.nums(),))
        self.true_strict_mtx = SMatrix(shape=(self.nums(),))
        # The total number of the effective elements present in the symmetric matrix
        elems = int(self.nums() * (self.nums() - 1) / 2)

        if self.options["parallel"] == 1:  # Serial execution
            MCS_map: dict[tuple[int, int], str] = {}
            self.compute_mtx(
                0, elems - 1, self.strict_mtx, self.loose_mtx, self.true_strict_mtx, MCS_map
            )
            for idx in MCS_map:
                self.set_MCSmap(idx[0], idx[1], MCS_map[idx])
        else:
            # Parallel execution
            logging.info("Parallel mode is on")

            # Number of selected processes
            num_proc = self.options["parallel"]

            delta = int(elems / num_proc)
            rem = elems % num_proc

            if delta < 1:
                kmax = elems
            else:
                kmax = num_proc
            proc = []

            with multiprocessing.Manager() as manager:
                # Shared memory array used by the different allocated processes
                # At the moment we're using a combination of Array and Manager, which is nasty
                # Note: Ignore the call-overload typing issues for the next three. SMatrix & multiprocessing has odd typing.
                strict_mtx = multiprocessing.Array("d", self.strict_mtx)  # type: ignore[call-overload]
                loose_mtx = multiprocessing.Array("d", self.loose_mtx)  # type: ignore[call-overload]
                true_strict_mtx = multiprocessing.Array("d", self.true_strict_mtx)  # type: ignore[call-overload]
                MCS_map = manager.dict()  # type: ignore[assignment]

                # Chopping the indexes redistributing the remainder
                for k in range(0, kmax):
                    spc = delta + int(int(rem / (k + 1)) > 0)

                    if k == 0:
                        i = 0
                    else:
                        # Note: this loop structure assumes j will be defined
                        # in the first iteration where k == 0, it's not ideal and linters
                        # don't like it, but it should work
                        i = j + 1  # type: ignore[has-type] # noqa: F821

                    if k != kmax - 1:
                        j = i + spc - 1
                    else:
                        j = elems - 1

                    # Python multiprocessing allocation
                    p = multiprocessing.Process(
                        target=self.compute_mtx,
                        args=(
                            i,
                            j,
                            strict_mtx,
                            loose_mtx,
                            true_strict_mtx,
                            MCS_map,
                        ),
                    )
                    p.start()
                    proc.append(p)
                # End parallel execution
                for p in proc:
                    p.join()

                # Copying back the results
                self.strict_mtx[:] = strict_mtx[:]
                self.loose_mtx[:] = loose_mtx[:]
                self.true_strict_mtx[:] = true_strict_mtx[:]
                for idx in MCS_map.keys():
                    self.set_MCSmap(idx[0], idx[1], MCS_map[idx])

        return self.strict_mtx, self.loose_mtx

    def build_graph(self) -> nx.Graph:
        """
        Build the perturbation network graph from the computed score matrices.

        Optionally writes output files and displays the graph.

        Returns
        -------
        nx.Graph
          The generated perturbation network graph.
        """
        logging.info("\nGenerating graph in progress....")

        # The Graph is build from an instance of the Class GraphGen by passing
        # the selected user options
        Gr = graphgen.GraphGen(
            score_matrix=self.strict_mtx.to_numpy_2D_array(),
            ids=[m.getID() for m in self],
            names=[m.getName() for m in self],
            actives=[m.isActive() for m in self],
            max_path_length=self.options["max"],
            max_dist_from_active=self.options["max_dist_from_actives"],
            similarity_cutoff=self.options["cutoff"],
            require_cycle_covering=not self.options["allow_tree"],
            radial=self.options["radial"],
            fast=self.options["fast"],
            hub=self.options["hub"],
        )

        # Writing the results to files
        if self.options["output"]:
            try:
                Gr.write_graph(
                    self, self.options["output_no_images"], self.options["output_no_graph"]
                )
                with open(self.options["name"] + ".pickle", "wb") as pickle_f:
                    pickle.dump(Gr, pickle_f)
            except Exception as e:
                logging.error(str(e))

        # Handle to the NetworkX generated graph
        self.Graph = Gr.resultGraph

        # print self.Graph.nodes(data=True)

        # Display the graph by using Matplotlib
        if self.options["display"]:
            Gr.draw(self)

        return self.Graph

    def write_dic(self) -> None:
        """
        Write the molecule index-to-filename mapping to a text file.

        The output file is named ``<name>.txt`` where ``name`` is the prefix
        set at class construction time.

        Raises
        ------
        OSError
          If the output file cannot be opened for writing.
        """

        try:
            file_txt = open(self.options["name"] + ".txt", "w")
        except Exception:
            raise OSError("It was not possible to write out the mapping file")
        file_txt.write("#ID\tFileName\n")
        for key in self.dic_mapping:
            file_txt.write(f"{key}\t{self.dic_mapping[key]}\n")

        file_txt.close()


class SMatrix(np.ndarray):
    """
    This class implements a "basic" interface for symmetric matrices
    subclassing ndarray. The class internally stores a bi-dimensional
    numpy array as a linear array A[k], however the user can still
    access to the matrix elements by using a two indices notation A[i,j]

    """

    def __new__(
        subtype: type[SMatrix],
        shape: tuple[int, ...],
        dtype: type = float,
        buffer: np.ndarray | None = None,
        offset: int = 0,
        strides: tuple[int, ...] | None = None,
        order: Literal['K', 'A', 'C', 'F'] | None = None,
    ) -> SMatrix:
        if len(shape) > 2:
            raise ValueError("The matrix shape is greater than two")

        elif len(shape) == 2:
            if shape[0] != shape[1]:
                raise ValueError("The matrix must be a square matrix")

        elems = int(shape[0] * (shape[0] - 1) / 2)

        shape = (elems,)

        obj = np.ndarray.__new__(subtype, shape, dtype, buffer, offset, strides, order)

        # Array initialization
        # Note: ignore assignment typing issue - SMatrix is oddly typed
        obj = obj * 0.0  # type: ignore[assignment]

        return obj

    def __getitem__(self, *kargs: Any) -> float | np.ndarray:  # type: ignore[override]
        """
        Retrieve one or more elements from the symmetric matrix.

        Supports three calling forms:

        * ``A[k]`` — linear index into the flat storage array.
        * ``A[i:j]`` — slice of the flat storage array.
        * ``A[i, j]`` — element at row ``i``, column ``j`` of the square
          matrix (returns 0.0 when ``i == j``).

        Parameters
        ----------
        *kargs : int or slice or tuple[int, int]
          Index, slice, or pair of row/column indices.

        Returns
        -------
        float or np.ndarray
          Scalar element or array slice from the flat storage array.

        Raises
        ------
        ValueError
          If more than two matrix indices are provided, or if an index is
          out of bounds.
        """

        if isinstance(kargs[0], int):
            k = kargs[0]
            return super().__getitem__(k)

        if isinstance(kargs[0], slice):
            k = kargs[0]  # type: ignore[assignment]
            return super().__getitem__(k)

        elif len(kargs[0]) > 2:
            raise ValueError("Two indices can be addressed")

        i = kargs[0][0]
        j = kargs[0][1]

        if i == j:
            return 0.0

        # Total number of elements in the corresponding bi-dimensional symmetric matrix
        # where self.size is the length of the linear array
        n = int((1 + math.sqrt(1 + 8 * self.size)) / 2)

        if i > n - 1:
            raise ValueError("First index out of bound")

        if j > n - 1:
            raise ValueError("Second index out of bound")

        if i < j:
            k = int((n * (n - 1) / 2) - (n - i) * ((n - i) - 1) / 2 + j - i - 1)
        else:
            k = int((n * (n - 1) / 2) - (n - j) * ((n - j) - 1) / 2 + i - j - 1)

        return super().__getitem__(k)

    def __setitem__(self, *kargs: Any) -> None:
        """
        Set one or more elements in the symmetric matrix.

        Supports three calling forms:

        * ``A[k] = v`` — set a single element by linear index.
        * ``A[i:j] = v`` — set a slice of the flat storage array.
        * ``A[i, j] = v`` — set the element at row ``i``, column ``j``.

        Parameters
        ----------
        *kargs : int or slice or tuple[int, int], followed by the value
          Index/slice and the value to assign.

        Raises
        ------
        ValueError
          If more than two matrix indices are provided, or if an index is
          out of bounds.
        """

        if isinstance(kargs[0], int):
            k = kargs[0]
            value = kargs[1]
            return super().__setitem__(k, value)

        elif isinstance(kargs[0], slice):
            start, stop, step = kargs[0].indices(len(self))
            value = kargs[1]
            return super().__setitem__(kargs[0], value)

        elif len(kargs[0]) > 2:
            raise ValueError("Two indices can be addressed")

        # Passed indexes and value to set
        i = kargs[0][0]
        j = kargs[0][1]
        value = kargs[1]

        # Total number of elements in the corresponding bi-dimensional symmetric matrix
        # where self.size is the length of the linear array
        n = int((1 + math.sqrt(1 + 8 * self.size)) / 2)

        if i > n - 1:
            raise ValueError("First index out of bound")
        if j > n - 1:
            raise ValueError("Second index out of bound")

        if i < j:
            k = int((n * (n - 1) / 2) - (n - i) * ((n - i) - 1) / 2 + j - i - 1)
        else:
            k = int((n * (n - 1) / 2) - (n - j) * ((n - j) - 1) / 2 + i - j - 1)
        super().__setitem__(k, value)

    def to_numpy_2D_array(self) -> np.ndarray:
        """
        Return a 2D numpy arrray of the symmetric similarity score from
        the flat storage array.

        Returns
        -------
        np_mat : np.ndarray
          Square symmetric matrix of shape ``(n, n)`` where ``n`` is the
          number of molecules.
        """

        # Total number of elements in the corresponding bi-dimensional symmetric matrix
        # where self.size is the length of the linear array
        n = int((1 + math.sqrt(1 + 8 * self.size)) / 2)

        np_mat = np.zeros((n, n))

        for i in range(0, n):
            for j in range(0, n):
                np_mat[i, j] = self[i, j]

        return np_mat

    def mat_size(self) -> int:
        """
        Return the side length of the equivalent square matrix.

        Returns
        -------
        int
          Number of rows (and columns) in the corresponding square matrix.
        """

        # Total number of elements in the corresponding bi-dimensional symmetric matrix
        # where self.size is the length of the linear array
        n = int((1 + math.sqrt(1 + 8 * self.size)) / 2)

        return n


class Molecule:
    """
    This Class stores the Rdkit molecule objects, their identification number
    and the total number of instantiated molecules

    """

    def __init__(self, molecule: Chem.rdchem.Mol, mol_id: int, molname: str) -> None:
        """
        Parameters
        ----------
        molecule : Chem.rdchem.Mol
          The RDKit molecule object to wrap.
        mol_id : int
          Unique integer identifier for this molecule.
        molname : str
          File name (basename) associated with this molecule.

        Raises
        ------
        ValueError
          If ``molecule`` is not an RDKit molecule or ``molname`` is not a
          string.
        """

        # Check Inputs
        if not isinstance(molecule, Chem.rdchem.Mol):
            raise ValueError("The passed molecule object is not a RdKit molecule")

        if not isinstance(molname, str):
            raise ValueError("The passed molecule name must be a string")

        # The variable __molecule saves the current RDkit molecule object
        # The variable is defined as private
        self.__molecule = molecule

        # The variable __ID saves the molecule identification number
        # The variable is defined as private
        self.__ID = mol_id

        # The variable __name saves the molecule identification name
        # The variable is defined as private
        self.__name = molname

        # The variable __active saves whether the molecule is a known active
        # The variable is defined as private
        self.__active = False

    def getID(self) -> int:
        """
        Return the molecule's integer identifier.

        Returns
        -------
        int
          Unique identifier assigned to this molecule.
        """
        return self.__ID

    def getMolecule(self) -> Chem.Mol:
        """
        Return a copy of the underlying RDKit molecule object.

        Returns
        -------
        Chem.Mol
          A copy of the stored RDKit molecule.
        """
        mol_copy = Chem.Mol(self.__molecule)
        return mol_copy

    def getName(self) -> str:
        """
        Return the molecule's file name.

        Returns
        -------
        str
          Basename of the file from which this molecule was loaded.
        """
        return self.__name

    def isActive(self) -> bool:
        """
        Return whether this molecule is a known active.

        Returns
        -------
        bool
          ``True`` if this molecule has been marked as a known active.
        """
        return self.__active

    def setActive(self, active: bool) -> None:
        """
        Mark or unmark this molecule as a known active.

        Parameters
        ----------
        active : bool
          Pass ``True`` to mark the molecule as a known active, ``False``
          to unmark it.
        """
        self.__active = active


class CheckDir(argparse.Action):
    # Classes used to check some of the passed user options in the main function
    # Class used to check the input directory
    @classmethod
    def _check_directory(cls, directory: str) -> None:
        if not os.path.isdir(directory):
            raise argparse.ArgumentTypeError(f"The directory name is not a valid path: {directory}")
        if not os.access(directory, os.R_OK):
            raise argparse.ArgumentTypeError(f"The directory name is not readable: {directory}")

    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        directory: str,  # type: ignore[override]
        option_string: str | None = None,
    ) -> None:
        self._check_directory(directory)
        setattr(namespace, self.dest, directory)


class CheckPos(argparse.Action):
    # Class used to check the parallel, time and max user options
    @classmethod
    def _check(cls, value: int) -> None:
        if value < 1:
            raise argparse.ArgumentTypeError(f"{value} is not a positive integer number")

    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        value: int,  # type: ignore[override]
        option_string: str | None = None,
    ) -> None:
        self._check(value)
        setattr(namespace, self.dest, value)


class CheckCutoff(argparse.Action):
    # Class used to check the cutoff user option
    @classmethod
    def _check(cls, value: float) -> None:
        if not isinstance(value, float) or value < 0.0:
            raise argparse.ArgumentTypeError(f"{value} is not a positive real number")

    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        value: float,  # type: ignore[override]
        option_string: str | None = None,
    ) -> None:
        self._check(value)
        setattr(namespace, self.dest, value)


class CheckEcrscore(argparse.Action):
    # Class used to check the handicap user option
    @classmethod
    def _check(cls, value: float) -> None:
        if not isinstance(value, float) or value < 0.0 or value > 1.0:
            raise argparse.ArgumentTypeError(
                f"{value} is not a real number in the range [0.0, 1.0]"
            )

    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        value: float,  # type: ignore[override]
        option_string: str | None = None,
    ) -> None:
        self._check(value)
        setattr(namespace, self.dest, value)


def startup() -> None:
    # Emit user-facing warning that the CLI is deprecated and will be removed
    warnings.warn(
        "The dbmol CLI is deprecated and will be removed in the next major release. "
        "Please let us know if you keeping this CLI entry point is important to you, "
        "see https://github.com/OpenFreeEnergy/Lomap/issues/138 for more details.",
        category=DeprecationWarning,
    )
    # This is the CLI entrypoint
    # Options and arguments passed by the user
    ops = parser.parse_args()

    _startup_inner(
        directory=ops.directory,
        parallel=ops.parallel,
        verbose=ops.verbose,
        time=ops.time,
        ecrscore=ops.ecrscore,
        threed=ops.threed,
        max3d=ops.max3d,
        element_change=ops.element_change,
        output=ops.output,
        name=ops.name,
        output_no_images=ops.output_no_images,
        output_no_graph=ops.output_no_graph,
        display=ops.display,
        allow_tree=ops.allow_tree,
        max=ops.max,
        max_dist_from_actives=ops.max_dist_from_actives,
        cutoff=ops.cutoff,
        radial=ops.radial,
        hub=ops.hub,
        fast=ops.fast,
        links_file=ops.links_file,
        known_actives_file=ops.known_actives_file,
        common_core=ops.common_core,
    )


def _startup_inner(
    directory: str,  # TODO: Should really constant out the CLI constants to keep this DRY
    parallel: int = 1,
    verbose: Literal['off', 'info', 'pedantic'] = "info",
    time: int = 20,
    ecrscore: float = 0.0,
    threed: bool = False,
    max3d: float = 1000,
    element_change: bool = True,
    output: bool = True,
    name: str = "out",
    output_no_images: bool = False,
    output_no_graph: bool = False,
    display: bool = False,
    allow_tree: bool = False,
    max: int = 6,
    max_dist_from_actives: int = 2,
    cutoff: float = 0.4,
    radial: bool = False,
    hub: str | None = None,
    fast: bool = False,
    links_file: str = "",
    known_actives_file: str = "",
    common_core: bool = True,
) -> None:
    # Inside function of CLI interface, for start of "library" like calling

    # Molecule DataBase initialized with the passed user options
    db_mol = DBMolecules(
        directory,
        parallel,
        verbose,
        time,
        ecrscore,
        threed,
        max3d,
        element_change,
        output,
        name,
        output_no_images,
        output_no_graph,
        display,
        allow_tree,
        max,
        cutoff,
        radial,
        hub,
        fast,
        links_file,
        known_actives_file,
        max_dist_from_actives,
        use_common_core=common_core,
    )
    # Similarity score linear array generation
    strict, loose = db_mol.build_matrices()

    # Get the 2D numpy matrices
    # strict.to_numpy_2D_array()
    # loose.to_numpy_2D_array()

    # Graph generation based on the similarity score matrix
    _ = db_mol.build_graph()

    # print db_mol.Graph.nodes(data=True)
    # print db_mol.Graph.edges(data=True)


# Command line user interface
# ----------------------------------------------------------------
parser = argparse.ArgumentParser(
    description="Lead Optimization Mapper 2. A program to plan alchemical relative "
    "binding affinity calculations",
    prog=f"LOMAP v. {lomap.__version__}",  # type: ignore[has-type]
)
parser.add_argument("directory", action=CheckDir, help="The mol2/sdf file directory")
parser.add_argument(
    "-p",
    "--parallel",
    default=1,
    action=CheckPos,
    type=int,
    help="Set the parallel mode. If an integer number N is specified, N processes will be executed to "
    "build the similarity matrices",
)
parser.add_argument(
    "-v",
    "--verbose",
    default="info",
    type=str,
    choices=["off", "info", "pedantic"],
    help="verbose mode selection",
)

mcs_group = parser.add_argument_group("MCS setting")
mcs_group.add_argument(
    "-t",
    "--time",
    default=20,
    action=CheckPos,
    type=int,
    help="Set the maximum time in seconds to perform the mcs search between pair of molecules",
)
mcs_group.add_argument(
    "-e",
    "--ecrscore",
    default=0.0,
    action=CheckEcrscore,
    type=float,
    help="If different from 0.0 the value is used to set the electrostatic score between two molecules with different charges",
)
mcs_group.add_argument(
    "-3",
    "--threed",
    default=False,
    action="store_true",
    help="Use the input 3D coordinates to guide the preferred MCS mappings",
)
mcs_group.add_argument(
    "-x",
    "--max3d",
    default=1000,
    type=float,
    help="The MCS is trimmed to remove atoms which are further apart than this distance",
)
mcs_group.add_argument(
    "-s",
    "--shift",
    default=True,
    action="store_true",
    help="Translate molecules to maximise overlap before checking real space alignment for symmetry resolving",
)
mcs_group.add_argument(
    "-L",
    "--element_change",
    default=True,
    type=bool,
    help="Whether to allow element changes in mappings",
)

out_group = parser.add_argument_group("Output setting")
out_group.add_argument(
    "-o", "--output", default=True, action="store_true", help="Generates output files"
)
out_group.add_argument(
    "-n",
    "--name",
    type=str,
    default="out",
    help="File name prefix used to generate the output files",
)
out_group.add_argument(
    "--output-no-images",
    default=False,
    action="store_true",
    help="Disable the generation of the image files, removed the dependency on Pillow",
)
out_group.add_argument(
    "--output-no-graph",
    default=False,
    action="store_true",
    help="Disable the generation of the graph (.dot) file, removed the dependency on pygraphviz",
)

parser.add_argument(
    "-d",
    "--display",
    default=False,
    action="store_true",
    help="Display the generated graph by using Matplotlib",
)

graph_group = parser.add_argument_group("Graph setting")
graph_group.add_argument(
    "-T",
    "--allow-tree",
    default=False,
    action="store_true",
    help="Remove the requirement that all molecules be in a cycle, so that the returned "
    "graph will be a tree instead.",
)
graph_group.add_argument(
    "-m", "--max", default=6, action=CheckPos, type=int, help="The maximum diameter of the graph"
)
graph_group.add_argument(
    "-A",
    "--max-dist-from-actives",
    default=2,
    action=CheckPos,
    type=int,
    help="The maximum distance of any molecule from an active (requires -k)",
)
graph_group.add_argument(
    "-c",
    "--cutoff",
    default=0.4,
    action=CheckCutoff,
    type=float,
    help="The Minimum Similarity Score (MSS) used to build the graph",
)
graph_group.add_argument(
    "-r",
    "--radial",
    default=False,
    action="store_true",
    help="Using the radial option to build the graph",
)
graph_group.add_argument(
    "-b",
    "--hub",
    default=None,
    type=str,
    help="Using a radial graph approach with a manually specified hub compound",
)
graph_group.add_argument(
    "-a",
    "--fast",
    default=False,
    action="store_true",
    help="Using the fast graphing when the lead compound is specified",
)
graph_group.add_argument(
    "-l",
    "--links-file",
    type=str,
    default="",
    help="Specify a filename listing the pairs of molecule files that should be initialised as linked."
    'Each line can be "mol1 mol2", which indicates that Lomap should compute the score and mapping '
    'but that this link must be used in the final graph, or "mol1 mol2 score", which indicates that '
    'Lomap should use the provided score, or "mol1 mol2 score force", which indicates that Lomap '
    "should use the provided score and force this link to be used in the final graph.",
)
graph_group.add_argument(
    "-k",
    "--known-actives-file",
    type=str,
    default="",
    help='Specify a filename listing the molecule files that should be initialised as "known actives", one per line',
)
graph_group.add_argument(
    "-C",
    "--common-core",
    type=bool,
    default=True,
    help="Calculate a common core among all input molecules before pairwise calculations",
)

# ------------------------------------------------------------------


# Main function
if "__main__" == __name__:
    startup()
