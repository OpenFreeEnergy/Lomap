# ******************
# MODULE DOCSTRING
# ******************

"""

LOMAP: Graph generation
=====

Alchemical free energy calculations hold increasing promise as an aid to drug
discovery efforts. However, applications of these techniques in discovery
projects have been relatively few, partly because of the difficulty of planning
and setting up calculations. The Lead Optimization Mapper (LOMAP) is an
automated algorithm to plan efficient relative free energy calculations between
potential molecules within a substantial set of molecules.

"""

# *****************************************************************************
# Lomap2: A toolkit to plan alchemical relative binding affinity calculations
# Copyright 2015 - 2016  UC Irvine and the Authors
#
# Authors: Dr Gaetano Calabro' and Dr David Mobley
#
# This part of the code has been originally made by Jonathan Redmann,
# and Christopher Summa at Summa Lab, Dept. of Computer Science,
# University of New Orleans and it has just been adapted to the new Lomap code
#
# *****************************************************************************
from __future__ import annotations

import copy
import logging
import os.path
import shutil
import subprocess
import tempfile
import traceback
import warnings
from operator import itemgetter
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from lomap.dbmol import DBMolecules

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

from lomap.utils import requires_package

__all__ = ["GraphGen"]

logger = logging.getLogger(__name__)


def find_non_cyclic_nodes(subgraph: nx.Graph) -> set[int]:
    """
    Generates a list of nodes of the subgraph that are not in a cycle

    Parameters
    ----------
    subgraph : nx.Graph
      The subgraph to check for non-cycle nodes.

    Returns
    -------
    missingNodesSet : set[int]
      The set of graph nodes that are not in a cycle.

    """
    cycleList = nx.cycle_basis(subgraph)

    cycleNodes = [node for cycle in cycleList for node in cycle]

    missingNodesSet = {node for node in subgraph.nodes() if node not in cycleNodes}

    return missingNodesSet


def find_non_cyclic_edges(subgraph: nx.Graph) -> set[tuple[int, int]]:
    """
    Generates a set of edges of the subgraph that are not in a cycle (called
    "bridges" in networkX terminology).

    Parameters
    ----------
    subgraph : nx.Graph
      The subgraph to check for non-cycle edges.

    Returns
    -------
    missingEdgesSet : set[tuple[int, int]]
      The set of edges that are not in a cycle.

    """
    missingEdgesSet = set(nx.bridges(subgraph))

    return missingEdgesSet


class GraphGen:
    """This class is used to set and generate the graph used to plan binding free energy calculation

    Notes
    -----
    To use the :class:`GraphGen` :meth:`draw`, :meth:`generate_depictions`,
    and :meth:`_write_output_dot_graph` methods, you will need to
    install the optional dependency ``pygraphviz``.
    """

    def __init__(
        self,
        score_matrix: np.ndarray,
        ids: list[int],
        names: list[str],
        max_path_length: int,
        actives: list[bool],
        max_dist_from_active: int,
        similarity_cutoff: float,
        require_cycle_covering: bool,
        radial: bool,
        fast: bool,
        hub: str | None = None,
    ):
        """

        Parameters
        ----------
        score_matrix : np.ndarray
          Array of scores between each molecule. Should be a symmetric (n x n) matrix.
        ids : list[int]
          Indices for each molecule. Should be the same length as ``score_matrix``.
          These ids are used as the ``'ID'`` attribute in the resulting graph.
        names : list[str]
          List of string identifiers for each molecule.
          These names are used as the ``'fname_comp'`` attribute in the resulting graph.
        max_path_length : int
          The maximum allowed shortest path length between any two molecules in the graph.
        actives : list[bool]
          For each molecule in input, whether it is considered active. Used in conjunction
          with the ``max_dist_from_active`` argument.
        max_dist_from_active : int
          The maximum allowed shortest path length from any molecule to the nearest active molecule.
        similarity_cutoff : float
          The value above which edges must be to be considered viable. ``0.0`` would allow all edges.
        require_cycle_covering : bool
          If ``True``, enforce that every node is part of a cycle in the graph.
        radial : bool
          Whether to construct a radial graph. Note that this radial graph will still include cycles.
        fast : bool
          If ``True``, use the faster radial-only algorithm (requires ``radial=True``).
        hub : str, optional
          The **name** of the molecule to use as the center of the hub.
        """
        self.score_matrix = score_matrix
        self.maxPathLength = max_path_length
        self.maxDistFromActive = max_dist_from_active
        self.similarityScoresLimit = similarity_cutoff

        if radial:
            self.lead_index = self.pick_lead(hub=hub, names=names, strict_mtx=score_matrix)
        else:
            self.lead_index = None

        # A set of nodes that will be used to save nodes that are not a cycle cover for a given subgraph
        self.nonCycleNodesSet: set[int] = set()

        # A set of edges that will be used to save edges that are acyclic for given subgraph
        self.nonCycleEdgesSet: set[tuple[int, int]] = set()

        # A count of the number of nodes that are not within self.maxDistFromActive edges
        # of an active
        self.distanceToActiveFailures = 0

        # The following Section has been strongly copied/adapted from the original implementation

        # Generate a list related to the disconnected graphs present in the initial graph
        fast_map = fast and radial
        self.initialSubgraphList = self.generate_initial_subgraph_list(
            fast_map=fast_map,
            strict_mtx=score_matrix,
            ids=ids,
            names=names,
            is_active=actives,
            lead_index=self.lead_index,
        )

        # A list of elements made of [edge, weights] for each subgraph
        self.subgraphScoresLists = self.generate_subgraph_scores_lists(self.initialSubgraphList)

        # Eliminates from each subgraph those edges whose weights are less than the hard limit
        self.remove_edges_below_hard_limit(
            subgraphlist=self.initialSubgraphList,
            scores=self.subgraphScoresLists,
            similarity_scores_limit=similarity_cutoff,
        )

        # Make a new master list of subgraphs now that there may be more disconnected components
        self.workingSubgraphsList = self.generate_working_subgraphs_list(self.initialSubgraphList)

        # Make a new sorted list of [edge, weights] for each subgraph now that there may be new subgraphs
        self.workingSubgraphScoresLists = self.generate_subgraph_scores_lists(
            self.workingSubgraphsList
        )

        # Remove edges, whose removal does not violate constraints, from the subgraphs,
        # starting with lowest similarity score first

        if fast and radial:
            # if we use the fast and radial option, just need to add the surrounding edges from the initial graph
            # assert type for lead_index to make linters happy
            # note: radial=True guarantees pick_lead sets this
            assert self.lead_index is not None
            self.resultGraph = self.add_surrounding_edges(
                subgraphs=self.workingSubgraphsList,
                score_matrix=score_matrix,
                lead_index=self.lead_index,
                similarity_score_limit=similarity_cutoff,
            )
            # after adding the surround edges, some subgraphs may merge into a larger graph and so need to update the
            # current subgraphs
            # self.resultingSubgraphsList = copy.deepcopy(self.workingSubgraphsList)
            # merge all Subgraphs together for layout
            # self.resultGraph = self.merge_all_subgraphs()
        else:
            # >>>>>>>>>>>>>>>>>>>>>>>>>>>ISSUE ORDER PROBLEM<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            self.minimize_edges(require_cycle_covering)
            # >>>>>>>>>>>>>>>>>>>>>>>>>>>ISSUE ORDER PROBLEM<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            # Collect together disjoint subgraphs of like charge into subgraphs
            self.resultingSubgraphsList = copy.deepcopy(self.workingSubgraphsList)

            # Combine separate subgraphs into a single resulting graph
            self.resultGraph = self.merge_all_subgraphs(self.workingSubgraphsList)

            # Make a copy of the resulting graph for later processing in connectResultingComponents()
            self.copyResultGraph = self.resultGraph.copy()

            # Holds list of edges that were added in the connect components phase
            self.edgesAddedInFirstTreePass: list[tuple[int, int, float]] = []

            # Add edges to the resultingGraph to connect its components
            self.connect_subgraphs()

    @staticmethod
    def pick_lead(hub: str | None, names: list[str], strict_mtx: np.ndarray) -> int | None:
        """Pick a lead molecule.

        Parameters
        ----------
        hub : str | None
          The name of the desired hub molecule, or ``None`` to auto-select.
        names : list[str]
          Names of each molecule.
        strict_mtx : np.ndarray
          Scoring matrix.

        Returns
        -------
        int | None
          The index of the lead molecule, or ``None`` if no suitable hub was found.
        """
        # TODO: remove support for "None" string for hub in next release
        if hub == "None":
            msg = (
                "Support for passing hub as str(None) instead of None directly "
                "to `GraphGen.pick_lead` will be removed in the next release of Lomap"
            )
            warnings.warn(msg, DeprecationWarning)
            hub = None

        if hub is not None:
            # hub radial option. Use the provided reference molecule as a hub
            hub_index = None
            for i, nm in enumerate(names):
                if os.path.basename(nm) == hub:
                    hub_index = i
            if hub_index is None:
                logging.info(
                    f"Warning: the specified center molecule {hub} is not in the "
                    "molecule database, will not use the radial option."
                )
            return hub_index
        else:
            # complete radial option.
            # Pick the molecule with the highest total similarity to all other molecules to use as a hub
            N = len(names)
            all_sum_i = []
            for i in range(N):
                sum_i = 0
                for j in range(N):
                    sum_i += strict_mtx[i, j]
                all_sum_i.append(sum_i)
            max_value = max(all_sum_i)
            max_index = [i for i, x in enumerate(all_sum_i) if x == max_value]
            max_index_final = max_index[0]
            return max_index_final

    @staticmethod
    def generate_initial_subgraph_list(
        fast_map: bool,
        strict_mtx: np.ndarray,
        ids: list[int],
        names: list[str],
        is_active: list[bool],
        lead_index: int | None,
    ) -> list[nx.Graph]:
        """
        Generate a starting graph with edges connecting all the
        molecules with a positive strict similarity score.

        Parameters
        ----------
        fast_map : bool
          If ``True``, use the faster radial-only algorithm.
        strict_mtx : np.ndarray
          Matrix of scores between molecules.
        ids : list[int]
          List of identifiers for each molecule.
        names : list[str]
          Names of each molecule.
        is_active : list[bool]
          For each molecule, whether it is active.
        lead_index : int | None
          The index of the lead molecule; cannot be ``None`` if ``fast_map`` is ``True``.

        Returns
        -------
        initialSubgraphList : list[nx.Graph]
          The list of connected component graphs.

        Raises
        ------
        ValueError
          If ``lead_index`` is not defined and ``fast_map`` is ``True``.
        """
        compound_graph: nx.Graph = nx.Graph()

        if fast_map and (lead_index is None):
            msg = "`lead_index` must be defined if using the fast map option"
            raise ValueError(msg)

        if not fast_map:
            # if not fast map option, connect all possible nodes to generate the initial graph
            for i in range(len(ids)):
                if i == 0:
                    compound_graph.add_node(
                        i, ID=ids[i], fname_comp=os.path.basename(names[i]), active=is_active[i]
                    )

                for j in range(i + 1, len(ids)):
                    if i == 0:
                        compound_graph.add_node(
                            j, ID=ids[j], fname_comp=os.path.basename(names[j]), active=is_active[j]
                        )

                    wgt = strict_mtx[i, j]

                    if wgt > 0.0:
                        compound_graph.add_edge(i, j, similarity=wgt, strict_flag=True)
        else:
            # if fast map option, then add all possible radial edges as the initial graph
            for i in range(len(ids)):
                # add the node for i
                compound_graph.add_node(i, ID=ids[i], fname_comp=os.path.basename(names[i]))
                if i != lead_index:
                    wgt = strict_mtx[i, lead_index]
                    if wgt > 0:
                        compound_graph.add_edge(i, lead_index, similarity=wgt, strict_flag=True)

        initialSubgraphGen = [
            compound_graph.subgraph(c).copy() for c in nx.connected_components(compound_graph)
        ]
        initialSubgraphList = [x for x in initialSubgraphGen]

        return initialSubgraphList

    @staticmethod
    def generate_subgraph_scores_lists(
        subgraphList: list[nx.Graph],
    ) -> list[list[tuple[int, int, float]]]:
        """
        This function generate a list of lists where each inner list is the
        weights of each edge in a given subgraph in the subgraphList,
        sorted from lowest to highest

        Parameters
        ----------
        subgraphList : list[nx.Graph]
          The list of subgraphs to generate scores for.

        Returns
        -------
        subgraphScoresLists : list[list[tuple[int, int, float]]]
          Each list contains a tuple with the graph node indices and their
          similarity as weight, sorted from lowest to highest.
        """

        subgraphScoresLists = []

        for subgraph in subgraphList:
            weightsDictionary = nx.get_edge_attributes(subgraph, "similarity")

            subgraphWeightsList = [
                (edge[0], edge[1], weightsDictionary[edge]) for edge in weightsDictionary.keys()
            ]

            subgraphWeightsList.sort(key=lambda entry: entry[2])

            subgraphScoresLists.append(subgraphWeightsList)

        return subgraphScoresLists

    @staticmethod
    def remove_edges_below_hard_limit(
        subgraphlist: list[nx.Graph],
        scores: list[list[tuple[int, int, float]]],
        similarity_scores_limit: float,
    ) -> None:
        """
        Remove edges below the set hard limit from each subGraph
        and from each weightsList.

        Operates on subgraphlist in-place!

        Parameters
        ----------
        subgraphlist : list[nx.Graph]
          The list of subgraphs to remove edges from.
        scores : list[list[tuple[int, int, float]]]
          The sorted edge scores for each subgraph.
        similarity_scores_limit : float
          The minimum similarity score threshold; edges below this are removed.
        """

        totalEdges = 0

        for subgraph in subgraphlist:
            weightsList = scores[subgraphlist.index(subgraph)]

            index = 0

            for edge in weightsList:
                if edge[2] < similarity_scores_limit:
                    subgraph.remove_edge(edge[0], edge[1])

                    index = weightsList.index(edge)

            del weightsList[: index + 1]

            totalEdges = totalEdges + subgraph.number_of_edges()

    @staticmethod
    def generate_working_subgraphs_list(subgraph_list: list[nx.Graph]) -> list[nx.Graph]:
        """
        Generate a list of connected subgraphs.

        This method is called after the deletion of the edges that have a
        weight less than the selected threshold. When this happens the
        subgraph maybe disconnected and so this method is called to
        create a new master list of connected subgraphs.

        Parameters
        ----------
        subgraph_list : list[nx.Graph]
          The list of subgraphs to regenerate connected components from.

        Returns
        -------
        workingSubgraphsList : list[nx.Graph]
          The new list of connected component subgraphs.
        """
        workingSubgraphsList = []

        for subgraph in subgraph_list:
            newSubgraphList = [
                subgraph.subgraph(c).copy() for c in nx.connected_components(subgraph)
            ]

            for newSubgraph in newSubgraphList:
                workingSubgraphsList.append(newSubgraph)

        return workingSubgraphsList

    def minimize_edges(self, require_cycle_covering: bool) -> None:
        """
        Minimize edges in each subgraph while ensuring constraints are met

        Parameters
        ----------
        require_cycle_covering : bool
          If ``True``, enforce that every node is part of a cycle in the graph.
        """

        for subgraph in self.workingSubgraphsList:
            weightsList = self.workingSubgraphScoresLists[self.workingSubgraphsList.index(subgraph)]

            # ISSUE ORDER IS ORIGINATED HERE
            # weightsList = sorted(weightsList, key = itemgetter(1))

            # This part has been copied from the original code
            self.nonCycleNodesSet = find_non_cyclic_nodes(subgraph)
            self.nonCycleEdgesSet = find_non_cyclic_edges(subgraph)
            numberOfComponents = nx.number_connected_components(subgraph)
            self.distanceToActiveFailures = self.count_distance_to_active_failures(
                subgraph, self.maxDistFromActive
            )

            if len(subgraph.edges()) > 2:  # Graphs must have at least 3 edges to be minimized
                for edge in weightsList:
                    if self.lead_index is not None:
                        # Here the radial option is applied, will check if the remove_edge is connect to
                        # the hub(lead) molecule, if the edge is connected to the lead molecule,
                        # then add it back into the graph.
                        if self.lead_index not in [edge[0], edge[1]]:
                            subgraph.remove_edge(edge[0], edge[1])
                            if not self.check_constraints(
                                subgraph, numberOfComponents, require_cycle_covering
                            ):
                                subgraph.add_edge(
                                    edge[0], edge[1], similarity=edge[2], strict_flag=True
                                )
                    elif edge[2] < 1.0:  # Don't remove edges with similarity 1
                        logging.info(
                            f"Trying to remove edge {edge[0]}-{edge[1]} with similarity {edge[2]}"
                        )
                        subgraph.remove_edge(edge[0], edge[1])
                        if not self.check_constraints(
                            subgraph, numberOfComponents, require_cycle_covering
                        ):
                            subgraph.add_edge(
                                edge[0], edge[1], similarity=edge[2], strict_flag=True
                            )
                        else:
                            logging.info(f"Removed edge {edge[0]}-{edge[1]}")
                    else:
                        logging.info(f"Skipping edge {edge[0]}-{edge[1]} as it has similarity 1")

    def add_surrounding_edges(
        self,
        subgraphs: list[nx.Graph],
        score_matrix: np.ndarray,
        lead_index: int,
        similarity_score_limit: float,
    ) -> nx.Graph:
        """
        Add surrounding edges in each subgraph to make sure all nodes are in cycle.

        Parameters
        ----------
        subgraphs : list[nx.Graph]
          The list of subgraphs to add surrounding edges to.
        score_matrix : np.ndarray
          Matrix of similarity scores between molecules.
        lead_index : int
          The index of the lead (hub) molecule.
        similarity_score_limit : float
          The minimum similarity score threshold for adding edges.

        Returns
        -------
        nx.Graph
          The subgraph containing the lead molecule, with surrounding edges added.

        Raises
        ------
        ValueError
          If the ``lead_index`` is not in any subgraph.
        """
        for subgraph in subgraphs:
            subgraph_nodes = subgraph.nodes()
            if lead_index in subgraph_nodes:
                # here we only consider the subgraph with lead molecule
                self.nonCycleNodesSet = find_non_cyclic_nodes(subgraph)
                self.nonCycleEdgesSet = find_non_cyclic_edges(subgraph)
                for node in self.nonCycleNodesSet:
                    # for each node in the noncyclenodeset:
                    # find the similarity compare to all other surrounding nodes
                    # and pick the one with the max score and connect them
                    node_score_list = []
                    for i in range(score_matrix.shape[0]):
                        if i != node and i != lead_index:
                            node_score_list.append(score_matrix[node, i])
                        else:
                            node_score_list.append(0.0)
                    max_value = max(node_score_list)
                    if max_value > similarity_score_limit:
                        max_index = [i for i, x in enumerate(node_score_list) if x == max_value]
                        max_index_final = max_index[0]
                        subgraph.add_edge(
                            node,
                            max_index_final,
                            similarity=score_matrix[node, max_index_final],
                            strict_flag=True,
                        )
                return subgraph
        raise ValueError(f"lead_index {lead_index} not found in any subgraph")

    def check_constraints(
        self, subgraph: nx.Graph, numComp: int, require_cycle_covering: bool
    ) -> bool:
        """
        Determine if the given subgraph still meets the constraints.

        Parameters
        ----------
        subgraph : nx.Graph
          The subgraph to check for the constraints.
        numComp : int
          The number of connected components.
        require_cycle_covering : bool
          If ``True``, enforce cycle covering.

        Returns
        -------
        constraintsMet : bool
          ``True`` if all the constraints are met, ``False`` otherwise.
        """

        constraintsMet = True

        if not self.remains_connected(subgraph, numComp):
            constraintsMet = False

        # The requirement to keep a cycle covering is now optional
        if constraintsMet and require_cycle_covering:
            if not self.check_cycle_covering(subgraph, self.nonCycleEdgesSet):
                constraintsMet = False

        if constraintsMet:
            if not self.check_max_distance(subgraph, max_path_length=self.maxPathLength):
                constraintsMet = False

        if constraintsMet:
            if not self.check_distance_to_active(
                subgraph, self.distanceToActiveFailures, self.maxDistFromActive
            ):
                constraintsMet = False

        return constraintsMet

    @staticmethod
    def remains_connected(subgraph: nx.Graph, numComponents: int) -> bool:
        """
        Determine if the subgraph remains connected after an edge has been
        removed.

        Parameters
        ----------
        subgraph : nx.Graph
          The subgraph to check for connectivity after the edge deletion.
        numComponents : int
          The number of connected components.

        Returns
        -------
        isConnected : bool
          ``True`` if the subgraph is connected, ``False`` otherwise.
        """
        is_connected = numComponents == nx.number_connected_components(subgraph)
        if not is_connected:
            logging.info("Rejecting edge deletion on graph connectivity")

        return is_connected

    @staticmethod
    def check_cycle_covering(subgraph: nx.Graph, non_cycle_edges_set: set[tuple[int, int]]) -> bool:
        """
        Checks if the subgraph has a cycle covering.

        Parameters
        ----------
        subgraph : nx.Graph
          The subgraph to check for a cycle covering after edge deletion.
        non_cycle_edges_set : set[tuple[int, int]]
          The set of non-cyclic edges before the attempted edge removal.

        Returns
        -------
        hasCovering : bool
          ``True`` if the subgraph has a cycle covering, ``False`` otherwise.

        Notes
        -----
        This has been extended from the original algorithm. We not only care
        if the number of acyclic nodes has increased, but also if the number
        of acyclic edges (bridges) has increased.
        If the number of acycclicc edges hasn't increased, then the number of
        acyclic nodes hasn't eitherrr, so that test is included in the edges test.
        """
        hasCovering = True

        # Have we increased the number of non-cyclic edges?
        if find_non_cyclic_edges(subgraph).difference(non_cycle_edges_set):
            hasCovering = False
            logging.info("Rejecting edge deletion on cycle covering")

        return hasCovering

    @staticmethod
    def check_max_distance(subgraph: nx.Graph, max_path_length: int) -> bool:
        """
        Check to see if the graph has paths from all molecules to all other
        molecules within the specified limit.

        Parameters
        ----------
        subgraph : nx.Graph
          The subgraph to check for the max distance between nodes.
        max_path_length : int
          The maximum allowed shortest path length between any two nodes.

        Returns
        -------
        withinMaxDistance : bool
          ``True`` if all nodes are within the specified max distance of each other.
        """
        withinMaxDistance = True

        for node in subgraph:
            eccentricity = nx.eccentricity(subgraph, node)
            # eccentricity is int but mypy thinks int | dict[Any, int]
            if eccentricity > max_path_length:  # type: ignore[operator]
                withinMaxDistance = False
                logging.info(f"Rejecting edge deletion on graph diameter for node {node}")

        return withinMaxDistance

    @staticmethod
    def count_distance_to_active_failures(subgraph: nx.Graph, max_dist_from_active: int) -> int:
        """
        Count the number of molecules that don't have a minimum-length path to an active
        within the specified limit.

        Parameters
        ----------
        subgraph : nx.Graph
          The subgraph to check for distance to active nodes.
        max_dist_from_active : int
          The maximum allowed shortest path length from any molecule to the nearest active.

        Returns
        -------
        failures : int
          Number of nodes that are not within ``max_dist_from_active`` to any active node.
        """

        failures = 0

        hasActives = False
        for node in subgraph.nodes():
            if subgraph.nodes[node]["active"]:
                hasActives = True
        if not hasActives:
            return 0  # No actives, so don't bother checking

        # Backwards compatible networkx support
        paths = dict(nx.shortest_path(subgraph))
        for node in subgraph.nodes():
            if not subgraph.nodes[node]["active"]:
                ok = False
                for node2 in subgraph.nodes():
                    if subgraph.nodes[node2]["active"]:
                        pathlen = len(paths[node][node2]) - 1  # No. edges is 1 less than no. nodes
                        if pathlen <= max_dist_from_active:
                            ok = True
                if not ok:
                    failures = failures + 1

        return failures

    def check_distance_to_active(
        self, subgraph: nx.Graph, distance_to_active_failures: int, max_distance_from_active: int
    ) -> bool:
        """
        Check to see if we have increased the number of distance-to-active failures.

        Parameters
        ----------
        subgraph : nx.Graph
          The subgraph to check for the distance to active nodes.
        distance_to_active_failures : int
          The current count of nodes that fail the distance-to-active check.
        max_distance_from_active : int
          The maximum allowed shortest path length from any molecule to the nearest active.

        Returns
        -------
        ok : bool
          ``True`` if the number of failed nodes has not increased.
        """
        count = self.count_distance_to_active_failures(subgraph, max_distance_from_active)
        failed = count > distance_to_active_failures
        if failed:
            logging.info(
                f"Rejecting edge deletion on distance-to-actives {count} vs {distance_to_active_failures}"
            )
        logging.info(
            f"Checking edge deletion on distance-to-actives {count} vs {distance_to_active_failures}"
        )
        return not failed

    @staticmethod
    def merge_all_subgraphs(working_subgraphs: list[nx.Graph]) -> nx.Graph:
        """Generates a single networkx graph object from the subgraphs that have
        been processed

        Parameters
        ----------
        working_subgraphs : list[nx.Graph]
          The list of subgraphs to merge.

        Returns
        -------
        finalGraph : nx.Graph
          The final graph produced by merging all the subgraphs. The produced
          graph may have disconnected parts.
        """

        finalGraph: nx.Graph = nx.Graph()

        for subgraph in working_subgraphs:
            finalGraph = nx.union(finalGraph, subgraph)

        return finalGraph

    def connect_subgraphs(self) -> None:
        """
        Adds edges to the resultGraph to connect as many components of the final
        graph possible.
        """

        connectSuccess = self.connect_graph_components_brute_force()

        while connectSuccess:
            connectSuccess = self.connect_graph_components_brute_force()

        # WARNING: The self.workingSubgraphsList at this point is different from
        # the copy self.resultingSubgraphsList made before

        connectSuccess = self.connect_graph_components_brute_force_2()

        while connectSuccess:
            connectSuccess = self.connect_graph_components_brute_force_2()

    def connect_graph_components_brute_force(self) -> bool:
        """
        Adds edges to the resultGraph to connect all components that can be
        connected, only one edge is added per component, to form a tree like
        structure between the different components of the resultGraph.

        Returns
        -------
        bool
          ``True`` if there are edges that could have been added.
        """

        generator_graph = [
            self.resultGraph.subgraph(c).copy() for c in nx.connected_components(self.resultGraph)
        ]

        self.workingSubgraphsList = [x for x in generator_graph]

        if len(self.workingSubgraphsList) == 1:
            return False

        edgesToCheck = []
        numzeros = 0

        for i in range(0, len(self.workingSubgraphsList)):
            nodesOfI = self.workingSubgraphsList[i].nodes()

            for j in range(i + 1, len(self.workingSubgraphsList)):
                nodesOfJ = self.workingSubgraphsList[j].nodes()

                # change the following lines to be compatible with networkx 2.0
                for k in nodesOfI.keys():
                    for l in nodesOfJ.keys():  # noqa: E741
                        # produce an edge from nodesOfI[k] and nodesofJ[l] if nonzero weights push
                        # this edge into possibleEdgeList """

                        # print 'Molecules (%d,%d)' % (nodesOfI[k],nodesOfJ[l])
                        # I assumed that the score matrix is symmetric. In the Graph part this
                        # does not seems to be true:

                        similarity = self.score_matrix[nodesOfI[k]["ID"], nodesOfJ[l]["ID"]]

                        if similarity > 0.0:
                            edgesToCheck.append((nodesOfI[k]["ID"], nodesOfJ[l]["ID"], similarity))
                        else:
                            numzeros = numzeros + 1

        if len(edgesToCheck) > 0:
            sortedList = sorted(edgesToCheck, key=itemgetter(2), reverse=True)
            edgeToAdd = sortedList[0]
            self.edgesAddedInFirstTreePass.append(edgeToAdd)
            self.resultGraph.add_edge(
                edgeToAdd[0], edgeToAdd[1], similarity=edgeToAdd[2], strict_flag=False
            )

            generator_graph = [
                self.resultGraph.subgraph(c).copy()
                for c in nx.connected_components(self.resultGraph)
            ]
            self.workingSubgraphsList = [x for x in generator_graph]

            return True

        else:
            return False

    def connect_graph_components_brute_force_2(self) -> bool:
        """
        Adds a second edge between each of the (former) components of the
        resultGraph to try to provide cycles between (former) components.

        Returns
        -------
        bool
          ``True`` if the addition of edges was possible.
        """

        if len(self.resultingSubgraphsList) == 1:
            return False

        edgesToCheck = []

        for idx_i in range(0, len(self.resultingSubgraphsList)):
            nodesOfI = self.resultingSubgraphsList[idx_i].nodes()

            for idx_j in range(idx_i + 1, len(self.resultingSubgraphsList)):
                nodesOfJ = self.resultingSubgraphsList[idx_j].nodes()

                for node_i in nodesOfI.keys():
                    for node_j in nodesOfJ.keys():
                        # produce an edge from nodesOfI[node_i] and nodesofJ[node_j] if
                        # nonzero weights push this edge into possibleEdgeList """

                        # print 'Molecules (%d,%d)' % (nodesOfI[node_i],nodesOfJ[node_j])
                        # I assumed that the score matrix is symmetric. In the Graph part
                        # this does not seems to be true: <<<<<<<<<<<<<DEBUG>>>>>>>>>>>>>>>
                        similarity = self.score_matrix[
                            nodesOfI[node_i]["ID"], nodesOfJ[node_j]["ID"]
                        ]

                        if similarity > 0.0:
                            edgesToCheck.append(
                                (nodesOfI[node_i]["ID"], nodesOfJ[node_j]["ID"], similarity)
                            )

        finalEdgesToCheck = [
            edge for edge in edgesToCheck if edge not in self.edgesAddedInFirstTreePass
        ]

        if len(finalEdgesToCheck) > 0:
            sortedList = sorted(finalEdgesToCheck, key=itemgetter(2), reverse=True)
            edgeToAdd = sortedList[0]

            self.resultGraph.add_edge(
                edgeToAdd[0], edgeToAdd[1], similarity=edgeToAdd[2], strict_flag=False
            )
            self.copyResultGraph.add_edge(
                edgeToAdd[0], edgeToAdd[1], similarity=edgeToAdd[2], strict_flag=False
            )

            generator_graph = [
                self.copyResultGraph.subgraph(c).copy()
                for c in nx.connected_components(self.copyResultGraph)
            ]
            self.resultingSubgraphsList = [x for x in generator_graph]

            return True

        else:
            return False

    @requires_package("pygraphviz")
    def generate_depictions(
        self,
        dbase: DBMolecules,
        max_images: int = 2000,
        max_mol_size: float = 50.0,
        edge_labels: bool = True,
    ) -> None:
        """
        Generate images of the network.

        Parameters
        ----------
        dbase : DBMolecules
          The molecule database object.
        max_images : int
          Max number of displayed chemical molecule images as graph nodes.
        max_mol_size : float
          The maximum threshold distance in angstroms used to select if a molecule is depicted.
        edge_labels : bool
          If ``True``, add labels on edges.
        """

        def max_dist_mol(mol):

            max_dist = 0.0
            conf = mol.GetConformer()

            for i in range(0, conf.GetNumAtoms()):
                crdi = np.array(
                    [
                        conf.GetAtomPosition(i).x,
                        conf.GetAtomPosition(i).y,
                        conf.GetAtomPosition(i).z,
                    ]
                )

                for j in range(i + 1, conf.GetNumAtoms()):
                    crdj = np.array(
                        [
                            conf.GetAtomPosition(j).x,
                            conf.GetAtomPosition(i).y,
                            conf.GetAtomPosition(j).z,
                        ]
                    )
                    dist = np.linalg.norm(crdi - crdj)

                    if dist > max_dist:
                        max_dist = dist

            return max_dist

        directory_name = tempfile.mkdtemp()

        temp_graph = self.resultGraph.copy()

        if nx.number_of_nodes(temp_graph) <= max_images:
            # Draw.DrawingOptions.atomLabelFontSize=30
            # Draw.DrawingOptions.dotsPerAngstrom=100

            for n in temp_graph:
                id_mol = temp_graph.nodes[n]["ID"]
                mol = dbase[id_mol].getMolecule()
                max_dist = max_dist_mol(mol)

                if max_dist < max_mol_size:
                    fname = os.path.join(directory_name, dbase[id_mol].getName() + ".png")
                    # 1, modify here to calculate the 2D structure for molecules cannot remove Hydrogens by rdkit
                    # 2, change the graph size to get better resolution
                    try:
                        mol = Chem.RemoveHs(mol)
                    except (AllChem.KekulizeException, ValueError):
                        # Note: newer versions of RDKit now use KekulizeException
                        # for backwards compatibility, we also include ValueError which used to be thrown.
                        ###### need to ask RDKit to fix this if possible, see the code
                        # issue tracker for more details######
                        logging.info(
                            f"Error attempting to remove hydrogens for molecule {dbase[id_mol].getName()} using RDKit. RDKit cannot kekulize the molecule"
                        )
                    Chem.rdDepictor.Compute2DCoords(mol)
                    from rdkit.Chem.Draw.MolDrawing import DrawingOptions

                    DrawingOptions.bondLineWidth = 2.5
                    Draw.MolToFile(
                        mol, fname, size=(200, 200), kekulize=False, fitimage=True, imageType="png"
                    )
                    temp_graph.nodes[n]["image"] = fname
                    # self.resultGraph.nodes[n]['label'] = ''
                    temp_graph.nodes[n]["labelloc"] = "t"
                    temp_graph.nodes[n]["penwidth"] = 2.5
                    # self.resultGraph.node[n]['xlabel'] =  self.resultGraph.nodes[n]['ID']
        for u, v, d in temp_graph.edges(data=True):
            if d["strict_flag"]:
                temp_graph[u][v]["color"] = "blue"
                temp_graph[u][v]["penwidth"] = 2.5
            else:
                temp_graph[u][v]["color"] = "red"
                temp_graph[u][v]["penwidth"] = 2.5
            if edge_labels:
                temp_graph[u][v]["label"] = round(d["similarity"], 2)

        nx.nx_agraph.write_dot(temp_graph, dbase.options["name"] + "_tmp.dot")

        cmd = "dot -Tpng " + dbase.options["name"] + "_tmp.dot -o " + dbase.options["name"] + ".png"

        os.system(cmd)
        cmd = "dot -Teps " + dbase.options["name"] + "_tmp.dot -o " + dbase.options["name"] + ".eps"

        os.system(cmd)
        cmd = "dot -Tpdf " + dbase.options["name"] + "_tmp.dot -o " + dbase.options["name"] + ".pdf"

        os.system(cmd)
        os.remove(dbase.options["name"] + "_tmp.dot")
        shutil.rmtree(directory_name, ignore_errors=True)

    def layout_info(self, dbase: DBMolecules) -> None:
        """
        Write layout information (connectivity) and scoring data to text files.

        Parameters
        ----------
        dbase : DBMolecules
          The molecule database object.
        """
        # pass the lead molecule index if the radial option is on and generate the
        # morph type of output required by FESetup
        if self.lead_index is not None:
            morph_txt = open(dbase.options["name"] + "_morph.txt", "w")
            morph_data = "morph_pairs = "
        with open(dbase.options["name"] + "_score_with_connection.txt", "w") as info_txt:
            all_key_id = dbase.dic_mapping.keys()
            data = [
                f"{'Index_1':<10},{'Index_2':<10},{'Filename_1':<25},{'Filename_2':<25},{'Str_sim':<15},{'Eff_sim':<15},{'Loose_sim':<15},{'Connect':<10}\n"
            ]
            for i in range(len(all_key_id) - 1):
                for j in range(i + 1, len(all_key_id)):
                    morph_string = None
                    connected = False
                    similarity = 0
                    try:
                        # Note: this portion of the code seems to be solely relying on the use of
                        # an IndexError not being raised at some point in the next two lines..
                        # it's less than ideal but we also won't remove it.
                        edgedata = [
                            d
                            for (u, v, d) in self.resultGraph.edges(data=True)
                            if ((u == i and v == j) or (u == j and v == i))
                        ]
                        similarity = edgedata[0]["similarity"]  # noqa: F841
                        connected = True
                    except IndexError:
                        pass
                    Filename_i = dbase.dic_mapping[i]
                    Filename_j = dbase.dic_mapping[j]
                    MCmap = dbase.get_MCSmap(i, j)
                    mapString = ""
                    if MCmap is not None:
                        mapString = MCmap
                    # print "Check the filename", Filename_i, Filename_j
                    strict_similarity = dbase.strict_mtx[i, j]
                    loose_similarity = dbase.loose_mtx[i, j]
                    true_strict_similarity = dbase.true_strict_mtx[i, j]
                    if connected:
                        new_line = f"{str(i):<10},{str(j):<10},{Filename_i:<25},{Filename_j:<25},{true_strict_similarity:<15.5f},{strict_similarity:<15.5f},{loose_similarity:<15.5f},{'Yes':<10},{mapString}\n"
                        # generate the morph type, and pick the start molecule based on the similarity
                        if self.lead_index is not None:
                            morph_i = Filename_i.split(".")[0]
                            morph_j = Filename_j.split(".")[0]
                            if i == self.lead_index:
                                morph_string = f"{morph_i} > {morph_j}, "
                            elif j == self.lead_index:
                                morph_string = f"{morph_j} > {morph_i}, "
                            else:
                                # compare i and j with the lead molecule, and
                                # pick the one with the higher similarity as the start molecule
                                similarity_i = dbase.strict_mtx[self.lead_index, i]
                                similarity_j = dbase.strict_mtx[self.lead_index, j]
                                if similarity_i > similarity_j:
                                    morph_string = f"{morph_i} > {morph_j}, "
                                else:
                                    morph_string = f"{morph_j} > {morph_i}, "
                            morph_data += morph_string
                    else:
                        new_line = f"{i:<10},{j:<10},{Filename_i:<25},{Filename_j:<25},{true_strict_similarity:<15.5f},{strict_similarity:<15.5f},{loose_similarity:<15.5f},{'No':<10},{mapString}\n"
                    data.append(new_line)
            info_txt.writelines(data)
            if self.lead_index is not None:
                morph_txt.write(morph_data)

    @requires_package("pygraphviz")
    def _write_output_dot_graph(self, filename: str) -> None:
        """
        Helper method to write the graph to GraphViz dot format.

        Parameters
        ----------
        filename : str
          Name of the file to write the dot file to.
        """
        nx.nx_agraph.write_dot(self.resultGraph, filename)

    def write_graph(
        self, dbase: DBMolecules, output_no_images: bool, output_no_graph: bool
    ) -> None:
        """
        Write the final generated NetworkX graph as ``.dot`` and the ``.ps`` files.
        The mapping between molecule IDs and molecules name is saved as text file.

        Parameters
        ----------
        dbase : DBMolecules
          The molecule database object.
        output_no_images : bool
          If ``True``, do not generate image files.
        output_no_graph : bool
          If ``True``, do not generate the dot graph file.

        Raises
        ------
        OSError
          If database or image writing fails.
        """

        try:
            dbase.write_dic()
            self.layout_info(dbase)
        except Exception as e:
            traceback.print_exc()
            raise OSError(f"{str(e)}: {dbase.options['name']}.txt")

        try:
            if not output_no_images:
                self.generate_depictions(dbase)
            if not output_no_graph:
                self._write_output_dot_graph(dbase.options["name"] + ".dot")
        except Exception as e:
            traceback.print_exc()
            raise OSError(f"Problems during the file generation: {str(e)}")

        logging.info(30 * "-")

        log = "The following files have been generated:"
        if not output_no_graph:
            log += f"\n{dbase.options['name']}.dot\tGraph file"
        if not output_no_images:
            log += f"\n{dbase.options['name']}.png\tPng file"
        log += f"\n{dbase.options['name']}.txt\tMapping Text file"
        logging.info(log)

        logging.info(30 * "-")

    @requires_package("pygraphviz")
    def draw(
        self,
        dbase: DBMolecules,
        max_images: int = 2000,
        max_nodes: int = 100,
        edge_labels: bool = True,
    ) -> None:
        """Plot the NetworkX graph using Matplotlib.

        Parameters
        ----------
        dbase : DBMolecules
          The molecule database object.
        max_images : int
          Max number of displayed chemical molecule images as graph nodes.
        max_nodes : int
          Max number of displayed nodes in the graph.
        edge_labels : bool
          If ``True``, add edge labels in the graph.

        Notes
        -----
        This requires the optional dependency ``pygraphviz``.
        """

        logging.info("\nDrawing....")

        if nx.number_of_nodes(self.resultGraph) > max_nodes:
            logging.info(
                f"The number of generated graph nodes {nx.number_of_nodes(self.resultGraph)} exceed the max number of drawable nodes {max_nodes}"
            )
            return

        def max_dist_mol(mol):

            max_dist = 0.0
            conf = mol.GetConformer()

            for i in range(0, conf.GetNumAtoms()):
                crdi = np.array(
                    [
                        conf.GetAtomPosition(i).x,
                        conf.GetAtomPosition(i).y,
                        conf.GetAtomPosition(i).z,
                    ]
                )

                for j in range(i + 1, conf.GetNumAtoms()):
                    crdj = np.array(
                        [
                            conf.GetAtomPosition(j).x,
                            conf.GetAtomPosition(i).y,
                            conf.GetAtomPosition(j).z,
                        ]
                    )
                    dist = np.linalg.norm(crdi - crdj)

                    if dist > max_dist:
                        max_dist = dist

            return max_dist

        # Determine the screen resolution by using xdpyinfo and removing massive qt dependency
        command = "xdpyinfo | grep dimensions"
        p = subprocess.run(command, stdout=subprocess.PIPE, shell=True, executable="/bin/bash")
        width = int(p.stdout.split()[1].split(b"x")[0])
        height = int(p.stdout.split()[1].split(b"x")[1])

        # Canvas scale factor
        scale_canvas = 0.75

        # Canvas resolution
        max_canvas_size = (int(width * scale_canvas), int(height * scale_canvas))

        fig = plt.figure(1, facecolor="white")

        fig.set_dpi(100)

        fig.set_size_inches(
            max_canvas_size[0] / fig.get_dpi(), max_canvas_size[1] / fig.get_dpi(), forward=True
        )

        ax = plt.subplot(111)
        plt.axis("off")

        pos = nx.nx_agraph.graphviz_layout(self.resultGraph, prog="neato")

        strict_edges = [
            (u, v) for (u, v, d) in self.resultGraph.edges(data=True) if d["strict_flag"]
        ]
        loose_edges = [
            (u, v) for (u, v, d) in self.resultGraph.edges(data=True) if not d["strict_flag"]
        ]

        node_labels = {u: d["ID"] for u, d in self.resultGraph.nodes(data=True)}

        # Draw nodes
        nx.draw_networkx_nodes(self.resultGraph, pos, node_size=500, node_color="r")
        # Draw node labels
        nx.draw_networkx_labels(self.resultGraph, pos, labels=node_labels, font_size=10)

        if edge_labels:
            # fmt: off
            edge_weight_strict = {
                    (u, v): d["similarity"]
                    for u, v, d in self.resultGraph.edges(data=True)
                    if d["strict_flag"]
            }
            edge_weight_loose = {
                    (u, v): d["similarity"]
                    for u, v, d in self.resultGraph.edges(data=True)
                    if not d["strict_flag"]
            }
            # fmt: on

            for key in edge_weight_strict:
                edge_weight_strict[key] = round(edge_weight_strict[key], 2)

            for key in edge_weight_loose:
                edge_weight_loose[key] = round(edge_weight_loose[key], 2)

            # edge strict
            nx.draw_networkx_edge_labels(
                self.resultGraph, pos, edge_labels=edge_weight_strict, font_color="g"
            )
            # edge loose
            nx.draw_networkx_edge_labels(
                self.resultGraph, pos, edge_labels=edge_weight_loose, font_color="r"
            )

        # edges strict
        nx.draw_networkx_edges(self.resultGraph, pos, edgelist=strict_edges, edge_color="g")
        # edges loose
        nx.draw_networkx_edges(self.resultGraph, pos, edgelist=loose_edges, edge_color="r")

        if nx.number_of_nodes(self.resultGraph) <= max_images:
            trans = ax.transData.transform
            trans2 = fig.transFigure.inverted().transform

            cut = 1.0

            frame = 10
            xmax = cut * max(xx for xx, yy in pos.values()) + frame
            ymax = cut * max(yy for xx, yy in pos.values()) + frame

            xmin = cut * min(xx for xx, yy in pos.values()) - frame
            ymin = cut * min(yy for xx, yy in pos.values()) - frame

            plt.xlim(xmin, xmax)
            plt.ylim(ymin, ymax)

            h = 20
            w = 20

            mol_size = (200, 200)

            for each_node in self.resultGraph:
                id_mol = self.resultGraph.nodes[each_node]["ID"]
                # skip remove Hs by rdkit if Hs cannot be removed
                try:
                    mol = Chem.RemoveHs(dbase[id_mol].getMolecule())
                except (AllChem.KekulizeException, ValueError):
                    # Note: newer versions of RDKit now use KekulizeException
                    # for backwards compatibility, we also include ValueError which used to be thrown.
                    ###### need to ask RDKit to fix this if possible, see the code
                    # issue tracker for more details######
                    mol = dbase[id_mol].getMolecule()
                    logging.info(
                        f"Error attempting to remove hydrogens for molecule {dbase[id_mol].getName()} using RDKit. RDKit cannot kekulize the molecule"
                    )

                # max_dist = max_dist_mol(mol)
                # if max_dist > 7.0:
                #     continue

                Chem.rdDepictor.Compute2DCoords(mol)
                # add try exception for cases cannot be draw
                try:
                    img_mol = Draw.MolToImage(mol, mol_size, kekulize=False)
                except Exception:
                    img_mol = None
                    logging.exception(
                        "This mol cannot be drawn using the RDKit Draw function, need to check for more details..."
                    )

                xx, yy = trans(pos[each_node])
                xa, ya = trans2((xx, yy))

                nodesize_1 = 300.0 / (h * 100)
                nodesize_2 = 300.0 / (w * 100)

                p2_2 = nodesize_2 / 2
                p2_1 = nodesize_1 / 2

                a = plt.axes([xa - p2_2, ya - p2_1, nodesize_2, nodesize_1])  # type: ignore[arg-type]
                # self.resultGraph.nodes[id_mol]['image'] = img_mol
                # a.imshow(self.resultGraph.node[each_node]['image'])
                a.imshow(img_mol)
                a.axis("off")

        # plt.savefig('graph.png', facecolor=fig.get_facecolor())
        # print 'Graph .png file has been generated...'

        plt.show()
