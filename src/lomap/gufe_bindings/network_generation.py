from __future__ import annotations

import itertools
import logging
from collections.abc import Callable

import networkx as nx
import numpy as np

try:
    from gufe import (
        AtomMapper,
        LigandAtomMapping,
        LigandNetwork,
        SmallMoleculeComponent,
    )
except ImportError:
    pass

from lomap._due import Doi, due
from lomap.graphgen import GraphGen
from lomap.utils import deprecated_kwargs, requires_package

logger = logging.getLogger(__name__)


@requires_package("gufe")
@deprecated_kwargs(name_mappings={"molecules": "ligands"})
@due.dcite(Doi("https://doi.org/10.1007/s10822-013-9678-y"), description="LOMAP")
def generate_lomap_network(
    ligands: list[SmallMoleculeComponent],
    mappers: AtomMapper | list[AtomMapper],
    scorer: Callable,
    distance_cutoff: float = 0.4,
    max_path_length: int = 6,
    actives: list[bool] | None = None,
    max_dist_from_active: int = 2,
    require_cycle_covering: bool = True,
    radial: bool = False,
    fast: bool = False,
    hub: SmallMoleculeComponent | None = None,
) -> LigandNetwork:
    """Generate a LigandNetwork according to Lomap's network creation rules

    Parameters
    ----------
    ligands : list[SmallMoleculeComponent]
       Molecules to include in the network.
    mappers : list[AtomMapper] or AtomMapper
       One or more Mapper functions to use to propose edge mappings.
    scorer : Callable
       Scoring function for edges. Should be a function which takes an
       AtomMapping and returns a value from 0.0 (worst) to 1.0 (best), inclusive.
       These values are used as the "distance" between two molecules,
       and compared against the 'distance_cutoff' parameter.
    distance_cutoff : float, default 0.4
       Edges with a score < 1 - distance_cutoff will be rejected.
    max_path_length : int, default 6
      Maximum edge distance between any two molecules in the resulting network.
    actives : list[bool] | None, default None
      If defined, a tag for each ligand which defines if it is an active molecule.
    max_dist_from_active : int, default 2
      When 'actives' is given, constrains the resulting map to be within this
      number of edges (e.g. distance) from an active molecule.
    require_cycle_covering : bool, default True
      If ``True``, attempt to ensure that every ligand has
      redundant paths to its neighbors, giving the network robustness against
      individual perturbation failures. This is achieved by rejecting edge
      removals that would leave a node outside a cycle or create a new bridge
      (an edge whose removal disconnects the graph). If ``False``, this constraint is
      relaxed and the resulting network may have no cycles.
    radial : bool, default False
      construct a radial (star) network. Note that the map will not necessarily
      be a true radial map; edges will still obey the ``distance_cutoff`` and if
      ``require_cycle_covering`` is ``True``, this radial map will still feature cycles.
    fast : bool, default False
      When both ``fast`` and ``radial`` are ``True``, switch the initial
      graph construction to only consider hub-spoke edges (every ligand
      connected to the hub/lead) rather than all pairwise edges. This
      makes network construction faster, at the potential cost of a
      less optimal network.
    hub : SmallMoleculeComponent | None, default None
      If radial is ``True``, force this ligand to be the center/hub of the radial graph.
    """
    if not mappers:
        raise ValueError("At least one Mapper must be provided")
    if isinstance(mappers, AtomMapper):
        mappers = [mappers]
    if actives is None:
        actives = [False] * len(ligands)

    # gen n x n mappings with scores
    # initially all zero scores, i.e. impossible
    mtx = np.zeros((len(ligands), len(ligands)), dtype=float)
    # np array of mappings
    mps = np.zeros_like(mtx, dtype=object)

    # for all combinations of ligands
    for i, j in itertools.combinations(range(len(ligands)), 2):
        mA, mB = ligands[i], ligands[j]

        # pick best score across all mappings from all mappings
        best_mp: LigandAtomMapping | None = None
        best_score = 0.0
        for mapper in mappers:
            mp: LigandAtomMapping
            score: float
            try:
                score, mp = max((scorer(mp), mp) for mp in (mapper.suggest_mappings(mA, mB)))
            except ValueError:
                # if mapper returned no mappings
                continue
            else:
                if score > best_score:
                    best_mp = mp
                    best_score = score

        if best_mp is None:
            logger.debug(f"Found no mapping for {mA} {mB}")
            continue

        logger.debug(f"Mapping for {mA} {mB} has score {best_score}")

        mtx[i, j] = mtx[j, i] = best_score
        mps[i, j] = mps[j, i] = best_mp.with_annotations({"score": best_score})

    gg = GraphGen(
        score_matrix=mtx,
        ids=list(range(mtx.shape[0])),
        names=[m.name for m in ligands],
        max_path_length=max_path_length,
        actives=actives,
        max_dist_from_active=max_dist_from_active,
        similarity_cutoff=1 - distance_cutoff,
        require_cycle_covering=require_cycle_covering,
        radial=radial,
        fast=fast,
        hub=hub.name if hub else None,
    )
    n: nx.Graph = gg.resultGraph

    ln = LigandNetwork(
        edges=[mps[i, j] for i, j in n.edges],
        nodes=ligands,
    )

    return ln
