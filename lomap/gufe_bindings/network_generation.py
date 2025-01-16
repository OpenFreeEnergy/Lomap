from gufe import(
    LigandAtomMapping,
    AtomMapper,
    LigandNetwork,
)
import gufe
import itertools
import functools
import logging
import networkx as nx
import numpy as np
from typing import Callable, Optional, Union, Any
import warnings

from ..graphgen import GraphGen
from .._due import due, Doi


logger = logging.getLogger(__name__)

def deprecated_kwarg(old_name:str, new_name:str):
  def decorator(func: Callable)->Callable:
      # old_name = 'molecules'
      # new_name = 'ligands'
      @functools.wraps(func)
      def wrapper(*args:Any, **kwargs:Any):
          if old_name in kwargs:
            kwargs[new_name] = kwargs.pop(old_name)
            warnings.warn(f"{func.__name__} argument '{old_name}' is deprecated. Please use '{new_name}' instead.",
                          DeprecationWarning)
            return func(*args, **kwargs)
      return wrapper
  return decorator

@deprecated_kwarg(old_name='molecules', new_name='ligands')
@due.dcite(Doi("https://doi.org/10.1007/s10822-013-9678-y"), description="LOMAP")
def generate_lomap_network(
        ligands: list[gufe.SmallMoleculeComponent],
        mappers: Union[AtomMapper, list[AtomMapper]],
        scorer: Callable,
        distance_cutoff: float = 0.4,
        max_path_length=6,
        actives: Optional[list[bool]] = None,
        max_dist_from_active=2,
        require_cycle_covering: bool = True,
        radial: bool = False,
        fast: bool = False,
        hub: Optional[gufe.SmallMoleculeComponent] = None,
    ) -> LigandNetwork:
    """Generate a LigandNetwork according to Lomap's network creation rules

    Parameters
    ----------
    ligands : list[SmallMoleculeComponent]
       molecules to map
    mappers : list[AtomMapper] or AtomMapper
       one or more Mapper functions to use to propose edges
    scorer: function
       scoring function for edges.  Should be a function which takes an AtomMapping and returns a value from 1.0 (best)
       to 0.0 (worst).  These values are use as the "distance" between two molecules, and compared against the
       'distance_cutoff' parameter
    distance_cutoff : float
       the maximum distance/dissimilarity between two molecules for an edge to be accepted
    max_path_length : int
      maximum distance between any two molecules in the resulting network
    actives : list[bool]
      for each molecule, if it is tagged as an active molecule
    max_dist_from_active
      when 'actives' is given, constrains the resulting map to
    require_cycle_covering : bool
      add cycles into the network
    radial : bool
      construct a radial/starmap network.  Note that this the map will not necessarily be a true radial map; edges
      will still obey the distance_cutoff and if 'require_cycle_covering' is true, this radial map will still
      feature cycles, default False
    fast : bool
      hmmm...
    hub : SmallMoleculeComponent, optional
      if radial, force this ligand to be the centre/hub of the radial graph
    """
    if not mappers:
        raise ValueError("At least one Mapper must be provided")
    if isinstance(mappers, gufe.AtomMapper):
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
        best_mp: Optional[LigandAtomMapping] = None
        best_score = 0.0
        for mapper in mappers:
            mp: LigandAtomMapping
            score: float
            try:
                score, mp = max((scorer(mp), mp)
                                for mp in (mapper.suggest_mappings(mA, mB)))
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
        mps[i, j] = mps[j, i] = best_mp.with_annotations({'score': best_score})

    gg = GraphGen(score_matrix=mtx,
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
