"""
The MCS class wrapped to provide a gufe interface

"""
import gufe
from gufe import AtomMapper, LigandAtomMapping
from typing import Iterable

from .. import mcs as lomap_mcs
from .._due import due, Doi


class LomapAtomMapper(AtomMapper):
    time: int
    threed: bool
    max3d: float
    element_change: bool
    seed: str
    shift: bool

    def __init__(self, *, time: int = 20, threed: bool = True,
                 max3d: float = 1.0, element_change: bool = True,
                 seed: str = '', shift: bool = False):
        """Wraps the MCS atom mapper from Lomap.

        Kwargs are passed directly to the MCS class from Lomap for each mapping
        created

        Parameters
        ----------
        time : int, optional
          timeout of MCS algorithm, passed to RDKit
          default 20
        threed : bool, optional
          if true, positional info is used to choose between symmetrically
          equivalent mappings and prune the mapping, default True
        max3d : float, optional
          maximum discrepancy in Angstroms between atoms before mapping is not
          allowed, default 1000.0, which effectively trims no atoms
        element_change: bool, optional
          whether to allow element changes in the mappings, default True
        seed: str, optional
          SMARTS string to use as seed for MCS searches.  When used across an
          entire set of ligands, this can speed up calculations considerably
        shift: bool, optional
          when determining 3D overlap, if to translate the two molecules MCS to minimise
          RMSD to boost potential alignment.
        """
        self.time = time
        self.threed = threed
        self.max3d = max3d
        self.element_change = element_change
        self.seed = seed
        self.shift = shift

    def __repr__(self):
        return (f"<LomapAtomMapper (time={self.time}, threed={self.threed}, "
                f"max3d={self.max3d}, element_change={self.element_change}, "
                f"seed='{self.seed}', shift={self.shift})>")

    @due.dcite(Doi("https://doi.org/10.1007/s10822-013-9678-y"), description="LOMAP")
    def suggest_mappings(self,
                         componentA: gufe.SmallMoleculeComponent,
                         componentB: gufe.SmallMoleculeComponent) -> Iterable[LigandAtomMapping]:
        """Generate one or more mappings between two small molecules
        
        Parameters
        ----------
        componentA, componentB: gufe.SmallMoleculeComponent

        Returns
        -------
        mapping : Iterable[LigandAtomMapping]
          potential mappings
        """
        try:
            mcs = lomap_mcs.MCS(componentA.to_rdkit(), componentB.to_rdkit(),
                                time=self.time,
                                threed=self.threed, max3d=self.max3d,
                                element_change=self.element_change,
                                seed=self.seed,
                                shift=self.shift)
        except ValueError:
            # if no match found, Lomap throws ValueError, so we just yield
            # generator with no contents
            return

        mapping_string = mcs.all_atom_match_list()
        # lomap spits out "1:1,2:2,...,x:y", so split around commas,
        # then colons and coerce to ints
        mapping_dict = dict((map(int, v.split(':'))
                             for v in mapping_string.split(',')))

        yield LigandAtomMapping(
            componentA=componentA,
            componentB=componentB,
            componentA_to_componentB=mapping_dict,
        )
        return
