import importlib.resources
import gufe
import pytest
from rdkit import Chem

import lomap


@pytest.fixture
def basic():
    data_path = importlib.resources.files('lomap.tests.basic')

    fns = [data_path / f for f in
           (
                '1,3,7-trimethylnaphthalene.mol2',
                '1-butyl-4-methylbenzene.mol2',
                '2,6-dimethylnaphthalene.mol2',
                '2-methyl-6-propylnaphthalene.mol2',
                '2-methylnaphthalene.mol2',
                '2-naftanol.mol2',
                'methylcyclohexane.mol2',
                'toluene.mol2',
            )]

    mols = [gufe.SmallMoleculeComponent(Chem.MolFromMol2File(str(f), removeHs=False))
            for f in fns]

    return mols


def test_generate_network_smoketest(basic):
    with pytest.deprecated_call():
        network = lomap.generate_lomap_network(
            molecules=basic,
            mappers=lomap.LomapAtomMapper(),
            scorer=lomap.default_lomap_score,
        )

        assert isinstance(network, gufe.LigandNetwork)

