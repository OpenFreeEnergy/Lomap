import importlib
import importlib.resources
import os

import pytest
from rdkit import Chem

try:
    import gufe

    HAS_GUFE = True
except ImportError:
    HAS_GUFE = False


LOMAP_DOI = "https://doi.org/10.1007/s10822-013-9678-y"


def _skip_if_duecredit_disabled():
    try:
        import duecredit  # noqa: F401
    except ImportError:
        pytest.skip("duecredit is not installed")

    if os.environ.get("DUECREDIT_ENABLE", "no").lower() in ("no", "0", "false"):
        pytest.skip("duecredit is disabled (set DUECREDIT_ENABLE=yes to enable)")


@pytest.fixture(autouse=True)
def require_duecredit():
    _skip_if_duecredit_disabled()


@pytest.fixture
def two_smcs():
    data_path = importlib.resources.files("lomap.tests.basic")
    mols = [
        gufe.SmallMoleculeComponent(Chem.MolFromMol2File(str(data_path / f), removeHs=False))
        for f in ("toluene.mol2", "2-methylnaphthalene.mol2")
    ]
    return mols


class TestDuecredit:
    def test_duecredit_active(self):
        from lomap._due import due

        assert due.active

    @pytest.mark.parametrize(
        "module, doi",
        [
            ("lomap.gufe_bindings.scorers", LOMAP_DOI),
        ],
    )
    def test_module_citations(self, module, doi):
        from lomap._due import due

        importlib.import_module(module)
        assert due.citations[(module, doi)].cites_module

    @pytest.mark.skipif(not HAS_GUFE, reason="requires gufe")
    def test_mapper_citation(self, two_smcs):
        from lomap import LomapAtomMapper
        from lomap._due import due

        mapper = LomapAtomMapper()
        list(mapper.suggest_mappings(*two_smcs))

        assert ("lomap.gufe_bindings.mapper:suggest_mappings", LOMAP_DOI) in due.citations

    @pytest.mark.skipif(not HAS_GUFE, reason="requires gufe")
    def test_network_generation_citation(self, two_smcs):
        from lomap import LomapAtomMapper, default_lomap_score, generate_lomap_network
        from lomap._due import due

        generate_lomap_network(
            ligands=two_smcs,
            mappers=LomapAtomMapper(),
            scorer=default_lomap_score,
        )

        assert (
            "lomap.gufe_bindings.network_generation:generate_lomap_network",
            LOMAP_DOI,
        ) in due.citations
