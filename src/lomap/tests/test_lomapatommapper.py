import pytest

from lomap import LomapAtomMapper

try:
    import gufe
    from gufe.tests.test_tokenization import GufeTokenizableTestsMixin

    HAS_GUFE = True
except ImportError:
    # Fake assign to object to avoid issues
    GufeTokenizableTestsMixin = object
    HAS_GUFE = False


@pytest.mark.skipif(HAS_GUFE, reason="requires not having gufe installed")
def test_lomap_atommaper_no_gufe_error():
    msg = "gufe is required to use `LomapAtomMapper` but is not installed."
    with pytest.raises(ImportError, match=msg):
        _ = LomapAtomMapper()


@pytest.mark.skipif(not HAS_GUFE, reason="requires gufe installed")
class TestLomapAtomMapper(GufeTokenizableTestsMixin):
    cls = LomapAtomMapper
    key = None
    repr = (
        "<LomapAtomMapper (time=20, threed=True, max3d=1.0, "
        "element_change=True, seed='None', shift=False)>"
    )

    @pytest.fixture
    def instance(self):
        return LomapAtomMapper()

    def test_to_dict_roundtrip_different(self):
        ref_vals = {
            "time": 19,
            "threed": False,
            "max3d": 999.0,
            "element_change": False,
            "seed": "CC",
            "shift": False,
        }

        m = LomapAtomMapper(**ref_vals)

        d = m.to_dict()

        m2 = LomapAtomMapper.from_dict(d)

        assert isinstance(m2, LomapAtomMapper)
        assert m2.time == ref_vals["time"]
        assert m2.threed == ref_vals["threed"]
        assert m2.max3d == ref_vals["max3d"]
        assert m2.element_change == ref_vals["element_change"]
        assert m2.seed == ref_vals["seed"]
        assert m2.shift == ref_vals["shift"]

    def test_repr_different(self):
        m = LomapAtomMapper(time=15, seed="c1ccccc1")

        assert repr(m) == (
            "<LomapAtomMapper (time=15, threed=True, max3d=1.0, "
            "element_change=True, seed='c1ccccc1', shift=False)>"
        )
