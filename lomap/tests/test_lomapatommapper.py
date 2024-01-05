from lomap import LomapAtomMapper


def test_repr():
    m = LomapAtomMapper()

    assert repr(m) == ("<LomapAtomMapper (time=20, threed=True, max3d=1000.0, "
                       "element_change=True, seed='', shift=True)>")

    m = LomapAtomMapper(time=15, seed='c1ccccc1')

    assert repr(m) == ("<LomapAtomMapper (time=15, threed=True, max3d=1000.0, "
                       "element_change=True, seed='c1ccccc1', shift=True)>")