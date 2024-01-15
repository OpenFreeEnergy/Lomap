from lomap import LomapAtomMapper


def test_to_dict_roundtrip():
    ref_vals = {
        'time': 19,
        'threed': False,
        'max3d': 999.0,
        'element_change': False,
        'seed': 'CC',
        'shift': False,
    }

    m = LomapAtomMapper(
        **ref_vals
    )

    d = m.to_dict()

    m2 = LomapAtomMapper.from_dict(d)

    assert m2
    assert m2.time == ref_vals['time']
    assert m2.threed == ref_vals['threed']
    assert m2.max3d == ref_vals['max3d']
    assert m2.element_change == ref_vals['element_change']
    assert m2.seed == ref_vals['seed']
    assert m2.shift == ref_vals['shift']
