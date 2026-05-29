from importlib import reload
import numpy as np
import pytest

from lomap.graphgen import GraphGen
from lomap import fp


def test_pick_lead_hub_None_deprecation():
    wmsg = "Support for passing hub as"
    with pytest.warns(DeprecationWarning, match=wmsg):
        _ = GraphGen.pick_lead(
            hub="None",
            names=["foo", "bar", "foobar"],
            strict_mtx=np.array([[1, 7, 3], [7, 4, 5], [3, 5, 2]]),
        )


def test_generate_initial_subgraph_list_lead_error():
    msg = "`lead_index` must be defined if using the fast map option"

    with pytest.raises(ValueError, match=msg):
        _ = GraphGen.generate_initial_subgraph_list(
            fast_map=True, strict_mtx=None, ids=None, names=None, is_active=None, lead_index=None
        )


def test_fp_deprecation():
    msg = "The fp module and associated Figureprint class are deprecated"
    with pytest.warns(DeprecationWarning, match=msg):
        reload(fp)
