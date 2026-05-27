import numpy as np
import pytest

from lomap.graphgen import GraphGen


def test_pick_lead_hub_None_deprecation():
    wmsg = "Support for passing hub as"
    with pytest.warns(DeprecationWarning, match=wmsg):
        GraphGen.pick_lead(hub="None", names=["foo", "bar", "foobar"], strict_mtx=np.array([[1, 7, 3], [7, 4, 5], [3, 5, 2]])) 
