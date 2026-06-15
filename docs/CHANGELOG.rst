================
Lomap Change Log
================

.. current developments

v3.3.0
====================

**Changed:**

* Changed the default ``charge_changes_score`` from 0.0 to 0.1 in the gufe bindings to enable connected networks for
  ligands of different net charge by default.
* By default ``generate_lomap_network`` (gufe bindings) now
  fails if it cannot create a network where all the nodes
  are at least indirectly connected to each other. This behaviour
  can be controlled using the ``allow_disconnected`` keyword (PR #154).
* The `gufe` package is now an optional dependency of Lomap (PR #134).
* Minimum tested Python version has been raised on Python 3.11 and gufe
  raised to v1.0 (PR #124).
* The ``seed`` argument to ``LomapAtomMapper`` and ``MCS`` is now
  ``None`` by default, which has the same behaviour as the the
  previous default ``""`` (an empty string - which means no seed).
  This is purely an asthetic change to make the code more Pythonic.
  PR #156.

**Deprecated:**

* The ``dbmol`` CLI is deprecated and will be removed in the next
  major release of Lomap (Issue #138).
* The ``fp`` module and associated ``Figureprint`` class are deprecated
  and will be removed in the next release of Lomap (Issue #129).
* Deprecated the use of `str(None)` as an input to the `hub` keyword argument
  in `graphgen.GraphGen.pick_lead`. This option will be removed in the next release
  of Lomap, please use `None` instead. PR #125.

**Removed:**

* Testing code which was available by directly calling `mcs.py` has
  been removed (PR #137)

**Fixed:**

* The gufe bindings for the MNCAR score, which encodes the "minimum number of
  common atoms rule", incorrectly compared the common atoms threadshol ``ths``
  as `>` instead of `>=`. This is now fixed. Issue #147.
* Improvements to the optional ``pygraphviz`` dependency.
* Various historical typing issues, including modernizing various
  `rdkit.Chem.AllChem` calls. PR #125.


