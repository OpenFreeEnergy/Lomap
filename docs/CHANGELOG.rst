================
Lomap Change Log
================

.. current developments

v3.3.0
====================

This release offers a series of usability improvements and prepares ``Lomap`` for
its next major release.

**Changed:**

* Changed the default ``charge_changes_score`` from 0.0 to 0.1 in the gufe
  bindings to enable connected networks for ligands of different net charge by default
  (`PR #83 <https://github.com/OpenFreeEnergy/Lomap/pull/83>`_).
* By default ``generate_lomap_network`` (gufe bindings) now
  fails if it cannot create a network where all the nodes
  are at least indirectly connected to each other. This behaviour
  can be controlled using the ``allow_disconnected`` keyword
  (`PR #154 <https://github.com/OpenFreeEnergy/Lomap/pull/154>`_).
* The `gufe` package is now an optional dependency of Lomap
  (`PR #134 <https://github.com/OpenFreeEnergy/Lomap/pull/134>`_).
* Minimum tested Python version has been raised on Python 3.11 and gufe
  raised to v1.0 (`PR #124 <https://github.com/OpenFreeEnergy/Lomap/pull/124>`_).
* The ``seed`` argument to ``LomapAtomMapper`` and ``MCS`` is now
  ``None`` by default, which has the same behaviour as the
  previous default ``""`` (an empty string - which means no seed).
  This is purely an aesthetic change to make the code more Pythonic.
  (`PR #156 <https://github.com/OpenFreeEnergy/Lomap/pull/156>`_).

**Deprecated:**

* The ``dbmol`` CLI is deprecated and will be removed in the next
  major release of Lomap (`Issue #138 <https://github.com/OpenFreeEnergy/Lomap/issues/138>`_).
* The ``fp`` module and associated ``Figureprint`` class are deprecated
  and will be removed in the next release of Lomap
  (`Issue #129 <https://github.com/OpenFreeEnergy/Lomap/issues/129>`_).
* Deprecated the use of `str(None)` as an input to the `hub` keyword argument
  in `graphgen.GraphGen.pick_lead`. This option will be removed in the next release
  of Lomap, please use `None` instead
  (`PR #125 <https://github.com/OpenFreeEnergy/Lomap/pull/125>`_).

**Removed:**

* Testing code which was available by directly calling `mcs.py` has
  been removed (`PR #137 <https://github.com/OpenFreeEnergy/Lomap/pull/137>`_).

**Fixed:**

* The gufe bindings for the MNCAR score, which encodes the "minimum number of
  common atoms rule", incorrectly compared the common atoms threshold ``ths``
  as ``>`` instead of ``>=``. This is now fixed.
  (`Issue #147 <https://github.com/OpenFreeEnergy/Lomap/issues/147>`_).
* Improvements to the optional ``pygraphviz`` dependency
  (`PR #141 <https://github.com/OpenFreeEnergy/Lomap/pull/141>`_).
* Various historical typing issues, including modernizing various
  `rdkit.Chem.AllChem` calls
  (`PR #125 <https://github.com/OpenFreeEnergy/Lomap/pull/125>`_).



v3.2.1
====================

**Fixed:**

* Fixed bug where ``generate_lomap_network()`` would throw an error if
  ``molecules`` was passed in instead of ``ligands``.



v3.2.0
====================

**Deprecated:**

* Replaced the ``molecules`` argument with ``ligands`` in ``generate_lomap_network()``.
  Argument name ``molecules`` will be deprecated.



v2.3.0
====================

**Added:**

* Added shift option to MCS and cli (-s or --shift), this defaults to ``True``. 
  In combination with ``threed``/``-3`` this option controls if the two
  structures are shifted on top of each other (using only translation and no rotation) before
  geometric mismatch is judged.  This was always used as default until now.
  In cases where bad symmetric options are being chosen, using ``shift=False`` might resolve issues
  (`Issue #26 <https://github.com/OpenFreeEnergy/Lomap/issues/26>`_).



v2.2.0
====================

**Added:**

* Added ``use_common_core`` option to ``DBMolecules`` and ``-C`` / ``--common-core`` option to CLI.
  This is on by default and speeds up network creation by 10-50x
  (`PR #25 <https://github.com/OpenFreeEnergy/Lomap/pull/25>`_).

**Fixed:**

* Fixed ``DBMolecules`` iteration in Python v3 (`PR #18 <https://github.com/OpenFreeEnergy/Lomap/pull/18>`_).



v2.1.0
====================

**Added:**

* Added ``element_change`` keyword to ``MCS`` class. This toggles if elemental
  changes are allowed within a mapping. Defaults to ``True`` (this was the previous behaviour).
* Added ``-L`` switch command line usage to toggle element changes.



v2.0.0
====================

**Added:**

* Added regression tests for ``MCS`` class.

**Changed:**

* Replaced usage of of ``argparse.Namespace`` as input to class ``inits``
  (e.g. ``MCS``) in favour of dicts.
* Renamed various dunderscore variables to single underscore.
* Rewrote tests in pytest, using fixtures to iterate over cases.

**Fixed:**

* Fixed bug in explicit hydrogen atom mappings.



v1.0.0
====================

**Added:**

* Added versioning.
* Support for Python 3.5 and ``networkx`` v2.
* Improved visualization of output networks.


**Fixed:**

* Better enforcement of PEP8 style.

**Removed:**

* Support for Python v2.7.



v0.0.x
====================

This is an internal Alpha version.

**Added:**

* Basic functionality for preparing a Perturbation map for alchemical free energy calculations.
