Getting started using Lomap
===========================

Lomap plans networks of relative free energy calculations across a set of
ligands. Its interface provides
`gufe <https://gufe.openfree.energy/en/latest/>`_ bindings, so Lomap's atom
mapping, scoring, and network planning interoperate with the rest of the
`Open Free Energy <https://openfree.energy>`_ ecosystem. The workflow is:
load your ligands, propose atom mappings between them, score
those mappings, and assemble them into a network. Each step is available on its
own, but most users only need :func:`~lomap.generate_lomap_network`, covered at
the end of this guide.

These bindings require the optional ``gufe`` dependency
(see :doc:`installation`).

All examples below use these two example ligands:

.. code-block:: python

    import importlib.resources

    from gufe import SmallMoleculeComponent

    data = importlib.resources.files("lomap.tests.data")
    ligands = [
        SmallMoleculeComponent.from_sdf_file(data / name)
        for name in ["lig_41.sdf", "lig_74.sdf"]
    ]

Generating mappings
-------------------

A :class:`~lomap.LomapAtomMapper` proposes an atom mapping between a pair of
ligands. ``suggest_mappings`` yields ``LigandAtomMapping`` objects, or nothing
if the ligands share no suitable common substructure.

.. code-block:: python

    from lomap import LomapAtomMapper

    mapper = LomapAtomMapper()
    # suggest_mappings is a generator, so wrap it in list() to pull out all the possible mappings
    mappings = list(mapper.suggest_mappings(ligands[0], ligands[1]))

    mapping = mappings[0]

    # atom indices in ligand A mapped onto the corresponding atoms in ligand B
    print(mapping.componentA_to_componentB)

The mapper can be tuned with options such as ``time``,
``threed``, ``max3d``, ``element_change``, ``seed``, and ``shift``. See the
:ref:`API reference <gufe bindings api>` for details.

Generating scores
-----------------

A scorer takes a ``LigandAtomMapping`` and returns a value from ``0.0`` (worst)
to ``1.0`` (best).
:func:`~lomap.default_lomap_score` is the standard Lomap
score, a product of sub-scores penalizing changes such as broken rings, altered
ring sizes, and net-charge differences.

.. code-block:: python

    from lomap import default_lomap_score

    score = default_lomap_score(mapping)
    print(f"{score:.3f}")

The score of ``0.095`` for this pair is expected. ``lig_41`` and ``lig_74``
differ at both ends; fused ring system is created on one end and the
heteroaryl head changes, leaving a relatively small shared core, along with
element and hybridization changes.

The sub-scores live in ``lomap.gufe_bindings.scorers`` if you need a custom
scorer.

Generating networks
-------------------

:func:`~lomap.generate_lomap_network` combines the previous steps. Given a set of ligands, a
mapper, and a scorer, it proposes and scores edges and assembles a connected
``LigandNetwork``.

.. code-block:: python

    from lomap import LomapAtomMapper, default_lomap_score, generate_lomap_network

    network = generate_lomap_network(
        ligands=ligands,
        mappers=LomapAtomMapper(),
        scorer=default_lomap_score,
    )

    print(f"{len(network.nodes)} ligands, {len(network.edges)} edges")

``network.nodes`` are the ligands and ``network.edges`` the scored mappings.
Options such as ``distance_cutoff``, ``require_cycle_covering``, and ``radial``
are documented in the :ref:`API reference <gufe bindings api>`.

Command line interface
----------------------

The legacy ``lomap`` command-line tool is deprecated and will be removed in the
next major release; use :func:`~lomap.generate_lomap_network` instead. See
:doc:`legacy`.
