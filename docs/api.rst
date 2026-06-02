.. _gufe bindings api:

Lomap `gufe <https://gufe.openfree.energy/en/latest/>`_ bindings API guide
==========================================================================

This section describes the optional `gufe <https://gufe.openfree.energy/en/latest/>`_
bindings API. These allow Lomap to be seemlessly used with other components
of the `Open Free Energy ecosystem <https://openfree.energy/projects>`_.

Lomap Atom Mapper
-----------------

.. autoclass:: lomap.LomapAtomMapper
   :members:
   :inherited-members:


Lomap Network Generator
-----------------------

.. autofunction:: lomap.generate_lomap_network

Lomap Atom Mapping scorers
--------------------------

For most users the :func:`default_lomap_score` is all you will need to use:

.. autofunction:: lomap.default_lomap_score

However, you can also use the underlying sub-scores that make up the :func:`default_lomap_score`:

.. autofunction:: lomap.gufe_bindings.scorers.ecr_score

.. autofunction:: lomap.gufe_bindings.scorers.mcsr_score

.. autofunction:: lomap.gufe_bindings.scorers.mncar_score

.. autofunction:: lomap.gufe_bindings.scorers.atomic_number_score

.. autofunction:: lomap.gufe_bindings.scorers.hybridization_score

.. autofunction:: lomap.gufe_bindings.scorers.sulfonamides_score

.. autofunction:: lomap.gufe_bindings.scorers.heterocycles_score

.. autofunction:: lomap.gufe_bindings.scorers.transmuting_methyl_into_ring_score

.. autofunction:: lomap.gufe_bindings.scorers.transmuting_ring_sizes_score
