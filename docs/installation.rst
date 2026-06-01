Installation
============

Latest Release
--------------

You can install lomap using conda (or mamba). Note the package name is ``lomap2``:

.. parsed-literal::
    conda install -c conda-forge lomap2

Development Version
-------------------
Alternatively, you can install the development version of ``lomap`` directly from the ``main`` branch of this repository.

First install the package dependencies using conda (or mamba) in a virtual environment with:

.. parsed-literal::
    conda env create -f environment.yaml
    conda activate lomap-env


Then install ``lomap`` locally with:

.. parsed-literal::
    pip install -e .


Optional Dependencies
---------------------

``lomap2`` has optional software dependencies which can be installed to extend
the package's capabilities:


gufe
~~~~

The ``lomap2`` package has optional bindings for the `gufe <https://gufe.openfree.energy/en/latest/`_
package. These bindings allow the atom mapping, network planning and edge
scoring functionality in ``lomap2`` to be used seemlessly with other components
of the `Open Free Energy ecosystem <https://openfree.energy/projects>`_.

This package can be installed via conda:

.. parsed-literal::
    conda install -c conda-forge gufe

Please see the :ref:`gufe bindings api documentation <gufe bindings api>` for
more information on how to use them.

pygraphviz
~~~~~~~~~~

The :class:`GraphGen` class can optionally plot the network graph using its
:meth:`draw` method. To do this, you will need to install
`pygraphviz <https://pygraphviz.github.io/>`_.

This can be done using conda:

.. parsed-literal::
    conda install -c conda-forge pygraphviz

or using `pip`:

.. parsed-literal::
    python -m pip install pygraphviz
