Installation
============

Latest Release
--------------

You can install lomap using conda (or mamba). Note the package name is ``lomap2``:

.. parsed-literal::
    conda install -c conda-forge lomap2

Development Version
-------------------
Alternatively, you can install the development version of ``lomap``` directly from the ``main`` branch of this repository.

First install the package dependencies using conda (or mamba) in a virtual environment with:

.. parsed-literal::
    conda env create -f environment.yaml
    conda activate lomap-env


Then install ``lomap`` locally with:

.. parsed-literal::
    pip install -e .
