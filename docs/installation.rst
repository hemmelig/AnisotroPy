Installation
============
:mod:`AnisotroPy` is a predominantly Python package, though it will in future have some routines written and optimised in C.
Supported operating systems
---------------------------
AnisotroPy was developed and tested on Ubuntu 20.04, with the intention of being "platform agnostic". As of June 2022, the package has been successfully built and run on:

- Ubuntu 20.04

Prerequisites
-------------
AnisotroPy supports Python 3.7 or newer. We recommend using Anaconda as a package manager and environment management system to isolate and install the specific dependencies of AnisotroPy. Instructions for downloading and installing Anaconda can be found `here <https://docs.anaconda.com/anaconda/install/>`_. If drive space is limited, consider using Miniconda instead, which ships with a minimal collection of useful packages.

Setting up an environment
*************************
Using conda, you can use our ``environment.yml`` file to create and activate a minimally complete environment:

.. code-block:: bash
    
    conda env create -f environment.yml
    conda activate anisotropy

This will install the explicit dependencies of AnisotroPy. The full list of dependencies (and versions, where relevant) is:

- matplotlib
- numpy
- obspy
- pandas

.. note:: Any version pins are subject to change. We defer to ObsPy to select suitable versions for NumPy/SciPy.

Installing
----------
There are several ways to get a copy of AnisotroPy:

From source
***********
`Clone the repository <https://help.github.com/en/github/creating-cloning-and-archiving-repositories/cloning-a-repository>`_ from our `GitHub <https://github.com/hemmelig/AnisotroPy>`_ (note: you will need ``git`` installed on your system), or alternatively download the source code directly through the GitHub web interface. Once you have a local copy, navigate to the new AnisotroPy directory and run (ensuring your environment is activated):

.. code-block:: bash
    
    pip install .

You can optionally pass a ``-e`` argument to install the package in 'editable' mode.

You should now be able to import AnisotroPy within a Python session (but not in the git repo root directory!):

.. code-block:: bash
    
    python
    >>> import anisotropy

pip install
***********
We will be linking the package to PyPI (the Python Package Index) soon, after which you will be able to use the following command to install the package:

.. code-block:: bash
    
    pip install anisotropy

conda install
*************
We hope to link the package with the conda forge soon, after which you will be able to use the following command to install the package:

.. code-block:: bash
    
    conda install -c conda-forge anisotropy

Testing your installation
-------------------------
Navigate to ``AnisotroPy/tests/.`` and run:

.. code-block:: bash
    
    test_materials.py

This should execute with no failed tests. More tests are on the way!
