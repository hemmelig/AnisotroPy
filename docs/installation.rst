Installation
============
:mod:`AnisotroPy` is a predominantly Python package with some routines written and optimised in C. These are built and linked to AnisotroPy at installation, which means you will need to ensure that there is a suitable compiler available (see :ref:`C compilers`).

Supported operating systems
---------------------------
AnisotroPy was developed and tested on Ubuntu 20.04, with the intention of being "platform agnostic". As of January 2022, the package has been successfully built and run on:

- Ubuntu 20.04

Prerequisites
-------------
AnisotroPy supports Python 3.6 or newer. We recommend using Anaconda as a package manager and environment management system to isolate and install the specific dependencies of AnisotroPy. Instructions for downloading and installing Anaconda can be found `here <https://docs.anaconda.com/anaconda/install/>`_. If drive space is limited, consider using Miniconda instead, which ships with a minimal collection of useful packages.

Setting up an environment
*************************
Using conda, you can use our ``anisotropy.yml`` file to create and activate a minimally complete environment:

.. code-block:: bash
    
    conda env create -f anisotropy.yml
    conda activate anisotropy

This will install the explicit dependencies of AnisotroPy (as well as some additional sub-dependencies/useful packages). The full list of dependencies (and versions, where relevant) is:

- matplotlib
- numpy
- obspy
- pandas
- pyproj
- scipy

.. note:: Any version pins are subject to change. We defer to ObsPy to select suitable versions for NumPy/SciPy.

C compilers
***********
In order to install and use AnisotroPy, you will need a C compiler that will build the migration extension library.

If you already have a suitable compiler (e.g. gcc, MSVC) at the OS level, then you can proceed to the Installing section.

If you do not, or to be sure, we recommend installing a compiler using conda. Instructions for doing this on :ref:`Linux` and :ref:`macOS` operating systems are given below.

Linux
#####
We recommend installing the GNU compiler collection (GCC, which previously stood for the GNU C Compiler) `through conda <https://docs.conda.io/projects/conda-build/en/latest/resources/compiler-tools.html>`_.

.. code-block:: bash
    
    conda install gcc_linux-64

It is generally useful to install compilers at the OS level, including a C++ compiler (e.g. gxx), which is required to build the scikit-fmm package.

Once installed, you can proceed with the AnisotroPy :ref:`installation <Installing>`.

macOS
#####
As with Linux, we recommend installing GCC through conda.

.. code-block:: bash
    
    conda install gcc

.. note:: We have not yet tested compiling and/or running AnisotroPy against the Clang compiler.

Alternatively, installation of compilers at an OS level can be done using ``Homebrew``, `a package manager for macOS <https://brew.sh/>`_. It is then as simple as:

.. code-block:: bash
    
    brew install gcc

.. note:: To install gcc-9, replace ``gcc`` with ``gcc@9``

Once installed, you can proceed with the AnisotroPy :ref:`installation <Installing>`.

Windows
#######
Compilation and linking of the C extensions has been successful using the Microsoft Visual C++ (MSVC) build tools. We strongly recommend that you download and install these tools in order to use AnisotroPy. You can either install Visual Studio in its entirety, or just the Build Tools - `available here <https://visualstudio.microsoft.com/downloads/>`_. You will need to restart your computer once the installation process has completed. We recommend using the anaconda command line interface (unix shell-like) to install AnisotroPy over command prompt.

.. warning:: AnisotroPy has been tested and validated on Windows, but there may yet remain some unknown issues. If you encounter an issue (and/or resolve it), please let us know!

Once installed, you can proceed with the AnisotroPy :ref:`installation <Installing>`.

Installing
----------
There are several ways to get a copy of AnisotroPy:

From source
***********
`Clone the repository <https://help.github.com/en/github/creating-cloning-and-archiving-repositories/cloning-a-repository>`_ from our `GitHub <https://github.com/hemmelig/AnisotroPy>`_ (note: you will need ``git`` installed on your system), or alternatively download the source code directly through the GitHub web interface. Once you have a local copy, navigate to the new AnisotroPy directory and run (ensuring your environment is activated):

.. code-block:: bash
    
    pip install .

You can optionally pass a ``-e`` argument to install the package in 'editable' mode.

You should now be able to import AnisotroPy within a Python session:

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
In order to test your installation, you will need to have cloned the GitHub repository. This will ensure you have all of the required benchmarked data (which is not included in pip/conda installs).

Iceland icequake test
*********************
Navigate to ``AnisotroPy/examples/.`` and run the example scripts in the following order:

.. code-block:: bash
    
    python iceland_lut.py
    python iceland_detect.py
    python iceland_trigger.py
    python iceland_locate.py

Once these have all run successfully, navigate to ``AnisotroPy/tests`` and run:

.. code-block:: bash
    
    python test_benchmarks.py

This should execute with no failed tests.
