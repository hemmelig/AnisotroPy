.. figure:: img/AnilogoBig.png
   :figwidth: 100 %
   :width: 90%
   :align: center
   :alt: AnisotroPy: a toolkit for the study of seismic anisotropy.

AnisotroPy
==========

:mod:`AnisotroPy` is a cross-platform Python package for the study of seismic anisotropy.

The goal of AnisotroPy is to provide straightforward access to a suite of routines and utilities for the study of seismic anisotropy. The package is primarily designed for both scripting purposes, but also includes a number of easy-to-use command-line utilities. Currently, we support:

- Effective media modelling
- Anisotropic layer fitting to observations
- Visualisation of shear-wave splitting measurements
- Shear-wave splitting analysis

A sibling project—|AnisotropIO|—provides a standardised file format for shear-wave splitting measurements, as well as a range of parsers for other existing codes. The aim here is to make such analyses reproducible and straightforward to ingest into, for example, meta-studies.

The source code for the project is hosted on |github|.

This package is written by the AnisotroPy developers, and is distributed under
the GPLv3 License, Copyright AnisotroPy developers 2023.

.. |github| raw:: html

    <a href="https://github.com/hemmelig/AnisotroPy" target="_blank">GitHub</a>

.. |AnisotropIO| raw:: html

   <a href="https://github.com/hemmelig/AnisotropIO" target="_blank">AnisotropIO</a>

Supported operating systems
---------------------------
AnisotroPy was developed and tested on Ubuntu 20.04, with the intention of being "platform agnostic". As of April 2023, the package has been successfully built and run on:

- Ubuntu 20.04/22.04
- macOS Monterey v12.5.1

Citation
--------
If you use AnisotroPy in your work, please cite the following:

AnisotroPy Developers (2022). AnisotroPy: v0.0.1 (v0.0.1). Zenodo. https://doi.org/10.5281/zenodo.5931586

We hope at some future point to publish a short companion paper in the Journal of Open Source Software.

Contact
-------
Any comments/questions can be directed to:

* **Conor Bacon** - cbacon [ at ] ldeo.columbia.edu

License
-------
AnisotroPy is **free** and **open-source**, distributed under the GPLv3 License. Please see the `LICENSE <https://www.gnu.org/licenses/gpl-3.0.html>`_ for a complete description of the rights and freedoms that this provides the user.

Contents:
---------

.. toctree::
   :numbered:
   :maxdepth: 1

   installation
   tutorials
   sourcecode
