<p align="center">
  <!-- DOI -->
  <a href="https://doi.org/10.5281/zenodo.5931586">
    <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.5931586.svg" />
  </a>
  <!-- ReadTheDocs -->
  <a href="https://seismicanisotropy.readthedocs.io/en/latest">
    <img src="https://readthedocs.org/projects/seismicanisotropy/badge/?version=latest" />
  </a>
  <!-- License -->
  <a href="https://www.gnu.org/licenses/gpl-3.0">
    <img src="https://img.shields.io/badge/License-GPLv3-blue.svg" />
  </a>
</p>

<p align="center">
  <a href="https://seismicanisotropy.readthedocs.io/en/latest/index.html">AnisotroPy</a> is a cross-platform Python package for the study of <a href="https://en.wikipedia.org/wiki/Seismic_anisotropy">seismic anisotropy</a>.
</p>

<p align="center">
  <img src="https://raw.githubusercontent.com/hemmelig/AnisotroPy/main/docs/img/AnilogoBig.png" width="80%" />
</p>

Key features
------------
The goal of AnisotroPy is to provide straightforward access to a suite of routines and utilities for the study of seismic anisotropy. The package is primarily designed for both scripting purposes, but also includes a number of easy-to-use command-line utilities. Currently, we support:

- Effective media modelling
- Anisotropic layer fitting to observations
- Visualisation of shear-wave splitting measurements
- Shear-wave splitting analysis

A sibling project—[AnisotropIO](https://github.com/hemmelig/AnisotropIO)—provides a standardised file format for shear-wave splitting measurements, as well as a range of parsers for other existing codes. The aim here is to make such analyses reproducible and straightforward to ingest into, for example, meta-studies.

Documentation
-------------
Documentation for AnisotroPy is hosted [here](https://seismicanisotropy.readthedocs.io/en/latest/index.html).

Installation
------------
AnisotroPy requires Python version 3.8 and above. Installation of AnisotroPy, including all dependencies, can be done using pip:

```console
pip install anisotropy
```

For further information regarding installation—including virtual environment management and installation from source—please consult [our documentation](https://seismicanisotropy.readthedocs.io/en/latest/installation.html).

Usage
-----
We are working on tutorials covering how each individual aspect of the package works.

This is a work in progress - [see our documentation for full details](https://seismicanisotropy.readthedocs.io/en/latest/tutorials.html).

For a more comprehensive demonstration of the options available, see the [template scripts](examples/template_scripts).

Citation
--------
If you use AnisotroPy in your work, please cite the following:

AnisotroPy Developers (2022). AnisotroPy: v0.0.1 (v0.0.1). Zenodo. https://doi.org/10.5281/zenodo.5931586

Contributing to AnisotroPy
--------------------------
Contributions to AnisotroPy are welcomed. The first stop should be to reach out, either directly or—preferably—via the GitHub Issues panel, to discuss the proposed changes. Next, simply fork the AnisotroPy repository, make your changes/add your new contribution, then make a [pull request](https://help.github.com/articles/about-pull-requests/). All contributors to AnisotroPy will be listed as authors on the releases.

Bug reports, suggestions for new features and enhancements, and even links to projects that have made use of AnisotroPy are most welcome.

See our [contributions page](https://github.com/hemmelig/AnisotroPy/blob/main/.github/CONTRIBUTING.md) for more information.

Contact
-------
Any comments/questions can be directed to:
* **Conor Bacon** - cbacon [ at ] ldeo.columbia.edu

License
-------
AnisotroPy is **free** and **open source**, distributed under the GPLv3 License. Please see the [LICENSE](LICENSE) file for a complete description of the rights and freedoms that this provides the user.
