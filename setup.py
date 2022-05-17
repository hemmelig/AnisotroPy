# -*- coding: utf-8 -*-
"""
AnisotroPy - a Python toolkit for the study of seismic anisotropy.

:copyright:
    2021--2022, AnisotroPy developers.
:license:
    GNU General Public License, Version 3
    (https://www.gnu.org/licenses/gpl-3.0.html)

"""

import os
import pathlib
import platform
# import re
import shutil
import sys

from distutils.ccompiler import get_default_compiler
# from pkg_resources import get_build_platform
from setuptools import Extension, find_packages, setup
from setuptools.command.build_ext import build_ext


__version__ = "0.0.1"

# The minimum python version which can be used to run AnisotroPy
MIN_PYTHON_VERSION = (3, 6)

# Fail fast if the user is on an unsupported version of Python.
if sys.version_info < MIN_PYTHON_VERSION:
    msg = (f"AnisotroPy requires python version >= {MIN_PYTHON_VERSION}"
           f" you are using python version {sys.version_info}")
    print(msg, file=sys.stderr)
    sys.exit(1)

# Check if we are on RTD and don't build extensions if we are.
READ_THE_DOCS = os.environ.get("READTHEDOCS", None) == "True"
if READ_THE_DOCS:
    try:
        environ = os.environb
    except AttributeError:
        environ = os.environ

    environ[b"CC"] = b"x86_64-linux-gnu-gcc"
    environ[b"LD"] = b"x86_64-linux-gnu-ld"
    environ[b"AR"] = b"x86_64-linux-gnu-ar"

# Directory of the current file
SETUP_DIRECTORY = pathlib.Path.cwd()
DOCSTRING = __doc__.split("\n")

# Check for MSVC (Windows)
if platform.system() == "Windows" and (
        "msvc" in sys.argv or
        "-c" not in sys.argv and
        get_default_compiler() == "msvc"):
    IS_MSVC = True
else:
    IS_MSVC = False

INSTALL_REQUIRES = [
    "matplotlib",
    "numpy",
    "obspy",
    "pandas",
    "pyproj",
    "scipy"
]

if READ_THE_DOCS:
    EXTRAS_REQUIRES = {
        "docs": [
            "Sphinx >= 1.8.1",
            "docutils"
        ]
    }
else:
    EXTRAS_REQUIRES = {}

KEYWORDS = [
    "array", "anisotropy", "seismic", "seismology", "earthquake", "splitting",
    "modelling", "ObsPy", "waveform", "seismic anisotropy", "processing"
]

# Monkey patch for MS Visual Studio
if IS_MSVC:
    # Remove 'init' entry in exported symbols
    def _get_export_symbols(self, ext):
        return ext.export_symbols
    from setuptools.command.build_ext import build_ext
    build_ext.get_export_symbols = _get_export_symbols


def export_symbols(path):
    """
    Required for Windows systems - functions defined in anisotropylib.def.
    """
    with (SETUP_DIRECTORY / path).open("r") as f:
        lines = f.readlines()[2:]
    return [s.strip() for s in lines if s.strip() != ""]


# def get_package_data():
#     package_data = {}
#     if not READ_THE_DOCS:
#         if IS_MSVC:
#             package_data["anisotropy.splitting.core"] = [
#                 "anisotropy/splitting/core/src/*.dll"
#             ]

#     return package_data


# def get_package_dir():
#     package_dir = {}
#     if IS_MSVC:
#         package_dir["anisotropy.splitting.core"] = str(
#             pathlib.Path("anisotropy") / "splitting" / "core"
#         )

#     return package_dir


# def get_include_dirs():
#     include_dirs = [
#         str(pathlib.Path.cwd() / "anisotropy" / "core" / "src"),
#         str(pathlib.Path(sys.prefix) / "include")
#     ]

#     if get_build_platform().startswith("freebsd"):
#         include_dirs.append("/usr/local/include")

#     return include_dirs


# def get_library_dirs():
#     library_dirs = []
#     if IS_MSVC:
#         library_dirs.append(str(pathlib.Path.cwd() / "anisotropy" / "core"))
#         library_dirs.append(str(pathlib.Path(sys.prefix) / "bin"))

#     library_dirs.append(str(pathlib.Path(sys.prefix) / "lib"))
    # if get_build_platform().startswith("freebsd"):
    #     library_dirs.append("/usr/local/lib")

    # return library_dirs


def get_extensions():
    """
    Config function used to compile C code into a Python extension.
    """
    extensions = []

    if READ_THE_DOCS:
        return extensions

    common_extension_args = {
        "include_dirs": [
            str(pathlib.Path.cwd() / "anisotropy" / "core" / "src"),
            str(pathlib.Path(sys.prefix) / "include")
        ],
        "library_dirs": [str(pathlib.Path(sys.prefix) / "lib")]
    }

    sources = [
        str(pathlib.Path("anisotropy") / "splitting/core/src/anisotropy.c")
    ]

    extra_link_args = []
    if IS_MSVC:
        extra_compile_args = ["/openmp", "/TP", "/O2"]
        common_extension_args["export_symbols"] = export_symbols(
            "anisotropy/splitting/core/src/anisotropylib.def"
        )
        common_extension_args["library_dirs"].extend(
            str(pathlib.Path.cwd() / "anisotropy" / "core"),
            str(pathlib.Path(sys.prefix) / "bin")
        )
    else:
        extra_link_args.extend(["-lm", "-lgsl", "-lgslcblas"])
        if platform.system() == "Darwin":
            extra_link_args.extend(["-lomp"])
        else:
            extra_link_args.extend(["-lgomp"])
        extra_compile_args = ["-fopenmp", "-fPIC"]#, "-Ofast"]

    common_extension_args["extra_link_args"] = extra_link_args
    common_extension_args["extra_compile_args"] = extra_compile_args

    extensions.append(Extension("anisotropy.splitting.core.src.anisotropylib",
                      sources=sources, **common_extension_args))

    return extensions


class CustomBuildExt(build_ext):
    def finalize_options(self):
        build_ext.finalize_options(self)

        if self.compiler is None:
            compiler = get_default_compiler()
        else:
            compiler = self.compiler

        if IS_MSVC:
            # Sort linking issues with init exported symbols
            def _get_export_symbols(self, ext):
                return ext.export_symbols

            build_ext.get_export_symbols = _get_export_symbols


def setup_package():
    """Setup package"""

    if READ_THE_DOCS:
        INSTALL_REQUIRES.append("mock")

    package_dir, package_data = {}, {}
    if IS_MSVC:
        package_dir["anisotropy.splitting.core"] = str(
            pathlib.Path("anisotropy") / "splitting" / "core"
        )
        package_data["anisotropy.splitting.core"] = [
            "anisotropy/splitting/core/src/*.dll"
        ]

    setup_args = {
        "name": "anisotropy",
        "version": __version__,
        "description": DOCSTRING[1],
        "long_description": "\n".join(DOCSTRING[3:]),
        "url": "https://github.com/hemmelig/AnisotroPy",
        "author": "The AnisotroPy Development Team",
        "author_email": "conor.bacon@esc.cam.ac.uk",
        "license": "GNU General Public License, Version 3 (GPLv3)",
        "classifiers": [
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Science/Research",
            "Topic :: Scientific/Engineering",
            "Natural Language :: English",
            "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
            "Operating System :: OS Independent",
            "Programming Language :: Python",
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Python :: 3.7",
            "Programming Language :: Python :: 3.8",
        ],
        "keywords": KEYWORDS,
        "install_requires": INSTALL_REQUIRES,
        "extras_require": EXTRAS_REQUIRES,
        "zip_safe": False,
        # "packages": find_packages(),
        "packages": ["anisotropy",
                     "anisotropy.effective_modelling",
                     "anisotropy.materials",
                     "anisotropy.splitting",
                     "anisotropy.splitting.core",
                     "anisotropy.utils"],
        "ext_modules": get_extensions(),
        "package_data": package_data,
        "package_dir": package_dir
    }

    shutil.rmtree(str(SETUP_DIRECTORY / "build"), ignore_errors=True)

    setup(**setup_args)


if __name__ == "__main__":
    # clean --all does not remove extensions automatically
    if "clean" in sys.argv and "--all" in sys.argv:
        # Delete complete build directory
        path = SETUP_DIRECTORY / "build"
        shutil.rmtree(str(path), ignore_errors=True)

        # Delete all shared libs from clib directory
        path = SETUP_DIRECTORY / "anisotropy" / "core" / "src"
        for filename in path.glob("*.pyd"):
            filename.unlink(missing_ok=True)
        for filename in path.glob("*.so"):
            filename.unlink(missing_ok=True)
        for filename in path.glob("*.dll"):
            filename.unlink(missing_ok=True)
    else:
        setup_package()
