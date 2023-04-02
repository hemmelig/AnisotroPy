# -*- coding: utf-8 -*-
"""
Module that supplies various utility functions and classes.

:copyright:
    2023, AnisotroPy developers.
:license:
    GNU General Public License, Version 3
    (https://www.gnu.org/licenses/gpl-3.0.html)

"""

import logging
import sys
import time
from datetime import datetime
from functools import wraps
from itertools import tee

import numpy as np


log_spacer = "=" * 110


def azinc2vec(ψ, θ):
    """
    Return a 3-vector with the direction defined by (ψ, θ)

    Parameters
    ----------
    ψ : float
        Azimuth from x1 towards x2. Units of degrees.
    θ : float
        Inclination from the x1-x2 plane towards x3. Units of degrees.

    Returns
    -------
    x : array
        Vector pointing along (ψ, θ).

    """

    ψ, θ = np.deg2rad([ψ, θ])

    return np.array([np.cos(ψ) * np.cos(θ), -np.sin(ψ) * np.cos(θ), np.sin(θ)])


def rotate2xy(vector, azimuth, inclination):
    """
    Rotate a 3-vector into the x-y plane.

    Parameters
    ----------
    azimuth : list of float or float
        Azimuth, ψ, from x1 towards x2, in degrees.
    inclination : list of float or float
        Inclination, θ, from the x1-x2 plane towards x3, in degrees.

    Returns
    -------
    x : array
        Vector within the x-y plane.

    """

    ψ, θ = np.deg2rad([azimuth, inclination])

    r1 = np.array([[np.cos(ψ), np.sin(ψ), 0], [-np.sin(ψ), np.cos(ψ), 0], [0, 0, 1]])
    r2 = np.array([[np.cos(θ), 0, -np.sin(θ)], [0, 1, 0], [np.sin(θ), 0, np.cos(θ)]])

    return np.dot(np.dot(vector, r1), r2)


def make_directories(run, subdir=None):
    """
    Make run directory, and optionally make subdirectories within it.

    Parameters
    ----------
    run : `pathlib.Path` object
        Location of parent output directory, named by run name.
    subdir : str, optional
        subdir to make beneath the run level.

    """

    run.mkdir(exist_ok=True)

    if subdir:
        new_dir = run / subdir
        new_dir.mkdir(exist_ok=True, parents=True)


def logger(logstem, log, loglevel="info"):
    """
    Simple logger that will output to both a log file and stdout.

    Parameters
    ----------
    logstem : str
        Filestem for log file.
    log : bool
        Toggle for logging - default is to only print information to stdout.
        If True, will also create a log file.
    loglevel : str, optional
        Toggle for logging level - default is to print only "info" messages to
        log. To print more detailed "debug" messages, set to "debug".

    """

    level = logging.DEBUG if loglevel == "debug" else logging.INFO

    if log:
        now = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        logfile = logstem.parent / f"{logstem.name}_{now}"
        logfile.parent.mkdir(exist_ok=True, parents=True)
        handlers = [
            logging.FileHandler(str(logfile.with_suffix(".log"))),
            logging.StreamHandler(sys.stdout),
        ]
    else:
        handlers = [logging.StreamHandler(sys.stdout)]

    logging.basicConfig(level=level, format="%(message)s", handlers=handlers)


def pairwise(iterable):
    """Utility to iterate over an iterable pairwise."""

    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


def timeit(*args_, **kwargs_):
    """Function wrapper that measures the time elapsed during its execution."""

    def inner_function(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            ts = time.time()
            result = func(*args, **kwargs)
            msg = " " * 21 + f"Elapsed time: {time.time() - ts:6f} seconds."
            try:
                if args_[0] == "info":
                    logging.info(msg)
            except IndexError:
                logging.debug(msg)
            return result

        return wrapper

    return inner_function


def sample_sphere(n):
    """n even samples over a half sphere"""
    # https://stackoverflow.com/questions/9600801/evenly-distributing-n-points-on-a-sphere
    idx = np.arange(0, n, dtype=float) + 0.5
    incs = ((np.pi - np.arccos(idx / n)) * 180 / np.pi) % 90
    azis = ((np.pi * (1 + 5**0.5) * idx) * 180 / np.pi) % 360
    return incs, azis
