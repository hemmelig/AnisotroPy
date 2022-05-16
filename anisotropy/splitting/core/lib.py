# -*- coding: utf-8 -*-
"""
Bindings for the C library functions.

:copyright:
    2021--2022, AnisotroPy developers.
:license:
    GNU General Public License, Version 3
    (https://www.gnu.org/licenses/gpl-3.0.html)

"""

import numpy as np
import numpy.ctypeslib as clib

from anisotropy.splitting.core.libnames import _load_cdll


anisotropylib = _load_cdll("anisotropylib")

run_parameters = np.dtype(
    [
     ("fast_inc", np.float64),
     ("dt_inc", np.float64),
     ("dt_max", np.float64),
    ],
    align=True
)

window_parameters = np.dtype(
    [
     ("start_i", clib.ctypes.c_int32),
     ("end_i", clib.ctypes.c_int32),
     ("n_start", clib.ctypes.c_int32),
     ("n_end", clib.ctypes.c_int32)
    ],
    align=True
)

anisotropylib.trial_windows_sc91.argtypes = [
    clib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
    clib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
    clib.ctypes.c_int32,
    clib.ctypes.c_double,
    clib.ndpointer(dtype=run_parameters, flags="C_CONTIGUOUS"),
    clib.ndpointer(dtype=window_parameters, flags="C_CONTIGUOUS"),
    clib.ctypes.c_int32,
]


def evaluate_splitting(stream, run_params, window_params, threads, method="sc91"):
    """
    Compute the pair of splitting parameters, (phi, dt), that best remove any
    shear-wave splitting for a given event-station pair.

    So far, just the eigenvalue minimisation method of Silver and Chan, 1991,
    has been implemented, but there is room to incorporate other methods later.

    Parameters
    ----------
    stream : `obspy.Stream` object
        Pre-processed waveform data in the form of a Stream containing the
        east-west (x) and north-south (y) components.
    run_params : dict
        Table of parameters specifying the run
    window_params : dict
        Table of parameters defining the 
    threads : int
        Number of threads with which to perform the scan.
    method : str, optional
        Specify the method to use - default Silver and Chan, 1991.

    """

    # Build C data struct for the run parameters
    run_struct = np.empty(1, dtype=run_parameters)
    run_struct[:] = (
        run_params["fast_inc"],
        run_params["dt_inc"],
        run_params["dt_max"]
    )

    # Build C data struct for the window parameters
    window_struct = np.empty(1, dtype=window_parameters)
    window_struct[:] = (
        window_params["start_i"],
        window_params["end_i"],
        window_params["n_start"],
        window_params["n_end"]
    )

    # Ensure input waveform data are contiguous in memory
    east = np.ascontiguousarray(
        stream.select(component="E")[0].data,
        dtype=np.double
    )
    north = np.ascontiguousarray(
        stream.select(component="N")[0].data,
        dtype=np.double
    )

    anisotropylib.trial_windows_sc91(
        north,
        east,
        clib.ctypes.c_int32(len(east)),
        clib.ctypes.c_double(stream[0].stats.delta),
        run_struct,
        window_struct,
        clib.ctypes.c_int32(threads)
    )

