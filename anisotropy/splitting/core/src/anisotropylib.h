/*
 * =============================================================================
 *
 *       Filename:  anisotropylib.h
 *
 *        Purpose:  Header file to bring together definitions used in the
 *                  anisotropy.c library.
 *
 *      Copyright:  2021--2022, AnisotroPy developers.
 *        License:  GNU General Public License, Version 3
 *                  (https://www.gnu.org/licenses/gpl-3.0.html)
 *
 * =============================================================================
 */

#include <gsl/gsl_vector.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>

#ifndef _OPENMP
    /* Generate a compiler error to stop the build */
    mustLinkOpenMP
#endif

#define ROUND(a) (floor((a) + 0.5))

typedef struct
{
    double fast_inc;
    double dt_inc;
    double dt_max;
} run_parameters;

typedef struct
{
    int32_t start;
    int32_t npts;
    double delta;
} window_parameters;

typedef struct
{
    int32_t start_i;
    int32_t end_i;
    int32_t n_start;
    int32_t n_end;
} window_trials;

void trial_windows_sc91(
    double *x,
    double *y,
    int32_t npts,
    double delta,
    run_parameters *params,
    window_trials *window,
    int32_t threads
);

void split_sc91(
    gsl_vector *x,
    gsl_vector *y,
    run_parameters *params,
    window_parameters *window
);
