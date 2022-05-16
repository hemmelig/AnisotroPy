/*
 * =============================================================================
 *
 *       Filename:  anisotropy.c
 *
 *        Purpose:  Routines for computing the shear-wave splitting parameters
 *                  for an event.
 *
 *      Copyright:  2021--2022, AnisotroPy developers.
 *        License:  GNU General Public License, Version 3
 *                  (https://www.gnu.org/licenses/gpl-3.0.html)
 *
 * =============================================================================
 */

#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_vector.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "anisotropylib.h"

#define PI 3.14159265358979323846

void trial_windows_sc91(
    double *x,
    double *y,
    int32_t npts,
    double delta,
    run_parameters *params,
    window_trials *window,
    int32_t threads)
{
    int32_t i, j, tmp; // Counters for window start/end trials
    window_parameters w_params;

    // Make vector views of the input data
    gsl_vector_view east = gsl_vector_view_array(x, npts);
    gsl_vector_view north = gsl_vector_view_array(y, npts);

    // Grab out window parameters and run over the n trials
    #pragma omp parallel for \
    private(j, tmp, w_params) \
    num_threads(threads)
    for (i=window->start_i; i<window->n_start+window->start_i; i++) {
        for (j=window->end_i; j<window->n_end+window->end_i; j++) {
            printf("%d %d\n", i, j);
            w_params = (window_parameters){.start = i, .npts = j-i+1, .delta = delta};
            split_sc91(&east.vector, &north.vector, params, &w_params);
        }
    }
}

void split_sc91(
    gsl_vector *x,
    gsl_vector *y,
    run_parameters *params,
    window_parameters *window)
{
    /*
    Perform a grid search over fast directions, phi, and delay time, dt.

    Args:
        x: Pointer to a GSL vector containing the east-west components of the
           seismogram.
    */
    int i, j;
    double fast, dt, current_min, dt_min, phi_min;
    int npts = x->size;

    // Build workspaces (pre-allocated matrices in memory) for analysis
    // Waveform rotation operators
    gsl_matrix *rotator = gsl_matrix_alloc(2, 2);
    gsl_matrix *rotated_data = gsl_matrix_alloc(2, npts);
    // Eigenvalue/vector analysis
    gsl_eigen_symmv_workspace *eigen_ws = gsl_eigen_symmv_alloc(2);
    gsl_matrix *evec = gsl_matrix_alloc(2, 2); /* eigenvectors */
    gsl_vector *eval = gsl_vector_alloc(2);    /* eigenvalues */
    // Covariance results
    gsl_matrix *C = gsl_matrix_alloc(2, 2);

    // Assign input data to a matrix
    gsl_matrix *data = gsl_matrix_alloc(2, npts);
    gsl_matrix_set_row(data, 0, x);
    gsl_matrix_set_row(data, 1, y);

    current_min = 10000.0;
    for (fast = -90, j = 0; fast <= 90; fast += params->fast_inc, j++)
    {
        // Build rotation operator
        double c = cos(fast * PI / 180.), s = sin(fast * PI / 180.);
        gsl_matrix_set(rotator, 0, 0, c);
        gsl_matrix_set(rotator, 0, 1, -s);
        gsl_matrix_set(rotator, 1, 0, s);
        gsl_matrix_set(rotator, 1, 1, c);

        // Rotate the input data
        gsl_blas_dgemm(
            CblasNoTrans,
            CblasNoTrans,
            1.0,
            rotator,
            data,
            0.0,
            rotated_data);

        // Added some tolerance to the upper limit as dealing with floating points
        for (dt = 0, i = 0; dt <= params->dt_max+0.0001; dt += params->dt_inc, i++)
        {
            // Calculate number of samples by which to shift data
            size_t dt_samples = (size_t)ROUND(dt / (2 * window->delta));

            // Window data (incorporating the 'time lag' and any window start time information)
            gsl_vector_view rot_x = gsl_matrix_subrow(
                rotated_data,
                0,
                window->start - dt_samples,
                window->npts
            );

            gsl_vector_view rot_y = gsl_matrix_subrow(
                rotated_data,
                1,
                window->start + dt_samples,
                window->npts
            );

            // Calculate covariance (no need to rotate back to input frame)
            double cov_xx = gsl_stats_covariance(
                rot_x.vector.data,
                rot_x.vector.stride,
                rot_x.vector.data,
                rot_x.vector.stride,
                rot_x.vector.size
            );
            double cov_xy = gsl_stats_covariance(
                rot_x.vector.data,
                rot_x.vector.stride,
                rot_y.vector.data,
                rot_y.vector.stride,
                rot_x.vector.size
            );
            double cov_yy = gsl_stats_covariance(
                rot_y.vector.data,
                rot_y.vector.stride,
                rot_y.vector.data,
                rot_y.vector.stride,
                rot_y.vector.size
            );

            gsl_matrix_set(C, 0, 0, cov_xx);
            gsl_matrix_set(C, 0, 1, cov_xy);
            gsl_matrix_set(C, 1, 0, cov_xy);
            gsl_matrix_set(C, 1, 1, cov_yy);

            gsl_eigen_symmv(C, eval, evec, eigen_ws);
            gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_ASC);

            if (gsl_vector_get(eval, 0) < current_min) {
                current_min = gsl_vector_get(eval, 0);
                dt_min = dt;
                phi_min = fast;
            }
        }
    }
    printf("The minimum eigenvalue was ");
    printf("%f at phi=%f and dt=%f s\n", current_min, phi_min, dt_min);

    /* Free memory allocated to data arrays. */
    gsl_matrix_free(rotator);
    gsl_matrix_free(rotated_data);
    gsl_eigen_symmv_free(eigen_ws);
    gsl_matrix_free(evec);
    gsl_vector_free(eval);
    gsl_matrix_free(data);
    gsl_matrix_free(C);
}
