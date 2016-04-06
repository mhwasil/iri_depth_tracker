/*
 * compute_prob1.h
 *
 * Code generation for function 'compute_prob1'
 *
 * C source code generated on: Tue Apr  1 12:30:12 2014
 *
 */

#ifndef __COMPUTE_PROB1_H__
#define __COMPUTE_PROB1_H__
/* Include files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_defines.h"
#include "rt_nonfinite.h"

#include "rtwtypes.h"
#include "Optimal_affine_tracking_3d16_fast_realtime_types.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern void compute_prob1(real_T X_par[2250000], real_T X_par_pred[2250000], real_T AR_velocity[2250000], const real_T dw_dp[15552], real_T t, real_T center_x, real_T center_y, const real_T Ixyz[504063], const real_T point_matrix[3888], const real_T mean_img[3888], real_T tracked_images[38880000], real_T Aff_matrix[9], real_T centroid[3]);
extern void eml_li_find(const boolean_T x[1296], int32_T y_data[1296], int32_T y_size[1]);
#endif
/* End of code generation (compute_prob1.h) */
