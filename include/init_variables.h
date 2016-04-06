/*
 * init_variables.h
 *
 * Code generation for function 'init_variables'
 *
 * C source code generated on: Tue Apr  1 12:30:12 2014
 *
 */

#ifndef __INIT_VARIABLES_H__
#define __INIT_VARIABLES_H__
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
extern void init_variables(const real_T Ixyz[504063], const real_T ptx[4], const real_T pty[4], real_T X_par_pred[2250000], real_T tracked_images[38880000], real_T dw_dp[15552], real_T X_par[2250000], real_T AR_velocity[2250000], real_T point_matrix[3888], real_T mean_img[3888], real_T corner_p[12], real_T *center_x, real_T *center_y);
#endif
/* End of code generation (init_variables.h) */
