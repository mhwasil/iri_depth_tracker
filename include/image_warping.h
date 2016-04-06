/*
 * image_warping.h
 *
 * Code generation for function 'image_warping'
 *
 * C source code generated on: Tue Apr  1 12:30:12 2014
 *
 */

#ifndef __IMAGE_WARPING_H__
#define __IMAGE_WARPING_H__
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
extern void b_image_warping(const real_T frame[336042], const real_T aff_matrix[9], const real_T point_matrix[3888], real_T warped_img[1296]);
extern void image_warping(const real_T frame[168021], const real_T aff_matrix[9], const real_T point_matrix[3888], real_T warped_img[1296]);
#endif
/* End of code generation (image_warping.h) */
