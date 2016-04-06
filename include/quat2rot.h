/*
 * quat2rot.h
 *
 * Code generation for function 'quat2rot'
 *
 * C source code generated on: Tue Apr  1 12:30:12 2014
 *
 */

#ifndef __QUAT2ROT_H__
#define __QUAT2ROT_H__
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
extern void quat2rot(const real_T Q[4], real_T R[9]);
#endif
/* End of code generation (quat2rot.h) */
