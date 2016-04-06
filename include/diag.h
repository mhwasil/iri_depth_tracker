/*
 * diag.h
 *
 * Code generation for function 'diag'
 *
 * C source code generated on: Tue Apr  1 12:30:12 2014
 *
 */

#ifndef __DIAG_H__
#define __DIAG_H__
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
extern void b_diag(const creal_T v[2], creal_T d[4]);
extern void diag(const creal_T v[3], creal_T d[9]);
#endif
/* End of code generation (diag.h) */
