/*
 * meshgrid.h
 *
 * Code generation for function 'meshgrid'
 *
 * C source code generated on: Tue Apr  1 12:30:12 2014
 *
 */

#ifndef __MESHGRID_H__
#define __MESHGRID_H__
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
extern void meshgrid(const real_T x[25], const real_T y[25], real_T xx[625], real_T yy[625]);
#endif
/* End of code generation (meshgrid.h) */
