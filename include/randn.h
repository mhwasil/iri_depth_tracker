/*
 * randn.h
 *
 * Code generation for function 'randn'
 *
 * C source code generated on: Tue Apr  1 12:30:12 2014
 *
 */

#ifndef __RANDN_H__
#define __RANDN_H__
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
extern real_T genrandu(uint32_T mt[625]);
extern void randn(real_T r[6]);
#endif
/* End of code generation (randn.h) */
