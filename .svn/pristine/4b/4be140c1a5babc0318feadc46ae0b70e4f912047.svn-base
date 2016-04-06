/*
 * cumulative_sum.cpp
 *
 * Code generation for function 'cumulative_sum'
 *
 * C source code generated on: Tue Apr  1 12:30:12 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Optimal_affine_tracking_3d16_fast_realtime.h"
#include "cumulative_sum.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
void cumulative_sum(const real_T x[25], real_T csum[25])
{
  int32_T i;
  memset(&csum[0], 0, 25U * sizeof(real_T));
  csum[0] = x[0];
  for (i = 0; i < 24; i++) {
    csum[1 + i] = x[1 + i] + csum[i];
  }
}

/* End of code generation (cumulative_sum.cpp) */
