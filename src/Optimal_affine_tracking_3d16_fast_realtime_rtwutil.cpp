/*
 * Optimal_affine_tracking_3d16_fast_realtime_rtwutil.cpp
 *
 * Code generation for function 'Optimal_affine_tracking_3d16_fast_realtime_rtwutil'
 *
 * C source code generated on: Tue Apr  1 12:30:12 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Optimal_affine_tracking_3d16_fast_realtime.h"
#include "Optimal_affine_tracking_3d16_fast_realtime_rtwutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
real_T rt_hypotd_snf(real_T u0, real_T u1)
{
  real_T y;
  real_T a;
  real_T b;
  a = fabs(u0);
  b = fabs(u1);
  if (a < b) {
    a /= b;
    y = b * sqrt(a * a + 1.0);
  } else if (a > b) {
    b /= a;
    y = a * sqrt(b * b + 1.0);
  } else if (rtIsNaN(b)) {
    y = b;
  } else {
    y = a * 1.4142135623730951;
  }

  return y;
}

/* End of code generation (Optimal_affine_tracking_3d16_fast_realtime_rtwutil.cpp) */
