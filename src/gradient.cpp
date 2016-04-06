/*
 * gradient.cpp
 *
 * Code generation for function 'gradient'
 *
 * C source code generated on: Tue Apr  1 12:30:12 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Optimal_affine_tracking_3d16_fast_realtime.h"
#include "gradient.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
void gradient(const real_T f[168021], real_T varargout_1[168021], real_T
              varargout_2[168021])
{
  int32_T i1;
  int32_T i2;
  int32_T j;
  int32_T k;
  i1 = -1;
  i2 = 167639;
  for (j = 0; j < 381; j++) {
    i1++;
    i2++;
    varargout_1[i1] = f[i1 + 381] - f[i1];
    for (k = 0; k < 439; k++) {
      varargout_1[i1 + (k + 1) * 381] = (f[i1 + (k + 2) * 381] - f[i1 + k * 381])
        / 2.0;
    }

    varargout_1[i2] = f[i2] - f[i2 - 381];
  }

  i2 = -1;
  for (j = 0; j < 441; j++) {
    i1 = i2 + 1;
    i2 += 381;
    varargout_2[i1] = f[i1 + 1] - f[i1];
    for (k = 0; k < 379; k++) {
      varargout_2[(i1 + k) + 1] = (f[(i1 + k) + 2] - f[i1 + k]) / 2.0;
    }

    varargout_2[i2] = f[i2] - f[i2 - 1];
  }
}

/* End of code generation (gradient.cpp) */
