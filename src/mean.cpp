/*
 * mean.cpp
 *
 * Code generation for function 'mean'
 *
 * C source code generated on: Tue Apr  1 12:30:12 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Optimal_affine_tracking_3d16_fast_realtime.h"
#include "mean.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
void b_mean(const real_T x[50], real_T y[2])
{
  int32_T iy;
  int32_T ixstart;
  int32_T j;
  int32_T ix;
  real_T s;
  int32_T k;
  iy = -1;
  ixstart = -1;
  for (j = 0; j < 2; j++) {
    ixstart++;
    ix = ixstart;
    s = x[ixstart];
    for (k = 0; k < 24; k++) {
      ix += 2;
      s += x[ix];
    }

    iy++;
    y[iy] = s;
  }

  for (iy = 0; iy < 2; iy++) {
    y[iy] /= 25.0;
  }
}

void c_mean(const real_T x[3888], real_T y[3])
{
  int32_T ix;
  int32_T iy;
  int32_T i;
  int32_T ixstart;
  real_T s;
  ix = -1;
  iy = -1;
  for (i = 0; i < 3; i++) {
    ixstart = ix + 1;
    ix++;
    s = x[ixstart];
    for (ixstart = 0; ixstart < 1295; ixstart++) {
      ix++;
      s += x[ix];
    }

    iy++;
    y[iy] = s;
  }

  for (ixstart = 0; ixstart < 3; ixstart++) {
    y[ixstart] /= 1296.0;
  }
}

real_T mean(const real_T x_data[1296], const int32_T x_size[1])
{
  real_T y;
  int32_T k;
  if (x_size[0] == 0) {
    y = 0.0;
  } else {
    y = x_data[0];
    for (k = 2; k <= x_size[0]; k++) {
      y += x_data[k - 1];
    }
  }

  y /= (real_T)x_size[0];
  return y;
}

/* End of code generation (mean.cpp) */
