/*
 * diag.cpp
 *
 * Code generation for function 'diag'
 *
 * C source code generated on: Tue Apr  1 12:30:12 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Optimal_affine_tracking_3d16_fast_realtime.h"
#include "diag.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
void b_diag(const creal_T v[2], creal_T d[4])
{
  int32_T j;
  for (j = 0; j < 4; j++) {
    d[j].re = 0.0;
    d[j].im = 0.0;
  }

  for (j = 0; j < 2; j++) {
    d[j + (j << 1)] = v[j];
  }
}

void diag(const creal_T v[3], creal_T d[9])
{
  int32_T j;
  for (j = 0; j < 9; j++) {
    d[j].re = 0.0;
    d[j].im = 0.0;
  }

  for (j = 0; j < 3; j++) {
    d[j + 3 * j] = v[j];
  }
}

/* End of code generation (diag.cpp) */
