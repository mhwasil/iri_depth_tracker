/*
 * mldivide.cpp
 *
 * Code generation for function 'mldivide'
 *
 * C source code generated on: Tue Apr  1 12:30:12 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Optimal_affine_tracking_3d16_fast_realtime.h"
#include "mldivide.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
void b_mldivide(const real_T A[4], const real_T B[4], real_T Y[4])
{
  int32_T r1;
  int32_T r2;
  real_T a21;
  real_T a22;
  int32_T k;
  if (fabs(A[1]) > fabs(A[0])) {
    r1 = 1;
    r2 = 0;
  } else {
    r1 = 0;
    r2 = 1;
  }

  a21 = A[r2] / A[r1];
  a22 = A[2 + r2] - a21 * A[2 + r1];
  for (k = 0; k < 2; k++) {
    Y[1 + (k << 1)] = (B[r2 + (k << 1)] - B[r1 + (k << 1)] * a21) / a22;
    Y[k << 1] = (B[r1 + (k << 1)] - Y[1 + (k << 1)] * A[2 + r1]) / A[r1];
  }
}

void mldivide(const real_T A[9], const real_T B[9], real_T Y[9])
{
  real_T b_A[9];
  int32_T r1;
  int32_T r2;
  int32_T r3;
  real_T maxval;
  real_T a21;
  int32_T rtemp;
  memcpy(&b_A[0], (real_T *)&A[0], 9U * sizeof(real_T));
  r1 = 0;
  r2 = 1;
  r3 = 2;
  maxval = fabs(A[0]);
  a21 = fabs(A[1]);
  if (a21 > maxval) {
    maxval = a21;
    r1 = 1;
    r2 = 0;
  }

  if (fabs(A[2]) > maxval) {
    r1 = 2;
    r2 = 1;
    r3 = 0;
  }

  b_A[r2] = A[r2] / A[r1];
  b_A[r3] /= b_A[r1];
  b_A[3 + r2] -= b_A[r2] * b_A[3 + r1];
  b_A[3 + r3] -= b_A[r3] * b_A[3 + r1];
  b_A[6 + r2] -= b_A[r2] * b_A[6 + r1];
  b_A[6 + r3] -= b_A[r3] * b_A[6 + r1];
  if (fabs(b_A[3 + r3]) > fabs(b_A[3 + r2])) {
    rtemp = r2;
    r2 = r3;
    r3 = rtemp;
  }

  b_A[3 + r3] /= b_A[3 + r2];
  b_A[6 + r3] -= b_A[3 + r3] * b_A[6 + r2];
  for (rtemp = 0; rtemp < 3; rtemp++) {
    Y[3 * rtemp] = B[r1 + 3 * rtemp];
    Y[1 + 3 * rtemp] = B[r2 + 3 * rtemp] - Y[3 * rtemp] * b_A[r2];
    Y[2 + 3 * rtemp] = (B[r3 + 3 * rtemp] - Y[3 * rtemp] * b_A[r3]) - Y[1 + 3 *
      rtemp] * b_A[3 + r3];
    Y[2 + 3 * rtemp] /= b_A[6 + r3];
    Y[3 * rtemp] -= Y[2 + 3 * rtemp] * b_A[6 + r1];
    Y[1 + 3 * rtemp] -= Y[2 + 3 * rtemp] * b_A[6 + r2];
    Y[1 + 3 * rtemp] /= b_A[3 + r2];
    Y[3 * rtemp] -= Y[1 + 3 * rtemp] * b_A[3 + r1];
    Y[3 * rtemp] /= b_A[r1];
  }
}

/* End of code generation (mldivide.cpp) */
