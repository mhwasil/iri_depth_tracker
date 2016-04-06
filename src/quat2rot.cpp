/*
 * quat2rot.cpp
 *
 * Code generation for function 'quat2rot'
 *
 * C source code generated on: Tue Apr  1 12:30:12 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Optimal_affine_tracking_3d16_fast_realtime.h"
#include "quat2rot.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
void quat2rot(const real_T Q[4], real_T R[9])
{
  /*  QUAT2ROT */
  /*    R = QUAT2ROT(Q) converts a quaternion (4x1 or 1x4) into a 3x3 rotation mattrix */
  /*  */
  /*    reference: google! */
  /*  */
  /*    Babak Taati, 2003 */
  /*    (revised 2009) */
  R[0] = ((Q[0] * Q[0] + Q[1] * Q[1]) - Q[2] * Q[2]) - Q[3] * Q[3];
  R[3] = 2.0 * (Q[1] * Q[2] - Q[0] * Q[3]);
  R[6] = 2.0 * (Q[1] * Q[3] + Q[0] * Q[2]);
  R[1] = 2.0 * (Q[1] * Q[2] + Q[0] * Q[3]);
  R[4] = ((Q[0] * Q[0] - Q[1] * Q[1]) + Q[2] * Q[2]) - Q[3] * Q[3];
  R[7] = 2.0 * (Q[2] * Q[3] - Q[0] * Q[1]);
  R[2] = 2.0 * (Q[1] * Q[3] - Q[0] * Q[2]);
  R[5] = 2.0 * (Q[2] * Q[3] + Q[0] * Q[1]);
  R[8] = ((Q[0] * Q[0] - Q[1] * Q[1]) - Q[2] * Q[2]) + Q[3] * Q[3];
}

/* End of code generation (quat2rot.cpp) */
