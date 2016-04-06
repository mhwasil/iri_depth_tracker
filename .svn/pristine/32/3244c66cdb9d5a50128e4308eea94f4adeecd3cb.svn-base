/*
 * estimateRigidTransform.cpp
 *
 * Code generation for function 'estimateRigidTransform'
 *
 * C source code generated on: Tue Apr  1 12:30:12 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Optimal_affine_tracking_3d16_fast_realtime.h"
#include "estimateRigidTransform.h"
#include "quat2rot.h"
#include "eye.h"
#include "svd.h"
#include "crossTimesMatrix.h"
#include "sum.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
void estimateRigidTransform(const real_T x[3888], const real_T y[3888], real_T
  T[16])
{
  static real_T b_x[3888];
  static real_T b_y[3888];
  int32_T i12;
  int32_T i13;
  real_T x_centroid[3];
  real_T y_centroid[3];
  real_T x_centrized[3888];
  real_T y_centrized[3888];
  real_T b_y_centrized[3888];
  static real_T R22[11664];
  real_T B[16];
  static real_T A[20736];
  int32_T ii;
  real_T d2;
  int32_T i14;
  real_T b_B[16];
  real_T V[16];
  real_T unusedU1[16];
  real_T dv6[9];
  real_T dv7[9];
  real_T dv8[9];
  static const int8_T iv6[4] = { 0, 0, 0, 1 };

  real_T dv9[16];

  /*  ESTIMATERIGIDTRANSFORM */
  /*    [T, EPS] = ESTIMATERIGIDTRANSFORM(X, Y) estimates the rigid transformation */
  /*    that best aligns x with y (in the least-squares sense). */
  /*    */
  /*    Reference: "Estimating Rigid Transformations" in  */
  /*    "Computer Vision, a modern approach" by Forsyth and Ponce (1993), page 480 */
  /*    (page 717(?) of the newer edition) */
  /*  */
  /*    Input: */
  /*        X: 3xN, N 3-D points (N>=3) */
  /*        Y: 3xN, N 3-D points (N>=3) */
  /*  */
  /*    Output */
  /*        T: the rigid transformation that aligns x and y as:  xh = T * yh */
  /*           (h denotes homogenous coordinates)   */
  /*           (corrspondence between points x(:,i) and y(:,i) is assumed) */
  /*         */
  /*        EPS: the smallest singular value. The closer this value it is  */
  /*           to 0, the better the estimate is. (large values mean that the  */
  /*           transform between the two data sets cannot be approximated */
  /*           well with a rigid transform. */
  /*  */
  /*    Babak Taati, 2003 */
  /*    (revised 2009) */
  for (i12 = 0; i12 < 1296; i12++) {
    for (i13 = 0; i13 < 3; i13++) {
      b_x[i13 + 3 * i12] = x[i12 + 1296 * i13];
      b_y[i13 + 3 * i12] = y[i12 + 1296 * i13];
    }
  }

  /*  if nargin ~= 2 */
  /*      error('Requires two input arguments.') */
  /*  end */
  /*   */
  /*  if size(x,1)~=3 || size(y,1)~=3 */
  /*      error('Input point clouds must be a 3xN matrix.'); */
  /*  end */
  /*   */
  /*  if size(x, 2) ~= size(y,2) */
  /*      error('Input point clouds must be of the same size'); */
  /*  end                             */
  /*   */
  /*  if size(x,2)<3 || size(y,2)<3 */
  /*      error('At least 3 point matches are needed'); */
  /*  end                             */
  /*  since x has N=3+ points, length shows the number of points */
  sum(b_x, x_centroid);
  sum(b_y, y_centroid);
  for (i12 = 0; i12 < 3; i12++) {
    y_centroid[i12] /= 1296.0;
    x_centroid[i12] /= 1296.0;
  }

  for (i12 = 0; i12 < 1296; i12++) {
    x_centrized[3 * i12] = b_x[3 * i12] - x_centroid[0];
    x_centrized[1 + 3 * i12] = b_x[1 + 3 * i12] - x_centroid[1];
    x_centrized[2 + 3 * i12] = b_x[2 + 3 * i12] - x_centroid[2];
    y_centrized[3 * i12] = b_y[3 * i12] - y_centroid[0];
    y_centrized[1 + 3 * i12] = b_y[1 + 3 * i12] - y_centroid[1];
    y_centrized[2 + 3 * i12] = b_y[2 + 3 * i12] - y_centroid[2];
  }

  for (i12 = 0; i12 < 3; i12++) {
    for (i13 = 0; i13 < 1296; i13++) {
      b_y[i13 + 1296 * i12] = y_centrized[i12 + 3 * i13] - x_centrized[i12 + 3 *
        i13];
    }
  }

  for (i12 = 0; i12 < 3888; i12++) {
    b_x[i12] = x_centrized[i12] - y_centrized[i12];
    b_y_centrized[i12] = y_centrized[i12] + x_centrized[i12];
  }

  crossTimesMatrix(b_y_centrized, R22);
  memset(&B[0], 0, sizeof(real_T) << 4);
  memset(&A[0], 0, 20736U * sizeof(real_T));
  for (ii = 0; ii < 1296; ii++) {
    A[ii << 4] = 0.0;
    for (i12 = 0; i12 < 3; i12++) {
      A[((i12 + 1) << 2) + (ii << 4)] = b_y[ii + 1296 * i12];
    }

    for (i12 = 0; i12 < 3; i12++) {
      A[(i12 + (ii << 4)) + 1] = b_x[i12 + 3 * ii];
    }

    for (i12 = 0; i12 < 3; i12++) {
      for (i13 = 0; i13 < 3; i13++) {
        A[((i13 + ((i12 + 1) << 2)) + (ii << 4)) + 1] = R22[(i13 + 3 * i12) + 9 *
          ii];
      }
    }

    for (i12 = 0; i12 < 4; i12++) {
      for (i13 = 0; i13 < 4; i13++) {
        d2 = 0.0;
        for (i14 = 0; i14 < 4; i14++) {
          d2 += A[(i14 + (i12 << 2)) + (ii << 4)] * A[(i14 + (i13 << 2)) + (ii <<
            4)];
        }

        B[i12 + (i13 << 2)] += d2;
      }
    }
  }

  memcpy(&b_B[0], &B[0], sizeof(real_T) << 4);
  svd(b_B, B, unusedU1, V);
  b_eye(dv6);
  b_eye(dv7);
  quat2rot(*(real_T (*)[4])&V[12], dv8);
  for (i12 = 0; i12 < 3; i12++) {
    for (i13 = 0; i13 < 3; i13++) {
      unusedU1[i13 + (i12 << 2)] = dv7[i13 + 3 * i12];
    }
  }

  for (i12 = 0; i12 < 3; i12++) {
    unusedU1[12 + i12] = x_centroid[i12];
  }

  for (i12 = 0; i12 < 4; i12++) {
    unusedU1[3 + (i12 << 2)] = (real_T)iv6[i12];
  }

  for (i12 = 0; i12 < 3; i12++) {
    for (i13 = 0; i13 < 3; i13++) {
      V[i13 + (i12 << 2)] = dv8[i13 + 3 * i12];
    }
  }

  for (i12 = 0; i12 < 3; i12++) {
    V[12 + i12] = 0.0;
  }

  for (i12 = 0; i12 < 4; i12++) {
    V[3 + (i12 << 2)] = (real_T)iv6[i12];
  }

  for (i12 = 0; i12 < 4; i12++) {
    for (i13 = 0; i13 < 4; i13++) {
      dv9[i12 + (i13 << 2)] = 0.0;
      for (i14 = 0; i14 < 4; i14++) {
        dv9[i12 + (i13 << 2)] += unusedU1[i12 + (i14 << 2)] * V[i14 + (i13 << 2)];
      }
    }
  }

  for (i12 = 0; i12 < 3; i12++) {
    for (i13 = 0; i13 < 3; i13++) {
      unusedU1[i13 + (i12 << 2)] = dv6[i13 + 3 * i12];
    }
  }

  for (i12 = 0; i12 < 3; i12++) {
    unusedU1[12 + i12] = -y_centroid[i12];
  }

  for (i12 = 0; i12 < 4; i12++) {
    unusedU1[3 + (i12 << 2)] = (real_T)iv6[i12];
  }

  for (i12 = 0; i12 < 4; i12++) {
    for (i13 = 0; i13 < 4; i13++) {
      T[i12 + (i13 << 2)] = 0.0;
      for (i14 = 0; i14 < 4; i14++) {
        T[i12 + (i13 << 2)] += dv9[i12 + (i14 << 2)] * unusedU1[i14 + (i13 << 2)];
      }
    }
  }

  /*  Eps = S(4,4); */
}

/* End of code generation (estimateRigidTransform.cpp) */
