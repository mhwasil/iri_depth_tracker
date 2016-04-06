/*
 * crossTimesMatrix.cpp
 *
 * Code generation for function 'crossTimesMatrix'
 *
 * C source code generated on: Tue Apr  1 12:30:12 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Optimal_affine_tracking_3d16_fast_realtime.h"
#include "crossTimesMatrix.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
void crossTimesMatrix(const real_T V[3888], real_T V_times[11664])
{
  int32_T i15;

  /*  CROSSTIMESMATRIX */
  /*    V_TIMES = CROSSTIMESMATRIX(V) returns a 3x3 (or a series of 3x3) cross times matrices of input vector(s) V */
  /*   */
  /*    Input: */
  /*        V a 3xN matrix, rpresenting a series of 3x1 vectors */
  /*   */
  /*    Output:    */
  /*        V_TIMES (Vx) a series of 3x3 matrices where V_times(:,:,i) is the Vx matrix for the vector V(:,i) */
  /*   */
  /*  	Babak Taati, 2003 */
  /*    (revised 2009) */
  memset(&V_times[0], 0, 11664U * sizeof(real_T));

  /*  V_times(1,1,:) = 0; */
  for (i15 = 0; i15 < 1296; i15++) {
    V_times[3 + 9 * i15] = -V[2 + 3 * i15];
    V_times[6 + 9 * i15] = V[1 + 3 * i15];
    V_times[1 + 9 * i15] = V[2 + 3 * i15];

    /*  V_times(2,2,:) = 0; */
    V_times[7 + 9 * i15] = -V[3 * i15];
    V_times[2 + 9 * i15] = -V[1 + 3 * i15];
    V_times[5 + 9 * i15] = V[3 * i15];
  }

  /*  V_times(3,3,:) = 0; */
}

/* End of code generation (crossTimesMatrix.cpp) */
