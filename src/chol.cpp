/*
 * chol.cpp
 *
 * Code generation for function 'chol'
 *
 * C source code generated on: Tue Apr  1 12:30:12 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Optimal_affine_tracking_3d16_fast_realtime.h"
#include "chol.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
void chol(real_T A[36])
{
  int32_T info;
  int32_T colj;
  int32_T j;
  boolean_T exitg1;
  int32_T jj;
  real_T ajj;
  int32_T ix;
  int32_T iy;
  int32_T jmax;
  int32_T i24;
  real_T c;
  int32_T i;
  int32_T ia;
  info = 0;
  colj = 1;
  j = 1;
  exitg1 = FALSE;
  while ((exitg1 == FALSE) && (j < 7)) {
    jj = (colj + j) - 2;
    ajj = 0.0;
    if (j - 1 < 1) {
    } else {
      ix = colj;
      iy = colj;
      for (jmax = 1; jmax < j; jmax++) {
        ajj += A[ix - 1] * A[iy - 1];
        ix++;
        iy++;
      }
    }

    ajj = A[jj] - ajj;
    if (ajj > 0.0) {
      ajj = sqrt(ajj);
      A[jj] = ajj;
      if (j < 6) {
        if (j - 1 == 0) {
        } else {
          iy = jj + 6;
          i24 = (colj + 6 * (5 - j)) + 6;
          for (jmax = colj + 6; jmax <= i24; jmax += 6) {
            ix = colj;
            c = 0.0;
            i = (jmax + j) - 2;
            for (ia = jmax; ia <= i; ia++) {
              c += A[ia - 1] * A[ix - 1];
              ix++;
            }

            A[iy] += -c;
            iy += 6;
          }
        }

        ajj = 1.0 / ajj;
        i24 = jj + 6 * (5 - j);
        for (jmax = jj + 6; jmax + 1 <= i24 + 7; jmax += 6) {
          A[jmax] *= ajj;
        }

        colj += 6;
      }

      j++;
    } else {
      A[jj] = ajj;
      info = j;
      exitg1 = TRUE;
    }
  }

  if (info == 0) {
    jmax = 6;
  } else {
    jmax = info - 1;
  }

  for (j = 0; j + 1 <= jmax; j++) {
    for (i = j + 1; i + 1 <= jmax; i++) {
      A[i + 6 * j] = 0.0;
    }
  }
}

/* End of code generation (chol.cpp) */
