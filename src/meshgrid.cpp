/*
 * meshgrid.cpp
 *
 * Code generation for function 'meshgrid'
 *
 * C source code generated on: Tue Apr  1 12:30:12 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Optimal_affine_tracking_3d16_fast_realtime.h"
#include "meshgrid.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
void meshgrid(const real_T x[25], const real_T y[25], real_T xx[625], real_T yy
              [625])
{
  int32_T ia;
  int32_T ib;
  int32_T iacol;
  int32_T jcol;
  int32_T itilerow;
  ia = 1;
  ib = 0;
  iacol = 1;
  for (jcol = 0; jcol < 25; jcol++) {
    for (itilerow = 0; itilerow < 25; itilerow++) {
      xx[ib] = x[iacol - 1];
      ia = iacol + 1;
      ib++;
    }

    iacol = ia;
  }

  ib = 0;
  for (iacol = 0; iacol < 25; iacol++) {
    ia = 0;
    for (jcol = 0; jcol < 25; jcol++) {
      yy[ib] = y[ia];
      ia++;
      ib++;
    }
  }
}

/* End of code generation (meshgrid.cpp) */
