/*
 * resampling.cpp
 *
 * Code generation for function 'resampling'
 *
 * C source code generated on: Tue Apr  1 12:30:12 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Optimal_affine_tracking_3d16_fast_realtime.h"
#include "resampling.h"
#include "meshgrid.h"
#include "cumulative_sum.h"
#include "rand.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
void resampling(const real_T w[25], real_T outindex[25])
{
  real_T ndbl;
  real_T y;
  int32_T xtmp;
  real_T apnd;
  real_T absa;
  real_T b_absa;
  real_T U_data[26];
  int32_T nm1d2;
  int32_T k;
  real_T dv10[25];
  real_T cum[26];
  real_T dv11[25];
  real_T dv12[25];
  int32_T i;
  real_T jj[625];
  real_T ii[625];
  int32_T i16;
  boolean_T cond[625];
  int32_T ind;
  int32_T b_ind;
  int16_T ii_data[1];
  boolean_T exitg1;
  int16_T tmp_data[625];
  int32_T loop_ub;
  int32_T i17;
  int16_T b_tmp_data[625];
  int16_T b_ii_data[1];
  int16_T j_data[1];
  ndbl = c_rand();
  y = ndbl / 25.0;
  if (rtIsNaN(y)) {
    xtmp = 0;
    y = rtNaN;
    apnd = 1.0;
  } else if (rtIsInf(y)) {
    xtmp = 0;
    y = rtNaN;
    apnd = 1.0;
  } else {
    ndbl = floor((1.0 - y) / 0.04 + 0.5);
    apnd = y + ndbl * 0.04;
    absa = fabs(y);
    if (absa > 1.0) {
      b_absa = absa;
    } else {
      b_absa = 1.0;
    }

    if (fabs(apnd - 1.0) < 4.4408920985006262E-16 * b_absa) {
      ndbl++;
      apnd = 1.0;
    } else if (apnd - 1.0 > 0.0) {
      apnd = y + (ndbl - 1.0) * 0.04;
    } else {
      ndbl++;
    }

    xtmp = (int32_T)ndbl - 1;
  }

  U_data[0] = y;
  if (xtmp + 1 > 1) {
    U_data[xtmp] = apnd;
    nm1d2 = xtmp / 2;
    for (k = 1; k < nm1d2; k++) {
      ndbl = (real_T)k * 0.04;
      U_data[k] = y + ndbl;
      U_data[xtmp - k] = apnd - ndbl;
    }

    if (nm1d2 << 1 == xtmp) {
      U_data[nm1d2] = (y + apnd) / 2.0;
    } else {
      ndbl = (real_T)nm1d2 * 0.04;
      U_data[nm1d2] = y + ndbl;
      U_data[nm1d2 + 1] = apnd - ndbl;
    }
  }

  cumulative_sum(w, dv10);
  cum[0] = 0.0;
  for (i = 0; i < 25; i++) {
    cum[i + 1] = dv10[i];
    outindex[i] = 0.0;
    dv11[i] = 1.0 + (real_T)i;
    dv12[i] = 2.0 + (real_T)i;
  }

  meshgrid(dv11, dv12, ii, jj);
  for (i16 = 0; i16 < 625; i16++) {
    cond[i16] = ((U_data[(int32_T)ii[i16] - 1] > cum[(int32_T)(jj[i16] - 1.0) -
                  1]) && (U_data[(int32_T)ii[i16] - 1] <= cum[(int32_T)jj[i16] -
                          1]));
  }

  i = 0;
  for (ind = 0; ind < 25; ind++) {
    b_ind = 1 + ind * 25;
    if (b_ind > b_ind + 24) {
      b_ind = 1;
      i16 = 0;
    } else {
      i16 = b_ind + 24;
    }

    xtmp = (i16 - b_ind) + 1;
    if (1 <= xtmp) {
      k = 1;
    } else {
      k = xtmp;
    }

    nm1d2 = 0;
    xtmp = (i16 - b_ind) + 1;
    exitg1 = FALSE;
    while ((exitg1 == FALSE) && (xtmp > 0)) {
      loop_ub = i16 - b_ind;
      for (i17 = 0; i17 <= loop_ub; i17++) {
        tmp_data[i17] = (int16_T)(b_ind + i17);
      }

      loop_ub = (i16 - b_ind) + 1;
      for (i17 = 0; i17 < loop_ub; i17++) {
        b_tmp_data[i17] = tmp_data[i17];
      }

      if (cond[b_tmp_data[xtmp - 1] - 1]) {
        nm1d2 = 1;
        ii_data[0] = (int16_T)xtmp;
        exitg1 = TRUE;
      } else {
        xtmp--;
      }
    }

    if (k == 1) {
      if (nm1d2 == 0) {
        k = 0;
      }
    } else {
      if (1 > nm1d2) {
        loop_ub = -1;
      } else {
        loop_ub = 0;
      }

      xtmp = loop_ub + 1;
      i16 = 0;
      while (i16 <= xtmp - 1) {
        i16 = 0;
        while (i16 <= 0) {
          b_ii_data[0] = ii_data[0];
          i16 = 1;
        }

        i16 = 1;
      }

      k = loop_ub + 1;
      xtmp = loop_ub + 1;
      i16 = 0;
      while (i16 <= xtmp - 1) {
        ii_data[0] = b_ii_data[0];
        i16 = 1;
      }

      i16 = loop_ub + 1;
      nm1d2 = i16 / 2;
      xtmp = 1;
      while (xtmp <= nm1d2) {
        xtmp = ii_data[0];
        ii_data[0] = ii_data[loop_ub];
        ii_data[loop_ub] = (int16_T)xtmp;
        xtmp = 2;
      }
    }

    loop_ub = k;
    i16 = 0;
    while (i16 <= loop_ub - 1) {
      j_data[0] = ii_data[0];
      i16 = 1;
    }

    outindex[i] = (real_T)j_data[0];
    i++;
  }
}

/* End of code generation (resampling.cpp) */
