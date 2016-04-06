/*
 * repmat.cpp
 *
 * Code generation for function 'repmat'
 *
 * C source code generated on: Tue Apr  1 12:30:12 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Optimal_affine_tracking_3d16_fast_realtime.h"
#include "repmat.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
static int32_T div_s32(int32_T numerator, int32_T denominator);

/* Function Definitions */
static int32_T div_s32(int32_T numerator, int32_T denominator)
{
  int32_T quotient;
  uint32_T absNumerator;
  uint32_T absDenominator;
  int32_T quotientNeedsNegation;
  if (denominator == 0) {
    if (numerator >= 0) {
      quotient = MAX_int32_T;
    } else {
      quotient = MIN_int32_T;
    }
  } else {
    if (numerator >= 0) {
      absNumerator = (uint32_T)numerator;
    } else {
      absNumerator = (uint32_T)-numerator;
    }

    if (denominator >= 0) {
      absDenominator = (uint32_T)denominator;
    } else {
      absDenominator = (uint32_T)-denominator;
    }

    quotientNeedsNegation = (int32_T)((int32_T)(numerator < 0) != (int32_T)
      (denominator < 0));
    absNumerator /= absDenominator;
    if ((uint32_T)quotientNeedsNegation) {
      quotient = -(int32_T)absNumerator;
    } else {
      quotient = (int32_T)absNumerator;
    }
  }

  return quotient;
}

void b_repmat(const real_T a[168021], real_T b[336042])
{
  int32_T ib;
  int32_T jtilecol;
  int32_T iacol;
  int32_T jcol;
  int32_T k;
  ib = 0;
  for (jtilecol = 0; jtilecol < 2; jtilecol++) {
    iacol = 1;
    for (jcol = 0; jcol < 441; jcol++) {
      for (k = 0; k < 381; k++) {
        b[ib] = a[iacol - 1];
        iacol++;
        ib++;
      }
    }
  }
}

void repmat(const real_T a[504063], real_T b[1008126])
{
  int32_T db[3];
  int32_T da[3];
  int32_T ibtmp;
  int32_T k;
  static const int16_T iv2[3] = { 381, 441, 3 };

  static const int16_T iv3[3] = { 381, 882, 3 };

  int32_T ib;
  int32_T ia;
  int32_T r;
  int32_T u1;
  for (ibtmp = 0; ibtmp < 3; ibtmp++) {
    db[ibtmp] = 1;
    da[ibtmp] = 1;
  }

  for (k = 0; k < 2; k++) {
    da[k + 1] = da[k] * iv2[k];
    db[k + 1] = db[k] * iv3[k];
  }

  for (ib = 0; ib < 1008126; ib++) {
    ia = 0;
    ibtmp = ib;
    for (k = 2; k > -1; k += -1) {
      r = ibtmp - div_s32(ibtmp, db[k]) * db[k];
      ibtmp = div_s32(ibtmp - r, db[k]);
      u1 = iv2[k];
      ibtmp -= u1 * div_s32(ibtmp, u1);
      ia += da[k] * ibtmp;
      ibtmp = r;
    }

    b[ib] = a[ia];
  }
}

/* End of code generation (repmat.cpp) */
