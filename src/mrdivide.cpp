/*
 * mrdivide.cpp
 *
 * Code generation for function 'mrdivide'
 *
 * C source code generated on: Tue Apr  1 12:30:12 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Optimal_affine_tracking_3d16_fast_realtime.h"
#include "mrdivide.h"
#include "log.h"
#include "Optimal_affine_tracking_3d16_fast_realtime_rtwutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
static real_T b_eml_xnrm2(int32_T n, const real_T x[12], int32_T ix0);
static real_T eml_xnrm2(int32_T n, const real_T x[12], int32_T ix0);

/* Function Definitions */
static real_T b_eml_xnrm2(int32_T n, const real_T x[12], int32_T ix0)
{
  real_T y;
  real_T scale;
  int32_T kend;
  int32_T k;
  real_T absxk;
  real_T t;
  y = 0.0;
  if (n == 1) {
    y = fabs(x[ix0 - 1]);
  } else {
    scale = 2.2250738585072014E-308;
    kend = (ix0 + n) - 1;
    for (k = ix0; k <= kend; k++) {
      absxk = fabs(x[k - 1]);
      if (absxk > scale) {
        t = scale / absxk;
        y = 1.0 + y * t * t;
        scale = absxk;
      } else {
        t = absxk / scale;
        y += t * t;
      }
    }

    y = scale * sqrt(y);
  }

  return y;
}

static real_T eml_xnrm2(int32_T n, const real_T x[12], int32_T ix0)
{
  real_T y;
  real_T scale;
  int32_T kend;
  int32_T k;
  real_T absxk;
  real_T t;
  y = 0.0;
  if (n == 1) {
    y = fabs(x[ix0 - 1]);
  } else {
    scale = 2.2250738585072014E-308;
    kend = (ix0 + n) - 1;
    for (k = ix0; k <= kend; k++) {
      absxk = fabs(x[k - 1]);
      if (absxk > scale) {
        t = scale / absxk;
        y = 1.0 + y * t * t;
        scale = absxk;
      } else {
        t = absxk / scale;
        y += t * t;
      }
    }

    y = scale * sqrt(y);
  }

  return y;
}

void b_mrdivide(const creal_T A[9], const creal_T B[9], creal_T y[9])
{
  creal_T b_A[9];
  creal_T b_B[9];
  int32_T rtemp;
  int32_T k;
  int32_T r1;
  int32_T r2;
  int32_T r3;
  real_T maxval;
  real_T a21;
  creal_T Y[9];
  for (rtemp = 0; rtemp < 3; rtemp++) {
    for (k = 0; k < 3; k++) {
      b_A[k + 3 * rtemp].re = B[rtemp + 3 * k].re;
      b_A[k + 3 * rtemp].im = -B[rtemp + 3 * k].im;
      b_B[k + 3 * rtemp].re = A[rtemp + 3 * k].re;
      b_B[k + 3 * rtemp].im = -A[rtemp + 3 * k].im;
    }
  }

  r1 = 0;
  r2 = 1;
  r3 = 2;
  maxval = fabs(b_A[0].re) + fabs(b_A[0].im);
  a21 = fabs(b_A[1].re) + fabs(b_A[1].im);
  if (a21 > maxval) {
    maxval = a21;
    r1 = 1;
    r2 = 0;
  }

  if (fabs(b_A[2].re) + fabs(b_A[2].im) > maxval) {
    r1 = 2;
    r2 = 1;
    r3 = 0;
  }

  b_A[r2] = c_eml_div(b_A[r2], b_A[r1]);
  b_A[r3] = c_eml_div(b_A[r3], b_A[r1]);
  maxval = b_A[r2].re * b_A[3 + r1].im + b_A[r2].im * b_A[3 + r1].re;
  b_A[3 + r2].re -= b_A[r2].re * b_A[3 + r1].re - b_A[r2].im * b_A[3 + r1].im;
  b_A[3 + r2].im -= maxval;
  maxval = b_A[r3].re * b_A[3 + r1].im + b_A[r3].im * b_A[3 + r1].re;
  b_A[3 + r3].re -= b_A[r3].re * b_A[3 + r1].re - b_A[r3].im * b_A[3 + r1].im;
  b_A[3 + r3].im -= maxval;
  maxval = b_A[r2].re * b_A[6 + r1].im + b_A[r2].im * b_A[6 + r1].re;
  b_A[6 + r2].re -= b_A[r2].re * b_A[6 + r1].re - b_A[r2].im * b_A[6 + r1].im;
  b_A[6 + r2].im -= maxval;
  maxval = b_A[r3].re * b_A[6 + r1].im + b_A[r3].im * b_A[6 + r1].re;
  b_A[6 + r3].re -= b_A[r3].re * b_A[6 + r1].re - b_A[r3].im * b_A[6 + r1].im;
  b_A[6 + r3].im -= maxval;
  if (fabs(b_A[3 + r3].re) + fabs(b_A[3 + r3].im) > fabs(b_A[3 + r2].re) + fabs
      (b_A[3 + r2].im)) {
    rtemp = r2;
    r2 = r3;
    r3 = rtemp;
  }

  b_A[3 + r3] = c_eml_div(b_A[3 + r3], b_A[3 + r2]);
  maxval = b_A[3 + r3].re * b_A[6 + r2].im + b_A[3 + r3].im * b_A[6 + r2].re;
  b_A[6 + r3].re -= b_A[3 + r3].re * b_A[6 + r2].re - b_A[3 + r3].im * b_A[6 +
    r2].im;
  b_A[6 + r3].im -= maxval;
  for (k = 0; k < 3; k++) {
    Y[3 * k] = b_B[r1 + 3 * k];
    maxval = Y[3 * k].re * b_A[r2].im + Y[3 * k].im * b_A[r2].re;
    Y[1 + 3 * k].re = b_B[r2 + 3 * k].re - (Y[3 * k].re * b_A[r2].re - Y[3 * k].
      im * b_A[r2].im);
    Y[1 + 3 * k].im = b_B[r2 + 3 * k].im - maxval;
    maxval = Y[3 * k].re * b_A[r3].im + Y[3 * k].im * b_A[r3].re;
    a21 = Y[1 + 3 * k].re * b_A[3 + r3].im + Y[1 + 3 * k].im * b_A[3 + r3].re;
    Y[2 + 3 * k].re = (b_B[r3 + 3 * k].re - (Y[3 * k].re * b_A[r3].re - Y[3 * k]
      .im * b_A[r3].im)) - (Y[1 + 3 * k].re * b_A[3 + r3].re - Y[1 + 3 * k].im *
      b_A[3 + r3].im);
    Y[2 + 3 * k].im = (b_B[r3 + 3 * k].im - maxval) - a21;
    Y[2 + 3 * k] = c_eml_div(Y[2 + 3 * k], b_A[6 + r3]);
    maxval = Y[2 + 3 * k].re * b_A[6 + r1].im + Y[2 + 3 * k].im * b_A[6 + r1].re;
    Y[3 * k].re -= Y[2 + 3 * k].re * b_A[6 + r1].re - Y[2 + 3 * k].im * b_A[6 +
      r1].im;
    Y[3 * k].im -= maxval;
    maxval = Y[2 + 3 * k].re * b_A[6 + r2].im + Y[2 + 3 * k].im * b_A[6 + r2].re;
    Y[1 + 3 * k].re -= Y[2 + 3 * k].re * b_A[6 + r2].re - Y[2 + 3 * k].im * b_A
      [6 + r2].im;
    Y[1 + 3 * k].im -= maxval;
    Y[1 + 3 * k] = c_eml_div(Y[1 + 3 * k], b_A[3 + r2]);
    maxval = Y[1 + 3 * k].re * b_A[3 + r1].im + Y[1 + 3 * k].im * b_A[3 + r1].re;
    Y[3 * k].re -= Y[1 + 3 * k].re * b_A[3 + r1].re - Y[1 + 3 * k].im * b_A[3 +
      r1].im;
    Y[3 * k].im -= maxval;
    Y[3 * k] = c_eml_div(Y[3 * k], b_A[r1]);
  }

  for (rtemp = 0; rtemp < 3; rtemp++) {
    for (k = 0; k < 3; k++) {
      y[k + 3 * rtemp].re = Y[rtemp + 3 * k].re;
      y[k + 3 * rtemp].im = -Y[rtemp + 3 * k].im;
    }
  }
}

creal_T c_eml_div(const creal_T x, const creal_T y)
{
  creal_T z;
  real_T brm;
  real_T bim;
  real_T d;
  if (y.im == 0.0) {
    if (x.im == 0.0) {
      z.re = x.re / y.re;
      z.im = 0.0;
    } else if (x.re == 0.0) {
      z.re = 0.0;
      z.im = x.im / y.re;
    } else {
      z.re = x.re / y.re;
      z.im = x.im / y.re;
    }
  } else if (y.re == 0.0) {
    if (x.re == 0.0) {
      z.re = x.im / y.im;
      z.im = 0.0;
    } else if (x.im == 0.0) {
      z.re = 0.0;
      z.im = -(x.re / y.im);
    } else {
      z.re = x.im / y.im;
      z.im = -(x.re / y.im);
    }
  } else {
    brm = fabs(y.re);
    bim = fabs(y.im);
    if (brm > bim) {
      bim = y.im / y.re;
      d = y.re + bim * y.im;
      z.re = (x.re + bim * x.im) / d;
      z.im = (x.im - bim * x.re) / d;
    } else if (bim == brm) {
      if (y.re > 0.0) {
        bim = 0.5;
      } else {
        bim = -0.5;
      }

      if (y.im > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }

      z.re = (x.re * bim + x.im * d) / brm;
      z.im = (x.im * bim - x.re * d) / brm;
    } else {
      bim = y.re / y.im;
      d = y.im + bim * y.re;
      z.re = (bim * x.re + x.im) / d;
      z.im = (bim * x.im - x.re) / d;
    }
  }

  return z;
}

void c_mrdivide(const creal_T A[4], const creal_T B[4], creal_T y[4])
{
  creal_T b_A[4];
  creal_T b_B[4];
  int32_T r1;
  int32_T r2;
  creal_T a21;
  creal_T a22;
  creal_T Y[4];
  int32_T k;
  creal_T c_B;
  creal_T d_B;
  for (r1 = 0; r1 < 2; r1++) {
    for (r2 = 0; r2 < 2; r2++) {
      b_A[r2 + (r1 << 1)].re = B[r1 + (r2 << 1)].re;
      b_A[r2 + (r1 << 1)].im = -B[r1 + (r2 << 1)].im;
      b_B[r2 + (r1 << 1)].re = A[r1 + (r2 << 1)].re;
      b_B[r2 + (r1 << 1)].im = -A[r1 + (r2 << 1)].im;
    }
  }

  if (fabs(b_A[1].re) + fabs(b_A[1].im) > fabs(b_A[0].re) + fabs(b_A[0].im)) {
    r1 = 1;
    r2 = 0;
  } else {
    r1 = 0;
    r2 = 1;
  }

  a21 = c_eml_div(b_A[r2], b_A[r1]);
  a22.re = b_A[2 + r2].re - (a21.re * b_A[2 + r1].re - a21.im * b_A[2 + r1].im);
  a22.im = b_A[2 + r2].im - (a21.re * b_A[2 + r1].im + a21.im * b_A[2 + r1].re);
  for (k = 0; k < 2; k++) {
    c_B.re = b_B[r2 + (k << 1)].re - (b_B[r1 + (k << 1)].re * a21.re - b_B[r1 +
                                      (k << 1)].im * a21.im);
    c_B.im = b_B[r2 + (k << 1)].im - (b_B[r1 + (k << 1)].re * a21.im + b_B[r1 +
                                      (k << 1)].im * a21.re);
    Y[1 + (k << 1)] = c_eml_div(c_B, a22);
    d_B.re = b_B[r1 + (k << 1)].re - (Y[1 + (k << 1)].re * b_A[2 + r1].re - Y[1
      + (k << 1)].im * b_A[2 + r1].im);
    d_B.im = b_B[r1 + (k << 1)].im - (Y[1 + (k << 1)].re * b_A[2 + r1].im + Y[1
      + (k << 1)].im * b_A[2 + r1].re);
    Y[k << 1] = c_eml_div(d_B, b_A[r1]);
  }

  for (r1 = 0; r1 < 2; r1++) {
    for (r2 = 0; r2 < 2; r2++) {
      y[r2 + (r1 << 1)].re = Y[r1 + (r2 << 1)].re;
      y[r2 + (r1 << 1)].im = -Y[r1 + (r2 << 1)].im;
    }
  }
}

void mrdivide(const real_T A[12], const real_T B[12], real_T y[9])
{
  real_T tau[3];
  real_T b_A[12];
  real_T b_B[12];
  int8_T jpvt[3];
  real_T work[3];
  int32_T i4;
  int32_T itemp;
  real_T vn1[3];
  real_T vn2[3];
  int32_T k;
  int32_T iy;
  real_T b_y;
  real_T wj;
  real_T rankR;
  real_T t;
  int32_T i;
  int32_T i_i;
  int32_T ix;
  int32_T pvt;
  int32_T i_ip1;
  int32_T lastv;
  int32_T lastc;
  boolean_T exitg2;
  int32_T exitg1;
  real_T Y[9];
  for (i4 = 0; i4 < 3; i4++) {
    for (itemp = 0; itemp < 4; itemp++) {
      b_A[itemp + (i4 << 2)] = B[i4 + 3 * itemp];
      b_B[itemp + (i4 << 2)] = A[i4 + 3 * itemp];
    }

    jpvt[i4] = (int8_T)(1 + i4);
    work[i4] = 0.0;
  }

  k = 1;
  for (iy = 0; iy < 3; iy++) {
    b_y = 0.0;
    wj = 2.2250738585072014E-308;
    for (itemp = k; itemp <= k + 3; itemp++) {
      rankR = fabs(b_A[itemp - 1]);
      if (rankR > wj) {
        t = wj / rankR;
        b_y = 1.0 + b_y * t * t;
        wj = rankR;
      } else {
        t = rankR / wj;
        b_y += t * t;
      }
    }

    b_y = wj * sqrt(b_y);
    vn1[iy] = b_y;
    vn2[iy] = vn1[iy];
    k += 4;
  }

  for (i = 0; i < 3; i++) {
    i_i = i + (i << 2);
    itemp = 0;
    if (3 - i > 1) {
      ix = i;
      wj = fabs(vn1[i]);
      for (k = 1; k + 1 <= 3 - i; k++) {
        ix++;
        rankR = fabs(vn1[ix]);
        if (rankR > wj) {
          itemp = k;
          wj = rankR;
        }
      }
    }

    pvt = i + itemp;
    if (pvt + 1 != i + 1) {
      ix = pvt << 2;
      iy = i << 2;
      for (k = 0; k < 4; k++) {
        wj = b_A[ix];
        b_A[ix] = b_A[iy];
        b_A[iy] = wj;
        ix++;
        iy++;
      }

      itemp = jpvt[pvt];
      jpvt[pvt] = jpvt[i];
      jpvt[i] = (int8_T)itemp;
      vn1[pvt] = vn1[i];
      vn2[pvt] = vn2[i];
    }

    t = b_A[i_i];
    rankR = 0.0;
    wj = eml_xnrm2(3 - i, b_A, i_i + 2);
    if (wj != 0.0) {
      wj = rt_hypotd_snf(fabs(b_A[i_i]), fabs(wj));
      if (b_A[i_i] >= 0.0) {
        wj = -wj;
      }

      if (fabs(wj) < 1.0020841800044864E-292) {
        itemp = 0;
        do {
          itemp++;
          i4 = i_i - i;
          for (k = i_i + 1; k + 1 <= i4 + 4; k++) {
            b_A[k] *= 9.9792015476736E+291;
          }

          wj *= 9.9792015476736E+291;
          t *= 9.9792015476736E+291;
        } while (!(fabs(wj) >= 1.0020841800044864E-292));

        wj = rt_hypotd_snf(fabs(t), fabs(eml_xnrm2(3 - i, b_A, i_i + 2)));
        if (t >= 0.0) {
          wj = -wj;
        }

        rankR = (wj - t) / wj;
        t = 1.0 / (t - wj);
        i4 = i_i - i;
        for (k = i_i + 1; k + 1 <= i4 + 4; k++) {
          b_A[k] *= t;
        }

        for (k = 1; k <= itemp; k++) {
          wj *= 1.0020841800044864E-292;
        }

        t = wj;
      } else {
        rankR = (wj - b_A[i_i]) / wj;
        t = 1.0 / (b_A[i_i] - wj);
        i4 = i_i - i;
        for (k = i_i + 1; k + 1 <= i4 + 4; k++) {
          b_A[k] *= t;
        }

        t = wj;
      }
    }

    tau[i] = rankR;
    b_A[i_i] = t;
    if (i + 1 < 3) {
      t = b_A[i_i];
      b_A[i_i] = 1.0;
      i_ip1 = (i + ((i + 1) << 2)) + 1;
      if (tau[i] != 0.0) {
        lastv = 4 - i;
        itemp = i_i - i;
        while ((lastv > 0) && (b_A[itemp + 3] == 0.0)) {
          lastv--;
          itemp--;
        }

        lastc = 2 - i;
        exitg2 = FALSE;
        while ((exitg2 == FALSE) && (lastc > 0)) {
          itemp = i_ip1 + ((lastc - 1) << 2);
          k = itemp;
          do {
            exitg1 = 0;
            if (k <= (itemp + lastv) - 1) {
              if (b_A[k - 1] != 0.0) {
                exitg1 = 1;
              } else {
                k++;
              }
            } else {
              lastc--;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = TRUE;
          }
        }
      } else {
        lastv = 0;
        lastc = 0;
      }

      if (lastv > 0) {
        if (lastc == 0) {
        } else {
          for (iy = 1; iy <= lastc; iy++) {
            work[iy - 1] = 0.0;
          }

          iy = 0;
          i4 = i_ip1 + ((lastc - 1) << 2);
          for (pvt = i_ip1; pvt <= i4; pvt += 4) {
            ix = i_i;
            wj = 0.0;
            itemp = (pvt + lastv) - 1;
            for (k = pvt; k <= itemp; k++) {
              wj += b_A[k - 1] * b_A[ix];
              ix++;
            }

            work[iy] += wj;
            iy++;
          }
        }

        if (-tau[i] == 0.0) {
        } else {
          itemp = i_ip1 - 1;
          pvt = 0;
          for (iy = 1; iy <= lastc; iy++) {
            if (work[pvt] != 0.0) {
              wj = work[pvt] * -tau[i];
              ix = i_i;
              i4 = lastv + itemp;
              for (k = itemp; k + 1 <= i4; k++) {
                b_A[k] += b_A[ix] * wj;
                ix++;
              }
            }

            pvt++;
            itemp += 4;
          }
        }
      }

      b_A[i_i] = t;
    }

    for (iy = i + 1; iy + 1 < 4; iy++) {
      if (vn1[iy] != 0.0) {
        rankR = fabs(b_A[i + (iy << 2)]) / vn1[iy];
        b_y = rankR * rankR;
        rankR = 1.0 - rankR * rankR;
        if (1.0 - b_y < 0.0) {
          rankR = 0.0;
        }

        wj = vn1[iy] / vn2[iy];
        if (rankR * (wj * wj) <= 1.4901161193847656E-8) {
          vn1[iy] = b_eml_xnrm2(3 - i, b_A, (i + (iy << 2)) + 2);
          vn2[iy] = vn1[iy];
        } else {
          vn1[iy] *= sqrt(rankR);
        }
      }
    }
  }

  rankR = 0.0;
  k = 0;
  while ((k < 3) && (!(fabs(b_A[k + (k << 2)]) <= 4.0 * fabs(b_A[0]) *
                       2.2204460492503131E-16))) {
    rankR++;
    k++;
  }

  memset(&Y[0], 0, 9U * sizeof(real_T));
  for (iy = 0; iy < 3; iy++) {
    if (tau[iy] != 0.0) {
      for (k = 0; k < 3; k++) {
        wj = b_B[iy + (k << 2)];
        for (i = 0; i <= 2 - iy; i++) {
          itemp = (iy + i) + 1;
          wj += b_A[itemp + (iy << 2)] * b_B[itemp + (k << 2)];
        }

        wj *= tau[iy];
        if (wj != 0.0) {
          b_B[iy + (k << 2)] -= wj;
          for (i = 0; i <= 2 - iy; i++) {
            itemp = (iy + i) + 1;
            b_B[itemp + (k << 2)] -= b_A[itemp + (iy << 2)] * wj;
          }
        }
      }
    }
  }

  for (k = 0; k < 3; k++) {
    for (i = 0; i < (int32_T)rankR; i++) {
      Y[(jpvt[i] + 3 * k) - 1] = b_B[i + (k << 2)];
    }

    for (iy = 0; iy < (int32_T)-(1.0 + (-1.0 - rankR)); iy++) {
      wj = rankR + -(real_T)iy;
      Y[(jpvt[(int32_T)wj - 1] + 3 * k) - 1] /= b_A[((int32_T)wj + (((int32_T)wj
        - 1) << 2)) - 1];
      for (i = 0; i <= (int32_T)wj - 2; i++) {
        Y[(jpvt[i] + 3 * k) - 1] -= Y[(jpvt[(int32_T)wj - 1] + 3 * k) - 1] *
          b_A[i + (((int32_T)wj - 1) << 2)];
      }
    }
  }

  for (i4 = 0; i4 < 3; i4++) {
    for (itemp = 0; itemp < 3; itemp++) {
      y[itemp + 3 * i4] = Y[i4 + 3 * itemp];
    }
  }
}

/* End of code generation (mrdivide.cpp) */
