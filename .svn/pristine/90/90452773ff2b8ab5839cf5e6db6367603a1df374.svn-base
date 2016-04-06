/*
 * svd.cpp
 *
 * Code generation for function 'svd'
 *
 * C source code generated on: Tue Apr  1 12:30:12 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Optimal_affine_tracking_3d16_fast_realtime.h"
#include "svd.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
static void b_eml_xaxpy(int32_T n, real_T a, const real_T x[16], int32_T ix0,
  real_T y[4], int32_T iy0);
static void b_eml_xscal(int32_T n, real_T a, real_T x[4], int32_T ix0);
static void c_eml_xaxpy(int32_T n, real_T a, const real_T x[4], int32_T ix0,
  real_T y[16], int32_T iy0);
static real_T c_eml_xnrm2(int32_T n, const real_T x[16], int32_T ix0);
static real_T d_eml_xnrm2(int32_T n, const real_T x[4], int32_T ix0);
static real_T eml_div(real_T x, real_T y);
static void eml_xaxpy(int32_T n, real_T a, int32_T ix0, real_T y[16], int32_T
                      iy0);
static real_T eml_xdotc(int32_T n, const real_T x[16], int32_T ix0, const real_T
  y[16], int32_T iy0);
static void eml_xgesvd(const real_T A[16], real_T U[16], real_T S[4], real_T V
  [16]);
static void eml_xrot(real_T x[16], int32_T ix0, int32_T iy0, real_T c, real_T s);
static void eml_xrotg(real_T *a, real_T *b, real_T *c, real_T *s);
static void eml_xscal(int32_T n, real_T a, real_T x[16], int32_T ix0);
static void eml_xswap(real_T x[16], int32_T ix0, int32_T iy0);

/* Function Definitions */
static void b_eml_xaxpy(int32_T n, real_T a, const real_T x[16], int32_T ix0,
  real_T y[4], int32_T iy0)
{
  int32_T ix;
  int32_T iy;
  int32_T k;
  if (a == 0.0) {
  } else {
    ix = ix0 - 1;
    iy = iy0 - 1;
    for (k = 0; k < n; k++) {
      y[iy] += a * x[ix];
      ix++;
      iy++;
    }
  }
}

static void b_eml_xscal(int32_T n, real_T a, real_T x[4], int32_T ix0)
{
  int32_T i23;
  int32_T k;
  i23 = (ix0 + n) - 1;
  for (k = ix0; k <= i23; k++) {
    x[k - 1] *= a;
  }
}

static void c_eml_xaxpy(int32_T n, real_T a, const real_T x[4], int32_T ix0,
  real_T y[16], int32_T iy0)
{
  int32_T ix;
  int32_T iy;
  int32_T k;
  if (a == 0.0) {
  } else {
    ix = ix0 - 1;
    iy = iy0 - 1;
    for (k = 0; k < n; k++) {
      y[iy] += a * x[ix];
      ix++;
      iy++;
    }
  }
}

static real_T c_eml_xnrm2(int32_T n, const real_T x[16], int32_T ix0)
{
  real_T y;
  real_T scale;
  int32_T kend;
  int32_T k;
  real_T absxk;
  real_T t;
  y = 0.0;
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

  return scale * sqrt(y);
}

static real_T d_eml_xnrm2(int32_T n, const real_T x[4], int32_T ix0)
{
  real_T y;
  real_T scale;
  int32_T kend;
  int32_T k;
  real_T absxk;
  real_T t;
  y = 0.0;
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

  return scale * sqrt(y);
}

static real_T eml_div(real_T x, real_T y)
{
  return x / y;
}

static void eml_xaxpy(int32_T n, real_T a, int32_T ix0, real_T y[16], int32_T
                      iy0)
{
  int32_T ix;
  int32_T iy;
  int32_T k;
  if ((n < 1) || (a == 0.0)) {
  } else {
    ix = ix0 - 1;
    iy = iy0 - 1;
    for (k = 0; k < n; k++) {
      y[iy] += a * y[ix];
      ix++;
      iy++;
    }
  }
}

static real_T eml_xdotc(int32_T n, const real_T x[16], int32_T ix0, const real_T
  y[16], int32_T iy0)
{
  real_T d;
  int32_T ix;
  int32_T iy;
  int32_T k;
  d = 0.0;
  if (n < 1) {
  } else {
    ix = ix0;
    iy = iy0;
    for (k = 1; k <= n; k++) {
      d += x[ix - 1] * y[iy - 1];
      ix++;
      iy++;
    }
  }

  return d;
}

static void eml_xgesvd(const real_T A[16], real_T U[16], real_T S[4], real_T V
  [16])
{
  real_T b_A[16];
  real_T s[4];
  real_T e[4];
  real_T work[4];
  int32_T i;
  int32_T q;
  int32_T qs;
  real_T ztest0;
  int32_T ii;
  int32_T m;
  real_T rt;
  real_T ztest;
  int32_T iter;
  real_T tiny;
  real_T snorm;
  int32_T exitg3;
  boolean_T exitg2;
  real_T sn;
  real_T varargin_1[5];
  boolean_T exitg1;
  real_T sqds;
  real_T b;
  memcpy(&b_A[0], (real_T *)&A[0], sizeof(real_T) << 4);
  for (i = 0; i < 4; i++) {
    s[i] = 0.0;
    e[i] = 0.0;
    work[i] = 0.0;
  }

  for (i = 0; i < 16; i++) {
    U[i] = 0.0;
    V[i] = 0.0;
  }

  for (q = 0; q < 3; q++) {
    qs = q + (q << 2);
    ztest0 = c_eml_xnrm2(4 - q, b_A, qs + 1);
    if (ztest0 > 0.0) {
      if (b_A[qs] < 0.0) {
        s[q] = -ztest0;
      } else {
        s[q] = ztest0;
      }

      eml_xscal(4 - q, eml_div(1.0, s[q]), b_A, qs + 1);
      b_A[qs]++;
      s[q] = -s[q];
    } else {
      s[q] = 0.0;
    }

    for (ii = q + 1; ii + 1 < 5; ii++) {
      i = q + (ii << 2);
      if (s[q] != 0.0) {
        ztest0 = -eml_div(eml_xdotc(4 - q, b_A, qs + 1, b_A, i + 1), b_A[q + (q <<
          2)]);
        eml_xaxpy(4 - q, ztest0, qs + 1, b_A, i + 1);
      }

      e[ii] = b_A[i];
    }

    for (ii = q; ii + 1 < 5; ii++) {
      U[ii + (q << 2)] = b_A[ii + (q << 2)];
    }

    if (q + 1 <= 2) {
      ztest0 = d_eml_xnrm2(3 - q, e, q + 2);
      if (ztest0 == 0.0) {
        e[q] = 0.0;
      } else {
        if (e[q + 1] < 0.0) {
          e[q] = -ztest0;
        } else {
          e[q] = ztest0;
        }

        ztest0 = eml_div(1.0, e[q]);
        b_eml_xscal(3 - q, ztest0, e, q + 2);
        e[q + 1]++;
      }

      e[q] = -e[q];
      if (e[q] != 0.0) {
        for (ii = q + 1; ii + 1 < 5; ii++) {
          work[ii] = 0.0;
        }

        for (ii = q + 1; ii + 1 < 5; ii++) {
          b_eml_xaxpy(3 - q, e[ii], b_A, (q + (ii << 2)) + 2, work, q + 2);
        }

        for (ii = q + 1; ii + 1 < 5; ii++) {
          c_eml_xaxpy(3 - q, eml_div(-e[ii], e[q + 1]), work, q + 2, b_A, (q +
            (ii << 2)) + 2);
        }
      }

      for (ii = q + 1; ii + 1 < 5; ii++) {
        V[ii + (q << 2)] = e[ii];
      }
    }
  }

  m = 2;
  s[3] = b_A[15];
  e[2] = b_A[14];
  e[3] = 0.0;
  for (ii = 0; ii < 4; ii++) {
    U[12 + ii] = 0.0;
  }

  U[15] = 1.0;
  for (q = 2; q > -1; q += -1) {
    qs = q + (q << 2);
    if (s[q] != 0.0) {
      for (ii = q + 1; ii + 1 < 5; ii++) {
        i = (q + (ii << 2)) + 1;
        ztest0 = -eml_div(eml_xdotc(4 - q, U, qs + 1, U, i), U[qs]);
        eml_xaxpy(4 - q, ztest0, qs + 1, U, i);
      }

      for (ii = q; ii + 1 < 5; ii++) {
        U[ii + (q << 2)] = -U[ii + (q << 2)];
      }

      U[qs]++;
      for (ii = 1; ii <= q; ii++) {
        U[(ii + (q << 2)) - 1] = 0.0;
      }
    } else {
      for (ii = 0; ii < 4; ii++) {
        U[ii + (q << 2)] = 0.0;
      }

      U[qs] = 1.0;
    }
  }

  for (q = 3; q > -1; q += -1) {
    if ((q + 1 <= 2) && (e[q] != 0.0)) {
      i = (q + (q << 2)) + 2;
      for (ii = q + 1; ii + 1 < 5; ii++) {
        qs = (q + (ii << 2)) + 2;
        ztest0 = -eml_div(eml_xdotc(3 - q, V, i, V, qs), V[i - 1]);
        eml_xaxpy(3 - q, ztest0, i, V, qs);
      }
    }

    for (ii = 0; ii < 4; ii++) {
      V[ii + (q << 2)] = 0.0;
    }

    V[q + (q << 2)] = 1.0;
  }

  for (q = 0; q < 4; q++) {
    ztest0 = e[q];
    if (s[q] != 0.0) {
      rt = fabs(s[q]);
      ztest = eml_div(s[q], rt);
      s[q] = rt;
      if (q + 1 < 4) {
        ztest0 = eml_div(e[q], ztest);
      }

      eml_xscal(4, ztest, U, (q << 2) + 1);
    }

    if ((q + 1 < 4) && (ztest0 != 0.0)) {
      rt = fabs(ztest0);
      ztest = eml_div(rt, ztest0);
      ztest0 = rt;
      s[q + 1] *= ztest;
      eml_xscal(4, ztest, V, ((q + 1) << 2) + 1);
    }

    e[q] = ztest0;
  }

  iter = 0;
  tiny = eml_div(2.2250738585072014E-308, 2.2204460492503131E-16);
  snorm = 0.0;
  for (ii = 0; ii < 4; ii++) {
    ztest0 = fabs(s[ii]);
    ztest = fabs(e[ii]);
    if ((ztest0 >= ztest) || rtIsNaN(ztest)) {
    } else {
      ztest0 = ztest;
    }

    if ((snorm >= ztest0) || rtIsNaN(ztest0)) {
    } else {
      snorm = ztest0;
    }
  }

  while ((m + 2 > 0) && (!(iter >= 75))) {
    ii = m;
    do {
      exitg3 = 0;
      q = ii + 1;
      if (ii + 1 == 0) {
        exitg3 = 1;
      } else {
        ztest0 = fabs(e[ii]);
        if ((ztest0 <= 2.2204460492503131E-16 * (fabs(s[ii]) + fabs(s[ii + 1])))
            || (ztest0 <= tiny) || ((iter > 20) && (ztest0 <=
              2.2204460492503131E-16 * snorm))) {
          e[ii] = 0.0;
          exitg3 = 1;
        } else {
          ii--;
        }
      }
    } while (exitg3 == 0);

    if (ii + 1 == m + 1) {
      i = 4;
    } else {
      qs = m + 2;
      i = m + 2;
      exitg2 = FALSE;
      while ((exitg2 == FALSE) && (i >= ii + 1)) {
        qs = i;
        if (i == ii + 1) {
          exitg2 = TRUE;
        } else {
          ztest0 = 0.0;
          if (i < m + 2) {
            ztest0 = fabs(e[i - 1]);
          }

          if (i > ii + 2) {
            ztest0 += fabs(e[i - 2]);
          }

          ztest = fabs(s[i - 1]);
          if ((ztest <= 2.2204460492503131E-16 * ztest0) || (ztest <= tiny)) {
            s[i - 1] = 0.0;
            exitg2 = TRUE;
          } else {
            i--;
          }
        }
      }

      if (qs == ii + 1) {
        i = 3;
      } else if (qs == m + 2) {
        i = 1;
      } else {
        i = 2;
        q = qs;
      }
    }

    switch (i) {
     case 1:
      ztest = e[m];
      e[m] = 0.0;
      for (i = m; i + 1 >= q + 1; i--) {
        ztest0 = s[i];
        eml_xrotg(&ztest0, &ztest, &rt, &sn);
        s[i] = ztest0;
        if (i + 1 > q + 1) {
          ztest = -sn * e[i - 1];
          e[i - 1] *= rt;
        }

        eml_xrot(V, (i << 2) + 1, ((m + 1) << 2) + 1, rt, sn);
      }
      break;

     case 2:
      ztest = e[q - 1];
      e[q - 1] = 0.0;
      for (i = q; i + 1 <= m + 2; i++) {
        eml_xrotg(&s[i], &ztest, &rt, &sn);
        ztest = -sn * e[i];
        e[i] *= rt;
        eml_xrot(U, (i << 2) + 1, ((q - 1) << 2) + 1, rt, sn);
      }
      break;

     case 3:
      varargin_1[0] = fabs(s[m + 1]);
      varargin_1[1] = fabs(s[m]);
      varargin_1[2] = fabs(e[m]);
      varargin_1[3] = fabs(s[q]);
      varargin_1[4] = fabs(e[q]);
      i = 1;
      sn = varargin_1[0];
      if (rtIsNaN(varargin_1[0])) {
        qs = 2;
        exitg1 = FALSE;
        while ((exitg1 == FALSE) && (qs < 6)) {
          i = qs;
          if (!rtIsNaN(varargin_1[qs - 1])) {
            sn = varargin_1[qs - 1];
            exitg1 = TRUE;
          } else {
            qs++;
          }
        }
      }

      if (i < 5) {
        while (i + 1 < 6) {
          if (varargin_1[i] > sn) {
            sn = varargin_1[i];
          }

          i++;
        }
      }

      rt = eml_div(s[m + 1], sn);
      ztest0 = eml_div(s[m], sn);
      ztest = eml_div(e[m], sn);
      sqds = eml_div(s[q], sn);
      b = eml_div((ztest0 + rt) * (ztest0 - rt) + ztest * ztest, 2.0);
      ztest0 = rt * ztest;
      ztest0 *= ztest0;
      ztest = 0.0;
      if ((b != 0.0) || (ztest0 != 0.0)) {
        ztest = sqrt(b * b + ztest0);
        if (b < 0.0) {
          ztest = -ztest;
        }

        ztest = eml_div(ztest0, b + ztest);
      }

      ztest += (sqds + rt) * (sqds - rt);
      ztest0 = sqds * eml_div(e[q], sn);
      for (i = q + 1; i <= m + 1; i++) {
        eml_xrotg(&ztest, &ztest0, &rt, &sn);
        if (i > q + 1) {
          e[i - 2] = ztest;
        }

        ztest0 = rt * s[i - 1];
        ztest = sn * e[i - 1];
        e[i - 1] = rt * e[i - 1] - sn * s[i - 1];
        b = s[i];
        s[i] *= rt;
        eml_xrot(V, ((i - 1) << 2) + 1, (i << 2) + 1, rt, sn);
        s[i - 1] = ztest0 + ztest;
        ztest0 = sn * b;
        eml_xrotg(&s[i - 1], &ztest0, &rt, &sn);
        ztest = rt * e[i - 1] + sn * s[i];
        s[i] = -sn * e[i - 1] + rt * s[i];
        ztest0 = sn * e[i];
        e[i] *= rt;
        eml_xrot(U, ((i - 1) << 2) + 1, (i << 2) + 1, rt, sn);
      }

      e[m] = ztest;
      iter++;
      break;

     default:
      if (s[q] < 0.0) {
        s[q] = -s[q];
        eml_xscal(4, -1.0, V, (q << 2) + 1);
      }

      i = q + 1;
      while ((q + 1 < 4) && (s[q] < s[i])) {
        rt = s[q];
        s[q] = s[i];
        s[i] = rt;
        eml_xswap(V, (q << 2) + 1, ((q + 1) << 2) + 1);
        eml_xswap(U, (q << 2) + 1, ((q + 1) << 2) + 1);
        q = i;
        i++;
      }

      iter = 0;
      m--;
      break;
    }
  }

  for (i = 0; i < 4; i++) {
    S[i] = s[i];
  }
}

static void eml_xrot(real_T x[16], int32_T ix0, int32_T iy0, real_T c, real_T s)
{
  int32_T ix;
  int32_T iy;
  int32_T k;
  real_T y;
  real_T b_y;
  ix = ix0 - 1;
  iy = iy0 - 1;
  for (k = 0; k < 4; k++) {
    y = c * x[ix];
    b_y = s * x[iy];
    x[iy] = c * x[iy] - s * x[ix];
    x[ix] = y + b_y;
    iy++;
    ix++;
  }
}

static void eml_xrotg(real_T *a, real_T *b, real_T *c, real_T *s)
{
  real_T roe;
  real_T absa;
  real_T absb;
  real_T scale;
  real_T ads;
  real_T bds;
  roe = *b;
  absa = fabs(*a);
  absb = fabs(*b);
  if (absa > absb) {
    roe = *a;
  }

  scale = absa + absb;
  if (scale == 0.0) {
    *s = 0.0;
    *c = 1.0;
    ads = 0.0;
    scale = 0.0;
  } else {
    ads = absa / scale;
    bds = absb / scale;
    ads = scale * sqrt(ads * ads + bds * bds);
    if (roe < 0.0) {
      ads = -ads;
    }

    *c = *a / ads;
    *s = *b / ads;
    if (absa > absb) {
      scale = *s;
    } else if (*c != 0.0) {
      scale = 1.0 / *c;
    } else {
      scale = 1.0;
    }
  }

  *a = ads;
  *b = scale;
}

static void eml_xscal(int32_T n, real_T a, real_T x[16], int32_T ix0)
{
  int32_T i22;
  int32_T k;
  i22 = (ix0 + n) - 1;
  for (k = ix0; k <= i22; k++) {
    x[k - 1] *= a;
  }
}

static void eml_xswap(real_T x[16], int32_T ix0, int32_T iy0)
{
  int32_T ix;
  int32_T iy;
  int32_T k;
  real_T temp;
  ix = ix0 - 1;
  iy = iy0 - 1;
  for (k = 0; k < 4; k++) {
    temp = x[ix];
    x[ix] = x[iy];
    x[iy] = temp;
    ix++;
    iy++;
  }
}

void svd(const real_T A[16], real_T U[16], real_T S[16], real_T V[16])
{
  real_T s[4];
  int32_T k;
  eml_xgesvd(A, U, s, V);
  memset(&S[0], 0, sizeof(real_T) << 4);
  for (k = 0; k < 4; k++) {
    S[k + (k << 2)] = s[k];
  }
}

/* End of code generation (svd.cpp) */
