/*
 * eig.cpp
 *
 * Code generation for function 'eig'
 *
 * C source code generated on: Tue Apr  1 12:30:12 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Optimal_affine_tracking_3d16_fast_realtime.h"
#include "eig.h"
#include "mrdivide.h"
#include "log.h"
#include "sqrt.h"
#include "Optimal_affine_tracking_3d16_fast_realtime_rtwutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
static creal_T b_eml_div(const creal_T x, real_T y);
static void b_eml_matlab_zggev(creal_T A[4], real_T *info, creal_T alpha1[2],
  creal_T beta1[2], creal_T V[4]);
static void b_eml_matlab_zhgeqz(creal_T A[4], int32_T ilo, int32_T ihi, creal_T
  Z[4], real_T *info, creal_T alpha1[2], creal_T beta1[2]);
static void b_eml_matlab_zlartg(const creal_T f, const creal_T g, real_T *cs,
  creal_T *sn);
static void b_eml_matlab_ztgevc(const creal_T A[4], creal_T V[4]);
static void eml_matlab_zggev(creal_T A[9], real_T *info, creal_T alpha1[3],
  creal_T beta1[3], creal_T V[9]);
static void eml_matlab_zhgeqz(creal_T A[9], int32_T ilo, int32_T ihi, creal_T Z
  [9], real_T *info, creal_T alpha1[3], creal_T beta1[3]);
static void eml_matlab_zlartg(const creal_T f, const creal_T g, real_T *cs,
  creal_T *sn, creal_T *r);
static void eml_matlab_ztgevc(const creal_T A[9], creal_T V[9]);

/* Function Definitions */
static creal_T b_eml_div(const creal_T x, real_T y)
{
  creal_T z;
  if (x.im == 0.0) {
    z.re = x.re / y;
    z.im = 0.0;
  } else if (x.re == 0.0) {
    z.re = 0.0;
    z.im = x.im / y;
  } else {
    z.re = x.re / y;
    z.im = x.im / y;
  }

  return z;
}

static void b_eml_matlab_zggev(creal_T A[4], real_T *info, creal_T alpha1[2],
  creal_T beta1[2], creal_T V[4])
{
  real_T anrm;
  int32_T nzcount;
  boolean_T exitg5;
  real_T absxk;
  int32_T i;
  boolean_T ilascl;
  real_T anrmto;
  real_T ctoc;
  boolean_T notdone;
  real_T cfrom1;
  real_T cto1;
  real_T mul;
  creal_T b_A[4];
  int8_T rscale[2];
  int32_T ilo;
  int32_T ihi;
  int32_T j;
  int32_T ii;
  boolean_T exitg3;
  int32_T jj;
  boolean_T exitg4;
  boolean_T guard2 = FALSE;
  boolean_T exitg1;
  boolean_T exitg2;
  boolean_T guard1 = FALSE;
  static const int8_T iv7[4] = { 1, 0, 0, 1 };

  *info = 0.0;
  anrm = 0.0;
  nzcount = 0;
  exitg5 = FALSE;
  while ((exitg5 == FALSE) && (nzcount < 4)) {
    absxk = rt_hypotd_snf(fabs(A[nzcount].re), fabs(A[nzcount].im));
    if (rtIsNaN(absxk)) {
      anrm = rtNaN;
      exitg5 = TRUE;
    } else {
      if (absxk > anrm) {
        anrm = absxk;
      }

      nzcount++;
    }
  }

  if (!((!rtIsInf(anrm)) && (!rtIsNaN(anrm)))) {
    for (i = 0; i < 2; i++) {
      alpha1[i].re = rtNaN;
      alpha1[i].im = 0.0;
      beta1[i].re = rtNaN;
      beta1[i].im = 0.0;
    }

    for (nzcount = 0; nzcount < 4; nzcount++) {
      V[nzcount].re = rtNaN;
      V[nzcount].im = 0.0;
    }
  } else {
    ilascl = FALSE;
    anrmto = anrm;
    if ((anrm > 0.0) && (anrm < 6.7178761075670888E-139)) {
      anrmto = 6.7178761075670888E-139;
      ilascl = TRUE;
    } else {
      if (anrm > 1.4885657073574029E+138) {
        anrmto = 1.4885657073574029E+138;
        ilascl = TRUE;
      }
    }

    if (ilascl) {
      absxk = anrm;
      ctoc = anrmto;
      notdone = TRUE;
      while (notdone) {
        cfrom1 = absxk * 2.0041683600089728E-292;
        cto1 = ctoc / 4.9896007738368E+291;
        if ((cfrom1 > ctoc) && (ctoc != 0.0)) {
          mul = 2.0041683600089728E-292;
          absxk = cfrom1;
        } else if (cto1 > absxk) {
          mul = 4.9896007738368E+291;
          ctoc = cto1;
        } else {
          mul = ctoc / absxk;
          notdone = FALSE;
        }

        for (nzcount = 0; nzcount < 4; nzcount++) {
          A[nzcount].re *= mul;
          A[nzcount].im *= mul;
        }
      }
    }

    memcpy(&b_A[0], &A[0], sizeof(creal_T) << 2);
    for (i = 0; i < 2; i++) {
      rscale[i] = 0;
    }

    ilo = 1;
    ihi = 2;
    i = 0;
    j = 0;
    notdone = FALSE;
    ii = 2;
    exitg3 = FALSE;
    while ((exitg3 == FALSE) && (ii > 0)) {
      nzcount = 0;
      i = ii;
      j = 2;
      jj = 1;
      exitg4 = FALSE;
      while ((exitg4 == FALSE) && (jj <= 2)) {
        guard2 = FALSE;
        if ((A[(ii + ((jj - 1) << 1)) - 1].re != 0.0) || (A[(ii + ((jj - 1) << 1))
             - 1].im != 0.0) || (ii == jj)) {
          if (nzcount == 0) {
            j = jj;
            nzcount = 1;
            guard2 = TRUE;
          } else {
            nzcount = 2;
            exitg4 = TRUE;
          }
        } else {
          guard2 = TRUE;
        }

        if (guard2 == TRUE) {
          jj++;
        }
      }

      if (nzcount < 2) {
        notdone = TRUE;
        exitg3 = TRUE;
      } else {
        ii--;
      }
    }

    if (!notdone) {
      i = 0;
      j = 0;
      notdone = FALSE;
      jj = 1;
      exitg1 = FALSE;
      while ((exitg1 == FALSE) && (jj <= 2)) {
        nzcount = 0;
        i = 2;
        j = jj;
        ii = 1;
        exitg2 = FALSE;
        while ((exitg2 == FALSE) && (ii <= 2)) {
          guard1 = FALSE;
          if ((b_A[(ii + ((jj - 1) << 1)) - 1].re != 0.0) || (b_A[(ii + ((jj - 1)
                 << 1)) - 1].im != 0.0) || (ii == jj)) {
            if (nzcount == 0) {
              i = ii;
              nzcount = 1;
              guard1 = TRUE;
            } else {
              nzcount = 2;
              exitg2 = TRUE;
            }
          } else {
            guard1 = TRUE;
          }

          if (guard1 == TRUE) {
            ii++;
          }
        }

        if (nzcount < 2) {
          notdone = TRUE;
          exitg1 = TRUE;
        } else {
          jj++;
        }
      }

      if (!notdone) {
      } else {
        if (i != 1) {
          for (nzcount = 0; nzcount + 1 < 3; nzcount++) {
            absxk = b_A[(i + (nzcount << 1)) - 1].re;
            ctoc = b_A[(i + (nzcount << 1)) - 1].im;
            b_A[(i + (nzcount << 1)) - 1] = b_A[nzcount << 1];
            b_A[nzcount << 1].re = absxk;
            b_A[nzcount << 1].im = ctoc;
          }
        }

        if (j != 1) {
          for (nzcount = 0; nzcount < 2; nzcount++) {
            absxk = b_A[nzcount + ((j - 1) << 1)].re;
            ctoc = b_A[nzcount + ((j - 1) << 1)].im;
            b_A[nzcount + ((j - 1) << 1)] = b_A[nzcount];
            b_A[nzcount].re = absxk;
            b_A[nzcount].im = ctoc;
          }
        }

        rscale[0] = (int8_T)j;
        ilo = 2;
        rscale[1] = 2;
      }
    } else {
      if (i != 2) {
        for (nzcount = 0; nzcount < 2; nzcount++) {
          absxk = b_A[nzcount << 1].re;
          ctoc = b_A[nzcount << 1].im;
          b_A[nzcount << 1] = b_A[1 + (nzcount << 1)];
          b_A[1 + (nzcount << 1)].re = absxk;
          b_A[1 + (nzcount << 1)].im = ctoc;
        }
      }

      if (j != 2) {
        for (nzcount = 0; nzcount < 2; nzcount++) {
          absxk = b_A[nzcount].re;
          ctoc = b_A[nzcount].im;
          b_A[nzcount] = b_A[2 + nzcount];
          b_A[2 + nzcount].re = absxk;
          b_A[2 + nzcount].im = ctoc;
        }
      }

      rscale[1] = (int8_T)j;
      ihi = 1;
      rscale[0] = 1;
    }

    for (nzcount = 0; nzcount < 4; nzcount++) {
      V[nzcount].re = (real_T)iv7[nzcount];
      V[nzcount].im = 0.0;
    }

    b_eml_matlab_zhgeqz(b_A, ilo, ihi, V, info, alpha1, beta1);
    if (*info != 0.0) {
    } else {
      b_eml_matlab_ztgevc(b_A, V);
      if ((ilo > 1) && (rscale[0] != 1)) {
        for (j = 0; j < 2; j++) {
          absxk = V[j << 1].re;
          ctoc = V[j << 1].im;
          V[j << 1] = V[(rscale[0] + (j << 1)) - 1];
          V[(rscale[0] + (j << 1)) - 1].re = absxk;
          V[(rscale[0] + (j << 1)) - 1].im = ctoc;
        }
      }

      if ((ihi < 2) && (rscale[1] != 2)) {
        for (j = 0; j < 2; j++) {
          absxk = V[1 + (j << 1)].re;
          ctoc = V[1 + (j << 1)].im;
          V[1 + (j << 1)] = V[j << 1];
          V[(rscale[1] + (j << 1)) - 1].re = absxk;
          V[(rscale[1] + (j << 1)) - 1].im = ctoc;
        }
      }

      for (nzcount = 0; nzcount < 2; nzcount++) {
        absxk = fabs(V[1 + (nzcount << 1)].re) + fabs(V[1 + (nzcount << 1)].im);
        ctoc = fabs(V[nzcount << 1].re) + fabs(V[nzcount << 1].im);
        if (absxk > ctoc) {
          ctoc = absxk;
        }

        if (ctoc >= 6.7178761075670888E-139) {
          ctoc = 1.0 / ctoc;
          for (ii = 0; ii < 2; ii++) {
            V[ii + (nzcount << 1)].re *= ctoc;
            V[ii + (nzcount << 1)].im *= ctoc;
          }
        }
      }

      if (ilascl) {
        notdone = TRUE;
        while (notdone) {
          cfrom1 = anrmto * 2.0041683600089728E-292;
          cto1 = anrm / 4.9896007738368E+291;
          if ((cfrom1 > anrm) && (anrm != 0.0)) {
            mul = 2.0041683600089728E-292;
            anrmto = cfrom1;
          } else if (cto1 > anrmto) {
            mul = 4.9896007738368E+291;
            anrm = cto1;
          } else {
            mul = anrm / anrmto;
            notdone = FALSE;
          }

          for (nzcount = 0; nzcount < 2; nzcount++) {
            alpha1[nzcount].re *= mul;
            alpha1[nzcount].im *= mul;
          }
        }
      }
    }
  }
}

static void b_eml_matlab_zhgeqz(creal_T A[4], int32_T ilo, int32_T ihi, creal_T
  Z[4], real_T *info, creal_T alpha1[2], creal_T beta1[2])
{
  int32_T i;
  real_T eshift_re;
  real_T eshift_im;
  creal_T ctemp;
  real_T rho_re;
  real_T rho_im;
  real_T anorm;
  real_T scale;
  real_T sumsq;
  boolean_T firstNonZero;
  int32_T j;
  int32_T istart;
  real_T reAij;
  real_T imAij;
  real_T temp2;
  real_T b_atol;
  real_T ascale;
  boolean_T failed;
  boolean_T guard1 = FALSE;
  boolean_T guard2 = FALSE;
  int32_T ifirst;
  int32_T ilast;
  int32_T ilastm1;
  int32_T iiter;
  boolean_T goto60;
  boolean_T goto70;
  boolean_T goto90;
  int32_T jiter;
  int32_T exitg1;
  boolean_T exitg3;
  boolean_T b_guard1 = FALSE;
  creal_T b_A;
  creal_T t1;
  creal_T d;
  real_T sigma1_im;
  boolean_T exitg2;
  creal_T c_A;
  *info = 0.0;
  for (i = 0; i < 2; i++) {
    alpha1[i].re = 0.0;
    alpha1[i].im = 0.0;
    beta1[i].re = 1.0;
    beta1[i].im = 0.0;
  }

  eshift_re = 0.0;
  eshift_im = 0.0;
  ctemp.re = 0.0;
  ctemp.im = 0.0;
  rho_re = 0.0;
  rho_im = 0.0;
  anorm = 0.0;
  if (ilo > ihi) {
  } else {
    scale = 0.0;
    sumsq = 0.0;
    firstNonZero = TRUE;
    for (j = ilo; j <= ihi; j++) {
      istart = j + 1;
      if (ihi < j + 1) {
        istart = ihi;
      }

      for (i = ilo; i <= istart; i++) {
        reAij = A[(i + ((j - 1) << 1)) - 1].re;
        imAij = A[(i + ((j - 1) << 1)) - 1].im;
        if (reAij != 0.0) {
          anorm = fabs(reAij);
          if (firstNonZero) {
            sumsq = 1.0;
            scale = anorm;
            firstNonZero = FALSE;
          } else if (scale < anorm) {
            temp2 = scale / anorm;
            sumsq = 1.0 + sumsq * temp2 * temp2;
            scale = anorm;
          } else {
            temp2 = anorm / scale;
            sumsq += temp2 * temp2;
          }
        }

        if (imAij != 0.0) {
          anorm = fabs(imAij);
          if (firstNonZero) {
            sumsq = 1.0;
            scale = anorm;
            firstNonZero = FALSE;
          } else if (scale < anorm) {
            temp2 = scale / anorm;
            sumsq = 1.0 + sumsq * temp2 * temp2;
            scale = anorm;
          } else {
            temp2 = anorm / scale;
            sumsq += temp2 * temp2;
          }
        }
      }
    }

    anorm = scale * sqrt(sumsq);
  }

  reAij = 2.2204460492503131E-16 * anorm;
  b_atol = 2.2250738585072014E-308;
  if (reAij > 2.2250738585072014E-308) {
    b_atol = reAij;
  }

  reAij = 2.2250738585072014E-308;
  if (anorm > 2.2250738585072014E-308) {
    reAij = anorm;
  }

  ascale = 1.0 / reAij;
  failed = TRUE;
  j = ihi + 1;
  while (j < 3) {
    alpha1[1] = A[3];
    j = 3;
  }

  guard1 = FALSE;
  guard2 = FALSE;
  if (ihi >= ilo) {
    ifirst = ilo;
    istart = ilo;
    ilast = ihi - 1;
    ilastm1 = ihi - 2;
    iiter = 0;
    goto60 = FALSE;
    goto70 = FALSE;
    goto90 = FALSE;
    jiter = 1;
    do {
      exitg1 = 0;
      if (jiter <= 30 * ((ihi - ilo) + 1)) {
        if (ilast + 1 == ilo) {
          goto60 = TRUE;
        } else if (fabs(A[ilast + (ilastm1 << 1)].re) + fabs(A[ilast + (ilastm1 <<
          1)].im) <= b_atol) {
          A[ilast + (ilastm1 << 1)].re = 0.0;
          A[ilast + (ilastm1 << 1)].im = 0.0;
          goto60 = TRUE;
        } else {
          j = ilastm1;
          exitg3 = FALSE;
          while ((exitg3 == FALSE) && (j + 1 >= ilo)) {
            if (j + 1 == ilo) {
              firstNonZero = TRUE;
            } else if (fabs(A[j].re) + fabs(A[j].im) <= b_atol) {
              A[j].re = 0.0;
              A[j].im = 0.0;
              firstNonZero = TRUE;
            } else {
              firstNonZero = FALSE;
            }

            if (firstNonZero) {
              ifirst = j + 1;
              goto70 = TRUE;
              exitg3 = TRUE;
            } else {
              j--;
            }
          }
        }

        if (goto60 || goto70) {
          firstNonZero = TRUE;
        } else {
          firstNonZero = FALSE;
        }

        if (!firstNonZero) {
          for (i = 0; i < 2; i++) {
            alpha1[i].re = rtNaN;
            alpha1[i].im = 0.0;
            beta1[i].re = rtNaN;
            beta1[i].im = 0.0;
          }

          for (istart = 0; istart < 2; istart++) {
            for (ifirst = 0; ifirst < 2; ifirst++) {
              Z[ifirst + (istart << 1)].re = rtNaN;
              Z[ifirst + (istart << 1)].im = 0.0;
            }
          }

          *info = -1.0;
          exitg1 = 1;
        } else {
          b_guard1 = FALSE;
          if (goto60) {
            goto60 = FALSE;
            alpha1[ilast] = A[ilast + (ilast << 1)];
            ilast = ilastm1;
            ilastm1--;
            if (ilast + 1 < ilo) {
              failed = FALSE;
              guard2 = TRUE;
              exitg1 = 1;
            } else {
              iiter = 0;
              eshift_re = 0.0;
              eshift_im = 0.0;
              b_guard1 = TRUE;
            }
          } else {
            if (goto70) {
              goto70 = FALSE;
              iiter++;
              if (iiter - iiter / 10 * 10 != 0) {
                b_A.re = -(A[ilast + (ilast << 1)].re - A[ilastm1 + (ilastm1 <<
                            1)].re);
                b_A.im = -(A[ilast + (ilast << 1)].im - A[ilastm1 + (ilastm1 <<
                            1)].im);
                t1 = b_eml_div(b_A, 2.0);
                anorm = A[ilastm1 + (ilast << 1)].re * A[ilast + (ilastm1 << 1)]
                  .re - A[ilastm1 + (ilast << 1)].im * A[ilast + (ilastm1 << 1)]
                  .im;
                reAij = A[ilastm1 + (ilast << 1)].re * A[ilast + (ilastm1 << 1)]
                  .im + A[ilastm1 + (ilast << 1)].im * A[ilast + (ilastm1 << 1)]
                  .re;
                d.re = (t1.re * t1.re - t1.im * t1.im) + anorm;
                d.im = (t1.re * t1.im + t1.im * t1.re) + reAij;
                b_sqrt(&d);
                temp2 = A[ilastm1 + (ilastm1 << 1)].re - (t1.re - d.re);
                sigma1_im = A[ilastm1 + (ilastm1 << 1)].im - (t1.im - d.im);
                scale = A[ilastm1 + (ilastm1 << 1)].re - (t1.re + d.re);
                sumsq = A[ilastm1 + (ilastm1 << 1)].im - (t1.im + d.im);
                rho_re = temp2 - A[ilast + (ilast << 1)].re;
                rho_im = sigma1_im - A[ilast + (ilast << 1)].im;
                anorm = scale - A[ilast + (ilast << 1)].re;
                reAij = sumsq - A[ilast + (ilast << 1)].im;
                if (rt_hypotd_snf(fabs(rho_re), fabs(rho_im)) <= rt_hypotd_snf
                    (fabs(anorm), fabs(reAij))) {
                  scale = temp2;
                  sumsq = sigma1_im;
                  rho_re = t1.re - d.re;
                  rho_im = t1.im - d.im;
                } else {
                  rho_re = t1.re + d.re;
                  rho_im = t1.im + d.im;
                }
              } else {
                eshift_re += A[ilast + (ilastm1 << 1)].re;
                eshift_im += A[ilast + (ilastm1 << 1)].im;
                scale = eshift_re;
                sumsq = eshift_im;
              }

              j = ilastm1;
              exitg2 = FALSE;
              while ((exitg2 == FALSE) && (j + 1 > ifirst)) {
                istart = 2;
                ctemp.re = A[3].re - scale;
                ctemp.im = A[3].im - sumsq;
                anorm = ascale * (fabs(ctemp.re) + fabs(ctemp.im));
                temp2 = ascale * (fabs(A[3].re) + fabs(A[3].im));
                reAij = anorm;
                if (temp2 > anorm) {
                  reAij = temp2;
                }

                if ((reAij < 1.0) && (reAij != 0.0)) {
                  anorm /= reAij;
                  temp2 /= reAij;
                }

                if ((fabs(A[1].re) + fabs(A[1].im)) * temp2 <= anorm * b_atol) {
                  goto90 = TRUE;
                  exitg2 = TRUE;
                } else {
                  j = 0;
                }
              }

              if (!goto90) {
                istart = ifirst;
                if (ifirst == ilastm1 + 1) {
                  ctemp.re = rho_re;
                  ctemp.im = rho_im;
                } else {
                  ctemp.re = A[(ifirst + ((ifirst - 1) << 1)) - 1].re - scale;
                  ctemp.im = A[(ifirst + ((ifirst - 1) << 1)) - 1].im - sumsq;
                }

                goto90 = TRUE;
              }
            }

            if (goto90) {
              goto90 = FALSE;
              c_A = A[1 + ((istart - 1) << 1)];
              b_eml_matlab_zlartg(ctemp, c_A, &imAij, &t1);
              j = istart;
              while (j < ilast + 1) {
                for (j = 0; j < 2; j++) {
                  d.re = imAij * A[j << 1].re;
                  d.im = imAij * A[j << 1].im;
                  temp2 = t1.re * A[1 + (j << 1)].re - t1.im * A[1 + (j << 1)].
                    im;
                  sigma1_im = t1.re * A[1 + (j << 1)].im + t1.im * A[1 + (j << 1)]
                    .re;
                  anorm = A[j << 1].re;
                  reAij = A[j << 1].im;
                  scale = A[j << 1].im;
                  sumsq = A[j << 1].re;
                  A[1 + (j << 1)].re = imAij * A[1 + (j << 1)].re - (t1.re *
                    anorm + t1.im * reAij);
                  A[1 + (j << 1)].im = imAij * A[1 + (j << 1)].im - (t1.re *
                    scale - t1.im * sumsq);
                  A[j << 1].re = d.re + temp2;
                  A[j << 1].im = d.im + sigma1_im;
                }

                t1.re = -t1.re;
                t1.im = -t1.im;
                for (i = 0; i < 2; i++) {
                  d.re = imAij * A[2 + i].re;
                  d.im = imAij * A[2 + i].im;
                  temp2 = t1.re * A[i].re - t1.im * A[i].im;
                  sigma1_im = t1.re * A[i].im + t1.im * A[i].re;
                  anorm = A[2 + i].re;
                  reAij = A[2 + i].im;
                  scale = A[2 + i].im;
                  sumsq = A[2 + i].re;
                  A[i].re = imAij * A[i].re - (t1.re * anorm + t1.im * reAij);
                  A[i].im = imAij * A[i].im - (t1.re * scale - t1.im * sumsq);
                  A[2 + i].re = d.re + temp2;
                  A[2 + i].im = d.im + sigma1_im;
                }

                for (i = 0; i < 2; i++) {
                  d.re = imAij * Z[2 + i].re;
                  d.im = imAij * Z[2 + i].im;
                  temp2 = t1.re * Z[i].re - t1.im * Z[i].im;
                  sigma1_im = t1.re * Z[i].im + t1.im * Z[i].re;
                  anorm = Z[2 + i].re;
                  reAij = Z[2 + i].im;
                  scale = Z[2 + i].im;
                  sumsq = Z[2 + i].re;
                  Z[i].re = imAij * Z[i].re - (t1.re * anorm + t1.im * reAij);
                  Z[i].im = imAij * Z[i].im - (t1.re * scale - t1.im * sumsq);
                  Z[2 + i].re = d.re + temp2;
                  Z[2 + i].im = d.im + sigma1_im;
                }

                j = 2;
              }
            }

            b_guard1 = TRUE;
          }

          if (b_guard1 == TRUE) {
            jiter++;
          }
        }
      } else {
        guard2 = TRUE;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  } else {
    guard1 = TRUE;
  }

  if (guard2 == TRUE) {
    if (failed) {
      *info = (real_T)(ilast + 1);
      for (ifirst = 0; ifirst + 1 <= ilast + 1; ifirst++) {
        alpha1[ifirst].re = rtNaN;
        alpha1[ifirst].im = 0.0;
        beta1[ifirst].re = rtNaN;
        beta1[ifirst].im = 0.0;
      }

      for (istart = 0; istart < 2; istart++) {
        for (ifirst = 0; ifirst < 2; ifirst++) {
          Z[ifirst + (istart << 1)].re = rtNaN;
          Z[ifirst + (istart << 1)].im = 0.0;
        }
      }
    } else {
      guard1 = TRUE;
    }
  }

  if (guard1 == TRUE) {
    for (j = 0; j + 1 < ilo; j++) {
      alpha1[j] = A[j + (j << 1)];
    }
  }
}

static void b_eml_matlab_zlartg(const creal_T f, const creal_T g, real_T *cs,
  creal_T *sn)
{
  real_T scale;
  real_T f2s;
  real_T g2;
  real_T fs_re;
  real_T fs_im;
  real_T gs_re;
  real_T gs_im;
  boolean_T guard1 = FALSE;
  real_T g2s;
  scale = fabs(f.re);
  f2s = fabs(f.im);
  if (f2s > scale) {
    scale = f2s;
  }

  f2s = fabs(g.re);
  g2 = fabs(g.im);
  if (g2 > f2s) {
    f2s = g2;
  }

  if (f2s > scale) {
    scale = f2s;
  }

  fs_re = f.re;
  fs_im = f.im;
  gs_re = g.re;
  gs_im = g.im;
  guard1 = FALSE;
  if (scale >= 7.4428285367870146E+137) {
    do {
      fs_re *= 1.3435752215134178E-138;
      fs_im *= 1.3435752215134178E-138;
      gs_re *= 1.3435752215134178E-138;
      gs_im *= 1.3435752215134178E-138;
      scale *= 1.3435752215134178E-138;
    } while (!(scale < 7.4428285367870146E+137));

    guard1 = TRUE;
  } else if (scale <= 1.3435752215134178E-138) {
    if ((g.re == 0.0) && (g.im == 0.0)) {
      *cs = 1.0;
      sn->re = 0.0;
      sn->im = 0.0;
    } else {
      do {
        fs_re *= 7.4428285367870146E+137;
        fs_im *= 7.4428285367870146E+137;
        gs_re *= 7.4428285367870146E+137;
        gs_im *= 7.4428285367870146E+137;
        scale *= 7.4428285367870146E+137;
      } while (!(scale > 1.3435752215134178E-138));

      guard1 = TRUE;
    }
  } else {
    guard1 = TRUE;
  }

  if (guard1 == TRUE) {
    scale = fs_re * fs_re + fs_im * fs_im;
    g2 = gs_re * gs_re + gs_im * gs_im;
    f2s = g2;
    if (1.0 > g2) {
      f2s = 1.0;
    }

    if (scale <= f2s * 2.0041683600089728E-292) {
      if ((f.re == 0.0) && (f.im == 0.0)) {
        *cs = 0.0;
        scale = rt_hypotd_snf(fabs(gs_re), fabs(gs_im));
        sn->re = gs_re / scale;
        sn->im = -gs_im / scale;
      } else {
        g2s = sqrt(g2);
        *cs = rt_hypotd_snf(fabs(fs_re), fabs(fs_im)) / g2s;
        f2s = fabs(f.re);
        g2 = fabs(f.im);
        if (g2 > f2s) {
          f2s = g2;
        }

        if (f2s > 1.0) {
          scale = rt_hypotd_snf(fabs(f.re), fabs(f.im));
          fs_re = f.re / scale;
          fs_im = f.im / scale;
        } else {
          f2s = 7.4428285367870146E+137 * f.re;
          g2 = 7.4428285367870146E+137 * f.im;
          scale = rt_hypotd_snf(fabs(f2s), fabs(g2));
          fs_re = f2s / scale;
          fs_im = g2 / scale;
        }

        gs_re /= g2s;
        gs_im = -gs_im / g2s;
        sn->re = fs_re * gs_re - fs_im * gs_im;
        sn->im = fs_re * gs_im + fs_im * gs_re;
      }
    } else {
      f2s = sqrt(1.0 + g2 / scale);
      *cs = 1.0 / f2s;
      scale += g2;
      fs_re = f2s * fs_re / scale;
      fs_im = f2s * fs_im / scale;
      sn->re = fs_re * gs_re - fs_im * -gs_im;
      sn->im = fs_re * -gs_im + fs_im * gs_re;
    }
  }
}

static void b_eml_matlab_ztgevc(const creal_T A[4], creal_T V[4])
{
  creal_T work1[2];
  creal_T work2[2];
  int32_T i;
  real_T y;
  real_T anorm;
  real_T dmin;
  real_T ascale;
  int32_T je;
  real_T temp;
  real_T salpha_re;
  real_T salpha_im;
  real_T acoeff;
  boolean_T b2;
  boolean_T b3;
  real_T scale;
  real_T b_y;
  creal_T d;
  creal_T b_work1;
  int32_T jc;
  for (i = 0; i < 2; i++) {
    work1[i].re = 0.0;
    work1[i].im = 0.0;
    work2[i].re = 0.0;
    work2[i].im = 0.0;
  }

  y = (fabs(A[2].re) + fabs(A[2].im)) + (fabs(A[3].re) + fabs(A[3].im));
  anorm = fabs(A[0].re) + fabs(A[0].im);
  if (y > anorm) {
    anorm = y;
  }

  dmin = anorm;
  if (2.2250738585072014E-308 > anorm) {
    dmin = 2.2250738585072014E-308;
  }

  ascale = 1.0 / dmin;
  for (je = 0; je < 2; je++) {
    y = (fabs(A[(((1 - je) << 1) - je) + 1].re) + fabs(A[(((1 - je) << 1) - je)
          + 1].im)) * ascale;
    if (1.0 > y) {
      y = 1.0;
    }

    temp = 1.0 / y;
    salpha_re = ascale * (temp * A[(((1 - je) << 1) - je) + 1].re);
    salpha_im = ascale * (temp * A[(((1 - je) << 1) - je) + 1].im);
    acoeff = temp * ascale;
    if ((fabs(temp) >= 2.2250738585072014E-308) && (fabs(acoeff) <
         2.0041683600089728E-292)) {
      b2 = TRUE;
    } else {
      b2 = FALSE;
    }

    if ((fabs(salpha_re) + fabs(salpha_im) >= 2.2250738585072014E-308) && (fabs
         (salpha_re) + fabs(salpha_im) < 2.0041683600089728E-292)) {
      b3 = TRUE;
    } else {
      b3 = FALSE;
    }

    scale = 1.0;
    if (b2) {
      dmin = anorm;
      if (4.9896007738368E+291 < anorm) {
        dmin = 4.9896007738368E+291;
      }

      scale = 2.0041683600089728E-292 / fabs(temp) * dmin;
    }

    if (b3) {
      dmin = 2.0041683600089728E-292 / (fabs(salpha_re) + fabs(salpha_im));
      if (dmin > scale) {
        scale = dmin;
      }
    }

    if (b2 || b3) {
      y = fabs(acoeff);
      b_y = fabs(salpha_re) + fabs(salpha_im);
      if (1.0 > y) {
        y = 1.0;
      }

      if (b_y > y) {
        y = b_y;
      }

      y = 1.0 / (2.2250738585072014E-308 * y);
      if (y < scale) {
        scale = y;
      }

      if (b2) {
        acoeff = ascale * (scale * temp);
      } else {
        acoeff *= scale;
      }

      if (b3) {
        salpha_re *= scale;
        salpha_im *= scale;
      } else {
        salpha_re *= scale;
        salpha_im *= scale;
      }
    }

    for (i = 0; i < 2; i++) {
      work1[i].re = 0.0;
      work1[i].im = 0.0;
    }

    work1[1 - je].re = 1.0;
    work1[1 - je].im = 0.0;
    dmin = 2.2204460492503131E-16 * fabs(acoeff) * anorm;
    y = 2.2204460492503131E-16 * (fabs(salpha_re) + fabs(salpha_im));
    if (y > dmin) {
      dmin = y;
    }

    if (2.2250738585072014E-308 > dmin) {
      dmin = 2.2250738585072014E-308;
    }

    i = 0;
    while (i <= -je) {
      work1[0].re = acoeff * A[(1 - je) << 1].re;
      work1[0].im = acoeff * A[(1 - je) << 1].im;
      i = 1;
    }

    work1[1 - je].re = 1.0;
    work1[1 - je].im = 0.0;
    i = 0;
    while (i <= -je) {
      d.re = acoeff * A[0].re - salpha_re;
      d.im = acoeff * A[0].im - salpha_im;
      if (fabs(d.re) + fabs(d.im) <= dmin) {
        d.re = dmin;
        d.im = 0.0;
      }

      if ((fabs(d.re) + fabs(d.im) < 1.0) && (fabs(work1[0].re) + fabs(work1[0].
            im) >= 2.2471164185778949E+307 * (fabs(d.re) + fabs(d.im)))) {
        temp = 1.0 / (fabs(work1[0].re) + fabs(work1[0].im));
        for (i = 0; i <= 1 - je; i++) {
          work1[i].re *= temp;
          work1[i].im *= temp;
        }
      }

      b_work1.re = -work1[0].re;
      b_work1.im = -work1[0].im;
      work1[0] = c_eml_div(b_work1, d);
      i = 1;
    }

    for (i = 0; i < 2; i++) {
      work2[i].re = 0.0;
      work2[i].im = 0.0;
    }

    for (jc = 0; jc <= 1 - je; jc++) {
      for (i = 0; i < 2; i++) {
        dmin = V[i + (jc << 1)].re * work1[jc].re - V[i + (jc << 1)].im *
          work1[jc].im;
        b_y = V[i + (jc << 1)].re * work1[jc].im + V[i + (jc << 1)].im *
          work1[jc].re;
        work2[i].re += dmin;
        work2[i].im += b_y;
      }
    }

    y = fabs(work2[1].re) + fabs(work2[1].im);
    dmin = fabs(work2[0].re) + fabs(work2[0].im);
    if (y > dmin) {
      dmin = y;
    }

    if (dmin > 2.2250738585072014E-308) {
      temp = 1.0 / dmin;
      for (i = 0; i < 2; i++) {
        V[i + ((1 - je) << 1)].re = temp * work2[i].re;
        V[i + ((1 - je) << 1)].im = temp * work2[i].im;
      }
    } else {
      for (i = 0; i < 2; i++) {
        V[i + ((1 - je) << 1)].re = 0.0;
        V[i + ((1 - je) << 1)].im = 0.0;
      }
    }
  }
}

static void eml_matlab_zggev(creal_T A[9], real_T *info, creal_T alpha1[3],
  creal_T beta1[3], creal_T V[9])
{
  real_T anrm;
  int32_T nzcount;
  boolean_T exitg7;
  real_T absxk;
  int32_T i;
  int32_T ii;
  boolean_T ilascl;
  real_T anrmto;
  real_T ctoc;
  boolean_T notdone;
  real_T cfrom1;
  real_T cto1;
  real_T mul;
  creal_T b_A[9];
  int32_T rscale[3];
  int32_T ilo;
  int32_T ihi;
  int32_T exitg2;
  int32_T j;
  boolean_T exitg5;
  int32_T jj;
  boolean_T exitg6;
  boolean_T guard2 = FALSE;
  creal_T atmp;
  int32_T exitg1;
  boolean_T exitg3;
  boolean_T exitg4;
  boolean_T guard1 = FALSE;
  static const int8_T iv4[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };

  creal_T c_A;
  creal_T d_A;
  real_T y_re;
  real_T y_im;
  *info = 0.0;
  anrm = 0.0;
  nzcount = 0;
  exitg7 = FALSE;
  while ((exitg7 == FALSE) && (nzcount < 9)) {
    absxk = rt_hypotd_snf(fabs(A[nzcount].re), fabs(A[nzcount].im));
    if (rtIsNaN(absxk)) {
      anrm = rtNaN;
      exitg7 = TRUE;
    } else {
      if (absxk > anrm) {
        anrm = absxk;
      }

      nzcount++;
    }
  }

  if (!((!rtIsInf(anrm)) && (!rtIsNaN(anrm)))) {
    for (i = 0; i < 3; i++) {
      alpha1[i].re = rtNaN;
      alpha1[i].im = 0.0;
      beta1[i].re = rtNaN;
      beta1[i].im = 0.0;
    }

    for (ii = 0; ii < 9; ii++) {
      V[ii].re = rtNaN;
      V[ii].im = 0.0;
    }
  } else {
    ilascl = FALSE;
    anrmto = anrm;
    if ((anrm > 0.0) && (anrm < 6.7178761075670888E-139)) {
      anrmto = 6.7178761075670888E-139;
      ilascl = TRUE;
    } else {
      if (anrm > 1.4885657073574029E+138) {
        anrmto = 1.4885657073574029E+138;
        ilascl = TRUE;
      }
    }

    if (ilascl) {
      absxk = anrm;
      ctoc = anrmto;
      notdone = TRUE;
      while (notdone) {
        cfrom1 = absxk * 2.0041683600089728E-292;
        cto1 = ctoc / 4.9896007738368E+291;
        if ((cfrom1 > ctoc) && (ctoc != 0.0)) {
          mul = 2.0041683600089728E-292;
          absxk = cfrom1;
        } else if (cto1 > absxk) {
          mul = 4.9896007738368E+291;
          ctoc = cto1;
        } else {
          mul = ctoc / absxk;
          notdone = FALSE;
        }

        for (ii = 0; ii < 9; ii++) {
          A[ii].re *= mul;
          A[ii].im *= mul;
        }
      }
    }

    memcpy(&b_A[0], &A[0], 9U * sizeof(creal_T));
    for (i = 0; i < 3; i++) {
      rscale[i] = 0;
    }

    ilo = 1;
    ihi = 3;
    do {
      exitg2 = 0;
      i = 0;
      j = 0;
      notdone = FALSE;
      ii = ihi;
      exitg5 = FALSE;
      while ((exitg5 == FALSE) && (ii > 0)) {
        nzcount = 0;
        i = ii;
        j = ihi;
        jj = 1;
        exitg6 = FALSE;
        while ((exitg6 == FALSE) && (jj <= ihi)) {
          guard2 = FALSE;
          if ((b_A[(ii + 3 * (jj - 1)) - 1].re != 0.0) || (b_A[(ii + 3 * (jj - 1))
               - 1].im != 0.0) || (ii == jj)) {
            if (nzcount == 0) {
              j = jj;
              nzcount = 1;
              guard2 = TRUE;
            } else {
              nzcount = 2;
              exitg6 = TRUE;
            }
          } else {
            guard2 = TRUE;
          }

          if (guard2 == TRUE) {
            jj++;
          }
        }

        if (nzcount < 2) {
          notdone = TRUE;
          exitg5 = TRUE;
        } else {
          ii--;
        }
      }

      if (!notdone) {
        exitg2 = 2;
      } else {
        if (i != ihi) {
          for (nzcount = 0; nzcount < 3; nzcount++) {
            atmp = b_A[(i + 3 * nzcount) - 1];
            b_A[(i + 3 * nzcount) - 1] = b_A[(ihi + 3 * nzcount) - 1];
            b_A[(ihi + 3 * nzcount) - 1] = atmp;
          }
        }

        if (j != ihi) {
          for (nzcount = 0; nzcount + 1 <= ihi; nzcount++) {
            atmp = b_A[nzcount + 3 * (j - 1)];
            b_A[nzcount + 3 * (j - 1)] = b_A[nzcount + 3 * (ihi - 1)];
            b_A[nzcount + 3 * (ihi - 1)] = atmp;
          }
        }

        rscale[ihi - 1] = j;
        ihi--;
        if (ihi == 1) {
          rscale[0] = 1;
          exitg2 = 1;
        }
      }
    } while (exitg2 == 0);

    if (exitg2 == 1) {
    } else {
      do {
        exitg1 = 0;
        i = 0;
        j = 0;
        notdone = FALSE;
        jj = ilo;
        exitg3 = FALSE;
        while ((exitg3 == FALSE) && (jj <= ihi)) {
          nzcount = 0;
          i = ihi;
          j = jj;
          ii = ilo;
          exitg4 = FALSE;
          while ((exitg4 == FALSE) && (ii <= ihi)) {
            guard1 = FALSE;
            if ((b_A[(ii + 3 * (jj - 1)) - 1].re != 0.0) || (b_A[(ii + 3 * (jj -
                   1)) - 1].im != 0.0) || (ii == jj)) {
              if (nzcount == 0) {
                i = ii;
                nzcount = 1;
                guard1 = TRUE;
              } else {
                nzcount = 2;
                exitg4 = TRUE;
              }
            } else {
              guard1 = TRUE;
            }

            if (guard1 == TRUE) {
              ii++;
            }
          }

          if (nzcount < 2) {
            notdone = TRUE;
            exitg3 = TRUE;
          } else {
            jj++;
          }
        }

        if (!notdone) {
          exitg1 = 1;
        } else {
          if (i != ilo) {
            for (nzcount = ilo - 1; nzcount + 1 < 4; nzcount++) {
              atmp = b_A[(i + 3 * nzcount) - 1];
              b_A[(i + 3 * nzcount) - 1] = b_A[(ilo + 3 * nzcount) - 1];
              b_A[(ilo + 3 * nzcount) - 1] = atmp;
            }
          }

          if (j != ilo) {
            for (nzcount = 0; nzcount + 1 <= ihi; nzcount++) {
              atmp = b_A[nzcount + 3 * (j - 1)];
              b_A[nzcount + 3 * (j - 1)] = b_A[nzcount + 3 * (ilo - 1)];
              b_A[nzcount + 3 * (ilo - 1)] = atmp;
            }
          }

          rscale[ilo - 1] = j;
          ilo++;
          if (ilo == ihi) {
            rscale[ilo - 1] = ilo;
            exitg1 = 1;
          }
        }
      } while (exitg1 == 0);
    }

    for (ii = 0; ii < 9; ii++) {
      V[ii].re = (real_T)iv4[ii];
      V[ii].im = 0.0;
    }

    if (ihi < ilo + 2) {
    } else {
      ii = ilo;
      while (ii < 2) {
        c_A = b_A[1];
        d_A = b_A[2];
        eml_matlab_zlartg(c_A, d_A, &cfrom1, &atmp, &b_A[1]);
        b_A[2].re = 0.0;
        b_A[2].im = 0.0;
        for (j = 0; j < 2; j++) {
          cto1 = cfrom1 * b_A[1 + 3 * (j + 1)].re;
          mul = cfrom1 * b_A[1 + 3 * (j + 1)].im;
          y_re = atmp.re * b_A[2 + 3 * (j + 1)].re - atmp.im * b_A[2 + 3 * (j +
            1)].im;
          y_im = atmp.re * b_A[2 + 3 * (j + 1)].im + atmp.im * b_A[2 + 3 * (j +
            1)].re;
          absxk = b_A[1 + 3 * (j + 1)].im;
          ctoc = b_A[1 + 3 * (j + 1)].re;
          b_A[2 + 3 * (j + 1)].re = cfrom1 * b_A[2 + 3 * (j + 1)].re - (atmp.re *
            b_A[1 + 3 * (j + 1)].re + atmp.im * b_A[1 + 3 * (j + 1)].im);
          b_A[2 + 3 * (j + 1)].im = cfrom1 * b_A[2 + 3 * (j + 1)].im - (atmp.re *
            absxk - atmp.im * ctoc);
          b_A[1 + 3 * (j + 1)].re = cto1 + y_re;
          b_A[1 + 3 * (j + 1)].im = mul + y_im;
        }

        atmp.re = -atmp.re;
        atmp.im = -atmp.im;
        for (i = ilo - 1; i + 1 < 4; i++) {
          cto1 = cfrom1 * b_A[6 + i].re;
          mul = cfrom1 * b_A[6 + i].im;
          y_re = atmp.re * b_A[3 + i].re - atmp.im * b_A[3 + i].im;
          y_im = atmp.re * b_A[3 + i].im + atmp.im * b_A[3 + i].re;
          absxk = b_A[6 + i].im;
          ctoc = b_A[6 + i].re;
          b_A[3 + i].re = cfrom1 * b_A[3 + i].re - (atmp.re * b_A[6 + i].re +
            atmp.im * b_A[6 + i].im);
          b_A[3 + i].im = cfrom1 * b_A[3 + i].im - (atmp.re * absxk - atmp.im *
            ctoc);
          b_A[6 + i].re = cto1 + y_re;
          b_A[6 + i].im = mul + y_im;
        }

        for (i = 0; i < 3; i++) {
          cto1 = cfrom1 * V[6 + i].re;
          mul = cfrom1 * V[6 + i].im;
          y_re = atmp.re * V[3 + i].re - atmp.im * V[3 + i].im;
          y_im = atmp.re * V[3 + i].im + atmp.im * V[3 + i].re;
          absxk = V[6 + i].im;
          ctoc = V[6 + i].re;
          V[3 + i].re = cfrom1 * V[3 + i].re - (atmp.re * V[6 + i].re + atmp.im *
            V[6 + i].im);
          V[3 + i].im = cfrom1 * V[3 + i].im - (atmp.re * absxk - atmp.im * ctoc);
          V[6 + i].re = cto1 + y_re;
          V[6 + i].im = mul + y_im;
        }

        ii = 2;
      }
    }

    eml_matlab_zhgeqz(b_A, ilo, ihi, V, info, alpha1, beta1);
    if (*info != 0.0) {
    } else {
      eml_matlab_ztgevc(b_A, V);
      if (ilo > 1) {
        for (i = ilo - 2; i + 1 >= 1; i--) {
          if (rscale[i] != i + 1) {
            for (j = 0; j < 3; j++) {
              atmp = V[i + 3 * j];
              V[i + 3 * j] = V[(rscale[i] + 3 * j) - 1];
              V[(rscale[i] + 3 * j) - 1] = atmp;
            }
          }
        }
      }

      if (ihi < 3) {
        while (ihi + 1 < 4) {
          if (rscale[ihi] != ihi + 1) {
            for (j = 0; j < 3; j++) {
              atmp = V[ihi + 3 * j];
              V[ihi + 3 * j] = V[(rscale[ihi] + 3 * j) - 1];
              V[(rscale[ihi] + 3 * j) - 1] = atmp;
            }
          }

          ihi++;
        }
      }

      for (ii = 0; ii < 3; ii++) {
        absxk = fabs(V[3 * ii].re) + fabs(V[3 * ii].im);
        for (nzcount = 0; nzcount < 2; nzcount++) {
          ctoc = fabs(V[(nzcount + 3 * ii) + 1].re) + fabs(V[(nzcount + 3 * ii)
            + 1].im);
          if (ctoc > absxk) {
            absxk = ctoc;
          }
        }

        if (absxk >= 6.7178761075670888E-139) {
          absxk = 1.0 / absxk;
          for (nzcount = 0; nzcount < 3; nzcount++) {
            V[nzcount + 3 * ii].re *= absxk;
            V[nzcount + 3 * ii].im *= absxk;
          }
        }
      }

      if (ilascl) {
        notdone = TRUE;
        while (notdone) {
          cfrom1 = anrmto * 2.0041683600089728E-292;
          cto1 = anrm / 4.9896007738368E+291;
          if ((cfrom1 > anrm) && (anrm != 0.0)) {
            mul = 2.0041683600089728E-292;
            anrmto = cfrom1;
          } else if (cto1 > anrmto) {
            mul = 4.9896007738368E+291;
            anrm = cto1;
          } else {
            mul = anrm / anrmto;
            notdone = FALSE;
          }

          for (ii = 0; ii < 3; ii++) {
            alpha1[ii].re *= mul;
            alpha1[ii].im *= mul;
          }
        }
      }
    }
  }
}

static void eml_matlab_zhgeqz(creal_T A[9], int32_T ilo, int32_T ihi, creal_T Z
  [9], real_T *info, creal_T alpha1[3], creal_T beta1[3])
{
  int32_T i;
  real_T eshift_re;
  real_T eshift_im;
  creal_T ctemp;
  real_T rho_re;
  real_T rho_im;
  real_T anorm;
  real_T scale;
  real_T sumsq;
  boolean_T firstNonZero;
  int32_T j;
  int32_T jm1;
  real_T reAij;
  real_T imAij;
  real_T temp2;
  real_T b_atol;
  real_T ascale;
  boolean_T failed;
  boolean_T guard1 = FALSE;
  boolean_T guard2 = FALSE;
  int32_T ifirst;
  int32_T istart;
  int32_T ilast;
  int32_T ilastm1;
  int32_T iiter;
  boolean_T goto60;
  boolean_T goto70;
  boolean_T goto90;
  int32_T jiter;
  int32_T exitg1;
  boolean_T exitg3;
  int32_T jp1;
  boolean_T b_guard1 = FALSE;
  creal_T b_A;
  creal_T t1;
  creal_T d;
  real_T sigma1_im;
  boolean_T exitg2;
  creal_T c_A;
  creal_T d_A;
  creal_T dc0;
  *info = 0.0;
  for (i = 0; i < 3; i++) {
    alpha1[i].re = 0.0;
    alpha1[i].im = 0.0;
    beta1[i].re = 1.0;
    beta1[i].im = 0.0;
  }

  eshift_re = 0.0;
  eshift_im = 0.0;
  ctemp.re = 0.0;
  ctemp.im = 0.0;
  rho_re = 0.0;
  rho_im = 0.0;
  anorm = 0.0;
  if (ilo > ihi) {
  } else {
    scale = 0.0;
    sumsq = 0.0;
    firstNonZero = TRUE;
    for (j = ilo; j <= ihi; j++) {
      jm1 = j + 1;
      if (ihi < j + 1) {
        jm1 = ihi;
      }

      for (i = ilo; i <= jm1; i++) {
        reAij = A[(i + 3 * (j - 1)) - 1].re;
        imAij = A[(i + 3 * (j - 1)) - 1].im;
        if (reAij != 0.0) {
          anorm = fabs(reAij);
          if (firstNonZero) {
            sumsq = 1.0;
            scale = anorm;
            firstNonZero = FALSE;
          } else if (scale < anorm) {
            temp2 = scale / anorm;
            sumsq = 1.0 + sumsq * temp2 * temp2;
            scale = anorm;
          } else {
            temp2 = anorm / scale;
            sumsq += temp2 * temp2;
          }
        }

        if (imAij != 0.0) {
          anorm = fabs(imAij);
          if (firstNonZero) {
            sumsq = 1.0;
            scale = anorm;
            firstNonZero = FALSE;
          } else if (scale < anorm) {
            temp2 = scale / anorm;
            sumsq = 1.0 + sumsq * temp2 * temp2;
            scale = anorm;
          } else {
            temp2 = anorm / scale;
            sumsq += temp2 * temp2;
          }
        }
      }
    }

    anorm = scale * sqrt(sumsq);
  }

  reAij = 2.2204460492503131E-16 * anorm;
  b_atol = 2.2250738585072014E-308;
  if (reAij > 2.2250738585072014E-308) {
    b_atol = reAij;
  }

  reAij = 2.2250738585072014E-308;
  if (anorm > 2.2250738585072014E-308) {
    reAij = anorm;
  }

  ascale = 1.0 / reAij;
  failed = TRUE;
  for (j = ihi; j + 1 < 4; j++) {
    alpha1[j] = A[j + 3 * j];
  }

  guard1 = FALSE;
  guard2 = FALSE;
  if (ihi >= ilo) {
    ifirst = ilo;
    istart = ilo;
    ilast = ihi - 1;
    ilastm1 = ihi - 2;
    iiter = 0;
    goto60 = FALSE;
    goto70 = FALSE;
    goto90 = FALSE;
    jiter = 1;
    do {
      exitg1 = 0;
      if (jiter <= 30 * ((ihi - ilo) + 1)) {
        if (ilast + 1 == ilo) {
          goto60 = TRUE;
        } else if (fabs(A[ilast + 3 * ilastm1].re) + fabs(A[ilast + 3 * ilastm1]
                    .im) <= b_atol) {
          A[ilast + 3 * ilastm1].re = 0.0;
          A[ilast + 3 * ilastm1].im = 0.0;
          goto60 = TRUE;
        } else {
          j = ilastm1;
          exitg3 = FALSE;
          while ((exitg3 == FALSE) && (j + 1 >= ilo)) {
            if (j + 1 == ilo) {
              firstNonZero = TRUE;
            } else if (fabs(A[j + 3 * (j - 1)].re) + fabs(A[j + 3 * (j - 1)].im)
                       <= b_atol) {
              A[j + 3 * (j - 1)].re = 0.0;
              A[j + 3 * (j - 1)].im = 0.0;
              firstNonZero = TRUE;
            } else {
              firstNonZero = FALSE;
            }

            if (firstNonZero) {
              ifirst = j + 1;
              goto70 = TRUE;
              exitg3 = TRUE;
            } else {
              j--;
            }
          }
        }

        if (goto60 || goto70) {
          firstNonZero = TRUE;
        } else {
          firstNonZero = FALSE;
        }

        if (!firstNonZero) {
          for (i = 0; i < 3; i++) {
            alpha1[i].re = rtNaN;
            alpha1[i].im = 0.0;
            beta1[i].re = rtNaN;
            beta1[i].im = 0.0;
          }

          for (jm1 = 0; jm1 < 3; jm1++) {
            for (jp1 = 0; jp1 < 3; jp1++) {
              Z[jp1 + 3 * jm1].re = rtNaN;
              Z[jp1 + 3 * jm1].im = 0.0;
            }
          }

          *info = -1.0;
          exitg1 = 1;
        } else {
          b_guard1 = FALSE;
          if (goto60) {
            goto60 = FALSE;
            alpha1[ilast] = A[ilast + 3 * ilast];
            ilast = ilastm1;
            ilastm1--;
            if (ilast + 1 < ilo) {
              failed = FALSE;
              guard2 = TRUE;
              exitg1 = 1;
            } else {
              iiter = 0;
              eshift_re = 0.0;
              eshift_im = 0.0;
              b_guard1 = TRUE;
            }
          } else {
            if (goto70) {
              goto70 = FALSE;
              iiter++;
              if (iiter - iiter / 10 * 10 != 0) {
                b_A.re = -(A[ilast + 3 * ilast].re - A[ilastm1 + 3 * ilastm1].re);
                b_A.im = -(A[ilast + 3 * ilast].im - A[ilastm1 + 3 * ilastm1].im);
                t1 = b_eml_div(b_A, 2.0);
                anorm = A[ilastm1 + 3 * ilast].re * A[ilast + 3 * ilastm1].re -
                  A[ilastm1 + 3 * ilast].im * A[ilast + 3 * ilastm1].im;
                reAij = A[ilastm1 + 3 * ilast].re * A[ilast + 3 * ilastm1].im +
                  A[ilastm1 + 3 * ilast].im * A[ilast + 3 * ilastm1].re;
                d.re = (t1.re * t1.re - t1.im * t1.im) + anorm;
                d.im = (t1.re * t1.im + t1.im * t1.re) + reAij;
                b_sqrt(&d);
                temp2 = A[ilastm1 + 3 * ilastm1].re - (t1.re - d.re);
                sigma1_im = A[ilastm1 + 3 * ilastm1].im - (t1.im - d.im);
                scale = A[ilastm1 + 3 * ilastm1].re - (t1.re + d.re);
                sumsq = A[ilastm1 + 3 * ilastm1].im - (t1.im + d.im);
                rho_re = temp2 - A[ilast + 3 * ilast].re;
                rho_im = sigma1_im - A[ilast + 3 * ilast].im;
                anorm = scale - A[ilast + 3 * ilast].re;
                reAij = sumsq - A[ilast + 3 * ilast].im;
                if (rt_hypotd_snf(fabs(rho_re), fabs(rho_im)) <= rt_hypotd_snf
                    (fabs(anorm), fabs(reAij))) {
                  scale = temp2;
                  sumsq = sigma1_im;
                  rho_re = t1.re - d.re;
                  rho_im = t1.im - d.im;
                } else {
                  rho_re = t1.re + d.re;
                  rho_im = t1.im + d.im;
                }
              } else {
                eshift_re += A[ilast + 3 * ilastm1].re;
                eshift_im += A[ilast + 3 * ilastm1].im;
                scale = eshift_re;
                sumsq = eshift_im;
              }

              j = ilastm1;
              jp1 = ilastm1 + 1;
              exitg2 = FALSE;
              while ((exitg2 == FALSE) && (j + 1 > ifirst)) {
                istart = j + 1;
                ctemp.re = A[j + 3 * j].re - scale;
                ctemp.im = A[j + 3 * j].im - sumsq;
                anorm = ascale * (fabs(ctemp.re) + fabs(ctemp.im));
                temp2 = ascale * (fabs(A[jp1 + 3 * j].re) + fabs(A[jp1 + 3 * j].
                  im));
                reAij = anorm;
                if (temp2 > anorm) {
                  reAij = temp2;
                }

                if ((reAij < 1.0) && (reAij != 0.0)) {
                  anorm /= reAij;
                  temp2 /= reAij;
                }

                if ((fabs(A[j + 3 * (j - 1)].re) + fabs(A[j + 3 * (j - 1)].im)) *
                    temp2 <= anorm * b_atol) {
                  goto90 = TRUE;
                  exitg2 = TRUE;
                } else {
                  jp1 = j;
                  j--;
                }
              }

              if (!goto90) {
                istart = ifirst;
                if (ifirst == ilastm1 + 1) {
                  ctemp.re = rho_re;
                  ctemp.im = rho_im;
                } else {
                  ctemp.re = A[(ifirst + 3 * (ifirst - 1)) - 1].re - scale;
                  ctemp.im = A[(ifirst + 3 * (ifirst - 1)) - 1].im - sumsq;
                }

                goto90 = TRUE;
              }
            }

            if (goto90) {
              goto90 = FALSE;
              c_A = A[istart + 3 * (istart - 1)];
              b_eml_matlab_zlartg(ctemp, c_A, &imAij, &t1);
              j = istart;
              jm1 = istart - 2;
              while (j < ilast + 1) {
                if (j > istart) {
                  c_A = A[1 + 3 * jm1];
                  d_A = A[j + 3 * jm1];
                  eml_matlab_zlartg(c_A, d_A, &imAij, &t1, &dc0);
                  A[1 + 3 * jm1] = dc0;
                  A[j + 3 * jm1].re = 0.0;
                  A[j + 3 * jm1].im = 0.0;
                }

                for (jp1 = j - 1; jp1 + 1 < 4; jp1++) {
                  d.re = imAij * A[(j + 3 * jp1) - 1].re;
                  d.im = imAij * A[(j + 3 * jp1) - 1].im;
                  temp2 = t1.re * A[j + 3 * jp1].re - t1.im * A[j + 3 * jp1].im;
                  sigma1_im = t1.re * A[j + 3 * jp1].im + t1.im * A[j + 3 * jp1]
                    .re;
                  anorm = A[(j + 3 * jp1) - 1].re;
                  reAij = A[(j + 3 * jp1) - 1].im;
                  scale = A[(j + 3 * jp1) - 1].im;
                  sumsq = A[(j + 3 * jp1) - 1].re;
                  A[j + 3 * jp1].re = imAij * A[j + 3 * jp1].re - (t1.re * anorm
                    + t1.im * reAij);
                  A[j + 3 * jp1].im = imAij * A[j + 3 * jp1].im - (t1.re * scale
                    - t1.im * sumsq);
                  A[(j + 3 * jp1) - 1].re = d.re + temp2;
                  A[(j + 3 * jp1) - 1].im = d.im + sigma1_im;
                }

                t1.re = -t1.re;
                t1.im = -t1.im;
                jp1 = j;
                if (ilast + 1 < j + 2) {
                  jp1 = ilast - 1;
                }

                for (i = 0; i + 1 <= jp1 + 2; i++) {
                  d.re = imAij * A[i + 3 * j].re;
                  d.im = imAij * A[i + 3 * j].im;
                  temp2 = t1.re * A[i + 3 * (j - 1)].re - t1.im * A[i + 3 * (j -
                    1)].im;
                  sigma1_im = t1.re * A[i + 3 * (j - 1)].im + t1.im * A[i + 3 *
                    (j - 1)].re;
                  anorm = A[i + 3 * j].re;
                  reAij = A[i + 3 * j].im;
                  scale = A[i + 3 * j].im;
                  sumsq = A[i + 3 * j].re;
                  A[i + 3 * (j - 1)].re = imAij * A[i + 3 * (j - 1)].re - (t1.re
                    * anorm + t1.im * reAij);
                  A[i + 3 * (j - 1)].im = imAij * A[i + 3 * (j - 1)].im - (t1.re
                    * scale - t1.im * sumsq);
                  A[i + 3 * j].re = d.re + temp2;
                  A[i + 3 * j].im = d.im + sigma1_im;
                }

                for (i = 0; i < 3; i++) {
                  d.re = imAij * Z[i + 3 * j].re;
                  d.im = imAij * Z[i + 3 * j].im;
                  temp2 = t1.re * Z[i + 3 * (j - 1)].re - t1.im * Z[i + 3 * (j -
                    1)].im;
                  sigma1_im = t1.re * Z[i + 3 * (j - 1)].im + t1.im * Z[i + 3 *
                    (j - 1)].re;
                  anorm = Z[i + 3 * j].re;
                  reAij = Z[i + 3 * j].im;
                  scale = Z[i + 3 * j].im;
                  sumsq = Z[i + 3 * j].re;
                  Z[i + 3 * (j - 1)].re = imAij * Z[i + 3 * (j - 1)].re - (t1.re
                    * anorm + t1.im * reAij);
                  Z[i + 3 * (j - 1)].im = imAij * Z[i + 3 * (j - 1)].im - (t1.re
                    * scale - t1.im * sumsq);
                  Z[i + 3 * j].re = d.re + temp2;
                  Z[i + 3 * j].im = d.im + sigma1_im;
                }

                jm1 = j - 1;
                j++;
              }
            }

            b_guard1 = TRUE;
          }

          if (b_guard1 == TRUE) {
            jiter++;
          }
        }
      } else {
        guard2 = TRUE;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  } else {
    guard1 = TRUE;
  }

  if (guard2 == TRUE) {
    if (failed) {
      *info = (real_T)(ilast + 1);
      for (jp1 = 0; jp1 + 1 <= ilast + 1; jp1++) {
        alpha1[jp1].re = rtNaN;
        alpha1[jp1].im = 0.0;
        beta1[jp1].re = rtNaN;
        beta1[jp1].im = 0.0;
      }

      for (jm1 = 0; jm1 < 3; jm1++) {
        for (jp1 = 0; jp1 < 3; jp1++) {
          Z[jp1 + 3 * jm1].re = rtNaN;
          Z[jp1 + 3 * jm1].im = 0.0;
        }
      }
    } else {
      guard1 = TRUE;
    }
  }

  if (guard1 == TRUE) {
    for (j = 0; j + 1 < ilo; j++) {
      alpha1[j] = A[j + 3 * j];
    }
  }
}

static void eml_matlab_zlartg(const creal_T f, const creal_T g, real_T *cs,
  creal_T *sn, creal_T *r)
{
  real_T scale;
  real_T f2s;
  real_T g2;
  real_T fs_re;
  real_T fs_im;
  real_T gs_re;
  real_T gs_im;
  int32_T count;
  int32_T rescaledir;
  boolean_T guard1 = FALSE;
  real_T g2s;
  scale = fabs(f.re);
  f2s = fabs(f.im);
  if (f2s > scale) {
    scale = f2s;
  }

  f2s = fabs(g.re);
  g2 = fabs(g.im);
  if (g2 > f2s) {
    f2s = g2;
  }

  if (f2s > scale) {
    scale = f2s;
  }

  fs_re = f.re;
  fs_im = f.im;
  gs_re = g.re;
  gs_im = g.im;
  count = 0;
  rescaledir = 0;
  guard1 = FALSE;
  if (scale >= 7.4428285367870146E+137) {
    do {
      count++;
      fs_re *= 1.3435752215134178E-138;
      fs_im *= 1.3435752215134178E-138;
      gs_re *= 1.3435752215134178E-138;
      gs_im *= 1.3435752215134178E-138;
      scale *= 1.3435752215134178E-138;
    } while (!(scale < 7.4428285367870146E+137));

    rescaledir = 1;
    guard1 = TRUE;
  } else if (scale <= 1.3435752215134178E-138) {
    if ((g.re == 0.0) && (g.im == 0.0)) {
      *cs = 1.0;
      sn->re = 0.0;
      sn->im = 0.0;
      *r = f;
    } else {
      do {
        count++;
        fs_re *= 7.4428285367870146E+137;
        fs_im *= 7.4428285367870146E+137;
        gs_re *= 7.4428285367870146E+137;
        gs_im *= 7.4428285367870146E+137;
        scale *= 7.4428285367870146E+137;
      } while (!(scale > 1.3435752215134178E-138));

      rescaledir = -1;
      guard1 = TRUE;
    }
  } else {
    guard1 = TRUE;
  }

  if (guard1 == TRUE) {
    scale = fs_re * fs_re + fs_im * fs_im;
    g2 = gs_re * gs_re + gs_im * gs_im;
    f2s = g2;
    if (1.0 > g2) {
      f2s = 1.0;
    }

    if (scale <= f2s * 2.0041683600089728E-292) {
      if ((f.re == 0.0) && (f.im == 0.0)) {
        *cs = 0.0;
        r->re = rt_hypotd_snf(fabs(g.re), fabs(g.im));
        r->im = 0.0;
        f2s = rt_hypotd_snf(fabs(gs_re), fabs(gs_im));
        sn->re = gs_re / f2s;
        sn->im = -gs_im / f2s;
      } else {
        g2s = sqrt(g2);
        *cs = rt_hypotd_snf(fabs(fs_re), fabs(fs_im)) / g2s;
        f2s = fabs(f.re);
        g2 = fabs(f.im);
        if (g2 > f2s) {
          f2s = g2;
        }

        if (f2s > 1.0) {
          f2s = rt_hypotd_snf(fabs(f.re), fabs(f.im));
          fs_re = f.re / f2s;
          fs_im = f.im / f2s;
        } else {
          g2 = 7.4428285367870146E+137 * f.re;
          scale = 7.4428285367870146E+137 * f.im;
          f2s = rt_hypotd_snf(fabs(g2), fabs(scale));
          fs_re = g2 / f2s;
          fs_im = scale / f2s;
        }

        gs_re /= g2s;
        gs_im = -gs_im / g2s;
        sn->re = fs_re * gs_re - fs_im * gs_im;
        sn->im = fs_re * gs_im + fs_im * gs_re;
        r->re = *cs * f.re + (sn->re * g.re - sn->im * g.im);
        r->im = *cs * f.im + (sn->re * g.im + sn->im * g.re);
      }
    } else {
      f2s = sqrt(1.0 + g2 / scale);
      r->re = f2s * fs_re;
      r->im = f2s * fs_im;
      *cs = 1.0 / f2s;
      f2s = scale + g2;
      g2 = r->re / f2s;
      f2s = r->im / f2s;
      sn->re = g2 * gs_re - f2s * -gs_im;
      sn->im = g2 * -gs_im + f2s * gs_re;
      if (rescaledir > 0) {
        for (rescaledir = 1; rescaledir <= count; rescaledir++) {
          r->re *= 7.4428285367870146E+137;
          r->im *= 7.4428285367870146E+137;
        }
      } else {
        if (rescaledir < 0) {
          for (rescaledir = 1; rescaledir <= count; rescaledir++) {
            r->re *= 1.3435752215134178E-138;
            r->im *= 1.3435752215134178E-138;
          }
        }
      }
    }
  }
}

static void eml_matlab_ztgevc(const creal_T A[9], creal_T V[9])
{
  creal_T work1[3];
  creal_T work2[3];
  real_T rworka[3];
  int32_T i;
  real_T anorm;
  int32_T j;
  real_T y;
  real_T xmx;
  real_T ascale;
  int32_T je;
  real_T temp;
  real_T salpha_re;
  real_T salpha_im;
  real_T acoeff;
  boolean_T b0;
  boolean_T b1;
  real_T scale;
  real_T acoefa;
  int32_T jr;
  real_T dmin;
  creal_T d;
  creal_T b_work1;
  for (i = 0; i < 3; i++) {
    work1[i].re = 0.0;
    work1[i].im = 0.0;
    work2[i].re = 0.0;
    work2[i].im = 0.0;
    rworka[i] = 0.0;
  }

  anorm = fabs(A[0].re) + fabs(A[0].im);
  for (j = 0; j < 2; j++) {
    for (i = 0; i <= j; i++) {
      rworka[1 + j] += fabs(A[i + 3 * (1 + j)].re) + fabs(A[i + 3 * (1 + j)].im);
    }

    y = rworka[1 + j] + (fabs(A[(j + 3 * (1 + j)) + 1].re) + fabs(A[(j + 3 * (1
      + j)) + 1].im));
    if (y > anorm) {
      anorm = y;
    }
  }

  xmx = anorm;
  if (2.2250738585072014E-308 > anorm) {
    xmx = 2.2250738585072014E-308;
  }

  ascale = 1.0 / xmx;
  for (je = 0; je < 3; je++) {
    y = (fabs(A[(3 * (2 - je) - je) + 2].re) + fabs(A[(3 * (2 - je) - je) + 2].
          im)) * ascale;
    if (1.0 > y) {
      y = 1.0;
    }

    temp = 1.0 / y;
    salpha_re = ascale * (temp * A[(3 * (2 - je) - je) + 2].re);
    salpha_im = ascale * (temp * A[(3 * (2 - je) - je) + 2].im);
    acoeff = temp * ascale;
    if ((fabs(temp) >= 2.2250738585072014E-308) && (fabs(acoeff) <
         3.0062525400134592E-292)) {
      b0 = TRUE;
    } else {
      b0 = FALSE;
    }

    if ((fabs(salpha_re) + fabs(salpha_im) >= 2.2250738585072014E-308) && (fabs
         (salpha_re) + fabs(salpha_im) < 3.0062525400134592E-292)) {
      b1 = TRUE;
    } else {
      b1 = FALSE;
    }

    scale = 1.0;
    if (b0) {
      xmx = anorm;
      if (3.3264005158911995E+291 < anorm) {
        xmx = 3.3264005158911995E+291;
      }

      scale = 3.0062525400134592E-292 / fabs(temp) * xmx;
    }

    if (b1) {
      xmx = 3.0062525400134592E-292 / (fabs(salpha_re) + fabs(salpha_im));
      if (xmx > scale) {
        scale = xmx;
      }
    }

    if (b0 || b1) {
      y = fabs(acoeff);
      xmx = fabs(salpha_re) + fabs(salpha_im);
      if (1.0 > y) {
        y = 1.0;
      }

      if (xmx > y) {
        y = xmx;
      }

      y = 1.0 / (2.2250738585072014E-308 * y);
      if (y < scale) {
        scale = y;
      }

      if (b0) {
        acoeff = ascale * (scale * temp);
      } else {
        acoeff *= scale;
      }

      if (b1) {
        salpha_re *= scale;
        salpha_im *= scale;
      } else {
        salpha_re *= scale;
        salpha_im *= scale;
      }
    }

    acoefa = fabs(acoeff);
    for (jr = 0; jr < 3; jr++) {
      work1[jr].re = 0.0;
      work1[jr].im = 0.0;
    }

    work1[2 - je].re = 1.0;
    work1[2 - je].im = 0.0;
    dmin = 2.2204460492503131E-16 * acoefa * anorm;
    y = 2.2204460492503131E-16 * (fabs(salpha_re) + fabs(salpha_im));
    if (y > dmin) {
      dmin = y;
    }

    if (2.2250738585072014E-308 > dmin) {
      dmin = 2.2250738585072014E-308;
    }

    for (jr = 0; jr <= 1 - je; jr++) {
      work1[jr].re = acoeff * A[jr + 3 * (2 - je)].re;
      work1[jr].im = acoeff * A[jr + 3 * (2 - je)].im;
    }

    work1[2 - je].re = 1.0;
    work1[2 - je].im = 0.0;
    for (j = 0; j <= 1 - je; j++) {
      i = (-je - j) + 1;
      d.re = acoeff * A[i + 3 * i].re - salpha_re;
      d.im = acoeff * A[i + 3 * i].im - salpha_im;
      if (fabs(d.re) + fabs(d.im) <= dmin) {
        d.re = dmin;
        d.im = 0.0;
      }

      if ((fabs(d.re) + fabs(d.im) < 1.0) && (fabs(work1[i].re) + fabs(work1[i].
            im) >= 1.4980776123852632E+307 * (fabs(d.re) + fabs(d.im)))) {
        temp = 1.0 / (fabs(work1[i].re) + fabs(work1[i].im));
        for (jr = 0; jr <= 2 - je; jr++) {
          work1[jr].re *= temp;
          work1[jr].im *= temp;
        }
      }

      b_work1.re = -work1[i].re;
      b_work1.im = -work1[i].im;
      work1[i] = c_eml_div(b_work1, d);
      if (i + 1 > 1) {
        if (fabs(work1[1].re) + fabs(work1[1].im) > 1.0) {
          temp = 1.0 / (fabs(work1[1].re) + fabs(work1[1].im));
          if (acoefa * rworka[1] >= 1.4980776123852632E+307 * temp) {
            for (jr = 0; jr <= 2 - je; jr++) {
              work1[jr].re *= temp;
              work1[jr].im *= temp;
            }
          }
        }

        xmx = acoeff * work1[1].re;
        scale = acoeff * work1[1].im;
        work1[0].re += xmx * A[3].re - scale * A[3].im;
        work1[0].im += xmx * A[3].im + scale * A[3].re;
      }
    }

    for (jr = 0; jr < 3; jr++) {
      work2[jr].re = 0.0;
      work2[jr].im = 0.0;
    }

    for (i = 0; i <= 2 - je; i++) {
      for (jr = 0; jr < 3; jr++) {
        xmx = V[jr + 3 * i].re * work1[i].re - V[jr + 3 * i].im * work1[i].im;
        scale = V[jr + 3 * i].re * work1[i].im + V[jr + 3 * i].im * work1[i].re;
        work2[jr].re += xmx;
        work2[jr].im += scale;
      }
    }

    xmx = fabs(work2[0].re) + fabs(work2[0].im);
    for (jr = 0; jr < 2; jr++) {
      y = fabs(work2[1 + jr].re) + fabs(work2[1 + jr].im);
      if (y > xmx) {
        xmx = y;
      }
    }

    if (xmx > 2.2250738585072014E-308) {
      temp = 1.0 / xmx;
      for (jr = 0; jr < 3; jr++) {
        V[jr + 3 * (2 - je)].re = temp * work2[jr].re;
        V[jr + 3 * (2 - je)].im = temp * work2[jr].im;
      }
    } else {
      for (jr = 0; jr < 3; jr++) {
        V[jr + 3 * (2 - je)].re = 0.0;
        V[jr + 3 * (2 - je)].im = 0.0;
      }
    }
  }
}

void b_eig(const real_T A[4], creal_T V[4], creal_T D[4])
{
  creal_T b_A[4];
  int32_T j;
  creal_T beta1[2];
  creal_T alpha1[2];
  real_T colnorm;
  int32_T coltop;
  real_T scale;
  real_T absxk;
  real_T t;
  real_T alpha1_im;
  for (j = 0; j < 4; j++) {
    b_A[j].re = A[j];
    b_A[j].im = 0.0;
  }

  b_eml_matlab_zggev(b_A, &colnorm, alpha1, beta1, V);
  for (coltop = 0; coltop < 4; coltop += 2) {
    colnorm = 0.0;
    scale = 2.2250738585072014E-308;
    for (j = coltop; j + 1 <= coltop + 2; j++) {
      absxk = fabs(V[j].re);
      if (absxk > scale) {
        t = scale / absxk;
        colnorm = 1.0 + colnorm * t * t;
        scale = absxk;
      } else {
        t = absxk / scale;
        colnorm += t * t;
      }

      absxk = fabs(V[j].im);
      if (absxk > scale) {
        t = scale / absxk;
        colnorm = 1.0 + colnorm * t * t;
        scale = absxk;
      } else {
        t = absxk / scale;
        colnorm += t * t;
      }
    }

    colnorm = scale * sqrt(colnorm);
    for (j = coltop; j + 1 <= coltop + 2; j++) {
      V[j] = b_eml_div(V[j], colnorm);
    }
  }

  for (j = 0; j < 2; j++) {
    t = alpha1[j].re;
    alpha1_im = alpha1[j].im;
    if (beta1[j].im == 0.0) {
      if (alpha1[j].im == 0.0) {
        alpha1[j].re /= beta1[j].re;
        alpha1[j].im = 0.0;
      } else if (alpha1[j].re == 0.0) {
        alpha1[j].re = 0.0;
        alpha1[j].im = alpha1_im / beta1[j].re;
      } else {
        alpha1[j].re /= beta1[j].re;
        alpha1[j].im = alpha1_im / beta1[j].re;
      }
    } else if (beta1[j].re == 0.0) {
      if (alpha1[j].re == 0.0) {
        alpha1[j].re = alpha1[j].im / beta1[j].im;
        alpha1[j].im = 0.0;
      } else if (alpha1[j].im == 0.0) {
        alpha1[j].re = 0.0;
        alpha1[j].im = -(t / beta1[j].im);
      } else {
        alpha1[j].re = alpha1[j].im / beta1[j].im;
        alpha1[j].im = -(t / beta1[j].im);
      }
    } else {
      absxk = fabs(beta1[j].re);
      colnorm = fabs(beta1[j].im);
      if (absxk > colnorm) {
        colnorm = beta1[j].im / beta1[j].re;
        scale = beta1[j].re + colnorm * beta1[j].im;
        alpha1[j].re = (alpha1[j].re + colnorm * alpha1[j].im) / scale;
        alpha1[j].im = (alpha1_im - colnorm * t) / scale;
      } else if (colnorm == absxk) {
        if (beta1[j].re > 0.0) {
          colnorm = 0.5;
        } else {
          colnorm = -0.5;
        }

        if (beta1[j].im > 0.0) {
          scale = 0.5;
        } else {
          scale = -0.5;
        }

        alpha1[j].re = (alpha1[j].re * colnorm + alpha1[j].im * scale) / absxk;
        alpha1[j].im = (alpha1_im * colnorm - t * scale) / absxk;
      } else {
        colnorm = beta1[j].re / beta1[j].im;
        scale = beta1[j].im + colnorm * beta1[j].re;
        alpha1[j].re = (colnorm * alpha1[j].re + alpha1[j].im) / scale;
        alpha1[j].im = (colnorm * alpha1_im - t) / scale;
      }
    }
  }

  for (j = 0; j < 4; j++) {
    D[j].re = 0.0;
    D[j].im = 0.0;
  }

  for (j = 0; j < 2; j++) {
    D[j + (j << 1)] = alpha1[j];
  }
}

void eig(const real_T A[9], creal_T V[9], creal_T D[9])
{
  creal_T b_A[9];
  int32_T j;
  creal_T beta1[3];
  creal_T alpha1[3];
  real_T colnorm;
  int32_T coltop;
  real_T scale;
  real_T absxk;
  real_T t;
  real_T alpha1_im;
  for (j = 0; j < 9; j++) {
    b_A[j].re = A[j];
    b_A[j].im = 0.0;
  }

  eml_matlab_zggev(b_A, &colnorm, alpha1, beta1, V);
  for (coltop = 0; coltop < 8; coltop += 3) {
    colnorm = 0.0;
    scale = 2.2250738585072014E-308;
    for (j = coltop; j + 1 <= coltop + 3; j++) {
      absxk = fabs(V[j].re);
      if (absxk > scale) {
        t = scale / absxk;
        colnorm = 1.0 + colnorm * t * t;
        scale = absxk;
      } else {
        t = absxk / scale;
        colnorm += t * t;
      }

      absxk = fabs(V[j].im);
      if (absxk > scale) {
        t = scale / absxk;
        colnorm = 1.0 + colnorm * t * t;
        scale = absxk;
      } else {
        t = absxk / scale;
        colnorm += t * t;
      }
    }

    colnorm = scale * sqrt(colnorm);
    for (j = coltop; j + 1 <= coltop + 3; j++) {
      V[j] = b_eml_div(V[j], colnorm);
    }
  }

  for (j = 0; j < 3; j++) {
    t = alpha1[j].re;
    alpha1_im = alpha1[j].im;
    if (beta1[j].im == 0.0) {
      if (alpha1[j].im == 0.0) {
        alpha1[j].re /= beta1[j].re;
        alpha1[j].im = 0.0;
      } else if (alpha1[j].re == 0.0) {
        alpha1[j].re = 0.0;
        alpha1[j].im = alpha1_im / beta1[j].re;
      } else {
        alpha1[j].re /= beta1[j].re;
        alpha1[j].im = alpha1_im / beta1[j].re;
      }
    } else if (beta1[j].re == 0.0) {
      if (alpha1[j].re == 0.0) {
        alpha1[j].re = alpha1[j].im / beta1[j].im;
        alpha1[j].im = 0.0;
      } else if (alpha1[j].im == 0.0) {
        alpha1[j].re = 0.0;
        alpha1[j].im = -(t / beta1[j].im);
      } else {
        alpha1[j].re = alpha1[j].im / beta1[j].im;
        alpha1[j].im = -(t / beta1[j].im);
      }
    } else {
      absxk = fabs(beta1[j].re);
      colnorm = fabs(beta1[j].im);
      if (absxk > colnorm) {
        colnorm = beta1[j].im / beta1[j].re;
        scale = beta1[j].re + colnorm * beta1[j].im;
        alpha1[j].re = (alpha1[j].re + colnorm * alpha1[j].im) / scale;
        alpha1[j].im = (alpha1_im - colnorm * t) / scale;
      } else if (colnorm == absxk) {
        if (beta1[j].re > 0.0) {
          colnorm = 0.5;
        } else {
          colnorm = -0.5;
        }

        if (beta1[j].im > 0.0) {
          scale = 0.5;
        } else {
          scale = -0.5;
        }

        alpha1[j].re = (alpha1[j].re * colnorm + alpha1[j].im * scale) / absxk;
        alpha1[j].im = (alpha1_im * colnorm - t * scale) / absxk;
      } else {
        colnorm = beta1[j].re / beta1[j].im;
        scale = beta1[j].im + colnorm * beta1[j].re;
        alpha1[j].re = (colnorm * alpha1[j].re + alpha1[j].im) / scale;
        alpha1[j].im = (colnorm * alpha1_im - t) / scale;
      }
    }
  }

  for (j = 0; j < 9; j++) {
    D[j].re = 0.0;
    D[j].im = 0.0;
  }

  for (j = 0; j < 3; j++) {
    D[j + 3 * j] = alpha1[j];
  }
}

/* End of code generation (eig.cpp) */
