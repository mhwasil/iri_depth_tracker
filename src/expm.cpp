/*
 * expm.cpp
 *
 * Code generation for function 'expm'
 *
 * C source code generated on: Tue Apr  1 12:30:12 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Optimal_affine_tracking_3d16_fast_realtime.h"
#include "expm.h"
#include "mldivide.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
static void PadeApproximantOfDegree(const real_T A[9], real_T m, real_T F[9]);
static void b_PadeApproximantOfDegree(const real_T A[4], real_T m, real_T F[4]);
static real_T rt_powd_snf(real_T u0, real_T u1);

/* Function Definitions */
static void PadeApproximantOfDegree(const real_T A[9], real_T m, real_T F[9])
{
  int32_T i8;
  int32_T k;
  real_T A2[9];
  int32_T i9;
  real_T U[9];
  real_T A4[9];
  real_T V[9];
  real_T d;
  real_T A3[9];
  real_T b_A4[9];
  for (i8 = 0; i8 < 3; i8++) {
    for (k = 0; k < 3; k++) {
      A2[i8 + 3 * k] = 0.0;
      for (i9 = 0; i9 < 3; i9++) {
        A2[i8 + 3 * k] += A[i8 + 3 * i9] * A[i9 + 3 * k];
      }
    }
  }

  if (m == 3.0) {
    memcpy(&U[0], &A2[0], 9U * sizeof(real_T));
    for (k = 0; k < 3; k++) {
      U[k + 3 * k] += 60.0;
    }

    for (i8 = 0; i8 < 3; i8++) {
      for (k = 0; k < 3; k++) {
        A4[i8 + 3 * k] = 0.0;
        for (i9 = 0; i9 < 3; i9++) {
          A4[i8 + 3 * k] += A[i8 + 3 * i9] * U[i9 + 3 * k];
        }
      }
    }

    for (i8 = 0; i8 < 3; i8++) {
      for (k = 0; k < 3; k++) {
        U[k + 3 * i8] = A4[k + 3 * i8];
      }
    }

    for (i8 = 0; i8 < 9; i8++) {
      V[i8] = 12.0 * A2[i8];
    }

    d = 120.0;
  } else {
    for (i8 = 0; i8 < 3; i8++) {
      for (k = 0; k < 3; k++) {
        A3[i8 + 3 * k] = 0.0;
        for (i9 = 0; i9 < 3; i9++) {
          A3[i8 + 3 * k] += A2[i8 + 3 * i9] * A2[i9 + 3 * k];
        }
      }
    }

    if (m == 5.0) {
      for (i8 = 0; i8 < 9; i8++) {
        U[i8] = A3[i8] + 420.0 * A2[i8];
      }

      for (k = 0; k < 3; k++) {
        U[k + 3 * k] += 15120.0;
      }

      for (i8 = 0; i8 < 3; i8++) {
        for (k = 0; k < 3; k++) {
          A4[i8 + 3 * k] = 0.0;
          for (i9 = 0; i9 < 3; i9++) {
            A4[i8 + 3 * k] += A[i8 + 3 * i9] * U[i9 + 3 * k];
          }
        }
      }

      for (i8 = 0; i8 < 3; i8++) {
        for (k = 0; k < 3; k++) {
          U[k + 3 * i8] = A4[k + 3 * i8];
        }
      }

      for (i8 = 0; i8 < 9; i8++) {
        V[i8] = 30.0 * A3[i8] + 3360.0 * A2[i8];
      }

      d = 30240.0;
    } else {
      for (i8 = 0; i8 < 3; i8++) {
        for (k = 0; k < 3; k++) {
          b_A4[i8 + 3 * k] = 0.0;
          for (i9 = 0; i9 < 3; i9++) {
            b_A4[i8 + 3 * k] += A3[i8 + 3 * i9] * A2[i9 + 3 * k];
          }
        }
      }

      if (m == 7.0) {
        for (i8 = 0; i8 < 9; i8++) {
          U[i8] = (b_A4[i8] + 1512.0 * A3[i8]) + 277200.0 * A2[i8];
        }

        for (k = 0; k < 3; k++) {
          U[k + 3 * k] += 8.64864E+6;
        }

        for (i8 = 0; i8 < 3; i8++) {
          for (k = 0; k < 3; k++) {
            A4[i8 + 3 * k] = 0.0;
            for (i9 = 0; i9 < 3; i9++) {
              A4[i8 + 3 * k] += A[i8 + 3 * i9] * U[i9 + 3 * k];
            }
          }
        }

        for (i8 = 0; i8 < 3; i8++) {
          for (k = 0; k < 3; k++) {
            U[k + 3 * i8] = A4[k + 3 * i8];
          }
        }

        for (i8 = 0; i8 < 9; i8++) {
          V[i8] = (56.0 * b_A4[i8] + 25200.0 * A3[i8]) + 1.99584E+6 * A2[i8];
        }

        d = 1.729728E+7;
      } else if (m == 9.0) {
        for (i8 = 0; i8 < 3; i8++) {
          for (k = 0; k < 3; k++) {
            V[i8 + 3 * k] = 0.0;
            for (i9 = 0; i9 < 3; i9++) {
              V[i8 + 3 * k] += b_A4[i8 + 3 * i9] * A2[i9 + 3 * k];
            }
          }
        }

        for (i8 = 0; i8 < 9; i8++) {
          U[i8] = ((V[i8] + 3960.0 * b_A4[i8]) + 2.16216E+6 * A3[i8]) +
            3.027024E+8 * A2[i8];
        }

        for (k = 0; k < 3; k++) {
          U[k + 3 * k] += 8.8216128E+9;
        }

        for (i8 = 0; i8 < 3; i8++) {
          for (k = 0; k < 3; k++) {
            A4[i8 + 3 * k] = 0.0;
            for (i9 = 0; i9 < 3; i9++) {
              A4[i8 + 3 * k] += A[i8 + 3 * i9] * U[i9 + 3 * k];
            }
          }
        }

        for (i8 = 0; i8 < 3; i8++) {
          for (k = 0; k < 3; k++) {
            U[k + 3 * i8] = A4[k + 3 * i8];
          }
        }

        for (i8 = 0; i8 < 9; i8++) {
          V[i8] = ((90.0 * V[i8] + 110880.0 * b_A4[i8]) + 3.027024E+7 * A3[i8])
            + 2.0756736E+9 * A2[i8];
        }

        d = 1.76432256E+10;
      } else {
        for (i8 = 0; i8 < 9; i8++) {
          U[i8] = (3.352212864E+10 * b_A4[i8] + 1.05594705216E+13 * A3[i8]) +
            1.1873537964288E+15 * A2[i8];
        }

        for (k = 0; k < 3; k++) {
          U[k + 3 * k] += 3.238237626624E+16;
          for (i8 = 0; i8 < 3; i8++) {
            A4[i8 + 3 * k] = (b_A4[i8 + 3 * k] + 16380.0 * A3[i8 + 3 * k]) +
              4.08408E+7 * A2[i8 + 3 * k];
          }
        }

        for (i8 = 0; i8 < 3; i8++) {
          for (k = 0; k < 3; k++) {
            d = 0.0;
            for (i9 = 0; i9 < 3; i9++) {
              d += b_A4[i8 + 3 * i9] * A4[i9 + 3 * k];
            }

            V[i8 + 3 * k] = d + U[i8 + 3 * k];
          }
        }

        for (i8 = 0; i8 < 3; i8++) {
          for (k = 0; k < 3; k++) {
            U[i8 + 3 * k] = 0.0;
            for (i9 = 0; i9 < 3; i9++) {
              U[i8 + 3 * k] += A[i8 + 3 * i9] * V[i9 + 3 * k];
            }

            A4[k + 3 * i8] = (182.0 * b_A4[k + 3 * i8] + 960960.0 * A3[k + 3 *
                              i8]) + 1.32324192E+9 * A2[k + 3 * i8];
          }
        }

        for (i8 = 0; i8 < 3; i8++) {
          for (k = 0; k < 3; k++) {
            d = 0.0;
            for (i9 = 0; i9 < 3; i9++) {
              d += b_A4[i8 + 3 * i9] * A4[i9 + 3 * k];
            }

            V[i8 + 3 * k] = ((d + 6.704425728E+11 * b_A4[i8 + 3 * k]) +
                             1.29060195264E+14 * A3[i8 + 3 * k]) +
              7.7717703038976E+15 * A2[i8 + 3 * k];
          }
        }

        d = 6.476475253248E+16;
      }
    }
  }

  for (k = 0; k < 3; k++) {
    V[k + 3 * k] += d;
  }

  for (k = 0; k < 9; k++) {
    d = V[k] - U[k];
    V[k] += U[k];
    U[k] = d;
  }

  mldivide(U, V, F);
}

static void b_PadeApproximantOfDegree(const real_T A[4], real_T m, real_T F[4])
{
  int32_T i19;
  int32_T k;
  real_T A2[4];
  int32_T i20;
  real_T U[4];
  real_T A4[4];
  real_T V[4];
  real_T d;
  real_T A3[4];
  real_T b_A4[4];
  for (i19 = 0; i19 < 2; i19++) {
    for (k = 0; k < 2; k++) {
      A2[i19 + (k << 1)] = 0.0;
      for (i20 = 0; i20 < 2; i20++) {
        A2[i19 + (k << 1)] += A[i19 + (i20 << 1)] * A[i20 + (k << 1)];
      }
    }
  }

  if (m == 3.0) {
    for (i19 = 0; i19 < 4; i19++) {
      U[i19] = A2[i19];
    }

    for (k = 0; k < 2; k++) {
      U[k + (k << 1)] += 60.0;
    }

    for (i19 = 0; i19 < 2; i19++) {
      for (k = 0; k < 2; k++) {
        A4[i19 + (k << 1)] = 0.0;
        for (i20 = 0; i20 < 2; i20++) {
          A4[i19 + (k << 1)] += A[i19 + (i20 << 1)] * U[i20 + (k << 1)];
        }
      }
    }

    for (i19 = 0; i19 < 2; i19++) {
      for (k = 0; k < 2; k++) {
        U[k + (i19 << 1)] = A4[k + (i19 << 1)];
      }
    }

    for (i19 = 0; i19 < 4; i19++) {
      V[i19] = 12.0 * A2[i19];
    }

    d = 120.0;
  } else {
    for (i19 = 0; i19 < 2; i19++) {
      for (k = 0; k < 2; k++) {
        A3[i19 + (k << 1)] = 0.0;
        for (i20 = 0; i20 < 2; i20++) {
          A3[i19 + (k << 1)] += A2[i19 + (i20 << 1)] * A2[i20 + (k << 1)];
        }
      }
    }

    if (m == 5.0) {
      for (i19 = 0; i19 < 4; i19++) {
        U[i19] = A3[i19] + 420.0 * A2[i19];
      }

      for (k = 0; k < 2; k++) {
        U[k + (k << 1)] += 15120.0;
      }

      for (i19 = 0; i19 < 2; i19++) {
        for (k = 0; k < 2; k++) {
          A4[i19 + (k << 1)] = 0.0;
          for (i20 = 0; i20 < 2; i20++) {
            A4[i19 + (k << 1)] += A[i19 + (i20 << 1)] * U[i20 + (k << 1)];
          }
        }
      }

      for (i19 = 0; i19 < 2; i19++) {
        for (k = 0; k < 2; k++) {
          U[k + (i19 << 1)] = A4[k + (i19 << 1)];
        }
      }

      for (i19 = 0; i19 < 4; i19++) {
        V[i19] = 30.0 * A3[i19] + 3360.0 * A2[i19];
      }

      d = 30240.0;
    } else {
      for (i19 = 0; i19 < 2; i19++) {
        for (k = 0; k < 2; k++) {
          b_A4[i19 + (k << 1)] = 0.0;
          for (i20 = 0; i20 < 2; i20++) {
            b_A4[i19 + (k << 1)] += A3[i19 + (i20 << 1)] * A2[i20 + (k << 1)];
          }
        }
      }

      if (m == 7.0) {
        for (i19 = 0; i19 < 4; i19++) {
          U[i19] = (b_A4[i19] + 1512.0 * A3[i19]) + 277200.0 * A2[i19];
        }

        for (k = 0; k < 2; k++) {
          U[k + (k << 1)] += 8.64864E+6;
        }

        for (i19 = 0; i19 < 2; i19++) {
          for (k = 0; k < 2; k++) {
            A4[i19 + (k << 1)] = 0.0;
            for (i20 = 0; i20 < 2; i20++) {
              A4[i19 + (k << 1)] += A[i19 + (i20 << 1)] * U[i20 + (k << 1)];
            }
          }
        }

        for (i19 = 0; i19 < 2; i19++) {
          for (k = 0; k < 2; k++) {
            U[k + (i19 << 1)] = A4[k + (i19 << 1)];
          }
        }

        for (i19 = 0; i19 < 4; i19++) {
          V[i19] = (56.0 * b_A4[i19] + 25200.0 * A3[i19]) + 1.99584E+6 * A2[i19];
        }

        d = 1.729728E+7;
      } else if (m == 9.0) {
        for (i19 = 0; i19 < 2; i19++) {
          for (k = 0; k < 2; k++) {
            V[i19 + (k << 1)] = 0.0;
            for (i20 = 0; i20 < 2; i20++) {
              V[i19 + (k << 1)] += b_A4[i19 + (i20 << 1)] * A2[i20 + (k << 1)];
            }
          }
        }

        for (i19 = 0; i19 < 4; i19++) {
          U[i19] = ((V[i19] + 3960.0 * b_A4[i19]) + 2.16216E+6 * A3[i19]) +
            3.027024E+8 * A2[i19];
        }

        for (k = 0; k < 2; k++) {
          U[k + (k << 1)] += 8.8216128E+9;
        }

        for (i19 = 0; i19 < 2; i19++) {
          for (k = 0; k < 2; k++) {
            A4[i19 + (k << 1)] = 0.0;
            for (i20 = 0; i20 < 2; i20++) {
              A4[i19 + (k << 1)] += A[i19 + (i20 << 1)] * U[i20 + (k << 1)];
            }
          }
        }

        for (i19 = 0; i19 < 2; i19++) {
          for (k = 0; k < 2; k++) {
            U[k + (i19 << 1)] = A4[k + (i19 << 1)];
          }
        }

        for (i19 = 0; i19 < 4; i19++) {
          V[i19] = ((90.0 * V[i19] + 110880.0 * b_A4[i19]) + 3.027024E+7 *
                    A3[i19]) + 2.0756736E+9 * A2[i19];
        }

        d = 1.76432256E+10;
      } else {
        for (i19 = 0; i19 < 4; i19++) {
          U[i19] = (3.352212864E+10 * b_A4[i19] + 1.05594705216E+13 * A3[i19]) +
            1.1873537964288E+15 * A2[i19];
        }

        for (k = 0; k < 2; k++) {
          U[k + (k << 1)] += 3.238237626624E+16;
          for (i19 = 0; i19 < 2; i19++) {
            A4[i19 + (k << 1)] = (b_A4[i19 + (k << 1)] + 16380.0 * A3[i19 + (k <<
              1)]) + 4.08408E+7 * A2[i19 + (k << 1)];
          }
        }

        for (i19 = 0; i19 < 2; i19++) {
          for (k = 0; k < 2; k++) {
            d = 0.0;
            for (i20 = 0; i20 < 2; i20++) {
              d += b_A4[i19 + (i20 << 1)] * A4[i20 + (k << 1)];
            }

            V[i19 + (k << 1)] = d + U[i19 + (k << 1)];
          }
        }

        for (i19 = 0; i19 < 2; i19++) {
          for (k = 0; k < 2; k++) {
            U[i19 + (k << 1)] = 0.0;
            for (i20 = 0; i20 < 2; i20++) {
              U[i19 + (k << 1)] += A[i19 + (i20 << 1)] * V[i20 + (k << 1)];
            }

            A4[k + (i19 << 1)] = (182.0 * b_A4[k + (i19 << 1)] + 960960.0 * A3[k
                                  + (i19 << 1)]) + 1.32324192E+9 * A2[k + (i19 <<
              1)];
          }
        }

        for (i19 = 0; i19 < 2; i19++) {
          for (k = 0; k < 2; k++) {
            d = 0.0;
            for (i20 = 0; i20 < 2; i20++) {
              d += b_A4[i19 + (i20 << 1)] * A4[i20 + (k << 1)];
            }

            V[i19 + (k << 1)] = ((d + 6.704425728E+11 * b_A4[i19 + (k << 1)]) +
                                 1.29060195264E+14 * A3[i19 + (k << 1)]) +
              7.7717703038976E+15 * A2[i19 + (k << 1)];
          }
        }

        d = 6.476475253248E+16;
      }
    }
  }

  for (k = 0; k < 2; k++) {
    V[k + (k << 1)] += d;
  }

  for (k = 0; k < 4; k++) {
    d = V[k] - U[k];
    V[k] += U[k];
    U[k] = d;
  }

  b_mldivide(U, V, F);
}

static real_T rt_powd_snf(real_T u0, real_T u1)
{
  real_T y;
  real_T d3;
  real_T d4;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else {
    d3 = fabs(u0);
    d4 = fabs(u1);
    if (rtIsInf(u1)) {
      if (d3 == 1.0) {
        y = rtNaN;
      } else if (d3 > 1.0) {
        if (u1 > 0.0) {
          y = rtInf;
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = rtInf;
      }
    } else if (d4 == 0.0) {
      y = 1.0;
    } else if (d4 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > floor(u1))) {
      y = rtNaN;
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
}

void b_expm(real_T A[4], real_T F[4])
{
  real_T normA;
  int32_T j;
  boolean_T exitg2;
  real_T s;
  int32_T i;
  boolean_T exitg1;
  static const real_T theta[5] = { 0.01495585217958292, 0.253939833006323,
    0.95041789961629319, 2.097847961257068, 5.3719203511481517 };

  static const int8_T iv8[5] = { 3, 5, 7, 9, 13 };

  int32_T eint;
  real_T b_F[4];
  int32_T i18;
  normA = 0.0;
  j = 0;
  exitg2 = FALSE;
  while ((exitg2 == FALSE) && (j < 2)) {
    s = 0.0;
    for (i = 0; i < 2; i++) {
      s += fabs(A[i + (j << 1)]);
    }

    if (rtIsNaN(s)) {
      normA = rtNaN;
      exitg2 = TRUE;
    } else {
      if (s > normA) {
        normA = s;
      }

      j++;
    }
  }

  if (normA <= 5.3719203511481517) {
    i = 0;
    exitg1 = FALSE;
    while ((exitg1 == FALSE) && (i < 5)) {
      if (normA <= theta[i]) {
        b_PadeApproximantOfDegree(A, (real_T)iv8[i], F);
        exitg1 = TRUE;
      } else {
        i++;
      }
    }
  } else {
    normA /= 5.3719203511481517;
    if ((!rtIsInf(normA)) && (!rtIsNaN(normA))) {
      normA = frexp(normA, &eint);
    } else {
      eint = 0;
    }

    s = (real_T)eint;
    if (normA == 0.5) {
      s = (real_T)eint - 1.0;
    }

    normA = rt_powd_snf(2.0, s);
    for (i = 0; i < 4; i++) {
      A[i] /= normA;
    }

    b_PadeApproximantOfDegree(A, 13.0, F);
    for (j = 0; j < (int32_T)s; j++) {
      for (i = 0; i < 2; i++) {
        for (eint = 0; eint < 2; eint++) {
          b_F[i + (eint << 1)] = 0.0;
          for (i18 = 0; i18 < 2; i18++) {
            b_F[i + (eint << 1)] += F[i + (i18 << 1)] * F[i18 + (eint << 1)];
          }
        }
      }

      for (i = 0; i < 2; i++) {
        for (eint = 0; eint < 2; eint++) {
          F[eint + (i << 1)] = b_F[eint + (i << 1)];
        }
      }
    }
  }
}

void expm(real_T A[9], real_T F[9])
{
  real_T normA;
  int32_T j;
  boolean_T exitg2;
  real_T s;
  int32_T i;
  boolean_T exitg1;
  static const real_T theta[5] = { 0.01495585217958292, 0.253939833006323,
    0.95041789961629319, 2.097847961257068, 5.3719203511481517 };

  static const int8_T iv5[5] = { 3, 5, 7, 9, 13 };

  int32_T eint;
  real_T b_F[9];
  int32_T i7;
  normA = 0.0;
  j = 0;
  exitg2 = FALSE;
  while ((exitg2 == FALSE) && (j < 3)) {
    s = 0.0;
    for (i = 0; i < 3; i++) {
      s += fabs(A[i + 3 * j]);
    }

    if (rtIsNaN(s)) {
      normA = rtNaN;
      exitg2 = TRUE;
    } else {
      if (s > normA) {
        normA = s;
      }

      j++;
    }
  }

  if (normA <= 5.3719203511481517) {
    i = 0;
    exitg1 = FALSE;
    while ((exitg1 == FALSE) && (i < 5)) {
      if (normA <= theta[i]) {
        PadeApproximantOfDegree(A, (real_T)iv5[i], F);
        exitg1 = TRUE;
      } else {
        i++;
      }
    }
  } else {
    normA /= 5.3719203511481517;
    if ((!rtIsInf(normA)) && (!rtIsNaN(normA))) {
      normA = frexp(normA, &eint);
    } else {
      eint = 0;
    }

    s = (real_T)eint;
    if (normA == 0.5) {
      s = (real_T)eint - 1.0;
    }

    normA = rt_powd_snf(2.0, s);
    for (i = 0; i < 9; i++) {
      A[i] /= normA;
    }

    PadeApproximantOfDegree(A, 13.0, F);
    for (j = 0; j < (int32_T)s; j++) {
      for (i = 0; i < 3; i++) {
        for (eint = 0; eint < 3; eint++) {
          b_F[i + 3 * eint] = 0.0;
          for (i7 = 0; i7 < 3; i7++) {
            b_F[i + 3 * eint] += F[i + 3 * i7] * F[i7 + 3 * eint];
          }
        }
      }

      for (i = 0; i < 3; i++) {
        for (eint = 0; eint < 3; eint++) {
          F[eint + 3 * i] = b_F[eint + 3 * i];
        }
      }
    }
  }
}

/* End of code generation (expm.cpp) */
