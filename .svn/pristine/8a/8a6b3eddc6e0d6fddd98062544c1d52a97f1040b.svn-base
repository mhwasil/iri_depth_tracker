/*
 * compute_prob1.cpp
 *
 * Code generation for function 'compute_prob1'
 *
 * C source code generated on: Tue Apr  1 12:30:12 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Optimal_affine_tracking_3d16_fast_realtime.h"
#include "compute_prob1.h"
#include "norm.h"
#include "estimateRigidTransform.h"
#include "mean.h"
#include "image_warping.h"
#include "mldivide.h"
#include "expm.h"
#include "randn.h"
#include "chol.h"
#include "sum.h"
#include "mrdivide.h"
#include "diag.h"
#include "log.h"
#include "eig.h"
#include "resampling.h"
#include "exp.h"
#include "repmat.h"
#include "gradient.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
static void b_eml_li_find(const boolean_T x[25], int32_T y_data[25], int32_T
  y_size[1]);
static int32_T compute_nones(const boolean_T x[1296]);

/* Function Definitions */
static void b_eml_li_find(const boolean_T x[25], int32_T y_data[25], int32_T
  y_size[1])
{
  int32_T k;
  int32_T i;
  k = 0;
  for (i = 0; i < 25; i++) {
    if (x[i]) {
      k++;
    }
  }

  y_size[0] = k;
  k = 0;
  for (i = 0; i < 25; i++) {
    if (x[i]) {
      y_data[k] = i + 1;
      k++;
    }
  }
}

static int32_T compute_nones(const boolean_T x[1296])
{
  int32_T k;
  int32_T i;
  k = 0;
  for (i = 0; i < 1296; i++) {
    if (x[i]) {
      k++;
    }
  }

  return k;
}

void compute_prob1(real_T X_par[2250000], real_T X_par_pred[2250000], real_T
                   AR_velocity[2250000], const real_T dw_dp[15552], real_T t,
                   real_T center_x, real_T center_y, const real_T Ixyz[504063],
                   const real_T point_matrix[3888], const real_T mean_img[3888],
                   real_T tracked_images[38880000], real_T Aff_matrix[9], real_T
                   centroid[3])
{
  static real_T warped_img[3888];
  static real_T frame[504063];
  static real_T b_Ixyz[168021];
  static real_T c_Ixyz[168021];
  int32_T i21;
  int32_T i;
  static real_T ext_frame[1008126];
  static real_T ext_grad_x[336042];
  static real_T ext_grad_y[336042];
  real_T rigid_transform_par[400];
  real_T dist1Array[25];
  int32_T par;
  creal_T v[9];
  creal_T V[9];
  creal_T dcv0[3];
  creal_T b_V[9];
  int32_T ix;
  real_T LOGMX[9];
  real_T b_LOGMX[9];
  real_T X_par_temp[9];
  real_T mean_trans[2];
  static const int8_T iv9[3] = { 0, 0, 1 };

  boolean_T zero_ind_z[1296];
  boolean_T b_zero_ind_z;
  boolean_T c_zero_ind_z[1296];
  int32_T tmp_size[1];
  int32_T tmp_data[1296];
  int32_T b_tmp_size[1];
  int32_T b_tmp_data[1296];
  real_T warped_grad_x3[1296];
  int32_T warped_img_size[1];
  real_T d;
  int32_T b_warped_img_size[1];
  int32_T c_warped_img_size[1];
  real_T warped_grad_y3[1296];
  real_T img_grad[2592];
  real_T transform[16];
  static real_T b_warped_img[5184];
  static real_T b_transform[5184];
  static real_T a[7776];
  real_T b_X_par_temp[6];
  real_T c_X_par_temp[6];
  real_T d_X_par_temp[6];
  real_T e_X_par_temp[6];
  real_T dv13[6];
  real_T dv14[6];
  static real_T b_a[7776];
  real_T c_a[2592];
  real_T CH[36];
  real_T sum_jacobian[6];
  real_T sigma_12[6];
  static const real_T sigma_11[36] = { 0.00066666666666666675, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 2.6666666666666667E-5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.002666666666666667, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.6666666666666667E-5,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.2666666666666666, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 3.2666666666666666 };

  real_T y[6];
  real_T b;
  static const int8_T b_b[9] = { 0, 0, 0, 0, 0, 0, 0, 1, 0 };

  static const int8_T c_b[9] = { 0, 0, 0, 0, 0, 0, 1, 0, 0 };

  static const int8_T d_b[9] = { 0, 1, 0, 1, 0, 0, 0, 0, 0 };

  static const int8_T e_b[9] = { 0, 1, 0, -1, 0, 0, 0, 0, 0 };

  static const int8_T f_b[9] = { 1, 0, 0, 0, -1, 0, 0, 0, 0 };

  static const int8_T g_b[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 0 };

  real_T b_sigma_12[6];
  real_T h_b[9];
  real_T b_center_x[2];
  int32_T d_warped_img_size[1];
  int32_T e_warped_img_size[1];
  int32_T f_warped_img_size[1];
  static real_T c_transform[5184];
  static real_T warped_img_tr[3888];
  static real_T b_warped_img_tr[3888];
  boolean_T b_dist1Array[25];
  int32_T c_tmp_data[25];
  real_T outindex[25];
  real_T b_X_par_pred[225];
  real_T b_rigid_transform_par[400];
  int32_T itmp;
  boolean_T exitg1;
  real_T mean_log[4];
  real_T c_X_par_pred[4];
  real_T b_X_par[4];
  real_T c_LOGMX[4];
  creal_T b_v[4];
  creal_T c_V[4];
  creal_T dcv1[2];
  creal_T d_V[4];
  real_T c_X_par[50];
  real_T c_center_x[2];
  int32_T g_warped_img_size[1];
  int32_T h_warped_img_size[1];
  int32_T i_warped_img_size[1];
  static real_T c_rigid_transform_par[5184];
  memset(&warped_img[0], 0, 3888U * sizeof(real_T));
  memset(&frame[0], 0, 504063U * sizeof(real_T));
  for (i21 = 0; i21 < 441; i21++) {
    for (i = 0; i < 381; i++) {
      frame[336042 + (i + 381 * i21)] = Ixyz[336042 + (i + 381 * i21)];
      b_Ixyz[i + 381 * i21] = (Ixyz[i + 381 * i21] + 2.0) * (real_T)(frame
        [336042 + (i + 381 * i21)] != 0.0);
      frame[i + 381 * i21] = b_Ixyz[i + 381 * i21];
      c_Ixyz[i + 381 * i21] = (Ixyz[168021 + (i + 381 * i21)] + 2.0) * (real_T)
        (frame[336042 + (i + 381 * i21)] != 0.0);
      frame[168021 + (i + 381 * i21)] = c_Ixyz[i + 381 * i21];
    }
  }

  repmat(frame, ext_frame);
  gradient(*(real_T (*)[168021])&frame[336042], b_Ixyz, c_Ixyz);
  b_repmat(c_Ixyz, ext_grad_x);
  b_repmat(b_Ixyz, ext_grad_y);
  for (par = 0; par < 25; par++) {
    eig(*(real_T (*)[9])&AR_velocity[9 * par + 225 * ((int32_T)(t - 1.0) - 1)],
        V, v);
    for (i21 = 0; i21 < 3; i21++) {
      dcv0[i21] = v[i21 << 2];
    }

    b_log(dcv0);
    diag(dcv0, v);
    for (i21 = 0; i21 < 3; i21++) {
      for (i = 0; i < 3; i++) {
        b_V[i21 + 3 * i].re = 0.0;
        b_V[i21 + 3 * i].im = 0.0;
        for (ix = 0; ix < 3; ix++) {
          b_V[i21 + 3 * i].re += V[i21 + 3 * ix].re * v[ix + 3 * i].re - V[i21 +
            3 * ix].im * v[ix + 3 * i].im;
          b_V[i21 + 3 * i].im += V[i21 + 3 * ix].re * v[ix + 3 * i].im + V[i21 +
            3 * ix].im * v[ix + 3 * i].re;
        }
      }
    }

    b_mrdivide(b_V, V, v);
    for (i21 = 0; i21 < 9; i21++) {
      LOGMX[i21] = v[i21].re * 0.5;
    }

    expm(LOGMX, b_LOGMX);
    for (i21 = 0; i21 < 3; i21++) {
      for (i = 0; i < 3; i++) {
        X_par_temp[i21 + 3 * i] = 0.0;
        for (ix = 0; ix < 3; ix++) {
          X_par_temp[i21 + 3 * i] += X_par[((i21 + 3 * ix) + 9 * par) + 225 *
            ((int32_T)(t - 1.0) - 1)] * b_LOGMX[ix + 3 * i];
        }
      }
    }

    mean_trans[0] = center_x;
    mean_trans[1] = center_y;
    for (i21 = 0; i21 < 2; i21++) {
      for (i = 0; i < 2; i++) {
        Aff_matrix[i + 3 * i21] = X_par_temp[i + 3 * i21];
      }
    }

    for (i21 = 0; i21 < 2; i21++) {
      Aff_matrix[6 + i21] = X_par_temp[6 + i21] + mean_trans[i21];
    }

    for (i21 = 0; i21 < 3; i21++) {
      Aff_matrix[2 + 3 * i21] = (real_T)iv9[i21];
    }

    b_image_warping(*(real_T (*)[336042])&ext_frame[0], Aff_matrix, point_matrix,
                    *(real_T (*)[1296])&warped_img[0]);
    b_image_warping(*(real_T (*)[336042])&ext_frame[336042], Aff_matrix,
                    point_matrix, *(real_T (*)[1296])&warped_img[1296]);
    b_image_warping(*(real_T (*)[336042])&ext_frame[672084], Aff_matrix,
                    point_matrix, *(real_T (*)[1296])&warped_img[2592]);

    /*             %% replace zeros with average values */
    for (i21 = 0; i21 < 1296; i21++) {
      b_zero_ind_z = (warped_img[2592 + i21] == 0.0);
      zero_ind_z[i21] = !b_zero_ind_z;
      c_zero_ind_z[i21] = b_zero_ind_z;
    }

    eml_li_find(c_zero_ind_z, tmp_data, tmp_size);
    eml_li_find(zero_ind_z, b_tmp_data, b_tmp_size);
    warped_img_size[0] = b_tmp_size[0];
    i = b_tmp_size[0];
    for (i21 = 0; i21 < i; i21++) {
      warped_grad_x3[i21] = warped_img[b_tmp_data[i21] + 2591];
    }

    d = mean(warped_grad_x3, warped_img_size);
    i = tmp_size[0];
    for (i21 = 0; i21 < i; i21++) {
      warped_img[tmp_data[i21] + 2591] = d;
    }

    eml_li_find(c_zero_ind_z, tmp_data, tmp_size);
    for (i = 0; i < 1296; i++) {
      zero_ind_z[i] = !c_zero_ind_z[i];
    }

    eml_li_find(zero_ind_z, b_tmp_data, b_tmp_size);
    b_warped_img_size[0] = b_tmp_size[0];
    i = b_tmp_size[0];
    for (i21 = 0; i21 < i; i21++) {
      warped_grad_x3[i21] = warped_img[b_tmp_data[i21] - 1];
    }

    d = mean(warped_grad_x3, b_warped_img_size);
    i = tmp_size[0];
    for (i21 = 0; i21 < i; i21++) {
      warped_img[tmp_data[i21] - 1] = d;
    }

    eml_li_find(c_zero_ind_z, tmp_data, tmp_size);
    for (i = 0; i < 1296; i++) {
      zero_ind_z[i] = !c_zero_ind_z[i];
    }

    eml_li_find(zero_ind_z, b_tmp_data, b_tmp_size);
    c_warped_img_size[0] = b_tmp_size[0];
    i = b_tmp_size[0];
    for (i21 = 0; i21 < i; i21++) {
      warped_grad_x3[i21] = warped_img[b_tmp_data[i21] + 1295];
    }

    d = mean(warped_grad_x3, c_warped_img_size);
    i = tmp_size[0];
    for (i21 = 0; i21 < i; i21++) {
      warped_img[tmp_data[i21] + 1295] = d;
    }

    /*             %% */
    b_image_warping(ext_grad_x, Aff_matrix, point_matrix, warped_grad_x3);
    b_image_warping(ext_grad_y, Aff_matrix, point_matrix, warped_grad_y3);
    for (i21 = 0; i21 < 1296; i21++) {
      img_grad[i21] = warped_grad_x3[i21];
      img_grad[1296 + i21] = warped_grad_y3[i21];
    }

    /*             %% compute rigid transform */
    estimateRigidTransform(mean_img, warped_img, transform);
    for (i21 = 0; i21 < 1296; i21++) {
      for (i = 0; i < 3; i++) {
        b_warped_img[i + (i21 << 2)] = warped_img[i21 + 1296 * i];
      }

      b_warped_img[3 + (i21 << 2)] = 1.0;
      for (i = 0; i < 4; i++) {
        b_transform[i21 + 1296 * i] = 0.0;
        for (ix = 0; ix < 4; ix++) {
          b_transform[i21 + 1296 * i] += transform[i + (ix << 2)] *
            b_warped_img[ix + (i21 << 2)];
        }
      }
    }

    for (i21 = 0; i21 < 1296; i21++) {
      warped_grad_x3[i21] = b_transform[2592 + i21] - mean_img[2592 + i21];
    }

    b_X_par_temp[0] = X_par_temp[0];
    b_X_par_temp[1] = X_par_temp[1];
    b_X_par_temp[2] = X_par_temp[3];
    b_X_par_temp[3] = X_par_temp[4];
    b_X_par_temp[4] = 0.0;
    b_X_par_temp[5] = 0.0;
    c_X_par_temp[0] = X_par_temp[0];
    c_X_par_temp[1] = X_par_temp[1];
    c_X_par_temp[2] = -X_par_temp[3];
    c_X_par_temp[3] = -X_par_temp[4];
    c_X_par_temp[4] = 0.0;
    c_X_par_temp[5] = 0.0;
    d_X_par_temp[0] = X_par_temp[3];
    d_X_par_temp[1] = X_par_temp[4];
    d_X_par_temp[2] = -X_par_temp[0];
    d_X_par_temp[3] = -X_par_temp[1];
    d_X_par_temp[4] = 0.0;
    d_X_par_temp[5] = 0.0;
    e_X_par_temp[0] = X_par_temp[3];
    e_X_par_temp[1] = X_par_temp[4];
    e_X_par_temp[2] = X_par_temp[0];
    e_X_par_temp[3] = X_par_temp[1];
    e_X_par_temp[4] = 0.0;
    e_X_par_temp[5] = 0.0;
    dv13[0] = 0.0;
    dv13[1] = 0.0;
    dv13[2] = 0.0;
    dv13[3] = 0.0;
    dv13[4] = X_par_temp[0];
    dv13[5] = X_par_temp[1];
    dv14[0] = 0.0;
    dv14[1] = 0.0;
    dv14[2] = 0.0;
    dv14[3] = 0.0;
    dv14[4] = X_par_temp[3];
    dv14[5] = X_par_temp[4];
    for (i = 0; i < 1296; i++) {
      d = 2.0 * warped_grad_x3[i];
      c_a[i] = d * img_grad[i];
      c_a[1296 + i] = d * img_grad[1296 + i];
      b_a[i] = c_a[i] * dw_dp[12 * i] + c_a[1296 + i] * dw_dp[1 + 12 * i];
      b_a[1296 + i] = c_a[i] * dw_dp[2 + 12 * i] + c_a[1296 + i] * dw_dp[3 + 12 *
        i];
      b_a[2592 + i] = c_a[i] * dw_dp[4 + 12 * i] + c_a[1296 + i] * dw_dp[5 + 12 *
        i];
      b_a[3888 + i] = c_a[i] * dw_dp[6 + 12 * i] + c_a[1296 + i] * dw_dp[7 + 12 *
        i];
      b_a[5184 + i] = c_a[i] * dw_dp[8 + 12 * i] + c_a[1296 + i] * dw_dp[9 + 12 *
        i];
      b_a[6480 + i] = c_a[i] * dw_dp[10 + 12 * i] + c_a[1296 + i] * dw_dp[11 +
        12 * i];
    }

    for (i21 = 0; i21 < 6; i21++) {
      CH[i21] = b_X_par_temp[i21];
      CH[6 + i21] = c_X_par_temp[i21];
      CH[12 + i21] = d_X_par_temp[i21];
      CH[18 + i21] = e_X_par_temp[i21];
      CH[24 + i21] = dv13[i21];
      CH[30 + i21] = dv14[i21];
    }

    for (i21 = 0; i21 < 1296; i21++) {
      for (i = 0; i < 6; i++) {
        a[i21 + 1296 * i] = 0.0;
        for (ix = 0; ix < 6; ix++) {
          a[i21 + 1296 * i] += b_a[i21 + 1296 * ix] * CH[ix + 6 * i];
        }
      }
    }

    b_sum(a, sum_jacobian);

    /* eq 17 -> P */
    d = 0.0;
    for (i21 = 0; i21 < 6; i21++) {
      sigma_12[i21] = 0.0;
      for (i = 0; i < 6; i++) {
        sigma_12[i21] += sigma_11[i21 + 6 * i] * sum_jacobian[i];
      }

      /* eq 17 -> P*J' */
      y[i21] = 0.0;
      for (i = 0; i < 6; i++) {
        y[i21] += sum_jacobian[i] * sigma_11[i + 6 * i21];
      }

      d += y[i21] * sum_jacobian[i21];
    }

    /* eq 17 -> J*P*J'+R */
    b = -norm(warped_grad_x3);
    for (i = 0; i < 6; i++) {
      sum_jacobian[i] = sigma_12[i] / (d + 1.0) * b;
    }

    for (i21 = 0; i21 < 9; i21++) {
      LOGMX[i21] = ((((sum_jacobian[0] * (real_T)g_b[i21] + sum_jacobian[1] *
                       (real_T)f_b[i21]) + sum_jacobian[2] * (real_T)e_b[i21]) +
                     sum_jacobian[3] * (real_T)d_b[i21]) + sum_jacobian[4] *
                    (real_T)c_b[i21]) + sum_jacobian[5] * (real_T)b_b[i21];
    }

    expm(LOGMX, b_LOGMX);
    for (i21 = 0; i21 < 6; i21++) {
      b_sigma_12[i21] = sigma_12[i21] / (d + 1.0);
    }

    for (i21 = 0; i21 < 6; i21++) {
      for (i = 0; i < 6; i++) {
        CH[i21 + 6 * i] = sigma_11[i21 + 6 * i] - b_sigma_12[i21] * sigma_12[i];
      }
    }

    chol(CH);
    randn(sigma_12);
    for (i21 = 0; i21 < 6; i21++) {
      sum_jacobian[i21] = 0.0;
      for (i = 0; i < 6; i++) {
        sum_jacobian[i21] += CH[i + 6 * i21] * sigma_12[i];
      }
    }

    for (i21 = 0; i21 < 9; i21++) {
      LOGMX[i21] = ((((sum_jacobian[0] * (real_T)g_b[i21] / 2.0 + sum_jacobian[1]
                       * (real_T)f_b[i21]) + sum_jacobian[2] * (real_T)e_b[i21])
                     + sum_jacobian[3] * (real_T)d_b[i21]) + sum_jacobian[4] *
                    (real_T)c_b[i21]) + sum_jacobian[5] * (real_T)b_b[i21];
    }

    expm(LOGMX, h_b);
    for (i21 = 0; i21 < 3; i21++) {
      for (i = 0; i < 3; i++) {
        LOGMX[i21 + 3 * i] = 0.0;
        for (ix = 0; ix < 3; ix++) {
          LOGMX[i21 + 3 * i] += X_par_temp[i21 + 3 * ix] * b_LOGMX[ix + 3 * i];
        }
      }
    }

    for (i21 = 0; i21 < 3; i21++) {
      for (i = 0; i < 3; i++) {
        X_par_pred[((i21 + 3 * i) + 9 * par) + 225 * ((int32_T)t - 1)] = 0.0;
        for (ix = 0; ix < 3; ix++) {
          X_par_pred[((i21 + 3 * i) + 9 * par) + 225 * ((int32_T)t - 1)] +=
            LOGMX[i21 + 3 * ix] * h_b[ix + 3 * i];
        }
      }
    }

    mldivide(*(real_T (*)[9])&X_par[9 * par + 225 * ((int32_T)(t - 1.0) - 1)],
             *(real_T (*)[9])&X_par_pred[9 * par + 225 * ((int32_T)t - 1)],
             LOGMX);
    for (i21 = 0; i21 < 3; i21++) {
      for (i = 0; i < 3; i++) {
        AR_velocity[((i + 3 * i21) + 9 * par) + 225 * ((int32_T)t - 1)] =
          LOGMX[i + 3 * i21];
      }
    }

    b_center_x[0] = center_x;
    b_center_x[1] = center_y;
    for (i21 = 0; i21 < 2; i21++) {
      for (i = 0; i < 2; i++) {
        Aff_matrix[i + 3 * i21] = X_par_pred[((i + 3 * i21) + 9 * par) + 225 *
          ((int32_T)t - 1)];
      }
    }

    for (i21 = 0; i21 < 2; i21++) {
      Aff_matrix[6 + i21] = X_par_pred[6 + ((i21 + 9 * par) + 225 * ((int32_T)t
        - 1))] + b_center_x[i21];
    }

    for (i21 = 0; i21 < 3; i21++) {
      Aff_matrix[2 + 3 * i21] = (real_T)iv9[i21];
    }

    b_image_warping(*(real_T (*)[336042])&ext_frame[0], Aff_matrix, point_matrix,
                    *(real_T (*)[1296])&warped_img[0]);
    b_image_warping(*(real_T (*)[336042])&ext_frame[336042], Aff_matrix,
                    point_matrix, *(real_T (*)[1296])&warped_img[1296]);
    b_image_warping(*(real_T (*)[336042])&ext_frame[672084], Aff_matrix,
                    point_matrix, *(real_T (*)[1296])&warped_img[2592]);

    /*             %% replace zeros with average values */
    for (i21 = 0; i21 < 1296; i21++) {
      b_zero_ind_z = (warped_img[2592 + i21] == 0.0);
      zero_ind_z[i21] = !b_zero_ind_z;
      c_zero_ind_z[i21] = b_zero_ind_z;
    }

    eml_li_find(c_zero_ind_z, tmp_data, tmp_size);
    eml_li_find(zero_ind_z, b_tmp_data, b_tmp_size);
    d_warped_img_size[0] = b_tmp_size[0];
    i = b_tmp_size[0];
    for (i21 = 0; i21 < i; i21++) {
      warped_grad_x3[i21] = warped_img[b_tmp_data[i21] + 2591];
    }

    d = mean(warped_grad_x3, d_warped_img_size);
    i = tmp_size[0];
    for (i21 = 0; i21 < i; i21++) {
      warped_img[tmp_data[i21] + 2591] = d;
    }

    eml_li_find(c_zero_ind_z, tmp_data, tmp_size);
    for (i = 0; i < 1296; i++) {
      zero_ind_z[i] = !c_zero_ind_z[i];
    }

    eml_li_find(zero_ind_z, b_tmp_data, b_tmp_size);
    e_warped_img_size[0] = b_tmp_size[0];
    i = b_tmp_size[0];
    for (i21 = 0; i21 < i; i21++) {
      warped_grad_x3[i21] = warped_img[b_tmp_data[i21] - 1];
    }

    d = mean(warped_grad_x3, e_warped_img_size);
    i = tmp_size[0];
    for (i21 = 0; i21 < i; i21++) {
      warped_img[tmp_data[i21] - 1] = d;
    }

    eml_li_find(c_zero_ind_z, tmp_data, tmp_size);
    for (i = 0; i < 1296; i++) {
      zero_ind_z[i] = !c_zero_ind_z[i];
    }

    eml_li_find(zero_ind_z, b_tmp_data, b_tmp_size);
    f_warped_img_size[0] = b_tmp_size[0];
    i = b_tmp_size[0];
    for (i21 = 0; i21 < i; i21++) {
      warped_grad_x3[i21] = warped_img[b_tmp_data[i21] + 1295];
    }

    d = mean(warped_grad_x3, f_warped_img_size);
    i = tmp_size[0];
    for (i21 = 0; i21 < i; i21++) {
      warped_img[tmp_data[i21] + 1295] = d;
    }

    /*             %% compute rigid transform */
    estimateRigidTransform(mean_img, warped_img, transform);
    for (i21 = 0; i21 < 4; i21++) {
      for (i = 0; i < 4; i++) {
        rigid_transform_par[(i + (i21 << 2)) + (par << 4)] = transform[i + (i21 <<
          2)];
      }
    }

    for (i21 = 0; i21 < 1296; i21++) {
      for (i = 0; i < 3; i++) {
        b_warped_img[i + (i21 << 2)] = warped_img[i21 + 1296 * i];
      }

      b_warped_img[3 + (i21 << 2)] = 1.0;
      for (i = 0; i < 4; i++) {
        c_transform[i21 + 1296 * i] = 0.0;
        for (ix = 0; ix < 4; ix++) {
          c_transform[i21 + 1296 * i] += transform[i + (ix << 2)] *
            b_warped_img[ix + (i21 << 2)];
        }
      }
    }

    for (i21 = 0; i21 < 3; i21++) {
      memcpy(&warped_img_tr[1296 * i21], &c_transform[1296 * i21], 1296U *
             sizeof(real_T));
    }

    for (i21 = 0; i21 < 1296; i21++) {
      b_warped_img_tr[i21] = warped_img_tr[2592 + i21] - mean_img[2592 + i21];
    }

    for (i21 = 0; i21 < 1296; i21++) {
      b_warped_img_tr[i21 + 1296] = warped_img_tr[1296 + i21] - mean_img[1296 +
        i21];
    }

    for (i21 = 0; i21 < 1296; i21++) {
      b_warped_img_tr[i21 + 2592] = warped_img_tr[i21] - mean_img[i21];
    }

    d = b_norm(b_warped_img_tr) / 3.0;

    /*             %% */
    dist1Array[par] = d * d;
  }

  for (i = 0; i < 25; i++) {
    b_dist1Array[i] = (dist1Array[i] > 1000.0);
  }

  b_eml_li_find(b_dist1Array, c_tmp_data, tmp_size);
  i = tmp_size[0];
  for (i21 = 0; i21 < i; i21++) {
    dist1Array[c_tmp_data[i21] - 1] = 1000.0;
  }

  for (i21 = 0; i21 < 25; i21++) {
    dist1Array[i21] *= -0.5;
  }

  b_exp(dist1Array);

  /* % Particle resampling */
  d = c_sum(dist1Array);
  for (i21 = 0; i21 < 25; i21++) {
    dist1Array[i21] /= d;
  }

  resampling(dist1Array, outindex);
  for (i21 = 0; i21 < 25; i21++) {
    for (i = 0; i < 3; i++) {
      for (ix = 0; ix < 3; ix++) {
        b_X_par_pred[(ix + 3 * i) + 9 * i21] = X_par_pred[((ix + 3 * i) + 9 *
          ((int32_T)outindex[i21] - 1)) + 225 * ((int32_T)t - 1)];
      }
    }
  }

  for (i21 = 0; i21 < 25; i21++) {
    for (i = 0; i < 3; i++) {
      for (ix = 0; ix < 3; ix++) {
        X_par[((ix + 3 * i) + 9 * i21) + 225 * ((int32_T)t - 1)] = b_X_par_pred
          [(ix + 3 * i) + 9 * i21];
      }
    }
  }

  for (i21 = 0; i21 < 25; i21++) {
    for (i = 0; i < 3; i++) {
      for (ix = 0; ix < 3; ix++) {
        b_X_par_pred[(ix + 3 * i) + 9 * i21] = AR_velocity[((ix + 3 * i) + 9 *
          ((int32_T)outindex[i21] - 1)) + 225 * ((int32_T)t - 1)];
      }
    }
  }

  for (i21 = 0; i21 < 25; i21++) {
    for (i = 0; i < 3; i++) {
      for (ix = 0; ix < 3; ix++) {
        AR_velocity[((ix + 3 * i) + 9 * i21) + 225 * ((int32_T)t - 1)] =
          b_X_par_pred[(ix + 3 * i) + 9 * i21];
      }
    }
  }

  for (i21 = 0; i21 < 25; i21++) {
    for (i = 0; i < 4; i++) {
      for (ix = 0; ix < 4; ix++) {
        b_rigid_transform_par[(ix + (i << 2)) + (i21 << 4)] =
          rigid_transform_par[(ix + (i << 2)) + (((int32_T)outindex[i21] - 1) <<
          4)];
      }
    }
  }

  for (i21 = 0; i21 < 25; i21++) {
    for (i = 0; i < 4; i++) {
      for (ix = 0; ix < 4; ix++) {
        rigid_transform_par[(ix + (i << 2)) + (i21 << 4)] =
          b_rigid_transform_par[(ix + (i << 2)) + (i21 << 4)];
      }
    }
  }

  /* % Computing affine state particle mean  */
  i = 1;
  d = dist1Array[0];
  itmp = 0;
  if (rtIsNaN(dist1Array[0])) {
    ix = 1;
    exitg1 = FALSE;
    while ((exitg1 == FALSE) && (ix + 1 < 26)) {
      i = ix + 1;
      if (!rtIsNaN(dist1Array[ix])) {
        d = dist1Array[ix];
        itmp = ix;
        exitg1 = TRUE;
      } else {
        ix++;
      }
    }
  }

  if (i < 25) {
    while (i + 1 < 26) {
      if (dist1Array[i] > d) {
        d = dist1Array[i];
        itmp = i;
      }

      i++;
    }
  }

  for (i21 = 0; i21 < 4; i21++) {
    mean_log[i21] = 0.0;
  }

  for (par = 0; par < 25; par++) {
    for (i21 = 0; i21 < 2; i21++) {
      for (i = 0; i < 2; i++) {
        c_X_par_pred[i + (i21 << 1)] = X_par_pred[((i + 3 * i21) + 9 * itmp) +
          225 * ((int32_T)t - 1)];
      }
    }

    for (i21 = 0; i21 < 2; i21++) {
      for (i = 0; i < 2; i++) {
        b_X_par[i + (i21 << 1)] = X_par[((i + 3 * i21) + 9 * par) + 225 *
          ((int32_T)t - 1)];
      }
    }

    b_mldivide(c_X_par_pred, b_X_par, c_LOGMX);
    b_eig(c_LOGMX, c_V, b_v);
    for (i21 = 0; i21 < 2; i21++) {
      dcv1[i21] = b_v[3 * i21];
    }

    c_log(dcv1);
    b_diag(dcv1, b_v);
    for (i21 = 0; i21 < 2; i21++) {
      for (i = 0; i < 2; i++) {
        d_V[i21 + (i << 1)].re = 0.0;
        d_V[i21 + (i << 1)].im = 0.0;
        for (ix = 0; ix < 2; ix++) {
          d_V[i21 + (i << 1)].re += c_V[i21 + (ix << 1)].re * b_v[ix + (i << 1)]
            .re - c_V[i21 + (ix << 1)].im * b_v[ix + (i << 1)].im;
          d_V[i21 + (i << 1)].im += c_V[i21 + (ix << 1)].re * b_v[ix + (i << 1)]
            .im + c_V[i21 + (ix << 1)].im * b_v[ix + (i << 1)].re;
        }
      }
    }

    c_mrdivide(d_V, c_V, b_v);
    for (i21 = 0; i21 < 4; i21++) {
      mean_log[i21] += b_v[i21].re / 25.0;
    }
  }

  b_expm(mean_log, c_LOGMX);
  for (i21 = 0; i21 < 25; i21++) {
    for (i = 0; i < 2; i++) {
      c_X_par[i + (i21 << 1)] = X_par[6 + ((i + 9 * i21) + 225 * ((int32_T)t - 1))];
    }
  }

  b_mean(c_X_par, mean_trans);

  /*  mean_X_par(:,:,t) = [mean_gl mean_trans;0 0 1]; */
  /* % Tracked object image update */
  c_center_x[0] = center_x;
  c_center_x[1] = center_y;
  for (i21 = 0; i21 < 2; i21++) {
    for (i = 0; i < 2; i++) {
      c_X_par_pred[i21 + (i << 1)] = 0.0;
      for (ix = 0; ix < 2; ix++) {
        c_X_par_pred[i21 + (i << 1)] += X_par_pred[((i21 + 3 * ix) + 9 * itmp) +
          225 * ((int32_T)t - 1)] * c_LOGMX[ix + (i << 1)];
      }
    }
  }

  for (i21 = 0; i21 < 2; i21++) {
    for (i = 0; i < 2; i++) {
      Aff_matrix[i + 3 * i21] = c_X_par_pred[i + (i21 << 1)];
    }
  }

  for (i21 = 0; i21 < 2; i21++) {
    Aff_matrix[6 + i21] = mean_trans[i21] + c_center_x[i21];
  }

  for (i21 = 0; i21 < 3; i21++) {
    Aff_matrix[2 + 3 * i21] = (real_T)iv9[i21];
  }

  /*  %% apply rigid transform to the tracked image and replace zeros with average values */
  memset(&warped_img[0], 0, 3888U * sizeof(real_T));
  b_image_warping(*(real_T (*)[336042])&ext_frame[0], Aff_matrix, point_matrix, *
                  (real_T (*)[1296])&warped_img[0]);
  b_image_warping(*(real_T (*)[336042])&ext_frame[336042], Aff_matrix,
                  point_matrix, *(real_T (*)[1296])&warped_img[1296]);
  b_image_warping(*(real_T (*)[336042])&ext_frame[672084], Aff_matrix,
                  point_matrix, *(real_T (*)[1296])&warped_img[2592]);
  for (i21 = 0; i21 < 1296; i21++) {
    b_zero_ind_z = (warped_img[2592 + i21] == 0.0);
    zero_ind_z[i21] = !b_zero_ind_z;
    c_zero_ind_z[i21] = b_zero_ind_z;
  }

  eml_li_find(c_zero_ind_z, tmp_data, tmp_size);
  eml_li_find(zero_ind_z, b_tmp_data, b_tmp_size);
  g_warped_img_size[0] = b_tmp_size[0];
  i = b_tmp_size[0];
  for (i21 = 0; i21 < i; i21++) {
    warped_grad_x3[i21] = warped_img[b_tmp_data[i21] + 2591];
  }

  d = mean(warped_grad_x3, g_warped_img_size);
  i = tmp_size[0];
  for (i21 = 0; i21 < i; i21++) {
    warped_img[tmp_data[i21] + 2591] = d;
  }

  eml_li_find(c_zero_ind_z, tmp_data, tmp_size);
  for (i = 0; i < 1296; i++) {
    zero_ind_z[i] = !c_zero_ind_z[i];
  }

  eml_li_find(zero_ind_z, b_tmp_data, b_tmp_size);
  h_warped_img_size[0] = b_tmp_size[0];
  i = b_tmp_size[0];
  for (i21 = 0; i21 < i; i21++) {
    warped_grad_x3[i21] = warped_img[b_tmp_data[i21] - 1];
  }

  d = mean(warped_grad_x3, h_warped_img_size);
  i = tmp_size[0];
  for (i21 = 0; i21 < i; i21++) {
    warped_img[tmp_data[i21] - 1] = d;
  }

  eml_li_find(c_zero_ind_z, tmp_data, tmp_size);
  for (i = 0; i < 1296; i++) {
    zero_ind_z[i] = !c_zero_ind_z[i];
  }

  eml_li_find(zero_ind_z, b_tmp_data, b_tmp_size);
  i_warped_img_size[0] = b_tmp_size[0];
  i = b_tmp_size[0];
  for (i21 = 0; i21 < i; i21++) {
    warped_grad_x3[i21] = warped_img[b_tmp_data[i21] + 1295];
  }

  d = mean(warped_grad_x3, i_warped_img_size);
  i = tmp_size[0];
  for (i21 = 0; i21 < i; i21++) {
    warped_img[tmp_data[i21] + 1295] = d;
  }

  for (i21 = 0; i21 < 1296; i21++) {
    for (i = 0; i < 3; i++) {
      b_warped_img[i + (i21 << 2)] = warped_img[i21 + 1296 * i];
    }

    b_warped_img[3 + (i21 << 2)] = 1.0;
    for (i = 0; i < 4; i++) {
      c_rigid_transform_par[i21 + 1296 * i] = 0.0;
      for (ix = 0; ix < 4; ix++) {
        c_rigid_transform_par[i21 + 1296 * i] += rigid_transform_par[(i + (ix <<
          2)) + (itmp << 4)] * b_warped_img[ix + (i21 << 2)];
      }
    }
  }

  for (i21 = 0; i21 < 3; i21++) {
    memcpy(&warped_img_tr[1296 * i21], &c_rigid_transform_par[1296 * i21], 1296U
           * sizeof(real_T));
  }

  for (i21 = 0; i21 < 1296; i21++) {
    tracked_images[i21 + 1296 * ((int32_T)t - 1)] = warped_img_tr[i21];
  }

  for (i21 = 0; i21 < 1296; i21++) {
    tracked_images[12960000 + (i21 + 1296 * ((int32_T)t - 1))] = warped_img_tr
      [1296 + i21];
  }

  for (i21 = 0; i21 < 1296; i21++) {
    tracked_images[25920000 + (i21 + 1296 * ((int32_T)t - 1))] = warped_img_tr
      [2592 + i21];
  }

  c_mean(warped_img, centroid);

  /*  if t >= init_size && rem(t,update_period) == 0 */
  /*          mean_img(:,1) = mean(tracked_images(:,1:t,1),2); */
  /*          mean_img(:,2) = mean(tracked_images(:,1:t,2),2); */
  /*          mean_img(:,3) = mean(tracked_images(:,1:t,3),2);       */
  /*          %% replace non connected values with mean values in the mean img */
  /*          mean_img = removeOutliers(mean_img); */
  /*  end */
}

void eml_li_find(const boolean_T x[1296], int32_T y_data[1296], int32_T y_size[1])
{
  int32_T j;
  int32_T i;
  y_size[0] = compute_nones(x);
  j = 0;
  for (i = 0; i < 1296; i++) {
    if (x[i]) {
      y_data[j] = i + 1;
      j++;
    }
  }
}

/* End of code generation (compute_prob1.cpp) */
