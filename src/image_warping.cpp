/*
 * image_warping.cpp
 *
 * Code generation for function 'image_warping'
 *
 * C source code generated on: Tue Apr  1 12:30:12 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Optimal_affine_tracking_3d16_fast_realtime.h"
#include "image_warping.h"
#include "round.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
void b_image_warping(const real_T frame[336042], const real_T aff_matrix[9],
                     const real_T point_matrix[3888], real_T warped_img[1296])
{
  int32_T k;
  int32_T i10;
  real_T transformed_point[3888];
  int32_T i11;
  real_T maxval[1296];
  real_T b_maxval[1296];
  real_T u0;
  for (k = 0; k < 3; k++) {
    for (i10 = 0; i10 < 1296; i10++) {
      transformed_point[k + 3 * i10] = 0.0;
      for (i11 = 0; i11 < 3; i11++) {
        transformed_point[k + 3 * i10] += aff_matrix[k + 3 * i11] *
          point_matrix[i11 + 3 * i10];
      }
    }
  }

  for (k = 0; k < 1296; k++) {
    u0 = transformed_point[3 * k];
    if (u0 <= 381.0) {
    } else {
      u0 = 381.0;
    }

    if ((1.0 >= u0) || rtIsNaN(u0)) {
      u0 = 1.0;
    }

    maxval[k] = u0;
    u0 = transformed_point[1 + 3 * k];
    if (u0 <= 882.0) {
    } else {
      u0 = 882.0;
    }

    if ((1.0 >= u0) || rtIsNaN(u0)) {
      u0 = 1.0;
    }

    b_maxval[k] = u0;
  }

  b_round(maxval);
  b_round(b_maxval);
  for (k = 0; k < 1296; k++) {
    warped_img[k] = frame[(int32_T)((b_maxval[k] - 1.0) * 381.0 + maxval[k]) - 1];
  }

  /*  warped_img = reshape(frame((y-1)*size_x + x),[crop_x crop_y]); */
}

void image_warping(const real_T frame[168021], const real_T aff_matrix[9], const
                   real_T point_matrix[3888], real_T warped_img[1296])
{
  int32_T k;
  int32_T i5;
  real_T transformed_point[3888];
  int32_T i6;
  real_T maxval[1296];
  real_T b_maxval[1296];
  real_T u0;
  for (k = 0; k < 3; k++) {
    for (i5 = 0; i5 < 1296; i5++) {
      transformed_point[k + 3 * i5] = 0.0;
      for (i6 = 0; i6 < 3; i6++) {
        transformed_point[k + 3 * i5] += aff_matrix[k + 3 * i6] *
          point_matrix[i6 + 3 * i5];
      }
    }
  }

  for (k = 0; k < 1296; k++) {
    u0 = transformed_point[3 * k];
    if (u0 <= 381.0) {
    } else {
      u0 = 381.0;
    }

    if ((1.0 >= u0) || rtIsNaN(u0)) {
      u0 = 1.0;
    }

    maxval[k] = u0;
    u0 = transformed_point[1 + 3 * k];
    if (u0 <= 441.0) {
    } else {
      u0 = 441.0;
    }

    if ((1.0 >= u0) || rtIsNaN(u0)) {
      u0 = 1.0;
    }

    b_maxval[k] = u0;
  }

  b_round(maxval);
  b_round(b_maxval);
  for (k = 0; k < 1296; k++) {
    warped_img[k] = frame[(int32_T)((b_maxval[k] - 1.0) * 381.0 + maxval[k]) - 1];
  }

  /*  warped_img = reshape(frame((y-1)*size_x + x),[crop_x crop_y]); */
}

/* End of code generation (image_warping.cpp) */
