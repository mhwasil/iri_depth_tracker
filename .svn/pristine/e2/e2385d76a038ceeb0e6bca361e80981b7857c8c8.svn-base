/*
 * Optimal_affine_tracking_3d16_fast_realtime.cpp
 *
 * Code generation for function 'Optimal_affine_tracking_3d16_fast_realtime'
 *
 * C source code generated on: Tue Apr  1 12:30:12 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Optimal_affine_tracking_3d16_fast_realtime.h"
#include "compute_prob1.h"
#include "LoadKinectMesh_realtime.h"
#include "init_variables.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
void Optimal_affine_tracking_3d16_fast_realtime(const real_T pty[4], const
  real_T ptx[4], real_T p[12])
{
  static real_T Irgb[504063];
  static real_T Ixyz[504063];
  real_T center_y;
  real_T center_x;
  real_T corner_p[12];
  static real_T mean_img[3888];
  static real_T point_matrix[3888];
  static real_T AR_velocity[2250000];
  static real_T X_par[2250000];
  static real_T dw_dp[15552];
  static real_T tracked_images[38880000];
  static real_T X_par_pred[2250000];
  int32_T t;
  real_T centroid[3];
  real_T Aff_matrix[9];
  int32_T i0;
  int32_T i1;
  int32_T i2;

  /* % Object template initialization */
  /*  [Ixyz,Irgb] = LoadKinectMeshBinary_realtime; */
  LoadKinectMesh_realtime(Ixyz, Irgb);

  /*  h1=figure(1); */
  init_variables(Ixyz, ptx, pty, X_par_pred, tracked_images, dw_dp, X_par,
                 AR_velocity, point_matrix, mean_img, corner_p, &center_x,
                 &center_y);

  /* % Execute tracking */
  for (t = 0; t < 9999; t++) {
    LoadKinectMesh_realtime(Ixyz, Irgb);

    /*      [Ixyz,Irgb] = LoadKinectMeshBinary_realtime;     */
    /*     %% Particle propagation and likelihood computation  */
    compute_prob1(X_par, X_par_pred, AR_velocity, dw_dp, (real_T)t + 2.0,
                  center_x, center_y, Ixyz, point_matrix, mean_img,
                  tracked_images, Aff_matrix, centroid);

    /*     %% Display the tracking results */
    /*      set(0,'CurrentFigure',h1); */
    /*      imshow(Irgb/255); */
    for (i0 = 0; i0 < 3; i0++) {
      for (i1 = 0; i1 < 4; i1++) {
        p[i0 + 3 * i1] = 0.0;
        for (i2 = 0; i2 < 3; i2++) {
          p[i0 + 3 * i1] += Aff_matrix[i0 + 3 * i2] * corner_p[i2 + 3 * i1];
        }
      }
    }

    /*      hold on;line([p(2,1) p(2,2)],[p(1,1) p(1,2)],'Color','red','LineWidth',3); */
    /*      hold on;line([p(2,2) p(2,3)],[p(1,2) p(1,3)],'Color','red','LineWidth',3); */
    /*      hold on;line([p(2,3) p(2,4)],[p(1,3) p(1,4)],'Color','red','LineWidth',3); */
    /*      hold on;line([p(2,4) p(2,1)],[p(1,4) p(1,1)],'Color','red','LineWidth',3); */
    /*      hold off; */
    /*      drawnow; */
  }
}

/* End of code generation (Optimal_affine_tracking_3d16_fast_realtime.cpp) */
