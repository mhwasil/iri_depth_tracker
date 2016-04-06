/*
 * Optimal_affine_tracking_3d16_fast_realtime_initialize.cpp
 *
 * Code generation for function 'Optimal_affine_tracking_3d16_fast_realtime_initialize'
 *
 * C source code generated on: Tue Apr  1 12:30:12 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Optimal_affine_tracking_3d16_fast_realtime.h"
#include "Optimal_affine_tracking_3d16_fast_realtime_initialize.h"
#include "Optimal_affine_tracking_3d16_fast_realtime_data.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
void Optimal_affine_tracking_3d16_fast_realtime_initialize(void)
{
  uint32_T r;
  int32_T mti;
  rt_InitInfAndNaN(8U);
  method_not_empty = FALSE;
  memset(&state[0], 0, 625U * sizeof(uint32_T));
  r = 5489U;
  state[0] = 5489U;
  for (mti = 0; mti < 623; mti++) {
    r = (r ^ r >> 30U) * 1812433253U + (uint32_T)(1 + mti);
    state[1 + mti] = r;
  }

  memset(&state[0], 0, 625U * sizeof(uint32_T));
  r = 5489U;
  state[0] = 5489U;
  for (mti = 0; mti < 623; mti++) {
    r = (r ^ r >> 30U) * 1812433253U + (uint32_T)(1 + mti);
    state[1 + mti] = r;
  }

  memset(&state[0], 0, 625U * sizeof(uint32_T));
  r = 5489U;
  state[0] = 5489U;
  for (mti = 0; mti < 623; mti++) {
    r = (r ^ r >> 30U) * 1812433253U + (uint32_T)(1 + mti);
    state[1 + mti] = r;
  }

  state[624] = 624U;
}

/* End of code generation (Optimal_affine_tracking_3d16_fast_realtime_initialize.cpp) */
