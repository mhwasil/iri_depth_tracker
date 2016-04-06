/*
 * rand.cpp
 *
 * Code generation for function 'rand'
 *
 * C source code generated on: Tue Apr  1 12:30:12 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Optimal_affine_tracking_3d16_fast_realtime.h"
#include "rand.h"
#include "randn.h"
#include "Optimal_affine_tracking_3d16_fast_realtime_data.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
static void eml_rand_mt19937ar_stateful(real_T r[1843200]);

/* Function Definitions */
static void eml_rand_mt19937ar_stateful(real_T r[1843200])
{
  int32_T k;
  real_T d0;
  for (k = 0; k < 1843200; k++) {
    d0 = genrandu(state);
    r[k] = d0;
  }
}

void b_rand(real_T r[1843200])
{
  if (!method_not_empty) {
    method_not_empty = TRUE;
  }

  eml_rand_mt19937ar_stateful(r);
}

real_T c_rand(void)
{
  if (!method_not_empty) {
    method_not_empty = TRUE;
  }

  return genrandu(state);
}

/* End of code generation (rand.cpp) */
