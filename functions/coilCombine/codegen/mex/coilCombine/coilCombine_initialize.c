/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * coilCombine_initialize.c
 *
 * Code generation for function 'coilCombine_initialize'
 *
 */

/* Include files */
#include "coilCombine_initialize.h"
#include "_coder_coilCombine_mex.h"
#include "coilCombine_data.h"
#include "rt_nonfinite.h"

/* Function Declarations */
static void coilCombine_once(void);

/* Function Definitions */
static void coilCombine_once(void)
{
  mex_InitInfAndNan();
}

void coilCombine_initialize(void)
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtBreakCheckR2012bFlagVar = emlrtGetBreakCheckFlagAddressR2022b(&st);
  emlrtClearAllocCountR2012b(&st, false, 0U, NULL);
  emlrtEnterRtStackR2012b(&st);
  if (emlrtFirstTimeR2012b(emlrtRootTLSGlobal)) {
    coilCombine_once();
  }
}

/* End of code generation (coilCombine_initialize.c) */
