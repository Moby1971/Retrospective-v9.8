/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * coilCombine_terminate.c
 *
 * Code generation for function 'coilCombine_terminate'
 *
 */

/* Include files */
#include "coilCombine_terminate.h"
#include "_coder_coilCombine_mex.h"
#include "coilCombine_data.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void coilCombine_atexit(void)
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

void coilCombine_terminate(void)
{
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (coilCombine_terminate.c) */
