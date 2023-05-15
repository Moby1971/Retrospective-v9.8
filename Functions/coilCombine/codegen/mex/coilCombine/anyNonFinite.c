/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * anyNonFinite.c
 *
 * Code generation for function 'anyNonFinite'
 *
 */

/* Include files */
#include "anyNonFinite.h"
#include "coilCombine_data.h"
#include "coilCombine_types.h"
#include "eml_int_forloop_overflow_check.h"
#include "rt_nonfinite.h"
#include "mwmathutil.h"

/* Variable Definitions */
static emlrtRSInfo q_emlrtRSI = {
    29,             /* lineNo */
    "anyNonFinite", /* fcnName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/"
    "anyNonFinite.m" /* pathName */
};

/* Function Definitions */
boolean_T anyNonFinite(const emlrtStack *sp, const emxArray_creal_T *x)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack st;
  const creal_T *x_data;
  int32_T k;
  int32_T nx;
  boolean_T p;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  x_data = x->data;
  st.site = &q_emlrtRSI;
  b_st.site = &r_emlrtRSI;
  nx = x->size[0] * x->size[1];
  p = true;
  c_st.site = &s_emlrtRSI;
  if (nx > 2147483646) {
    d_st.site = &l_emlrtRSI;
    check_forloop_overflow_error(&d_st);
  }
  for (k = 0; k < nx; k++) {
    if ((!p) || (muDoubleScalarIsInf(x_data[k].re) ||
                 muDoubleScalarIsInf(x_data[k].im) ||
                 (muDoubleScalarIsNaN(x_data[k].re) ||
                  muDoubleScalarIsNaN(x_data[k].im)))) {
      p = false;
    }
  }
  return !p;
}

/* End of code generation (anyNonFinite.c) */
