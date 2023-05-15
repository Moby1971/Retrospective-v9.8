/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * squeeze.c
 *
 * Code generation for function 'squeeze'
 *
 */

/* Include files */
#include "squeeze.h"
#include "coilCombine_data.h"
#include "coilCombine_emxutil.h"
#include "coilCombine_types.h"
#include "rt_nonfinite.h"
#include "mwmathutil.h"

/* Variable Definitions */
static emlrtRSInfo
    qb_emlrtRSI =
        {
            38,        /* lineNo */
            "squeeze", /* fcnName */
            "/Applications/MATLAB_R2022b.app/toolbox/eml/lib/matlab/elmat/"
            "squeeze.m" /* pathName */
};

static emlrtRTEInfo
    jb_emlrtRTEI =
        {
            38,        /* lineNo */
            1,         /* colNo */
            "squeeze", /* fName */
            "/Applications/MATLAB_R2022b.app/toolbox/eml/lib/matlab/elmat/"
            "squeeze.m" /* pName */
};

/* Function Definitions */
void squeeze(const emlrtStack *sp, const emxArray_creal_T *a,
             emxArray_creal_T *b)
{
  emlrtStack st;
  const creal_T *a_data;
  creal_T *b_data;
  int32_T szb[2];
  int32_T j;
  int32_T nx;
  boolean_T p;
  st.prev = sp;
  st.tls = sp->tls;
  a_data = a->data;
  szb[0] = 1;
  szb[1] = 1;
  p = (a->size[2] == 1);
  if ((!p) || (a->size[3] != 1)) {
    p = false;
  }
  if (!p) {
    j = 0;
    if (a->size[2] != 1) {
      j = 1;
      szb[0] = a->size[2];
    }
    if (a->size[3] != 1) {
      szb[j] = a->size[3];
    }
  }
  st.site = &qb_emlrtRSI;
  nx = a->size[2] * a->size[3];
  j = 1;
  if (a->size[2] > 1) {
    j = a->size[2];
  }
  if (a->size[3] > j) {
    j = a->size[3];
  }
  j = muIntScalarMax_sint32(nx, j);
  if (szb[0] > j) {
    emlrtErrorWithMessageIdR2018a(&st, &e_emlrtRTEI,
                                  "Coder:toolbox:reshape_emptyReshapeLimit",
                                  "Coder:toolbox:reshape_emptyReshapeLimit", 0);
  }
  if (szb[1] > j) {
    emlrtErrorWithMessageIdR2018a(&st, &e_emlrtRTEI,
                                  "Coder:toolbox:reshape_emptyReshapeLimit",
                                  "Coder:toolbox:reshape_emptyReshapeLimit", 0);
  }
  j = szb[0] * szb[1];
  if (j != nx) {
    emlrtErrorWithMessageIdR2018a(
        &st, &d_emlrtRTEI, "Coder:MATLAB:getReshapeDims_notSameNumel",
        "Coder:MATLAB:getReshapeDims_notSameNumel", 0);
  }
  nx = b->size[0] * b->size[1];
  b->size[0] = szb[0];
  b->size[1] = szb[1];
  emxEnsureCapacity_creal_T(sp, b, nx, &jb_emlrtRTEI);
  b_data = b->data;
  for (nx = 0; nx < j; nx++) {
    b_data[nx] = a_data[nx];
  }
}

/* End of code generation (squeeze.c) */
