/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_coilCombine_api.c
 *
 * Code generation for function '_coder_coilCombine_api'
 *
 */

/* Include files */
#include "_coder_coilCombine_api.h"
#include "coilCombine.h"
#include "coilCombine_data.h"
#include "coilCombine_emxutil.h"
#include "coilCombine_types.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRTEInfo bc_emlrtRTEI = {
    1,                        /* lineNo */
    1,                        /* colNo */
    "_coder_coilCombine_api", /* fName */
    ""                        /* pName */
};

/* Function Declarations */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               emxArray_creal_T *y);

static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               emxArray_creal_T *ret);

static void emlrt_marshallIn(const emlrtStack *sp, const mxArray *im1,
                             const char_T *identifier, emxArray_creal_T *y);

static const mxArray *emlrt_marshallOut(const emlrtStack *sp,
                                        const emxArray_creal_T *u);

/* Function Definitions */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               emxArray_creal_T *y)
{
  c_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               emxArray_creal_T *ret)
{
  static const int32_T dims[5] = {-1, -1, -1, -1, -1};
  creal_T *ret_data;
  int32_T iv[5];
  int32_T i;
  const boolean_T bv[5] = {true, true, true, true, true};
  emlrtCheckVsBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", true, 5U,
                            (const void *)&dims[0], &bv[0], &iv[0]);
  i = ret->size[0] * ret->size[1] * ret->size[2] * ret->size[3] * ret->size[4];
  ret->size[0] = iv[0];
  ret->size[1] = iv[1];
  ret->size[2] = iv[2];
  ret->size[3] = iv[3];
  ret->size[4] = iv[4];
  emxEnsureCapacity_creal_T(sp, ret, i, (emlrtRTEInfo *)NULL);
  ret_data = ret->data;
  emlrtImportArrayR2015b((emlrtConstCTX)sp, src, &ret_data[0], 8, true);
  emlrtDestroyArray(&src);
}

static void emlrt_marshallIn(const emlrtStack *sp, const mxArray *im1,
                             const char_T *identifier, emxArray_creal_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  b_emlrt_marshallIn(sp, emlrtAlias(im1), &thisId, y);
  emlrtDestroyArray(&im1);
}

static const mxArray *emlrt_marshallOut(const emlrtStack *sp,
                                        const emxArray_creal_T *u)
{
  const mxArray *m;
  const mxArray *y;
  const creal_T *u_data;
  int32_T iv[4];
  u_data = u->data;
  y = NULL;
  iv[0] = u->size[0];
  iv[1] = u->size[1];
  iv[2] = u->size[2];
  iv[3] = u->size[3];
  m = emlrtCreateNumericArray(4, &iv[0], mxDOUBLE_CLASS, mxCOMPLEX);
  emlrtExportNumericArrayR2013b((emlrtConstCTX)sp, m, (const void *)&u_data[0],
                                8);
  emlrtAssign(&y, m);
  return y;
}

void coilCombine_api(const mxArray *prhs, const mxArray **plhs)
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  emxArray_creal_T *im1;
  emxArray_creal_T *im2;
  st.tls = emlrtRootTLSGlobal;
  emlrtHeapReferenceStackEnterFcnR2012b(&st);
  /* Marshall function inputs */
  emxInit_creal_T(&st, &im1, 5, &bc_emlrtRTEI);
  emlrt_marshallIn(&st, emlrtAliasP(prhs), "im1", im1);
  /* Invoke the target function */
  emxInit_creal_T(&st, &im2, 4, &bc_emlrtRTEI);
  coilCombine(&st, im1, im2);
  emxFree_creal_T(&st, &im1);
  /* Marshall function outputs */
  *plhs = emlrt_marshallOut(&st, im2);
  emxFree_creal_T(&st, &im2);
  emlrtHeapReferenceStackLeaveFcnR2012b(&st);
}

/* End of code generation (_coder_coilCombine_api.c) */
