/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * permute.c
 *
 * Code generation for function 'permute'
 *
 */

/* Include files */
#include "permute.h"
#include "coilCombine_data.h"
#include "coilCombine_emxutil.h"
#include "coilCombine_types.h"
#include "eml_int_forloop_overflow_check.h"
#include "reshapeSizeChecks.h"
#include "rt_nonfinite.h"
#include "mwmathutil.h"

/* Variable Definitions */
static emlrtRSInfo
    e_emlrtRSI =
        {
            44,        /* lineNo */
            "permute", /* fcnName */
            "/Applications/MATLAB_R2022b.app/toolbox/eml/lib/matlab/elmat/"
            "permute.m" /* pathName */
};

static emlrtRSInfo
    f_emlrtRSI =
        {
            47,        /* lineNo */
            "permute", /* fcnName */
            "/Applications/MATLAB_R2022b.app/toolbox/eml/lib/matlab/elmat/"
            "permute.m" /* pathName */
};

static emlrtRSInfo
    g_emlrtRSI =
        {
            53,        /* lineNo */
            "permute", /* fcnName */
            "/Applications/MATLAB_R2022b.app/toolbox/eml/lib/matlab/elmat/"
            "permute.m" /* pathName */
};

static emlrtRSInfo h_emlrtRSI = {
    51,                  /* lineNo */
    "reshapeSizeChecks", /* fcnName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/"
    "reshapeSizeChecks.m" /* pathName */
};

static emlrtRSInfo i_emlrtRSI = {
    119,               /* lineNo */
    "computeDimsData", /* fcnName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/"
    "reshapeSizeChecks.m" /* pathName */
};

static emlrtRSInfo
    j_emlrtRSI =
        {
            71,       /* lineNo */
            "looper", /* fcnName */
            "/Applications/MATLAB_R2022b.app/toolbox/eml/lib/matlab/elmat/"
            "permute.m" /* pathName */
};

static emlrtRSInfo
    k_emlrtRSI =
        {
            72,       /* lineNo */
            "looper", /* fcnName */
            "/Applications/MATLAB_R2022b.app/toolbox/eml/lib/matlab/elmat/"
            "permute.m" /* pathName */
};

static emlrtRTEInfo
    fb_emlrtRTEI =
        {
            47,        /* lineNo */
            20,        /* colNo */
            "permute", /* fName */
            "/Applications/MATLAB_R2022b.app/toolbox/eml/lib/matlab/elmat/"
            "permute.m" /* pName */
};

static emlrtRTEInfo
    gb_emlrtRTEI =
        {
            44,        /* lineNo */
            5,         /* colNo */
            "permute", /* fName */
            "/Applications/MATLAB_R2022b.app/toolbox/eml/lib/matlab/elmat/"
            "permute.m" /* pName */
};

/* Function Definitions */
void b_permute(const emlrtStack *sp, const emxArray_creal_T *a,
               emxArray_creal_T *b)
{
  static const int8_T iv[5] = {1, 2, 5, 4, 3};
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack g_st;
  emlrtStack st;
  const creal_T *a_data;
  creal_T *b_data;
  int32_T b_k;
  int32_T c_k;
  int32_T d_k;
  int32_T k;
  int32_T nx;
  int32_T plast;
  int32_T subsa_idx_1;
  int32_T subsa_idx_3;
  int32_T subsa_idx_4;
  boolean_T b_b;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  e_st.prev = &d_st;
  e_st.tls = d_st.tls;
  f_st.prev = &e_st;
  f_st.tls = e_st.tls;
  g_st.prev = &f_st;
  g_st.tls = f_st.tls;
  a_data = a->data;
  b_b = true;
  if ((a->size[0] != 0) && (a->size[1] != 0) && (a->size[3] != 0) &&
      (a->size[4] != 0)) {
    boolean_T exitg1;
    plast = 0;
    k = 0;
    exitg1 = false;
    while ((!exitg1) && (k < 5)) {
      if (a->size[iv[k] - 1] != 1) {
        if (plast > iv[k]) {
          b_b = false;
          exitg1 = true;
        } else {
          plast = iv[k];
          k++;
        }
      } else {
        k++;
      }
    }
  }
  if (b_b) {
    st.site = &e_emlrtRSI;
    nx = a->size[0] * a->size[1] * a->size[3] * a->size[4];
    b_st.site = &h_emlrtRSI;
    c_st.site = &i_emlrtRSI;
    c_st.site = &i_emlrtRSI;
    c_st.site = &i_emlrtRSI;
    c_st.site = &i_emlrtRSI;
    plast = a->size[0];
    if (a->size[1] > a->size[0]) {
      plast = a->size[1];
    }
    if (plast < 1) {
      plast = 1;
    }
    if (a->size[3] > plast) {
      plast = a->size[3];
    }
    if (a->size[4] > plast) {
      plast = a->size[4];
    }
    plast = muIntScalarMax_sint32(nx, plast);
    if (a->size[0] > plast) {
      emlrtErrorWithMessageIdR2018a(
          &st, &e_emlrtRTEI, "Coder:toolbox:reshape_emptyReshapeLimit",
          "Coder:toolbox:reshape_emptyReshapeLimit", 0);
    }
    if (a->size[1] > plast) {
      emlrtErrorWithMessageIdR2018a(
          &st, &e_emlrtRTEI, "Coder:toolbox:reshape_emptyReshapeLimit",
          "Coder:toolbox:reshape_emptyReshapeLimit", 0);
    }
    if (a->size[4] > plast) {
      emlrtErrorWithMessageIdR2018a(
          &st, &e_emlrtRTEI, "Coder:toolbox:reshape_emptyReshapeLimit",
          "Coder:toolbox:reshape_emptyReshapeLimit", 0);
    }
    if (a->size[3] > plast) {
      emlrtErrorWithMessageIdR2018a(
          &st, &e_emlrtRTEI, "Coder:toolbox:reshape_emptyReshapeLimit",
          "Coder:toolbox:reshape_emptyReshapeLimit", 0);
    }
    if (a->size[0] * a->size[1] * a->size[4] * a->size[3] != nx) {
      emlrtErrorWithMessageIdR2018a(
          &st, &d_emlrtRTEI, "Coder:MATLAB:getReshapeDims_notSameNumel",
          "Coder:MATLAB:getReshapeDims_notSameNumel", 0);
    }
    nx = b->size[0] * b->size[1] * b->size[2] * b->size[3];
    b->size[0] = a->size[0];
    b->size[1] = a->size[1];
    b->size[2] = a->size[4];
    b->size[3] = a->size[3];
    emxEnsureCapacity_creal_T(sp, b, nx, &gb_emlrtRTEI);
    b_data = b->data;
    plast = a->size[0] * a->size[1] * a->size[4] * a->size[3];
    for (nx = 0; nx < plast; nx++) {
      b_data[nx] = a_data[nx];
    }
  } else {
    st.site = &f_emlrtRSI;
    nx = a->size[0] * a->size[1] * a->size[3] * a->size[4];
    b_st.site = &h_emlrtRSI;
    c_st.site = &i_emlrtRSI;
    c_st.site = &i_emlrtRSI;
    c_st.site = &i_emlrtRSI;
    c_st.site = &i_emlrtRSI;
    plast = a->size[0];
    if (a->size[1] > a->size[0]) {
      plast = a->size[1];
    }
    if (plast < 1) {
      plast = 1;
    }
    if (a->size[3] > plast) {
      plast = a->size[3];
    }
    if (a->size[4] > plast) {
      plast = a->size[4];
    }
    plast = muIntScalarMax_sint32(nx, plast);
    if (a->size[0] > plast) {
      emlrtErrorWithMessageIdR2018a(
          &st, &e_emlrtRTEI, "Coder:toolbox:reshape_emptyReshapeLimit",
          "Coder:toolbox:reshape_emptyReshapeLimit", 0);
    }
    if (a->size[1] > plast) {
      emlrtErrorWithMessageIdR2018a(
          &st, &e_emlrtRTEI, "Coder:toolbox:reshape_emptyReshapeLimit",
          "Coder:toolbox:reshape_emptyReshapeLimit", 0);
    }
    if (a->size[4] > plast) {
      emlrtErrorWithMessageIdR2018a(
          &st, &e_emlrtRTEI, "Coder:toolbox:reshape_emptyReshapeLimit",
          "Coder:toolbox:reshape_emptyReshapeLimit", 0);
    }
    if (a->size[3] > plast) {
      emlrtErrorWithMessageIdR2018a(
          &st, &e_emlrtRTEI, "Coder:toolbox:reshape_emptyReshapeLimit",
          "Coder:toolbox:reshape_emptyReshapeLimit", 0);
    }
    if (a->size[0] * a->size[1] * a->size[4] * a->size[3] != nx) {
      emlrtErrorWithMessageIdR2018a(
          &st, &d_emlrtRTEI, "Coder:MATLAB:getReshapeDims_notSameNumel",
          "Coder:MATLAB:getReshapeDims_notSameNumel", 0);
    }
    nx = b->size[0] * b->size[1] * b->size[2] * b->size[3];
    b->size[0] = a->size[0];
    b->size[1] = a->size[1];
    b->size[2] = a->size[4];
    b->size[3] = a->size[3];
    emxEnsureCapacity_creal_T(sp, b, nx, &fb_emlrtRTEI);
    b_data = b->data;
    st.site = &g_emlrtRSI;
    plast = a->size[4];
    b_st.site = &j_emlrtRSI;
    if (a->size[4] > 2147483646) {
      c_st.site = &l_emlrtRSI;
      check_forloop_overflow_error(&c_st);
    }
    for (k = 0; k < plast; k++) {
      b_st.site = &k_emlrtRSI;
      nx = a->size[3];
      c_st.site = &j_emlrtRSI;
      if (a->size[3] > 2147483646) {
        d_st.site = &l_emlrtRSI;
        check_forloop_overflow_error(&d_st);
      }
      for (b_k = 0; b_k < nx; b_k++) {
        int32_T c_b;
        c_st.site = &k_emlrtRSI;
        d_st.site = &k_emlrtRSI;
        c_b = a->size[1];
        e_st.site = &j_emlrtRSI;
        if (a->size[1] > 2147483646) {
          f_st.site = &l_emlrtRSI;
          check_forloop_overflow_error(&f_st);
        }
        for (c_k = 0; c_k < c_b; c_k++) {
          int32_T d_b;
          e_st.site = &k_emlrtRSI;
          d_b = a->size[0];
          f_st.site = &j_emlrtRSI;
          if (a->size[0] > 2147483646) {
            g_st.site = &l_emlrtRSI;
            check_forloop_overflow_error(&g_st);
          }
          if (a->size[0] - 1 >= 0) {
            subsa_idx_1 = c_k + 1;
            subsa_idx_3 = b_k + 1;
            subsa_idx_4 = k + 1;
          }
          for (d_k = 0; d_k < d_b; d_k++) {
            b_data[((d_k + b->size[0] * (subsa_idx_1 - 1)) +
                    b->size[0] * b->size[1] * (subsa_idx_4 - 1)) +
                   b->size[0] * b->size[1] * b->size[2] * (subsa_idx_3 - 1)] =
                a_data[((d_k + a->size[0] * (subsa_idx_1 - 1)) +
                        a->size[0] * a->size[1] * (subsa_idx_3 - 1)) +
                       a->size[0] * a->size[1] * a->size[3] *
                           (subsa_idx_4 - 1)];
          }
        }
      }
    }
  }
}

void c_permute(const emlrtStack *sp, const emxArray_creal_T *a,
               emxArray_creal_T *b)
{
  static const int8_T iv[5] = {1, 2, 3, 5, 4};
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack g_st;
  emlrtStack st;
  const creal_T *a_data;
  creal_T *b_data;
  int32_T b_k;
  int32_T c_k;
  int32_T d_k;
  int32_T k;
  int32_T nx;
  int32_T plast;
  int32_T subsa_idx_1;
  int32_T subsa_idx_2;
  int32_T subsa_idx_4;
  boolean_T b_b;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  e_st.prev = &d_st;
  e_st.tls = d_st.tls;
  f_st.prev = &e_st;
  f_st.tls = e_st.tls;
  g_st.prev = &f_st;
  g_st.tls = f_st.tls;
  a_data = a->data;
  b_b = true;
  if ((a->size[0] != 0) && (a->size[1] != 0) && (a->size[2] != 0) &&
      (a->size[4] != 0)) {
    boolean_T exitg1;
    plast = 0;
    k = 0;
    exitg1 = false;
    while ((!exitg1) && (k < 5)) {
      if (a->size[iv[k] - 1] != 1) {
        if (plast > iv[k]) {
          b_b = false;
          exitg1 = true;
        } else {
          plast = iv[k];
          k++;
        }
      } else {
        k++;
      }
    }
  }
  if (b_b) {
    st.site = &e_emlrtRSI;
    nx = a->size[0] * a->size[1] * a->size[2] * a->size[4];
    b_st.site = &h_emlrtRSI;
    c_st.site = &i_emlrtRSI;
    c_st.site = &i_emlrtRSI;
    c_st.site = &i_emlrtRSI;
    c_st.site = &i_emlrtRSI;
    plast = a->size[0];
    if (a->size[1] > a->size[0]) {
      plast = a->size[1];
    }
    if (a->size[2] > plast) {
      plast = a->size[2];
    }
    if (plast < 1) {
      plast = 1;
    }
    if (a->size[4] > plast) {
      plast = a->size[4];
    }
    plast = muIntScalarMax_sint32(nx, plast);
    if (a->size[0] > plast) {
      emlrtErrorWithMessageIdR2018a(
          &st, &e_emlrtRTEI, "Coder:toolbox:reshape_emptyReshapeLimit",
          "Coder:toolbox:reshape_emptyReshapeLimit", 0);
    }
    if (a->size[1] > plast) {
      emlrtErrorWithMessageIdR2018a(
          &st, &e_emlrtRTEI, "Coder:toolbox:reshape_emptyReshapeLimit",
          "Coder:toolbox:reshape_emptyReshapeLimit", 0);
    }
    if (a->size[2] > plast) {
      emlrtErrorWithMessageIdR2018a(
          &st, &e_emlrtRTEI, "Coder:toolbox:reshape_emptyReshapeLimit",
          "Coder:toolbox:reshape_emptyReshapeLimit", 0);
    }
    if (a->size[4] > plast) {
      emlrtErrorWithMessageIdR2018a(
          &st, &e_emlrtRTEI, "Coder:toolbox:reshape_emptyReshapeLimit",
          "Coder:toolbox:reshape_emptyReshapeLimit", 0);
    }
    if (a->size[0] * a->size[1] * a->size[2] * a->size[4] != nx) {
      emlrtErrorWithMessageIdR2018a(
          &st, &d_emlrtRTEI, "Coder:MATLAB:getReshapeDims_notSameNumel",
          "Coder:MATLAB:getReshapeDims_notSameNumel", 0);
    }
    nx = b->size[0] * b->size[1] * b->size[2] * b->size[3];
    b->size[0] = a->size[0];
    b->size[1] = a->size[1];
    b->size[2] = a->size[2];
    b->size[3] = a->size[4];
    emxEnsureCapacity_creal_T(sp, b, nx, &gb_emlrtRTEI);
    b_data = b->data;
    plast = a->size[0] * a->size[1] * a->size[2] * a->size[4];
    for (nx = 0; nx < plast; nx++) {
      b_data[nx] = a_data[nx];
    }
  } else {
    st.site = &f_emlrtRSI;
    nx = a->size[0] * a->size[1] * a->size[2] * a->size[4];
    b_st.site = &h_emlrtRSI;
    c_st.site = &i_emlrtRSI;
    c_st.site = &i_emlrtRSI;
    c_st.site = &i_emlrtRSI;
    c_st.site = &i_emlrtRSI;
    plast = a->size[0];
    if (a->size[1] > a->size[0]) {
      plast = a->size[1];
    }
    if (a->size[2] > plast) {
      plast = a->size[2];
    }
    if (plast < 1) {
      plast = 1;
    }
    if (a->size[4] > plast) {
      plast = a->size[4];
    }
    plast = muIntScalarMax_sint32(nx, plast);
    if (a->size[0] > plast) {
      emlrtErrorWithMessageIdR2018a(
          &st, &e_emlrtRTEI, "Coder:toolbox:reshape_emptyReshapeLimit",
          "Coder:toolbox:reshape_emptyReshapeLimit", 0);
    }
    if (a->size[1] > plast) {
      emlrtErrorWithMessageIdR2018a(
          &st, &e_emlrtRTEI, "Coder:toolbox:reshape_emptyReshapeLimit",
          "Coder:toolbox:reshape_emptyReshapeLimit", 0);
    }
    if (a->size[2] > plast) {
      emlrtErrorWithMessageIdR2018a(
          &st, &e_emlrtRTEI, "Coder:toolbox:reshape_emptyReshapeLimit",
          "Coder:toolbox:reshape_emptyReshapeLimit", 0);
    }
    if (a->size[4] > plast) {
      emlrtErrorWithMessageIdR2018a(
          &st, &e_emlrtRTEI, "Coder:toolbox:reshape_emptyReshapeLimit",
          "Coder:toolbox:reshape_emptyReshapeLimit", 0);
    }
    if (a->size[0] * a->size[1] * a->size[2] * a->size[4] != nx) {
      emlrtErrorWithMessageIdR2018a(
          &st, &d_emlrtRTEI, "Coder:MATLAB:getReshapeDims_notSameNumel",
          "Coder:MATLAB:getReshapeDims_notSameNumel", 0);
    }
    nx = b->size[0] * b->size[1] * b->size[2] * b->size[3];
    b->size[0] = a->size[0];
    b->size[1] = a->size[1];
    b->size[2] = a->size[2];
    b->size[3] = a->size[4];
    emxEnsureCapacity_creal_T(sp, b, nx, &fb_emlrtRTEI);
    b_data = b->data;
    st.site = &g_emlrtRSI;
    plast = a->size[4];
    b_st.site = &j_emlrtRSI;
    if (a->size[4] > 2147483646) {
      c_st.site = &l_emlrtRSI;
      check_forloop_overflow_error(&c_st);
    }
    for (k = 0; k < plast; k++) {
      b_st.site = &k_emlrtRSI;
      c_st.site = &k_emlrtRSI;
      nx = a->size[2];
      d_st.site = &j_emlrtRSI;
      if (a->size[2] > 2147483646) {
        e_st.site = &l_emlrtRSI;
        check_forloop_overflow_error(&e_st);
      }
      for (b_k = 0; b_k < nx; b_k++) {
        int32_T c_b;
        d_st.site = &k_emlrtRSI;
        c_b = a->size[1];
        e_st.site = &j_emlrtRSI;
        if (a->size[1] > 2147483646) {
          f_st.site = &l_emlrtRSI;
          check_forloop_overflow_error(&f_st);
        }
        for (c_k = 0; c_k < c_b; c_k++) {
          int32_T d_b;
          e_st.site = &k_emlrtRSI;
          d_b = a->size[0];
          f_st.site = &j_emlrtRSI;
          if (a->size[0] > 2147483646) {
            g_st.site = &l_emlrtRSI;
            check_forloop_overflow_error(&g_st);
          }
          if (a->size[0] - 1 >= 0) {
            subsa_idx_1 = c_k + 1;
            subsa_idx_2 = b_k + 1;
            subsa_idx_4 = k + 1;
          }
          for (d_k = 0; d_k < d_b; d_k++) {
            b_data[((d_k + b->size[0] * (subsa_idx_1 - 1)) +
                    b->size[0] * b->size[1] * (subsa_idx_2 - 1)) +
                   b->size[0] * b->size[1] * b->size[2] * (subsa_idx_4 - 1)] =
                a_data[((d_k + a->size[0] * (subsa_idx_1 - 1)) +
                        a->size[0] * a->size[1] * (subsa_idx_2 - 1)) +
                       a->size[0] * a->size[1] * a->size[2] *
                           (subsa_idx_4 - 1)];
          }
        }
      }
    }
  }
}

void permute(const emlrtStack *sp, const emxArray_creal_T *a,
             emxArray_creal_T *b)
{
  static const int8_T iv[5] = {1, 2, 3, 5, 4};
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack g_st;
  emlrtStack st;
  const creal_T *a_data;
  creal_T *b_data;
  int32_T b_k;
  int32_T c_k;
  int32_T d_k;
  int32_T e_k;
  int32_T k;
  int32_T nx;
  int32_T plast;
  int32_T subsa_idx_1;
  int32_T subsa_idx_2;
  int32_T subsa_idx_3;
  int32_T subsa_idx_4;
  boolean_T b_b;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  e_st.prev = &d_st;
  e_st.tls = d_st.tls;
  f_st.prev = &e_st;
  f_st.tls = e_st.tls;
  g_st.prev = &f_st;
  g_st.tls = f_st.tls;
  a_data = a->data;
  b_b = true;
  if ((a->size[0] != 0) && (a->size[1] != 0) && (a->size[2] != 0) &&
      (a->size[3] != 0) && (a->size[4] != 0)) {
    boolean_T exitg1;
    plast = 0;
    k = 0;
    exitg1 = false;
    while ((!exitg1) && (k < 5)) {
      if (a->size[iv[k] - 1] != 1) {
        if (plast > iv[k]) {
          b_b = false;
          exitg1 = true;
        } else {
          plast = iv[k];
          k++;
        }
      } else {
        k++;
      }
    }
  }
  if (b_b) {
    st.site = &e_emlrtRSI;
    nx = a->size[0] * a->size[1] * a->size[2] * a->size[3] * a->size[4];
    b_st.site = &h_emlrtRSI;
    computeDimsData();
    plast = a->size[0];
    if (a->size[1] > a->size[0]) {
      plast = a->size[1];
    }
    if (a->size[2] > plast) {
      plast = a->size[2];
    }
    if (a->size[3] > plast) {
      plast = a->size[3];
    }
    if (a->size[4] > plast) {
      plast = a->size[4];
    }
    plast = muIntScalarMax_sint32(nx, plast);
    if (a->size[0] > plast) {
      emlrtErrorWithMessageIdR2018a(
          &st, &e_emlrtRTEI, "Coder:toolbox:reshape_emptyReshapeLimit",
          "Coder:toolbox:reshape_emptyReshapeLimit", 0);
    }
    if (a->size[1] > plast) {
      emlrtErrorWithMessageIdR2018a(
          &st, &e_emlrtRTEI, "Coder:toolbox:reshape_emptyReshapeLimit",
          "Coder:toolbox:reshape_emptyReshapeLimit", 0);
    }
    if (a->size[2] > plast) {
      emlrtErrorWithMessageIdR2018a(
          &st, &e_emlrtRTEI, "Coder:toolbox:reshape_emptyReshapeLimit",
          "Coder:toolbox:reshape_emptyReshapeLimit", 0);
    }
    if (a->size[4] > plast) {
      emlrtErrorWithMessageIdR2018a(
          &st, &e_emlrtRTEI, "Coder:toolbox:reshape_emptyReshapeLimit",
          "Coder:toolbox:reshape_emptyReshapeLimit", 0);
    }
    if (a->size[3] > plast) {
      emlrtErrorWithMessageIdR2018a(
          &st, &e_emlrtRTEI, "Coder:toolbox:reshape_emptyReshapeLimit",
          "Coder:toolbox:reshape_emptyReshapeLimit", 0);
    }
    if (a->size[0] * a->size[1] * a->size[2] * a->size[4] * a->size[3] != nx) {
      emlrtErrorWithMessageIdR2018a(
          &st, &d_emlrtRTEI, "Coder:MATLAB:getReshapeDims_notSameNumel",
          "Coder:MATLAB:getReshapeDims_notSameNumel", 0);
    }
    nx = b->size[0] * b->size[1] * b->size[2] * b->size[3] * b->size[4];
    b->size[0] = a->size[0];
    b->size[1] = a->size[1];
    b->size[2] = a->size[2];
    b->size[3] = a->size[4];
    b->size[4] = a->size[3];
    emxEnsureCapacity_creal_T(sp, b, nx, &gb_emlrtRTEI);
    b_data = b->data;
    plast = a->size[0] * a->size[1] * a->size[2] * a->size[4] * a->size[3];
    for (nx = 0; nx < plast; nx++) {
      b_data[nx] = a_data[nx];
    }
  } else {
    st.site = &f_emlrtRSI;
    nx = a->size[0] * a->size[1] * a->size[2] * a->size[3] * a->size[4];
    b_st.site = &h_emlrtRSI;
    computeDimsData();
    plast = a->size[0];
    if (a->size[1] > a->size[0]) {
      plast = a->size[1];
    }
    if (a->size[2] > plast) {
      plast = a->size[2];
    }
    if (a->size[3] > plast) {
      plast = a->size[3];
    }
    if (a->size[4] > plast) {
      plast = a->size[4];
    }
    plast = muIntScalarMax_sint32(nx, plast);
    if (a->size[0] > plast) {
      emlrtErrorWithMessageIdR2018a(
          &st, &e_emlrtRTEI, "Coder:toolbox:reshape_emptyReshapeLimit",
          "Coder:toolbox:reshape_emptyReshapeLimit", 0);
    }
    if (a->size[1] > plast) {
      emlrtErrorWithMessageIdR2018a(
          &st, &e_emlrtRTEI, "Coder:toolbox:reshape_emptyReshapeLimit",
          "Coder:toolbox:reshape_emptyReshapeLimit", 0);
    }
    if (a->size[2] > plast) {
      emlrtErrorWithMessageIdR2018a(
          &st, &e_emlrtRTEI, "Coder:toolbox:reshape_emptyReshapeLimit",
          "Coder:toolbox:reshape_emptyReshapeLimit", 0);
    }
    if (a->size[4] > plast) {
      emlrtErrorWithMessageIdR2018a(
          &st, &e_emlrtRTEI, "Coder:toolbox:reshape_emptyReshapeLimit",
          "Coder:toolbox:reshape_emptyReshapeLimit", 0);
    }
    if (a->size[3] > plast) {
      emlrtErrorWithMessageIdR2018a(
          &st, &e_emlrtRTEI, "Coder:toolbox:reshape_emptyReshapeLimit",
          "Coder:toolbox:reshape_emptyReshapeLimit", 0);
    }
    if (a->size[0] * a->size[1] * a->size[2] * a->size[4] * a->size[3] != nx) {
      emlrtErrorWithMessageIdR2018a(
          &st, &d_emlrtRTEI, "Coder:MATLAB:getReshapeDims_notSameNumel",
          "Coder:MATLAB:getReshapeDims_notSameNumel", 0);
    }
    nx = b->size[0] * b->size[1] * b->size[2] * b->size[3] * b->size[4];
    b->size[0] = a->size[0];
    b->size[1] = a->size[1];
    b->size[2] = a->size[2];
    b->size[3] = a->size[4];
    b->size[4] = a->size[3];
    emxEnsureCapacity_creal_T(sp, b, nx, &fb_emlrtRTEI);
    b_data = b->data;
    st.site = &g_emlrtRSI;
    plast = a->size[4];
    b_st.site = &j_emlrtRSI;
    if (a->size[4] > 2147483646) {
      c_st.site = &l_emlrtRSI;
      check_forloop_overflow_error(&c_st);
    }
    for (k = 0; k < plast; k++) {
      b_st.site = &k_emlrtRSI;
      nx = a->size[3];
      c_st.site = &j_emlrtRSI;
      if (a->size[3] > 2147483646) {
        d_st.site = &l_emlrtRSI;
        check_forloop_overflow_error(&d_st);
      }
      for (b_k = 0; b_k < nx; b_k++) {
        int32_T c_b;
        c_st.site = &k_emlrtRSI;
        c_b = a->size[2];
        d_st.site = &j_emlrtRSI;
        if (a->size[2] > 2147483646) {
          e_st.site = &l_emlrtRSI;
          check_forloop_overflow_error(&e_st);
        }
        for (c_k = 0; c_k < c_b; c_k++) {
          int32_T d_b;
          d_st.site = &k_emlrtRSI;
          d_b = a->size[1];
          e_st.site = &j_emlrtRSI;
          if (a->size[1] > 2147483646) {
            f_st.site = &l_emlrtRSI;
            check_forloop_overflow_error(&f_st);
          }
          for (d_k = 0; d_k < d_b; d_k++) {
            int32_T e_b;
            e_st.site = &k_emlrtRSI;
            e_b = a->size[0];
            f_st.site = &j_emlrtRSI;
            if (a->size[0] > 2147483646) {
              g_st.site = &l_emlrtRSI;
              check_forloop_overflow_error(&g_st);
            }
            if (a->size[0] - 1 >= 0) {
              subsa_idx_1 = d_k + 1;
              subsa_idx_2 = c_k + 1;
              subsa_idx_3 = b_k + 1;
              subsa_idx_4 = k + 1;
            }
            for (e_k = 0; e_k < e_b; e_k++) {
              b_data[(((e_k + b->size[0] * (subsa_idx_1 - 1)) +
                       b->size[0] * b->size[1] * (subsa_idx_2 - 1)) +
                      b->size[0] * b->size[1] * b->size[2] *
                          (subsa_idx_4 - 1)) +
                     b->size[0] * b->size[1] * b->size[2] * b->size[3] *
                         (subsa_idx_3 - 1)] =
                  a_data[(((e_k + a->size[0] * (subsa_idx_1 - 1)) +
                           a->size[0] * a->size[1] * (subsa_idx_2 - 1)) +
                          a->size[0] * a->size[1] * a->size[2] *
                              (subsa_idx_3 - 1)) +
                         a->size[0] * a->size[1] * a->size[2] * a->size[3] *
                             (subsa_idx_4 - 1)];
            }
          }
        }
      }
    }
  }
}

/* End of code generation (permute.c) */
