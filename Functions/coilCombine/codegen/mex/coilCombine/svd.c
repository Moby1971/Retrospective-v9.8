/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * svd.c
 *
 * Code generation for function 'svd'
 *
 */

/* Include files */
#include "svd.h"
#include "coilCombine_data.h"
#include "coilCombine_emxutil.h"
#include "coilCombine_types.h"
#include "eml_int_forloop_overflow_check.h"
#include "rt_nonfinite.h"
#include "lapacke.h"
#include "mwmathutil.h"
#include <stddef.h>

/* Variable Definitions */
static emlrtRSInfo
    ub_emlrtRSI =
        {
            23,    /* lineNo */
            "svd", /* fcnName */
            "/Applications/MATLAB_R2023b.app/toolbox/eml/eml/+coder/+internal/"
            "svd.m" /* pathName */
};

static emlrtRSInfo
    vb_emlrtRSI =
        {
            52,    /* lineNo */
            "svd", /* fcnName */
            "/Applications/MATLAB_R2023b.app/toolbox/eml/eml/+coder/+internal/"
            "svd.m" /* pathName */
};

static emlrtRSInfo
    wb_emlrtRSI =
        {
            163,              /* lineNo */
            "getUSVForEmpty", /* fcnName */
            "/Applications/MATLAB_R2023b.app/toolbox/eml/eml/+coder/+internal/"
            "svd.m" /* pathName */
};

static emlrtRSInfo
    xb_emlrtRSI =
        {
            171,              /* lineNo */
            "getUSVForEmpty", /* fcnName */
            "/Applications/MATLAB_R2023b.app/toolbox/eml/eml/+coder/+internal/"
            "svd.m" /* pathName */
};

static emlrtRSInfo yb_emlrtRSI = {
    50,    /* lineNo */
    "eye", /* fcnName */
    "/Applications/MATLAB_R2023b.app/toolbox/eml/lib/matlab/elmat/eye.m" /* pathName
                                                                          */
};

static emlrtRSInfo ac_emlrtRSI = {
    96,    /* lineNo */
    "eye", /* fcnName */
    "/Applications/MATLAB_R2023b.app/toolbox/eml/lib/matlab/elmat/eye.m" /* pathName
                                                                          */
};

static emlrtRSInfo
    cc_emlrtRSI =
        {
            89,           /* lineNo */
            "callLAPACK", /* fcnName */
            "/Applications/MATLAB_R2023b.app/toolbox/eml/eml/+coder/+internal/"
            "svd.m" /* pathName */
};

static emlrtRSInfo
    dc_emlrtRSI =
        {
            81,           /* lineNo */
            "callLAPACK", /* fcnName */
            "/Applications/MATLAB_R2023b.app/toolbox/eml/eml/+coder/+internal/"
            "svd.m" /* pathName */
};

static emlrtRSInfo ec_emlrtRSI = {
    209,      /* lineNo */
    "xgesdd", /* fcnName */
    "/Applications/MATLAB_R2023b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "xgesdd.m" /* pathName */
};

static emlrtRSInfo fc_emlrtRSI = {
    31,       /* lineNo */
    "xgesvd", /* fcnName */
    "/Applications/MATLAB_R2023b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "xgesvd.m" /* pathName */
};

static emlrtRSInfo gc_emlrtRSI = {
    197,            /* lineNo */
    "ceval_xgesvd", /* fcnName */
    "/Applications/MATLAB_R2023b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "xgesvd.m" /* pathName */
};

static emlrtRTEInfo
    f_emlrtRTEI =
        {
            111,          /* lineNo */
            5,            /* colNo */
            "callLAPACK", /* fName */
            "/Applications/MATLAB_R2023b.app/toolbox/eml/eml/+coder/+internal/"
            "svd.m" /* pName */
};

static emlrtRTEInfo g_emlrtRTEI = {
    45,          /* lineNo */
    13,          /* colNo */
    "infocheck", /* fName */
    "/Applications/MATLAB_R2023b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "infocheck.m" /* pName */
};

static emlrtRTEInfo h_emlrtRTEI = {
    48,          /* lineNo */
    13,          /* colNo */
    "infocheck", /* fName */
    "/Applications/MATLAB_R2023b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "infocheck.m" /* pName */
};

static emlrtRTEInfo
    mb_emlrtRTEI =
        {
            57,    /* lineNo */
            33,    /* colNo */
            "svd", /* fName */
            "/Applications/MATLAB_R2023b.app/toolbox/eml/eml/+coder/+internal/"
            "svd.m" /* pName */
};

static emlrtRTEInfo
    nb_emlrtRTEI =
        {
            162,   /* lineNo */
            1,     /* colNo */
            "svd", /* fName */
            "/Applications/MATLAB_R2023b.app/toolbox/eml/eml/+coder/+internal/"
            "svd.m" /* pName */
};

static emlrtRTEInfo
    ob_emlrtRTEI =
        {
            81,    /* lineNo */
            63,    /* colNo */
            "svd", /* fName */
            "/Applications/MATLAB_R2023b.app/toolbox/eml/eml/+coder/+internal/"
            "svd.m" /* pName */
};

static emlrtRTEInfo pb_emlrtRTEI = {
    45,       /* lineNo */
    24,       /* colNo */
    "xgesdd", /* fName */
    "/Applications/MATLAB_R2023b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "xgesdd.m" /* pName */
};

static emlrtRTEInfo qb_emlrtRTEI = {
    47,       /* lineNo */
    25,       /* colNo */
    "xgesdd", /* fName */
    "/Applications/MATLAB_R2023b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "xgesdd.m" /* pName */
};

static emlrtRTEInfo rb_emlrtRTEI = {
    57,       /* lineNo */
    20,       /* colNo */
    "xgesdd", /* fName */
    "/Applications/MATLAB_R2023b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "xgesdd.m" /* pName */
};

static emlrtRTEInfo
    sb_emlrtRTEI =
        {
            171,   /* lineNo */
            9,     /* colNo */
            "svd", /* fName */
            "/Applications/MATLAB_R2023b.app/toolbox/eml/eml/+coder/+internal/"
            "svd.m" /* pName */
};

static emlrtRTEInfo tb_emlrtRTEI = {
    218,      /* lineNo */
    5,        /* colNo */
    "xgesvd", /* fName */
    "/Applications/MATLAB_R2023b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "xgesvd.m" /* pName */
};

static emlrtRTEInfo ub_emlrtRTEI = {
    75,       /* lineNo */
    24,       /* colNo */
    "xgesvd", /* fName */
    "/Applications/MATLAB_R2023b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "xgesvd.m" /* pName */
};

static emlrtRTEInfo vb_emlrtRTEI = {
    82,       /* lineNo */
    25,       /* colNo */
    "xgesvd", /* fName */
    "/Applications/MATLAB_R2023b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "xgesvd.m" /* pName */
};

static emlrtRTEInfo wb_emlrtRTEI = {
    90,       /* lineNo */
    20,       /* colNo */
    "xgesvd", /* fName */
    "/Applications/MATLAB_R2023b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "xgesvd.m" /* pName */
};

static emlrtRTEInfo xb_emlrtRTEI = {
    123,      /* lineNo */
    9,        /* colNo */
    "xgesvd", /* fName */
    "/Applications/MATLAB_R2023b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "xgesvd.m" /* pName */
};

static emlrtRTEInfo yb_emlrtRTEI = {
    121,      /* lineNo */
    33,       /* colNo */
    "xgesvd", /* fName */
    "/Applications/MATLAB_R2023b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "xgesvd.m" /* pName */
};

static emlrtRTEInfo
    ac_emlrtRTEI =
        {
            1,     /* lineNo */
            20,    /* colNo */
            "svd", /* fName */
            "/Applications/MATLAB_R2023b.app/toolbox/eml/eml/+coder/+internal/"
            "svd.m" /* pName */
};

static emlrtRTEInfo bc_emlrtRTEI = {
    82,       /* lineNo */
    5,        /* colNo */
    "xgesvd", /* fName */
    "/Applications/MATLAB_R2023b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "xgesvd.m" /* pName */
};

static emlrtRTEInfo cc_emlrtRTEI = {
    121,      /* lineNo */
    9,        /* colNo */
    "xgesvd", /* fName */
    "/Applications/MATLAB_R2023b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "xgesvd.m" /* pName */
};

/* Function Definitions */
void svd(const emlrtStack *sp, const emxArray_creal_T *A, emxArray_creal_T *U,
         emxArray_real_T *s, emxArray_creal_T *V)
{
  static const char_T b_fname[14] = {'L', 'A', 'P', 'A', 'C', 'K', 'E',
                                     '_', 'z', 'g', 'e', 's', 'v', 'd'};
  static const char_T fname[14] = {'L', 'A', 'P', 'A', 'C', 'K', 'E',
                                   '_', 'z', 'g', 'e', 's', 'd', 'd'};
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack st;
  emxArray_creal_T *Vt;
  emxArray_creal_T *b_A;
  emxArray_creal_T *c_A;
  emxArray_real_T *superb;
  const creal_T *A_data;
  creal_T *U_data;
  creal_T *Vt_data;
  creal_T *b_A_data;
  creal_T *c_A_data;
  real_T *s_data;
  real_T *superb_data;
  int32_T i;
  int32_T m;
  int32_T n;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  A_data = A->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  if ((A->size[0] == 0) || (A->size[1] == 0)) {
    int32_T info;
    int32_T loop_ub;
    st.site = &ub_emlrtRSI;
    n = A->size[0];
    m = A->size[1];
    i = U->size[0] * U->size[1];
    U->size[0] = A->size[0];
    U->size[1] = A->size[0];
    emxEnsureCapacity_creal_T(&st, U, i, &nb_emlrtRTEI);
    U_data = U->data;
    loop_ub = A->size[0] * A->size[0];
    for (i = 0; i < loop_ub; i++) {
      U_data[i].re = 0.0;
      U_data[i].im = 0.0;
    }
    info = muIntScalarMin_sint32(n, n);
    b_st.site = &wb_emlrtRSI;
    if (info > 2147483646) {
      c_st.site = &l_emlrtRSI;
      check_forloop_overflow_error(&c_st);
    }
    for (n = 0; n < info; n++) {
      U_data[n + U->size[0] * n].re = 1.0;
      U_data[n + U->size[0] * n].im = 0.0;
    }
    b_st.site = &xb_emlrtRSI;
    c_st.site = &yb_emlrtRSI;
    i = V->size[0] * V->size[1];
    V->size[0] = A->size[1];
    V->size[1] = A->size[1];
    emxEnsureCapacity_creal_T(&b_st, V, i, &sb_emlrtRTEI);
    U_data = V->data;
    loop_ub = A->size[1] * A->size[1];
    for (i = 0; i < loop_ub; i++) {
      U_data[i].re = 0.0;
      U_data[i].im = 0.0;
    }
    if (A->size[1] > 0) {
      c_st.site = &ac_emlrtRSI;
      if (A->size[1] > 2147483646) {
        d_st.site = &l_emlrtRSI;
        check_forloop_overflow_error(&d_st);
      }
      for (n = 0; n < m; n++) {
        U_data[n + V->size[0] * n].re = 1.0;
        U_data[n + V->size[0] * n].im = 0.0;
      }
    }
    s->size[0] = 0;
  } else {
    ptrdiff_t info_t;
    int32_T info;
    st.site = &vb_emlrtRSI;
    emxInit_creal_T(&st, &b_A, 2, &ac_emlrtRTEI);
    i = b_A->size[0] * b_A->size[1];
    b_A->size[0] = A->size[0];
    b_A->size[1] = A->size[1];
    emxEnsureCapacity_creal_T(&st, b_A, i, &mb_emlrtRTEI);
    b_A_data = b_A->data;
    info = A->size[0] * A->size[1];
    for (i = 0; i < info; i++) {
      b_A_data[i] = A_data[i];
    }
    m = A->size[0];
    n = A->size[1];
    b_st.site = &dc_emlrtRSI;
    emxInit_creal_T(&b_st, &c_A, 2, &ob_emlrtRTEI);
    i = c_A->size[0] * c_A->size[1];
    c_A->size[0] = A->size[0];
    c_A->size[1] = A->size[1];
    emxEnsureCapacity_creal_T(&b_st, c_A, i, &ob_emlrtRTEI);
    c_A_data = c_A->data;
    for (i = 0; i < info; i++) {
      c_A_data[i] = A_data[i];
    }
    i = U->size[0] * U->size[1];
    U->size[0] = A->size[0];
    U->size[1] = A->size[0];
    emxEnsureCapacity_creal_T(&b_st, U, i, &pb_emlrtRTEI);
    U_data = U->data;
    emxInit_creal_T(&b_st, &Vt, 2, &bc_emlrtRTEI);
    i = Vt->size[0] * Vt->size[1];
    Vt->size[0] = A->size[1];
    Vt->size[1] = A->size[1];
    emxEnsureCapacity_creal_T(&b_st, Vt, i, &qb_emlrtRTEI);
    Vt_data = Vt->data;
    i = muIntScalarMin_sint32(n, m);
    m = s->size[0];
    s->size[0] = i;
    emxEnsureCapacity_real_T(&b_st, s, m, &rb_emlrtRTEI);
    s_data = s->data;
    info_t = LAPACKE_zgesdd(
        102, 'A', (ptrdiff_t)A->size[0], (ptrdiff_t)A->size[1],
        (lapack_complex_double *)&c_A_data[0], (ptrdiff_t)A->size[0],
        &s_data[0], (lapack_complex_double *)&U_data[0], (ptrdiff_t)A->size[0],
        (lapack_complex_double *)&Vt_data[0], (ptrdiff_t)A->size[1]);
    emxFree_creal_T(&b_st, &c_A);
    c_st.site = &ec_emlrtRSI;
    if ((int32_T)info_t < 0) {
      if ((int32_T)info_t == -1010) {
        emlrtErrorWithMessageIdR2018a(&c_st, &g_emlrtRTEI, "MATLAB:nomem",
                                      "MATLAB:nomem", 0);
      } else {
        emlrtErrorWithMessageIdR2018a(&c_st, &h_emlrtRTEI,
                                      "Coder:toolbox:LAPACKCallErrorInfo",
                                      "Coder:toolbox:LAPACKCallErrorInfo", 5, 4,
                                      14, &fname[0], 12, (int32_T)info_t);
      }
    }
    info = (int32_T)info_t;
    if ((int32_T)info_t > 0) {
      int32_T loop_ub;
      b_st.site = &cc_emlrtRSI;
      c_st.site = &fc_emlrtRSI;
      m = U->size[0] * U->size[1];
      U->size[0] = A->size[0];
      U->size[1] = A->size[0];
      emxEnsureCapacity_creal_T(&c_st, U, m, &ub_emlrtRTEI);
      U_data = U->data;
      m = Vt->size[0] * Vt->size[1];
      Vt->size[0] = A->size[1];
      Vt->size[1] = A->size[1];
      emxEnsureCapacity_creal_T(&c_st, Vt, m, &vb_emlrtRTEI);
      Vt_data = Vt->data;
      m = s->size[0];
      s->size[0] = i;
      emxEnsureCapacity_real_T(&c_st, s, m, &wb_emlrtRTEI);
      s_data = s->data;
      emxInit_real_T(&c_st, &superb, &cc_emlrtRTEI);
      if (i > 1) {
        m = superb->size[0];
        superb->size[0] = i - 1;
        emxEnsureCapacity_real_T(&c_st, superb, m, &yb_emlrtRTEI);
        superb_data = superb->data;
      } else {
        i = superb->size[0];
        superb->size[0] = 1;
        emxEnsureCapacity_real_T(&c_st, superb, i, &xb_emlrtRTEI);
        superb_data = superb->data;
      }
      info_t = LAPACKE_zgesvd(
          102, 'A', 'A', (ptrdiff_t)A->size[0], (ptrdiff_t)A->size[1],
          (lapack_complex_double *)&b_A_data[0], (ptrdiff_t)A->size[0],
          &s_data[0], (lapack_complex_double *)&U_data[0],
          (ptrdiff_t)A->size[0], (lapack_complex_double *)&Vt_data[0],
          (ptrdiff_t)A->size[1], &superb_data[0]);
      emxFree_real_T(&c_st, &superb);
      i = V->size[0] * V->size[1];
      V->size[0] = Vt->size[1];
      V->size[1] = Vt->size[0];
      emxEnsureCapacity_creal_T(&c_st, V, i, &tb_emlrtRTEI);
      U_data = V->data;
      loop_ub = Vt->size[0];
      for (i = 0; i < loop_ub; i++) {
        n = Vt->size[1];
        for (m = 0; m < n; m++) {
          U_data[m + V->size[0] * i].re = Vt_data[i + Vt->size[0] * m].re;
          U_data[m + V->size[0] * i].im = -Vt_data[i + Vt->size[0] * m].im;
        }
      }
      d_st.site = &gc_emlrtRSI;
      if ((int32_T)info_t < 0) {
        if ((int32_T)info_t == -1010) {
          emlrtErrorWithMessageIdR2018a(&d_st, &g_emlrtRTEI, "MATLAB:nomem",
                                        "MATLAB:nomem", 0);
        } else {
          emlrtErrorWithMessageIdR2018a(
              &d_st, &h_emlrtRTEI, "Coder:toolbox:LAPACKCallErrorInfo",
              "Coder:toolbox:LAPACKCallErrorInfo", 5, 4, 14, &b_fname[0], 12,
              (int32_T)info_t);
        }
      }
      info = (int32_T)info_t;
    } else {
      int32_T loop_ub;
      i = V->size[0] * V->size[1];
      V->size[0] = Vt->size[1];
      V->size[1] = Vt->size[0];
      emxEnsureCapacity_creal_T(&st, V, i, &tb_emlrtRTEI);
      U_data = V->data;
      loop_ub = Vt->size[0];
      for (i = 0; i < loop_ub; i++) {
        n = Vt->size[1];
        for (m = 0; m < n; m++) {
          U_data[m + V->size[0] * i].re = Vt_data[i + Vt->size[0] * m].re;
          U_data[m + V->size[0] * i].im = -Vt_data[i + Vt->size[0] * m].im;
        }
      }
    }
    emxFree_creal_T(&st, &Vt);
    emxFree_creal_T(&st, &b_A);
    if (info > 0) {
      emlrtErrorWithMessageIdR2018a(&st, &f_emlrtRTEI,
                                    "Coder:MATLAB:svd_NoConvergence",
                                    "Coder:MATLAB:svd_NoConvergence", 0);
    }
  }
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

/* End of code generation (svd.c) */
