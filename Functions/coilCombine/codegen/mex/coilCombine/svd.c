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
            "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/"
            "svd.m" /* pathName */
};

static emlrtRSInfo
    vb_emlrtRSI =
        {
            52,    /* lineNo */
            "svd", /* fcnName */
            "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/"
            "svd.m" /* pathName */
};

static emlrtRSInfo
    wb_emlrtRSI =
        {
            163,              /* lineNo */
            "getUSVForEmpty", /* fcnName */
            "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/"
            "svd.m" /* pathName */
};

static emlrtRSInfo
    xb_emlrtRSI =
        {
            171,              /* lineNo */
            "getUSVForEmpty", /* fcnName */
            "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/"
            "svd.m" /* pathName */
};

static emlrtRSInfo yb_emlrtRSI = {
    50,    /* lineNo */
    "eye", /* fcnName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/lib/matlab/elmat/eye.m" /* pathName
                                                                          */
};

static emlrtRSInfo ac_emlrtRSI = {
    96,    /* lineNo */
    "eye", /* fcnName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/lib/matlab/elmat/eye.m" /* pathName
                                                                          */
};

static emlrtRSInfo
    cc_emlrtRSI =
        {
            89,           /* lineNo */
            "callLAPACK", /* fcnName */
            "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/"
            "svd.m" /* pathName */
};

static emlrtRSInfo
    dc_emlrtRSI =
        {
            81,           /* lineNo */
            "callLAPACK", /* fcnName */
            "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/"
            "svd.m" /* pathName */
};

static emlrtRSInfo ec_emlrtRSI = {
    209,      /* lineNo */
    "xgesdd", /* fcnName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "xgesdd.m" /* pathName */
};

static emlrtRSInfo fc_emlrtRSI = {
    31,       /* lineNo */
    "xgesvd", /* fcnName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "xgesvd.m" /* pathName */
};

static emlrtRSInfo gc_emlrtRSI = {
    197,            /* lineNo */
    "ceval_xgesvd", /* fcnName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "xgesvd.m" /* pathName */
};

static emlrtRTEInfo
    f_emlrtRTEI =
        {
            111,          /* lineNo */
            5,            /* colNo */
            "callLAPACK", /* fName */
            "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/"
            "svd.m" /* pName */
};

static emlrtRTEInfo g_emlrtRTEI = {
    44,          /* lineNo */
    13,          /* colNo */
    "infocheck", /* fName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "infocheck.m" /* pName */
};

static emlrtRTEInfo h_emlrtRTEI = {
    47,          /* lineNo */
    13,          /* colNo */
    "infocheck", /* fName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "infocheck.m" /* pName */
};

static emlrtRTEInfo
    kb_emlrtRTEI =
        {
            57,    /* lineNo */
            33,    /* colNo */
            "svd", /* fName */
            "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/"
            "svd.m" /* pName */
};

static emlrtRTEInfo
    lb_emlrtRTEI =
        {
            162,   /* lineNo */
            1,     /* colNo */
            "svd", /* fName */
            "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/"
            "svd.m" /* pName */
};

static emlrtRTEInfo
    mb_emlrtRTEI =
        {
            81,    /* lineNo */
            63,    /* colNo */
            "svd", /* fName */
            "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/"
            "svd.m" /* pName */
};

static emlrtRTEInfo nb_emlrtRTEI = {
    45,       /* lineNo */
    24,       /* colNo */
    "xgesdd", /* fName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "xgesdd.m" /* pName */
};

static emlrtRTEInfo ob_emlrtRTEI = {
    47,       /* lineNo */
    25,       /* colNo */
    "xgesdd", /* fName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "xgesdd.m" /* pName */
};

static emlrtRTEInfo pb_emlrtRTEI = {
    57,       /* lineNo */
    20,       /* colNo */
    "xgesdd", /* fName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "xgesdd.m" /* pName */
};

static emlrtRTEInfo
    qb_emlrtRTEI =
        {
            171,   /* lineNo */
            9,     /* colNo */
            "svd", /* fName */
            "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/"
            "svd.m" /* pName */
};

static emlrtRTEInfo rb_emlrtRTEI = {
    218,      /* lineNo */
    5,        /* colNo */
    "xgesvd", /* fName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "xgesvd.m" /* pName */
};

static emlrtRTEInfo sb_emlrtRTEI = {
    75,       /* lineNo */
    24,       /* colNo */
    "xgesvd", /* fName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "xgesvd.m" /* pName */
};

static emlrtRTEInfo tb_emlrtRTEI = {
    82,       /* lineNo */
    25,       /* colNo */
    "xgesvd", /* fName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "xgesvd.m" /* pName */
};

static emlrtRTEInfo ub_emlrtRTEI = {
    90,       /* lineNo */
    20,       /* colNo */
    "xgesvd", /* fName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "xgesvd.m" /* pName */
};

static emlrtRTEInfo vb_emlrtRTEI = {
    123,      /* lineNo */
    9,        /* colNo */
    "xgesvd", /* fName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "xgesvd.m" /* pName */
};

static emlrtRTEInfo wb_emlrtRTEI = {
    121,      /* lineNo */
    33,       /* colNo */
    "xgesvd", /* fName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "xgesvd.m" /* pName */
};

static emlrtRTEInfo
    xb_emlrtRTEI =
        {
            1,     /* lineNo */
            20,    /* colNo */
            "svd", /* fName */
            "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/"
            "svd.m" /* pName */
};

static emlrtRTEInfo yb_emlrtRTEI = {
    82,       /* lineNo */
    5,        /* colNo */
    "xgesvd", /* fName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
    "xgesvd.m" /* pName */
};

static emlrtRTEInfo ac_emlrtRTEI = {
    121,      /* lineNo */
    9,        /* colNo */
    "xgesvd", /* fName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/+lapack/"
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
  int32_T i1;
  int32_T info;
  int32_T loop_ub;
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
    int32_T m;
    st.site = &ub_emlrtRSI;
    info = A->size[0];
    m = A->size[1];
    i = U->size[0] * U->size[1];
    U->size[0] = A->size[0];
    U->size[1] = A->size[0];
    emxEnsureCapacity_creal_T(&st, U, i, &lb_emlrtRTEI);
    U_data = U->data;
    loop_ub = A->size[0] * A->size[0];
    for (i = 0; i < loop_ub; i++) {
      U_data[i].re = 0.0;
      U_data[i].im = 0.0;
    }
    info = muIntScalarMin_sint32(info, info);
    b_st.site = &wb_emlrtRSI;
    if (info > 2147483646) {
      c_st.site = &l_emlrtRSI;
      check_forloop_overflow_error(&c_st);
    }
    for (loop_ub = 0; loop_ub < info; loop_ub++) {
      U_data[loop_ub + U->size[0] * loop_ub].re = 1.0;
      U_data[loop_ub + U->size[0] * loop_ub].im = 0.0;
    }
    b_st.site = &xb_emlrtRSI;
    c_st.site = &yb_emlrtRSI;
    i = V->size[0] * V->size[1];
    V->size[0] = A->size[1];
    V->size[1] = A->size[1];
    emxEnsureCapacity_creal_T(&b_st, V, i, &qb_emlrtRTEI);
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
      for (info = 0; info < m; info++) {
        U_data[info + V->size[0] * info].re = 1.0;
        U_data[info + V->size[0] * info].im = 0.0;
      }
    }
    s->size[0] = 0;
  } else {
    ptrdiff_t info_t;
    int32_T m;
    st.site = &vb_emlrtRSI;
    emxInit_creal_T(&st, &b_A, 2, &xb_emlrtRTEI);
    i = b_A->size[0] * b_A->size[1];
    b_A->size[0] = A->size[0];
    b_A->size[1] = A->size[1];
    emxEnsureCapacity_creal_T(&st, b_A, i, &kb_emlrtRTEI);
    b_A_data = b_A->data;
    loop_ub = A->size[0] * A->size[1];
    for (i = 0; i < loop_ub; i++) {
      b_A_data[i] = A_data[i];
    }
    m = A->size[0];
    info = A->size[1];
    b_st.site = &dc_emlrtRSI;
    emxInit_creal_T(&b_st, &c_A, 2, &mb_emlrtRTEI);
    i = c_A->size[0] * c_A->size[1];
    c_A->size[0] = A->size[0];
    c_A->size[1] = A->size[1];
    emxEnsureCapacity_creal_T(&b_st, c_A, i, &mb_emlrtRTEI);
    c_A_data = c_A->data;
    loop_ub = A->size[0] * A->size[1];
    for (i = 0; i < loop_ub; i++) {
      c_A_data[i] = A_data[i];
    }
    i = U->size[0] * U->size[1];
    U->size[0] = A->size[0];
    U->size[1] = A->size[0];
    emxEnsureCapacity_creal_T(&b_st, U, i, &nb_emlrtRTEI);
    U_data = U->data;
    emxInit_creal_T(&b_st, &Vt, 2, &yb_emlrtRTEI);
    i = Vt->size[0] * Vt->size[1];
    Vt->size[0] = A->size[1];
    Vt->size[1] = A->size[1];
    emxEnsureCapacity_creal_T(&b_st, Vt, i, &ob_emlrtRTEI);
    Vt_data = Vt->data;
    i = s->size[0];
    s->size[0] = muIntScalarMin_sint32(info, m);
    emxEnsureCapacity_real_T(&b_st, s, i, &pb_emlrtRTEI);
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
      b_st.site = &cc_emlrtRSI;
      c_st.site = &fc_emlrtRSI;
      m = A->size[0];
      info = A->size[1];
      info = muIntScalarMin_sint32(info, m);
      i = U->size[0] * U->size[1];
      U->size[0] = A->size[0];
      U->size[1] = A->size[0];
      emxEnsureCapacity_creal_T(&c_st, U, i, &sb_emlrtRTEI);
      U_data = U->data;
      i = Vt->size[0] * Vt->size[1];
      Vt->size[0] = A->size[1];
      Vt->size[1] = A->size[1];
      emxEnsureCapacity_creal_T(&c_st, Vt, i, &tb_emlrtRTEI);
      Vt_data = Vt->data;
      i = s->size[0];
      s->size[0] = info;
      emxEnsureCapacity_real_T(&c_st, s, i, &ub_emlrtRTEI);
      s_data = s->data;
      emxInit_real_T(&c_st, &superb, &ac_emlrtRTEI);
      if (info > 1) {
        i = superb->size[0];
        superb->size[0] = info - 1;
        emxEnsureCapacity_real_T(&c_st, superb, i, &wb_emlrtRTEI);
        superb_data = superb->data;
      } else {
        i = superb->size[0];
        superb->size[0] = 1;
        emxEnsureCapacity_real_T(&c_st, superb, i, &vb_emlrtRTEI);
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
      emxEnsureCapacity_creal_T(&c_st, V, i, &rb_emlrtRTEI);
      U_data = V->data;
      loop_ub = Vt->size[0];
      for (i = 0; i < loop_ub; i++) {
        m = Vt->size[1];
        for (i1 = 0; i1 < m; i1++) {
          U_data[i1 + V->size[0] * i].re = Vt_data[i + Vt->size[0] * i1].re;
          U_data[i1 + V->size[0] * i].im = -Vt_data[i + Vt->size[0] * i1].im;
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
      i = V->size[0] * V->size[1];
      V->size[0] = Vt->size[1];
      V->size[1] = Vt->size[0];
      emxEnsureCapacity_creal_T(&st, V, i, &rb_emlrtRTEI);
      U_data = V->data;
      loop_ub = Vt->size[0];
      for (i = 0; i < loop_ub; i++) {
        m = Vt->size[1];
        for (i1 = 0; i1 < m; i1++) {
          U_data[i1 + V->size[0] * i].re = Vt_data[i + Vt->size[0] * i1].re;
          U_data[i1 + V->size[0] * i].im = -Vt_data[i + Vt->size[0] * i1].im;
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
