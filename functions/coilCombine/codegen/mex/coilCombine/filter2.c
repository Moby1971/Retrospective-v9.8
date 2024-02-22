/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * filter2.c
 *
 * Code generation for function 'filter2'
 *
 */

/* Include files */
#include "filter2.h"
#include "anyNonFinite.h"
#include "coilCombine_data.h"
#include "coilCombine_emxutil.h"
#include "coilCombine_types.h"
#include "conv2AXPYSameCMP.h"
#include "eml_int_forloop_overflow_check.h"
#include "rt_nonfinite.h"
#include "blas.h"
#include "mwmathutil.h"
#include <stddef.h>

/* Variable Definitions */
static emlrtRSInfo q_emlrtRSI =
    {
        55,        /* lineNo */
        "filter2", /* fcnName */
        "/Applications/MATLAB_R2023b.app/toolbox/eml/lib/matlab/datafun/"
        "filter2.m" /* pathName */
};

static emlrtRSInfo r_emlrtRSI =
    {
        61,        /* lineNo */
        "filter2", /* fcnName */
        "/Applications/MATLAB_R2023b.app/toolbox/eml/lib/matlab/datafun/"
        "filter2.m" /* pathName */
};

static emlrtRSInfo s_emlrtRSI =
    {
        62,        /* lineNo */
        "filter2", /* fcnName */
        "/Applications/MATLAB_R2023b.app/toolbox/eml/lib/matlab/datafun/"
        "filter2.m" /* pathName */
};

static emlrtRSInfo t_emlrtRSI =
    {
        69,        /* lineNo */
        "filter2", /* fcnName */
        "/Applications/MATLAB_R2023b.app/toolbox/eml/lib/matlab/datafun/"
        "filter2.m" /* pathName */
};

static emlrtRSInfo
    u_emlrtRSI =
        {
            42,      /* lineNo */
            "conv2", /* fcnName */
            "/Applications/MATLAB_R2023b.app/toolbox/eml/lib/matlab/datafun/"
            "conv2.m" /* pathName */
};

static emlrtRSInfo
    v_emlrtRSI =
        {
            265,              /* lineNo */
            "conv2Separable", /* fcnName */
            "/Applications/MATLAB_R2023b.app/toolbox/eml/lib/matlab/datafun/"
            "conv2.m" /* pathName */
};

static emlrtRSInfo
    w_emlrtRSI =
        {
            234,              /* lineNo */
            "conv2Separable", /* fcnName */
            "/Applications/MATLAB_R2023b.app/toolbox/eml/lib/matlab/datafun/"
            "conv2.m" /* pathName */
};

static emlrtRSInfo
    x_emlrtRSI =
        {
            229,              /* lineNo */
            "conv2Separable", /* fcnName */
            "/Applications/MATLAB_R2023b.app/toolbox/eml/lib/matlab/datafun/"
            "conv2.m" /* pathName */
};

static emlrtRSInfo
    y_emlrtRSI =
        {
            228,              /* lineNo */
            "conv2Separable", /* fcnName */
            "/Applications/MATLAB_R2023b.app/toolbox/eml/lib/matlab/datafun/"
            "conv2.m" /* pathName */
};

static emlrtRSInfo
    ab_emlrtRSI =
        {
            199,              /* lineNo */
            "conv2Separable", /* fcnName */
            "/Applications/MATLAB_R2023b.app/toolbox/eml/lib/matlab/datafun/"
            "conv2.m" /* pathName */
};

static emlrtRSInfo
    bb_emlrtRSI =
        {
            119,              /* lineNo */
            "conv2Separable", /* fcnName */
            "/Applications/MATLAB_R2023b.app/toolbox/eml/lib/matlab/datafun/"
            "conv2.m" /* pathName */
};

static emlrtRSInfo fb_emlrtRSI = {
    53,      /* lineNo */
    "xaxpy", /* fcnName */
    "/Applications/MATLAB_R2023b.app/toolbox/eml/eml/+coder/+internal/+blas/"
    "xaxpy.m" /* pathName */
};

static emlrtRSInfo hb_emlrtRSI = {
    16,      /* lineNo */
    "round", /* fcnName */
    "/Applications/MATLAB_R2023b.app/toolbox/eml/lib/matlab/elfun/round.m" /* pathName
                                                                            */
};

static emlrtRSInfo ib_emlrtRSI = {
    33,                           /* lineNo */
    "applyScalarFunctionInPlace", /* fcnName */
    "/Applications/MATLAB_R2023b.app/toolbox/eml/eml/+coder/+internal/"
    "applyScalarFunctionInPlace.m" /* pathName */
};

static emlrtRSInfo
    jb_emlrtRSI =
        {
            58,      /* lineNo */
            "conv2", /* fcnName */
            "/Applications/MATLAB_R2023b.app/toolbox/eml/lib/matlab/datafun/"
            "conv2.m" /* pathName */
};

static emlrtRSInfo
    kb_emlrtRSI =
        {
            75,                  /* lineNo */
            "conv2NonSeparable", /* fcnName */
            "/Applications/MATLAB_R2023b.app/toolbox/eml/lib/matlab/datafun/"
            "conv2.m" /* pathName */
};

static emlrtRSInfo lb_emlrtRSI =
    {
        105,     /* lineNo */
        "conv2", /* fcnName */
        "/Applications/MATLAB_R2023b.app/toolbox/eml/eml/+coder/+internal/"
        "conv2.m" /* pathName */
};

static emlrtRSInfo mb_emlrtRSI = {
    32,          /* lineNo */
    "conv2AXPY", /* fcnName */
    "/Applications/MATLAB_R2023b.app/toolbox/eml/eml/+coder/+internal/private/"
    "conv2AXPY.m" /* pathName */
};

static emlrtRTEInfo jb_emlrtRTEI =
    {
        55,        /* lineNo */
        9,         /* colNo */
        "filter2", /* fName */
        "/Applications/MATLAB_R2023b.app/toolbox/eml/lib/matlab/datafun/"
        "filter2.m" /* pName */
};

static emlrtRTEInfo
    kb_emlrtRTEI =
        {
            188,     /* lineNo */
            1,       /* colNo */
            "conv2", /* fName */
            "/Applications/MATLAB_R2023b.app/toolbox/eml/lib/matlab/datafun/"
            "conv2.m" /* pName */
};

/* Function Definitions */
void filter2(const emlrtStack *sp, const emxArray_creal_T *x,
             emxArray_creal_T *y)
{
  static const real_T hcol[7] = {
      1.0000000000000002, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  static const real_T hrow[7] = {1.0,
                                 1.0000000000000004,
                                 1.0000000000000002,
                                 1.0000000000000002,
                                 1.0000000000000002,
                                 1.0000000000000002,
                                 1.0000000000000002};
  ptrdiff_t incx_t;
  ptrdiff_t incy_t;
  ptrdiff_t n_t;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack st;
  emxArray_creal_T *work;
  const creal_T *x_data;
  creal_T *work_data;
  creal_T *y_data;
  int32_T i;
  int32_T j;
  int32_T k;
  int32_T loop_ub;
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
  x_data = x->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  loop_ub = x->size[0] * x->size[1];
  if (loop_ub >= 49) {
    real_T re_tmp;
    int32_T ihi;
    int32_T mb;
    int32_T nb;
    boolean_T p;
    st.site = &q_emlrtRSI;
    b_st.site = &u_emlrtRSI;
    mb = x->size[0];
    nb = x->size[1];
    c_st.site = &bb_emlrtRSI;
    anyNonFinite(&c_st, x);
    i = y->size[0] * y->size[1];
    y->size[0] = x->size[0];
    y->size[1] = x->size[1];
    emxEnsureCapacity_creal_T(&b_st, y, i, &jb_emlrtRTEI);
    y_data = y->data;
    for (i = 0; i < loop_ub; i++) {
      y_data[i].re = 0.0;
      y_data[i].im = 0.0;
    }
    emxInit_creal_T(&b_st, &work, 2, &kb_emlrtRTEI);
    i = work->size[0] * work->size[1];
    work->size[0] = x->size[0];
    work->size[1] = x->size[1];
    emxEnsureCapacity_creal_T(&b_st, work, i, &kb_emlrtRTEI);
    work_data = work->data;
    for (i = 0; i < loop_ub; i++) {
      work_data[i].re = 0.0;
      work_data[i].im = 0.0;
    }
    if ((x->size[0] != 0) && (x->size[1] != 0)) {
      creal_T a;
      int32_T i1;
      int32_T jlo;
      i = 5 - x->size[0];
      i = muIntScalarMax_sint32(1, i);
      i1 = x->size[0] + 3;
      i1 = muIntScalarMin_sint32(7, i1);
      c_st.site = &ab_emlrtRSI;
      for (k = i; k <= i1; k++) {
        int32_T ilo;
        if (k - 3 > 0) {
          ilo = k;
        } else {
          ilo = 4;
        }
        ihi = (mb + k) - 4;
        if (ihi > mb) {
          ihi = mb;
        }
        ihi = (ihi - ilo) + 4;
        c_st.site = &y_emlrtRSI;
        if (nb > 2147483646) {
          d_st.site = &l_emlrtRSI;
          check_forloop_overflow_error(&d_st);
        }
        for (j = 0; j < nb; j++) {
          c_st.site = &x_emlrtRSI;
          if (ihi >= 1) {
            d_st.site = &fb_emlrtRSI;
            a.re = hcol[k - 1];
            a.im = 0.0;
            n_t = (ptrdiff_t)ihi;
            incx_t = (ptrdiff_t)1;
            incy_t = (ptrdiff_t)1;
            jlo = j * mb + ilo;
            zaxpy(&n_t, (real_T *)&a, (real_T *)&x_data[jlo - k], &incx_t,
                  (real_T *)&work_data[jlo - 4], &incy_t);
          }
        }
      }
      i = 5 - x->size[1];
      i = muIntScalarMax_sint32(1, i);
      i1 = x->size[1] + 3;
      i1 = muIntScalarMin_sint32(7, i1);
      c_st.site = &w_emlrtRSI;
      for (k = i; k <= i1; k++) {
        if (k - 3 > 0) {
          jlo = k;
        } else {
          jlo = 4;
        }
        ihi = (nb + k) - 4;
        if (ihi > nb) {
          ihi = nb;
        }
        ihi = mb * ((ihi - jlo) + 4);
        c_st.site = &v_emlrtRSI;
        if (ihi >= 1) {
          d_st.site = &fb_emlrtRSI;
          a.re = hrow[k - 1];
          a.im = 0.0;
          n_t = (ptrdiff_t)ihi;
          incx_t = (ptrdiff_t)1;
          incy_t = (ptrdiff_t)1;
          zaxpy(&n_t, (real_T *)&a, (real_T *)&work_data[(jlo - k) * mb],
                &incx_t, (real_T *)&y_data[(jlo - 4) * mb], &incy_t);
        }
      }
    }
    emxFree_creal_T(&b_st, &work);
    st.site = &r_emlrtRSI;
    b_st.site = &db_emlrtRSI;
    p = true;
    c_st.site = &eb_emlrtRSI;
    if (loop_ub > 2147483646) {
      d_st.site = &l_emlrtRSI;
      check_forloop_overflow_error(&d_st);
    }
    k = 0;
    while ((k <= loop_ub - 1) && p) {
      real_T im_tmp;
      re_tmp = x_data[k].re;
      im_tmp = x_data[k].im;
      p = ((muDoubleScalarFloor(re_tmp) == re_tmp) &&
           (muDoubleScalarFloor(im_tmp) == im_tmp));
      k++;
    }
    if (p) {
      st.site = &s_emlrtRSI;
      b_st.site = &hb_emlrtRSI;
      ihi = y->size[0] * y->size[1];
      c_st.site = &ib_emlrtRSI;
      if (ihi > 2147483646) {
        d_st.site = &l_emlrtRSI;
        check_forloop_overflow_error(&d_st);
      }
      for (k = 0; k < ihi; k++) {
        re_tmp = y_data[k].im;
        y_data[k].re = muDoubleScalarRound(y_data[k].re);
        y_data[k].im = muDoubleScalarRound(re_tmp);
      }
    }
  } else {
    st.site = &t_emlrtRSI;
    b_st.site = &jb_emlrtRSI;
    c_st.site = &kb_emlrtRSI;
    d_st.site = &lb_emlrtRSI;
    e_st.site = &mb_emlrtRSI;
    conv2AXPYSameCMP(&e_st, x, y);
  }
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

/* End of code generation (filter2.c) */
