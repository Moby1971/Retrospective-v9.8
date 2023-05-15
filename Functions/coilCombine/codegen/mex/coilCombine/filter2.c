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
static emlrtRSInfo t_emlrtRSI =
    {
        55,        /* lineNo */
        "filter2", /* fcnName */
        "/Applications/MATLAB_R2022b.app/toolbox/eml/lib/matlab/datafun/"
        "filter2.m" /* pathName */
};

static emlrtRSInfo u_emlrtRSI =
    {
        61,        /* lineNo */
        "filter2", /* fcnName */
        "/Applications/MATLAB_R2022b.app/toolbox/eml/lib/matlab/datafun/"
        "filter2.m" /* pathName */
};

static emlrtRSInfo v_emlrtRSI =
    {
        62,        /* lineNo */
        "filter2", /* fcnName */
        "/Applications/MATLAB_R2022b.app/toolbox/eml/lib/matlab/datafun/"
        "filter2.m" /* pathName */
};

static emlrtRSInfo w_emlrtRSI =
    {
        69,        /* lineNo */
        "filter2", /* fcnName */
        "/Applications/MATLAB_R2022b.app/toolbox/eml/lib/matlab/datafun/"
        "filter2.m" /* pathName */
};

static emlrtRSInfo
    x_emlrtRSI =
        {
            42,      /* lineNo */
            "conv2", /* fcnName */
            "/Applications/MATLAB_R2022b.app/toolbox/eml/lib/matlab/datafun/"
            "conv2.m" /* pathName */
};

static emlrtRSInfo
    y_emlrtRSI =
        {
            265,              /* lineNo */
            "conv2Separable", /* fcnName */
            "/Applications/MATLAB_R2022b.app/toolbox/eml/lib/matlab/datafun/"
            "conv2.m" /* pathName */
};

static emlrtRSInfo
    ab_emlrtRSI =
        {
            234,              /* lineNo */
            "conv2Separable", /* fcnName */
            "/Applications/MATLAB_R2022b.app/toolbox/eml/lib/matlab/datafun/"
            "conv2.m" /* pathName */
};

static emlrtRSInfo
    bb_emlrtRSI =
        {
            229,              /* lineNo */
            "conv2Separable", /* fcnName */
            "/Applications/MATLAB_R2022b.app/toolbox/eml/lib/matlab/datafun/"
            "conv2.m" /* pathName */
};

static emlrtRSInfo
    cb_emlrtRSI =
        {
            228,              /* lineNo */
            "conv2Separable", /* fcnName */
            "/Applications/MATLAB_R2022b.app/toolbox/eml/lib/matlab/datafun/"
            "conv2.m" /* pathName */
};

static emlrtRSInfo
    db_emlrtRSI =
        {
            199,              /* lineNo */
            "conv2Separable", /* fcnName */
            "/Applications/MATLAB_R2022b.app/toolbox/eml/lib/matlab/datafun/"
            "conv2.m" /* pathName */
};

static emlrtRSInfo
    eb_emlrtRSI =
        {
            119,              /* lineNo */
            "conv2Separable", /* fcnName */
            "/Applications/MATLAB_R2022b.app/toolbox/eml/lib/matlab/datafun/"
            "conv2.m" /* pathName */
};

static emlrtRSInfo fb_emlrtRSI = {
    53,      /* lineNo */
    "xaxpy", /* fcnName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/+blas/"
    "xaxpy.m" /* pathName */
};

static emlrtRSInfo hb_emlrtRSI = {
    10,      /* lineNo */
    "round", /* fcnName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/lib/matlab/elfun/round.m" /* pathName
                                                                            */
};

static emlrtRSInfo ib_emlrtRSI = {
    33,                           /* lineNo */
    "applyScalarFunctionInPlace", /* fcnName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/"
    "applyScalarFunctionInPlace.m" /* pathName */
};

static emlrtRSInfo
    jb_emlrtRSI =
        {
            58,      /* lineNo */
            "conv2", /* fcnName */
            "/Applications/MATLAB_R2022b.app/toolbox/eml/lib/matlab/datafun/"
            "conv2.m" /* pathName */
};

static emlrtRSInfo
    kb_emlrtRSI =
        {
            75,                  /* lineNo */
            "conv2NonSeparable", /* fcnName */
            "/Applications/MATLAB_R2022b.app/toolbox/eml/lib/matlab/datafun/"
            "conv2.m" /* pathName */
};

static emlrtRSInfo lb_emlrtRSI =
    {
        105,     /* lineNo */
        "conv2", /* fcnName */
        "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/"
        "conv2.m" /* pathName */
};

static emlrtRSInfo mb_emlrtRSI = {
    32,          /* lineNo */
    "conv2AXPY", /* fcnName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/private/"
    "conv2AXPY.m" /* pathName */
};

static emlrtRTEInfo hb_emlrtRTEI =
    {
        55,        /* lineNo */
        9,         /* colNo */
        "filter2", /* fName */
        "/Applications/MATLAB_R2022b.app/toolbox/eml/lib/matlab/datafun/"
        "filter2.m" /* pName */
};

static emlrtRTEInfo
    ib_emlrtRTEI =
        {
            188,     /* lineNo */
            1,       /* colNo */
            "conv2", /* fName */
            "/Applications/MATLAB_R2022b.app/toolbox/eml/lib/matlab/datafun/"
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
  if (x->size[0] * x->size[1] >= 49) {
    int32_T ihi;
    int32_T mb;
    int32_T nb;
    boolean_T p;
    st.site = &t_emlrtRSI;
    b_st.site = &x_emlrtRSI;
    mb = x->size[0];
    nb = x->size[1];
    c_st.site = &eb_emlrtRSI;
    anyNonFinite(&c_st, x);
    i = y->size[0] * y->size[1];
    y->size[0] = x->size[0];
    y->size[1] = x->size[1];
    emxEnsureCapacity_creal_T(&b_st, y, i, &hb_emlrtRTEI);
    y_data = y->data;
    ihi = x->size[0] * x->size[1];
    for (i = 0; i < ihi; i++) {
      y_data[i].re = 0.0;
      y_data[i].im = 0.0;
    }
    emxInit_creal_T(&b_st, &work, 2, &ib_emlrtRTEI);
    i = work->size[0] * work->size[1];
    work->size[0] = x->size[0];
    work->size[1] = x->size[1];
    emxEnsureCapacity_creal_T(&b_st, work, i, &ib_emlrtRTEI);
    work_data = work->data;
    ihi = x->size[0] * x->size[1];
    for (i = 0; i < ihi; i++) {
      work_data[i].re = 0.0;
      work_data[i].im = 0.0;
    }
    if ((x->size[0] != 0) && (x->size[1] != 0)) {
      creal_T a;
      int32_T i1;
      int32_T ilo;
      i = 5 - x->size[0];
      i = muIntScalarMax_sint32(1, i);
      i1 = x->size[0] + 3;
      i1 = muIntScalarMin_sint32(7, i1);
      c_st.site = &db_emlrtRSI;
      for (k = i; k <= i1; k++) {
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
        c_st.site = &cb_emlrtRSI;
        if (nb > 2147483646) {
          d_st.site = &l_emlrtRSI;
          check_forloop_overflow_error(&d_st);
        }
        for (j = 0; j < nb; j++) {
          c_st.site = &bb_emlrtRSI;
          if (ihi >= 1) {
            d_st.site = &fb_emlrtRSI;
            a.re = hcol[k - 1];
            a.im = 0.0;
            n_t = (ptrdiff_t)ihi;
            incx_t = (ptrdiff_t)1;
            incy_t = (ptrdiff_t)1;
            zaxpy(&n_t, (real_T *)&a, (real_T *)&x_data[(j * mb + ilo) - k],
                  &incx_t, (real_T *)&work_data[(j * mb + ilo) - 4], &incy_t);
          }
        }
      }
      i = 5 - x->size[1];
      i = muIntScalarMax_sint32(1, i);
      i1 = x->size[1] + 3;
      i1 = muIntScalarMin_sint32(7, i1);
      c_st.site = &ab_emlrtRSI;
      for (k = i; k <= i1; k++) {
        if (k - 3 > 0) {
          ilo = k;
        } else {
          ilo = 4;
        }
        ihi = (nb + k) - 4;
        if (ihi > nb) {
          ihi = nb;
        }
        ihi = mb * ((ihi - ilo) + 4);
        c_st.site = &y_emlrtRSI;
        if (ihi >= 1) {
          d_st.site = &fb_emlrtRSI;
          a.re = hrow[k - 1];
          a.im = 0.0;
          n_t = (ptrdiff_t)ihi;
          incx_t = (ptrdiff_t)1;
          incy_t = (ptrdiff_t)1;
          zaxpy(&n_t, (real_T *)&a, (real_T *)&work_data[(ilo - k) * mb],
                &incx_t, (real_T *)&y_data[(ilo - 4) * mb], &incy_t);
        }
      }
    }
    emxFree_creal_T(&b_st, &work);
    st.site = &u_emlrtRSI;
    b_st.site = &r_emlrtRSI;
    ihi = x->size[0] * x->size[1];
    p = true;
    c_st.site = &s_emlrtRSI;
    if (ihi > 2147483646) {
      d_st.site = &l_emlrtRSI;
      check_forloop_overflow_error(&d_st);
    }
    k = 0;
    while ((k <= ihi - 1) && p) {
      p = ((muDoubleScalarFloor(x_data[k].re) == x_data[k].re) &&
           (muDoubleScalarFloor(x_data[k].im) == x_data[k].im));
      k++;
    }
    if (p) {
      st.site = &v_emlrtRSI;
      b_st.site = &hb_emlrtRSI;
      ihi = y->size[0] * y->size[1];
      c_st.site = &ib_emlrtRSI;
      if (ihi > 2147483646) {
        d_st.site = &l_emlrtRSI;
        check_forloop_overflow_error(&d_st);
      }
      for (k = 0; k < ihi; k++) {
        real_T d;
        d = y_data[k].im;
        y_data[k].re = muDoubleScalarRound(y_data[k].re);
        y_data[k].im = muDoubleScalarRound(d);
      }
    }
  } else {
    st.site = &w_emlrtRSI;
    b_st.site = &jb_emlrtRSI;
    c_st.site = &kb_emlrtRSI;
    d_st.site = &lb_emlrtRSI;
    e_st.site = &mb_emlrtRSI;
    conv2AXPYSameCMP(&e_st, x, y);
  }
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

/* End of code generation (filter2.c) */
