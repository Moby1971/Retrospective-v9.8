/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * conv2AXPYSameCMP.c
 *
 * Code generation for function 'conv2AXPYSameCMP'
 *
 */

/* Include files */
#include "conv2AXPYSameCMP.h"
#include "coilCombine_data.h"
#include "coilCombine_emxutil.h"
#include "coilCombine_types.h"
#include "eml_int_forloop_overflow_check.h"
#include "rt_nonfinite.h"
#include "blas.h"
#include "omp.h"
#include <stddef.h>

/* Variable Definitions */
static emlrtRSInfo nb_emlrtRSI = {
    63,                 /* lineNo */
    "conv2AXPYSameCMP", /* fcnName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/private/"
    "conv2AXPYSameCMP.m" /* pathName */
};

static emlrtRSInfo ob_emlrtRSI = {
    48,                 /* lineNo */
    "conv2AXPYSameCMP", /* fcnName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/private/"
    "conv2AXPYSameCMP.m" /* pathName */
};

static emlrtRSInfo pb_emlrtRSI = {
    36,                 /* lineNo */
    "conv2AXPYSameCMP", /* fcnName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/private/"
    "conv2AXPYSameCMP.m" /* pathName */
};

static emlrtRTEInfo cc_emlrtRTEI = {
    26,                 /* lineNo */
    5,                  /* colNo */
    "conv2AXPYSameCMP", /* fName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/private/"
    "conv2AXPYSameCMP.m" /* pName */
};

static emlrtRTEInfo dc_emlrtRTEI = {
    34,                 /* lineNo */
    20,                 /* colNo */
    "conv2AXPYSameCMP", /* fName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/private/"
    "conv2AXPYSameCMP.m" /* pName */
};

static emlrtRTEInfo ec_emlrtRTEI = {
    37,                 /* lineNo */
    5,                  /* colNo */
    "conv2AXPYSameCMP", /* fName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/private/"
    "conv2AXPYSameCMP.m" /* pName */
};

/* Function Definitions */
void conv2AXPYSameCMP(const emlrtStack *sp, const emxArray_creal_T *a,
                      emxArray_creal_T *c)
{
  jmp_buf emlrtJBEnviron;
  ptrdiff_t incx_t;
  ptrdiff_t incy_t;
  ptrdiff_t n_t;
  jmp_buf *volatile emlrtJBStack;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack st;
  emxArray_creal_T *cj;
  const creal_T *a_data;
  creal_T *c_data;
  creal_T *cj_data;
  int32_T conv2AXPYSameCMP_numThreads;
  int32_T i;
  int32_T ia0;
  int32_T ib;
  int32_T ic0;
  int32_T j;
  int32_T jb;
  int32_T jbmax;
  int32_T jbmin;
  int32_T jcut;
  int32_T ma;
  int32_T n;
  int32_T na;
  int32_T ub_loop;
  boolean_T emlrtHadParallelError = false;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  a_data = a->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  ma = a->size[0];
  na = a->size[1];
  if ((a->size[0] == 0) || (a->size[1] == 0)) {
    i = c->size[0] * c->size[1];
    c->size[0] = a->size[0];
    c->size[1] = a->size[1];
    emxEnsureCapacity_creal_T(sp, c, i, &cc_emlrtRTEI);
    c_data = c->data;
    ub_loop = a->size[0] * a->size[1];
    for (i = 0; i < ub_loop; i++) {
      c_data[i].re = 0.0;
      c_data[i].im = 0.0;
    }
  } else {
    jcut = a->size[1] - 3;
    i = c->size[0] * c->size[1];
    c->size[0] = a->size[0];
    c->size[1] = a->size[1];
    emxEnsureCapacity_creal_T(sp, c, i, &dc_emlrtRTEI);
    c_data = c->data;
    st.site = &pb_emlrtRSI;
    if (a->size[1] > 2147483646) {
      b_st.site = &l_emlrtRSI;
      check_forloop_overflow_error(&b_st);
    }
    ub_loop = a->size[1] - 1;
    emlrtEnterParallelRegion((emlrtCTX)sp, omp_in_parallel());
    emlrtPushJmpBuf((emlrtCTX)sp, &emlrtJBStack);
    conv2AXPYSameCMP_numThreads = emlrtAllocRegionTLSs(
        sp->tls, omp_in_parallel(), omp_get_max_threads(), omp_get_num_procs());
#pragma omp parallel num_threads(conv2AXPYSameCMP_numThreads) private(         \
    cj_data, cj, ia0, n, ic0, jbmax, jbmin, emlrtJBEnviron, e_st, jb, ib, n_t, \
    incx_t, incy_t) firstprivate(c_st, d_st, emlrtHadParallelError)
    {
      if (setjmp(emlrtJBEnviron) == 0) {
        c_st.prev = sp;
        c_st.tls = emlrtAllocTLS((emlrtCTX)sp, omp_get_thread_num());
        c_st.site = NULL;
        emlrtSetJmpBuf(&c_st, &emlrtJBEnviron);
        d_st.prev = &c_st;
        d_st.tls = c_st.tls;
        e_st.prev = &d_st;
        e_st.tls = d_st.tls;
        emxInit_creal_T(&c_st, &cj, 1, &ec_emlrtRTEI);
      } else {
        emlrtHadParallelError = true;
      }
#pragma omp for nowait
      for (j = 0; j <= ub_loop; j++) {
        if (emlrtHadParallelError) {
          continue;
        }
        if (setjmp(emlrtJBEnviron) == 0) {
          n = cj->size[0];
          cj->size[0] = ma;
          emxEnsureCapacity_creal_T(&c_st, cj, n, &ec_emlrtRTEI);
          cj_data = cj->data;
          for (n = 0; n < ma; n++) {
            cj_data[n].re = 0.0;
            cj_data[n].im = 0.0;
          }
          if (j + 1 > 3) {
            jbmin = 0;
          } else {
            jbmin = 3 - j;
          }
          if (j + 1 < jcut) {
            jbmax = 6;
          } else {
            jbmax = (na - j) + 2;
          }
          d_st.site = &ob_emlrtRSI;
          if ((jbmin <= jbmax) && (jbmax > 2147483646)) {
            e_st.site = &l_emlrtRSI;
            check_forloop_overflow_error(&e_st);
          }
          for (jb = jbmin; jb <= jbmax; jb++) {
            for (ib = 0; ib < 7; ib++) {
              if (ib < 3) {
                ic0 = 4 - ib;
                n = (ma + ib) - 3;
              } else {
                ic0 = 1;
                n = (ma - ib) + 3;
              }
              ia0 = ((ic0 + ib) + ((jb + j) - 3) * ma) - 3;
              d_st.site = &nb_emlrtRSI;
              if (n >= 1) {
                n_t = (ptrdiff_t)n;
                incx_t = (ptrdiff_t)1;
                incy_t = (ptrdiff_t)1;
                zaxpy(&n_t, (real_T *)&dc, (real_T *)&a_data[ia0 - 1], &incx_t,
                      (real_T *)&cj_data[ic0 - 1], &incy_t);
              }
            }
          }
          ia0 = cj->size[0];
          for (n = 0; n < ia0; n++) {
            c_data[n + c->size[0] * j] = cj_data[n];
          }
        } else {
          emlrtHadParallelError = true;
        }
      }
      if (!emlrtHadParallelError) {
        emlrtHeapReferenceStackLeaveScope(&c_st, 1);
        emxFree_creal_T(&c_st, &cj);
      }
    }
    emlrtPopJmpBuf((emlrtCTX)sp, &emlrtJBStack);
    emlrtExitParallelRegion((emlrtCTX)sp, omp_in_parallel());
  }
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

/* End of code generation (conv2AXPYSameCMP.c) */
