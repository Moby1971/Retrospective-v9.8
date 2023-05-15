/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * coilCombine.c
 *
 * Code generation for function 'coilCombine'
 *
 */

/* Include files */
#include "coilCombine.h"
#include "anyNonFinite.h"
#include "coilCombine_data.h"
#include "coilCombine_emxutil.h"
#include "coilCombine_types.h"
#include "filter2.h"
#include "permute.h"
#include "rt_nonfinite.h"
#include "squeeze.h"
#include "svd.h"
#include "blas.h"
#include "mwmathutil.h"
#include "omp.h"
#include <stddef.h>

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = {
    6,             /* lineNo */
    "coilCombine", /* fcnName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m" /* pathName */
};

static emlrtRSInfo b_emlrtRSI = {
    12,            /* lineNo */
    "coilCombine", /* fcnName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m" /* pathName */
};

static emlrtRSInfo c_emlrtRSI = {
    15,            /* lineNo */
    "coilCombine", /* fcnName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m" /* pathName */
};

static emlrtRSInfo
    d_emlrtRSI =
        {
            7,         /* lineNo */
            "ref/ref", /* fcnName */
            "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/"
            "ref.m" /* pathName */
};

static emlrtRSInfo m_emlrtRSI = {
    20,                        /* lineNo */
    "coilCombine/coilCombine", /* fcnName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m" /* pathName */
};

static emlrtRSInfo n_emlrtRSI = {
    34,                        /* lineNo */
    "coilCombine/coilCombine", /* fcnName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m" /* pathName */
};

static emlrtRSInfo o_emlrtRSI = {
    42,                        /* lineNo */
    "coilCombine/coilCombine", /* fcnName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m" /* pathName */
};

static emlrtRSInfo p_emlrtRSI = {
    44,                        /* lineNo */
    "coilCombine/coilCombine", /* fcnName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m" /* pathName */
};

static emlrtRSInfo rb_emlrtRSI = {
    14,    /* lineNo */
    "svd", /* fcnName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/lib/matlab/matfun/svd.m" /* pathName
                                                                           */
};

static emlrtRSInfo sb_emlrtRSI = {
    36,    /* lineNo */
    "svd", /* fcnName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/lib/matlab/matfun/svd.m" /* pathName
                                                                           */
};

static emlrtRSInfo tb_emlrtRSI = {
    42,    /* lineNo */
    "svd", /* fcnName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/lib/matlab/matfun/svd.m" /* pathName
                                                                           */
};

static emlrtRSInfo hc_emlrtRSI = {
    40,                  /* lineNo */
    "reshapeSizeChecks", /* fcnName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/"
    "reshapeSizeChecks.m" /* pathName */
};

static emlrtRSInfo ic_emlrtRSI = {
    94,                  /* lineNo */
    "eml_mtimes_helper", /* fcnName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/lib/matlab/ops/"
    "eml_mtimes_helper.m" /* pathName */
};

static emlrtRSInfo jc_emlrtRSI = {
    69,                  /* lineNo */
    "eml_mtimes_helper", /* fcnName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/lib/matlab/ops/"
    "eml_mtimes_helper.m" /* pathName */
};

static emlrtRSInfo kc_emlrtRSI = {
    142,      /* lineNo */
    "mtimes", /* fcnName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/+blas/"
    "mtimes.m" /* pathName */
};

static emlrtRSInfo lc_emlrtRSI = {
    178,           /* lineNo */
    "mtimes_blas", /* fcnName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/+blas/"
    "mtimes.m" /* pathName */
};

static emlrtRTEInfo emlrtRTEI = {
    133,                   /* lineNo */
    23,                    /* colNo */
    "dynamic_size_checks", /* fName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/lib/matlab/ops/"
    "eml_mtimes_helper.m" /* pName */
};

static emlrtRTEInfo b_emlrtRTEI = {
    138,                   /* lineNo */
    23,                    /* colNo */
    "dynamic_size_checks", /* fName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/lib/matlab/ops/"
    "eml_mtimes_helper.m" /* pName */
};

static emlrtRTEInfo c_emlrtRTEI = {
    64,                   /* lineNo */
    15,                   /* colNo */
    "assertValidSizeArg", /* fName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/"
    "assertValidSizeArg.m" /* pName */
};

static emlrtECInfo emlrtECI = {
    -1,                        /* nDims */
    44,                        /* lineNo */
    17,                        /* colNo */
    "coilCombine/coilCombine", /* fName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m" /* pName */
};

static emlrtBCInfo emlrtBCI = {
    -1,                        /* iFirst */
    -1,                        /* iLast */
    44,                        /* lineNo */
    24,                        /* colNo */
    "im2",                     /* aName */
    "coilCombine/coilCombine", /* fName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m", /* pName */
    0          /* checkKind */
};

static emlrtBCInfo b_emlrtBCI = {
    -1,                        /* iFirst */
    -1,                        /* iLast */
    44,                        /* lineNo */
    21,                        /* colNo */
    "im2",                     /* aName */
    "coilCombine/coilCombine", /* fName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m", /* pName */
    0          /* checkKind */
};

static emlrtBCInfo c_emlrtBCI = {
    -1,                        /* iFirst */
    -1,                        /* iLast */
    44,                        /* lineNo */
    67,                        /* colNo */
    "im1",                     /* aName */
    "coilCombine/coilCombine", /* fName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m", /* pName */
    0          /* checkKind */
};

static emlrtBCInfo d_emlrtBCI = {
    -1,                        /* iFirst */
    -1,                        /* iLast */
    44,                        /* lineNo */
    64,                        /* colNo */
    "im1",                     /* aName */
    "coilCombine/coilCombine", /* fName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m", /* pName */
    0          /* checkKind */
};

static emlrtBCInfo e_emlrtBCI = {
    -1,                        /* iFirst */
    -1,                        /* iLast */
    43,                        /* lineNo */
    30,                        /* colNo */
    "U",                       /* aName */
    "coilCombine/coilCombine", /* fName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m", /* pName */
    0          /* checkKind */
};

static emlrtBCInfo f_emlrtBCI = {
    -1,                        /* iFirst */
    -1,                        /* iLast */
    42,                        /* lineNo */
    43,                        /* colNo */
    "Rs",                      /* aName */
    "coilCombine/coilCombine", /* fName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m", /* pName */
    0          /* checkKind */
};

static emlrtBCInfo g_emlrtBCI = {
    -1,                        /* iFirst */
    -1,                        /* iLast */
    42,                        /* lineNo */
    40,                        /* colNo */
    "Rs",                      /* aName */
    "coilCombine/coilCombine", /* fName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m", /* pName */
    0          /* checkKind */
};

static emlrtECInfo b_emlrtECI = {
    -1,                        /* nDims */
    34,                        /* lineNo */
    21,                        /* colNo */
    "coilCombine/coilCombine", /* fName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m" /* pName */
};

static emlrtBCInfo h_emlrtBCI = {
    -1,                        /* iFirst */
    -1,                        /* iLast */
    34,                        /* lineNo */
    32,                        /* colNo */
    "Rs",                      /* aName */
    "coilCombine/coilCombine", /* fName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m", /* pName */
    0          /* checkKind */
};

static emlrtBCInfo i_emlrtBCI = {
    -1,                        /* iFirst */
    -1,                        /* iLast */
    34,                        /* lineNo */
    28,                        /* colNo */
    "Rs",                      /* aName */
    "coilCombine/coilCombine", /* fName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m", /* pName */
    0          /* checkKind */
};

static emlrtECInfo c_emlrtECI = {
    2,                         /* nDims */
    34,                        /* lineNo */
    39,                        /* colNo */
    "coilCombine/coilCombine", /* fName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m" /* pName */
};

static emlrtECInfo d_emlrtECI = {
    1,                         /* nDims */
    34,                        /* lineNo */
    39,                        /* colNo */
    "coilCombine/coilCombine", /* fName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m" /* pName */
};

static emlrtECInfo e_emlrtECI = {
    2,                         /* nDims */
    34,                        /* lineNo */
    80,                        /* colNo */
    "coilCombine/coilCombine", /* fName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m" /* pName */
};

static emlrtBCInfo j_emlrtBCI = {
    -1,                        /* iFirst */
    -1,                        /* iLast */
    34,                        /* lineNo */
    113,                       /* colNo */
    "im1",                     /* aName */
    "coilCombine/coilCombine", /* fName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m", /* pName */
    0          /* checkKind */
};

static emlrtBCInfo k_emlrtBCI = {
    -1,                        /* iFirst */
    -1,                        /* iLast */
    34,                        /* lineNo */
    110,                       /* colNo */
    "im1",                     /* aName */
    "coilCombine/coilCombine", /* fName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m", /* pName */
    0          /* checkKind */
};

static emlrtBCInfo l_emlrtBCI = {
    -1,                        /* iFirst */
    -1,                        /* iLast */
    34,                        /* lineNo */
    91,                        /* colNo */
    "im1",                     /* aName */
    "coilCombine/coilCombine", /* fName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m", /* pName */
    0          /* checkKind */
};

static emlrtBCInfo m_emlrtBCI = {
    -1,                        /* iFirst */
    -1,                        /* iLast */
    34,                        /* lineNo */
    88,                        /* colNo */
    "im1",                     /* aName */
    "coilCombine/coilCombine", /* fName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m", /* pName */
    0          /* checkKind */
};

static emlrtBCInfo n_emlrtBCI = {
    -1,                        /* iFirst */
    -1,                        /* iLast */
    34,                        /* lineNo */
    50,                        /* colNo */
    "Rs",                      /* aName */
    "coilCombine/coilCombine", /* fName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m", /* pName */
    0          /* checkKind */
};

static emlrtBCInfo o_emlrtBCI = {
    -1,                        /* iFirst */
    -1,                        /* iLast */
    34,                        /* lineNo */
    46,                        /* colNo */
    "Rs",                      /* aName */
    "coilCombine/coilCombine", /* fName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m", /* pName */
    0          /* checkKind */
};

static emlrtECInfo f_emlrtECI = {
    -1,            /* nDims */
    12,            /* lineNo */
    5,             /* colNo */
    "coilCombine", /* fName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m" /* pName */
};

static emlrtBCInfo p_emlrtBCI = {
    -1,            /* iFirst */
    -1,            /* iLast */
    12,            /* lineNo */
    13,            /* colNo */
    "im2",         /* aName */
    "coilCombine", /* fName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m", /* pName */
    0          /* checkKind */
};

static emlrtBCInfo q_emlrtBCI = {
    -1,            /* iFirst */
    -1,            /* iLast */
    12,            /* lineNo */
    43,            /* colNo */
    "im1",         /* aName */
    "coilCombine", /* fName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m", /* pName */
    0          /* checkKind */
};

static emlrtRTEInfo j_emlrtRTEI = {
    6,             /* lineNo */
    15,            /* colNo */
    "coilCombine", /* fName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m" /* pName */
};

static emlrtRTEInfo k_emlrtRTEI = {
    10,            /* lineNo */
    1,             /* colNo */
    "coilCombine", /* fName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m" /* pName */
};

static emlrtRTEInfo l_emlrtRTEI = {
    12,            /* lineNo */
    35,            /* colNo */
    "coilCombine", /* fName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m" /* pName */
};

static emlrtRTEInfo m_emlrtRTEI = {
    12,            /* lineNo */
    23,            /* colNo */
    "coilCombine", /* fName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m" /* pName */
};

static emlrtRTEInfo n_emlrtRTEI = {
    28,            /* lineNo */
    9,             /* colNo */
    "coilCombine", /* fName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m" /* pName */
};

static emlrtRTEInfo o_emlrtRTEI = {
    34,            /* lineNo */
    97,            /* colNo */
    "coilCombine", /* fName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m" /* pName */
};

static emlrtRTEInfo p_emlrtRTEI = {
    42,            /* lineNo */
    37,            /* colNo */
    "coilCombine", /* fName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m" /* pName */
};

static emlrtRTEInfo q_emlrtRTEI = {
    41,    /* lineNo */
    14,    /* colNo */
    "svd", /* fName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/lib/matlab/matfun/svd.m" /* pName
                                                                           */
};

static emlrtRTEInfo r_emlrtRTEI = {
    43,            /* lineNo */
    17,            /* colNo */
    "coilCombine", /* fName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m" /* pName */
};

static emlrtRTEInfo s_emlrtRTEI = {
    34,            /* lineNo */
    80,            /* colNo */
    "coilCombine", /* fName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m" /* pName */
};

static emlrtRTEInfo t_emlrtRTEI = {
    43,    /* lineNo */
    9,     /* colNo */
    "svd", /* fName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/lib/matlab/matfun/svd.m" /* pName
                                                                           */
};

static emlrtRTEInfo u_emlrtRTEI = {
    44,            /* lineNo */
    60,            /* colNo */
    "coilCombine", /* fName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m" /* pName */
};

static emlrtRTEInfo v_emlrtRTEI = {
    34,            /* lineNo */
    39,            /* colNo */
    "coilCombine", /* fName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m" /* pName */
};

static emlrtRTEInfo w_emlrtRTEI = {
    44,            /* lineNo */
    52,            /* colNo */
    "coilCombine", /* fName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m" /* pName */
};

static emlrtRTEInfo x_emlrtRTEI = {
    44,            /* lineNo */
    44,            /* colNo */
    "coilCombine", /* fName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m" /* pName */
};

static emlrtRTEInfo y_emlrtRTEI = {
    44,            /* lineNo */
    36,            /* colNo */
    "coilCombine", /* fName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m" /* pName */
};

static emlrtRTEInfo ab_emlrtRTEI = {
    218,      /* lineNo */
    20,       /* colNo */
    "mtimes", /* fName */
    "/Applications/MATLAB_R2022b.app/toolbox/eml/eml/+coder/+internal/+blas/"
    "mtimes.m" /* pName */
};

static emlrtRTEInfo bb_emlrtRTEI = {
    4,             /* lineNo */
    10,            /* colNo */
    "coilCombine", /* fName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m" /* pName */
};

static emlrtRTEInfo cb_emlrtRTEI = {
    4,             /* lineNo */
    16,            /* colNo */
    "coilCombine", /* fName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m" /* pName */
};

static emlrtRTEInfo db_emlrtRTEI = {
    18,            /* lineNo */
    32,            /* colNo */
    "coilCombine", /* fName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m" /* pName */
};

static emlrtRTEInfo eb_emlrtRTEI = {
    42,            /* lineNo */
    29,            /* colNo */
    "coilCombine", /* fName */
    "/Users/gustav/Library/CloudStorage/Dropbox/"
    "Reconstruction_and_analysis_tools/Matlab/MRI-apps/AI4MRI/Functions/"
    "coilCombine/coilCo"
    "mbine.m" /* pName */
};

/* Function Declarations */
static void binary_expand_op(const emlrtStack *sp, emxArray_creal_T *in1,
                             const emxArray_creal_T *in2, int32_T in3,
                             int32_T in4);

/* Function Definitions */
static void binary_expand_op(const emlrtStack *sp, emxArray_creal_T *in1,
                             const emxArray_creal_T *in2, int32_T in3,
                             int32_T in4)
{
  emxArray_creal_T *b_in2;
  const creal_T *in2_data;
  creal_T *b_in2_data;
  creal_T *in1_data;
  int32_T aux_0_1;
  int32_T aux_1_1;
  int32_T b_loop_ub;
  int32_T i;
  int32_T i1;
  int32_T loop_ub;
  int32_T stride_0_0;
  int32_T stride_0_1;
  int32_T stride_1_0;
  int32_T stride_1_1;
  in2_data = in2->data;
  in1_data = in1->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  emxInit_creal_T(sp, &b_in2, 2, &v_emlrtRTEI);
  i = b_in2->size[0] * b_in2->size[1];
  if (in1->size[0] == 1) {
    b_in2->size[0] = in2->size[0];
  } else {
    b_in2->size[0] = in1->size[0];
  }
  if (in1->size[1] == 1) {
    b_in2->size[1] = in2->size[1];
  } else {
    b_in2->size[1] = in1->size[1];
  }
  emxEnsureCapacity_creal_T(sp, b_in2, i, &v_emlrtRTEI);
  b_in2_data = b_in2->data;
  stride_0_0 = (in2->size[0] != 1);
  stride_0_1 = (in2->size[1] != 1);
  stride_1_0 = (in1->size[0] != 1);
  stride_1_1 = (in1->size[1] != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  if (in1->size[1] == 1) {
    loop_ub = in2->size[1];
  } else {
    loop_ub = in1->size[1];
  }
  for (i = 0; i < loop_ub; i++) {
    i1 = in1->size[0];
    if (i1 == 1) {
      b_loop_ub = in2->size[0];
    } else {
      b_loop_ub = i1;
    }
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      int32_T i2;
      i2 = i1 * stride_1_0;
      b_in2_data[i1 + b_in2->size[0] * i].re =
          in2_data[((i1 * stride_0_0 + in2->size[0] * aux_0_1) +
                    in2->size[0] * in2->size[1] * in3) +
                   in2->size[0] * in2->size[1] * in2->size[2] * in4]
              .re +
          in1_data[i2 + in1->size[0] * aux_1_1].re;
      b_in2_data[i1 + b_in2->size[0] * i].im =
          in2_data[((i1 * stride_0_0 + in2->size[0] * aux_0_1) +
                    in2->size[0] * in2->size[1] * in3) +
                   in2->size[0] * in2->size[1] * in2->size[2] * in4]
              .im +
          in1_data[i2 + in1->size[0] * aux_1_1].im;
    }
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }
  i = in1->size[0] * in1->size[1];
  in1->size[0] = b_in2->size[0];
  in1->size[1] = b_in2->size[1];
  emxEnsureCapacity_creal_T(sp, in1, i, &v_emlrtRTEI);
  in1_data = in1->data;
  loop_ub = b_in2->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_loop_ub = b_in2->size[0];
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      in1_data[i1 + in1->size[0] * i] = b_in2_data[i1 + b_in2->size[0] * i];
    }
  }
  emxFree_creal_T(sp, &b_in2);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

void coilCombine(const emlrtStack *sp, emxArray_creal_T *im1,
                 emxArray_creal_T *im2)
{
  static const creal_T beta1 = {
      0.0, /* re */
      0.0  /* im */
  };
  ptrdiff_t k_t;
  ptrdiff_t lda_t;
  ptrdiff_t ldb_t;
  ptrdiff_t ldc_t;
  ptrdiff_t m_t;
  ptrdiff_t n_t;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack st;
  emxArray_creal_T *A;
  emxArray_creal_T *Rs;
  emxArray_creal_T *V1;
  emxArray_creal_T *b_Rs;
  emxArray_creal_T *b_im1;
  emxArray_creal_T *b_im2;
  emxArray_creal_T *c_im1;
  emxArray_creal_T *d_im1;
  emxArray_creal_T *e_im1;
  emxArray_creal_T *myfilt;
  emxArray_creal_T *r;
  emxArray_creal_T *r1;
  emxArray_real_T *s1;
  creal_T *A_data;
  creal_T *Rs_data;
  creal_T *V1_data;
  creal_T *b_im1_data;
  creal_T *c_im1_data;
  creal_T *im1_data;
  creal_T *im2_data;
  creal_T *myfilt_data;
  creal_T *r2;
  int32_T c_im2[5];
  int32_T f_im1[2];
  int32_T C;
  int32_T c_loop_ub;
  int32_T d_loop_ub;
  int32_T i;
  int32_T i1;
  int32_T i2;
  int32_T i3;
  int32_T i4;
  int32_T kn;
  int32_T ky;
  int32_T kz;
  int32_T loop_ub;
  int32_T maxdimlen;
  int32_T nx;
  char_T TRANSA1;
  char_T TRANSB1;
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
  im1_data = im1->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  st.site = &emlrtRSI;
  b_st.site = &d_emlrtRSI;
  st.site = &emlrtRSI;
  b_st.site = &d_emlrtRSI;
  st.site = &emlrtRSI;
  b_st.site = &d_emlrtRSI;
  /*  ---------------------------------------------------------------------------------
   */
  /*  Combine images from different coils */
  /*  ---------------------------------------------------------------------------------
   */
  emxInit_creal_T(sp, &b_im1, 5, &j_emlrtRTEI);
  i = b_im1->size[0] * b_im1->size[1] * b_im1->size[2] * b_im1->size[3] *
      b_im1->size[4];
  b_im1->size[0] = im1->size[0];
  b_im1->size[1] = im1->size[1];
  b_im1->size[2] = im1->size[2];
  b_im1->size[3] = im1->size[3];
  b_im1->size[4] = im1->size[4];
  emxEnsureCapacity_creal_T(sp, b_im1, i, &j_emlrtRTEI);
  b_im1_data = b_im1->data;
  loop_ub =
      im1->size[0] * im1->size[1] * im1->size[2] * im1->size[3] * im1->size[4] -
      1;
  for (i = 0; i <= loop_ub; i++) {
    b_im1_data[i] = im1_data[i];
  }
  st.site = &emlrtRSI;
  permute(&st, b_im1, im1);
  im1_data = im1->data;
  emxFree_creal_T(sp, &b_im1);
  /*  Loop over slices */
  emxInit_creal_T(sp, &b_im2, 5, &bb_emlrtRTEI);
  i = b_im2->size[0] * b_im2->size[1] * b_im2->size[2] * b_im2->size[3] *
      b_im2->size[4];
  b_im2->size[0] = im1->size[0];
  b_im2->size[1] = im1->size[1];
  b_im2->size[2] = im1->size[2];
  b_im2->size[3] = 1;
  b_im2->size[4] = im1->size[4];
  emxEnsureCapacity_creal_T(sp, b_im2, i, &k_emlrtRTEI);
  im2_data = b_im2->data;
  loop_ub = im1->size[0] * im1->size[1] * im1->size[2] * im1->size[4];
  for (i = 0; i < loop_ub; i++) {
    im2_data[i].re = 0.0;
    im2_data[i].im = 0.0;
  }
  i = im1->size[2];
  emxInit_creal_T(sp, &r, 5, &cb_emlrtRTEI);
  emxInit_creal_T(sp, &Rs, 4, &n_emlrtRTEI);
  emxInit_creal_T(sp, &myfilt, 1, &r_emlrtRTEI);
  emxInit_creal_T(sp, &c_im1, 4, &db_emlrtRTEI);
  emxInit_creal_T(sp, &r1, 2, &cb_emlrtRTEI);
  emxInit_creal_T(sp, &A, 2, &eb_emlrtRTEI);
  emxInit_real_T(sp, &s1, &cb_emlrtRTEI);
  emxInit_creal_T(sp, &V1, 2, &x_emlrtRTEI);
  emxInit_creal_T(sp, &d_im1, 5, &l_emlrtRTEI);
  emxInit_creal_T(sp, &b_Rs, 4, &p_emlrtRTEI);
  emxInit_creal_T(sp, &e_im1, 2, &s_emlrtRTEI);
  for (kz = 0; kz < i; kz++) {
    real_T d;
    int32_T N_contents;
    int32_T b_loop_ub;
    if (kz + 1 > b_im2->size[2]) {
      emlrtDynamicBoundsCheckR2012b(kz + 1, 1, b_im2->size[2], &p_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    st.site = &b_emlrtRSI;
    if (kz + 1 > im1->size[2]) {
      emlrtDynamicBoundsCheckR2012b(kz + 1, 1, im1->size[2], &q_emlrtBCI, &st);
    }
    /*  Single slice coil combine */
    i1 = d_im1->size[0] * d_im1->size[1] * d_im1->size[2] * d_im1->size[3] *
         d_im1->size[4];
    d_im1->size[0] = im1->size[0];
    d_im1->size[1] = im1->size[1];
    d_im1->size[2] = 1;
    d_im1->size[3] = im1->size[3];
    loop_ub = im1->size[4];
    d_im1->size[4] = im1->size[4];
    emxEnsureCapacity_creal_T(&st, d_im1, i1, &l_emlrtRTEI);
    b_im1_data = d_im1->data;
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_loop_ub = im1->size[3];
      for (i2 = 0; i2 < b_loop_ub; i2++) {
        maxdimlen = im1->size[1];
        for (i3 = 0; i3 < maxdimlen; i3++) {
          nx = im1->size[0];
          for (C = 0; C < nx; C++) {
            b_im1_data[((C + d_im1->size[0] * i3) +
                        d_im1->size[0] * d_im1->size[1] * i2) +
                       d_im1->size[0] * d_im1->size[1] * d_im1->size[3] * i1] =
                im1_data[(((C + im1->size[0] * i3) +
                           im1->size[0] * im1->size[1] * kz) +
                          im1->size[0] * im1->size[1] * im1->size[2] * i2) +
                         im1->size[0] * im1->size[1] * im1->size[2] *
                             im1->size[3] * i1];
          }
        }
      }
    }
    b_st.site = &m_emlrtRSI;
    b_permute(&b_st, d_im1, c_im1);
    b_im1_data = c_im1->data;
    /*  Get image dimensions and set filter size */
    C = c_im1->size[3];
    N_contents = c_im1->size[2];
    /*  Initialize */
    i1 = r->size[0] * r->size[1] * r->size[2] * r->size[3] * r->size[4];
    r->size[0] = c_im1->size[0];
    r->size[1] = c_im1->size[1];
    r->size[2] = 1;
    r->size[3] = 1;
    r->size[4] = c_im1->size[2];
    emxEnsureCapacity_creal_T(&st, r, i1, &m_emlrtRTEI);
    r2 = r->data;
    loop_ub = c_im1->size[0] * c_im1->size[1] * c_im1->size[2];
    for (i1 = 0; i1 < loop_ub; i1++) {
      r2[i1].re = 0.0;
      r2[i1].im = 0.0;
    }
    i1 = Rs->size[0] * Rs->size[1] * Rs->size[2] * Rs->size[3];
    Rs->size[0] = c_im1->size[0];
    Rs->size[1] = c_im1->size[1];
    Rs->size[2] = c_im1->size[3];
    Rs->size[3] = c_im1->size[3];
    emxEnsureCapacity_creal_T(&st, Rs, i1, &n_emlrtRTEI);
    Rs_data = Rs->data;
    loop_ub = c_im1->size[0] * c_im1->size[1] * c_im1->size[3] * c_im1->size[3];
    for (i1 = 0; i1 < loop_ub; i1++) {
      Rs_data[i1].re = 0.0;
      Rs_data[i1].im = 0.0;
    }
    /*  Get correlation matrices */
    i1 = c_im1->size[3];
    for (maxdimlen = 0; maxdimlen < i1; maxdimlen++) {
      for (nx = 0; nx < C; nx++) {
        for (kn = 0; kn < N_contents; kn++) {
          if (kn + 1 > c_im1->size[2]) {
            emlrtDynamicBoundsCheckR2012b(kn + 1, 1, c_im1->size[2],
                                          &m_emlrtBCI, &st);
          }
          if (maxdimlen + 1 > c_im1->size[3]) {
            emlrtDynamicBoundsCheckR2012b(maxdimlen + 1, 1, c_im1->size[3],
                                          &l_emlrtBCI, &st);
          }
          if (kn + 1 > c_im1->size[2]) {
            emlrtDynamicBoundsCheckR2012b(kn + 1, 1, c_im1->size[2],
                                          &k_emlrtBCI, &st);
          }
          if (nx + 1 > c_im1->size[3]) {
            emlrtDynamicBoundsCheckR2012b(nx + 1, 1, c_im1->size[3],
                                          &j_emlrtBCI, &st);
          }
          i2 = A->size[0] * A->size[1];
          A->size[0] = c_im1->size[0];
          loop_ub = c_im1->size[1];
          A->size[1] = c_im1->size[1];
          emxEnsureCapacity_creal_T(&st, A, i2, &o_emlrtRTEI);
          A_data = A->data;
          for (i2 = 0; i2 < loop_ub; i2++) {
            b_loop_ub = c_im1->size[0];
            for (i3 = 0; i3 < b_loop_ub; i3++) {
              A_data[i3 + A->size[0] * i2].re =
                  b_im1_data[((i3 + c_im1->size[0] * i2) +
                              c_im1->size[0] * c_im1->size[1] * kn) +
                             c_im1->size[0] * c_im1->size[1] * c_im1->size[2] *
                                 nx]
                      .re;
              A_data[i3 + A->size[0] * i2].im =
                  -b_im1_data[((i3 + c_im1->size[0] * i2) +
                               c_im1->size[0] * c_im1->size[1] * kn) +
                              c_im1->size[0] * c_im1->size[1] * c_im1->size[2] *
                                  nx]
                       .im;
            }
          }
          f_im1[0] = c_im1->size[0];
          loop_ub = c_im1->size[1];
          f_im1[1] = c_im1->size[1];
          if ((c_im1->size[0] != A->size[0]) ||
              (c_im1->size[1] != A->size[1])) {
            emlrtSizeEqCheckNDErrorR2021b(&f_im1[0], &A->size[0], &e_emlrtECI,
                                          &st);
          }
          if (maxdimlen + 1 > Rs->size[2]) {
            emlrtDynamicBoundsCheckR2012b(maxdimlen + 1, 1, Rs->size[2],
                                          &o_emlrtBCI, &st);
          }
          if (nx + 1 > Rs->size[3]) {
            emlrtDynamicBoundsCheckR2012b(nx + 1, 1, Rs->size[3], &n_emlrtBCI,
                                          &st);
          }
          i2 = e_im1->size[0] * e_im1->size[1];
          e_im1->size[0] = c_im1->size[0];
          e_im1->size[1] = c_im1->size[1];
          emxEnsureCapacity_creal_T(&st, e_im1, i2, &s_emlrtRTEI);
          c_im1_data = e_im1->data;
          for (i2 = 0; i2 < loop_ub; i2++) {
            b_loop_ub = c_im1->size[0];
            for (i3 = 0; i3 < b_loop_ub; i3++) {
              real_T d1;
              real_T d2;
              real_T d3;
              d = b_im1_data[((i3 + c_im1->size[0] * i2) +
                              c_im1->size[0] * c_im1->size[1] * kn) +
                             c_im1->size[0] * c_im1->size[1] * c_im1->size[2] *
                                 maxdimlen]
                      .re;
              d1 = A_data[i3 + A->size[0] * i2].im;
              d2 = b_im1_data[((i3 + c_im1->size[0] * i2) +
                               c_im1->size[0] * c_im1->size[1] * kn) +
                              c_im1->size[0] * c_im1->size[1] * c_im1->size[2] *
                                  maxdimlen]
                       .im;
              d3 = A_data[i3 + A->size[0] * i2].re;
              c_im1_data[i3 + e_im1->size[0] * i2].re = d * d3 - d2 * d1;
              c_im1_data[i3 + e_im1->size[0] * i2].im = d * d1 + d2 * d3;
            }
          }
          b_st.site = &n_emlrtRSI;
          filter2(&b_st, e_im1, A);
          if ((Rs->size[0] != A->size[0]) &&
              ((Rs->size[0] != 1) && (A->size[0] != 1))) {
            emlrtDimSizeImpxCheckR2021b(Rs->size[0], A->size[0], &d_emlrtECI,
                                        &st);
          }
          loop_ub = Rs->size[1];
          if ((Rs->size[1] != A->size[1]) &&
              ((Rs->size[1] != 1) && (A->size[1] != 1))) {
            emlrtDimSizeImpxCheckR2021b(Rs->size[1], A->size[1], &c_emlrtECI,
                                        &st);
          }
          if (maxdimlen + 1 > Rs->size[2]) {
            emlrtDynamicBoundsCheckR2012b(maxdimlen + 1, 1, Rs->size[2],
                                          &i_emlrtBCI, &st);
          }
          if (nx + 1 > Rs->size[3]) {
            emlrtDynamicBoundsCheckR2012b(nx + 1, 1, Rs->size[3], &h_emlrtBCI,
                                          &st);
          }
          if ((Rs->size[0] == A->size[0]) && (Rs->size[1] == A->size[1])) {
            i2 = A->size[0] * A->size[1];
            A->size[0] = Rs->size[0];
            A->size[1] = Rs->size[1];
            emxEnsureCapacity_creal_T(&st, A, i2, &v_emlrtRTEI);
            A_data = A->data;
            for (i2 = 0; i2 < loop_ub; i2++) {
              b_loop_ub = Rs->size[0];
              for (i3 = 0; i3 < b_loop_ub; i3++) {
                A_data[i3 + A->size[0] * i2].re +=
                    Rs_data[((i3 + Rs->size[0] * i2) +
                             Rs->size[0] * Rs->size[1] * maxdimlen) +
                            Rs->size[0] * Rs->size[1] * Rs->size[2] * nx]
                        .re;
                A_data[i3 + A->size[0] * i2].im +=
                    Rs_data[((i3 + Rs->size[0] * i2) +
                             Rs->size[0] * Rs->size[1] * maxdimlen) +
                            Rs->size[0] * Rs->size[1] * Rs->size[2] * nx]
                        .im;
              }
            }
          } else {
            b_st.site = &n_emlrtRSI;
            binary_expand_op(&b_st, A, Rs, maxdimlen, nx);
            A_data = A->data;
          }
          f_im1[0] = Rs->size[0];
          f_im1[1] = Rs->size[1];
          emlrtSubAssignSizeCheckR2012b(&f_im1[0], 2, &A->size[0], 2,
                                        &b_emlrtECI, &st);
          loop_ub = A->size[1];
          for (i2 = 0; i2 < loop_ub; i2++) {
            b_loop_ub = A->size[0];
            for (i3 = 0; i3 < b_loop_ub; i3++) {
              Rs_data[((i3 + Rs->size[0] * i2) +
                       Rs->size[0] * Rs->size[1] * maxdimlen) +
                      Rs->size[0] * Rs->size[1] * Rs->size[2] * nx] =
                  A_data[i3 + A->size[0] * i2];
            }
          }
          if (*emlrtBreakCheckR2012bFlagVar != 0) {
            emlrtBreakCheckR2012b(&st);
          }
        }
        if (*emlrtBreakCheckR2012bFlagVar != 0) {
          emlrtBreakCheckR2012b(&st);
        }
      }
      if (*emlrtBreakCheckR2012bFlagVar != 0) {
        emlrtBreakCheckR2012b(&st);
      }
    }
    /*  Compute and apply filter at each voxel */
    i1 = c_im1->size[0];
    if (c_im1->size[0] - 1 >= 0) {
      i4 = c_im1->size[1] - 1;
      if (c_im1->size[1] - 1 >= 0) {
        c_loop_ub = Rs->size[3];
        d_loop_ub = c_im1->size[3];
      }
    }
    for (kn = 0; kn < i1; kn++) {
      for (ky = 0; ky <= i4; ky++) {
        boolean_T allFiniteA;
        b_st.site = &o_emlrtRSI;
        if (kn + 1 > Rs->size[0]) {
          emlrtDynamicBoundsCheckR2012b(kn + 1, 1, Rs->size[0], &g_emlrtBCI,
                                        &b_st);
        }
        if (ky + 1 > Rs->size[1]) {
          emlrtDynamicBoundsCheckR2012b(ky + 1, 1, Rs->size[1], &f_emlrtBCI,
                                        &b_st);
        }
        i2 = b_Rs->size[0] * b_Rs->size[1] * b_Rs->size[2] * b_Rs->size[3];
        b_Rs->size[0] = 1;
        b_Rs->size[1] = 1;
        b_Rs->size[2] = Rs->size[2];
        b_Rs->size[3] = Rs->size[3];
        emxEnsureCapacity_creal_T(&b_st, b_Rs, i2, &p_emlrtRTEI);
        c_im1_data = b_Rs->data;
        for (i2 = 0; i2 < c_loop_ub; i2++) {
          loop_ub = Rs->size[2];
          for (i3 = 0; i3 < loop_ub; i3++) {
            c_im1_data[i3 + b_Rs->size[2] * i2] =
                Rs_data[((kn + Rs->size[0] * ky) +
                         Rs->size[0] * Rs->size[1] * i3) +
                        Rs->size[0] * Rs->size[1] * Rs->size[2] * i2];
          }
        }
        c_st.site = &o_emlrtRSI;
        squeeze(&c_st, b_Rs, A);
        c_st.site = &rb_emlrtRSI;
        allFiniteA = !anyNonFinite(&c_st, A);
        if (allFiniteA) {
          c_st.site = &sb_emlrtRSI;
          svd(&c_st, A, e_im1, s1, V1);
          c_im1_data = e_im1->data;
        } else {
          i2 = e_im1->size[0] * e_im1->size[1];
          e_im1->size[0] = A->size[0];
          e_im1->size[1] = A->size[1];
          emxEnsureCapacity_creal_T(&b_st, e_im1, i2, &q_emlrtRTEI);
          c_im1_data = e_im1->data;
          loop_ub = A->size[0] * A->size[1];
          for (i2 = 0; i2 < loop_ub; i2++) {
            c_im1_data[i2] = beta1;
          }
          c_st.site = &tb_emlrtRSI;
          svd(&c_st, e_im1, A, s1, V1);
          i2 = e_im1->size[0] * e_im1->size[1];
          e_im1->size[0] = A->size[0];
          e_im1->size[1] = A->size[1];
          emxEnsureCapacity_creal_T(&b_st, e_im1, i2, &t_emlrtRTEI);
          c_im1_data = e_im1->data;
          loop_ub = A->size[0] * A->size[1];
          for (i2 = 0; i2 < loop_ub; i2++) {
            c_im1_data[i2].re = rtNaN;
            c_im1_data[i2].im = 0.0;
          }
        }
        if (e_im1->size[1] < 1) {
          emlrtDynamicBoundsCheckR2012b(1, 1, e_im1->size[1], &e_emlrtBCI, &st);
        }
        i2 = myfilt->size[0];
        myfilt->size[0] = e_im1->size[0];
        emxEnsureCapacity_creal_T(&st, myfilt, i2, &r_emlrtRTEI);
        myfilt_data = myfilt->data;
        loop_ub = e_im1->size[0];
        for (i2 = 0; i2 < loop_ub; i2++) {
          myfilt_data[i2] = c_im1_data[i2];
        }
        if (kn + 1 > c_im1->size[0]) {
          emlrtDynamicBoundsCheckR2012b(kn + 1, 1, c_im1->size[0], &d_emlrtBCI,
                                        &st);
        }
        if (ky + 1 > c_im1->size[1]) {
          emlrtDynamicBoundsCheckR2012b(ky + 1, 1, c_im1->size[1], &c_emlrtBCI,
                                        &st);
        }
        i2 = b_Rs->size[0] * b_Rs->size[1] * b_Rs->size[2] * b_Rs->size[3];
        b_Rs->size[0] = 1;
        b_Rs->size[1] = 1;
        b_Rs->size[2] = c_im1->size[2];
        b_Rs->size[3] = c_im1->size[3];
        emxEnsureCapacity_creal_T(&st, b_Rs, i2, &u_emlrtRTEI);
        c_im1_data = b_Rs->data;
        for (i2 = 0; i2 < d_loop_ub; i2++) {
          loop_ub = c_im1->size[2];
          for (i3 = 0; i3 < loop_ub; i3++) {
            c_im1_data[i3 + b_Rs->size[2] * i2] =
                b_im1_data[((kn + c_im1->size[0] * ky) +
                            c_im1->size[0] * c_im1->size[1] * i3) +
                           c_im1->size[0] * c_im1->size[1] * c_im1->size[2] *
                               i2];
          }
        }
        b_st.site = &p_emlrtRSI;
        squeeze(&b_st, b_Rs, A);
        A_data = A->data;
        i2 = V1->size[0] * V1->size[1];
        V1->size[0] = A->size[1];
        V1->size[1] = A->size[0];
        emxEnsureCapacity_creal_T(&st, V1, i2, &w_emlrtRTEI);
        V1_data = V1->data;
        loop_ub = A->size[0];
        for (i2 = 0; i2 < loop_ub; i2++) {
          b_loop_ub = A->size[1];
          for (i3 = 0; i3 < b_loop_ub; i3++) {
            V1_data[i3 + V1->size[0] * i2] = A_data[i2 + A->size[0] * i3];
          }
        }
        b_st.site = &p_emlrtRSI;
        nx = V1->size[0] * V1->size[1];
        c_st.site = &hc_emlrtRSI;
        if (C <= 0) {
          d = 0.0;
        } else {
          d = (uint32_T)C;
        }
        if (N_contents <= 0) {
          d = 0.0;
        } else {
          d *= (real_T)(uint32_T)N_contents;
        }
        if (!(d <= 2.147483647E+9)) {
          emlrtErrorWithMessageIdR2018a(&c_st, &c_emlrtRTEI,
                                        "Coder:MATLAB:pmaxsize",
                                        "Coder:MATLAB:pmaxsize", 0);
        }
        maxdimlen = V1->size[0];
        if (V1->size[1] > V1->size[0]) {
          maxdimlen = V1->size[1];
        }
        maxdimlen = muIntScalarMax_sint32(nx, maxdimlen);
        if (C > maxdimlen) {
          emlrtErrorWithMessageIdR2018a(
              &b_st, &e_emlrtRTEI, "Coder:toolbox:reshape_emptyReshapeLimit",
              "Coder:toolbox:reshape_emptyReshapeLimit", 0);
        }
        if (N_contents > maxdimlen) {
          emlrtErrorWithMessageIdR2018a(
              &b_st, &e_emlrtRTEI, "Coder:toolbox:reshape_emptyReshapeLimit",
              "Coder:toolbox:reshape_emptyReshapeLimit", 0);
        }
        if (C * N_contents != nx) {
          emlrtErrorWithMessageIdR2018a(
              &b_st, &d_emlrtRTEI, "Coder:MATLAB:getReshapeDims_notSameNumel",
              "Coder:MATLAB:getReshapeDims_notSameNumel", 0);
        }
        if (kn + 1 > r->size[0]) {
          emlrtDynamicBoundsCheckR2012b(kn + 1, 1, r->size[0], &b_emlrtBCI,
                                        &st);
        }
        if (ky + 1 > r->size[1]) {
          emlrtDynamicBoundsCheckR2012b(ky + 1, 1, r->size[1], &emlrtBCI, &st);
        }
        b_st.site = &p_emlrtRSI;
        i2 = V1->size[0] * V1->size[1];
        V1->size[0] = C;
        V1->size[1] = N_contents;
        emxEnsureCapacity_creal_T(&b_st, V1, i2, &x_emlrtRTEI);
        V1_data = V1->data;
        c_st.site = &jc_emlrtRSI;
        if (e_im1->size[0] != C) {
          if ((e_im1->size[0] == 1) || ((C == 1) && (N_contents == 1))) {
            emlrtErrorWithMessageIdR2018a(
                &c_st, &emlrtRTEI,
                "Coder:toolbox:mtimes_noDynamicScalarExpansion",
                "Coder:toolbox:mtimes_noDynamicScalarExpansion", 0);
          } else {
            emlrtErrorWithMessageIdR2018a(
                &c_st, &b_emlrtRTEI, "MATLAB:innerdim", "MATLAB:innerdim", 0);
          }
        }
        c_st.site = &ic_emlrtRSI;
        if ((e_im1->size[0] == 0) || (C == 0) || (N_contents == 0)) {
          i2 = r1->size[0] * r1->size[1];
          r1->size[0] = 1;
          r1->size[1] = N_contents;
          emxEnsureCapacity_creal_T(&c_st, r1, i2, &y_emlrtRTEI);
          c_im1_data = r1->data;
          for (i2 = 0; i2 < N_contents; i2++) {
            c_im1_data[i2] = beta1;
          }
        } else {
          d_st.site = &kc_emlrtRSI;
          e_st.site = &lc_emlrtRSI;
          TRANSB1 = 'N';
          TRANSA1 = 'C';
          m_t = (ptrdiff_t)1;
          n_t = (ptrdiff_t)(uint32_T)N_contents;
          k_t = (ptrdiff_t)e_im1->size[0];
          lda_t = (ptrdiff_t)e_im1->size[0];
          ldb_t = (ptrdiff_t)(uint32_T)C;
          ldc_t = (ptrdiff_t)1;
          i2 = r1->size[0] * r1->size[1];
          r1->size[0] = 1;
          r1->size[1] = N_contents;
          emxEnsureCapacity_creal_T(&e_st, r1, i2, &ab_emlrtRTEI);
          c_im1_data = r1->data;
          zgemm(&TRANSA1, &TRANSB1, &m_t, &n_t, &k_t, (real_T *)&dc,
                (real_T *)&myfilt_data[0], &lda_t, (real_T *)&V1_data[0],
                &ldb_t, (real_T *)&beta1, (real_T *)&c_im1_data[0], &ldc_t);
        }
        c_im2[0] = 1;
        c_im2[1] = 1;
        c_im2[2] = 1;
        c_im2[3] = 1;
        c_im2[4] = r->size[4];
        emlrtSubAssignSizeCheckR2012b(&c_im2[0], 5, &r1->size[0], 2, &emlrtECI,
                                      &st);
        loop_ub = r->size[4];
        for (i2 = 0; i2 < loop_ub; i2++) {
          r2[(kn + r->size[0] * ky) + r->size[0] * r->size[1] * i2] =
              c_im1_data[i2];
        }
        if (*emlrtBreakCheckR2012bFlagVar != 0) {
          emlrtBreakCheckR2012b(&st);
        }
      }
      if (*emlrtBreakCheckR2012bFlagVar != 0) {
        emlrtBreakCheckR2012b(&st);
      }
    }
    c_im2[0] = b_im2->size[0];
    c_im2[1] = b_im2->size[1];
    c_im2[2] = 1;
    c_im2[3] = 1;
    c_im2[4] = b_im2->size[4];
    emlrtSubAssignSizeCheckR2012b(&c_im2[0], 5, &r->size[0], 5, &f_emlrtECI,
                                  (emlrtCTX)sp);
    loop_ub = r->size[4];
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_loop_ub = r->size[1];
      for (i2 = 0; i2 < b_loop_ub; i2++) {
        maxdimlen = r->size[0];
        for (i3 = 0; i3 < maxdimlen; i3++) {
          im2_data[((i3 + b_im2->size[0] * i2) +
                    b_im2->size[0] * b_im2->size[1] * kz) +
                   b_im2->size[0] * b_im2->size[1] * b_im2->size[2] * i1] =
              r2[(i3 + r->size[0] * i2) + r->size[0] * r->size[1] * i1];
        }
      }
    }
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b((emlrtConstCTX)sp);
    }
  }
  emxFree_creal_T(sp, &e_im1);
  emxFree_creal_T(sp, &b_Rs);
  emxFree_creal_T(sp, &d_im1);
  emxFree_creal_T(sp, &V1);
  emxFree_real_T(sp, &s1);
  emxFree_creal_T(sp, &A);
  emxFree_creal_T(sp, &r1);
  emxFree_creal_T(sp, &c_im1);
  emxFree_creal_T(sp, &myfilt);
  emxFree_creal_T(sp, &Rs);
  emxFree_creal_T(sp, &r);
  st.site = &c_emlrtRSI;
  c_permute(&st, b_im2, im2);
  emxFree_creal_T(sp, &b_im2);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

emlrtCTX emlrtGetRootTLSGlobal(void)
{
  return emlrtRootTLSGlobal;
}

void emlrtLockerFunction(EmlrtLockeeFunction aLockee, emlrtConstCTX aTLS,
                         void *aData)
{
  omp_set_lock(&emlrtLockGlobal);
  emlrtCallLockeeFunction(aLockee, aTLS, aData);
  omp_unset_lock(&emlrtLockGlobal);
}

/* End of code generation (coilCombine.c) */
