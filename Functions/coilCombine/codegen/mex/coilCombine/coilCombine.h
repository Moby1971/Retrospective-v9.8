/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * coilCombine.h
 *
 * Code generation for function 'coilCombine'
 *
 */

#pragma once

/* Include files */
#include "coilCombine_types.h"
#include "rtwtypes.h"
#include "emlrt.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
void coilCombine(const emlrtStack *sp, emxArray_creal_T *im1,
                 emxArray_creal_T *im2);

emlrtCTX emlrtGetRootTLSGlobal(void);

void emlrtLockerFunction(EmlrtLockeeFunction aLockee, emlrtConstCTX aTLS,
                         void *aData);

/* End of code generation (coilCombine.h) */
