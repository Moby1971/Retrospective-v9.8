/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * permute.h
 *
 * Code generation for function 'permute'
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
void b_permute(const emlrtStack *sp, const emxArray_creal_T *a,
               emxArray_creal_T *b);

void c_permute(const emlrtStack *sp, const emxArray_creal_T *a,
               emxArray_creal_T *b);

void permute(const emlrtStack *sp, const emxArray_creal_T *a,
             emxArray_creal_T *b);

/* End of code generation (permute.h) */
