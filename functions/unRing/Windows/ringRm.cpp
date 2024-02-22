
/*
LICENCE

The Software remains the property of the University Medical Center
Freiburg ("the University").

The Software is distributed "AS IS" under this Licence solely for
non-commercial use in the hope that it will be useful, but in order
that the University as a charitable foundation protects its assets for
the benefit of its educational and research purposes, the University
makes clear that no condition is made or to be implied, nor is any
warranty given or to be implied, as to the accuracy of the Software,
or that it will be suitable for any particular purpose or for use
under any specific conditions. Furthermore, the University disclaims
all responsibility for the use which is made of the Software. It
further disclaims any liability for the outcomes arising from using
the Software.

The Licensee agrees to indemnify the University and hold the
University harmless from and against any and all claims, damages and
liabilities asserted by third parties (including claims for
negligence) which arise directly or indirectly from the use of the
Software or the sale of any products based on the Software.

No part of the Software may be reproduced, modified, transmitted or
transferred in any form or by any means, electronic or mechanical,
without the express permission of the University. The permission of
the University is not required if the said reproduction, modification,
transmission or transference is done without financial return, the
conditions of this Licence are imposed upon the receiver of the
product, and all original and amended source code is included in any
transmitted product. You may be held legally responsible for any
copyright infringement that is caused or encouraged by your failure to
abide by these terms and conditions.

You are not permitted under this Licence to use this Software
commercially. Use for which any financial return is received shall be
defined as commercial use, and includes (1) integration of all or part
of the source code or the Software into a product for sale or license
by or on behalf of Licensee to third parties or (2) use of the
Software or any derivative of it for research with the final aim of
developing software products for sale or license to a third party or
(3) use of the Software or any derivative of it for research with the
final aim of developing non-software products for sale or license to a
third party, or (4) use of the Software to provide any service to an
external organisation for which payment is received. */


#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits>

#include <mex.h>
#include <matrix.h>
#include <math.h>
#include <map>
#include <string.h>

// for windows
#include <iostream> 
#include <algorithm>

// #include "/usr/local/include/fftw3.h"

#include "fftw3.h"

using namespace std;

#define PI  3.1416

// mex -v -L/usr/local/lib -lfftw3 -I/usr/local/include/ ringRm.cpp
// mex -compatibleArrayDims -lfftw3-3 ringRm.cpp

void unring_1D(fftw_complex *data,int n, int numlines,int nsh,int minW, int maxW)
{


    fftw_complex *in, *out;
    fftw_plan p,pinv;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
    p = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    pinv = fftw_plan_dft_1d(n, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_complex *sh = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n *(2*nsh+1));
    fftw_complex *sh2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n *(2*nsh+1));


    double nfac = 1/double(n);

    int *shifts = (int*) malloc(sizeof(int)*(2*nsh+1));
    shifts[0] = 0;
    for (int j = 0; j < nsh;j++)
    {
        shifts[j+1] = j+1;
        shifts[1+nsh+j] = -(j+1);
    }
   
//    double TV1arr[2*nsh+1];
//    double TV2arr[2*nsh+1];
  
    double *TV1arr = (double*) malloc(sizeof(double)*(2*nsh+1));
    double *TV2arr = (double*) malloc(sizeof(double)*(2*nsh+1));

    for (int k = 0; k < numlines; k++)
    {


        fftw_execute_dft(p,&(data[n*k]),sh);

        int maxn = (n%2 == 1)? (n-1)/2 : n/2 -1;

        for (int j = 1; j < 2*nsh+1; j++)
        {
            double phi = PI/double(n) * double(shifts[j])/double(nsh);
            fftw_complex u = {cos(phi),sin(phi)};
            fftw_complex e = {1,0};

            sh[j*n ][0] = sh[0][0];
            sh[j*n ][1] = sh[0][1];

            if (n%2 == 0)
            {
                sh[j*n + n/2][0] = 0;
                sh[j*n + n/2][1] = 0;
            }

            for (int l = 0; l < maxn; l++)
            {

                double tmp = e[0];
                e[0] = u[0]*e[0] - u[1]*e[1];
                e[1] = tmp*u[1] + u[0]*e[1];

                int L ;
                L = l+1;
                sh[j*n +L][0] = (e[0]*sh[L][0] - e[1]*sh[L][1]);
                sh[j*n +L][1] = (e[0]*sh[L][1] + e[1]*sh[L][0]);
                L = n-1-l;
                sh[j*n +L][0] = (e[0]*sh[L][0] + e[1]*sh[L][1]);
                sh[j*n +L][1] = (e[0]*sh[L][1] - e[1]*sh[L][0]);

            }
        }


        for (int j = 0; j < 2*nsh+1; j++)
        {
            fftw_execute_dft(pinv,&(sh[j*n]),&sh2[j*n]);
        }

        for (int j=0;j < 2*nsh+1;j++)
        {
            TV1arr[j] = 0;
            TV2arr[j] = 0;
            const int l = 0;
            for (int t = minW; t <= maxW;t++)
            {
                TV1arr[j] += fabs(sh2[j*n + (l-t+n)%n ][0] - sh2[j*n + (l-(t+1)+n)%n ][0]);
                TV1arr[j] += fabs(sh2[j*n + (l-t+n)%n ][1] - sh2[j*n + (l-(t+1)+n)%n ][1]);
                TV2arr[j] += fabs(sh2[j*n + (l+t+n)%n ][0] - sh2[j*n + (l+(t+1)+n)%n ][0]);
                TV2arr[j] += fabs(sh2[j*n + (l+t+n)%n ][1] - sh2[j*n + (l+(t+1)+n)%n ][1]);
            }
        }




        for(int l=0; l < n; l++)
        {
            double minTV = 999999999999;
            int minidx= 0;
            for (int j=0;j < 2*nsh+1;j++)
            {

                if (TV1arr[j] < minTV)
                {
                    minTV = TV1arr[j];
                    minidx = j;
                }
                if (TV2arr[j] < minTV)
                {
                    minTV = TV2arr[j];
                    minidx = j;
                }

                TV1arr[j] += fabs(sh2[j*n + (l-minW+1+n)%n ][0] - sh2[j*n + (l-(minW)+n)%n ][0]);
                TV1arr[j] -= fabs(sh2[j*n + (l-maxW+n)%n ][0] - sh2[j*n + (l-(maxW+1)+n)%n ][0]);
                TV2arr[j] += fabs(sh2[j*n + (l+maxW+1+n)%n ][0] - sh2[j*n + (l+(maxW+2)+n)%n ][0]);
                TV2arr[j] -= fabs(sh2[j*n + (l+minW+n)%n ][0] - sh2[j*n + (l+(minW+1)+n)%n ][0]);

                TV1arr[j] += fabs(sh2[j*n + (l-minW+1+n)%n ][1] - sh2[j*n + (l-(minW)+n)%n ][1]);
                TV1arr[j] -= fabs(sh2[j*n + (l-maxW+n)%n ][1] - sh2[j*n + (l-(maxW+1)+n)%n ][1]);
                TV2arr[j] += fabs(sh2[j*n + (l+maxW+1+n)%n ][1] - sh2[j*n + (l+(maxW+2)+n)%n ][1]);
                TV2arr[j] -= fabs(sh2[j*n + (l+minW+n)%n ][1] - sh2[j*n + (l+(minW+1)+n)%n ][1]);

            }


            double a0r = sh2[minidx*n + (l-1+n)%n ][0];
            double a1r = sh2[minidx*n + l][0];
            double a2r = sh2[minidx*n + (l+1+n)%n ][0];
            double a0i = sh2[minidx*n + (l-1+n)%n ][1];
            double a1i = sh2[minidx*n + l][1];
            double a2i = sh2[minidx*n + (l+1+n)%n ][1];
            double s = double(shifts[minidx])/nsh/2;

            //data[k*n + l][0] =  (a1r - 0.5*(a2r-a0r)*s + (0.5*(a2r+a0r) - a1r)*s*s)*nfac;
            //data[k*n + l][1] =  (a1i - 0.5*(a2i-a0i)*s + (0.5*(a2i+a0i) - a1i)*s*s)*nfac;


            if (s>0)
            {
                data[k*n + l][0] =  (a1r*(1-s) + a0r*s)*nfac;
                data[k*n + l][1] =  (a1i*(1-s) + a0i*s)*nfac;
            }
            else
            {
                s = -s;
                data[k*n + l][0] =  (a1r*(1-s) + a2r*s)*nfac;
                data[k*n + l][1] =  (a1i*(1-s) + a2i*s)*nfac;
            }

        }



    }



    free(shifts);
    fftw_destroy_plan(p);
    fftw_destroy_plan(pinv);
    fftw_free(in);
    fftw_free(out);
    fftw_free(sh);
    fftw_free(sh2);
    free(TV1arr);
    free(TV2arr);




}


// void unring_2D(fftw_complex *data1,fftw_complex *tmp2, const int *dim_sz, int nsh, int minW, int maxW)

// void unring_2D(fftw_complex *data1,fftw_complex *tmp2, const unsigned long *dim_sz, int nsh, int minW, int maxW)
void unring_2D(fftw_complex *data1,fftw_complex *tmp2, const mwSize *dim_sz, int nsh, int minW, int maxW)
{


    double eps = 0;
    fftw_complex *tmp1 =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim_sz[0]*dim_sz[1]);
    fftw_complex *data2 =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim_sz[0]*dim_sz[1]);

    fftw_plan p,pinv,p_tr,pinv_tr;
    p = fftw_plan_dft_2d(dim_sz[1],dim_sz[0], data1, tmp1, FFTW_FORWARD, FFTW_ESTIMATE);
    pinv = fftw_plan_dft_2d(dim_sz[1],dim_sz[0], data1, tmp1, FFTW_BACKWARD, FFTW_ESTIMATE);
    p_tr = fftw_plan_dft_2d(dim_sz[0],dim_sz[1], data2, tmp2, FFTW_FORWARD, FFTW_ESTIMATE);
    pinv_tr = fftw_plan_dft_2d(dim_sz[0],dim_sz[1], data2, tmp2, FFTW_BACKWARD, FFTW_ESTIMATE);
    double nfac = 1/double(dim_sz[0]*dim_sz[1]);

    for (int k = 0 ; k < dim_sz[1];k++)
        for (int j = 0 ; j < dim_sz[0];j++)
        {
            data2[j*dim_sz[1]+k][0] = data1[k*dim_sz[0]+j][0];
            data2[j*dim_sz[1]+k][1] = data1[k*dim_sz[0]+j][1];
        }

    fftw_execute_dft(p,data1,tmp1);
    fftw_execute_dft(p_tr,data2,tmp2);

    for (int k = 0 ; k < dim_sz[1];k++)
    {
        double ck = (1+cos(2*PI*(double(k)/dim_sz[1])))*0.5 +eps;
        for (int j = 0 ; j < dim_sz[0];j++)
        {
            double cj = (1+cos(2*PI*(double(j)/dim_sz[0])))*0.5 +eps;
            tmp1[k*dim_sz[0]+j][0] = nfac*(tmp1[k*dim_sz[0]+j][0] * ck) / (ck+cj);
            tmp1[k*dim_sz[0]+j][1] = nfac*(tmp1[k*dim_sz[0]+j][1] * ck) / (ck+cj);
            tmp2[j*dim_sz[1]+k][0] = nfac*(tmp2[j*dim_sz[1]+k][0] * cj) / (ck+cj);
            tmp2[j*dim_sz[1]+k][1] = nfac*(tmp2[j*dim_sz[1]+k][1] * cj) / (ck+cj);
        }
    }

    fftw_execute_dft(pinv,tmp1,data1);
    fftw_execute_dft(pinv_tr,tmp2,data2);

    unring_1D(data1,dim_sz[0],dim_sz[1],nsh,minW,maxW);
    unring_1D(data2,dim_sz[1],dim_sz[0],nsh,minW,maxW);


    fftw_execute_dft(p,data1,tmp1);
    fftw_execute_dft(p_tr,data2,tmp2);


    for (int k = 0 ; k < dim_sz[1];k++)
    {
        double ck = (1+cos(2*PI*(double(k)/dim_sz[1])))*0.5 +eps;
        for (int j = 0 ; j < dim_sz[0];j++)
        {
            double cj = (1+cos(2*PI*(double(j)/dim_sz[0])))*0.5 +eps;
            tmp1[k*dim_sz[0]+j][0] = nfac*(tmp1[k*dim_sz[0]+j][0]  + tmp2[j*dim_sz[1]+k][0] ) ;
            tmp1[k*dim_sz[0]+j][1] = nfac*(tmp1[k*dim_sz[0]+j][1]  + tmp2[j*dim_sz[1]+k][1] ) ;
            //           tmp1[k*dim_sz[0]+j][0] = nfac*(tmp1[k*dim_sz[0]+j][0]  + tmp2[j*dim_sz[1]+k][0] ) /(ck+cj);
            //           tmp1[k*dim_sz[0]+j][1] = nfac*(tmp1[k*dim_sz[0]+j][1]  + tmp2[j*dim_sz[1]+k][1] ) /(ck+cj);
            //                 tmp1[k*dim_sz[0]+j][0] = nfac*(tmp1[k*dim_sz[0]+j][0]*ck  + tmp2[j*dim_sz[1]+k][0]*cj ) /(ck+cj);
            //                 tmp1[k*dim_sz[0]+j][1] = nfac*(tmp1[k*dim_sz[0]+j][1]*ck  + tmp2[j*dim_sz[1]+k][1]*cj ) /(ck+cj);
        }
    }

    fftw_execute_dft(pinv,tmp1,tmp2);

    fftw_free(data2);
    fftw_free(tmp1);
}






void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{


    int pcnt = 0;

    //////////////////////////// get all the arguments from MATLAB

    const mxArray *Img = prhs[pcnt++];
    const int numdim = mxGetNumberOfDimensions(Img);
    //    const int *dim_sz = mxGetDimensions(Img);
    // const unsigned long *dim_sz = mxGetDimensions(Img);
    const mwSize *dim_sz = mxGetDimensions(Img);
    double *data = (double*) mxGetData(Img);
    double *data_i = (double*) mxGetImagData(Img);

    const mxArray *Params = prhs[pcnt++];
    double *params = mxGetPr(Params);


    int nsh = int(params[2]);
    int minW = int(params[0]);
    int maxW = int(params[1]);

    if (numdim == 2)
    {

        fftw_complex *data_complex =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim_sz[0]*dim_sz[1]);
        fftw_complex *res_complex =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim_sz[0]*dim_sz[1]);
        if (data_i == 0)
        {
            plhs[0] = mxCreateNumericArray(numdim,dim_sz,mxGetClassID(Img),mxREAL);
            for (int k = 0 ; k < dim_sz[1];k++)
                for (int j = 0 ; j < dim_sz[0];j++)
                {
                    data_complex[k*dim_sz[0]+j][0] = data[k*dim_sz[0]+j];
                    data_complex[k*dim_sz[0]+j][1] = 0;
                }
        }
        else
        {
            plhs[0] = mxCreateNumericArray(numdim,dim_sz,mxGetClassID(Img),mxCOMPLEX);
            for (int k = 0 ; k < dim_sz[1];k++)
                for (int j = 0 ; j < dim_sz[0];j++)
                {
                    data_complex[k*dim_sz[0]+j][0] = data[k*dim_sz[0]+j];
                    data_complex[k*dim_sz[0]+j][1] = data_i[k*dim_sz[0]+j];
                }
        }

        unring_2D(data_complex,res_complex, dim_sz,nsh,minW,maxW);

        double *res =  (double*) mxGetData(plhs[0]);
        double *res_i =  (double*) mxGetImagData(plhs[0]);

        if (res_i == 0)
        {
            for (int k = 0 ; k < dim_sz[1];k++)
                for (int j = 0 ; j < dim_sz[0];j++)
                    res[k*dim_sz[0]+j] = res_complex[k*dim_sz[0]+j][0];
        }
        else
        {
            for (int k = 0 ; k < dim_sz[1];k++)
                for (int j = 0 ; j < dim_sz[0];j++)
                {
                    res[k*dim_sz[0]+j] = res_complex[k*dim_sz[0]+j][0];
                    res_i[k*dim_sz[0]+j] = res_complex[k*dim_sz[0]+j][1];
                }
        }

        fftw_free(data_complex);
        fftw_free(res_complex);
    }
    else if (numdim >= 3)
    {
        int prodsz = 1;
        for (int k = 2; k < numdim;k++)
            prodsz *= dim_sz[k];


        fftw_complex *data_complex =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim_sz[0]*dim_sz[1]*prodsz);
        fftw_complex *res_complex  =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dim_sz[0]*dim_sz[1]*prodsz);
        if (data_i == 0)
        {
            plhs[0] = mxCreateNumericArray(numdim,dim_sz,mxGetClassID(Img),mxREAL);
            for (int k = 0 ; k < dim_sz[1]*prodsz;k++)
                for (int j = 0 ; j < dim_sz[0];j++)
                {
                    data_complex[k*dim_sz[0]+j][0] = data[k*dim_sz[0]+j];
                    data_complex[k*dim_sz[0]+j][1] = 0;
                }
        }
        else
        {
            plhs[0] = mxCreateNumericArray(numdim,dim_sz,mxGetClassID(Img),mxCOMPLEX);
            for (int k = 0 ; k < dim_sz[1]*prodsz;k++)
                for (int j = 0 ; j < dim_sz[0];j++)
                {
                    data_complex[k*dim_sz[0]+j][0] = data[k*dim_sz[0]+j];
                    data_complex[k*dim_sz[0]+j][1] = data_i[k*dim_sz[0]+j];
                }
        }

        for (int k = 0; k < prodsz; k++)
            unring_2D(&(data_complex[k*dim_sz[0]*dim_sz[1]]),&(res_complex[k*dim_sz[0]*dim_sz[1]]), dim_sz,nsh,minW,maxW);

        double *res =  (double*) mxGetData(plhs[0]);
        double *res_i =  (double*) mxGetImagData(plhs[0]);

        if (res_i == 0)
        {
            for (int k = 0 ; k < dim_sz[1]*prodsz;k++)
                for (int j = 0 ; j < dim_sz[0];j++)
                    res[k*dim_sz[0]+j] = res_complex[k*dim_sz[0]+j][0];
        }
        else
        {
            for (int k = 0 ; k < dim_sz[1]*prodsz;k++)
                for (int j = 0 ; j < dim_sz[0];j++)
                {
                    res[k*dim_sz[0]+j] = res_complex[k*dim_sz[0]+j][0];
                    res_i[k*dim_sz[0]+j] = res_complex[k*dim_sz[0]+j][1];
                }
        }

        fftw_free(data_complex);
        fftw_free(res_complex);
    }


   


}


