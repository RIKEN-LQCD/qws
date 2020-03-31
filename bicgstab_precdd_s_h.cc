//****************************************************************************************
//
//  Copyright (c) 2015-2020, Yoshifumi Nakamura <nakamura@riken.jp>
//  Copyright (c) 2015-2020, Yuta Mukai         <mukai.yuta@fujitsu.com>
//  Copyright (c) 2018-2020, Ken-Ichi Ishikawa  <ishikawa@theo.phys.sci.hirosima-u.ac.jp>
//  Copyright (c) 2019-2020, Issaku Kanamori    <kanamori-i@riken.jp>
//
//
//  All rights reserved.
//
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions are
//  met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer. 
//
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer listed
//    in this license in the documentation and/or other materials
//    provided with the distribution.
//
//  * Neither the name of the copyright holders nor the names of its
//    contributors may be used to endorse or promote products derived from
//    this software without specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT  
//  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
//  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
//  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
//  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
//  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
//  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT  
//  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
//  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
//
//----------------------------------------------------------------------------------------
//  ACKNOWLEDGMENT
//
//  This software has been developed in a co-design working group for the lattice QCD
//  supported by MEXT's programs for the Development and Improvement for the Next
//  Generation Ultra High-Speed Computer System, under its Subsidies for Operating the
//  Specific Advanced Large Research Facilities, and Priority Issue 9 
//  (Elucidation of the Fundamental Laws and Evolution of the Universe) to be tackled by
//  using the Supercomputer Fugaku.
//
//****************************************************************************************
#include "qws.h"
#include "qws_h.h"

#include <complex>
#include <math.h>

//#undef _PREC_
#define _PREC_

//#undef _DEBUG_
#define _DEBUG_

//#undef _CHECK_
#define _CHECK_

#ifdef __cplusplus
extern "C"{
#endif

  using namespace std;

  extern int vols;

  static inline void qws_allreduce(float *array, const int count)
  {
#ifdef _MPI_
    MPI_Allreduce(MPI_IN_PLACE,(void *)array,count,MPI_REAL,MPI_SUM,MPI_COMM_WORLD);
#endif
  }

  void check_residual_s_h_(const scs_t *x, const scs_t *b, const int *iter, const float *bnorm, const int *nsap, const int *nm);


  void bicgstab_precdd_s_h_(scs_t *x, const scs_t *b,
                            const double *tol, int *conviter, const int *maxiter,
                            const int *nsap, const int *nm)
  //
  // Flexible BiCGStab solver for Wilson/Clover quark
  //
  //  D x = b
  //
  // is solved in single precision with the SAP preconsitioner in half precision
  //
  //         x : quark field solution (output)
  //         b : quark field source (input)
  //       tol : tolerance for convergence
  //  conviter : count for multiplication of D for convergence
  //   maxiter : maximum intertion count of BiCGStab iteration
  //      nsap : fixed iteration count for SAP preconditioner
  //        nm : fixed iteration count for Jacobbi iteration for approximate inverse of D in a block
  //
  {

    static scs_t *q  = nullptr;
    static scs_t *r  = nullptr;
    static scs_t *p  = nullptr;
    static scs_t *t  = nullptr;
    static scs_t *w  = nullptr;
    static scs_t *r0 = nullptr;

    float bnorm, rnorm, rtmp0, rtmp1, rtmp2, redu[3];
    float ar, ai, br, bi;
    complex< float > rho0, rho, beta, omega, alpha, ctmp;
    int mult_count;
    extern int rank;

    if(nullptr ==  q)  q = (scs_t *)malloc( sizeof(scs_t) * vols*2);
    if(nullptr ==  r)  r = (scs_t *)malloc( sizeof(scs_t) * vols*2);
    if(nullptr ==  p)  p = (scs_t *)malloc( sizeof(scs_t) * vols*2);
    if(nullptr ==  t)  t = (scs_t *)malloc( sizeof(scs_t) * vols*2);
    if(nullptr ==  w)  w = (scs_t *)malloc( sizeof(scs_t) * vols*2);
    if(nullptr == r0) r0 = (scs_t *)malloc( sizeof(scs_t) * vols*2);

    /////////////////////////////
    //
    //  x = 0 (zero start)
    //  r = b
    //  p = b
    // r0 = b
    // rho0 = <r0|r> = |r|^2
    // bnorm = |b|^2 = |r|^2
    //
    /////////////////////////////
    rtmp0 = 0.0f;
#pragma omp parallel for reduction(+:rtmp0)
    for(int i=0; i<vols*2; i++){
      for(int j=0; j<24; j++){
        for(int v=0; v<VLENS; v++){
           x[i].ccs[j].v[v] = 0.0f;
           r[i].ccs[j].v[v] = b[i].ccs[j].v[v];
           p[i].ccs[j].v[v] = b[i].ccs[j].v[v];
          r0[i].ccs[j].v[v] = b[i].ccs[j].v[v];
          rtmp0 += b[i].ccs[j].v[v] * b[i].ccs[j].v[v];
        }
      }
    }
    qws_allreduce(&rtmp0,1);
    bnorm = rtmp0;
    rho0 = complex<float>(rtmp0,0.0f);

#ifdef _DEBUG_
    mult_count = 0;
#ifdef _CHECK_
    check_residual_s_h_(x,b,&mult_count,&bnorm,nsap,nm);
#endif
    double err = 1.0e0;
    if (0==rank) printf("@@ %5d %12.7e\n",mult_count,err);
#endif

    for (int iter=0; iter<(*maxiter);iter++){

      ////////////////////
#ifdef _PREC_
      // w = Msap p
      // q = D w
      ////////////////////
      assign_mult_mSAP_s_h_( w, p, nsap, nm);
      assign_mult_wd_s_(     q, w);
#else
      // q = D p
      assign_mult_wd_s_(     q, p);
#endif

      //////////////////////////
      // alpha = rho0 / <r0,q>
      //////////////////////////
      rtmp0 = 0.0f;
      rtmp1 = 0.0f;
#pragma omp parallel for reduction(+:rtmp0, rtmp1)
      for(int i=0; i<vols*2; i++){
        for(int j=0; j<12; j++){
          for(int v=0; v<VLENS; v++){
            rtmp0 += r0[i].cs[j][0].v[v] * q[i].cs[j][0].v[v] + r0[i].cs[j][1].v[v] * q[i].cs[j][1].v[v];
            rtmp1 += r0[i].cs[j][0].v[v] * q[i].cs[j][1].v[v] - r0[i].cs[j][1].v[v] * q[i].cs[j][0].v[v];
          }
        }
      }
      redu[0] = rtmp0;
      redu[1] = rtmp1;
      qws_allreduce(redu,2);
      ctmp = complex<float>(redu[0], redu[1]);
      alpha = rho0 / ctmp;
      if(abs(ctmp)<=1e-15)alpha=0;
      ////////////////////////////
      // x = x + alpha w
      // r = r - alpha q
      // rnorm = |r|^2
      ////////////////////////////
      ar = alpha.real();
      ai = alpha.imag();
      rtmp0 = 0.0f;
#pragma omp parallel for reduction(+:rtmp0)
      for(int i=0; i<vols*2; i++){
        for(int j=0; j<12; j++){
          for(int v=0; v<VLENS; v++){

#ifdef _PREC_
            x[i].cs[j][0].v[v] += + ar * w[i].cs[j][0].v[v] - ai * w[i].cs[j][1].v[v];
            x[i].cs[j][1].v[v] += + ar * w[i].cs[j][1].v[v] + ai * w[i].cs[j][0].v[v];
#else
            x[i].cs[j][0].v[v] += + ar * p[i].cs[j][0].v[v] - ai * p[i].cs[j][1].v[v];
            x[i].cs[j][1].v[v] += + ar * p[i].cs[j][1].v[v] + ai * p[i].cs[j][0].v[v];
#endif

            r[i].cs[j][0].v[v] += - ar * q[i].cs[j][0].v[v] + ai * q[i].cs[j][1].v[v];
            r[i].cs[j][1].v[v] += - ar * q[i].cs[j][1].v[v] - ai * q[i].cs[j][0].v[v];

            rtmp0 += r[i].cs[j][0].v[v] * r[i].cs[j][0].v[v]
                   + r[i].cs[j][1].v[v] * r[i].cs[j][1].v[v];

          }
        }
      }
      qws_allreduce(&rtmp0,1);
      rnorm = rtmp0;

#ifdef _DEBUG_
      mult_count++;
#ifdef _CHECK_
      check_residual_s_h_(x,b,&mult_count,&bnorm,nsap,nm);
#endif
      err = (double)sqrt(rnorm/bnorm);
      if (0==rank) printf("@@ %5d %12.7e\n",mult_count,err);
#endif
      if (sqrt(rnorm/bnorm) < *tol){
        *conviter = mult_count;
	break;
      }

      ///////////////////////////
      // w = Msap r
      // t = D w
      ///////////////////////////
#ifdef _PREC_
      assign_mult_mSAP_s_h_( w, r, nsap, nm);
      assign_mult_wd_s_(     t, w);
#else
      assign_mult_wd_s_(     t, r);
#endif

      ///////////////////////////
      // omega = <t,r> / |t|^2
      ///////////////////////////
      rtmp0 = 0.0f;
      rtmp1 = 0.0f;
      rtmp2 = 0.0f;
#pragma omp parallel for reduction(+:rtmp0, rtmp1, rtmp2)
      for(int i=0; i<vols*2; i++){
        for(int j=0; j<12; j++){
          for(int v=0; v<VLENS; v++){
            rtmp0 +=  t[i].cs[j][0].v[v] * t[i].cs[j][0].v[v] 
                    + t[i].cs[j][1].v[v] * t[i].cs[j][1].v[v];

            rtmp1 +=  t[i].cs[j][0].v[v] * r[i].cs[j][0].v[v]
                    + t[i].cs[j][1].v[v] * r[i].cs[j][1].v[v];
            rtmp2 +=  t[i].cs[j][0].v[v] * r[i].cs[j][1].v[v] 
                    - t[i].cs[j][1].v[v] * r[i].cs[j][0].v[v];
          }
        }
      }
      redu[0] = rtmp0;
      redu[1] = rtmp1;
      redu[2] = rtmp2;
      qws_allreduce(redu,3);
      omega = complex<float>(redu[1]/redu[0], redu[2]/redu[0] );
      if(redu[0]<=1e-30)omega=0;
      /////////////////////////
      // x = x + omega w
      // r = r - omega t
      // rnorm = |r|^2
      /////////////////////////
      ar = omega.real();
      ai = omega.imag();
      rtmp0 = 0.0f;
#pragma omp parallel for reduction(+:rtmp0)
      for(int i=0; i<vols*2; i++){
        for(int j=0; j<12; j++){
          for(int v=0; v<VLENS; v++){

#ifdef _PREC_
            x[i].cs[j][0].v[v] += + ar * w[i].cs[j][0].v[v] - ai * w[i].cs[j][1].v[v];
            x[i].cs[j][1].v[v] += + ar * w[i].cs[j][1].v[v] + ai * w[i].cs[j][0].v[v];
#else
            x[i].cs[j][0].v[v] += + ar * r[i].cs[j][0].v[v] - ai * r[i].cs[j][1].v[v];
            x[i].cs[j][1].v[v] += + ar * r[i].cs[j][1].v[v] + ai * r[i].cs[j][0].v[v];
#endif

            r[i].cs[j][0].v[v] += - ar * t[i].cs[j][0].v[v] + ai * t[i].cs[j][1].v[v];
            r[i].cs[j][1].v[v] += - ar * t[i].cs[j][1].v[v] - ai * t[i].cs[j][0].v[v];

            rtmp0 += r[i].cs[j][0].v[v] * r[i].cs[j][0].v[v]
                   + r[i].cs[j][1].v[v] * r[i].cs[j][1].v[v];
          }
        }
      }
      qws_allreduce(&rtmp0,1);
      rnorm = rtmp0;

#ifdef _DEBUG_
      mult_count++;
#ifdef _CHECK_
      check_residual_s_h_(x,b,&mult_count,&bnorm,nsap,nm);
#endif
      err = (double)sqrt(rnorm/bnorm);
      if (0==rank) printf("@@ %5d %12.7e\n",mult_count,err);
#endif
      if (sqrt(rnorm/bnorm) < *tol){
        *conviter = mult_count;
	break;
      }

      /////////////////////
      // rho = <r0|r>
      /////////////////////
      rtmp0 = 0.0f;
      rtmp1 = 0.0f;
#pragma omp parallel for reduction(+:rtmp0, rtmp1)
      for(int i=0; i<vols*2; i++){
        for(int j=0; j<12; j++){
          for(int v=0; v<VLENS; v++){
            rtmp0 += r0[i].cs[j][0].v[v] * r[i].cs[j][0].v[v] + r0[i].cs[j][1].v[v] * r[i].cs[j][1].v[v];
            rtmp1 += r0[i].cs[j][0].v[v] * r[i].cs[j][1].v[v] - r0[i].cs[j][1].v[v] * r[i].cs[j][0].v[v];
          }
        }
      }
      redu[0] = rtmp0;
      redu[1] = rtmp1;
      qws_allreduce(redu,2);
      rho = complex<float>(redu[0],redu[1]);

      beta = alpha*rho/( rho0 * omega);
      rho0 = rho;
      if(abs(omega)==0||abs(rho0)==0)beta=0;
      //////////////////////
      // p = p - omega q
      // p = r + beta p
      //////////////////////
      ar = omega.real();
      ai = omega.imag();
      br =  beta.real();
      bi =  beta.imag();
#pragma omp parallel for
      for(int i=0; i<vols*2; i++){
        for(int j=0; j<12; j++){
          for(int v=0; v<VLENS; v++){
            float cr = p[i].cs[j][0].v[v] - ar * q[i].cs[j][0].v[v] + ai * q[i].cs[j][1].v[v];
            float ci = p[i].cs[j][1].v[v] - ar * q[i].cs[j][1].v[v] - ai * q[i].cs[j][0].v[v];
            p[i].cs[j][0].v[v] = br * cr - bi * ci + r[i].cs[j][0].v[v];
            p[i].cs[j][1].v[v] = br * ci + bi * cr + r[i].cs[j][1].v[v];
          }
        }
      }

    } // end for iter

#ifdef _DEBUG_
    check_residual_s_h_(x,b,&mult_count,&bnorm,nsap,nm);
#endif

  } // end of bicgstab


void check_residual_s_h_(const scs_t *x, const scs_t *b, const int *iter, const float *bnorm, const int *nsap, const int *nm)
//
// check |D x - b|/|b|
//
{
  extern int rank;

  scs_t *w = (scs_t*)malloc( sizeof(scs_t) * vols*2);
  scs_t *s = (scs_t*)malloc( sizeof(scs_t) * vols*2);

  //
  // check residual
  //
  // w = D x
  //
  assign_mult_wd_s_( w, x);

  //
  // s = b - w
  //
  // err = |s|/|b|
  //
  float rnorm = 0.0f;
#pragma omp parallel for reduction(+:rnorm)
  for(int i=0; i<vols*2; i++){
    for(int j=0; j<24; j++){
      for(int v=0; v<VLENS; v++){
        s[i].ccs[j].v[v] = b[i].ccs[j].v[v] - w[i].ccs[j].v[v];
        rnorm += s[i].ccs[j].v[v] * s[i].ccs[j].v[v];
      }
    }
  }
  qws_allreduce(&rnorm,1);

  float err = sqrt(rnorm/(*bnorm));
  if (0 == rank) printf("@C %5d %12.7e\n",*iter,err);

  free(w);
  free(s);
  w = nullptr;
  s = nullptr;

}

#ifdef __cplusplus
}
#endif
