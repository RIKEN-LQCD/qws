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
#include <complex>
#include <math.h>

#define _CHECK_


#ifdef __cplusplus
extern "C"{
#endif

  using namespace std;
  extern int vold;

  static inline void qws_allreduce(double *array, const int count)
  {
#ifdef _MPI_
    MPI_Allreduce(MPI_IN_PLACE,(void *)array,count,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD);
#endif
  }

  void ddd_d_(scd_t* out, scd_t* in);

  void bicgstab_precdd_s_h_(scs_t *x, const scs_t *b,
                            const double *tol, int *conviter, const int *maxiter,
                            const int *nsap, const int *nm);

  void prec_s_(scs_t* out, scs_t* in, int* nsap, int* nm);

  void fermi_reorder_s2d_dd_(scd_t *out, const scs_t *in);
  void fermi_reorder_d2s_dd_(scs_t *out, const scd_t *in);
  void fermi_reorder_s2d_dd_renorm_(scd_t *out, const scs_t *in, double *dnorm);
  void fermi_reorder_d2s_dd_renorm_(scs_t *out, const scd_t *in, const double *dnorm);


  void check_residual_dd(scd_t *x, const scd_t *b, const int *iter, const double *bnorm);

  //----------------------------------------------------------------------------------------
  void approx_inv_dirac_op_hf(scd_t *x, const scd_t *b, 
                              const double *tol, int *iter, const int *maxiter,
                              const int *nsap, const int *nm)
  {
    extern int vols;

    scs_t *b_s = (scs_t *)malloc( sizeof(scs_t) * vols*2);
    scs_t *x_s = (scs_t *)malloc( sizeof(scs_t) * vols*2);

    double bnorm;

    fermi_reorder_d2s_dd_renorm_(b_s, b, &bnorm);

    bicgstab_precdd_s_h_(x_s, b_s, tol, iter, maxiter, nsap, nm);

    fermi_reorder_s2d_dd_renorm_(x, x_s, &bnorm);

    free(b_s);
    free(x_s);
  }

  //----------------------------------------------------------------------------------------
  void bicgstab_dd_mix2_hf_(scd_t *x, const scd_t *b,
                            const double *tol, int *conviter, const int *maxiter,
                            const double *tol_s, const int *maxiter_s,
                            const int *nsap, const int *nm)
  {

    scd_t *q  = (scd_t *)malloc( sizeof(scd_t) * vold*2);
    scd_t *r  = (scd_t *)malloc( sizeof(scd_t) * vold*2);
    scd_t *p  = (scd_t *)malloc( sizeof(scd_t) * vold*2);
    scd_t *t  = (scd_t *)malloc( sizeof(scd_t) * vold*2);
    scd_t *u  = (scd_t *)malloc( sizeof(scd_t) * vold*2);
    scd_t *r0 = (scd_t *)malloc( sizeof(scd_t) * vold*2);

    double bnorm, rnorm, xnorm, rtmp0, rtmp1, rtmp2, redu[3];
    complex< double > rho0, rho, beta, omega, alpha, ctmp;

    double ar, ai, br, bi;
    int mult_count_all;
    extern int rank;

    mult_count_all = 0;
    
    ///////////////////////////////
    // compute initial residual 
    //
    // q = Dx
    // r = b - q = b - Dx
    // r0 = r
    // bnorm = |b|^2
    // rho0 = <r0|r> = |r|^2
    //
    ///////////////////////////////
    ddd_d_(q, x);

    rtmp0 = 0.0;
    rtmp1 = 0.0;
#pragma omp parallel for reduction(+:rtmp0, rtmp1)
    for(int i=0; i<vold*2; i++){
      for(int j=0; j<24; j++){
        for(int v=0; v<VLEND; v++){
           r[i].ccs[j].v[v] = b[i].ccs[j].v[v] - q[i].ccs[j].v[v];
           p[i].ccs[j].v[v] = r[i].ccs[j].v[v];
          r0[i].ccs[j].v[v] = r[i].ccs[j].v[v];
          rtmp0 += b[i].ccs[j].v[v]*b[i].ccs[j].v[v];
          rtmp1 += r[i].ccs[j].v[v]*r[i].ccs[j].v[v];
        }
      }
    }
    redu[0] = rtmp0;
    redu[1] = rtmp1;
    qws_allreduce(redu, 2);
    bnorm = redu[0];
    rho0  = complex<double>(redu[1],0);

#ifdef _CHECK_
    rnorm = redu[1];
    if (rank == 0)printf("@D %5d %24.15e\n", mult_count_all, sqrt(rnorm/bnorm));
#endif

    for (int iter=0; iter < (*maxiter);iter++) {

      /////////////////
      // u = Mp
      // q = Au
      /////////////////
      int iter_s = mult_count_all;
      approx_inv_dirac_op_hf(u, p, tol_s, &iter_s, maxiter_s, nsap, nm);
      ddd_d_(q, u);
      mult_count_all += iter_s + 1;

      // alpha = rho0 / <r0,q>
      rtmp0 = 0;
      rtmp1 = 0;
#pragma omp parallel for reduction(+:rtmp0, rtmp1)
      for(int i=0; i<vold*2; i++){
        for(int j=0; j<12; j++){
          for(int v=0; v<VLEND; v++){
            rtmp0 += r0[i].cs[j][0].v[v] * q[i].cs[j][0].v[v] + r0[i].cs[j][1].v[v] * q[i].cs[j][1].v[v];
            rtmp1 += r0[i].cs[j][0].v[v] * q[i].cs[j][1].v[v] - r0[i].cs[j][1].v[v] * q[i].cs[j][0].v[v];
          }
        }
      }
      redu[0] = rtmp0;
      redu[1] = rtmp1;
      qws_allreduce(redu,2);
      ctmp = complex<double>(redu[0], redu[1]);
      alpha = rho0 / ctmp;
      if(abs(ctmp)==0)alpha=0;
      ///////////////////////////
      // x = x + alpha u
      // r = r - alpha q
      // rnorm = |r|^2
      ///////////////////////////
      ar = alpha.real();
      ai = alpha.imag();
      rtmp0 = 0.0;
#pragma omp parallel for reduction(+:rtmp0)
      for(int i=0; i<vold*2; i++){
        for(int j=0; j<12; j++){
          for(int v=0; v<VLEND; v++){
            x[i].cs[j][0].v[v] = x[i].cs[j][0].v[v] + ar *  u[i].cs[j][0].v[v] - ai *  u[i].cs[j][1].v[v];
            x[i].cs[j][1].v[v] = x[i].cs[j][1].v[v] + ar *  u[i].cs[j][1].v[v] + ai *  u[i].cs[j][0].v[v];
            r[i].cs[j][0].v[v] = r[i].cs[j][0].v[v] - ar *  q[i].cs[j][0].v[v] + ai *  q[i].cs[j][1].v[v];
            r[i].cs[j][1].v[v] = r[i].cs[j][1].v[v] - ar *  q[i].cs[j][1].v[v] - ai *  q[i].cs[j][0].v[v];
            rtmp0 += r[i].cs[j][0].v[v] * r[i].cs[j][0].v[v] 
                   + r[i].cs[j][1].v[v] * r[i].cs[j][1].v[v];
          }
        }
      }
      qws_allreduce(&rtmp0,1);
      rnorm = rtmp0;


#ifdef _CHECK_
      // |x|  for check 
      rtmp0 = 0.0;
#pragma omp parallel for reduction(+:rtmp0)
      for(int i=0; i<vold*2; i++){
        for(int j=0; j<24; j++){
          for(int v=0; v<VLEND; v++){
            rtmp0 += x[i].ccs[j].v[v] * x[i].ccs[j].v[v];
          }
        }
      }
      qws_allreduce(&rtmp0,1);
      xnorm = rtmp0;

      if (rank == 0)printf("@D %5d %24.15e %24.15e\n", mult_count_all, sqrt(rnorm/bnorm), xnorm);

      check_residual_dd(x,b,&mult_count_all,&bnorm);
#endif

      if (sqrt(rnorm/bnorm) < *tol) {
        *conviter = iter;
        break;
      }

      //////////////////////////
      // u = Mr
      // t = Au
      //////////////////////////
      iter_s = mult_count_all;
      approx_inv_dirac_op_hf(u, r, tol_s, &iter_s, maxiter_s, nsap, nm);
      ddd_d_(t, u);
      mult_count_all += iter_s + 1;

      //////////////////////////
      // omega = <t,r> / |t|^2
      //////////////////////////
      rtmp0 = 0.0;
      rtmp1 = 0.0;
      rtmp2 = 0.0;
#pragma omp parallel for reduction(+:rtmp0, rtmp1, rtmp2)
      for(int i=0; i<vold*2; i++){
        for(int j=0; j<12; j++){
          for(int v=0; v<VLEND; v++){
            rtmp0 += (t[i].cs[j][0].v[v] * t[i].cs[j][0].v[v] + t[i].cs[j][1].v[v] *  t[i].cs[j][1].v[v]);
            rtmp1 += (t[i].cs[j][0].v[v] * r[i].cs[j][0].v[v] + t[i].cs[j][1].v[v] *  r[i].cs[j][1].v[v]);
            rtmp2 += (t[i].cs[j][0].v[v] * r[i].cs[j][1].v[v] - t[i].cs[j][1].v[v] *  r[i].cs[j][0].v[v]);
          }
        }
      }
      redu[0] = rtmp0;
      redu[1] = rtmp1;
      redu[2] = rtmp2;
      qws_allreduce(redu,3);

      omega = complex<double>(redu[1]/redu[0], redu[2]/redu[0] );
      if(redu[0]==0)omega=0;
      ////////////////////////////
      // x = x + omega u
      // r = r - omega t
      // rnorm = |r|^2
      ////////////////////////////
      ar = omega.real();
      ai = omega.imag();
      rtmp0 = 0.0;
#pragma omp parallel for reduction(+:rtmp0)
      for(int i=0; i<vold*2; i++){
        for(int j=0; j<12; j++){
          for(int v=0; v<VLEND; v++){
            x[i].cs[j][0].v[v] = x[i].cs[j][0].v[v] + ar *  u[i].cs[j][0].v[v] - ai *  u[i].cs[j][1].v[v];
            x[i].cs[j][1].v[v] = x[i].cs[j][1].v[v] + ar *  u[i].cs[j][1].v[v] + ai *  u[i].cs[j][0].v[v];
            r[i].cs[j][0].v[v] = r[i].cs[j][0].v[v] - ar *  t[i].cs[j][0].v[v] + ai *  t[i].cs[j][1].v[v];
            r[i].cs[j][1].v[v] = r[i].cs[j][1].v[v] - ar *  t[i].cs[j][1].v[v] - ai *  t[i].cs[j][0].v[v];
            rtmp0 += r[i].cs[j][0].v[v] * r[i].cs[j][0].v[v] 
                   + r[i].cs[j][1].v[v] * r[i].cs[j][1].v[v];
          }
        }
      }
      qws_allreduce(&rtmp0,1);
      rnorm = rtmp0;

#ifdef _CHECK_
      // |x|  for check 
      rtmp0 = 0.0;
#pragma omp parallel for reduction(+:rtmp0)
      for(int i=0; i<vold*2; i++){
        for(int j=0; j<24; j++){
          for(int v=0; v<VLEND; v++){
            rtmp0 += x[i].ccs[j].v[v] * x[i].ccs[j].v[v];
          }
        }
      }
      qws_allreduce(&rtmp0,1);
      xnorm = rtmp0;

      if (rank==0)printf("@D %5d %24.15e %24.15e\n", mult_count_all, sqrt(rnorm/bnorm), xnorm);

      check_residual_dd(x,b,&mult_count_all,&bnorm);
#endif

      if (sqrt(rnorm/bnorm) < *tol) {
        *conviter = iter;
        break;
      }

      //////////////////
      // rho = <r0,r>
      //////////////////
      rtmp0 = 0.0;
      rtmp1 = 0.0;
#pragma omp parallel for reduction(+:rtmp0, rtmp1)
      for(int i=0; i<vold*2; i++){
        for(int j=0; j<12; j++){
          for(int v=0; v<VLEND; v++){
            rtmp0 += r0[i].cs[j][0].v[v] * r[i].cs[j][0].v[v] + r0[i].cs[j][1].v[v] * r[i].cs[j][1].v[v];
            rtmp1 += r0[i].cs[j][0].v[v] * r[i].cs[j][1].v[v] - r0[i].cs[j][1].v[v] * r[i].cs[j][0].v[v];
          }
        }
      }
      redu[0] = rtmp0;
      redu[1] = rtmp1;
      qws_allreduce(redu,2);
      rho = complex<double>(redu[0],redu[1]);

      beta = alpha*rho/( rho0 * omega);
      rho0 = rho;
      if(abs(omega)==0||abs(rho0)==0)beta=0;
      //////////////////////////
      // p = p - omega q
      // p = r + beta p
      //////////////////////////
      ar = omega.real();
      ai = omega.imag();
      br =  beta.real();
      bi =  beta.imag();
#pragma omp parallel for
      for(int i=0; i<vold*2; i++){
        for(int j=0; j<12; j++){
          for(int v=0; v<VLEND; v++){
            double cr = p[i].cs[j][0].v[v] - ar *  q[i].cs[j][0].v[v] + ai *  q[i].cs[j][1].v[v];
            double ci = p[i].cs[j][1].v[v] - ar *  q[i].cs[j][1].v[v] - ai *  q[i].cs[j][0].v[v];
            p[i].cs[j][0].v[v] = br * cr - bi * ci + r[i].cs[j][0].v[v];
            p[i].cs[j][1].v[v] = br * ci + bi * cr + r[i].cs[j][1].v[v];
          }
        }
      }


    } // end for iter


  } // end of bicgstab


#ifdef __cplusplus
}
#endif
