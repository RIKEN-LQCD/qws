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
#include "qwsintrin.h"
#include <complex>
#include <math.h>
//#define _DEBUG_
//#define _CHECK_


#ifdef __cplusplus
extern "C"{
#endif

  using std::complex;
  extern int rank, vold, vols;
  extern void ddd_d_(scd_t* out, scd_t* in);
  extern void bicgstab_precdd_s_(scs_t* x, scs_t* b, double* tol, int* iter, int* maxiter, int* nsap, int* nm);
  extern void prec_s_(scs_t* out, scs_t* in, int* nsap, int* nm);
  extern void fermi_reorder_d2s_dd_(scs_t *out, scd_t *in);
  extern void fermi_reorder_s2d_dd_(scd_t *out, scs_t *in);
  extern void fermi_reorder_s2d_dd_renorm_(scd_t *out, const scs_t *in, double *dnorm);
  extern void fermi_reorder_d2s_dd_renorm_(scs_t *out, const scd_t *in, const double *dnorm);


  extern void check_residual_dd(scd_t *x, const scd_t *b, const int *iter, const double *bnorm);

  //----------------------------------------------------------------------------------------
  void approx_inv_dirac_op(scd_t* x, scd_t* b, double *tol, int* iter, int* maxiter, int* nsap, int* nm){
    __attribute__((aligned(64))) scs_t* b_s = (scs_t*)malloc( sizeof(scs_t) * vols*2);
    __attribute__((aligned(64))) scs_t* x_s = (scs_t*)malloc( sizeof(scs_t) * vols*2);
//    __attribute__((aligned(64))) scs_t* t_s = (scs_t*)malloc( sizeof(scs_t) * vols*2);
    scs_t* t_s;
    posix_memalign((void **)&t_s, CLS, sizeof(scs_t)*vols*2);

    //fermi_reorder_d2s_dd_(b_s, b);
    //fermi_reorder_d2s_dd_(t_s, b);

    double bnorm;

    fermi_reorder_d2s_dd_renorm_(b_s, b, &bnorm);

    bicgstab_precdd_s_(t_s, b_s, tol, iter, maxiter, nsap, nm);
    prec_s_(x_s, t_s, nsap, nm);

    fermi_reorder_s2d_dd_renorm_(x, x_s, &bnorm);

    //fermi_reorder_s2d_dd_(x, x_s);
    free(b_s);
    free(x_s);
    free(t_s);
    //if (rank==0)printf("%50s %10d\n", "approx_inv_dirac_op: converged iteration",*iter);
  }

  //----------------------------------------------------------------------------------------
  void bicgstab_dd_mix_(scd_t* x, scd_t* b, double *tol, int* conviter, int* maxiter, double *tol_s, int* maxiter_s, int* nsap, int* nm){



    scd_t* q = (scd_t*)malloc( sizeof(scd_t) * vold*2);
    scd_t* r = (scd_t*)malloc( sizeof(scd_t) * vold*2);
    scd_t* p = (scd_t*)malloc( sizeof(scd_t) * vold*2);
    scd_t* t = (scd_t*)malloc( sizeof(scd_t) * vold*2);
    scd_t* u = (scd_t*)malloc( sizeof(scd_t) * vold*2);
    scd_t* r0 = (scd_t*)malloc( sizeof(scd_t) * vold*2);
    double bnorm, rnorm, xnorm, rtmp0, rtmp1, rtmp2, redu[3];
    complex< double > rho0, rho, beta, omega, alpha, ctmp;
    //rvecd_t rvd0, rvd1, rvd2, rvd3, rvd4, rvd5;
    //rvecd_t xr, xi, rr, ri, pr, pi;
    //rvecd_t y0, y1, z0, z1;
    double ar, ai, br, bi;
    int iter, iter_s;
    int mult_count_all;

    mult_count_all = 0;










    
    ddd_d_(q, x);

    rtmp0 = 0;
    rtmp1 = 0;
#pragma omp parallel for reduction(+:rtmp0, rtmp1)
    for(int i=0; i<vold*2; i++){
      for(int j=0; j<24; j++){
        for(int v=0; v<VLEND; v++){
          r[i].ccs[j].v[v]  = b[i].ccs[j].v[v] - q[i].ccs[j].v[v];
          p[i].ccs[j].v[v]  = r[i].ccs[j].v[v];
          r0[i].ccs[j].v[v] = r[i].ccs[j].v[v];
          rtmp0 += b[i].ccs[j].v[v]*b[i].ccs[j].v[v];
          rtmp1 += r[i].ccs[j].v[v]*r[i].ccs[j].v[v];
        }
      }
    }
    redu[0] = rtmp0;
    redu[1] = rtmp1;
#ifdef _MPI_
    MPI_Allreduce(MPI_IN_PLACE,(void *)redu,2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD);
#endif
    bnorm= redu[0];
    rho0 = complex<double>(redu[1],0);

#ifdef _CHECK_
    rnorm = redu[1];
    if (rank==0)printf("@D %5d %24.15e\n", mult_count_all, sqrt(rnorm/bnorm));
#endif

    for (iter=1; iter<(*maxiter);iter++){
      // u = Mp
      // q = Au
      iter_s = mult_count_all;
      approx_inv_dirac_op(u, p, tol_s, &iter_s, maxiter_s, nsap, nm);
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
#ifdef _MPI_
      MPI_Allreduce(MPI_IN_PLACE,(void *)redu,2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD);
#endif
      ctmp = complex<double>(redu[0], redu[1]);
      alpha = rho0 / ctmp;
      if(abs(ctmp)==0)alpha=0;
      ///////////////////////////
      // x = x + alpha u
      // r = r - alpha q
      // rnorm = |r|^2
      ///////////////////////////
      ar=alpha.real();
      ai=alpha.imag();
#pragma omp parallel for
      for(int i=0; i<vold*2; i++){
        for(int j=0; j<12; j++){
          for(int v=0; v<VLEND; v++){
            x[i].cs[j][0].v[v] = x[i].cs[j][0].v[v] + ar *  u[i].cs[j][0].v[v] - ai *  u[i].cs[j][1].v[v];
            x[i].cs[j][1].v[v] = x[i].cs[j][1].v[v] + ar *  u[i].cs[j][1].v[v] + ai *  u[i].cs[j][0].v[v];
            r[i].cs[j][0].v[v] = r[i].cs[j][0].v[v] - ar *  q[i].cs[j][0].v[v] + ai *  q[i].cs[j][1].v[v];
            r[i].cs[j][1].v[v] = r[i].cs[j][1].v[v] - ar *  q[i].cs[j][1].v[v] - ai *  q[i].cs[j][0].v[v];
          }
        }
      }

      // |r|
      rtmp0 = 0;
#pragma omp parallel for reduction(+:rtmp0)
      for(int i=0; i<vold*2; i++){
        for(int j=0; j<24; j++){
          for(int v=0; v<VLEND; v++){
            rtmp0 += r[i].ccs[j].v[v] * r[i].ccs[j].v[v];
          }
        }
      }
      redu[0] = rtmp0;
#ifdef _MPI_
      MPI_Allreduce(MPI_IN_PLACE,(void *)redu,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD);
#endif
      rnorm = redu[0];

      // |x|  for check 
      rtmp0 = 0;
#pragma omp parallel for reduction(+:rtmp0)
      for(int i=0; i<vold*2; i++){
        for(int j=0; j<24; j++){
          for(int v=0; v<VLEND; v++){
            rtmp0 += x[i].ccs[j].v[v] * x[i].ccs[j].v[v];
          }
        }
      }
      redu[0] = rtmp0;
#ifdef _MPI_
      MPI_Allreduce(MPI_IN_PLACE,(void *)redu,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD);
#endif
      xnorm = redu[0];

#ifdef _CHECK_
      if (rank==0)printf("@D %5d %24.15e %24.15e\n", mult_count_all, sqrt(rnorm/bnorm), xnorm);
      check_residual_dd(x,b,&mult_count_all,&bnorm);
#endif
      if (rank==0)printf("iter = %d, rnorm = %24.14e, sqrt(rnorm/bnorm) = %24.14e, xnorm = %24.14e\n", iter, rnorm, sqrt(rnorm/bnorm), xnorm);
      if (sqrt(rnorm/bnorm) < *tol){break;}
      // u = Mr
      // t = Au
      iter_s = mult_count_all;
      approx_inv_dirac_op(u, r, tol_s, &iter_s, maxiter_s, nsap, nm);
      ddd_d_(t, u);
      mult_count_all += iter_s + 1;

      // omega = <t,r> / |t|^2
      rtmp0 = 0;
      rtmp1 = 0;
      rtmp2 = 0;
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
#ifdef _MPI_
      MPI_Allreduce(MPI_IN_PLACE,(void *)redu,3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD);
#endif

      omega = complex<double>(redu[1]/redu[0], redu[2]/redu[0] );
      if(redu[0]==0)omega=0;
      ////////////////////////////
      // x = x + omega u
      // r = r - omega t
      // rnorm = |r|^2
      ////////////////////////////
      ar=omega.real();
      ai=omega.imag();
#pragma omp parallel for
      for(int i=0; i<vold*2; i++){
        for(int j=0; j<12; j++){
          for(int v=0; v<VLEND; v++){
            x[i].cs[j][0].v[v] = x[i].cs[j][0].v[v] + ar *  u[i].cs[j][0].v[v] - ai *  u[i].cs[j][1].v[v];
            x[i].cs[j][1].v[v] = x[i].cs[j][1].v[v] + ar *  u[i].cs[j][1].v[v] + ai *  u[i].cs[j][0].v[v];
            r[i].cs[j][0].v[v] = r[i].cs[j][0].v[v] - ar *  t[i].cs[j][0].v[v] + ai *  t[i].cs[j][1].v[v];
            r[i].cs[j][1].v[v] = r[i].cs[j][1].v[v] - ar *  t[i].cs[j][1].v[v] - ai *  t[i].cs[j][0].v[v];
          }
        }
      }

      // |r|
      rtmp0 = 0;
#pragma omp parallel for reduction(+:rtmp0)
      for(int i=0; i<vold*2; i++){
        for(int j=0; j<24; j++){
          for(int v=0; v<VLEND; v++){
            rtmp0 += r[i].ccs[j].v[v] * r[i].ccs[j].v[v];
          }
        }
      }
      redu[0] = rtmp0;
#ifdef _MPI_
      MPI_Allreduce(MPI_IN_PLACE,(void *)redu,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD);
#endif
      rnorm = redu[0];

      // |x|  for check 
      rtmp0 = 0;
#pragma omp parallel for reduction(+:rtmp0)
      for(int i=0; i<vold*2; i++){
        for(int j=0; j<24; j++){
          for(int v=0; v<VLEND; v++){
            rtmp0 += x[i].ccs[j].v[v] * x[i].ccs[j].v[v];
          }
        }
      }
      redu[0] = rtmp0;
#ifdef _MPI_
      MPI_Allreduce(MPI_IN_PLACE,(void *)redu,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD);
#endif
      xnorm = redu[0];

#ifdef _CHECK_
      if (rank==0)printf("@D %5d %24.15e %24.15e\n", mult_count_all, sqrt(rnorm/bnorm), xnorm);
      check_residual_dd(x,b,&mult_count_all,&bnorm);
#endif
      if (rank==0)printf("iter = %d, rnorm = %24.14e, sqrt(rnorm/bnorm) = %24.14e, xnorm = %24.14e\n", iter, rnorm, sqrt(rnorm/bnorm), xnorm);
      if (sqrt(rnorm/bnorm) < *tol){break;}

      // rho = <r0,r>
      rtmp0 = 0;
      rtmp1 = 0;
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
#ifdef _MPI_
      MPI_Allreduce(MPI_IN_PLACE,(void *)redu,2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD);
#endif
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
      br = beta.real();
      bi = beta.imag();
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
    }//iter
    *conviter = iter;
  }//bicgstab


  void check_residual_dd(scd_t *x, const scd_t *b, const int *iter, const double *bnorm)
  {
    extern int rank, vold;

    __attribute__((aligned(64))) static scd_t *w,*s;

    if(0 == w) w = (scd_t*)malloc( sizeof(scd_t) * vold*2);
    if(0 == s) s = (scd_t*)malloc( sizeof(scd_t) * vold*2);

#ifdef _DEBUG_
    //
    // check residual
    //
    //  w = A x
    //
    ddd_d_(w, x);
    //
    //  s = b - w
    //  err = |s|/|b|
    //
    double rnorm = 0.0e0;
#pragma omp parallel for reduction(+:rnorm)
    for(int i=0; i<vold*2; i++){
      for(int j=0; j<24; j++){
	for(int v=0; v<VLEND; v++){
	  s[i].ccs[j].v[v] = b[i].ccs[j].v[v] - w[i].ccs[j].v[v];
	  rnorm += s[i].ccs[j].v[v] * s[i].ccs[j].v[v];
	}
      }
    }
#ifdef _MPI_
    MPI_Allreduce(MPI_IN_PLACE,(void *)&rnorm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
    double err = sqrt(rnorm/(*bnorm));
    if (0 == rank) printf("CD %5d %24.15e\n",*iter,err);
#endif
  }


#ifdef __cplusplus
}
#endif
