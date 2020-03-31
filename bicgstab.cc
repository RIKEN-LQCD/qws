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


#ifdef __cplusplus
extern "C"{
#endif

  using std::complex;
  extern int vold;
  extern void mtilde_vm_(scd_t* out, const scd_t* in);
  void check_residual(const scd_t *x,const scd_t *b,const int *iter,const double *bnorm);

#include "timing.h"
  //  extern void check_timing_ (const char *);

  void bicgstab_vm_(scd_t* x, scd_t* b, int* conviter, int* maxiter){
    __attribute__((aligned(64))) scd_t* q = (scd_t*)malloc( sizeof(scd_t) * vold);
    __attribute__((aligned(64))) scd_t* r = (scd_t*)malloc( sizeof(scd_t) * vold);
    __attribute__((aligned(64))) scd_t* p = (scd_t*)malloc( sizeof(scd_t) * vold);
    __attribute__((aligned(64))) scd_t* t = (scd_t*)malloc( sizeof(scd_t) * vold);
    __attribute__((aligned(64))) scd_t* r0 = (scd_t*)malloc( sizeof(scd_t) * vold);
    double bnorm, rnorm, rtmp0, rtmp1, rtmp2, redu[3];
    complex< double > rho0, rho, beta, omega, alpha, ctmp;
    //rvecd_t rvd0, rvd1, rvd2, rvd3, rvd4, rvd5;
    //rvecd_t xr, xi, rr, ri, pr, pi;
    double ar, ai, br, bi ,cr, ci;
    int i, j, iter;
    extern int rank;
    
    //    _MTILDE_TIC_;
    mtilde_vm_(q, x);
    //    _MTILDE_TOC_;

    rtmp0 = 0;
    rtmp1 = 0;
#pragma omp parallel for private(i, j) reduction(+:rtmp0, rtmp1)
    for(i=0; i<vold; i++){
      for(j=0; j<24; j++){
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
    MPI_Allreduce(MPI_IN_PLACE,(void *)redu,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
    bnorm= redu[0];
    rho0 = complex<double>(redu[1],0);

#ifdef _DEBUG_
    iter = 0;
    double err = 1.0e0;
    check_residual(x,b,&iter,&bnorm);
    if (0 == rank) printf("%% %5d %24.15e\n",iter,err);
#endif


    for (iter=1; iter<(*maxiter);iter++){
      _BCG_ITER_TIC_;

      // q = Ap
      //    _MTILDE_TIC_;
      mtilde_vm_(q, p);
      //    _MTILDE_TOC_;

      // alpha = rho0 / <r0,q>
      rtmp0 = 0;
      rtmp1 = 0;
#pragma omp parallel for private(i, j) reduction(+:rtmp0, rtmp1)
      for(i=0; i<vold; i++){
        for(j=0; j<12; j++){
          for(int v=0; v<VLEND; v++){
            rtmp0 += r0[i].cs[j][0].v[v] * q[i].cs[j][0].v[v] + r0[i].cs[j][1].v[v] * q[i].cs[j][1].v[v];
            rtmp1 += r0[i].cs[j][0].v[v] * q[i].cs[j][1].v[v] - r0[i].cs[j][1].v[v] * q[i].cs[j][0].v[v];
          }
        }
      }
      redu[0] = rtmp0;
      redu[1] = rtmp1;
#ifdef _MPI_
      MPI_Allreduce(MPI_IN_PLACE,(void *)redu,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
      ctmp = complex<double>(redu[0], redu[1]);
      alpha = rho0 / ctmp;

      // x = x + alpha p
      // r = r - alpha q
      ar=alpha.real();
      ai=alpha.imag();
#pragma omp parallel for private(i, j)
      for(i=0; i<vold; i++){
        for(j=0; j<12; j++){
          for(int v=0; v<VLEND; v++){
            x[i].cs[j][0].v[v] = x[i].cs[j][0].v[v] + ar *  p[i].cs[j][0].v[v] - ai *  p[i].cs[j][1].v[v];
            x[i].cs[j][1].v[v] = x[i].cs[j][1].v[v] + ar *  p[i].cs[j][1].v[v] + ai *  p[i].cs[j][0].v[v];
            r[i].cs[j][0].v[v] = r[i].cs[j][0].v[v] - ar *  q[i].cs[j][0].v[v] + ai *  q[i].cs[j][1].v[v];
            r[i].cs[j][1].v[v] = r[i].cs[j][1].v[v] - ar *  q[i].cs[j][1].v[v] - ai *  q[i].cs[j][0].v[v];
          }
        }
      }

      // |r|
      rtmp0 = 0;
#pragma omp parallel for private(i, j) reduction(+:rtmp0)
      for(i=0; i<vold; i++){
        for(j=0; j<24; j++){
          for(int v=0; v<VLEND; v++){
            rtmp0 += r[i].ccs[j].v[v] * r[i].ccs[j].v[v];
          }
        }
      }
      redu[0] = rtmp0;
#ifdef _MPI_
      MPI_Allreduce(MPI_IN_PLACE,(void *)redu,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
      rnorm = redu[0];

      // Check
      //printf("iter = %d, rnorm = %24.14e, sqrt(rnorm/bnorm) = %24.14e\n", iter, rnorm, sqrt(rnorm/bnorm));
#ifdef _DEBUG_
      check_residual(x,b,&iter,&bnorm);
      err = sqrt(rnorm/bnorm);
      if (0 == rank) printf("%% %5d %24.15e\n",iter,err);
#endif
      if (sqrt(rnorm/bnorm) < 1e-14){
	_BCG_ITER_TOC_;
	break;
      }

      // t = Ar
      //      _MTILDE_TIC_;
      mtilde_vm_(t, r);
      //      _MTILDE_TOC_;

      // omega = <t,r> / |t|^2
      rtmp0 = 0;
      rtmp1 = 0;
      rtmp2 = 0;
#pragma omp parallel for private(i, j) reduction(+:rtmp0, rtmp1, rtmp2)
      for(i=0; i<vold; i++){
        for(j=0; j<12; j++){
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
      MPI_Allreduce(MPI_IN_PLACE,(void *)redu,3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
      omega = complex<double>(redu[1]/redu[0], redu[2]/redu[0] );

      // x = x + omega r
      // r = r - omega t
      ar=omega.real();
      ai=omega.imag();
#pragma omp parallel for private(i, j)
      for(i=0; i<vold; i++){
        for(j=0; j<12; j++){
          for(int v=0; v<VLEND; v++){
            x[i].cs[j][0].v[v] = x[i].cs[j][0].v[v] + ar *  r[i].cs[j][0].v[v] - ai *  r[i].cs[j][1].v[v];
            x[i].cs[j][1].v[v] = x[i].cs[j][1].v[v] + ar *  r[i].cs[j][1].v[v] + ai *  r[i].cs[j][0].v[v];
            r[i].cs[j][0].v[v] = r[i].cs[j][0].v[v] - ar *  t[i].cs[j][0].v[v] + ai *  t[i].cs[j][1].v[v];
            r[i].cs[j][1].v[v] = r[i].cs[j][1].v[v] - ar *  t[i].cs[j][1].v[v] - ai *  t[i].cs[j][0].v[v];
          }
        }
      }

      // |r|      
      rtmp0 = 0;
#pragma omp parallel for private(i, j) reduction(+:rtmp0)
      for(i=0; i<vold; i++){
        for(j=0; j<24; j++){
          for(int v=0; v<VLEND; v++){
            rtmp0 += r[i].ccs[j].v[v] * r[i].ccs[j].v[v];
          }
        }
      }
      redu[0] = rtmp0;
#ifdef _MPI_
      MPI_Allreduce(MPI_IN_PLACE,(void *)redu,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
      rnorm = redu[0];

      // Check
      //printf("iter = %d, rnorm = %24.14e, sqrt(rnorm/bnorm) = %24.14e\n", iter, rnorm, sqrt(rnorm/bnorm));
#ifdef _DEBUG_
      check_residual(x,b,&iter,&bnorm);
      err = sqrt(rnorm/bnorm);
      if (0 == rank) printf("%% %5d %24.15e\n",iter,err);
#endif
      if (sqrt(rnorm/bnorm) < 1e-14){
	_BCG_ITER_TOC_;
	break;
      }

      // rho = <r0,r>
      rtmp0 = 0;
      rtmp1 = 0;
#pragma omp parallel for private(i, j) reduction(+:rtmp0, rtmp1)
      for(i=0; i<vold; i++){
        for(j=0; j<12; j++){
          for(int v=0; v<VLEND; v++){
            rtmp0 += r0[i].cs[j][0].v[v] * r[i].cs[j][0].v[v] + r0[i].cs[j][1].v[v] * r[i].cs[j][1].v[v];
            rtmp1 += r0[i].cs[j][0].v[v] * r[i].cs[j][1].v[v] - r0[i].cs[j][1].v[v] * r[i].cs[j][0].v[v];
          }
        }
      }
      redu[0] = rtmp0;
      redu[1] = rtmp1;
#ifdef _MPI_
      MPI_Allreduce(MPI_IN_PLACE,(void *)redu,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
      rho = complex<double>(redu[0],redu[1]);

      beta = alpha*rho/( rho0 * omega);
      rho0 = rho;


      // p = p - omega q
      // p = r + beta p
      ar = omega.real();
      ai = omega.imag();
      br = beta.real();
      bi = beta.imag();
#pragma omp parallel for private(i, j, cr, ci)
      for(i=0; i<vold; i++){
        for(j=0; j<12; j++){
          for(int v=0; v<VLEND; v++){
            cr = p[i].cs[j][0].v[v] - ar *  q[i].cs[j][0].v[v] + ai *  q[i].cs[j][1].v[v];
            ci = p[i].cs[j][1].v[v] - ar *  q[i].cs[j][1].v[v] - ai *  q[i].cs[j][0].v[v];
            p[i].cs[j][0].v[v] = br * cr - bi * ci + r[i].cs[j][0].v[v];
            p[i].cs[j][1].v[v] = br * ci + bi * cr + r[i].cs[j][1].v[v];
          }
        }
      }
      _BCG_ITER_TOC_;
    }//iter


  }//bicgstab

  void check_residual(const scd_t *x, const scd_t *b, const int *iter, const double *bnorm)
  {
    extern int rank;
    __attribute__((aligned(64))) static scd_t *w,*s;

    if(0 == w) w = (scd_t*)malloc( sizeof(scd_t) * vold);
    if(0 == s) s = (scd_t*)malloc( sizeof(scd_t) * vold);

#ifdef _DEBUG_
    extern int rank;
    //
    // check residual
    //
    //  w = A x
    //
    mtilde_vm_(w, x);
    //
    //  s = b - w
    //  err = |s|/|b|
    //
    double rnorm = 0.0e0;
#pragma omp parallel for reduction(+:rnorm)
    for(int i=0; i<vold; i++){
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
    if (0 == rank) printf("C %5d %24.15e\n",*iter,err);
#endif
  }


#ifdef __cplusplus
}
#endif
