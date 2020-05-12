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
#include "eml_lib.h"
#include "util.hh"
#include "addressing.hh"
#include "prefetch.h"
#include <omp.h>



#undef _DEBUG_
//#define _DEBUG_
#undef _CHECK_
//#define _CHECK_
#ifdef __cplusplus
extern "C"{
#endif
  using std::complex;
  extern int vols;
  extern int volse;
#pragma omp threadprivate(volse)
  extern void prec_ddd_s_(scs_t* out, const scs_t* in, const int* nsap, const int* nm);
  void check_residual_s_(const scs_t *x, const scs_t *b, const int *iter, const float *bnorm, const int *nsap, const int *nm);
  extern pglus_t glus;
  extern pclvs_t clvs;

#include "timing.h"
  //  extern void check_timing_ (const char *);

#ifdef _POWER_API_
  extern void power_measure_power();
#endif

  void bicgstab_precdd_s_(scs_t* x, scs_t* b, double* tol, int* conviter, int* maxiter, int* nsap, int* nm){
    __attribute__((aligned(64))) static scs_t *q, *r, *p , *t, *r0;
    if(q==0) {
      void *pp;
      posix_memalign(&pp, CLS, sizeof(scs_t)*vols*2);
      q = (scs_t*)pp;
    }
    if(r==0) {
      void *pp;
      posix_memalign(&pp, CLS, sizeof(scs_t)*vols*2);
      r = (scs_t*)pp;
    }
    if(p==0) {
      void *pp;
      posix_memalign(&pp, CLS, sizeof(scs_t)*vols*2);
      p = (scs_t*)pp;
    }
    if(t==0) {
      void *pp;
      posix_memalign(&pp, CLS, sizeof(scs_t)*vols*2);
      t = (scs_t*)pp;
    }
    if(r0==0) {
      void *pp;
      posix_memalign(&pp, CLS, sizeof(scs_t)*vols*2);
      r0 = (scs_t*)pp;
    }
    float bnorm, rnorm, rtmp0, rtmp1, rtmp2, redu[3];
    complex< float > rho0, rho, beta, omega, alpha, ctmp;
    //rvecs_t rvd0, rvd1, rvd2, rvd3, rvd4, rvd5;
    //rvecs_t xr, xi, rr, ri, pr, pi;
    float ar, ai, br, bi ,cr, ci;
    int i, j, iter, mult_count;
    extern int rank;

    mult_count = 0;


    _BCG_PRECDDS_TIC_;
#ifdef _DEBUG_
    if (0==rank) printf("%s : tol=%24.15e\n",__func__,*tol);
#endif
    static int tid, nthrd;
#pragma omp threadprivate(tid, nthrd)
#pragma omp parallel
    {
#ifdef _OPENMP
      tid = omp_get_thread_num();
      nthrd = omp_get_num_threads();
#else
      tid = 0;
      nthrd = 1;
#endif
    }
#ifdef COMPILE_TIME_DIM_SIZE
    const int vols = VOLS;
#else
    const int vols = ::vols;
#endif

#if 0
    //------------------------------------------------------------------ flexible start
    //    _PREC_DDD_S_TIC_;
    prec_ddd_s_(q, x, nsap, nm);
    //    _PREC_DDD_S_TOC_;
    rtmp0 = 0;
    rtmp1 = 0;
#pragma omp parallel for private(i, j) reduction(+:rtmp0, rtmp1)
    for(i=0; i<vols*2; i++){
      for(j=0; j<24; j++){
        for(int v=0; v<VLENS; v++){
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
    MPI_Allreduce(MPI_IN_PLACE,(void *)redu,2,MPI_REAL,MPI_SUM,MPI_COMM_WORLD);
#endif
    bnorm= redu[0];
    rho0 = complex<float>(redu[1],0);
#else
    //------------------------------------------------------------------ x=0 start
    rtmp0 = (float)0;
#pragma omp parallel for private(i, j) reduction(+:rtmp0)
    for(i=0; i<vols*2; i++){
      for(j=0; j<24; j++){
        for(int v=0; v<VLENS; v++){
          x[i].ccs[j].v[v]  = 0;
          r[i].ccs[j].v[v]  = b[i].ccs[j].v[v];
          p[i].ccs[j].v[v]  = b[i].ccs[j].v[v];
          r0[i].ccs[j].v[v] = b[i].ccs[j].v[v];
          rtmp0 += b[i].ccs[j].v[v]*b[i].ccs[j].v[v];
        }
      }
    }
    redu[0] = rtmp0;
#ifdef _MPI_
    MPI_Allreduce(MPI_IN_PLACE,(void *)redu,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD);
#endif
    bnorm= redu[0];
    rho0 = complex<float>(redu[0],0);
#endif

#ifdef _DEBUG_
    iter = 0;
    mult_count = 0;
#ifdef _CHECK_
    check_residual_s_(x,b,&mult_count,&bnorm,nsap,nm);
#endif
    double err = 1.0e0;
    if (0==rank) printf("@@ %5d %12.7e\n",mult_count,err);
#endif


#ifdef _MPI_
    MPI_Barrier(MPI_COMM_WORLD);
#endif

#ifdef _POWER_API_
    power_measure_power();
#endif //_POWER_API_

    for (iter=0; iter<(*maxiter);iter++){
      _BCG_PRECDDS_ITER_TIC_;
      // q = Ap
      //    _PREC_DDD_S_TIC_;
      prec_ddd_s_(q, p, nsap, nm);
      //    _PREC_DDD_S_TOC_;

      // alpha = rho0 / <r0,q>
      rtmp0 = 0;
      rtmp1 = 0;
#pragma omp parallel reduction(+:rtmp0, rtmp1)
      {
        vecs_t v_rtmp0_0 = fdup_s(0);
        vecs_t v_rtmp0_1 = fdup_s(0);
        vecs_t v_rtmp0_2 = fdup_s(0);
        vecs_t v_rtmp0_3 = fdup_s(0);
        vecs_t v_rtmp1_0 = fdup_s(0);
        vecs_t v_rtmp1_1 = fdup_s(0);
        vecs_t v_rtmp1_2 = fdup_s(0);
        vecs_t v_rtmp1_3 = fdup_s(0);
        pred_t pt = pred_true_all();
        int start, end;
        std::tie(start, end) = static_sched(vols*2, nthrd, tid);
        for(int i=start; i<end; i++){
          __prefetch_inp(r0+i, 2);
          __prefetch_inp(q+i, 2);
#define S(J, K)								\
          v_rtmp0_##K =  fmadd_s(pt, fload1_s(pt, &r0[i], dims_cs, J*4+K, 0), fload1_s(pt, &q[i], dims_cs, J*4+K, 0), v_rtmp0_##K); \
          v_rtmp0_##K =  fmadd_s(pt, fload1_s(pt, &r0[i], dims_cs, J*4+K, 1), fload1_s(pt, &q[i], dims_cs, J*4+K, 1), v_rtmp0_##K); \
          v_rtmp1_##K =  fmadd_s(pt, fload1_s(pt, &r0[i], dims_cs, J*4+K, 0), fload1_s(pt, &q[i], dims_cs, J*4+K, 1), v_rtmp1_##K); \
          v_rtmp1_##K = fnmadd_s(pt, fload1_s(pt, &r0[i], dims_cs, J*4+K, 1), fload1_s(pt, &q[i], dims_cs, J*4+K, 0), v_rtmp1_##K);
          LOOP_3(LOOP_4, S);
#undef S
        }
        v_rtmp0_0 = fadd_s(pt, v_rtmp0_0, v_rtmp0_1);
        v_rtmp0_2 = fadd_s(pt, v_rtmp0_2, v_rtmp0_3);
        v_rtmp0_0 = fadd_s(pt, v_rtmp0_0, v_rtmp0_2);
        v_rtmp1_0 = fadd_s(pt, v_rtmp1_0, v_rtmp1_1);
        v_rtmp1_2 = fadd_s(pt, v_rtmp1_2, v_rtmp1_3);
        v_rtmp1_0 = fadd_s(pt, v_rtmp1_0, v_rtmp1_2);
        rtmp0 = fsum_s(pt, v_rtmp0_0);
        rtmp1 = fsum_s(pt, v_rtmp1_0);
      }
      redu[0] = rtmp0;
      redu[1] = rtmp1;
#ifdef _MPI_
#ifdef _MPI_BARRIER_BEFORE_REDUC_
      _BCG_PRECDDS_ITER_BARRIER_BEFORE_REDUC2_TIC_;
      MPI_Barrier(MPI_COMM_WORLD);
      _BCG_PRECDDS_ITER_BARRIER_BEFORE_REDUC2_TOC_;
#endif
      _BCG_PRECDDS_ITER_REDUC2_TIC_;
      MPI_Allreduce(MPI_IN_PLACE,(void *)redu,2,MPI_REAL,MPI_SUM,MPI_COMM_WORLD);
      _BCG_PRECDDS_ITER_REDUC2_TOC_;
#endif
      ctmp = complex<float>(redu[0], redu[1]);
      alpha = rho0 / ctmp;
      if(abs(ctmp)==0)alpha=0;
      ////////////////////////////
      // x = x + alpha p
      // r = r - alpha q
      // rnorm = |r|^2
      ////////////////////////////
      ar=alpha.real();
      ai=alpha.imag();
#ifdef _CHECK_
      printf("iter = %d, alpha = %24.14e %24.14e, ctmp = %24.14e %24.14e\n", iter, ar, ai, redu[0], redu[1]);
#endif
#pragma omp parallel for private(i, j)
      for(i=0; i<vols*2; i++){
        __prefetch_inp(p+i, 2);
        __prefetch_inp(q+i, 2);
        __prefetch_out(x+i, 2);
        __prefetch_out(r+i, 2);
        for(j=0; j<12; j++){
          for(int v=0; v<VLENS; v++){
            x[i].cs[j][0].v[v] = x[i].cs[j][0].v[v] + ar *  p[i].cs[j][0].v[v] - ai *  p[i].cs[j][1].v[v];
            x[i].cs[j][1].v[v] = x[i].cs[j][1].v[v] + ar *  p[i].cs[j][1].v[v] + ai *  p[i].cs[j][0].v[v];
            r[i].cs[j][0].v[v] = r[i].cs[j][0].v[v] - ar *  q[i].cs[j][0].v[v] + ai *  q[i].cs[j][1].v[v];
            r[i].cs[j][1].v[v] = r[i].cs[j][1].v[v] - ar *  q[i].cs[j][1].v[v] - ai *  q[i].cs[j][0].v[v];
          }
        }
      }

      // |r|
      rtmp0 = 0;
#pragma omp parallel reduction(+:rtmp0)
      {
        vecs_t v_rtmp0_0 = fdup_s(0);
        vecs_t v_rtmp0_1 = fdup_s(0);
        vecs_t v_rtmp0_2 = fdup_s(0);
        vecs_t v_rtmp0_3 = fdup_s(0);
        vecs_t v_rtmp0_4 = fdup_s(0);
        vecs_t v_rtmp0_5 = fdup_s(0);
        vecs_t v_rtmp0_6 = fdup_s(0);
        vecs_t v_rtmp0_7 = fdup_s(0);
        pred_t pt = pred_true_all();
        int start, end;
        std::tie(start, end) = static_sched(vols*2, nthrd, tid);
        for(int i=start; i<end; i++){
          __prefetch_inp(r+i, 2);
#define S(J, K) \
          v_rtmp0_##K = fmadd_s(pt, fload1_s(pt, &r[i], dims_1d, J*8+K), fload1_s(pt, &r[i], dims_1d, J*8+K, 0), v_rtmp0_##K);
          LOOP_3(LOOP_8, S);
#undef S
        }
        v_rtmp0_0 = fadd_s(pt, v_rtmp0_0, v_rtmp0_1);
        v_rtmp0_2 = fadd_s(pt, v_rtmp0_2, v_rtmp0_3);
        v_rtmp0_4 = fadd_s(pt, v_rtmp0_4, v_rtmp0_5);
        v_rtmp0_6 = fadd_s(pt, v_rtmp0_6, v_rtmp0_7);
        v_rtmp0_0 = fadd_s(pt, v_rtmp0_0, v_rtmp0_2);
        v_rtmp0_4 = fadd_s(pt, v_rtmp0_4, v_rtmp0_6);
        v_rtmp0_0 = fadd_s(pt, v_rtmp0_0, v_rtmp0_4);
        rtmp0 = fsum_s(pt, v_rtmp0_0);
      }
      redu[0] = rtmp0;
#ifdef _MPI_
#ifdef _MPI_BARRIER_BEFORE_REDUC_
      _BCG_PRECDDS_ITER_BARRIER_BEFORE_REDUC1_TIC_;
      MPI_Barrier(MPI_COMM_WORLD);
      _BCG_PRECDDS_ITER_BARRIER_BEFORE_REDUC1_TOC_;
#endif
      _BCG_PRECDDS_ITER_REDUC1_TIC_;
      MPI_Allreduce(MPI_IN_PLACE,(void *)redu,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD);
      _BCG_PRECDDS_ITER_REDUC1_TOC_;
#endif
      rnorm = redu[0];

      // Check
      //printf("iter = %d, rnorm = %24.14e, sqrt(rnorm/bnorm) = %24.14e\n", iter, rnorm, sqrt(rnorm/bnorm));
#ifdef _DEBUG_
      mult_count++;
#ifdef _CHECK_
      check_residual_s_(x,b,&mult_count,&bnorm,nsap,nm);
#endif
      err = (double)sqrt(rnorm/bnorm);
      if (0==rank) printf("@@ %5d %12.7e\n",mult_count,err);
#endif
      if (sqrt(rnorm/bnorm) < *tol){
	_BCG_PRECDDS_ITER_TOC_;
	break;
      }

      // t = Ar
      //      _PREC_DDD_S_TIC_;
      prec_ddd_s_(t, r, nsap, nm);
      //      _PREC_DDD_S_TOC_;

      // omega = <t,r> / |t|^2
      rtmp0 = 0;
      rtmp1 = 0;
      rtmp2 = 0;
#pragma omp parallel reduction(+:rtmp0, rtmp1, rtmp2)
      {
        vecs_t v_rtmp0_0 = fdup_s(0);
        vecs_t v_rtmp0_1 = fdup_s(0);
        vecs_t v_rtmp0_2 = fdup_s(0);
        vecs_t v_rtmp1_0 = fdup_s(0);
        vecs_t v_rtmp1_1 = fdup_s(0);
        vecs_t v_rtmp1_2 = fdup_s(0);
        vecs_t v_rtmp2_0 = fdup_s(0);
        vecs_t v_rtmp2_1 = fdup_s(0);
        vecs_t v_rtmp2_2 = fdup_s(0);
        pred_t pt = pred_true_all();
        int start, end;
        std::tie(start, end) = static_sched(vols*2, nthrd, tid);
        for(int i=start; i<end; i++){
          __prefetch_inp(r+i, 2);
          __prefetch_inp(t+i, 2);
#define S(J, K)                                                         \
          v_rtmp0_##K =  fmadd_s(pt, fload1_s(pt, &t[i], dims_cs, J*3+K, 0), fload1_s(pt, &t[i], dims_cs, J*3+K, 0), v_rtmp0_##K); \
          v_rtmp0_##K =  fmadd_s(pt, fload1_s(pt, &t[i], dims_cs, J*3+K, 1), fload1_s(pt, &t[i], dims_cs, J*3+K, 1), v_rtmp0_##K); \
          v_rtmp1_##K =  fmadd_s(pt, fload1_s(pt, &t[i], dims_cs, J*3+K, 0), fload1_s(pt, &r[i], dims_cs, J*3+K, 0), v_rtmp1_##K); \
          v_rtmp1_##K =  fmadd_s(pt, fload1_s(pt, &t[i], dims_cs, J*3+K, 1), fload1_s(pt, &r[i], dims_cs, J*3+K, 1), v_rtmp1_##K); \
          v_rtmp2_##K =  fmadd_s(pt, fload1_s(pt, &t[i], dims_cs, J*3+K, 0), fload1_s(pt, &r[i], dims_cs, J*3+K, 1), v_rtmp2_##K); \
          v_rtmp2_##K = fnmadd_s(pt, fload1_s(pt, &t[i], dims_cs, J*3+K, 1), fload1_s(pt, &r[i], dims_cs, J*3+K, 0), v_rtmp2_##K);
          LOOP_4(LOOP_3, S);
#undef S
        }
        v_rtmp0_0 = fadd_s(pt, v_rtmp0_0, v_rtmp0_1);
        v_rtmp0_0 = fadd_s(pt, v_rtmp0_0, v_rtmp0_2);
        v_rtmp1_0 = fadd_s(pt, v_rtmp1_0, v_rtmp1_1);
        v_rtmp1_0 = fadd_s(pt, v_rtmp1_0, v_rtmp1_2);
        v_rtmp2_0 = fadd_s(pt, v_rtmp2_0, v_rtmp2_1);
        v_rtmp2_0 = fadd_s(pt, v_rtmp2_0, v_rtmp2_2);
        rtmp0 = fsum_s(pt, v_rtmp0_0);
        rtmp1 = fsum_s(pt, v_rtmp1_0);
        rtmp2 = fsum_s(pt, v_rtmp2_0);
      }
      redu[0] = rtmp0;
      redu[1] = rtmp1;
      redu[2] = rtmp2;
#ifdef _MPI_
#ifdef _MPI_BARRIER_BEFORE_REDUC_
      _BCG_PRECDDS_ITER_BARRIER_BEFORE_REDUC3_1_TIC_;
      MPI_Barrier(MPI_COMM_WORLD);
      _BCG_PRECDDS_ITER_BARRIER_BEFORE_REDUC3_1_TOC_;
#endif
      _BCG_PRECDDS_ITER_REDUC3_TIC_;
      MPI_Allreduce(MPI_IN_PLACE,(void *)redu,3,MPI_REAL,MPI_SUM,MPI_COMM_WORLD);
      _BCG_PRECDDS_ITER_REDUC3_TOC_;
#endif
      omega = complex<float>(redu[1]/redu[0], redu[2]/redu[0] );
#ifdef _CHECK_
      printf("iter = %d, omega = %24.14e %24.14e %24.14e\n", iter, redu[0], redu[1], redu[2]);
#endif
      if(redu[0]==0)omega=0;
      /////////////////////////
      // x = x + omega r
      // r = r - omega t
      // rnorm = |r|^2
      /////////////////////////
      ar=omega.real();
      ai=omega.imag();
#pragma omp parallel for private(i, j)
      for(i=0; i<vols*2; i++){
        __prefetch_inp(t+i, 2);
        __prefetch_out(x+i, 2);
        __prefetch_out(r+i, 2);
        for(j=0; j<12; j++){
          for(int v=0; v<VLENS; v++){
            x[i].cs[j][0].v[v] = x[i].cs[j][0].v[v] + ar *  r[i].cs[j][0].v[v] - ai *  r[i].cs[j][1].v[v];
            x[i].cs[j][1].v[v] = x[i].cs[j][1].v[v] + ar *  r[i].cs[j][1].v[v] + ai *  r[i].cs[j][0].v[v];
            r[i].cs[j][0].v[v] = r[i].cs[j][0].v[v] - ar *  t[i].cs[j][0].v[v] + ai *  t[i].cs[j][1].v[v];
            r[i].cs[j][1].v[v] = r[i].cs[j][1].v[v] - ar *  t[i].cs[j][1].v[v] - ai *  t[i].cs[j][0].v[v];
          }
        }
      }

      // |r|
      // rho = <r0,r>
      rtmp0 = 0;
      rtmp1 = 0;
      rtmp2 = 0;
#pragma omp parallel reduction(+:rtmp0, rtmp1, rtmp2)
      {
        vecs_t v_rtmp0_0 = fdup_s(0);
        vecs_t v_rtmp0_1 = fdup_s(0);
        vecs_t v_rtmp0_2 = fdup_s(0);
        vecs_t v_rtmp1_0 = fdup_s(0);
        vecs_t v_rtmp1_1 = fdup_s(0);
        vecs_t v_rtmp1_2 = fdup_s(0);
        vecs_t v_rtmp2_0 = fdup_s(0);
        vecs_t v_rtmp2_1 = fdup_s(0);
        vecs_t v_rtmp2_2 = fdup_s(0);
        pred_t pt = pred_true_all();
        int start, end;
        std::tie(start, end) = static_sched(vols*2, nthrd, tid);
        for(int i=start; i<end; i++){
          __prefetch_inp(r0+i, 2);
          __prefetch_inp(r+i, 2);
#define S(J, K)                                                         \
          v_rtmp0_##K =  fmadd_s(pt, fload1_s(pt, &r[i], dims_cs, J*3+K, 0), fload1_s(pt, &r[i], dims_cs, J*3+K, 0), v_rtmp0_##K); \
          v_rtmp0_##K =  fmadd_s(pt, fload1_s(pt, &r[i], dims_cs, J*3+K, 1), fload1_s(pt, &r[i], dims_cs, J*3+K, 1), v_rtmp0_##K); \
          v_rtmp1_##K =  fmadd_s(pt, fload1_s(pt, &r0[i], dims_cs, J*3+K, 0), fload1_s(pt, &r[i], dims_cs, J*3+K, 0), v_rtmp1_##K); \
          v_rtmp1_##K =  fmadd_s(pt, fload1_s(pt, &r0[i], dims_cs, J*3+K, 1), fload1_s(pt, &r[i], dims_cs, J*3+K, 1), v_rtmp1_##K); \
          v_rtmp2_##K =  fmadd_s(pt, fload1_s(pt, &r0[i], dims_cs, J*3+K, 0), fload1_s(pt, &r[i], dims_cs, J*3+K, 1), v_rtmp2_##K); \
          v_rtmp2_##K = fnmadd_s(pt, fload1_s(pt, &r0[i], dims_cs, J*3+K, 1), fload1_s(pt, &r[i], dims_cs, J*3+K, 0), v_rtmp2_##K);
          LOOP_4(LOOP_3, S);
#undef S
        }
        v_rtmp0_0 = fadd_s(pt, v_rtmp0_0, v_rtmp0_1);
        v_rtmp0_0 = fadd_s(pt, v_rtmp0_0, v_rtmp0_2);
        v_rtmp1_0 = fadd_s(pt, v_rtmp1_0, v_rtmp1_1);
        v_rtmp1_0 = fadd_s(pt, v_rtmp1_0, v_rtmp1_2);
        v_rtmp2_0 = fadd_s(pt, v_rtmp2_0, v_rtmp2_1);
        v_rtmp2_0 = fadd_s(pt, v_rtmp2_0, v_rtmp2_2);
        rtmp0 = fsum_s(pt, v_rtmp0_0);
        rtmp1 = fsum_s(pt, v_rtmp1_0);
        rtmp2 = fsum_s(pt, v_rtmp2_0);
      }
      redu[0] = rtmp0;
      redu[1] = rtmp1;
      redu[2] = rtmp2;
#ifdef _MPI_
#ifdef _MPI_BARRIER_BEFORE_REDUC_
      _BCG_PRECDDS_ITER_BARRIER_BEFORE_REDUC3_2_TIC_;
      MPI_Barrier(MPI_COMM_WORLD);
      _BCG_PRECDDS_ITER_BARRIER_BEFORE_REDUC3_2_TOC_;
#endif
      _BCG_PRECDDS_ITER_REDUC3_TIC_;
      MPI_Allreduce(MPI_IN_PLACE,(void *)redu,3,MPI_REAL,MPI_SUM,MPI_COMM_WORLD);
      _BCG_PRECDDS_ITER_REDUC3_TOC_;
#endif
      rnorm = redu[0];
      rho = complex<float>(redu[1],redu[2]);
      // Check
      //printf("iter = %d, rnorm = %24.14e, sqrt(rnorm/bnorm) = %24.14e\n", iter, rnorm, sqrt(rnorm/bnorm));
#ifdef _DEBUG_
      mult_count++;
#ifdef _CHECK_
      check_residual_s_(x,b,&mult_count,&bnorm,nsap,nm);
#endif
      err = (double)sqrt(rnorm/bnorm);
      if (0==rank) printf("@@ %5d %12.7e\n",mult_count,err);
#endif
      if (sqrt(rnorm/bnorm) < *tol){
	_BCG_PRECDDS_ITER_TOC_;
	break;
      }
#ifdef _CHECK_
      printf("iter = %d, rho = %24.14e %24.14e\n", iter, redu[1], redu[2]);
#endif
      beta = alpha*rho/( rho0 * omega);
      rho0 = rho;
      if(abs(omega)==0||abs(rho0)==0)beta=0;
#ifdef _CHECK_
      printf("iter = %d, beta = %24.14e %24.14e\n", iter, beta.real(), beta.imag());
#endif
      //////////////////////
      // p = p - omega q
      // p = r + beta p
      //////////////////////
      ar = omega.real();
      ai = omega.imag();
      br = beta.real();
      bi = beta.imag();
#pragma omp parallel for private(i, j, cr, ci)
      for(i=0; i<vols*2; i++){
        __prefetch_inp(q+i, 2);
        __prefetch_inp(r+i, 2);
        __prefetch_out(p+i, 2);
        for(j=0; j<12; j++){
          for(int v=0; v<VLENS; v++){
            cr = p[i].cs[j][0].v[v] - ar *  q[i].cs[j][0].v[v] + ai *  q[i].cs[j][1].v[v];
            ci = p[i].cs[j][1].v[v] - ar *  q[i].cs[j][1].v[v] - ai *  q[i].cs[j][0].v[v];
            p[i].cs[j][0].v[v] = br * cr - bi * ci + r[i].cs[j][0].v[v];
            p[i].cs[j][1].v[v] = br * ci + bi * cr + r[i].cs[j][1].v[v];
          }
        }
      }
      _BCG_PRECDDS_ITER_TOC_;
    }//iter

#ifdef _POWER_API_
    power_measure_power();
#endif //_POWER_API_



    *conviter = mult_count;
    //    *conviter = iter;
    _BCG_PRECDDS_TOC_;
#ifdef _DEBUG_
    check_residual_s_(x,b,&mult_count,&bnorm,nsap,nm);
#endif
  }//bicgstab

  void check_residual_s_(const scs_t *x, const scs_t *b, const int *iter, const float *bnorm, const int *nsap, const int *nm){
    extern int rank;
    __attribute__((aligned(64))) static scs_t *w,*s;
    
    if(0 == w) w = (scs_t*)malloc( sizeof(scs_t) * vols*2);
    if(0 == s) s = (scs_t*)malloc( sizeof(scs_t) * vols*2);
    
    //
    // check residual
    //
    // w = A x  (A=DMsap)
    //
    prec_ddd_s_(w, x, nsap, nm);
    
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
#ifdef _MPI_
    MPI_Allreduce(MPI_IN_PLACE,(void *)&rnorm,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD);
#endif
    float err = (double)sqrt(rnorm/(*bnorm));
    if (0==rank) printf("@C %5d %12.7e\n",*iter,err);
  }


#ifdef __cplusplus
}
#endif
