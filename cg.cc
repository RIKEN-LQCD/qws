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
#include <math.h>


#ifdef __cplusplus
extern "C"{
#endif

  extern int vold;
  extern void mtdagmt_vm_(scd_t* out, scd_t* in);

#include "timing.h"
  //  extern void check_timing_ (const char *);

  void cg_vm_(scd_t* x, scd_t* b, int* conviter, int* maxiter){
    scd_t* q = (scd_t*)malloc( sizeof(scd_t) * vold);
    scd_t* r = (scd_t*)malloc( sizeof(scd_t) * vold);
    scd_t* p = (scd_t*)malloc( sizeof(scd_t) * vold);
    double bnorm, rnorm, rtmp0, rtmp1, redu[3];
    double alpha, beta, rho, rho0;
    rvecd_t rvd0, rvd1;
    rvecd_t xx, rr;
    int i, j, iter;
    
    mtdagmt_vm_(q, x);
    rtmp0 = 0;
    rtmp1 = 0;
#pragma omp parallel for private(j, rvd0, rvd1) reduction(+:rtmp0, rtmp1)
    for(i=0; i<vold; i++){
      for(j=0; j<24; j++){
	r[i].ccs[j] = fsub_d(b[i].ccs[j], q[i].ccs[j]);
	p[i].ccs[j] = fcopy_d(r[i].ccs[j]);
	rvd0 = fmul_d(b[i].ccs[j],b[i].ccs[j]);
	rvd1 = fmul_d(r[i].ccs[j],r[i].ccs[j]);
	rtmp0 += fsum_d(rvd0);
	rtmp1 += fsum_d(rvd1);
      }
    }
    redu[0] = rtmp0;
    redu[1] = rtmp1;
#ifdef _MPI_
    MPI_Allreduce(MPI_IN_PLACE,(void *)redu,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
    bnorm= redu[0];
    rho0 = redu[1];

    for (iter=1; iter<(*maxiter);iter++){
      _CG_ITER_TIC_;
      
      // q = Ap
      mtdagmt_vm_(q, p);

      // alpha = rho0 / <p,q>
      rtmp0 = 0;
#pragma omp parallel for private(j, rvd0) reduction(+:rtmp0)
      for(i=0; i<vold; i++){
	for(j=0; j<12; j++){
	  rvd0 = fmadd_d(p[i].cs[j][0], q[i].cs[j][0], fmul_d(p[i].cs[j][1], q[i].cs[j][1]));
	  rtmp0 += fsum_d(rvd0);
	}
      }
      redu[0] = rtmp0;
#ifdef _MPI_
      MPI_Allreduce(MPI_IN_PLACE,(void *)redu,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
      alpha = rho0 / redu[0];

      // x = x + alpha p
      // r = r - alpha q
      rvd0 = fload1_d(alpha);
#pragma omp parallel for private(j, xx, rr)
      for(i=0; i<vold; i++){
	for(j=0; j<24; j++){
	  xx = fcopy_d(x[i].ccs[j]);
	  rr = fcopy_d(r[i].ccs[j]);
	  x[i].ccs[j] =  fmadd_d(rvd0, p[i].ccs[j], xx);
	  r[i].ccs[j] = fnmadd_d(rvd0, q[i].ccs[j], rr);
	}
      }

      // |r|      
      rtmp0 = 0;
#pragma omp parallel for private(j, rvd0) reduction(+:rtmp0)
      for(i=0; i<vold; i++){
	for(j=0; j<24; j++){
	  rvd0 = fmul_d(r[i].ccs[j],r[i].ccs[j]);
	  rtmp0 += fsum_d(rvd0);
	}
      }
      redu[0] = rtmp0;
#ifdef _MPI_
      MPI_Allreduce(MPI_IN_PLACE,(void *)redu,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
      rnorm = redu[0];

      // Check
      //printf("iter = %d, rnorm = %24.14e, sqrt(rnorm/bnorm) = %24.14e\n", iter, rnorm, sqrt(rnorm/bnorm));
      if (sqrt(rnorm/bnorm) < 1e-14){
	_CG_ITER_TOC_;
	break;
      }

      // rho = <r0,r>
      rho = rnorm;

      beta = rho/rho0;
      rho0 = rho;

      // p = r + beta p
      rvd0 = fload1_d(beta);
#pragma omp parallel for private(j, rr)
      for(i=0; i<vold; i++){
	for(j=0; j<24; j++){
	  rr =  fcopy_d(r[i].ccs[j]);
	  p[i].ccs[j] =  fmadd_d(rvd0, p[i].ccs[j], rr);
	}
      }
      
      _CG_ITER_TOC_;
    }//iter

  }//cg


#ifdef __cplusplus
}
#endif
