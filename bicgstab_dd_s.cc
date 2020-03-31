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
  extern int rank, vols;
  extern void ddd_s_(scs_t* out, scs_t* in);

#include "timing.h"
  //  extern void check_timing_ (const char *);

  void bicgstab_dd_s_(scs_t* x, scs_t* b, int* conviter, int* maxiter){
    __attribute__((aligned(64))) static scs_t *q, *r, *p , *t, *r0;
    if( q==0) q = (scs_t*)malloc( sizeof(scs_t) * vols*2);
    if( r==0) r = (scs_t*)malloc( sizeof(scs_t) * vols*2);
    if( p==0) p = (scs_t*)malloc( sizeof(scs_t) * vols*2);
    if( t==0) t = (scs_t*)malloc( sizeof(scs_t) * vols*2);
    if(r0==0) r0= (scs_t*)malloc( sizeof(scs_t) * vols*2);
    float bnorm, rnorm, rtmp0, rtmp1, rtmp2, redu[3];
    complex< float > rho0, rho, beta, omega, alpha, ctmp;
    //rvecs_t rvd0, rvd1, rvd2, rvd3, rvd4, rvd5;
    //rvecs_t xr, xi, rr, ri, pr, pi;
    float ar, ai, br, bi, cr, ci;
    int i, j, iter;
    
    //    _DDD_S_TIC_;
    ddd_s_(q, x);
    //    _DDD_S_TOC_;

    rtmp0 = 0;
    rtmp1 = 0;
#pragma omp parallel for private(i, j) reduction(+:rtmp0, rtmp1)
    for(i=0; i<vols*2; i++){
      for(j=0; j<24; j++){
	for(int v=0; v<VLENS; v++){
	  r[i].ccs[j].v[v]  = b[i].ccs[j].v[v] - q[i].ccs[j].v[v];
	  p[i].ccs[j].v[v]  = r[i].ccs[j].v[v];
	  r0[i].ccs[j].v[v] = r[i].ccs[j].v[v];
	  //rvd0 = fmul_s(b[i].ccs[j],b[i].ccs[j]);
	  //rvd1 = fmul_s(r[i].ccs[j],r[i].ccs[j]);
	  //rtmp0 += fsum_s(rvd0);
	  //rtmp1 += fsum_s(rvd1);
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

    for (iter=1; iter<(*maxiter);iter++){
      _BCG_DDS_ITER_TIC_;

      // q = Ap
      //    _DDD_S_TIC_;
      ddd_s_(q, p);
      //    _DDD_S_TOC_;

      // alpha = rho0 / <r0,q>
      rtmp0 = 0;
      rtmp1 = 0;
#pragma omp parallel for private(i, j) reduction(+:rtmp0, rtmp1)
      for(i=0; i<vols*2; i++){
	for(j=0; j<12; j++){
	  for(int v=0; v<VLENS; v++){
	    //rvd0 = fmadd_s(r0[i].cs[j][0], q[i].cs[j][0], fmul_s(r0[i].cs[j][1], q[i].cs[j][1]));
	    //rvd1 = fmsub_s(r0[i].cs[j][0], q[i].cs[j][1], fmul_s(r0[i].cs[j][1], q[i].cs[j][0]));
	    //rtmp0 += fsum_s(rvd0);
	    //rtmp1 += fsum_s(rvd1);
	    rtmp0 += r0[i].cs[j][0].v[v] * q[i].cs[j][0].v[v] + r0[i].cs[j][1].v[v] * q[i].cs[j][1].v[v];
	    rtmp1 += r0[i].cs[j][0].v[v] * q[i].cs[j][1].v[v] - r0[i].cs[j][1].v[v] * q[i].cs[j][0].v[v];
	  }
	}
      }
      redu[0] = rtmp0;
      redu[1] = rtmp1;
#ifdef _MPI_
      MPI_Allreduce(MPI_IN_PLACE,(void *)redu,2,MPI_REAL,MPI_SUM,MPI_COMM_WORLD);
#endif
      ctmp = complex<float>(redu[0], redu[1]);
      alpha = rho0 / ctmp;

      // x = x + alpha p
      // r = r - alpha q
      //rvd0 = fload1_s(alpha.real());
      //rvd1 = fload1_s(alpha.imag());
      ar=alpha.real();
      ai=alpha.imag();
#pragma omp parallel for private(i, j)
      for(i=0; i<vols*2; i++){
	for(j=0; j<12; j++){
	  for(int v=0; v<VLENS; v++){
	    //xr = fcopy_s(x[i].cs[j][0]);
	    //xi = fcopy_s(x[i].cs[j][1]);
	    //rr = fcopy_s(r[i].cs[j][0]);
	    //ri = fcopy_s(r[i].cs[j][1]);
	    //x[i].cs[j][0] =  fmadd_s(rvd0, p[i].cs[j][0], fnmadd_s(rvd1, p[i].cs[j][1], xr));
	    //x[i].cs[j][1] =  fmadd_s(rvd0, p[i].cs[j][1],  fmadd_s(rvd1, p[i].cs[j][0], xi));
	    //r[i].cs[j][0] = fnmadd_s(rvd0, q[i].cs[j][0],  fmadd_s(rvd1, q[i].cs[j][1], rr));
	    //r[i].cs[j][1] = fnmadd_s(rvd0, q[i].cs[j][1], fnmadd_s(rvd1, q[i].cs[j][0], ri));
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
      for(i=0; i<vols*2; i++){
	for(j=0; j<24; j++){
	  for(int v=0; v<VLENS; v++){
	    //rvd0 = fmul_s(r[i].ccs[j],r[i].ccs[j]);
	    //rtmp0 += fsum_s(rvd0);
	    rtmp0 += r[i].ccs[j].v[v] * r[i].ccs[j].v[v];
	  }
	}
      }
      redu[0] = rtmp0;
#ifdef _MPI_
      MPI_Allreduce(MPI_IN_PLACE,(void *)redu,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD);
#endif
      rnorm = redu[0];

      // Check
      //if(rank==0)printf("iter = %d, rnorm = %24.14e, sqrt(rnorm/bnorm) = %24.14e\n", iter, rnorm, sqrt(rnorm/bnorm));
      if (sqrt(rnorm/bnorm) < 1e-6){
	_BCG_DDS_ITER_TOC_;
	break;
      }

      // t = Ar
      //      _DDD_S_TIC_;
      ddd_s_(t, r);
      //      _DDD_S_TOC_;

      // omega = <t,r> / |t|^2
      rtmp0 = 0;
      rtmp1 = 0;
      rtmp2 = 0;
#pragma omp parallel for private(i, j) reduction(+:rtmp0, rtmp1, rtmp2)
      for(i=0; i<vols*2; i++){
	for(j=0; j<12; j++){
	  for(int v=0; v<VLENS; v++){
	    //rvd0 = fmadd_s(t[i].cs[j][0], t[i].cs[j][0], fmul_s(t[i].cs[j][1], t[i].cs[j][1]));
	    //rvd1 = fmadd_s(t[i].cs[j][0], r[i].cs[j][0], fmul_s(t[i].cs[j][1], r[i].cs[j][1]));
	    //rvd2 = fmsub_s(t[i].cs[j][0], r[i].cs[j][1], fmul_s(t[i].cs[j][1], r[i].cs[j][0]));
	    //rtmp0 += fsum_s(rvd0);
	    //rtmp1 += fsum_s(rvd1);
	    //rtmp2 += fsum_s(rvd2);
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
      MPI_Allreduce(MPI_IN_PLACE,(void *)redu,3,MPI_REAL,MPI_SUM,MPI_COMM_WORLD);
#endif
      omega = complex<float>(redu[1]/redu[0], redu[2]/redu[0] );

      // x = x + omega r
      // r = r - omega t
//      rvd0 = fload1_s(omega.real());
//      rvd1 = fload1_s(omega.imag());
//      //#pragma omp parallel for private(i, j, xr, xi, rr, ri)
//      for(i=0; i<vols*2; i++){
//	for(j=0; j<12; j++){
//	  xr = fcopy_s(x[i].cs[j][0]);
//	  xi = fcopy_s(x[i].cs[j][1]);
//	  rr = fcopy_s(r[i].cs[j][0]);
//	  ri = fcopy_s(r[i].cs[j][1]);
//	  x[i].cs[j][0] =  fmadd_s(rvd0, rr, fnmadd_s(rvd1, ri, xr));
//	  x[i].cs[j][1] =  fmadd_s(rvd0, ri,  fmadd_s(rvd1, rr, xi));
//	  r[i].cs[j][0] = fnmadd_s(rvd0, t[i].cs[j][0],  fmadd_s(rvd1, t[i].cs[j][1], rr));
//	  r[i].cs[j][1] = fnmadd_s(rvd0, t[i].cs[j][1], fnmadd_s(rvd1, t[i].cs[j][0], ri));
//	}
//      }
      ar=omega.real();
      ai=omega.imag();
#pragma omp parallel for private(i, j)
      for(i=0; i<vols*2; i++){
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
      rtmp0 = 0;
#pragma omp parallel for private(i, j) reduction(+:rtmp0)
      for(i=0; i<vols*2; i++){
	for(j=0; j<24; j++){
	  for(int v=0; v<VLENS; v++){
	    rtmp0 += r[i].ccs[j].v[v] * r[i].ccs[j].v[v];
	  }
	}
      }
      redu[0] = rtmp0;
#ifdef _MPI_
      MPI_Allreduce(MPI_IN_PLACE,(void *)redu,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD);
#endif
      rnorm = redu[0];

      // Check
      //if(rank==0)printf("iter = %d, rnorm = %24.14e, sqrt(rnorm/bnorm) = %24.14e\n", iter, rnorm, sqrt(rnorm/bnorm));
      if (sqrt(rnorm/bnorm) < 1e-6){
	_BCG_DDS_ITER_TOC_;
	break;
      }


      // rho = <r0,r>
      rtmp0 = 0;
      rtmp1 = 0;
#pragma omp parallel for private(i, j) reduction(+:rtmp0, rtmp1)
      for(i=0; i<vols*2; i++){
	for(j=0; j<12; j++){
	  for(int v=0; v<VLENS; v++){
	    //rvd0 = fmadd_s(r0[i].cs[j][0], r[i].cs[j][0], fmul_s(r0[i].cs[j][1], r[i].cs[j][1]));
	    //rvd1 = fmsub_s(r0[i].cs[j][0], r[i].cs[j][1], fmul_s(r0[i].cs[j][1], r[i].cs[j][0]));
	    //rtmp0 += fsum_s(rvd0);
	    //rtmp1 += fsum_s(rvd1);
	    rtmp0 += r0[i].cs[j][0].v[v] * r[i].cs[j][0].v[v] + r0[i].cs[j][1].v[v] * r[i].cs[j][1].v[v];
	    rtmp1 += r0[i].cs[j][0].v[v] * r[i].cs[j][1].v[v] - r0[i].cs[j][1].v[v] * r[i].cs[j][0].v[v];
	  }
	}
      }
      redu[0] = rtmp0;
      redu[1] = rtmp1;
#ifdef _MPI_
      MPI_Allreduce(MPI_IN_PLACE,(void *)redu,2,MPI_REAL,MPI_SUM,MPI_COMM_WORLD);
#endif
      rho = complex<float>(redu[0],redu[1]);

      beta = alpha*rho/( rho0 * omega);
      rho0 = rho;


      // p = p - omega q
      // p = r + beta p
      //rvd0 = fload1_s(omega.real());
      //rvd1 = fload1_s(omega.imag());
      //rvd2 = fload1_s(beta.real());
      //rvd3 = fload1_s(beta.imag());
      ar = omega.real();
      ai = omega.imag();
      br = beta.real();
      bi = beta.imag();
#pragma omp parallel for private(i, j, cr, ci)
      for(i=0; i<vols*2; i++){
	for(j=0; j<12; j++){
	  for(int v=0; v<VLENS; v++){
	    cr = p[i].cs[j][0].v[v] - ar *  q[i].cs[j][0].v[v] + ai *  q[i].cs[j][1].v[v];
	    ci = p[i].cs[j][1].v[v] - ar *  q[i].cs[j][1].v[v] - ai *  q[i].cs[j][0].v[v];
	    p[i].cs[j][0].v[v] = br * cr - bi * ci + r[i].cs[j][0].v[v];
	    p[i].cs[j][1].v[v] = br * ci + bi * cr + r[i].cs[j][1].v[v];
	  }
	}
      }
      
      _BCG_DDS_ITER_TOC_;
    }//iter
    *conviter = iter;
  }//bicgstab


#ifdef __cplusplus
}
#endif
