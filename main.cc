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
/*
!===============================================================================
!
! Copyright (C) 2014,2015 Yoshifumi Nakamura
!
! _CODENAME_ = ???
!
! This file is part of _CODENAME_
!
! _CODENAME_ is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! _CODENAME_ is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with _CODENAME_.  If not, see <http://www.gnu.org/licenses/>.
!
!-----------------------------------------------------------------------------
*/
#define MAIN

#ifdef _MPI_
#include <mpi.h>
#endif

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "qws.h"
#ifdef HALF_PREC
#include "qws_h.h"
#endif
#include "qwsintrin.h"
#include "timing.h"
#include <random>
static std::mt19937 mt;
static std::uniform_real_distribution<double> rand01(0.0,1.0);

extern int vold, vols, rank, nx, ny, nz, nt, nxh, nxd, nxs;
extern double kappa2, kappa, mkappa;

extern "C" void qws_init_(int* lx,  int* ly, int* lz, int* lt, int* npe_f, int* fbc_f, int* pce_f, int* pco_f, int* block_size);
extern "C" void qws_finalize_();
extern "C" void bicgstab_vm_(scd_t* x, scd_t* b, int* conviter, int* maxiter);
extern "C" void mtilde_vm_(scd_t* out, scd_t* in);

extern "C" void bicgstab_dd_d_(scd_t* x, scd_t* b, int* conviter, int* maxiter);
extern "C" void bicgstab_dd_mix_(scd_t* x, scd_t* b, double *tol, int* conviter, int* maxiter, double* tol_s, int* maxiter_s, int* nsap, int* nm);
extern "C" void bicgstab_dd_mix2_hf_(scd_t *x, const scd_t *b, const double *tol, int *conviter, const int *maxiter, const double *tol_s, const int *maxiter_s, const int *nsap, const int *nm);

extern "C" void ddd_d_(scd_t* x, scd_t* b);

extern "C" double get_clock (void);
extern "C" void print_timing_(void);

extern "C" void print_scdnorm2(const char* a, scd_t* in, int size);
extern "C" void print_scsnorm2(const char* a, scs_t* in, int size);
extern "C" void  zero_scd_field(scd_t *out,int i);
extern "C" void  zero_scs_field(scs_t *out,int i);

extern __attribute__((aligned(64))) pglud_t glud;
extern __attribute__((aligned(64))) pclvd_t clvd;
extern __attribute__((aligned(64))) pglus_t glus;
extern __attribute__((aligned(64))) pclvs_t clvs;

extern FILE *para_outputfile;

void generate_random_fields(void);
void generate_random_sc(scd_t *outd, scs_t *outs);
void test_double_prec_functions(scd_t* in);
void test_single_prec_functions(scs_t* in);
//----------------------------------------------------------------------------
int main( int argc, char *argv[] ){
  mt = std::mt19937(12345);

#ifdef _MPI_
  MPI_Init(&argc, &argv);
#endif
  int i, j;

  int lx=atoi(argv[1]);
  int ly=atoi(argv[2]);
  int lz=atoi(argv[3]);
  int lt=atoi(argv[4]);
  int npe_f[4];
  npe_f[0]=atoi(argv[5]);
  npe_f[1]=atoi(argv[6]);
  npe_f[2]=atoi(argv[7]);
  npe_f[3]=atoi(argv[8]);

  //double tol=1e-14;//-1;
  //double tol_s=1e-6;//-1;

  double tol   = atof(argv[9]);
  double tol_s = atof(argv[10]);
  int dd_maxiter= atoi(argv[11]);
  int dd_maxiter_s= atoi(argv[12]);

  int fbc_f[4];
  fbc_f[0]=1;
  fbc_f[1]=1;
  fbc_f[2]=1;
  fbc_f[3]=-1;
  int pce_f=0;
  int pco_f=1;
  
  int block_size[4];
  block_size[0]=1;
  block_size[1]=2;
  block_size[2]=2;
  block_size[3]=2;
  
  qws_init_(&lx, &ly, &lz, &lt, npe_f, fbc_f, &pce_f, &pco_f, block_size);
  //#ifdef _CHECK_PA
  PROF_INIT;
  //#endif

  kappa = (double)0.05;
  kappa2 = kappa * kappa;
  mkappa = - kappa;

  // ------------------------- test case
  generate_random_fields();
#ifdef HALF_PREC
  qws_h_init_(glus,clvs);
#endif
#if 1
  __attribute__((aligned(64))) scd_t* b  = (scd_t*)malloc( sizeof(scd_t) * vold*2);
  __attribute__((aligned(64))) scd_t* x  = (scd_t*)malloc( sizeof(scd_t) * vold*2);
  __attribute__((aligned(64))) scd_t* r  = (scd_t*)malloc( sizeof(scd_t) * vold*2);
  __attribute__((aligned(64))) scs_t* bs = (scs_t*)malloc( sizeof(scs_t) * vols*2);
#else
  scd_t* b;
  scd_t* x;
  scd_t* r;
  scs_t* bs;
  posix_memalign ((void**)&b,  CLS, sizeof(scd_t) * vold*2);
  posix_memalign ((void**)&x,  CLS, sizeof(scd_t) * vold*2);
  posix_memalign ((void**)&r,  CLS, sizeof(scd_t) * vold*2);
  posix_memalign ((void**)&bs, CLS, sizeof(scs_t) * vols*2);
#endif

  // double prec
  generate_random_sc(b,  bs);  // generate double precision random input b and single bs
  test_double_prec_functions(b );
  test_single_prec_functions(bs);
  // solver
  if(true){
  if (rank==0)printf("test bicgstab for double precision mtilde\n");
  int maxiter, iter;
  maxiter =100;
  //  zero_scd_field(x, vold);
#ifdef _MPI_
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  double start_time=get_clock();
  bicgstab_vm_(x, b, &iter, &maxiter);
#ifdef _MPI_
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  double end_time=get_clock();
  mtilde_vm_(r, x);
  rvecd_t rvd0, rvd1, rvd2, tmp;
  double redu[3];
  double rtmp0 = 0;
  double rtmp1 = 0;
  double rtmp2 = 0;
#pragma omp parallel for private(j, rvd0, rvd1, rvd2, tmp) reduction(+:rtmp0, rtmp1, rtmp2)
  for(i=0; i<vold; i++){
    for(j=0; j<24; j++){
      tmp = fsub_d(b[i].ccs[j], r[i].ccs[j]);
      rvd0 = fmul_d(b[i].ccs[j],b[i].ccs[j]);
      rvd1 = fmul_d(tmp,tmp);
      rvd2 = fmul_d(x[i].ccs[j],x[i].ccs[j]);
      rtmp0 += fsum_d(rvd0);
      rtmp1 += fsum_d(rvd1);
      rtmp2 += fsum_d(rvd2);
    }
  }
  redu[0] = rtmp0;
  redu[1] = rtmp1;
  redu[2] = rtmp2;
#ifdef _MPI_
  MPI_Allreduce(MPI_IN_PLACE,(void *)redu,3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
  if(rank==0)printf("bnorm^2     = %24.14e\n",redu[0]);
  if(rank==0)printf("rnorm^2     = %24.14e\n",redu[1]);
  if(rank==0)printf("xnorm^2     = %24.14e\n",redu[2]);
  if(rank==0)printf("rnorm/bnorm = %24.14e\n",sqrt(redu[1]/redu[0]));
  if(rank==0)printf("etime for sovler = %24.14e sec.\n", end_time-start_time);
  }
  rvecd_t rvd0, rvd1, rvd2, tmp;
  double redu[3];
  double rtmp0 = 0;
  double rtmp1 = 0;
  double rtmp2 = 0;
  int iter;
  // ------------------------- test case End
  if(rank==0)printf("\n");
  if(rank==0)printf("\n");
  // ------------------------- LDDHMC solver with Jacobi method instead of SSOR
  if(rank==0)printf("LDDHMC solver with Jacobi method instead of SSOR\n");
  //maxiter  =6;
  //int maxiter_s=50;
  int nsap, nm;
  nsap = 4;
  nm = 2;
  generate_random_sc(x,  bs);

  ////////////////////
  // keep x in xx
  ////////////////////
  scd_t *xx = (scd_t *)malloc( sizeof(scd_t) * vold*2);
#pragma omp parallel for
  for(i=0; i<vold*2; i++)
    for(j=0; j<24; j++)
      for(int v=0; v<VLEND; v++)
        xx[i].ccs[j].v[v] = x[i].ccs[j].v[v];

#ifdef _MPI_
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  double start_time=get_clock();
  bicgstab_dd_mix_(x, b, &tol, &iter, &dd_maxiter, &tol_s, &dd_maxiter_s, &nsap, &nm);
#ifdef _MPI_
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  double end_time=get_clock();

  ddd_d_(r, x);
  rtmp0 = 0;
  rtmp1 = 0;
  rtmp2 = 0;
#pragma omp parallel for private(j, rvd0, rvd1, rvd2, tmp) reduction(+:rtmp0, rtmp1, rtmp2)
  for(i=0; i<vold*2; i++){
    for(j=0; j<24; j++){
      tmp = fsub_d(b[i].ccs[j], r[i].ccs[j]);
      rvd0 = fmul_d(b[i].ccs[j],b[i].ccs[j]);
      rvd1 = fmul_d(tmp,tmp);
      rvd2 = fmul_d(x[i].ccs[j],x[i].ccs[j]);
      rtmp0 += fsum_d(rvd0);
      rtmp1 += fsum_d(rvd1);
      rtmp2 += fsum_d(rvd2);
    }
  }
  redu[0] = rtmp0;
  redu[1] = rtmp1;
  redu[2] = rtmp2;
#ifdef _MPI_
  MPI_Allreduce(MPI_IN_PLACE,(void *)redu,3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
  if(rank==0)printf("bnorm^2     = %24.14e\n",redu[0]);
  if(rank==0)printf("rnorm^2     = %24.14e\n",redu[1]);
  if(rank==0)printf("xnorm^2     = %24.14e\n",redu[2]);
  if(rank==0)printf("rnorm/bnorm = %24.14e\n",sqrt(redu[1]/redu[0]));
  if(rank==0)printf("etime for sovler = %24.14e sec.\n", end_time-start_time);
  // ------------------------- LDDHMC solver End
 

#ifdef HALF_PREC
  ///////////////////////////////////
  // test half precision
  ///////////////////////////////////
  if(rank==0)printf("\n");
  if(rank==0)printf("LDDHMC solver with Jacobi method instead of SSOR using Half Precision\n");

  ///////////////////////
  // restore x from xx
  ///////////////////////
#pragma omp parallel for
  for(i=0; i<vold*2; i++)
    for(j=0; j<24; j++)
      for(int v=0; v<VLEND; v++)
        x[i].ccs[j].v[v] = xx[i].ccs[j].v[v];
  free(xx);

#ifdef _MPI_
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  start_time=get_clock();
  bicgstab_dd_mix2_hf_(x, b, &tol, &iter, &dd_maxiter, &tol_s, &dd_maxiter_s, &nsap, &nm);
#ifdef _MPI_
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  double end_time=get_clock();


  ddd_d_(r, x);
  rtmp0 = 0;
  rtmp1 = 0;
  rtmp2 = 0;
#pragma omp parallel for private(j, rvd0, rvd1, rvd2, tmp) reduction(+:rtmp0, rtmp1, rtmp2)
  for(i=0; i<vold*2; i++){
    for(j=0; j<24; j++){
      tmp  = fsub_d(b[i].ccs[j],r[i].ccs[j]);
      rvd0 = fmul_d(b[i].ccs[j],b[i].ccs[j]);
      rvd1 = fmul_d(tmp,tmp);
      rvd2 = fmul_d(x[i].ccs[j],x[i].ccs[j]);
      rtmp0 += fsum_d(rvd0);
      rtmp1 += fsum_d(rvd1);
      rtmp2 += fsum_d(rvd2);
    }
  }
  redu[0] = rtmp0;
  redu[1] = rtmp1;
  redu[2] = rtmp2;
#ifdef _MPI_
  MPI_Allreduce(MPI_IN_PLACE,(void *)redu,3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
  if(rank==0)printf("bnorm^2     = %24.14e\n",redu[0]);
  if(rank==0)printf("rnorm^2     = %24.14e\n",redu[1]);
  if(rank==0)printf("xnorm^2     = %24.14e\n",redu[2]);
  if(rank==0)printf("rnorm/bnorm = %24.14e\n",sqrt(redu[1]/redu[0]));
  if(rank==0)printf("etime for sovler = %24.14e sec.\n", end_time-start_time);
  ///////////////////////////////////
  // test half precision end
  ///////////////////////////////////
#endif

#ifdef _CHECK_TIMING2
  fprintf(para_outputfile, "\n");
  fprintf(para_outputfile, "\n");
  // ------------------------- print timing
  fprintf(para_outputfile, "print timing\n");
  print_timing_();
#elif defined(_CHECK_TIMING)
  print_timing_();
#endif

  //#ifdef _CHECK_PA
  PROF_FINALIZE;
  //#endif
  if(rank==0)printf("end\n");
  if(rank==0)fflush(0);
  qws_finalize_();
#ifdef _MPI_
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
#endif
  exit(EXIT_SUCCESS);
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
void normal_rand(double *r0, double *r1){
  //double r = sqrt( - log(  (double)rand() / (double)RAND_MAX   ));
  //double t = (double)2.0 * (double)3.1415926535897932384626433832795028841971 * ((double)rand() / (double)RAND_MAX);
  double r = sqrt( - log(  rand01(mt)   ));
  double t = (double)2.0 * (double)3.1415926535897932384626433832795028841971 * rand01(mt);
  *r0 = r * cos(t);
  *r1 = r * sin(t);
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
void generate_random_fields(){
  int j, k, c1, c2, x, y, z, t, mu;
  double r0, r1;

  for(t=0; t<nt; t++){
    for(z=0; z<nz; z++){
      for(y=0; y<ny; y++){
	for(x=0; x<nx; x++){
          int id=(x%nxh)/VLEND + nxd*y + nxd*ny*z + nxd*ny*nz*t + NDIM*vold*(x/nxh);
          int jd=(x%nxh)%VLEND;
          int is=(x%nxh)/VLENS + nxs*y + nxs*ny*z + nxs*ny*nz*t + NDIM*vols*(x/nxh);
          int js=(x%nxh)%VLENS;
//	  int eo=(x+y+z+t)%2;
//          int id=(x/2)/VLEND + nxd*y + nxd*ny*z + nxd*ny*nz*t + NDIM*vold*eo;
//          int jd=(x/2)%VLEND;
//          int is=(x/2)/VLENS + nxs*y + nxs*ny*z + nxs*ny*nz*t + NDIM*vols*eo;
//          int js=(x/2)%VLENS;
          for (mu=0;mu<NDIM;mu++) {
	    for (c1=0;c1<3;c1++) {
	      for (c2=0;c2<2;c2++) {
		normal_rand(&r0, &r1);
		glud[id+ vold*mu].c[c2][c1][0][jd] = r0;
		glud[id+ vold*mu].c[c2][c1][1][jd] = r1;
		glus[is+ vols*mu].c[c2][c1][0][js] = (float)r0;
		glus[is+ vols*mu].c[c2][c1][1][js] = (float)r1;
	      }
	    }
#if SU3_RECONSTRUCT_D == 18 || SU3_RECONSTRUCT_S == 18
	    double u00r = glud[id+ vold*mu].c[0][0][0][jd];
	    double u00i = glud[id+ vold*mu].c[0][0][1][jd];
	    double u01r = glud[id+ vold*mu].c[0][1][0][jd];
	    double u01i = glud[id+ vold*mu].c[0][1][1][jd];
	    double u02r = glud[id+ vold*mu].c[0][2][0][jd];
	    double u02i = glud[id+ vold*mu].c[0][2][1][jd];
	    double u10r = glud[id+ vold*mu].c[1][0][0][jd];
	    double u10i = glud[id+ vold*mu].c[1][0][1][jd];
	    double u11r = glud[id+ vold*mu].c[1][1][0][jd];
	    double u11i = glud[id+ vold*mu].c[1][1][1][jd];
	    double u12r = glud[id+ vold*mu].c[1][2][0][jd];
	    double u12i = glud[id+ vold*mu].c[1][2][1][jd];
#endif
#if SU3_RECONSTRUCT_D == 18
	    glud[id+ vold*mu].c[2][0][0][jd] =  u01r*u12r - u01i*u12i - u02r*u11r + u02i*u11i;
	    glud[id+ vold*mu].c[2][1][0][jd] =  u02r*u10r - u02i*u10i - u00r*u12r + u00i*u12i;
	    glud[id+ vold*mu].c[2][2][0][jd] =  u00r*u11r - u00i*u11i - u01r*u10r + u01i*u10i;
	    glud[id+ vold*mu].c[2][0][1][jd] =- u01r*u12i - u01i+u12r + u02r*u11i + u02i*u11r;
	    glud[id+ vold*mu].c[2][1][1][jd] =- u02r*u10i - u02i+u10r + u00r*u12i + u00i*u12r;
	    glud[id+ vold*mu].c[2][2][1][jd] =- u00r*u11i - u00i+u11r + u01r*u10i + u01i*u10r;
#endif
#if SU3_RECONSTRUCT_S == 18
	    glus[is+ vols*mu].c[2][0][0][js] =(float)(  u01r*u12r - u01i*u12i - u02r*u11r + u02i*u11i );
	    glus[is+ vols*mu].c[2][1][0][js] =(float)(  u02r*u10r - u02i*u10i - u00r*u12r + u00i*u12i );
	    glus[is+ vols*mu].c[2][2][0][js] =(float)(  u00r*u11r - u00i*u11i - u01r*u10r + u01i*u10i );
	    glus[is+ vols*mu].c[2][0][1][js] =(float)( -u01r*u12i - u01i+u12r + u02r*u11i + u02i*u11r );
	    glus[is+ vols*mu].c[2][1][1][js] =(float)( -u02r*u10i - u02i+u10r + u00r*u12i + u00i*u12r );
	    glus[is+ vols*mu].c[2][2][1][js] =(float)( -u00r*u11i - u00i+u11r + u01r*u10i + u01i*u10r );
#endif
	  }
	}
      }
    }
  }

  for(t=0; t<nt; t++){
    for(z=0; z<nz; z++){
      for(y=0; y<ny; y++){
	for(x=0; x<nx; x++){
          int id = (x%nxh)/VLEND + nxd*y + nxd*ny*z + nxd*ny*nz*t + vold*(x/nxh);
          int jd = (x%nxh)%VLEND;
	  int is = (x%nxh)/VLENS + nxs*y + nxs*ny*z + nxs*ny*nz*t + vols*(x/nxh);
	  int js = (x%nxh)%VLENS;

//	  int eo=(x+y+z+t)%2;
//          int id=(x/2)/VLEND + nxd*y + nxd*ny*z + nxd*ny*nz*t + vold*eo;
//          int jd=(x/2)%VLEND;
//          int is=(x/2)/VLENS + nxs*y + nxs*ny*z + nxs*ny*nz*t + vols*eo;
//          int js=(x/2)%VLENS;
	  for (j=0;j<2;j++) {
	    for (k=0;k<18;k++) {
	      normal_rand(&r0, &r1);
	      clvd[id].c[j][2*k+0][jd] = r0/(double)10;
	      clvd[id].c[j][2*k+1][jd] = r1/(double)10;
	    }
	    clvd[id].c[j][0][jd] += 1;
	    clvd[id].c[j][1][jd] += 1;
	    clvd[id].c[j][2][jd] += 1;
	    clvd[id].c[j][3][jd] += 1;
	    clvd[id].c[j][4][jd] += 1;
	    clvd[id].c[j][5][jd] += 1;
	  }
	  for(j=0;j<2;j++){
	    for (k=0;k<36;k++){
	      clvs[is].c[j][k][js] = (float)clvd[id].c[j][k][jd];
	    }
	  }
	}
      }
    }
  }

}
//----------------------------------------------------------------------------
void generate_random_sc(scd_t *outd, scs_t *outs){
  double r0, r1;

  for(int t=0; t<nt; t++){
    for(int z=0; z<nz; z++){
      for(int y=0; y<ny; y++){
        for(int x=0; x<nx; x++){
	  int id = (x%nxh)/VLEND + nxd*y + nxd*ny*z + nxd*ny*nz*t + vold*(x/nxh);
	  int jd = (x%nxh)%VLEND;
	  int is = (x%nxh)/VLENS + nxs*y + nxs*ny*z + nxs*ny*nz*t + vols*(x/nxh);
	  int js = (x%nxh)%VLENS;

//	  int eo=(x+y+z+t)%2;
//          int id=(x/2)/VLEND + nxd*y + nxd*ny*z + nxd*ny*nz*t + NDIM*vold*eo;
//          int jd=(x/2)%VLEND;
//          int is=(x/2)/VLENS + nxs*y + nxs*ny*z + nxs*ny*nz*t + NDIM*vols*eo;
//          int js=(x/2)%VLENS;

	  for (int c=0;c<3;c++) {
	    for (int s=0;s<4;s++) {
	      normal_rand(&r0, &r1);
	      outd[id].c[c][s][0][jd]=r0;
	      outd[id].c[c][s][1][jd]=r1;
	      outs[is].c[c][s][0][js]=(float)r0;
	      outs[is].c[c][s][1][js]=(float)r1;
	    }
	  }
	}
      }
    }
  }
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
extern "C" void deo_in_vm_(int* pe, int* po, scd_t* out, scd_t* in);
extern "C" void dee_deo_in_vm_(int* pe, int* po, scd_t* out, scd_t* in);
extern "C" void deo_dag_in_vm_(int* pe, int* po, scd_t* out, scd_t* in);
extern "C" void dee_deo_dag_in_vm_(int* pe, int* po, scd_t* out, scd_t* in);
extern "C" void one_minus_dee_deo_in_vm_(int* pe, int* po, scd_t* out, scd_t* in, scd_t* in0, double factor);
extern "C" void one_minus_deo_dag_in_vm_(int* pe, int* po, scd_t* out, scd_t* in, scd_t* in0, double factor);
extern "C" void deo_out_pre_vm_(int* pe, int* po, scd_t* out, scd_t* in);
extern "C" void deo_dag_out_pre_vm_(int* pe, int* po, scd_t* out, scd_t* in);
extern "C" void deo_out_pos_vm_(int* pe, int* po, scd_t* out, scd_t* in, double factor);
extern "C" void deo_dag_out_pos_vm_(int* pe, int* po, scd_t* out, scd_t* in, double factor);
extern "C" void dee_deo_out_pos_vm_(int* pe, int* po, scd_t* out, scd_t* in, double factor);
extern "C" void dee_deo_dag_out_pos_vm_(int* pe, int* po, scd_t* out, scd_t* in, double factor);

extern "C" void clvd_vm_(int* pe, scd_t* inout);
extern "C" void clvd_vm2_(int* pe, scd_t* out, scd_t* in);
extern "C" void deo_vm_(int* pe, int* po, scd_t* out, scd_t* in);
extern "C" void dee_deo_vm_(int* pe, int* po, scd_t* out, scd_t* in);
extern "C" void one_minus_dee_deo_vm_(int* pe, int* po, scd_t* out, scd_t* in, scd_t* in0);
extern "C" void deo_dag_vm_(int* pe, int* po, scd_t* out, scd_t* in);
extern "C" void dee_deo_dag_vm_(int* pe, int* po, scd_t* out, scd_t* in);
extern "C" void one_minus_deo_dag_vm_(int* pe, int* po, scd_t* out, scd_t* in, scd_t* in0);
extern "C" void mtilde_vm_(scd_t* out, scd_t* in);
extern "C" void mtilde_dag_vm_(scd_t* out, scd_t* in);
extern "C" void mtdagmt_vm_(scd_t* out, scd_t* in);

extern "C" void ddd_out_pre_d_(scd_t* in, int* DEO);
extern "C" void ddd_out_pos_d_(scd_t* out, scd_t* in, int* DEO);
extern "C" void ddd_in_d_(scd_t* out, scd_t* in, int* DEO);
extern "C" void ddd_d_(scd_t* out, scd_t* in);
//----------------------------------------------------------------------------
void test_double_prec_functions(scd_t* in){
  __attribute__((aligned(64))) scd_t* out = (scd_t*)malloc( sizeof(scd_t) * vold*2);
  int pe=0;
  int po=1;
  double fac=-kappa2;
  if(rank==0)printf("test_double_prec_functions\n");
  print_scdnorm2("input",in, vold*2);

  zero_scd_field(out,vold);
  deo_in_vm_(&po, &pe,out,in);
  print_scdnorm2("deo_in_vm_",out,vold);

  zero_scd_field(out,vold);
  dee_deo_in_vm_(&po, &pe, out, in);
  print_scdnorm2("dee_deo_in_vm_",out,vold);

  zero_scd_field(out,vold);
  deo_dag_in_vm_(&po, &pe, out, in);
  print_scdnorm2("deo_dag_in_vm_",out,vold);

  zero_scd_field(out,vold);
  dee_deo_dag_in_vm_(&po, &pe, out, in);
  print_scdnorm2("dee_deo_dag_in_vm_",out,vold);
  
  zero_scd_field(out,vold);
  one_minus_dee_deo_in_vm_(&po, &pe, out, in, in, fac);
  print_scdnorm2("one_minus_dee_deo_in_vm_",out,vold);

  zero_scd_field(out,vold);
  one_minus_deo_dag_in_vm_(&po, &pe, out, in, in, fac);
  print_scdnorm2("one_minus_deo_dag_in_vm_",out,vold);

  zero_scd_field(out,vold);
  deo_out_pre_vm_(&po, &pe, out, in);
  deo_out_pos_vm_(&po, &pe, out, in, fac);
  print_scdnorm2("deo_out_pre|pos_vm_",out,vold);

  zero_scd_field(out,vold);
  deo_dag_out_pre_vm_(&po, &pe, out, in);
  deo_dag_out_pos_vm_(&po, &pe, out, in, fac);
  print_scdnorm2("deo_dag_out_pre|pos_vm_",out,vold);

  zero_scd_field(out,vold);
      deo_out_pre_vm_(&po, &pe, out, in);
  dee_deo_out_pos_vm_(&po, &pe, out, in, fac);
  print_scdnorm2("dee_deo_out_pre|pos_vm_",out,vold);

  zero_scd_field(out,vold);
      deo_dag_out_pre_vm_(&po, &pe, out, in);
  dee_deo_dag_out_pos_vm_(&po, &pe, out, in, fac);
  print_scdnorm2("dee_deo_dag_out_pre|pos_vm_",out,vold);


  memcpy(out, in, sizeof(scd_t)*vold);
  clvd_vm_(&pe, out);
  print_scdnorm2("clvd_vm_",out,vold);

  zero_scd_field(out,vold);
  clvd_vm2_(&pe, out, in);
  print_scdnorm2("clvd_vm2_",out,vold);

  zero_scd_field(out,vold);
  deo_vm_(&po, &pe, out, in);
  print_scdnorm2("deo_vm_",out,vold);
  
  zero_scd_field(out,vold);
  dee_deo_vm_(&po, &pe, out, in);
  print_scdnorm2("dee_deo_vm_",out,vold);
  
  zero_scd_field(out,vold);
  one_minus_dee_deo_vm_(&po, &pe, out, in, in);
  print_scdnorm2("one_minus_dee_deo_vm_",out,vold);
  
  zero_scd_field(out,vold);
  deo_dag_vm_(&po, &pe, out, in);
  print_scdnorm2("deo_dag_vm_",out,vold);

  zero_scd_field(out,vold);
  dee_deo_dag_vm_(&po, &pe, out, in);
  print_scdnorm2("dee_deo_dag_vm_",out,vold);

  zero_scd_field(out,vold);
  one_minus_deo_dag_vm_(&po, &pe, out, in, in);
  print_scdnorm2("one_minus_deo_dag_vm_",out,vold);

  zero_scd_field(out,vold*2);
  mtilde_vm_(out, in);
  print_scdnorm2("mtilde_vm_",out,vold);

  zero_scd_field(out,vold);
  mtilde_dag_vm_(out, in);
  print_scdnorm2("mtilde_dag_vm_",out,vold);

  zero_scd_field(out,vold);
  mtdagmt_vm_(out, in);
  print_scdnorm2("mtdagmt_vm_",out,vold);

  int DEO=0;

  zero_scd_field(out,vold*2);
  ddd_out_pre_d_(in, &DEO);
  ddd_out_pos_d_(out, in, &DEO);
  print_scdnorm2("ddd_out_pre|pos_d_",out,vold*2);

  zero_scd_field(out,vold*2);
  ddd_in_d_(out, in, &DEO);
  print_scdnorm2("ddd_in_d_",out,vold*2);

  zero_scd_field(out,vold*2);
  ddd_d_(out, in);
  print_scdnorm2("ddd_d_",out,vold*2);


  free(out);
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

extern "C" void deo_in_s_(int* pe, int* po, scs_t* out, scs_t* in);
extern "C" void dee_deo_in_s_(int* pe, int* po, scs_t* out, scs_t* in);
extern "C" void one_minus_dee_deo_in_s_(int* pe, int* po, scs_t* out, scs_t* in, scs_t* in0, float factor);
extern "C" void deo_out_pre_s_(int* pe, int* po, scs_t* out, scs_t* in);
extern "C" void deo_out_pos_s_(int* pe, int* po, scs_t* out, scs_t* in, float factor);
extern "C" void dee_deo_out_pos_s_(int* pe, int* po, scs_t* out, scs_t* in, float factor);
extern "C" void clv_s_(int* pe, scs_t* inout);
extern "C" void deo_s_(int* pe, int* po, scs_t* out, scs_t* in);
extern "C" void dee_deo_s_(int* pe, int* po, scs_t* out, scs_t* in);
extern "C" void one_minus_dee_deo_s_(int* pe, int* po, scs_t* out, scs_t* in, scs_t* in0);
extern "C" void mtilde_s_(scs_t* out, scs_t* in);
extern "C" void ddd_out_pre_s_(scs_t* in, int* DEO);
extern "C" void ddd_out_pre_s_no_timer_(scs_t* in, int* DEO);
extern "C" void ddd_out_pos_s_(scs_t* out, scs_t* in, int* DEO, float factor);
extern "C" void ddd_out_pos_s_no_timer_(scs_t* out, scs_t* in, int* DEO, float factor);
extern "C" void ddd_in_s_(scs_t* out, scs_t* in, int* DEO);
extern "C" void ddd_eo_s_(scs_t* out, scs_t* in, int* DEO);
extern "C" void ddd_s_(scs_t* out, scs_t* in);

//----------------------------------------------------------------------------
void test_single_prec_functions(scs_t* in){
  __attribute__((aligned(64))) scs_t* out = (scs_t*)malloc( sizeof(scs_t) * vols*2);
  int pe=0;
  int po=1;
  float fac = (float) -kappa2;

  if(rank==0)printf("test_single_prec_functions\n");
  print_scsnorm2("input",in, vols*2);

  deo_in_s_(&po, &pe, out, in);
  print_scsnorm2("deo_in_s_",out,vols);

  zero_scs_field(out,vols);
  dee_deo_in_s_(&po, &pe, out, in);
  print_scsnorm2("dee_deo_in_s_",out,vols);

  zero_scs_field(out,vols);
  one_minus_dee_deo_in_s_(&po, &pe, out, in, in, fac);
  print_scsnorm2("one_minus_dee_deo_in_s_",out,vols);

  zero_scs_field(out,vols);
  deo_out_pre_s_(&po, &pe, out, in);
  deo_out_pos_s_(&po, &pe, out, in, fac);
  print_scsnorm2("deo_out_pre|pos_s_",out,vols);

  zero_scs_field(out,vols);
      deo_out_pre_s_(&po, &pe, out, in);
  dee_deo_out_pos_s_(&po, &pe, out, in, fac);
  print_scsnorm2("dee_deo_out_pre|pos_s_",out,vols);

  memcpy(out, in, sizeof(scs_t)*vols);
  clv_s_(&pe, out);
  print_scsnorm2("clv_s_",out,vols);

  zero_scs_field(out,vols);
  deo_s_(&po, &pe, out, in);
  print_scsnorm2("deo_s_",out,vols);

  zero_scs_field(out,vols);
  dee_deo_s_(&po, &pe, out, in);
  print_scsnorm2("dee_deo_s_",out,vols);

  zero_scs_field(out,vols);
  one_minus_dee_deo_s_(&po, &pe, out, in, in);
  print_scsnorm2("one_minus_dee_deo_s_",out,vols);

  zero_scs_field(out,vols);
  mtilde_s_(out, in);
  print_scsnorm2("mtilde_s_",out,vols);

  int DEO=0;

  zero_scs_field(out,vols*2);
  ddd_out_pre_s_no_timer_(in, &DEO);
  ddd_out_pos_s_no_timer_(out, in, &DEO, (float)mkappa);
  print_scsnorm2("ddd_out_pre|pos_s_",out,vols*2);

  zero_scs_field(out,vols*2);
  ddd_in_s_(out, in, &DEO);
  print_scsnorm2("ddd_in_s_",out,vols*2);

  //zero_scs_field(out,vols*2);
  //ddd_eo_s_(out, in, &DEO);
  //print_scsnorm2("ddd_eo_s_",out,vols*2);

  zero_scs_field(out,vols*2);
  ddd_s_(out, in);
  print_scsnorm2("ddd_s_",out,vols*2);

  free(out);
}
