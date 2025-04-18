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
#define _PFI1 2
#define _PFI2 1
#define _PFI3 1
#define _PFI4 2
#include "wilson_d.h"
#include "wilson_s.h"
#include "clover_d.h"
#include "clover_s.h"
#include "clover_def.h"
#include "prefetch.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#define min(a,b) (a)>(b)?(b):(a) 

#include "eml_lib.h"

#ifdef __cplusplus
extern "C"{
#endif
#include "timing.h"

  __attribute__((aligned(64))) pglud_t glud;
  __attribute__((aligned(64))) pglus_t glus;
  __attribute__((aligned(64))) pclvd_t clvd;
  __attribute__((aligned(64))) pclvs_t clvs;
  //projscd_t *xfd_buff;
  //projscd_t *xbd_buff;
  projscd1_t *xfd_send, *xfd_recv;
  projscd1_t *xbd_send, *xbd_recv;
  projscd_t *yfd_send, *yfd_recv;
  projscd_t *ybd_send, *ybd_recv;
  projscd_t *zfd_send, *zfd_recv;
  projscd_t *zbd_send, *zbd_recv;
  projscd_t *tfd_send, *tfd_recv;
  projscd_t *tbd_send, *tbd_recv;

  projscs1_t *xfs_send, *xfs_recv;
  projscs1_t *xbs_send, *xbs_recv;
  projscs_t *yfs_send, *yfs_recv;
  projscs_t *ybs_send, *ybs_recv;
  projscs_t *zfs_send, *zfs_recv;
  projscs_t *zbs_send, *zbs_recv;
  projscs_t *tfs_send, *tfs_recv;
  projscs_t *tbs_send, *tbs_recv;


  scs_t *qs;

  double kappa, kappa2, mkappa;

  //  nh=nx/2, nd=nx/2/8, nxs=nx/2/16
  int nx, ny, nz, nt, nxh, nxd, nxs;
  int vold;// = nxd * ny * nz * nt;
  int vols;// = nxs * ny * nz * nt;
  int npe[4];// number of proccesses
  double fbc[4][2];// fermion boundary condiion [mu][0=fwd, 1=bwd]
  int rank, size, px, py, pz, pt;// process coordinate
  int pxf, pyf, pzf, ptf;
  int pxb, pyb, pzb, ptb;
  int domain_e, domain_o;
  int thmax,volsize;
  int iam;
#pragma omp threadprivate(iam)
  int volse;
#pragma omp threadprivate(volse)

  int *pce, *pco;//even-odd precondition notation

  block_map_t* block_map;
  int num_blocks;

  FILE *para_outputfile=0;
#ifdef _POWER_API_
  extern void power_api_init();
  extern void power_api_finalize();
#endif
  //---------------------------------------------------------------------------------------- qws init
  void qws_init_(int* lx,  int* ly, int* lz, int* lt, 
		 int* npe_f, int* fbc_f, int* pce_f, int* pco_f, int* block_size){

    int i;
    pce = pce_f;
    pco = pco_f;

    for (i=0;i<4;i++){
      npe[i]=npe_f[i];
    }

    nx=*lx;  ny=*ly;  nz=*lz;  nt=*lt;
    nxh=nx/2;
    nxd=nx/2/VLEND;
    nxs=nx/2/VLENS;
    if (nxd == 0) {printf("ERROR at %s : nx/2/VLEND must be >=1\n", __func__); exit(1);}
    if (nxs == 0) {printf("ERROR at %s : nx/2/VLENS must be >=1\n", __func__); exit(1);}
    vold = nxd * ny * nz * nt;
    vols = nxs * ny * nz * nt;
#ifdef _OPENMP
#pragma omp parallel 
    {
      iam = omp_get_thread_num();
#pragma omp master
      {
        thmax = omp_get_num_threads();
        volsize = vols / thmax;

      }
#pragma omp barrier
      volse = volsize * (iam+1);
    }
#else
    iam=0;
    thmax=1;
    volsize=vols;
    volse=vols;

#endif

    glud = (pglud_t)malloc( sizeof(glud_t) * vold*NDIM*NEO); /////////////////
    clvd = (pclvd_t)malloc( sizeof(clvd_t) * vold*NEO);

    void *tmp;

    posix_memalign(&tmp, CLS, sizeof(glus_t) * vols*NDIM*NEO);
    glus = (pglus_t)tmp;
    posix_memalign(&tmp, CLS, sizeof(clvs_t) * vols*NEO);
    clvs = (pclvs_t)tmp;
    posix_memalign(&tmp, CLS, sizeof(scs_t)*vols);
    qs = (scs_t*)tmp;
    for(i=0; i<vols; i++){
      for(int j=0; j<24; j++){
	for(int v=0; v<VLENS; v++){
	  qs[i].ccs[j].v[v]  = 0.0f;
	}
      }
    }

#ifdef _MPI_
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#else
    rank = 0;
    size = 1;
#endif

    if(rank==0)printf("QWS git commit hash : %s\n",VERSION);

    if (size != npe[0]*npe[1]*npe[2]*npe[3]) {
      printf("MPI size /= npe[0]*npe[1]*npe[2]*npe[3]\n size=%d, npe[0]*npe[1]*npe[2]*npe[3]=%d\nStop.\n",size,npe[0]*npe[1]*npe[2]*npe[3]);
#ifdef _MPI_
      MPI_Finalize();
#endif
    }
#ifdef _CHECK_TIMING
#ifdef _TIMING_EACH_RANK
    char filename[21];
    sprintf(filename, "output_rank%06d.txt", rank);
    para_outputfile = fopen(filename, "w");
#endif
#endif

    int map_id=-1;
    int rank_coords[4];
    int rank_size[4];
    int rankf[4];
    int rankb[4];
    map_id=xbound_get_rankmap(rank, rank_coords, rank_size, rankf, rankb);

    if (size == 1){
      printf("single %d\n", rank);
      fflush(stdout);
      px=0;    py=0;    pz=0;    pt=0;
      domain_e = ( py + pz + pt)%2;
      domain_o = 1 - domain_e;
    } else {
      if(map_id>=0){ //
	  if(rank==0){
	    printf("rankmap is found: map_id=%d\n", map_id);
            fflush(stdout);
	  }

	  if(npe[0] != rank_size[0] ||npe[1] != rank_size[1] ||
	     npe[2] != rank_size[2] ||npe[3] != rank_size[3] ){
	    if(rank==0){
	      printf("  WARNING!!! new rank size: %d %d %d %d\n", rank_size[0], rank_size[1], rank_size[2], rank_size[3]);
              fflush(stdout);
	    }
	  }
	for(int i=0; i<4; i++){
	  npe[i]=rank_size[i];
	}
	px=rank_coords[0];
	py=rank_coords[1];
	pz=rank_coords[2];
	pt=rank_coords[3];
	pxf=rankf[0];
	pyf=rankf[1];
	pzf=rankf[2];
	ptf=rankf[3];
	pxb=rankb[0];
	pyb=rankb[1];
	pzb=rankb[2];
	ptb=rankb[3];

      } else { // no rankmpa is found

	px = rank;
	pt = px / (npe[0]*npe[1]*npe[2]);
	px-= npe[0]*npe[1]*npe[2]*pt;
	pz = px / (npe[0]*npe[1]);
	px-= npe[0]*npe[1]*pz;
	py = px / npe[0];
	px-= npe[0]*py;
	pxf = ((       px+1)%npe[0]) + npe[0]*py + npe[0]*npe[1]*pz + npe[0]*npe[1]*npe[2]*pt;
	pxb = ((npe[0]+px-1)%npe[0]) + npe[0]*py + npe[0]*npe[1]*pz + npe[0]*npe[1]*npe[2]*pt;
	pyf = px + ((       py+1)%npe[1])*npe[0] + npe[0]*npe[1]*pz + npe[0]*npe[1]*npe[2]*pt;
	pyb = px + ((npe[1]+py-1)%npe[1])*npe[0] + npe[0]*npe[1]*pz + npe[0]*npe[1]*npe[2]*pt;
	pzf = px + npe[0]*py + ((       pz+1)%npe[2])*npe[0]*npe[1] + npe[0]*npe[1]*npe[2]*pt;
	pzb = px + npe[0]*py + ((npe[2]+pz-1)%npe[2])*npe[0]*npe[1] + npe[0]*npe[1]*npe[2]*pt;
	ptf = px + npe[0]*py + npe[0]*npe[1]*pz + ((       pt+1)%npe[3])*npe[0]*npe[1]*npe[2];
	ptb = px + npe[0]*py + npe[0]*npe[1]*pz + ((npe[3]+pt-1)%npe[3])*npe[0]*npe[1]*npe[2];
      }

      domain_e = ( py + pz + pt)%2;
      domain_o = 1 - domain_e;

#ifdef _DEBUG_
      printf("MPI %d\n", px);
      printf("MPI rank=%d px=%d py=%d pz=%d pt=%d\n", rank, px, py, pz, pt);
      printf("MPI rank=%d pxf=%d pxb=%d (forward and backward ranks in X-dir)\n", rank, pxf, pxb);
      printf("MPI rank=%d pyf=%d pyb=%d (forward and backward ranks in Y-dir)\n", rank, pyf, pyb);
      printf("MPI rank=%d pzf=%d pzb=%d (forward and backward ranks in Z-dir)\n", rank, pzf, pzb);
      printf("MPI rank=%d ptf=%d ptb=%d (forward and backward ranks in T-dir)\n", rank, ptf, ptb);
      printf("MPI rank=%d domain_e=%d\n", rank, domain_e);
      printf("MPI rank=%d domain_o=%d\n", rank, domain_o);
      fflush(stdout);
#endif

    }
    // initialize communicaiton
    xbound_init(rank,
                pxf,pyf,pzf,ptf,pxb,pyb,pzb,ptb);
#ifdef _DEBUG_
    printf("qws_init end: MPI rank=%d px=%d py=%d pz=%d pt=%d\n", rank, px, py, pz, pt);
#endif

    if (pt == npe[3]-1 && pt != 0 && fbc_f[3] != 1){
      fbc[3][0]=(double)(fbc_f[3]*2);
      fbc[3][1]=(double)2;
    } else if (pt != npe[3]-1 && pt == 0 && fbc_f[3] != 1){
      fbc[3][0]=(double)2;
      fbc[3][1]=(double)(fbc_f[3]*2);
    } else if (pt == npe[3]-1 && pt == 0 && fbc_f[3] != 1){
      fbc[3][0]=(double)(fbc_f[3]*2);
      fbc[3][1]=(double)(fbc_f[3]*2);
    } else {
      fbc[3][0]=(double)2;
      fbc[3][1]=(double)2;
    }


#ifdef _POWER_API_
    power_api_init();
#endif
  }
  //---------------------------------------------------------------------------------------- qws init end

  void qws_finalize_() {

    xbound_finalize();

    if(glud) free(glud);
    if(glus) free(glus);
    if(clvd) free(clvd);
    if(clvs) free(clvs);

#ifdef _POWER_API_
    power_api_finalize();
#endif
#ifdef _CHECK_TIMING
#ifdef _TIMING_EACH_RANK
    if(!para_outputfile)fclose(para_outputfile);
#endif 
#endif
  }

  //---------------------------------------------------------------------------------------- for testing, not using
  void deo_test_(int* pe, int* po, scd_t* out, scd_t* in){
    __attribute__((aligned(64))) scd_t tmp;
    int x, y, z, t;

    pglud_t gxf = &glud[vold*0 + NDIM*vold*(*pe)];
    pglud_t gxb = &glud[vold*0 + NDIM*vold*(*po)];
    pglud_t gyf = &glud[vold*1 + NDIM*vold*(*pe)];
    pglud_t gyb = &glud[vold*1 + NDIM*vold*(*po)];
    pglud_t gzf = &glud[vold*2 + NDIM*vold*(*pe)];
    pglud_t gzb = &glud[vold*2 + NDIM*vold*(*po)];
    pglud_t gtf = &glud[vold*3 + NDIM*vold*(*pe)];
    pglud_t gtb = &glud[vold*3 + NDIM*vold*(*po)];

    for(t=0; t<nt; t++){
      for(z=0; z<nz; z++){
	for(y=0; y<ny; y++){
	  int ick=  ( (*pe) + y     + z + t)%2;
	  //printf("ick e y z t %d %d %d %d %d\n",ick,e,y,z,t);
	  for(x=0; x<nxd; x++){
	    //printf("iiiii x=%d y=%d z=%d t=%d iyf=%d\n",x, y, z, t, iyf);
	    __zero_sc(tmp.c);
	    int i0  =      x        + nxd*y + nxd*ny*z + nxd*ny*nz*t;
	    int ixf = (    x+1)%nxd + nxd*y + nxd*ny*z + nxd*ny*nz*t;
	    int ixb = (nxd+x-1)%nxd + nxd*y + nxd*ny*z + nxd*ny*nz*t;
	    int iyf = x + ((   y+1)%ny)*nxd + nxd*ny*z + nxd*ny*nz*t;
	    int iyb = x + ((ny+y-1)%ny)*nxd + nxd*ny*z + nxd*ny*nz*t;
	    int izf = x + nxd*y + ((   z+1)%nz)*nxd*ny + nxd*ny*nz*t;
	    int izb = x + nxd*y + ((nz+z-1)%nz)*nxd*ny + nxd*ny*nz*t;
	    int itf = x + nxd*y + nxd*ny*z + ((   t+1)%nt)*nxd*ny*nz;
	    int itb = x + nxd*y + nxd*ny*z + ((nt+t-1)%nt)*nxd*ny*nz;
	    if (ick == 0) {
	      __mult_xfwd0(tmp.c, gxf, in, i0);
	      __mult_xbwd0(tmp.c, gxb, in, i0, ixb);
	    } else {
	      __mult_xfwd1(tmp.c, gxf, in, i0, ixf);
	      __mult_xbwd1(tmp.c, gxb, in, i0);
	    }
	    __mult_yfwd(tmp.c, gyf, in, i0, iyf);
	    __mult_ybwd(tmp.c, gyb, in,     iyb);
	    __mult_zfwd(tmp.c, gzf, in, i0, izf);
	    __mult_zbwd(tmp.c, gzb, in,     izb);
	    __mult_tfwd(tmp.c, gtf, in, i0, itf);
	    __mult_tbwd(tmp.c, gtb, in,     itb);
	    __store_sc(out[i0].c, tmp.c);
	  }
	}
      }
    }
  }

  //----------------------------------------------------------------------------------------
  void clvd_vm_(int* pe, scd_t* inout){
    __attribute__((aligned(64))) scd_t tmp;
    int x, y, z, t;
#pragma omp for private(x,y,z,t,tmp) collapse(3) schedule(static)
    for(t=0; t<nt; t++){
      for(z=0; z<nz; z++){
	for(y=0; y<ny; y++){
	  for(x=0; x<nxd; x++){
	    int i0  = x + nxd*y + nxd*ny*z + nxd*ny*nz*t;
	    __load_sc(tmp.c, inout[i0].c);
	    //__mult_clvd( tmp.cv, clv[i0 + vold*(*pe)].cv);
	    __mult_clvd_def( tmp.cv, clvd[i0 + vold*(*pe)].cv);
	    __store_sc(inout[i0].c, tmp.c);
	  }
	}
      }
    }
  }
  //----------------------------------------------------------------------------------------
  void clv_s_(int* pe, scs_t* inout){
    __attribute__((aligned(64))) scs_t tmp;
    int x, y, z, t;
#pragma omp for private(x,y,z,t,tmp) collapse(3) schedule(static)
    for(t=0; t<nt; t++){
      for(z=0; z<nz; z++){
	for(y=0; y<ny; y++){
	  for(x=0; x<nxs; x++){
	    int i0  = x + nxs*y + nxs*ny*z + nxs*ny*nz*t;
	    __load_sc_s(tmp.c, inout[i0].c);
	    __mult_clvs( tmp.cv, clvs[i0 + vols*(*pe)].cv);
	    __store_sc_s(inout[i0].c, tmp.c);
	  }
	}
      }
    }
  }
  //----------------------------------------------------------------------------------------
  void clvd_vm2_(int* pe, scd_t* out, scd_t* in){
    __attribute__((aligned(64))) scd_t tmp;
    int x, y, z, t;
#pragma omp for private(x,y,z,t,tmp) collapse(3) schedule(static)
    for(t=0; t<nt; t++){
      for(z=0; z<nz; z++){
	for(y=0; y<ny; y++){
	  for(x=0; x<nxd; x++){
	    int i0  = x + nxd*y + nxd*ny*z + nxd*ny*nz*t;
	    __load_sc(tmp.c, in[i0].c);
	    __mult_clvd( tmp.cv, clvd[i0 + vold*(*pe)].cv);
	    __store_sc(out[i0].c, tmp.c);
	  }
	}
      }
    }
  }
  //----------------------------------------------------------------------------
  void zero_scd_field(scd_t* inout, int size){
#pragma omp parallel for
    for (int i=0;i<size;i++) {
      for (int j=0;j<24;j++) {
	for (int iv=0;iv<VLEND;iv++) {
	  inout[i].ccs[j].v[iv]=0;
	}
      }
    }
  }
  //----------------------------------------------------------------------------
  void zero_scs_field(scs_t* inout, int size){
#pragma omp parallel for
    for (int i=0;i<size;i++) {
      for (int j=0;j<24;j++) {
	for (int iv=0;iv<VLENS;iv++) {
	  inout[i].ccs[j].v[iv]=0;
	}
      }
    }
  }
  //----------------------------------------------------------------------------
  void print_scdnorm2(const char* a, scd_t* in, int size){
    double rtmp0 = 0;
#pragma omp parallel for reduction(+:rtmp0)
    for(int i=0; i<size; i++){
      for(int j=0; j<24; j++){
	rvecd_t rvd0 = fmul_d(in[i].ccs[j],in[i].ccs[j]);
	rtmp0 += fsum_d(rvd0);
      }
    }
#ifdef _MPI_
    MPI_Allreduce(MPI_IN_PLACE, &rtmp0,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
    if(rank==0){
      printf("%30s : %24.14e\n", a, rtmp0);
      fflush(stdout);
    }
  }
  //----------------------------------------------------------------------------
  void print_scsnorm2(char* a, scs_t* in, int size){
    float rtmp0 = 0;
#pragma omp parallel for reduction(+:rtmp0)
    for(int i=0; i<size; i++){
      for(int j=0; j<24; j++){
        vecs_t v = fload1_s(pred_true_all(), &(in[i].ccs[j]), dims_1d, 0);
        vecs_t rvd0 = fmul_s(pred_true_all(), v, v);
        rtmp0 += fsum_s(pred_true_all(), rvd0);
      }
    }
#ifdef _MPI_
    MPI_Allreduce(MPI_IN_PLACE, &rtmp0,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD);
#endif
    if(rank==0){
      printf("%30s : %24.14e\n", a, rtmp0);
      fflush(stdout);
    }
  }
  //----------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------
  void deo_in_vm_(int* pe, int* po, scd_t* out, scd_t* in);
  void dee_deo_in_vm_(int* pe, int* po, scd_t* out, scd_t* in);
  void deo_dag_in_vm_(int* pe, int* po, scd_t* out, scd_t* in);
  void dee_deo_dag_in_vm_(int* pe, int* po, scd_t* out, scd_t* in);
  void one_minus_dee_deo_in_vm_(int* pe, int* po, scd_t* out, scd_t* in, scd_t* in0, double factor);
  void one_minus_deo_dag_in_vm_(int* pe, int* po, scd_t* out, scd_t* in, scd_t* in0, double factor);
  void deo_out_pre_vm_(int* pe, int* po, scd_t* out, scd_t* in);
  void deo_dag_out_pre_vm_(int* pe, int* po, scd_t* out, scd_t* in);
  void deo_out_pos_vm_(int* pe, int* po, scd_t* out, scd_t* in, double factor);
  void deo_dag_out_pos_vm_(int* pe, int* po, scd_t* out, scd_t* in, double factor);
  void dee_deo_out_pos_vm_(int* pe, int* po, scd_t* out, scd_t* in, double factor);
  void dee_deo_dag_out_pos_vm_(int* pe, int* po, scd_t* out, scd_t* in, double factor);
  //---------------------------------------------------------------------------------------- double precision
  void deo_vm_(int* pe, int* po, scd_t* out, scd_t* in){
    deo_out_pre_vm_(pe, po, out, in);
    deo_in_vm_(pe, po, out, in);
    deo_out_pos_vm_(pe, po, out, in, 1);
  }
  void dee_deo_vm_(int* pe, int* po, scd_t* out, scd_t* in){
    deo_out_pre_vm_(pe, po, out, in);
    dee_deo_in_vm_(pe, po, out, in);
    dee_deo_out_pos_vm_(pe, po, out, in, 1);
  }
  void one_minus_dee_deo_vm_(int* pe, int* po, scd_t* out, scd_t* in, scd_t* in0){
    double fac = -kappa2;
    deo_out_pre_vm_(pe, po, out, in);
    one_minus_dee_deo_in_vm_(pe, po, out, in, in0, fac);
    dee_deo_out_pos_vm_(pe, po, out, in, fac );
  }
  void deo_dag_vm_(int* pe, int* po, scd_t* out, scd_t* in){
    deo_dag_out_pre_vm_(pe, po, out, in);
    deo_dag_in_vm_(pe, po, out, in);
    deo_dag_out_pos_vm_(pe, po, out, in, 1);
  }
  void dee_deo_dag_vm_(int* pe, int* po, scd_t* out, scd_t* in){
    deo_dag_out_pre_vm_(pe, po, out, in);
    dee_deo_dag_in_vm_(pe, po, out, in);
    dee_deo_dag_out_pos_vm_(pe, po, out, in, 1);
  }
  void one_minus_deo_dag_vm_(int* pe, int* po, scd_t* out, scd_t* in, scd_t* in0){
    double fac = -kappa2;
    deo_dag_out_pre_vm_(pe, po, out, in);
    one_minus_deo_dag_in_vm_(pe, po, out, in, in0, fac);
    deo_dag_out_pos_vm_(pe, po, out, in, fac );
  }
  //----------------------------------------------------------------------------------------
  void mtilde_vm_(scd_t* out, scd_t* in){
    scd_t* tmp = (scd_t*)malloc( sizeof(scd_t) * nxd*ny*nz*nt);
    dee_deo_vm_(pco, pce, tmp, in);
    one_minus_dee_deo_vm_(pce, pco, out, tmp, in);
    free(tmp);
  }
  void mtilde_dag_vm_(scd_t* out, scd_t* in){
    scd_t* tmp0 = (scd_t*)malloc( sizeof(scd_t) * nxd*ny*nz*nt);
    scd_t* tmp1 = (scd_t*)malloc( sizeof(scd_t) * nxd*ny*nz*nt);
    clvd_vm2_(pce, tmp0, in);
    dee_deo_dag_vm_(pco, pce, tmp1, tmp0);
    one_minus_deo_dag_vm_(pce, pco, out, tmp1, in);
    free(tmp0);
    free(tmp1);
  }
  void mtdagmt_vm_(scd_t* out, scd_t* in){
    scd_t* tmp = (scd_t*)malloc( sizeof(scd_t) * nxd*ny*nz*nt);
    mtilde_vm_(tmp, in);
    mtilde_dag_vm_(out, tmp);
    free(tmp);
  }

  //---------------------------------------------------------------------------------------- single precision
  void deo_in_s_(int* pe, int* po, scs_t* out, scs_t* in);
  void dee_deo_in_s_(int* pe, int* po, scs_t* out, scs_t* in);
  void one_minus_dee_deo_in_s_(int* pe, int* po, scs_t* out, scs_t* in, scs_t* in0, float factor);
  void deo_out_pre_s_(const int* pe, const int* po, const scs_t* out, const scs_t* in);
  void deo_out_pos_s_(const int* pe, const int* po, const scs_t* out, const scs_t* in, float factor);
  void dee_deo_out_pos_s_(int* pe, int* po, scs_t* out, scs_t* in, float factor);
  //----------------------------------------------------------------------------------------
  void deo_s_(int* pe, int* po, scs_t* out, scs_t* in){
    float fac = 1;
    deo_out_pre_s_(pe, po, out, in);
    deo_in_s_(pe, po, out, in);
    deo_out_pos_s_(pe, po, out, in, fac);
  }
  void dee_deo_s_(int* pe, int* po, scs_t* out, scs_t* in){
    deo_out_pre_s_(pe, po, out, in);
    dee_deo_in_s_(pe, po, out, in);
    dee_deo_out_pos_s_(pe, po, out, in, 1);
  }
  void one_minus_dee_deo_s_(int* pe, int* po, scs_t* out, scs_t* in, scs_t* in0){
    float fac = -(float)kappa2;
    deo_out_pre_s_(pe, po, out, in);
    one_minus_dee_deo_in_s_(pe, po, out, in, in0, fac);
    dee_deo_out_pos_s_(pe, po, out, in, fac );
  }
  //----------------------------------------------------------------------------------------
  void mtilde_s_(scs_t* out, scs_t* in){
    scs_t* tmp = (scs_t*)malloc( sizeof(scs_t) * nxs*ny*nz*nt);
    dee_deo_s_(pco, pce, tmp, in);
    one_minus_dee_deo_s_(pce, pco, out, tmp, in);
    free(tmp);
  }

  //---------------------------------------------------------------------------------------- for DD solver
  void ddd_out_pre_d_(scd_t* in, int* domain);
  void ddd_out_pos_d_(scd_t* out, scd_t* in, int* domain);
  void ddd_in_d_(scd_t* out, scd_t* in, int* domain);

  void ddd_out_pre_s_(scs_t* in, int* domain);
  void ddd_out_pre_s_noprl_(scs_t* in, int* domain);
  void ddd_out_pre_s_no_timer_(scs_t* in, int* domain);
  void ddd_out_pre_s_noprl_no_timer_(scs_t* in, int* domain);
  void ddd_out_pos_s_(scs_t* out, scs_t* in, int* domain, float factor);
  void ddd_out_pos_s_no_timer_(scs_t* out, scs_t* in, int* domain, float factor);
  void ddd_in_s_(scs_t* out, scs_t* in, int* domain);
  void ddd_out_pos_s_noprl_(scs_t* out, scs_t* in, int* domain, float factor);
  void ddd_out_pos_s_noprl_no_timer_(scs_t* out, scs_t* in, int* domain, float factor);
  void ddd_in_s_noprl(scs_t* out, scs_t* in, int* domain);
  int domain0 = 0;
  int domain1 = 1;
  //----------------------------------------------------------------------------------------
  void ddd_d_(scd_t* out, scd_t* in){
    ddd_out_pre_d_(in, &domain0);
    ddd_in_d_(     &out[vold*0], &in[vold*0], &domain0);
    ddd_out_pos_d_(&out[vold*0], &in[vold*0], &domain0);
    ddd_out_pre_d_(in, &domain1);
    ddd_in_d_(     &out[vold*1], &in[vold*1], &domain1);
    ddd_out_pos_d_(&out[vold*1], &in[vold*0], &domain1);
  }
  //----------------------------------------------------------------------------------------
  void ddd_eo_s_(scs_t* out, scs_t* in, int* domain){
    
    ddd_out_pre_s_no_timer_(in, domain );
    ddd_in_s_(     &out[vols*(*domain)], &in[vols*(*domain)], domain);
    ddd_out_pos_s_no_timer_(&out[vols*(*domain)], in, domain, (float)mkappa);
  }
  //----------------------------------------------------------------------------------------
  void ddd_s_(scs_t* out, scs_t* in){
    ddd_eo_s_(out, in, &domain_e);
    ddd_eo_s_(out, in, &domain_o);
  }
  void jinv_ddd_in_s_(scs_t* x, scs_t* b, int *domain, int* maxiter);
  void jinv_ddd_in_s_noprl_(scs_t* x, scs_t* b, int *domain, int* maxiter);
  //----------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------
  inline void copy_vols2_s_noprl_(scs_t* out, scs_t* in){
#ifdef COMPILE_TIME_DIM_SIZE
    const int vols = VOLS;
#else
    const int vols = ::vols;
#endif
    int i, j;
#pragma omp for private(i,j)
    for(i=0; i<vols*2; i++){
      if (i+_PFI1 < volse*2) {
	__builtin_prefetch(&(((in+i+_PFI1 )->c_prefetch)[0]),0,3);
	__builtin_prefetch(&(((in+i+_PFI1 )->c_prefetch)[64]),0,3);
	__builtin_prefetch(&(((in+i+_PFI1 )->c_prefetch)[128]),0,3);
	__builtin_prefetch(&(((out+i+_PFI1 )->c_prefetch)[0]),1,3);
	__builtin_prefetch(&(((out+i+_PFI1 )->c_prefetch)[64]),1,3);
	__builtin_prefetch(&(((out+i+_PFI1 )->c_prefetch)[128]),1,3);
      }
      for(j=0; j<24; j++){
	for(int v=0; v<VLENS; v++){
	  out[i].ccs[j].v[v]  = in[i].ccs[j].v[v];
	}
      }
    }
  }
  //----------------------------------------------------------------------------------------
  inline void zero_vols_s_noprl_(scs_t* out){
#ifdef COMPILE_TIME_DIM_SIZE
    const int vols = VOLS;
#else
    const int vols = ::vols;
#endif
    int i, j;
#pragma omp  for private(i,j)
    for(i=0; i<vols; i++){
      if (i+_PFI4 < volse) {
	__builtin_prefetch(&((( out+i+_PFI4 )->c_prefetch)[0]),1,3);
	__builtin_prefetch(&((( out+i+_PFI4 )->c_prefetch)[64]),1,3);
	__builtin_prefetch(&((( out+i+_PFI4 )->c_prefetch)[128]),1,3);
      }
      for(j=0; j<24; j++){
	for(int v=0; v<VLENS; v++){
	  out[i].ccs[j].v[v] = 0;
	}
      }
    }
  }
  //----------------------------------------------------------------------------------------
  inline void accum_add_vols_s_noprl_(scs_t* out, scs_t* in){
#ifdef COMPILE_TIME_DIM_SIZE
    const int vols = VOLS;
#else
    const int vols = ::vols;
#endif
    int i, j;
#pragma omp for private(i,j)
    for(i=0; i<vols; i++){
      if (i+_PFI3 < volse) {
	__builtin_prefetch(&((( out+i+_PFI3 )->c_prefetch)[0]),1,3);
	__builtin_prefetch(&((( out+i+_PFI3 )->c_prefetch)[64]),1,3);
	__builtin_prefetch(&((( out+i+_PFI3 )->c_prefetch)[128]),1,3);
	__builtin_prefetch(&(((  in+i+_PFI3 )->c_prefetch)[0]),0,3);
	__builtin_prefetch(&(((  in+i+_PFI3 )->c_prefetch)[64]),0,3);
	__builtin_prefetch(&(((  in+i+_PFI3 )->c_prefetch)[128]),0,3);
      }
      for(j=0; j<24; j++){
	for(int v=0; v<VLENS; v++){
	  out[i].ccs[j].v[v] += in[i].ccs[j].v[v];
	}
      }
    }
  }
  //----------------------------------------------------------------------------------------
  inline void accum_sub_vols_s_noprl_(scs_t* out, scs_t* in){
#ifdef COMPILE_TIME_DIM_SIZE
    const int vols = VOLS;
#else
    const int vols = ::vols;
#endif
    int i, j;
#pragma omp for private(i,j)
    for(i=0; i<vols; i++){
      if (i+_PFI3 < volse) {
	__builtin_prefetch(&((( out+i+_PFI3 )->c_prefetch)[0]),1,3);
	__builtin_prefetch(&((( out+i+_PFI3 )->c_prefetch)[64]),1,3);
	__builtin_prefetch(&((( out+i+_PFI3 )->c_prefetch)[128]),1,3);
	__builtin_prefetch(&(((  in+i+_PFI3 )->c_prefetch)[0]),0,3);
	__builtin_prefetch(&(((  in+i+_PFI3 )->c_prefetch)[64]),0,3);
	__builtin_prefetch(&(((  in+i+_PFI3 )->c_prefetch)[128]),0,3);
      }
      for(j=0; j<24; j++){
	for(int v=0; v<VLENS; v++){
	  out[i].ccs[j].v[v] -= in[i].ccs[j].v[v];
	}
      }
    }
  }
  //----------------------------------------------------------------------------------------
  inline void accum_addsub_vols_s_noprl_(scs_t* out, scs_t* in1, scs_t* in2){
#ifdef COMPILE_TIME_DIM_SIZE
    const int vols = VOLS;
#else
    const int vols = ::vols;
#endif
    int i, j;
#pragma omp for private(i,j)
    for(i=0; i<vols; i++){
      if (i+_PFI2 < volse) {
	__builtin_prefetch(&(((in1+i+_PFI2 )->c_prefetch)[0]),0,3);
	__builtin_prefetch(&(((in1+i+_PFI2 )->c_prefetch)[64]),0,3);
	__builtin_prefetch(&(((in1+i+_PFI2 )->c_prefetch)[128]),0,3);
	__builtin_prefetch(&(((in2+i+_PFI2 )->c_prefetch)[0]),0,3);
	__builtin_prefetch(&(((in2+i+_PFI2 )->c_prefetch)[64]),0,3);
	__builtin_prefetch(&(((in2+i+_PFI2 )->c_prefetch)[128]),0,3);
	__builtin_prefetch(&(((out+i+_PFI2 )->c_prefetch)[0]),1,3);
	__builtin_prefetch(&(((out+i+_PFI2 )->c_prefetch)[64]),1,3);
	__builtin_prefetch(&(((out+i+_PFI2 )->c_prefetch)[128]),1,3);
      }
      for(j=0; j<24; j++){
	for(int v=0; v<VLENS; v++){
	  out[i].ccs[j].v[v] += in1[i].ccs[j].v[v] - in2[i].ccs[j].v[v];
	}
      }
    }
  }
  //----------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------
  void prec_ddd_s_(scs_t* out, scs_t* in, int* nsap, int* nm){
    static scs_t *s, *q;
    int ret;
    ret=0;
    if(s==0) {
      void *p;
      ret = posix_memalign(&p, CLS, sizeof(scs_t)*vols*2);
      s = (scs_t*)p;
    }
    if(q==0) {
      void *p;
      ret = posix_memalign(&p, CLS, sizeof(scs_t)*vols*2);
      q = (scs_t*)p;
    }

#ifdef COMPILE_TIME_DIM_SIZE
    const int vols = VOLS;
#else
    const int vols = ::vols;
#endif
    float kappa = ::kappa, mkappa = ::mkappa;

    if (*nsap > 10) { printf("prec_ddd_s_ nsap > 10\n"); exit(1); }
    if (*nm   > 10) { printf("prec_ddd_s_ nn   > 10\n"); exit(1); }

    _PREC_DDD_S_TIC_;
#pragma omp parallel 
    {
      copy_vols2_s_noprl_(s, in);
      if (npe[1]==1 || npe[2]==1 || npe[3]==1){
      	zero_vols_s_noprl_(&q[vols*domain_o]);
      }
      _OTHER_CALC_TOC_;
      for (int isap=0; isap < *nsap; isap++) {
	_JINV_DDD_IN_S_TIC_;
	jinv_ddd_in_s_noprl_(&q[vols*domain_e], &s[vols*domain_e], &domain_e, nm);
	_JINV_DDD_IN_S_TOC_;
	_DDD_OUT_PRE_S_TIC_;
	ddd_out_pre_s_noprl_(   q, &domain_o);
	_DDD_OUT_PRE_S_TOC_;
	_DDD_IN_S_TIC_;
	ddd_in_s_noprl(&q[vols*domain_o], &q[vols*domain_e], &domain_e);
	_DDD_IN_S_TOC_;
	_OTHER_CALC_TIC_;
	accum_addsub_vols_s_noprl_(&s[vols*domain_e], &in[vols*domain_e],  &q[vols*domain_o]);
	_OTHER_CALC_TOC_;
	_DDD_OUT_POS_S_TIC_;
	ddd_out_pos_s_noprl_(&s[vols*domain_o], q, &domain_o, (float)kappa);
	_DDD_OUT_POS_S_TOC_;

	_JINV_DDD_IN_S_TIC_;
	jinv_ddd_in_s_noprl_(&q[vols*domain_o], &s[vols*domain_o], &domain_o, nm);
	_JINV_DDD_IN_S_TOC_;
	_DDD_OUT_PRE_S_TIC_;
	ddd_out_pre_s_noprl_(   q, &domain_e);
	_DDD_OUT_PRE_S_TOC_;
	_DDD_IN_S_TIC_;
	ddd_in_s_noprl(&q[vols*domain_e], &q[vols*domain_o], &domain_o);
	_DDD_IN_S_TOC_;
	_OTHER_CALC_TIC_;
	accum_addsub_vols_s_noprl_(&s[vols*domain_o], &in[vols*domain_o],  &q[vols*domain_e]);
	_OTHER_CALC_TOC_;
	_DDD_OUT_POS_S_TIC_;
	ddd_out_pos_s_noprl_(&s[vols*domain_e], q, &domain_e, (float)kappa);
	_DDD_OUT_POS_S_TOC_;
      }
      if (npe[1]==1 || npe[2]==1 || npe[3]==1){
      	_OTHER_CALC_TIC_;
      	zero_vols_s_noprl_(&q[vols*domain_o]);
      	_OTHER_CALC_TOC_;
      }

      _JINV_DDD_IN_S_TIC_;
      jinv_ddd_in_s_noprl_(&q[vols*domain_e], &s[vols*domain_e], &domain_e, nm);
      _JINV_DDD_IN_S_TOC_;

      _DDD_OUT_PRE_S_TIC_;
      ddd_out_pre_s_noprl_(   q, &domain_o);
      _DDD_OUT_PRE_S_TOC_;

      _DDD_IN_S_TIC_;
      ddd_in_s_noprl(        &out[vols*domain_e], &q[vols*domain_e], &domain_e);
      _DDD_IN_S_TOC_;

      if (npe[1]==1 || npe[2]==1 || npe[3]==1){
      	_DDD_OUT_POS_S_TIC_;
      	ddd_out_pos_s_noprl_(   &s[vols*domain_o], q,  &domain_o, (float)kappa);
      	_DDD_OUT_POS_S_TOC_;
      
      	_JINV_DDD_IN_S_TIC_;
      	jinv_ddd_in_s_noprl_(&q[vols*domain_o], &s[vols*domain_o], &domain_o, nm);
      	_JINV_DDD_IN_S_TOC_;
      
      	_DDD_OUT_PRE_S_TIC_;
      	ddd_out_pre_s_noprl_(  q , &domain_e);
      	_DDD_OUT_PRE_S_TOC_;
      
      	_DDD_IN_S_TIC_;
      	ddd_in_s_noprl( &out[vols*domain_o], &q[vols*domain_o], &domain_o);
      	_DDD_IN_S_TOC_;
      
      	_DDD_OUT_POS_S_TIC_;
      	ddd_out_pos_s_noprl_(&out[vols*domain_e], q, &domain_e, (float)mkappa);
      	_DDD_OUT_POS_S_TOC_;
      
      	_DDD_OUT_PRE_S_TIC_;
      	ddd_out_pre_s_noprl_(  q , &domain_o);
      	_DDD_OUT_PRE_S_TOC_;
      
      	_DDD_OUT_POS_S_TIC_;
      	ddd_out_pos_s_noprl_(&out[vols*domain_o], q, &domain_o, (float)mkappa);
      	_DDD_OUT_POS_S_TOC_;
      
      }else{
	_OTHER_CALC_TIC_;
	zero_vols_s_noprl_(&out[vols*domain_o]);
	_OTHER_CALC_TOC_;

	_DDD_OUT_POS_S_TIC_;
	ddd_out_pos_s_noprl_(   &out[vols*domain_o], q,  &domain_o, (float)mkappa);
	_DDD_OUT_POS_S_TOC_;

	_OTHER_CALC_TIC_;
	accum_sub_vols_s_noprl_(&s[vols*domain_o], &out[vols*domain_o]);
	_OTHER_CALC_TOC_;

	_JINV_DDD_IN_S_TIC_;
	jinv_ddd_in_s_noprl_(&q[vols*domain_o], &s[vols*domain_o], &domain_o, nm);
	_JINV_DDD_IN_S_TOC_;

	_DDD_OUT_PRE_S_TIC_;
	ddd_out_pre_s_noprl_(  q , &domain_e);
	_DDD_OUT_PRE_S_TOC_;

	_DDD_IN_S_TIC_;
	ddd_in_s_noprl( &s[vols*domain_o], &q[vols*domain_o], &domain_o);
	_DDD_IN_S_TOC_;

	_OTHER_CALC_TIC_;
	accum_add_vols_s_noprl_(&out[vols*domain_o], &s[vols*domain_o]);
	_OTHER_CALC_TOC_;
	
	_DDD_OUT_POS_S_TIC_;
	ddd_out_pos_s_noprl_(&out[vols*domain_e], q, &domain_e, (float)mkappa);
	_DDD_OUT_POS_S_TOC_;
      }

      _OTHER_CALC_TIC_;
    } //end parallel
    _PREC_DDD_S_TOC_;
  }

  //----------------------------------------------------------------------------------------
  void prec_s_(scs_t* out, scs_t* in, int* nsap, int* nm){
    static scs_t *s, *q;
    int ret;
    ret=0;
    if(s==0) {
      void *p;
      ret = posix_memalign(&p, CLS, sizeof(scs_t)*vols*2);
      s = (scs_t*)p;
    }
    if(q==0) {
      void *p;
      ret = posix_memalign(&p, CLS, sizeof(scs_t)*vols*2);
      q = (scs_t*)p;
    }

#ifdef COMPILE_TIME_DIM_SIZE
    const int vols = VOLS;
#else
    const int vols = ::vols;
#endif
    float kappa = ::kappa;

    _PREC_S_TIC_;
#pragma omp parallel 
    {
      copy_vols2_s_noprl_(s, in);
      if (npe[1]==1 || npe[2]==1 || npe[3]==1)zero_vols_s_noprl_(&q[vols*domain_o]);
      for (int isap=0; isap < *nsap; isap++) {
	jinv_ddd_in_s_noprl_(&q[vols*domain_e], &s[vols*domain_e], &domain_e, nm);
	ddd_out_pre_s_noprl_no_timer_(   q, &domain_o);
	ddd_in_s_noprl(&q[vols*domain_o], &q[vols*domain_e], &domain_e);
        accum_addsub_vols_s_noprl_(&s[vols*domain_e], &in[vols*domain_e],  &q[vols*domain_o]);
	ddd_out_pos_s_noprl_no_timer_(&s[vols*domain_o], q, &domain_o, (float)kappa);
	jinv_ddd_in_s_noprl_(&q[vols*domain_o], &s[vols*domain_o], &domain_o, nm);
	ddd_out_pre_s_noprl_no_timer_(   q, &domain_e);
	ddd_in_s_noprl(&q[vols*domain_e], &q[vols*domain_o], &domain_o);
	accum_addsub_vols_s_noprl_(&s[vols*domain_o], &in[vols*domain_o],  &q[vols*domain_e]);
	ddd_out_pos_s_noprl_no_timer_(&s[vols*domain_e], q, &domain_e, (float)kappa);
      }
      if (npe[1]==1 || npe[2]==1 || npe[3]==1)zero_vols_s_noprl_(&out[vols*domain_o]);
      jinv_ddd_in_s_noprl_(&out[vols*domain_e], &s[vols*domain_e], &domain_e, nm);
      ddd_out_pre_s_noprl_no_timer_(out, &domain_o);
      ddd_out_pos_s_noprl_no_timer_(&s[vols*domain_o], out, &domain_o, (float)kappa);
      jinv_ddd_in_s_noprl_(&out[vols*domain_o], &s[vols*domain_o], &domain_o, nm);
    } //end parallel
    _PREC_S_TOC_;
  }

//  //----------------------------------------------------------------------------------------
//  void copy_s(float* out, float* in, int n){
//    int i;
//#pragma omp parallel for private(i)
//    for (i=0; i<n; i++){
//      out[i] = in[i];
//    }
//  }
//  //----------------------------------------------------------------------------------------
//  void zero_s(float* out, int n){
//    int i;
//#pragma omp parallel for private(i)
//    for (i=0; i<n; i++){
//      out[i] = 0;
//    }
//  }
//  //----------------------------------------------------------------------------------------
//  void accum_addsub_s(float* out, float* in1, float* in2, int n){
//    int i;
//#pragma omp parallel for private(i)
//    for (i=0; i<n; i++){
//      out[i] = in1[i] - in2[i] + out[i];
//    }
//  }
//  //----------------------------------------------------------------------------------------
//  void prec_part_s(scs_t* tmp0, scs_t* tmp1, scs_t* in, int* nm, int e){
//    int o = 1 -e;
//    jinv_ddd_in_s_(&tmp1[vols*e], &tmp0[vols*e], &e, nm);
//    ddd_out_pre_s_no_timer_(   tmp1, &o);
//    ddd_in_s_(&tmp1[vols*o], &tmp1[vols*e], &e); // tmp1_o can be used
//    accum_addsub_s((float*)&tmp0[vols*e], (float*)&in[vols*e],  (float*)&tmp1[vols*o], vols*24*VLENS);
//    ddd_out_pos_s_no_timer_(&tmp0[vols*o], tmp1, &o, (float)kappa);
//  }
//  //----------------------------------------------------------------------------------------
//  void prec_ddd_s_1(scs_t* out, scs_t* in, int* nsap, int* nm){
//    static scs_t *tmp;
//    if(tmp==0) tmp = (scs_t*)malloc( sizeof(scs_t) * vols*2);
//    prec_s_(tmp, in, nsap, nm);
//    ddd_s_(out, tmp);
//  }
//  //----------------------------------------------------------------------------------------
//  void prec_ddd_s_2(scs_t* out, scs_t* in, int* nsap, int* nm){
//    int i, j;
//    static scs_t *s, *q;
//    if(s==0) s = (scs_t*)malloc( sizeof(scs_t) * vols*2);
//    if(q==0) q = (scs_t*)malloc( sizeof(scs_t) * vols*2);
//
//    if (*nsap > 10) { printf("prec_ddd_s_ nsap > 10\n"); exit(1); }
//    if (*nm   > 10) { printf("prec_ddd_s_ nn   > 10\n"); exit(1); }
//
//    copy_s((float*)s, (float*)in, vols*2*24*VLENS);
//
//    if (npe[1]==1 || npe[2]==1 || npe[3]==1)zero_s((float*)&q[vols*domain_o], vols*24*VLENS);
//    for (int isap=0; isap < *nsap; isap++) {
//      prec_part_s(s, q, in, nm, domain_e);
//      prec_part_s(s, q, in, nm, domain_o);
//    }
//    if (npe[1]==1 || npe[2]==1 || npe[3]==1)zero_s((float*)&q[vols*domain_o], vols*24*VLENS);
//    jinv_ddd_in_s_(&q[vols*domain_e], &s[vols*domain_e], &domain_e, nm); // qe = Aee se
//    ddd_out_pre_s_no_timer_(   q, &domain_o);                            //      send Doe qe
//    ddd_in_s_(     &out[vols*domain_e], &q[vols*domain_e], &domain_e);   // xe = Dee qe
//
//    if (npe[1]==1 || npe[2]==1 || npe[3]==1){
//      ddd_out_pos_s_no_timer_(&s[vols*domain_o], q, &domain_o, (float)kappa);    // so = so - Doe qe recv Doe qe
//      jinv_ddd_in_s_(&q[vols*domain_o], &s[vols*domain_o], &domain_o, nm);       // qo = Aoo so
//      ddd_out_pre_s_no_timer_(q, &domain_e );                                    //                  send Deo qo
//      ddd_in_s_(     &out[vols*domain_o], &q[vols*domain_o], &domain_o);         // xo = Doo qo
//      ddd_out_pos_s_no_timer_(&out[vols*domain_e], q, &domain_e, (float)mkappa); // xe = xe + Deo qo recv Deo qo
//      ddd_out_pre_s_no_timer_(q, &domain_o );                                    //                  send Doe qe
//      ddd_out_pos_s_no_timer_(&out[vols*domain_o], q, &domain_o, (float)mkappa); // xo = xo + Doe qe recv Doe qe
//    }else{
//      zero_s((float*)&out[vols*domain_o], vols*24*VLENS);
//      ddd_out_pos_s_no_timer_(&out[vols*domain_o], q, &domain_o, (float)mkappa); // q = - Doe qe recv Doe qe
//#pragma omp parallel for private(i, j)
//      for(i=0; i<vols; i++){
//	for(j=0; j<24; j++){
//	  for(int v=0; v<VLENS; v++){
//	    s[i+vols*domain_o].ccs[j].v[v] -=  out[i+vols*domain_o].ccs[j].v[v];
//	  }
//	}
//      }
//      jinv_ddd_in_s_(&q[vols*domain_o], &s[vols*domain_o], &domain_o, nm);    // qo = Aoo so
//      ddd_out_pre_s_no_timer_(q, &domain_e );                                 // send Deo qo
//      ddd_in_s_(     &s[vols*domain_o], &q[vols*domain_o], &domain_o);        // xo = Doo qo
//#pragma omp parallel for private(i, j)
//      for(i=0; i<vols; i++){
//	for(j=0; j<24; j++){
//	  for(int v=0; v<VLENS; v++){
//	    out[i+vols*domain_o].ccs[j].v[v] += s[i+vols*domain_o].ccs[j].v[v];
//	  }
//	}
//      }
//      ddd_out_pos_s_no_timer_(&out[vols*domain_e], q, &domain_e, (float)mkappa); // xe = xe + Deo qo
//    }
//  }

#ifdef __cplusplus
}
#endif
