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
#include "wilson_s.h"
#include "clover_s.h"

#ifdef __cplusplus
extern "C"{
#endif


  extern __attribute__((aligned(64))) pglus_t glus;
  extern __attribute__((aligned(64))) pclvs_t clvs;
  extern   projscs1_t *xfs_send, *xfs_recv;
  extern   projscs1_t *xbs_send, *xbs_recv;
  extern   projscs_t *yfs_send, *yfs_recv;
  extern   projscs_t *ybs_send, *ybs_recv;
  extern   projscs_t *zfs_send, *zfs_recv;
  extern   projscs_t *zbs_send, *zbs_recv;
  extern   projscs_t *tfs_send, *tfs_recv;
  extern   projscs_t *tbs_send, *tbs_recv;
#ifdef _MPI_
  extern  MPI_Request ss_req[8];
  extern  MPI_Request rs_req[8];
#endif
  
  //  nh=nx/2, nd=nx/2/8, nxs=nx/2/16
  extern int nxs, ny, nz, nt, px, py, pz, pt;
  extern int vols;// = nxs * ny * nz * nt;
  extern int pt;// process coordinate
  extern int npe[4];// number of proccesses
  extern double fbc[4][2];// fermion boundary condiion
  
  void xbound(int req, int prec);
  void xbound_wait(int req, int prec);
  void xbound_send_waitall(int prec);
  void xbound_recv_updateall(int prec);
  void xbound_recv_update(int req, int prec);
  void xbound_recv_waitall(int prec);
  void xbound_recv_okall(int prec);
  void xbound_flip_parity(int prec);


  
  //---------------------------------------------------------------------------------------- preprocess mult D for boundary
  void deo_out_pre_s_(int* pe, int* po, scs_t* out, scs_t* in){
    int x, y, z, t, i0, i1, i2, i3, i4, i5, i6, i7;
    pglus_t gxb = &glus[vols*0 + NDIM*vols*(*po)];
    pglus_t gyb = &glus[vols*1 + NDIM*vols*(*po)];
    pglus_t gzb = &glus[vols*2 + NDIM*vols*(*po)];
    pglus_t gtb = &glus[vols*3 + NDIM*vols*(*po)];
    int jjj= ((*pe) + ny*py + nz*pz + nt*pt)%2;
#pragma omp for private(t,z,y,x,i0,i1,i2,i3,i4,i5,i6,i7) collapse(3) 
    for(t=0; t<nt; t++){
      for(z=0; z<nz; z++){
	for(y=0; y<ny; y++){
          int ick= ( jjj + y     + z + t)%2;
	  int it = nxs*ny*nz*t;
	  int iz = nxs*ny*z;
	  int iy = nxs*y;
	  for(x=0; x<nxs; x++){
            int ixf = 0     + iy + iz + it;
            int ixb = nxs-1 + iy + iz + it;
	    int iyf = x + iz + it;
	    int iyb = x + iz + it + (ny-1)*nxs;
	    int izf = x + iy + it;
	    int izb = x + iy + it + (nz-1)*nxs*ny;
	    int itf = x + iy + iz;
	    int itb = x + iy + iz + (nt-1)*nxs*ny*nz;
            if(x == nxs-1 && ick == 1) {i0 = y + ny*z + ny*nz*t;__l_xf_s(xfs_send[i0].c, in, ixf);}
	    if(y == ny-1             ) {i2 = x +nxs*z +nxs*nz*t;__l_yf_s(yfs_send[i2].c, in, iyf);}
	    if(z == nz-1             ) {i4 = x +nxs*y +nxs*ny*t;__l_zf_s(zfs_send[i4].c, in, izf);}
	    if(t == nt-1             ) {i6 = x +nxs*y +nxs*ny*z;__l_tf_s(tfs_send[i6].c, in, itf);}
            if(x == 0     && ick == 0) {i1 = y + ny*z + ny*nz*t;__pre_xbwd_s(xbs_send[i1].c, gxb, in, ixb);}
	    if(y == 0                ) {i3 = x +nxs*z +nxs*nz*t;__pre_ybwd_s(ybs_send[i3].c, gyb, in, iyb);}
	    if(z == 0                ) {i5 = x +nxs*y +nxs*ny*t;__pre_zbwd_s(zbs_send[i5].c, gzb, in, izb);}
	    if(t == 0                ) {i7 = x +nxs*y +nxs*ny*z;__pre_tbwd_s(tbs_send[i7].c, gtb, in, itb);}
	  }
	}
      }
    }
#ifdef _THREADED_RDMA
#pragma omp for nowait
  for(int req=0; req<8; req++) {
    xbound_recv_update(req, 4); // update the pointer to buffer before memcpy
    xbound(req,4);
  }
#else // _THREADED_RDMA
#ifdef _NO_OMP_SINGLE
#pragma omp master
#else
#pragma omp single nowait
#endif // _NO_OMP_SINGLE
    {
      xbound_recv_updateall(4);  // update the pointers to recv. buffer
      xbound(0,4);
      xbound(1,4);
      xbound(2,4);
      xbound(3,4);
      xbound(4,4);
      xbound(5,4);
      xbound(6,4);
      xbound(7,4);
    }
#endif    // _THREADED_RDMA
  }
  //---------------------------------------------------------------------------------------- postprocess mult D for boundary
  void deo_out_pos_s_(int* pe, int* po, scs_t* out, scs_t* in, float factor){
    __attribute__((aligned(64)))  scs_t tmp;
    int x, y, z, t, i2, i3, i4, i5, i6, i7;
    pglus_t gxf = &glus[vols*0 + NDIM*vols*(*pe)];
    pglus_t gyf = &glus[vols*1 + NDIM*vols*(*pe)];
    pglus_t gzf = &glus[vols*2 + NDIM*vols*(*pe)];
    pglus_t gtf = &glus[vols*3 + NDIM*vols*(*pe)];
    int jjj= ((*pe) + ny*py + nz*pz + nt*pt)%2;

#if 0 //#ifdef _THREADED_RDMA
#pragma omp for
  for(int req=0; req<8; req++) {
    xbound_wait(req,4);
  }
#else // #ifdef _THREADED_RDMA
#ifdef _NO_OMP_SINGLE
#pragma omp master
#else
#pragma omp single
#endif // _NO_OMP_SINGLE
    {
      xbound_wait(0,4);
      xbound_wait(1,4);
      xbound_wait(2,4);
      xbound_wait(3,4);
      xbound_wait(4,4);
      xbound_wait(5,4);
      xbound_wait(6,4);
      xbound_wait(7,4);
    }
#ifdef _NO_OMP_SINGLE
#pragma omp barrier
#endif
#endif // #ifdef _THREADED_RDMA

#pragma omp for private(t,z,y,x,tmp,i2,i3,i4,i5,i6,i7) collapse(3) 
    for(t=0; t<nt; t++){
      for(z=0; z<nz; z++){
	for(y=0; y<ny; y++){
	  int ick= ( jjj + y     + z + t)%2;
	  int ix = y + ny*z + ny*nz*t;
	  for(x=0; x<nxs; x++){
	    if(x == nxs-1 || y == ny-1 || z == nz-1 || t == nt-1 || x == 0 || y == 0 || z == 0 || t == 0 ){
	      int i0  = x + nxs*y + nxs*ny*z + nxs*ny*nz*t;
	      __zero_sc_s(tmp.c);
	      if(x == nxs-1 && ick == 1) {                        __pos_xfwd_s(tmp.c, gxf, xfs_recv[ix], i0);}
	      if(y == ny-1             ) {i2 = x +nxs*z +nxs*nz*t;__pos_yfwd_s(tmp.c, gyf, yfs_recv[i2], i0);}
	      if(z == nz-1             ) {i4 = x +nxs*y +nxs*ny*t;__pos_zfwd_s(tmp.c, gzf, zfs_recv[i4], i0);}
	      if(t == nt-1             ) {i6 = x +nxs*y +nxs*ny*z;__pos_tfwd_bc_s(tmp.c, gtf, tfs_recv[i6], i0, fbc[3][0]);}
	      if(x == 0     && ick == 0) {                        __s_xb1_s(tmp.c, xbs_recv[ix].c);}
	      if(y == 0                ) {i3 = x +nxs*z +nxs*nz*t;__s_yb_s(tmp.c, ybs_recv[i3].c);}
	      if(z == 0                ) {i5 = x +nxs*y +nxs*ny*t;__s_zb_s(tmp.c, zbs_recv[i5].c);}
	      if(t == 0                ) {i7 = x +nxs*y +nxs*ny*z;__s_tb_bc_s(tmp.c, tbs_recv[i7].c, fbc[3][1]);}
	      __store_add_sc_s(out[i0].c, tmp.c, factor);
	    }
	  }
	}
      }
    }
#ifdef _NO_OMP_SINGLE
#pragma omp master
#else
#pragma omp single
#endif // _NO_OMP_SINGLE
    {
      xbound_send_waitall(4);
      xbound_recv_okall(4);
      xbound_flip_parity(4);
    }
#ifdef _NO_OMP_SINGLE
#pragma omp barrier
#endif
  }
  //---------------------------------------------------------------------------------------- postprocess mult D for boundary with clv mult
  void dee_deo_out_pos_s_(int* pe, int* po, scs_t* out, scs_t* in, float factor){
    __attribute__((aligned(64)))  scs_t tmp;
    int x, y, z, t, i2, i3, i4, i5, i6, i7;
    pglus_t gxf = &glus[vols*0 + NDIM*vols*(*pe)];
    pglus_t gyf = &glus[vols*1 + NDIM*vols*(*pe)];
    pglus_t gzf = &glus[vols*2 + NDIM*vols*(*pe)];
    pglus_t gtf = &glus[vols*3 + NDIM*vols*(*pe)];
    int jjj= ((*pe) + ny*py + nz*pz + nt*pt)%2;
#ifdef _NO_OMP_SINGLE
#pragma omp master
#else
#pragma omp single
#endif // _NO_OMP_SINGLE
    {
      xbound_wait(0,4);
      xbound_wait(1,4);
      xbound_wait(2,4);
      xbound_wait(3,4);
      xbound_wait(4,4);
      xbound_wait(5,4);
      xbound_wait(6,4);
      xbound_wait(7,4);
    }
#ifdef _NO_OMP_SINGLE
#pragma omp barrier
#endif // _NO_OMP_SINGLE

#pragma omp for private(t,z,y,x,tmp,i2,i3,i4,i5,i6,i7) collapse(3) 
    for(t=0; t<nt; t++){
      for(z=0; z<nz; z++){
	for(y=0; y<ny; y++){
	  int ick= ( jjj + y     + z + t)%2;
	  int ix = y + ny*z + ny*nz*t;
	  for(x=0; x<nxs; x++){
	    if(x == nxs-1 || y == ny-1 || z == nz-1 || t == nt-1 || x == 0 || y == 0 || z == 0 || t == 0 ){
	      int i0  = x + nxs*y + nxs*ny*z + nxs*ny*nz*t;
	      __zero_sc_s(tmp.c);
	      if(x == nxs-1 && ick == 1) {                        __pos_xfwd_s(tmp.c, gxf, xfs_recv[ix], i0);}
	      if(y == ny-1             ) {i2 = x +nxs*z +nxs*nz*t;__pos_yfwd_s(tmp.c, gyf, yfs_recv[i2], i0);}
	      if(z == nz-1             ) {i4 = x +nxs*y +nxs*ny*t;__pos_zfwd_s(tmp.c, gzf, zfs_recv[i4], i0);}
	      if(t == nt-1             ) {i6 = x +nxs*y +nxs*ny*z;__pos_tfwd_bc_s(tmp.c, gtf, tfs_recv[i6], i0, fbc[3][0]);}
	      if(x == 0     && ick == 0) {                        __s_xb1_s(tmp.c, xbs_recv[ix].c);}
	      if(y == 0                ) {i3 = x +nxs*z +nxs*nz*t;__s_yb_s(tmp.c, ybs_recv[i3].c);}
	      if(z == 0                ) {i5 = x +nxs*y +nxs*ny*t;__s_zb_s(tmp.c, zbs_recv[i5].c);}
	      if(t == 0                ) {i7 = x +nxs*y +nxs*ny*z;__s_tb_bc_s(tmp.c, tbs_recv[i7].c, fbc[3][1]);}
	      __mult_clvs( tmp.cv, clvs[i0 + vols*(*pe)].cv);
	      __store_add_sc_s(out[i0].c, tmp.c, factor);
	    }
	  }
	}
      }
    }
#ifdef _NO_OMP_SINGLE
#pragma omp master
#else
#pragma omp single
#endif // _NO_OMP_SINGLE
    {
      xbound_send_waitall(4);
      xbound_recv_okall(4);
      xbound_flip_parity(4);
    }
#ifdef _NO_OMP_SINGLE
#pragma omp barrier
#endif
  }
#ifdef __cplusplus
}
#endif
