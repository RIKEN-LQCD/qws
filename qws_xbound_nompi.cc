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
//#define min(a,b) (a)>(b)?(b):(a) 

//#include "timing.h"
//#include "eml_lib.h"

#ifdef __cplusplus
extern "C"{
#endif

  extern   projscd1_t *xfd_send;
  extern   projscd1_t *xfd_recv;
  extern   projscd1_t *xbd_send;
  extern   projscd1_t *xbd_recv;
  extern   projscd_t *yfd_send;
  extern   projscd_t *yfd_recv;
  extern   projscd_t *ybd_send;
  extern   projscd_t *ybd_recv;
  extern   projscd_t *zfd_send;
  extern   projscd_t *zfd_recv;
  extern   projscd_t *zbd_send;
  extern   projscd_t *zbd_recv;
  extern   projscd_t *tfd_send;
  extern   projscd_t *tfd_recv;
  extern   projscd_t *tbd_send;
  extern   projscd_t *tbd_recv;

  extern   projscs1_t *xfs_send;
  extern   projscs1_t *xfs_recv;
  extern   projscs1_t *xbs_send;
  extern   projscs1_t *xbs_recv;
  extern   projscs_t *yfs_send;
  extern   projscs_t *yfs_recv;
  extern   projscs_t *ybs_send;
  extern   projscs_t *ybs_recv;
  extern   projscs_t *zfs_send;
  extern   projscs_t *zfs_recv;
  extern   projscs_t *zbs_send;
  extern   projscs_t *zbs_recv;
  extern   projscs_t *tfs_send;
  extern   projscs_t *tfs_recv;
  extern   projscs_t *tbs_send;
  extern   projscs_t *tbs_recv;

  extern   int npe[4];
  extern   int nx;
  extern   int ny;
  extern   int nz;
  extern   int nt;
  extern   int nxh;
  extern   int nxd;
  extern   int nxs;

  projscs1_t *xfs_recv_array[2];
  projscs1_t *xbs_recv_array[2];
  projscs_t *yfs_recv_array[2];
  projscs_t *ybs_recv_array[2];
  projscs_t *zfs_recv_array[2];
  projscs_t *zbs_recv_array[2];
  projscs_t *tfs_recv_array[2];
  projscs_t *tbs_recv_array[2];
  int bufs_parity;


  void xbound_set_parity(int parity, int prec){
    if(prec == 8){
      return;
    }
    bufs_parity=parity;
  }

  void xbound_flip_parity(int prec){
    if(prec == 8){
      return;
    }
    xbound_set_parity(1-bufs_parity, prec);
  }

  void xbound_recv_updateall(int prec){
    if(prec==4){
      xfs_recv=xfs_recv_array[bufs_parity];
      xbs_recv=xbs_recv_array[bufs_parity];
      yfs_recv=yfs_recv_array[bufs_parity];
      ybs_recv=ybs_recv_array[bufs_parity];
      zfs_recv=zfs_recv_array[bufs_parity];
      zbs_recv=zbs_recv_array[bufs_parity];
      tfs_recv=tfs_recv_array[bufs_parity];
      tbs_recv=tbs_recv_array[bufs_parity];
    }
  }

  //---------------------------------------------------------------------------------------- init communication
  void xbound_init(int rank,
		   int pxf, int pyf, int pzf, int ptf,
		   int pxb, int pyb, int pzb, int ptb)
  {


    // allocate communication buffers (double prec.)
    xfd_send = (projscd1_t*)malloc( sizeof(projscd1_t) * ny *nz*nt);
    xbd_send = (projscd1_t*)malloc( sizeof(projscd1_t) * ny *nz*nt);
    yfd_send = (projscd_t*)malloc( sizeof(projscd_t) * nxd*nz*nt);
    ybd_send = (projscd_t*)malloc( sizeof(projscd_t) * nxd*nz*nt);
    zfd_send = (projscd_t*)malloc( sizeof(projscd_t) * nxd*ny*nt);
    zbd_send = (projscd_t*)malloc( sizeof(projscd_t) * nxd*ny*nt);
    tfd_send = (projscd_t*)malloc( sizeof(projscd_t) * nxd*ny*nz);
    tbd_send = (projscd_t*)malloc( sizeof(projscd_t) * nxd*ny*nz);
    xfd_recv = (projscd1_t*)malloc( sizeof(projscd1_t) * ny *nz*nt);
    xbd_recv = (projscd1_t*)malloc( sizeof(projscd1_t) * ny *nz*nt);
    yfd_recv = (projscd_t*)malloc( sizeof(projscd_t) * nxd*nz*nt);
    ybd_recv = (projscd_t*)malloc( sizeof(projscd_t) * nxd*nz*nt);
    zfd_recv = (projscd_t*)malloc( sizeof(projscd_t) * nxd*ny*nt);
    zbd_recv = (projscd_t*)malloc( sizeof(projscd_t) * nxd*ny*nt);
    tfd_recv = (projscd_t*)malloc( sizeof(projscd_t) * nxd*ny*nz);
    tbd_recv = (projscd_t*)malloc( sizeof(projscd_t) * nxd*ny*nz);

    // allocate communication buffers (single prec.)
    void *tmp;
    posix_memalign(&tmp, CLS, sizeof(projscs1_t)*ny*nz*nt);
    xfs_send = (projscs1_t*)tmp;
    posix_memalign(&tmp, CLS, sizeof(projscs1_t)*ny*nz*nt);
    xbs_send = (projscs1_t*)tmp;
    posix_memalign(&tmp, CLS, sizeof(projscs_t)*nxs*nz*nt);
    yfs_send = (projscs_t*)tmp;
    posix_memalign(&tmp, CLS, sizeof(projscs_t)*nxs*nz*nt);
    ybs_send = (projscs_t*)tmp;
    posix_memalign(&tmp, CLS, sizeof(projscs_t)*nxs*ny*nt);
    zfs_send = (projscs_t*)tmp;
    posix_memalign(&tmp, CLS, sizeof(projscs_t)*nxs*ny*nt);
    zbs_send = (projscs_t*)tmp;
    posix_memalign(&tmp, CLS, sizeof(projscs_t)*nxs*ny*nz);
    tfs_send = (projscs_t*)tmp;
    posix_memalign(&tmp, CLS, sizeof(projscs_t)*nxs*ny*nz);
    tbs_send = (projscs_t*)tmp;
    for(int parity=0; parity<2; parity++){
      posix_memalign(&tmp, CLS, sizeof(projscs1_t)*ny*nz*nt);
      xfs_recv_array[parity] = (projscs1_t*)tmp;
      posix_memalign(&tmp, CLS, sizeof(projscs1_t)*ny*nz*nt);
      xbs_recv_array[parity] = (projscs1_t*)tmp;
      posix_memalign(&tmp, CLS, sizeof(projscs_t)*nxs*nz*nt);
      yfs_recv_array[parity] = (projscs_t*)tmp;
      posix_memalign(&tmp, CLS, sizeof(projscs_t)*nxs*nz*nt);
      ybs_recv_array[parity] = (projscs_t*)tmp;
      posix_memalign(&tmp, CLS, sizeof(projscs_t)*nxs*ny*nt);
      zfs_recv_array[parity] = (projscs_t*)tmp;
      posix_memalign(&tmp, CLS, sizeof(projscs_t)*nxs*ny*nt);
      zbs_recv_array[parity] = (projscs_t*)tmp;
      posix_memalign(&tmp, CLS, sizeof(projscs_t)*nxs*ny*nz);
      tfs_recv_array[parity] = (projscs_t*)tmp;
      posix_memalign(&tmp, CLS, sizeof(projscs_t)*nxs*ny*nz);
      tbs_recv_array[parity] = (projscs_t*)tmp;
    }

    xbound_set_parity(0, 4);
    xbound_recv_updateall(4);
  }

  void xbound_finalize() { 
    return;
  }


  //---------------------------------------------------------------------------------------- COMM
  void xbound_start(int req, int prec) {
    return;
  }

  //---------------------------------------------------------------------------------------- COMM
  void xbound(int req, int prec) {
#ifdef _DEBUG
    if (npe[req/2] != 1) {
      printf("cannot happen at %s : ivalid[req/2] without mpi\n", __func__); exit(1);
    }
#endif

    switch (req) {
    case 0 :
      if (prec == 8) {
	memcpy(xfd_recv, xfd_send, sizeof(double)*12*ny*nz*nt);
      } else {
	memcpy(xfs_recv, xfs_send, sizeof(float)*12*ny*nz*nt);
      }
      break;
    case 1 :
      if (prec == 8) {
	memcpy(xbd_recv, xbd_send, sizeof(double)*12*ny*nz*nt);
      } else {
	memcpy(xbs_recv, xbs_send, sizeof(float)*12*ny*nz*nt);
      }
      break;
    case 2 :
      if (prec == 8) {
	memcpy(yfd_recv, yfd_send, sizeof(projscd_t)*nxd*nz*nt);
      } else {
	memcpy(yfs_recv, yfs_send, sizeof(projscs_t)*nxs*nz*nt);
      }
      break;
    case 3 :
      if (prec == 8) {
	memcpy(ybd_recv, ybd_send, sizeof(projscd_t)*nxd*nz*nt);
      } else {
	memcpy(ybs_recv, ybs_send, sizeof(projscs_t)*nxs*nz*nt);
      }
      break;
    case 4 :
      if (prec == 8) {
	memcpy(zfd_recv, zfd_send, sizeof(projscd_t)*nxd*ny*nt);
      } else {
	memcpy(zfs_recv, zfs_send, sizeof(projscs_t)*nxs*ny*nt);
      }
      break;
    case 5 :
      if (prec == 8) {
	memcpy(zbd_recv, zbd_send, sizeof(projscd_t)*nxd*ny*nt);
      } else {
	memcpy(zbs_recv, zbs_send, sizeof(projscs_t)*nxs*ny*nt);
      }
      break;
    case 6 :
      if (prec == 8) {
	memcpy(tfd_recv, tfd_send, sizeof(projscd_t)*nxd*ny*nz);
      } else {
	memcpy(tfs_recv, tfs_send, sizeof(projscs_t)*nxs*ny*nz);
      }
      break;
    case 7 :
      if (prec == 8) {
	memcpy(tbd_recv, tbd_send, sizeof(projscd_t)*nxd*ny*nz);
      } else {
	memcpy(tbs_recv, tbs_send, sizeof(projscs_t)*nxs*ny*nz);
      }
      break;
    }
  }
  //---------------------------------------------------------------------------------------- wait
  void xbound_wait(int req, int prec) {
    return;
  }
  //---------------------------------------------------------------------------------------- wait
  void xbound_reset_comm(int req, int prec) {
    return;
  }
  //---------------------------------------------------------------------------------------- send
  void xbound_send_waitall(int prec) {
    return;
  }
  //---------------------------------------------------------------------------------------- recv
  void xbound_recv_waitall(int prec) {
    return;
  }

  void xbound_recv_okall(int prec) {
    return;
  }

#ifdef __cplusplus
}
#endif
