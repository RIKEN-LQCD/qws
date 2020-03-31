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

#ifdef _MPI_


#ifdef __cplusplus
extern "C"{
#endif

#ifdef _USE_RANKMAP
#include "rankmap_lib.h"
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


  MPI_Request sd_req[8];
  MPI_Request rd_req[8];
  MPI_Request *ss_req;
  MPI_Request *rs_req;
  MPI_Request *rs_req_opposite;
  MPI_Request ss_req_array[2][8];
  MPI_Request rs_req_array[2][8];

  int recv_started_array[2][8];
  int *recv_started;
  int *recv_started_opposite;

#ifdef _USE_RANKMAP
  int rankmap_id;
#endif

 void xbound_set_parity(int parity, int prec){
    if(prec == 8){
      return;
    }
    bufs_parity=parity;
    ss_req = ss_req_array[parity];
    rs_req = rs_req_array[parity];
    rs_req_opposite = rs_req_array[1-parity];
    recv_started = &recv_started_array[parity][0];
    recv_started_opposite = &recv_started_array[1-parity][0];
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

 //---------------------------------------------------------------------------------------- get rankmap
  int xbound_get_rankmap(int myrank, int *rank_coords, int *rank_size,
			 int *rankf, int *rankb) {
#ifdef _USE_RANKMAP
    static int dim=4;
    int mapid=rankmap_lib_set_rankmap4d();
    rankmap_id=mapid;
    if(mapid<0){ // no rank map is available
      return mapid;
    }
    int neighbor_ranks[2*dim];
    rankmap_lib_get_rankmap(rank_coords, neighbor_ranks, rank_size);
    for(int dir=0; dir<dim; dir++){
      rankf[dir]=neighbor_ranks[2*dir];
      rankb[dir]=neighbor_ranks[2*dir+1];
    }
    return mapid;
#else
    return -1;  // empty
#endif

  }

  //---------------------------------------------------------------------------------------- init communication
  void xbound_init(int rank,
		   int pxf, int pyf, int pzf, int ptf,
		   int pxb, int pyb, int pzb, int ptb)
  {

    int nxh=nx/2;
    int nxd=nx/2/VLEND;
    int nxs=nx/2/VLENS;

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

    // initializing communications: double prec.
    MPI_Send_init(xfd_send, 12*ny *nz*nt, MPI_DOUBLE_PRECISION, pxb, 0, MPI_COMM_WORLD, &sd_req[0]);
    MPI_Send_init(xbd_send, 12*ny *nz*nt, MPI_DOUBLE_PRECISION, pxf, 1, MPI_COMM_WORLD, &sd_req[1]);
    MPI_Send_init(yfd_send, 12*nxh*nz*nt, MPI_DOUBLE_PRECISION, pyb, 2, MPI_COMM_WORLD, &sd_req[2]);
    MPI_Send_init(ybd_send, 12*nxh*nz*nt, MPI_DOUBLE_PRECISION, pyf, 3, MPI_COMM_WORLD, &sd_req[3]);
    MPI_Send_init(zfd_send, 12*nxh*ny*nt, MPI_DOUBLE_PRECISION, pzb, 4, MPI_COMM_WORLD, &sd_req[4]);
    MPI_Send_init(zbd_send, 12*nxh*ny*nt, MPI_DOUBLE_PRECISION, pzf, 5, MPI_COMM_WORLD, &sd_req[5]);
    MPI_Send_init(tfd_send, 12*nxh*ny*nz, MPI_DOUBLE_PRECISION, ptb, 6, MPI_COMM_WORLD, &sd_req[6]);
    MPI_Send_init(tbd_send, 12*nxh*ny*nz, MPI_DOUBLE_PRECISION, ptf, 7, MPI_COMM_WORLD, &sd_req[7]);
    MPI_Recv_init(xfd_recv, 12*ny *nz*nt, MPI_DOUBLE_PRECISION, pxf, 0, MPI_COMM_WORLD, &rd_req[0]);
    MPI_Recv_init(xbd_recv, 12*ny *nz*nt, MPI_DOUBLE_PRECISION, pxb, 1, MPI_COMM_WORLD, &rd_req[1]);
    MPI_Recv_init(yfd_recv, 12*nxh*nz*nt, MPI_DOUBLE_PRECISION, pyf, 2, MPI_COMM_WORLD, &rd_req[2]);
    MPI_Recv_init(ybd_recv, 12*nxh*nz*nt, MPI_DOUBLE_PRECISION, pyb, 3, MPI_COMM_WORLD, &rd_req[3]);
    MPI_Recv_init(zfd_recv, 12*nxh*ny*nt, MPI_DOUBLE_PRECISION, pzf, 4, MPI_COMM_WORLD, &rd_req[4]);
    MPI_Recv_init(zbd_recv, 12*nxh*ny*nt, MPI_DOUBLE_PRECISION, pzb, 5, MPI_COMM_WORLD, &rd_req[5]);
    MPI_Recv_init(tfd_recv, 12*nxh*ny*nz, MPI_DOUBLE_PRECISION, ptf, 6, MPI_COMM_WORLD, &rd_req[6]);
    MPI_Recv_init(tbd_recv, 12*nxh*ny*nz, MPI_DOUBLE_PRECISION, ptb, 7, MPI_COMM_WORLD, &rd_req[7]);

    // initializing communications: single prec.
    for(int parity=0; parity<2; parity++){
      // send buffers are common for both parity
      int tag0=parity*8+8;
      MPI_Send_init(xfs_send, 12*ny *nz*nt, MPI_REAL, pxb, tag0+0, MPI_COMM_WORLD, &ss_req_array[parity][0]);
      MPI_Send_init(xbs_send, 12*ny *nz*nt, MPI_REAL, pxf, tag0+1, MPI_COMM_WORLD, &ss_req_array[parity][1]);
      MPI_Send_init(yfs_send, 12*nxh*nz*nt, MPI_REAL, pyb, tag0+2, MPI_COMM_WORLD, &ss_req_array[parity][2]);
      MPI_Send_init(ybs_send, 12*nxh*nz*nt, MPI_REAL, pyf, tag0+3, MPI_COMM_WORLD, &ss_req_array[parity][3]);
      MPI_Send_init(zfs_send, 12*nxh*ny*nt, MPI_REAL, pzb, tag0+4, MPI_COMM_WORLD, &ss_req_array[parity][4]);
      MPI_Send_init(zbs_send, 12*nxh*ny*nt, MPI_REAL, pzf, tag0+5, MPI_COMM_WORLD, &ss_req_array[parity][5]);
      MPI_Send_init(tfs_send, 12*nxh*ny*nz, MPI_REAL, ptb, tag0+6, MPI_COMM_WORLD, &ss_req_array[parity][6]);
      MPI_Send_init(tbs_send, 12*nxh*ny*nz, MPI_REAL, ptf, tag0+7, MPI_COMM_WORLD, &ss_req_array[parity][7]);
      MPI_Recv_init(xfs_recv_array[parity], 12*ny *nz*nt, MPI_REAL, pxf, tag0+0, MPI_COMM_WORLD, &rs_req_array[parity][0]);
      MPI_Recv_init(xbs_recv_array[parity], 12*ny *nz*nt, MPI_REAL, pxb, tag0+1, MPI_COMM_WORLD, &rs_req_array[parity][1]);
      MPI_Recv_init(yfs_recv_array[parity], 12*nxh*nz*nt, MPI_REAL, pyf, tag0+2, MPI_COMM_WORLD, &rs_req_array[parity][2]);
      MPI_Recv_init(ybs_recv_array[parity], 12*nxh*nz*nt, MPI_REAL, pyb, tag0+3, MPI_COMM_WORLD, &rs_req_array[parity][3]);
      MPI_Recv_init(zfs_recv_array[parity], 12*nxh*ny*nt, MPI_REAL, pzf, tag0+4, MPI_COMM_WORLD, &rs_req_array[parity][4]);
      MPI_Recv_init(zbs_recv_array[parity], 12*nxh*ny*nt, MPI_REAL, pzb, tag0+5, MPI_COMM_WORLD, &rs_req_array[parity][5]);
      MPI_Recv_init(tfs_recv_array[parity], 12*nxh*ny*nz, MPI_REAL, ptf, tag0+6, MPI_COMM_WORLD, &rs_req_array[parity][6]);
      MPI_Recv_init(tbs_recv_array[parity], 12*nxh*ny*nz, MPI_REAL, ptb, tag0+7, MPI_COMM_WORLD, &rs_req_array[parity][7]);
      for(int i=0; i<8; i++){
	recv_started_array[parity][i]=0;
      }
    }

    xbound_set_parity(0, 4);
    xbound_recv_updateall(4);
    for(int req=0; req<8; req++){ // start recieving
      if (npe[req/2] != 1) {
	MPI_Start(&rs_req[req]);
	recv_started[req]=1;
      }
    }
  }

  void xbound_finalize() { 

    // make sure all the communications are finished
    MPI_Barrier(MPI_COMM_WORLD);

    // stop recieving
    MPI_Status status;
    for(int parity=0; parity<2; parity++){
      for(int req=0; req<8; req++){
	if (npe[req/2] != 1) {
	  if(recv_started_array[parity][req]) {
	    MPI_Cancel(&rs_req_array[parity][req]);
	    //	    MPI_Wait(&rs_req_array[parity][req], &status);
	    MPI_Wait(&rs_req_array[parity][req], MPI_STATUSES_IGNORE);
	  }
	}
      }
    }

    // send/recv buffers: double precision
    if(xfd_send) free(xfd_send);
    if(xfd_recv) free(xfd_recv);
    if(xbd_send) free(xbd_send);
    if(xbd_recv) free(xbd_recv);

    if(yfd_send) free(yfd_send);
    if(yfd_recv) free(yfd_recv);
    if(ybd_send) free(ybd_send);
    if(ybd_recv) free(ybd_recv);

    if(zfd_send) free(zfd_send);
    if(zfd_recv) free(zfd_recv);
    if(zbd_send) free(zbd_send);
    if(zbd_recv) free(zbd_recv);

    if(tfd_send) free(tfd_send);
    if(tfd_recv) free(tfd_recv);
    if(tbd_send) free(tbd_send);
    if(tbd_recv) free(tbd_recv);

  // for double buffering (single precision)
    if(xfs_send) free(xfs_send);
    if(xbs_send) free(xbs_send);
    if(yfs_send) free(yfs_send);
    if(ybs_send) free(ybs_send);
    if(zfs_send) free(zfs_send);
    if(zbs_send) free(zbs_send);
    if(tfs_send) free(tfs_send);
    if(tbs_send) free(tbs_send);

    for(int parity=0; parity<2; parity++){
      if(xfs_recv_array[parity]) free(xfs_recv_array[parity]);
      if(xbs_recv_array[parity]) free(xbs_recv_array[parity]);
      if(yfs_recv_array[parity]) free(yfs_recv_array[parity]);
      if(ybs_recv_array[parity]) free(ybs_recv_array[parity]);
      if(zfs_recv_array[parity]) free(zfs_recv_array[parity]);
      if(zbs_recv_array[parity]) free(zbs_recv_array[parity]);
      if(tfs_recv_array[parity]) free(tfs_recv_array[parity]);
      if(tbs_recv_array[parity]) free(tbs_recv_array[parity]);
    }

    // wait for all cleaning is done
    MPI_Barrier(MPI_COMM_WORLD);
  }


  //---------------------------------------------------------------------------------------- COMM
  void xbound_start(int req, int prec) {
    if (npe[req/2] != 1) {
      if (prec == 8) {
	MPI_Start(&rd_req[req]);
      }
    }
  }

  //---------------------------------------------------------------------------------------- COMM
  void xbound(int req, int prec) {
    if (npe[req/2] != 1) {
      if (prec == 8) {
	MPI_Start(&sd_req[req]);
      } else {
	MPI_Start(&ss_req[req]);
      }
    } else {
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
  }
  //---------------------------------------------------------------------------------------- wait
  void xbound_wait(int req, int prec) {
    MPI_Status status;
    if (npe[req/2] != 1) {
      if (prec == 8) {
	MPI_Wait(&rd_req[req], &status);
      } else {
	MPI_Wait(&rs_req[req], &status);
	recv_started[req]=0;
	if(!recv_started_opposite[req]){
	  MPI_Start(&rs_req_opposite[req]);
	  recv_started_opposite[req]=1;
	}
      }
    }
  }
  //---------------------------------------------------------------------------------------- reset
  void xbound_reset_comm(int req, int prec) {
    if (prec == 4) {
      if (npe[req/2] != 1){
	if(!recv_started_opposite[req]){
	  MPI_Start(&rs_req_opposite[req]);
	  recv_started_opposite[req]=1;
	}
      }
    }
  }
  //---------------------------------------------------------------------------------------- send
  void xbound_send_waitall(int prec) {
    int i;
    MPI_Status status;
    if (prec == 8) {
      for (i=0;i<4;i++) {
	if (npe[i] != 1) {
	  MPI_Wait(&sd_req[0+2*i], &status);
	  MPI_Wait(&sd_req[1+2*i], &status);
	}
      }
    } else {
      for (i=0;i<4;i++) {
	if (npe[i] != 1) {
	  MPI_Wait(&ss_req[0+2*i], &status);
	  MPI_Wait(&ss_req[1+2*i], &status);
	}
      }
    }
  }

  void xbound_recv_okall(int prec) {
    return;
  }

#ifdef __cplusplus
}
#endif

#else // _MPI_
#error _MPI_ is not defeind

#endif // _MPI_
