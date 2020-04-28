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
#include "rdma_comlib_2buf.h"
#ifdef _UTOFU_RDMA
//#include "get_tni.h"
#endif

#ifdef _MPI_

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

  int bufs_parity;

  MPI_Request sd_req[8];
  MPI_Request rd_req[8];
  rdma_comlib_2buf buff_rdma[8];

  int rankmap_id;

  projscd1_t *xfd_send0;
  projscd1_t *xfd_recv0;
  projscd1_t *xbd_send0;
  projscd1_t *xbd_recv0;

  void xbound_set_parity(int parity, int prec){
    if(prec == 8){
      return;
    }
    if(bufs_parity == parity){
      return;
    }
    bufs_parity=parity;
    for(int dir=0; dir<8; dir++){
      buff_rdma[dir].flip_parity();
    }
  }

  void xbound_flip_parity(int prec){
    if(prec == 8){
      return;
    }
    xbound_set_parity(1-bufs_parity, prec);
  }

  void xbound_recv_updateall(int prec){
    if(prec==4){
      xfs_recv=(projscs1_t*)buff_rdma[0].rbuff();
      xbs_recv=(projscs1_t*)buff_rdma[1].rbuff();
      yfs_recv=(projscs_t*)buff_rdma[2].rbuff();
      ybs_recv=(projscs_t*)buff_rdma[3].rbuff();
      zfs_recv=(projscs_t*)buff_rdma[4].rbuff();
      zbs_recv=(projscs_t*)buff_rdma[5].rbuff();
      tfs_recv=(projscs_t*)buff_rdma[6].rbuff();
      tbs_recv=(projscs_t*)buff_rdma[7].rbuff();
    }
  }

  void xbound_recv_update(int req, int prec){
    if(prec==4){
      // buff_rdma[0]:  send to -x (xfs_send), recv from -x (xfb_recv)
      // buff_rdma[1]:  send to +x (xfs_send), recv from +x (xfb_recv)
      // ...
      switch(req){
      case 0:
        xfs_recv=(projscs1_t*)buff_rdma[0].rbuff();
        //#pragma omp flush(xfs_recv)
        break;
      case 1:
        xbs_recv=(projscs1_t*)buff_rdma[1].rbuff();
        //#pragma omp flush(xbs_recv)
        break;
      case 2:
        yfs_recv=(projscs_t*)buff_rdma[2].rbuff();
        //#pragma omp flush(yfs_recv)
        break;
      case 3:
        ybs_recv=(projscs_t*)buff_rdma[3].rbuff();
        //#pragma omp flush(ybs_recv)
        break;
      case 4:
        zfs_recv=(projscs_t*)buff_rdma[4].rbuff();
        //#pragma omp flush(zfs_recv)
        break;
      case 5:
        zbs_recv=(projscs_t*)buff_rdma[5].rbuff();
        //#pragma omp flush(zbs_recv)
        break;
      case 6:
        tfs_recv=(projscs_t*)buff_rdma[6].rbuff();
        //#pragma omp flush(tfs_recv)
        break;
      case 7:
        tbs_recv=(projscs_t*)buff_rdma[7].rbuff();
        //#pragma omp flush(zbs_recv)
        break;
      }
    }
  }

 //---------------------------------------------------------------------------------------- get rankmap
  int xbound_get_rankmap(int myrank, int *rank_coords, int *rank_size,
                         int *rankf, int *rankb) {

    static int dim=4;
    int mapid=rdma_comlib_2buf::set_rankmap(dim);
    rankmap_id=mapid;
    if(myrank==0){
      fprintf(stderr, "set_rankmap, done: rankmap_id=%d\n", rankmap_id);
    }
    if(mapid<0){ // no rank map is available
      return mapid;
    }
    int neighbor_ranks[2*dim];
    rdma_comlib_2buf::get_rank_map(rank_coords, neighbor_ranks, rank_size);
    for(int dir=0; dir<dim; dir++){
      rankf[dir]=neighbor_ranks[2*dir];
      rankb[dir]=neighbor_ranks[2*dir+1];
    }
    return mapid;
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
    xfd_send0 = xfd_send;
    xbd_send0 = xbd_send;
    xfd_recv0 = xfd_recv;
    xbd_recv0 = xbd_recv;

    // allocate communication buffers (single prec.)
    // no manual allcoation is needed for RDMA

    // initializing communications: double prec.
    MPI_Send_init(xfd_send0, 12*ny *nz*nt, MPI_DOUBLE_PRECISION, pxb, 0, MPI_COMM_WORLD, &sd_req[0]);
    MPI_Send_init(xbd_send0, 12*ny *nz*nt, MPI_DOUBLE_PRECISION, pxf, 1, MPI_COMM_WORLD, &sd_req[1]);
    MPI_Send_init(yfd_send, 12*nxh*nz*nt, MPI_DOUBLE_PRECISION, pyb, 2, MPI_COMM_WORLD, &sd_req[2]);
    MPI_Send_init(ybd_send, 12*nxh*nz*nt, MPI_DOUBLE_PRECISION, pyf, 3, MPI_COMM_WORLD, &sd_req[3]);
    MPI_Send_init(zfd_send, 12*nxh*ny*nt, MPI_DOUBLE_PRECISION, pzb, 4, MPI_COMM_WORLD, &sd_req[4]);
    MPI_Send_init(zbd_send, 12*nxh*ny*nt, MPI_DOUBLE_PRECISION, pzf, 5, MPI_COMM_WORLD, &sd_req[5]);
    MPI_Send_init(tfd_send, 12*nxh*ny*nz, MPI_DOUBLE_PRECISION, ptb, 6, MPI_COMM_WORLD, &sd_req[6]);
    MPI_Send_init(tbd_send, 12*nxh*ny*nz, MPI_DOUBLE_PRECISION, ptf, 7, MPI_COMM_WORLD, &sd_req[7]);
    MPI_Recv_init(xfd_recv0, 12*ny *nz*nt, MPI_DOUBLE_PRECISION, pxf, 0, MPI_COMM_WORLD, &rd_req[0]);
    MPI_Recv_init(xbd_recv0, 12*ny *nz*nt, MPI_DOUBLE_PRECISION, pxb, 1, MPI_COMM_WORLD, &rd_req[1]);
    MPI_Recv_init(yfd_recv, 12*nxh*nz*nt, MPI_DOUBLE_PRECISION, pyf, 2, MPI_COMM_WORLD, &rd_req[2]);
    MPI_Recv_init(ybd_recv, 12*nxh*nz*nt, MPI_DOUBLE_PRECISION, pyb, 3, MPI_COMM_WORLD, &rd_req[3]);
    MPI_Recv_init(zfd_recv, 12*nxh*ny*nt, MPI_DOUBLE_PRECISION, pzf, 4, MPI_COMM_WORLD, &rd_req[4]);
    MPI_Recv_init(zbd_recv, 12*nxh*ny*nt, MPI_DOUBLE_PRECISION, pzb, 5, MPI_COMM_WORLD, &rd_req[5]);
    MPI_Recv_init(tfd_recv, 12*nxh*ny*nz, MPI_DOUBLE_PRECISION, ptf, 6, MPI_COMM_WORLD, &rd_req[6]);
    MPI_Recv_init(tbd_recv, 12*nxh*ny*nz, MPI_DOUBLE_PRECISION, ptb, 7, MPI_COMM_WORLD, &rd_req[7]);

    // initializing communications: single prec.
    rdma_comlib_2buf::comlib_init();
    int tni_list[8];
    if(rank==0){
      fprintf(stderr, "calling get_tni_ids: rankmap_id=%d\n", rankmap_id);
    }
    rdma_comlib_2buf::get_tni_ids(tni_list, rankmap_id);

    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0){
      for(int i=0; i<8; i++){
        // fprintf(stderr, "rank=%d: tni_list[%d] = %d\n",rank, i,tni_list[i]);
        printf("rank=%d: tni_list[%d] = %d\n",rank, i,tni_list[i]);
      }
      printf("rank=%d: pxf=%d, pxb=%d, pyf=%d, pyb=%d, pzf=%d, pzb=%d, ptf=%d, ptb=%d\n", rank, pxf, pxb, pyf, pyb, pzf, pzb, ptf, ptb);
      fflush(0);
    }

    size_t size=sizeof(projscs1_t)*ny*nz*nt;
    buff_rdma[0].init(tni_list[0], pxb, pxf, size);  // put to pxb, recv from pxf
    buff_rdma[1].init(tni_list[1], pxf, pxb, size);  // put to pxf, recv from pxb
    rdma_comlib_2buf::swap_vcq_for_sending(buff_rdma[0], buff_rdma[1]);
    // after swapping:
    //   to/from pxf: tni_list[0]
    //   to/from pxb: tni_list[1]
    if(rank==0){
      printf("    swap_vcq_for_sending(buff_rdma[0], buff_rdma[1]), done\n");
      fflush(0);
    }

    size=sizeof(projscs_t)*nxs*nz*nt;
    buff_rdma[2].init(tni_list[2], pyb, pyf, size);  // put to pyb, recv from pyf
    buff_rdma[3].init(tni_list[3], pyf, pyb, size);  // put to pyf, recv from pyb
    rdma_comlib_2buf::swap_vcq_for_sending(buff_rdma[2], buff_rdma[3]);
    // after swapping:
    //   to/from pyf: tni_list[2]
    //   to/from pyb: tni_list[3]
    if(rank==0){
      printf("    swap_vcq_for_sending(buff_rdma[2], buff_rdma[3]), done\n");
      fflush(0);
    }

    size=sizeof(projscs_t)*nxs*ny*nt;
    buff_rdma[4].init(tni_list[4], pzb, pzf, size);  // put to pzb, recv from pzf
    buff_rdma[5].init(tni_list[5], pzf, pzb, size);  // put to pzf, recv from pzb
    rdma_comlib_2buf::swap_vcq_for_sending(buff_rdma[4], buff_rdma[5]);
    // after swapping:
    //   to/from pzf: tni_list[4]
    //   to/from pzb: tni_list[5]
    if(rank==0){
      printf("    swap_vcq_for_sending(buff_rdma[4], buff_rdma[5]), done\n");
      fflush(0);
    }

    size=sizeof(projscs_t)*nxs*ny*nz;
    buff_rdma[6].init(tni_list[6], ptb, ptf, size);  // put to ptb, recv from ptf

    buff_rdma[7].init(tni_list[7], ptf, ptb, size);  // put to ptf, recv from ptb

    rdma_comlib_2buf::swap_vcq_for_sending(buff_rdma[6], buff_rdma[7]);
    // after swapping:
    //   to/from ptf: tni_list[6]
    //   to/from ptb: tni_list[7]
    if(rank==0){
      printf("    swap_vcq_for_sending(buff_rdma[6], buff_rdma[7]), done\n");
      fflush(0);
    }

    xfs_send=(projscs1_t*)buff_rdma[0].sbuff();
    xbs_send=(projscs1_t*)buff_rdma[1].sbuff();
    yfs_send=(projscs_t*)buff_rdma[2].sbuff();
    ybs_send=(projscs_t*)buff_rdma[3].sbuff();
    zfs_send=(projscs_t*)buff_rdma[4].sbuff();
    zbs_send=(projscs_t*)buff_rdma[5].sbuff();
    tfs_send=(projscs_t*)buff_rdma[6].sbuff();
    tbs_send=(projscs_t*)buff_rdma[7].sbuff();
    xbound_set_parity(0, 4);
    xbound_recv_updateall(4);

    if(rank==0){
      printf("    set_parity, recv_updateall, done\n");
      fflush(0);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0){
      printf("qws_xbound is ready (rdma version).\n");
      fflush(0);
    }
  }


  void xbound_finalize() {

    // make sure all the communications are finished
    MPI_Barrier(MPI_COMM_WORLD);

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
    for(int i=0; i<8; i++){
      buff_rdma[i].finalize();
    }
    rdma_comlib_2buf::comlib_finalize();

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
        buff_rdma[req].isendrecv();
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
        buff_rdma[req].irecv_check();
      }
    }
  }
  //---------------------------------------------------------------------------------------- wait
  void xbound_reset_comm(int req, int prec) {
    if (npe[req/2] != 1) {
      if (prec == 4) {
        buff_rdma[req].reset_comm();
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
          buff_rdma[0+2*i].isend_check();
          buff_rdma[1+2*i].isend_check();
        }
      }
    }
  }
  //---------------------------------------------------------------------------------------- recv
  void xbound_recv_waitall(int prec) {
    int i;
    MPI_Status status;
    if (prec == 8) {
      for (i=0;i<4;i++) {
        if (npe[i] != 1) {
          MPI_Wait(&rd_req[0+2*i], &status);
          MPI_Wait(&rd_req[1+2*i], &status);
        }
      }
    } else {
      for (i=0;i<4;i++) {
        if (npe[i] != 1) {
          buff_rdma[0+2*i].irecv_check();
          buff_rdma[1+2*i].irecv_check();
        }
      }
    }
  }

  void xbound_recv_okall(int prec) {
    if(prec==4){
      for(int dir=0; dir<8; dir++){
	buff_rdma[dir].irecv_ok();
      }
      // reset the send buffers: might be pointing to the rbuff
      xfs_send=(projscs1_t*)buff_rdma[0].sbuff();
      xbs_send=(projscs1_t*)buff_rdma[1].sbuff();
    } else {
      xfd_send = xfd_send0;
      xfd_recv = xfd_recv0;
      xbd_send = xbd_send0;
      xbd_recv = xbd_recv0;
    }
  }

#ifdef __cplusplus
}
#endif

#else // _MPI_
#error _MPI_ is not defeind

#endif // _MPI_
