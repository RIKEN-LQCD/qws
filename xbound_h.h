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
#ifndef _XBOUND_H_H
#define _XBOUND_H_H

  projsch1_t * __restrict__ xfh_send;
  projsch1_t * __restrict__ xfh_recv;
  projsch1_t * __restrict__ xbh_send;
  projsch1_t * __restrict__ xbh_recv;
                            
  projsch_t  * __restrict__ yfh_send;
  projsch_t  * __restrict__ yfh_recv;
  projsch_t  * __restrict__ ybh_send;
  projsch_t  * __restrict__ ybh_recv;
                            
  projsch_t  * __restrict__ zfh_send;
  projsch_t  * __restrict__ zfh_recv;
  projsch_t  * __restrict__ zbh_send;
  projsch_t  * __restrict__ zbh_recv;
                           
  projsch_t  * __restrict__ tfh_send;
  projsch_t  * __restrict__ tfh_recv;
  projsch_t  * __restrict__ tbh_send;
  projsch_t  * __restrict__ tbh_recv;

#ifdef _MPI_
  MPI_Request sh_req[8];
  MPI_Request rh_req[8];
#endif

  //---------------------------------------------------------------------------------------- init
  void xbound_h_init_()
  {

//printf("sizeof(char)       = %d\n",sizeof(char));
//printf("sizeof(short)      = %d\n",sizeof(short));
//printf("sizeof(half)       = %d\n",sizeof(half));
//printf("sizeof(projsch_t)  = %d\n",sizeof(projsch_t));

    ////////////////////////////
    // allocate ghost/halo buffers
    ////////////////////////////
    xfh_send = (projsch1_t *)malloc( sizeof(projsch1_t) * ny *nz*nt);
    xbh_send = (projsch1_t *)malloc( sizeof(projsch1_t) * ny *nz*nt);
    yfh_send = (projsch_t  *)malloc( sizeof(projsch_t)  * nxs*nz*nt);
    ybh_send = (projsch_t  *)malloc( sizeof(projsch_t)  * nxs*nz*nt);
    zfh_send = (projsch_t  *)malloc( sizeof(projsch_t)  * nxs*ny*nt);
    zbh_send = (projsch_t  *)malloc( sizeof(projsch_t)  * nxs*ny*nt);
    tfh_send = (projsch_t  *)malloc( sizeof(projsch_t)  * nxs*ny*nz);
    tbh_send = (projsch_t  *)malloc( sizeof(projsch_t)  * nxs*ny*nz);
    xfh_recv = (projsch1_t *)malloc( sizeof(projsch1_t) * ny *nz*nt);
    xbh_recv = (projsch1_t *)malloc( sizeof(projsch1_t) * ny *nz*nt);
    yfh_recv = (projsch_t  *)malloc( sizeof(projsch_t)  * nxs*nz*nt);
    ybh_recv = (projsch_t  *)malloc( sizeof(projsch_t)  * nxs*nz*nt);
    zfh_recv = (projsch_t  *)malloc( sizeof(projsch_t)  * nxs*ny*nt);
    zbh_recv = (projsch_t  *)malloc( sizeof(projsch_t)  * nxs*ny*nt);
    tfh_recv = (projsch_t  *)malloc( sizeof(projsch_t)  * nxs*ny*nz);
    tbh_recv = (projsch_t  *)malloc( sizeof(projsch_t)  * nxs*ny*nz);

#ifdef _MPI_
    //////////////////////////////////////////////////////////////////////////////
    // Initialize Persistent Send/Recv requests for ghost/halo sites exchange
    //////////////////////////////////////////////////////////////////////////////
    MPI_Send_init(xfh_send, 12*ny *nz*nt, MPI_UNSIGNED_SHORT, pxb, 0, MPI_COMM_WORLD, &sh_req[0]);
    MPI_Send_init(xbh_send, 12*ny *nz*nt, MPI_UNSIGNED_SHORT, pxf, 1, MPI_COMM_WORLD, &sh_req[1]);
    MPI_Send_init(yfh_send, 12*nxh*nz*nt, MPI_UNSIGNED_SHORT, pyb, 2, MPI_COMM_WORLD, &sh_req[2]);
    MPI_Send_init(ybh_send, 12*nxh*nz*nt, MPI_UNSIGNED_SHORT, pyf, 3, MPI_COMM_WORLD, &sh_req[3]);
    MPI_Send_init(zfh_send, 12*nxh*ny*nt, MPI_UNSIGNED_SHORT, pzb, 4, MPI_COMM_WORLD, &sh_req[4]);
    MPI_Send_init(zbh_send, 12*nxh*ny*nt, MPI_UNSIGNED_SHORT, pzf, 5, MPI_COMM_WORLD, &sh_req[5]);
    MPI_Send_init(tfh_send, 12*nxh*ny*nz, MPI_UNSIGNED_SHORT, ptb, 6, MPI_COMM_WORLD, &sh_req[6]);
    MPI_Send_init(tbh_send, 12*nxh*ny*nz, MPI_UNSIGNED_SHORT, ptf, 7, MPI_COMM_WORLD, &sh_req[7]);
    MPI_Recv_init(xfh_recv, 12*ny *nz*nt, MPI_UNSIGNED_SHORT, pxf, 0, MPI_COMM_WORLD, &rh_req[0]);
    MPI_Recv_init(xbh_recv, 12*ny *nz*nt, MPI_UNSIGNED_SHORT, pxb, 1, MPI_COMM_WORLD, &rh_req[1]);
    MPI_Recv_init(yfh_recv, 12*nxh*nz*nt, MPI_UNSIGNED_SHORT, pyf, 2, MPI_COMM_WORLD, &rh_req[2]);
    MPI_Recv_init(ybh_recv, 12*nxh*nz*nt, MPI_UNSIGNED_SHORT, pyb, 3, MPI_COMM_WORLD, &rh_req[3]);
    MPI_Recv_init(zfh_recv, 12*nxh*ny*nt, MPI_UNSIGNED_SHORT, pzf, 4, MPI_COMM_WORLD, &rh_req[4]);
    MPI_Recv_init(zbh_recv, 12*nxh*ny*nt, MPI_UNSIGNED_SHORT, pzb, 5, MPI_COMM_WORLD, &rh_req[5]);
    MPI_Recv_init(tfh_recv, 12*nxh*ny*nz, MPI_UNSIGNED_SHORT, ptf, 6, MPI_COMM_WORLD, &rh_req[6]);
    MPI_Recv_init(tbh_recv, 12*nxh*ny*nz, MPI_UNSIGNED_SHORT, ptb, 7, MPI_COMM_WORLD, &rh_req[7]);
#endif

  }


  //---------------------------------------------------------------------------------------- COMM
  void xbound_h_start(int req)
  {
    if (npe[req/2] != 1) {
#ifdef _MPI_
      //
      // it is better to use MPI_Startall by packing all requests having (npe !=1) into a array.
      //
      MPI_Start(&rh_req[req]);
      MPI_Start(&sh_req[req]);
#endif
    } else {
      switch (req) {
      case 0 :   // X-forward
        memcpy(xfh_recv, xfh_send, sizeof(half)*12*ny*nz*nt);
	break;
      case 1 :   // X-backward
        memcpy(xbh_recv, xbh_send, sizeof(half)*12*ny*nz*nt);
	break;
      case 2 :   // Y-forward
        memcpy(yfh_recv, yfh_send, sizeof(projsch_t)*nxs*nz*nt);
	break;
      case 3 :   // Y-backward
        memcpy(ybh_recv, ybh_send, sizeof(projsch_t)*nxs*nz*nt);
	break;
      case 4 :   // Z-forward
        memcpy(zfh_recv, zfh_send, sizeof(projsch_t)*nxs*ny*nt);
	break;
      case 5 :   // Z-backward
        memcpy(zbh_recv, zbh_send, sizeof(projsch_t)*nxs*ny*nt);
	break;
      case 6 :   // T-forward
        memcpy(tfh_recv, tfh_send, sizeof(projsch_t)*nxs*ny*nz);
	break;
      case 7 :   // T-backward
        memcpy(tbh_recv, tbh_send, sizeof(projsch_t)*nxs*ny*nz);
	break;
      }
    }
  }
  //---------------------------------------------------------------------------------------- wait
  void xbound_h_wait(int req)
  {
#if _MPI_
    MPI_Status status;
    if (npe[req/2] != 1) {
      MPI_Wait(&rh_req[req], &status);
    }
#endif
    return;
  }
  //---------------------------------------------------------------------------------------- send
  void xbound_h_send_waitall() 
  {
#if _MPI_
    MPI_Status status;
    //
    // it is better to use MPI_Waitall by packing all requests with (npe !=1) into a array.
    //
    for (int i=0;i<4;i++) {
      if (npe[i] != 1) {
        MPI_Wait(&sh_req[0+2*i], &status);
        MPI_Wait(&sh_req[1+2*i], &status);
      }
    }
#endif
    return;
  }
  //---------------------------------------------------------------------------------------- recv
  void xbound_h_recv_waitall() {
#if _MPI_
    MPI_Status status;
    //
    // it is better to use MPI_Waitall by packing all requests with (npe !=1) into a array.
    //
    for (int i=0;i<4;i++) {
      if (npe[i] != 1) {
        MPI_Wait(&rh_req[0+2*i], &status);
        MPI_Wait(&rh_req[1+2*i], &status);
      }
    }
#endif
    return;
  }

#endif
