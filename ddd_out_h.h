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
#ifndef _DDD_OUT_H_H
#define _DDD_OUT_H_H

#include "ddd_out_h_0_inline.h"

  //---------------------------------------------------------------------------------------- preprocess mult D for boundary

  void ddd_out_pre_h_(const sch_t * __restrict__ in, const int *idomain) 
  //
  //  Send for Wilson/Clover operator, Deo or Doe, in a even/odd domain block
  //
  //   oute = oute + factor Meo ine 
  //
  //  or
  //
  //   outo = outo + factor Moe ino 
  //
  // Note: Deo = -kappa Meo, Doe = -kappa Moe
  //
  //       in : quark field in an even/odd domain (odd for idomain = 0, even for idomain = 1). 
  //            input domain is opposite to idomain
  //  idomain : domain index,  0 for even output,  1 for odd output
  //
  {

    //
    // pack data for send
    //

    _S_MULT_WD_DEO_OUT_SEND_HPC_TIC_;

    int xdmn, ydmn, zdmn, tdmn;

    xdmn = 1 - *idomain;
    if (npe[1] == 1) { ydmn = *idomain; } else { ydmn = 1 - *idomain; }
    if (npe[2] == 1) { zdmn = *idomain; } else { zdmn = 1 - *idomain; }
    if (npe[3] == 1) { tdmn = *idomain; } else { tdmn = 1 - *idomain; }

    pgluh_t __restrict__ gx = &gluh[vols*0 + NDIM*vols*xdmn];
    pgluh_t __restrict__ gy = &gluh[vols*1 + NDIM*vols*ydmn];
    pgluh_t __restrict__ gz = &gluh[vols*2 + NDIM*vols*zdmn];
    pgluh_t __restrict__ gt = &gluh[vols*3 + NDIM*vols*tdmn];


#pragma omp parallel
    {

#pragma omp for collapse(3) nowait
    for (int t = 0; t < nt; ++t) {
    for (int z = 0; z < nz; ++z) {
    for (int y = 0; y < ny; ++y) {

      int it = nxs*ny*nz*t;
      int iz = nxs*ny*z;
      int iy = nxs*y;

      // X-forward
      {
        //int x = nxs-1;

        int ixf = 0 + iy + iz + it;
        int i0  = y + ny*z + ny*nz*t;

        __mult_x_forw_pre_3_(xfh_send[i0],in[vols*xdmn+ixf]);

      }

      // X-backward
      { 
        //int x = 0;

        int ixb = nxs-1 + iy + iz + it;
        int i1  = y + ny*z + ny*nz*t;
        projsch1_t a __attribute__((aligned(_ALIGN_SIZE)));

        __mult_x_back_pre_3_(a,in[vols*xdmn+ixb]);
        __mult_udag_y_3_(xbh_send[i1],a,(*(gx + ixb)));

      }
    }
    }
    }

#pragma omp for collapse(3) nowait
    for (int t = 0; t < nt; ++t) {
    for (int z = 0; z < nz; ++z) {
      for (int x = 0; x < nxs; ++x) {

        int it = nxs*ny*nz*t;
        int iz = nxs*ny*z;

        //
        // Y-forward send
        //
        {
          //int y = ny-1;

          int iyf = x + iz + it;
          int i2  = x + nxs*z + nxs*nz*t;

          __mult_y_forw_pre_(yfh_send[i2],in[vols*ydmn+iyf]);

        }

        //
        // Y-backward send
        //
        {
          //int y = 0;

          int iyb = x + iz + it + (ny-1)*nxs;
          int i3  = x + nxs*z + nxs*nz*t;

          __attribute__((aligned(_ALIGN_SIZE))) projsch_t a;

          __mult_y_back_pre_(a,in[vols*ydmn+iyb]);
          __mult_udag_y_(ybh_send[i3],a,(*(gy + iyb)));

        }
      }
    }
    }

#pragma omp for collapse(3) nowait
    for (int t = 0; t < nt; ++t) {
    for (int y = 0; y < ny; ++y) {
      for (int x = 0; x < nxs; ++x) {

        int it = nxs*ny*nz*t;
        int iy = nxs*y;

        //
        // Z-forward send
        //
        { 
          //int z = nz-1;

          int izf = x + iy + it;
          int i4  = x + nxs*y + nxs*ny*t;

          __mult_z_forw_pre_(zfh_send[i4],in[vols*zdmn+izf]);


        }

        //
        // Z-backward send
        //
        {
          //int z = 0;

          int izb = x + iy + it + (nz-1)*nxs*ny;
          int i5  = x + nxs*y + nxs*ny*t;

          __attribute__((aligned(_ALIGN_SIZE))) projsch_t a;

          __mult_z_back_pre_(a,in[vols*zdmn+izb]);
          __mult_udag_y_(zbh_send[i5],a,(*(gz + izb)));


        }
      }
    }
    }

#pragma omp for collapse(3) nowait
    for (int z = 0; z < nz; ++z) {
    for (int y = 0; y < ny; ++y) {
      for (int x = 0; x < nxs; ++x) {

       int iz = nxs*ny*z;
       int iy = nxs*y;

        //
        // T-forward send
        //
        {
          //int t = nt-1;

          int itf = x + iy + iz;
          int i6  = x + nxs*y + nxs*ny*z;

          __mult_t_forw_pre_(tfh_send[i6],in[vols*tdmn+itf]);

        }

        //
        // T-backward send
        //
        {
          //int t = 0;

          int itb = x + iy + iz + (nt-1)*nxs*ny*nz;
          int i7  = x + nxs*y + nxs*ny*z;

          __attribute__((aligned(_ALIGN_SIZE))) projsch_t a;

          __mult_t_back_pre_(a,in[vols*tdmn+itb]);
          __mult_udag_y_(tbh_send[i7],a,(*(gt + itb)));

        }

      }
    }
    }

#pragma omp barrier
    } // end of omp parallel

    _S_MULT_WD_DEO_OUT_SEND_HPC_TOC_;
    _COMM_TIC_;
#pragma omp single nowait
    {
      if (*idomain == 0) {
        memcpy(xfh_recv, xfh_send, sizeof(half)*12*ny*nz*nt);
      } else {
        xbound_h_start(0);
      }
      if (*idomain == 1) {
        memcpy(xbh_recv, xbh_send, sizeof(half)*12*ny*nz*nt);
      } else {
        xbound_h_start(1);
      }
      xbound_h_start(2);
      xbound_h_start(3);
      xbound_h_start(4);
      xbound_h_start(5);
      xbound_h_start(6);
      xbound_h_start(7);
    }
    _COMM_TOC_;
  }

  //---------------------------------------------------------------------------------------- postprocess mult D for boundary

  void ddd_out_pos_h_(sch_t * __restrict__ out, const sch_t * __restrict__ in, const int *idomain, float factor) 
  //
  // Multiply Wilson/Clover operator from odd/even to even/odd domain
  //
  //   oute = oute + factor Meo ino
  // or
  //   outo = outo + factor Moe ine
  //
  // Note: Deo = -kappa Meo, Doe = -kappa Moe
  //
  //      out : quark field in an even/odd domain (output)
  //       in : quark field in an odd/even domain (input)
  //  idomain : domain index, 0 for even output, 1 for odd output
  //   factor : a factor for hopping parameter kappa or -kappa
  //
  {

    pgluh_t __restrict__ gx = &gluh[vols*0 + NDIM*vols*(*idomain)];
    pgluh_t __restrict__ gy = &gluh[vols*1 + NDIM*vols*(*idomain)];
    pgluh_t __restrict__ gz = &gluh[vols*2 + NDIM*vols*(*idomain)];
    pgluh_t __restrict__ gt = &gluh[vols*3 + NDIM*vols*(*idomain)];
    half hfactor = half(factor);

    _COMM_TIC_;
#pragma omp single
    {
      xbound_h_wait(0);
      xbound_h_wait(1);
      xbound_h_wait(2);
      xbound_h_wait(3);
      xbound_h_wait(4);
      xbound_h_wait(5);
      xbound_h_wait(6);
      xbound_h_wait(7);
    }
    _COMM_TOC_;

    _S_MULT_WD_DEO_OUT_RECV_HPC_CALC_TIC_;

#pragma omp parallel for collapse(4)
    for (int t = 0; t < nt; t++) {
    for (int z = 0; z < nz; z++) {
    for (int y = 0; y < ny; y++) {
      for (int x = 0; x < nxs; x++) {

        if (x == nxs-1 || y == ny-1 || z == nz-1 || t == nt-1 || x==0 || y == 0 || z == 0 || t == 0 ) {

          int it = nxs*ny*nz*t;
          int iz = nxs*ny*z;
          int iy = nxs*y;
          int i0 = x + iy + iz + it;

          sch_t tmp __attribute__((aligned(_ALIGN_SIZE))) ;

          for(int c = 0; c < 3; c++) {
          for(int s = 0; s < 4; s++) {
            for(int ri = 0; ri < 2; ri++) {
              for(int j = 0; j < VLENS; j++) {
                tmp.c[c][s][ri][j] = half(0.0f);
              }
            }
          }
          }

          //
          // X-forward
          //
          if ( x == nxs-1 ) { _mult_forw_x_recv_; }

          //
          // X-backward
          //
          if ( x == 0 ) { _mult_back_x_recv_; }

          //
          // Y-forward 
          //
          if ( y == ny-1 ) { _mult_forw_y_recv_; }

          //
          // Y-backward
          //
          if ( y == 0 ) { _mult_back_y_recv_; }

          //
          // Z-forward
          //
          if ( z == nz-1 ) { _mult_forw_z_recv_; }

          //
          // Z-backward
          //
          if ( z == 0 ) { _mult_back_z_recv_; }

          //
          // T-forward
          //
          if ( t == nt-1 ) { _mult_forw_t_recv_; }

          //
          // T-backward
          //
          if ( t == 0 ) { _mult_back_t_recv_; }

          __mult_clvh( tmp.cv, clvh[i0 + vols*(*idomain)].cv);

          sch_t *outi = out + i0;
//#pragma loop norecurrence
          for (int c = 0; c < 3; c++) {
          for (int s = 0; s < 4; s++) {
            for (int j = 0; j < VLENS; j++) {
              outi->c[c][s][0][j] += tmp.c[c][s][0][j] * hfactor;
              outi->c[c][s][1][j] += tmp.c[c][s][1][j] * hfactor;
            }
          }
          }

        } // if edge sites

      }
    }
    }
    }

    _S_MULT_WD_DEO_OUT_RECV_HPC_CALC_TOC_;
    _COMM_TIC_;
#pragma omp single
    {
      xbound_h_send_waitall();
    }
    _COMM_TOC_;
  }


#endif
