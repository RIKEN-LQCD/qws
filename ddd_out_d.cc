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
#include "wilson_d.h"
#include "clover_d.h"

#ifdef __cplusplus
extern "C"{
#endif


  extern __attribute__((aligned(64))) pglud_t glud;
  extern __attribute__((aligned(64))) pclvd_t clvd;
  extern   projscd1_t *xfd_send, *xfd_recv;
  extern   projscd1_t *xbd_send, *xbd_recv;
  extern   projscd_t *yfd_send, *yfd_recv;
  extern   projscd_t *ybd_send, *ybd_recv;
  extern   projscd_t *zfd_send, *zfd_recv;
  extern   projscd_t *zbd_send, *zbd_recv;
  extern   projscd_t *tfd_send, *tfd_recv;
  extern   projscd_t *tbd_send, *tbd_recv;

  extern double mkappa;
  extern int nxd, ny, nz, nt;
  extern int vold;
  extern int pt;
  extern int npe[4];
  extern double fbc[4][2];

  void xbound_start(int req, int prec);
  void xbound(int req, int prec);
  void xbound_wait(int req, int prec);
  void xbound_send_waitall(int prec);
  void xbound_recv_okall(int prec);

  //---------------------------------------------------------------------------------------- preprocess mult D for boundary
  void ddd_out_pre_d_(scd_t* in, int* idomain){
    int x, y, z, t, i0, i1, i2, i3, i4, i5, i6, i7;
    pglud_t gx = &glud[vold*0 + NDIM*vold*(1-(*idomain))];
    pglud_t gy = &glud[vold*1 + NDIM*vold*(*idomain)];
    pglud_t gz = &glud[vold*2 + NDIM*vold*(*idomain)];
    pglud_t gt = &glud[vold*3 + NDIM*vold*(*idomain)];

#ifdef _NO_OMP_SINGLE
#pragma omp master
#else
#pragma omp single nowait
#endif // _NO_OMP_SINGLE
    {
      if (*idomain == 0) {
	xbound_start(1,8);
      } else {
	xbound_start(0,8);
      }
      xbound_start(2,8);
      xbound_start(3,8);
      xbound_start(4,8);
      xbound_start(5,8);
      xbound_start(6,8);
      xbound_start(7,8);
    }

#pragma omp parallel for private(t,z,y,x,i0,i1,i2,i3,i4,i5,i6,i7) collapse(3)
    for(t=0; t<nt; t++){
      for(z=0; z<nz; z++){
	for(y=0; y<ny; y++){
	  int it = nxd*ny*nz*t;
	  int iz = nxd*ny*z;
	  int iy = nxd*y;
	  for(x=0; x<nxd; x++){
	    int ixf = 0     + iy + iz + it;
	    int ixb = nxd-1 + iy + iz + it;
	    int iyf = x + iz + it;
	    int iyb = x + iz + it + (ny-1)*nxd;
	    int izf = x + iy + it;
	    int izb = x + iy + it + (nz-1)*nxd*ny;
	    int itf = x + iy + iz;
	    int itb = x + iy + iz + (nt-1)*nxd*ny*nz;
	    if(x == nxd-1) {i0 = y + ny*z + ny*nz*t;__l_xf(xfd_send[i0].c, &in[vold*(1-(*idomain))], ixf);}
	    if(y == ny-1 ) {i2 = x +nxd*z +nxd*nz*t;__l_yf(yfd_send[i2].c, &in[vold*((*idomain))], iyf);}
	    if(z == nz-1 ) {i4 = x +nxd*y +nxd*ny*t;__l_zf(zfd_send[i4].c, &in[vold*((*idomain))], izf);}
	    if(t == nt-1 ) {i6 = x +nxd*y +nxd*ny*z;__l_tf(tfd_send[i6].c, &in[vold*((*idomain))], itf);}
	    if(x == 0    ) {i1 = y + ny*z + ny*nz*t;__pre_xbwd(xbd_send[i1].c, gx, &in[vold*(1-(*idomain))], ixb);}
	    if(y == 0    ) {i3 = x +nxd*z +nxd*nz*t;__pre_ybwd(ybd_send[i3].c, gy, &in[vold*((*idomain))], iyb);}
	    if(z == 0    ) {i5 = x +nxd*y +nxd*ny*t;__pre_zbwd(zbd_send[i5].c, gz, &in[vold*((*idomain))], izb);}
	    if(t == 0    ) {i7 = x +nxd*y +nxd*ny*z;__pre_tbwd(tbd_send[i7].c, gt, &in[vold*((*idomain))], itb);}
	  }
	}
      }
    }
#ifdef _NO_OMP_SINGLE
#pragma omp master
#else
#pragma omp single nowait
#endif // _NO_OMP_SINGLE
    {
      if (*idomain == 0) {
	// memcpy(xfd_recv, xfd_send, sizeof(double)*12*ny*nz*nt);
	projscd1_t *tmp = xfd_recv;
	xfd_recv = xfd_send;
	xfd_send = tmp;
	xbound(1,8);
      } else {
	//memcpy(xbd_recv, xbd_send, sizeof(double)*12*ny*nz*nt);
	projscd1_t *tmp = xbd_recv;
	xbd_recv = xbd_send;
	xbd_send = tmp;
	xbound(0,8);
      }
      xbound(2,8);
      xbound(3,8);
      xbound(4,8);
      xbound(5,8);
      xbound(6,8);
      xbound(7,8);
    }
  }

  //---------------------------------------------------------------------------------------- postprocess mult D for boundary
  void ddd_out_pos_d_(scd_t* out, scd_t* in, int* idomain){
    __attribute__((aligned(64)))  scd_t tmp;
    int x, y, z, t, i2, i3, i4, i5, i6, i7;
    pglud_t gxf = &glud[vold*0 + NDIM*vold*   (*idomain)];
    //pglud_t gxb = &glud[vold*0 + NDIM*vold*(1-(*idomain))];
    pglud_t gy = &glud[vold*1 + NDIM*vold*(*idomain)];
    pglud_t gz = &glud[vold*2 + NDIM*vold*(*idomain)];
    pglud_t gt = &glud[vold*3 + NDIM*vold*(*idomain)];
    //i2 = 0 ;  i3 = 0 ;  i4 = 0 ;  i5 = 0 ;  i6 = 0 ;  i7 = 0 ;
#ifdef _NO_OMP_SINGLE
#pragma omp master
#else
#pragma omp single
#endif // _NO_OMP_SINGLE
    {
      if (*idomain == 0) {
	xbound_wait(1,8);
      } else {
	xbound_wait(0,8);
      }
      //      xbound_wait(0,8);
      //      xbound_wait(1,8);
      xbound_wait(2,8);
      xbound_wait(3,8);
      xbound_wait(4,8);
      xbound_wait(5,8);
      xbound_wait(6,8);
      xbound_wait(7,8);
    }
#ifdef _NO_OMP_SINGLE
#pragma omp barrier
#endif // _NO_OMP_SINGLE

#pragma omp parallel for private(t,z,y,x,tmp,i2,i3,i4,i5,i6,i7) collapse(3)
    for(t=0; t<nt; t++){
      for(z=0; z<nz; z++){
	for(y=0; y<ny; y++){
	  //int it = nxd*ny*nz*t;
	  //int iz = nxd*ny*z;
	  //int iy = nxd*y;
	  int ix = y + ny*z + ny*nz*t;
	  for(x=0; x<nxd; x++){
	    if(x == nxd-1 || y == ny-1 || z == nz-1 || t == nt-1 || x==0 || y == 0 || z == 0 || t == 0 ){
	      int i0  = x + nxd*y + nxd*ny*z + nxd*ny*nz*t;
	      __zero_sc(tmp.c);
	      if(x == nxd-1) {	                      __pos_xfwd(tmp.c, gxf,xfd_recv[ix], i0);}
	      if(y == ny-1 ) {i2 = x +nxd*z +nxd*nz*t;__pos_yfwd(tmp.c, gy, yfd_recv[i2], i0);}
	      if(z == nz-1 ) {i4 = x +nxd*y +nxd*ny*t;__pos_zfwd(tmp.c, gz, zfd_recv[i4], i0);}
	      if(t == nt-1 ) {i6 = x +nxd*y +nxd*ny*z;__pos_tfwd_bc(tmp.c, gt, tfd_recv[i6], i0, fbc[3][0]);}
	      if(x == 0    ) {             	      __s_xb1(tmp.c, xbd_recv[ix].c);}
	      if(y == 0    ) {i3 = x +nxd*z +nxd*nz*t;__s_yb(tmp.c, ybd_recv[i3].c);}
	      if(z == 0    ) {i5 = x +nxd*y +nxd*ny*t;__s_zb(tmp.c, zbd_recv[i5].c);}
	      if(t == 0    ) {i7 = x +nxd*y +nxd*ny*z;__s_tb_bc(tmp.c, tbd_recv[i7].c, fbc[3][1]);}
	      __mult_clvd( tmp.cv, clvd[i0 + vold*(*idomain)].cv);
	      __store_add_sc(out[i0].c, tmp.c, mkappa);
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
      xbound_send_waitall(8);
      xbound_recv_okall(8);
    }
#ifdef _NO_OMP_SINGLE
#pragma omp barrier
#endif // _NO_OMP_SINGLE
  }
#ifdef __cplusplus
}
#endif
