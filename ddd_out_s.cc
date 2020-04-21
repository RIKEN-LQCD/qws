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
#include "mult_all.h"
#include "util.hh"
#include "prefetch.h"
#include <omp.h>

#ifdef __cplusplus
extern "C"{
#endif

#include "timing.h"

  extern pglus_t __restrict__ glus __attribute__((aligned(_ALIGN_SIZE)));
  extern pclvs_t __restrict__ clvs __attribute__((aligned(_ALIGN_SIZE)));
  extern projscs1_t * __restrict__ xfs_send;
  extern projscs1_t * __restrict__ xfs_recv;
  extern projscs1_t * __restrict__ xbs_send, * __restrict__ xbs_recv;
  extern projscs_t  * __restrict__ yfs_send, * __restrict__ yfs_recv;
  extern projscs_t  * __restrict__ ybs_send, * __restrict__ ybs_recv;
  extern projscs_t  * __restrict__ zfs_send, * __restrict__ zfs_recv;
  extern projscs_t  * __restrict__ zbs_send, * __restrict__ zbs_recv;
  extern projscs_t  * __restrict__ tfs_send, * __restrict__ tfs_recv;
  extern projscs_t  * __restrict__ tbs_send, * __restrict__ tbs_recv;

  extern double kappa, mkappa;
  extern int nxs, ny, nz, nt;
  extern int vols;
  extern int thmax;
  extern int iam;
#pragma omp threadprivate(iam)
  extern int pt;
  extern int npe[4];
  extern double fbc[4][2];

  void xbound(int req, int prec);
  void xbound_wait(int req, int prec);
  void xbound_reset_comm(int req, int prec);
  void xbound_send_waitall(int prec);
  void xbound_recv_okall(int prec);
  void xbound_flip_parity(int prec);
  void xbound_recv_updateall(int prec);
  void xbound_recv_update(int req, int prec);

  void ddd_out_pre_s_(scs_t* __restrict__  in, int* __restrict__  idomain);
  void ddd_out_pre_s_no_timer_(scs_t* __restrict__  in, int* __restrict__  idomain);
  void ddd_out_pre_s_noprl_(scs_t* __restrict__  in, int* __restrict__  idomain);
  void ddd_out_pre_s_noprl_no_timer_(scs_t* __restrict__  in, int* __restrict__  idomain);
  void ddd_out_pos_s_(scs_t* __restrict__ out, scs_t* __restrict__ in, int* idomain, float factor);
  void ddd_out_pos_s_no_timer_(scs_t* __restrict__ out, scs_t* __restrict__ in, int* idomain, float factor);
  void ddd_out_pos_s_noprl_(scs_t* __restrict__ out, scs_t* __restrict__ in, int* idomain, float factor);
  void ddd_out_pos_s_noprl_no_timer_(scs_t* __restrict__ out, scs_t* __restrict__ in, int* idomain, float factor);
#ifdef __cplusplus
}
#endif

//---------------------------------------------------------------------------------------- preprocess mult D for boundary
//
// pack data for send
//
template<bool is_timer_enabled>
static inline void ddd_out_pre_s_noprl(scs_t* __restrict__  in, int* __restrict__  idomain) {

  if(is_timer_enabled) {
    _S_MULT_WD_DEO_OUT_SEND_HPC_TIC_;
  }

#ifdef COMPILE_TIME_DIM_SIZE
    const long int nxs = NXS, ny = NY, nz = NZ, nt = NT;
    const long int vols = VOLS;
#else
    const long int nxs = ::nxs, ny = ::ny, nz = ::nz, nt = ::nt;
    const long int vols = ::vols;
#endif

  long int xdmn, ydmn, zdmn, tdmn;

  xdmn = 1 - *idomain;
  if (npe[1] == 1) { ydmn = *idomain; } else { ydmn = 1 - *idomain; }
  if (npe[2] == 1) { zdmn = *idomain; } else { zdmn = 1 - *idomain; }
  if (npe[3] == 1) { tdmn = *idomain; } else { tdmn = 1 - *idomain; }

  scs_t *inxdmn = in + vols*xdmn;
  scs_t *inydmn = in + vols*ydmn;
  scs_t *inzdmn = in + vols*zdmn;
  scs_t *intdmn = in + vols*tdmn;

  glus_t *gx = &glus[vols*0 + NDIM*vols*xdmn];
  glus_t *gy = &glus[vols*1 + NDIM*vols*ydmn];
  glus_t *gz = &glus[vols*2 + NDIM*vols*zdmn];
  glus_t *gt = &glus[vols*3 + NDIM*vols*tdmn];

  projscs1_t *xfs_send = ::xfs_send;
  projscs1_t *xbs_send = ::xbs_send;
  projscs_t  *yfs_send = ::yfs_send;
  projscs_t  *ybs_send = ::ybs_send;
  projscs_t  *zfs_send = ::zfs_send;
  projscs_t  *zbs_send = ::zbs_send;
  projscs_t  *tfs_send = ::tfs_send;
  projscs_t  *tbs_send = ::tbs_send;

#pragma omp for collapse(3)
  for(long int t=0; t<nt; t++){
    for(long int z=0; z<nz; z++){
      for(long int y=0; y<ny; y++){
        long int ix = nxs*y + nxs*ny*z + nxs*ny*nz*t;
        long int iy =         nxs*ny*z + nxs*ny*nz*t;
        long int iz = nxs*y +            nxs*ny*nz*t;
        long int it = nxs*y + nxs*ny*z;
        long int iwox =     y +     ny*z +     ny*nz*t;
        long int iwoy =         nxs*   z + nxs*   nz*t;
        long int iwoz = nxs*y +            nxs*ny*   t;
        long int iwot = nxs*y + nxs*ny*z              ;
        scs_t *inxdmn_i = inxdmn + ix;
        scs_t *inydmn_i = inydmn + iy;
        scs_t *inzdmn_i = inzdmn + iz;
        scs_t *intdmn_i = intdmn + it;
        glus_t *gx_i = gx + ix + nxs-1;
        glus_t *gy_i = gy + iy + (ny-1)*nxs;
        glus_t *gz_i = gz + iz + (nz-1)*nxs*ny;
        glus_t *gt_i = gt + it + (nt-1)*nxs*ny*nz;
        projscs1_t *xfs_send_i = xfs_send + iwox;
        projscs1_t *xbs_send_i = xbs_send + iwox;
        projscs_t  *yfs_send_i = yfs_send + iwoy;
        projscs_t  *ybs_send_i = ybs_send + iwoy;
        projscs_t  *zfs_send_i = zfs_send + iwoz;
        projscs_t  *zbs_send_i = zbs_send + iwoz;
        projscs_t  *tfs_send_i = tfs_send + iwot;
        projscs_t  *tbs_send_i = tbs_send + iwot;

        __define_preds;

#define S(C, X, Y)                                                      \
        vecs_t p_ ## C ## _ ## X ## _ ## Y, p2_ ## C ## _ ## X ## _ ## Y;
        LOOP_3(LOOP_2, LOOP_2a, S);
#undef S
              
        for(long int x=0; x<nxs; x++){
          if(x == nxs-1) {
            __mult_x_forw_pre_3_(*xfs_send_i,*inxdmn_i);
            __prefetch_inp(inxdmn_i, 1);
          }
          if(y == ny-1 ) {
            __mult_y_forw_pre_vec_(p,*inydmn_i);
            __store_projscs_vec_(yfs_send_i->c,p);
          }
          if(z == nz-1 ) {
            __prefetch_inp(inzdmn_i, 1);
            __mult_z_forw_pre_vec_(p,*inzdmn_i);
            __store_projscs_vec_(zfs_send_i->c,p);
            __prefetch_sendbuf(zfs_send_i, 1);
          }
          if(t == nt-1 ) {
            __prefetch_inp(intdmn_i, 1);
            __mult_t_forw_pre_vec_(p,*intdmn_i);
            __store_projscs_vec_(tfs_send_i->c,p);
            __prefetch_sendbuf(tfs_send_i, 1);
          }
          ASSUME_MODIFIED(inxdmn_i); 
          ASSUME_MODIFIED(inydmn_i); 
          ASSUME_MODIFIED(inzdmn_i); 
          ASSUME_MODIFIED(intdmn_i); 
          ASSUME_MODIFIED(gx_i);
          ASSUME_MODIFIED(gy_i);
          ASSUME_MODIFIED(gz_i);
          ASSUME_MODIFIED(gt_i);
          ASSUME_MODIFIED(xfs_send_i);
          ASSUME_MODIFIED(xbs_send_i);
          ASSUME_MODIFIED(yfs_send_i);
          ASSUME_MODIFIED(ybs_send_i);
          ASSUME_MODIFIED(zfs_send_i);
          ASSUME_MODIFIED(zbs_send_i);
          ASSUME_MODIFIED(tfs_send_i);
          ASSUME_MODIFIED(tbs_send_i);
          if(x == 0    ) {
            scs_t *intmp = inxdmn_i + nxs-1;
            projscs1_t a;
            __prefetch_inp(intmp, 1);
            __mult_x_back_pre_3_(a,*intmp);
            __mult_udag_y_3_(*xbs_send_i,a,*gx_i);
            __prefetch_su3(gx_i, 1);
          }
          if(y == 0    ) {
            scs_t *intmp = inydmn_i + (ny-1)*nxs;
            __mult_y_back_pre_vec_(p,*intmp);
            __mult_udag_y_vec_(p2,p,*gy_i);
            __store_projscs_vec_(ybs_send_i->c,p2);
          }
          if(z == 0    ) {
            scs_t *intmp = inzdmn_i + (nz-1)*nxs*ny;
            __prefetch_inp(intmp, 1);
            __mult_z_back_pre_vec_(p,*intmp);
            __mult_udag_y_vec_(p2,p,*gz_i);
            __store_projscs_vec_(zbs_send_i->c,p2);
            __prefetch_su3(gz_i, 1);
            __prefetch_sendbuf(zbs_send_i, 1);
          }
          if(t == 0    ) {
            scs_t *intmp = intdmn_i + (nt-1)*nxs*ny*nz;
            __prefetch_inp(intmp, 1);
            __mult_t_back_pre_vec_(p,*intmp);
            __mult_udag_y_vec_(p2,p,*gt_i);
            __store_projscs_vec_(tbs_send_i->c,p2);
            __prefetch_su3(gt_i, 1);
            __prefetch_sendbuf(tbs_send_i, 1);
          }
          inzdmn_i++;
          intdmn_i++;
          inydmn_i++;
          gy_i++;
          gz_i++;
          gt_i++;
          yfs_send_i++;
          ybs_send_i++;
          zfs_send_i++;
          zbs_send_i++;
          tfs_send_i++;
          tbs_send_i++;
        }
      }
    }
  }

  if(is_timer_enabled) {
    _S_MULT_WD_DEO_OUT_SEND_HPC_TOC_;
    _COMM_TIC_;
    _COMLIB_ISEND_ALL_C_TIC_;
  }
#ifdef _THREADED_RDMA
#pragma omp for nowait
  for(int req=0; req<8; req++) {
  xbound_recv_update(req, 4); // update the pointer to buffer before memcpy
    if(req==0){
      if (*idomain == 0) {
        xbound_reset_comm(0,4); // set a flag to supress send_check
	memcpy(xfs_recv, xfs_send, sizeof(float)*12*ny*nz*nt);
      } else {
	xbound(0,4);
      }
    } else if(req==1){
      if (*idomain == 0) {
	xbound(1,4);
      } else {
  	xbound_reset_comm(1,4); // set a flag to supress send_check
	memcpy(xbs_recv, xbs_send, sizeof(float)*12*ny*nz*nt);
      }
    } else { // req=2,3,...,7
      xbound(req,4);
    }
  }
#else // _THREADED_RDMA
#ifdef _NO_OMP_SINGLE
#pragma omp master
#else
#pragma omp single nowait
#endif // _NO_OMP_SINGLE
  {
    xbound_recv_updateall(4); // update the pointer to buffer before memcpy
    if (*idomain == 0) {
      xbound_reset_comm(0,4); // set a flag to supress send_check
      memcpy(xfs_recv, xfs_send, sizeof(float)*12*ny*nz*nt);
      xbound(1,4);
    } else {
      xbound_reset_comm(1,4); // set a flag to supress send_check
      memcpy(xbs_recv, xbs_send, sizeof(float)*12*ny*nz*nt);
      xbound(0,4);
    }
    xbound(2,4);
    xbound(3,4);
    xbound(4,4);
    xbound(5,4);
    xbound(6,4);
    xbound(7,4);
  }
#endif // _THREADED_RDMA
  if(is_timer_enabled) {
#pragma omp barrier
    _COMLIB_ISEND_ALL_C_TOC_;
    _COMM_TOC_;
    _OVERLAPPED_CALC_TIC_;
    _COMLIB_IRECV_ALL_C_TIC_;
  }
}

//---------------------------------------------------------------------------------------- postprocess mult D for boundary
template<bool is_timer_enabled>
static inline void ddd_out_pos_s_noprl(scs_t* __restrict__ out, scs_t* __restrict__ in, int* idomain, float factor) {

#ifdef COMPILE_TIME_DIM_SIZE
    const long int nxs = NXS, ny = NY, nz = NZ, nt = NT;
    const long int vols = VOLS;
#else
    const long int nxs = ::nxs, ny = ::ny, nz = ::nz, nt = ::nt;
    const long int vols = ::vols;
#endif
  if(is_timer_enabled) {
    _OVERLAPPED_CALC_TOC_;
    _COMM_TIC_;
  }

#if 0 //#ifdef _THREADED_RDMA
#pragma omp for
  for(int req=0; req<8; req++) {
    if(req==0){
      if (*idomain == 1) {
	xbound_wait(req,4);
      }
    } else if(req==1){
      if (*idomain == 0) {
	xbound_wait(req,4);
      }
    } else { // req=2,3,...,7
      xbound_wait(req,4);
    }
  }
#else // #ifdef _THREADED_RDMA
#ifdef _NO_OMP_SINGLE
#pragma omp master
#else
#pragma omp single
#endif // _NO_OMP_SINGLE
  {

    if (*idomain == 0) {
      xbound_wait(1,4);
    } else {
      xbound_wait(0,4);
    }
    //    xbound_wait(0,4);
    //    xbound_wait(1,4);
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
#endif // _THREADED_RDMA
  if(is_timer_enabled) {
    _COMLIB_IRECV_ALL_C_TOC_;
    _COMM_TOC_;
    _S_MULT_WD_DEO_OUT_RECV_HPC_CALC_TIC_;
  }

  glus_t *gx = &glus[vols*0 + NDIM*vols*(*idomain)];
  glus_t *gy = &glus[vols*1 + NDIM*vols*(*idomain)];
  glus_t *gz = &glus[vols*2 + NDIM*vols*(*idomain)];
  glus_t *gt = &glus[vols*3 + NDIM*vols*(*idomain)];
  clvs_t *clvs_idmn = clvs + vols*(*idomain);
  projscs1_t *xfs_recv = ::xfs_recv;
  projscs1_t *xbs_recv = ::xbs_recv;
  projscs_t  *yfs_recv = ::yfs_recv;
  projscs_t  *ybs_recv = ::ybs_recv;
  projscs_t  *zfs_recv = ::zfs_recv;
  projscs_t  *zbs_recv = ::zbs_recv;
  projscs_t  *tfs_recv = ::tfs_recv;
  projscs_t  *tbs_recv = ::tbs_recv;

  float tbc_fwd = ((float)fbc[3][0])*0.5f;
  float tbc_bwd = ((float)fbc[3][1])*0.5f;

#pragma omp for collapse(3)
  for (long int t = 0; t < nt; t++) {
    for (long int z = 0; z < nz; z++) {
      for (long int y = 0; y < ny; y++) {
        long int iwox =     y +     ny*z +     ny*nz*t;
        long int iwoy =         nxs*   z + nxs*   nz*t;
        long int iwoz = nxs*y +            nxs*ny*   t;
        long int iwot = nxs*y + nxs*ny*z              ;
        long int i0 =   nxs*y + nxs*ny*z + nxs*ny*nz*t;
        glus_t *gx_i = gx + i0;
        glus_t *gy_i = gy + i0;
        glus_t *gz_i = gz + i0;
        glus_t *gt_i = gt + i0;
        clvs_t *clvs_i = clvs_idmn + i0;
        scs_t *out_i = out + i0;
        projscs1_t *xfs_recv_i = xfs_recv + iwox;
        projscs1_t *xbs_recv_i = xbs_recv + iwox;
        projscs_t  *yfs_recv_i = yfs_recv + iwoy;
        projscs_t  *ybs_recv_i = ybs_recv + iwoy;
        projscs_t  *zfs_recv_i = zfs_recv + iwoz;
        projscs_t  *zbs_recv_i = zbs_recv + iwoz;
        projscs_t  *tfs_recv_i = tfs_recv + iwot;
        projscs_t  *tbs_recv_i = tbs_recv + iwot;

        __define_preds;

        for (long int x = 0; x < nxs; x++) {

          if (x == nxs-1 || y == ny-1 || z == nz-1 || t == nt-1 || x==0 || y == 0 || z == 0 || t == 0 ) {

#define S(C, X, Y)                                                      \
            vecs_t p_ ## C ## _ ## X ## _ ## Y, p2_ ## C ## _ ## X ## _ ## Y;
            LOOP_3(LOOP_2, LOOP_2a, S);
#undef S

#define S(C, S, R)                                      \
            vecs_t tmpv_ ## C ## _ ## S ## _ ## R;
            LOOP_3(LOOP_4, LOOP_2, S);
#undef S

            scs_t *tmp;
            alloca_aligned(tmp, sizeof(scs_t), CLS);

            if ( x == nxs-1 ) {
              __prefetch_su3(gx_i, 0);
            }

#define S(C, S, RI)                                             \
            fstore1_s(pt, fdup_s(0), tmp, dims_scs, C, S, RI);
            LOOP_3(LOOP_4, LOOP_2, S);
#undef S

            if ( y == ny-1 ) {
              __prefetch_su3(gy_i, 0);
              __prefetch_recvbuf(yfs_recv_i, 0);
            }

            //
            // X-forward
            //
            if ( x == nxs-1 ) {
              __attribute__((aligned(_ALIGN_SIZE))) projscs1_t ua;
              __mult_u_y_3_(ua,*xfs_recv_i,*gx_i);
              __mult_x_forw_pst_3_(*tmp,ua);
            }

            if ( z == nz-1 ) {
              __prefetch_su3(gz_i, 0);
              __prefetch_recvbuf(zfs_recv_i, 0);
            }

            //
            // Y-forward 
            //
            if ( y == ny-1 ) {
              __load_projscs_vec_(p, yfs_recv_i->c);
              __mult_u_y_vec_(p2,p,*gy_i);
              __mult_y_forw_pst_vec_((*tmp),p2);
            }

            if ( t == nt-1 ) {
              __prefetch_su3(gt_i, 0);
              __prefetch_recvbuf(tfs_recv_i, 0);
            }

            //
            // Z-forward
            //
            if ( z == nz-1 ) {
              __load_projscs_vec_(p, zfs_recv_i->c);
              __mult_u_y_vec_(p2,p,*gz_i);
              __mult_z_forw_pst_vec_((*tmp),p2);
            }

            //
            // T-forward
            //
            if ( t == nt-1 ) {
              __load_projscs_vec_(p, tfs_recv_i->c);
              __mult_u_y_vec_(p2,p,*gt_i);
              __mult_t_forw_pst_bc_vec_((*tmp),p2,tbc_fwd);
            }

            ASSUME_MODIFIED(gx_i);
            ASSUME_MODIFIED(gy_i);
            ASSUME_MODIFIED(gz_i);
            ASSUME_MODIFIED(gt_i);
            ASSUME_MODIFIED(clvs_i);
            ASSUME_MODIFIED(out_i);
            ASSUME_MODIFIED(xfs_recv_i);
            ASSUME_MODIFIED(xbs_recv_i);
            ASSUME_MODIFIED(yfs_recv_i);
            ASSUME_MODIFIED(ybs_recv_i);
            ASSUME_MODIFIED(zfs_recv_i);
            ASSUME_MODIFIED(zbs_recv_i);
            ASSUME_MODIFIED(tfs_recv_i);
            ASSUME_MODIFIED(tbs_recv_i);
            ASSUME_MODIFIED(tmp);

            __load_scs_vec_(tmpv, tmp);
            ASSUME_MODIFIED(tmp);

            if ( y == 0 ) {
              __prefetch_recvbuf(ybs_recv_i, 0);
            }

            //
            // X-backward
            //
            if ( x == 0 ) {
              __mult_x_back_pst_3_vec_(tmpv,*xbs_recv_i);
            }

            if ( z == 0 ) {
              __prefetch_recvbuf(zbs_recv_i, 0);
            }

            //
            // Y-backward
            //
            if ( y == 0 ) {
              __mult_y_back_pst_vec_reg_(tmpv, *ybs_recv_i);
            }

            if ( t == 0 ) {
              __prefetch_recvbuf(tbs_recv_i, 0);
            }

            __prefetch_clv(clvs_i, 0);
            
            //
            // Z-backward
            //
            if ( z == 0 ) {
              __mult_z_back_pst_vec_reg_(tmpv, *zbs_recv_i);
            }

            //
            // T-backward
            //
            if ( t == 0 ) {
              __mult_t_back_pst_bc_vec_reg_(tmpv,*tbs_recv_i,tbc_bwd);
            }

            __mult_clvs_vec_m_reg( (*tmp).cv, tmpv, clvs_i->cv);
            ASSUME_MODIFIED(tmp);

            __prefetch_out(out_i, 0);

#define S(C, S, RI)                                                     \
            fstore1_s(pt,                                               \
              fmadd_s(pt, fdup_s(factor), fload1_s(pt, tmp, dims_scs, C, S, RI), fload1_s(pt, out_i, dims_scs, C, S, RI)), \
              out_i, dims_scs, C, S, RI);
            LOOP_3(LOOP_4, LOOP_2, S);
#undef S
          } // if edge sites
          gx_i++;
          gy_i++;
          gz_i++;
          gt_i++;
          clvs_i++;
          out_i++;
          yfs_recv_i++;
          ybs_recv_i++;
          zfs_recv_i++;
          zbs_recv_i++;
          tfs_recv_i++;
          tbs_recv_i++;
        }
      }
    }
  }
  if(is_timer_enabled) {
    _S_MULT_WD_DEO_OUT_RECV_HPC_CALC_TOC_;
    _COMM_TIC_;
    _COMLIB_SEND_WAIT_ALL_C_TIC_;
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
#endif // _NO_OMP_SINGLE
  if(is_timer_enabled) {
    _COMLIB_SEND_WAIT_ALL_C_TOC_;
    _COMM_TOC_;
  }
}

void ddd_out_pre_s_noprl_(scs_t* __restrict__  in, int* __restrict__  idomain) {
  ddd_out_pre_s_noprl<true>(in, idomain);
}

void ddd_out_pre_s_noprl_no_timer_(scs_t* __restrict__  in, int* __restrict__  idomain) {
  ddd_out_pre_s_noprl<false>(in, idomain);
}
  
void ddd_out_pre_s_(scs_t* __restrict__  in, int* __restrict__  idomain) {
#pragma omp parallel
  {
    ddd_out_pre_s_noprl_(in, idomain);
  }
}

void ddd_out_pre_s_no_timer_(scs_t* __restrict__  in, int* __restrict__  idomain) {
#pragma omp parallel
  {
    ddd_out_pre_s_noprl_no_timer_(in, idomain);
  }
}

void ddd_out_pos_s_noprl_(scs_t* __restrict__ out, scs_t* __restrict__ in, int* idomain, float factor) {
  ddd_out_pos_s_noprl<true>(out, in, idomain, factor);
}

void ddd_out_pos_s_noprl_no_timer_(scs_t* __restrict__ out, scs_t* __restrict__ in, int* idomain, float factor) {
  ddd_out_pos_s_noprl<false>(out, in, idomain, factor);
}

void ddd_out_pos_s_(scs_t* __restrict__ out, scs_t* __restrict__ in, int* idomain, float factor) {
#pragma omp parallel
  {
    ddd_out_pos_s_noprl_(out, in, idomain, factor);
  }
}

void ddd_out_pos_s_no_timer_(scs_t* __restrict__ out, scs_t* __restrict__ in, int* idomain, float factor) {
#pragma omp parallel
  {
    ddd_out_pos_s_noprl_no_timer_(out, in, idomain, factor);
  }
}
