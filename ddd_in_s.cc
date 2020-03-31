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
#include "clover_s.h"
#include "mult_all.h"
#include "prefetch.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#include "util.hh"
#include "addressing.hh"

#ifdef __cplusplus
extern "C"{
#endif

  extern __attribute__((aligned(_ALIGN_SIZE))) pglus_t glus;
  extern __attribute__((aligned(_ALIGN_SIZE))) pclvs_t clvs;
  extern double kappa, mkappa;
  extern int nxs, ny, nz, nt;
  extern int vols;

  //---------------------------------------------------------------------------------------- mult D in bulk
  void ddd_in_s_noprl(scs_t* __restrict__ out, const scs_t* __restrict__ in, const int* DEO) {

    // for not to reload after function call
    const float mkappa = ::mkappa;
#ifdef COMPILE_TIME_DIM_SIZE
    const int nxs = NXS, ny = NY, nz = NZ, nt = NT;
    const int vols = VOLS;
#else
    const int nxs = ::nxs, ny = ::ny, nz = ::nz, nt = ::nt;
    const int vols = ::vols;
#endif

    const glus_t* __restrict__ gx = &glus[vols*0 + NDIM*vols*(*DEO)];
    const glus_t* __restrict__ gy = &glus[vols*1 + NDIM*vols*(*DEO)];
    const glus_t* __restrict__ gz = &glus[vols*2 + NDIM*vols*(*DEO)];
    const glus_t* __restrict__ gt = &glus[vols*3 + NDIM*vols*(*DEO)];
    const clvs_t* __restrict__ cl = &clvs[              vols*(*DEO)];

#pragma omp for collapse(3)
    for(int t = 0; t < nt; ++t) {
      for(int z = 0; z < nz; ++z) {
        for(int y = 0; y < ny; ++y) {
          int i0  = nxs*y + nxs*ny*z + nxs*ny*nz*t;
          int izf = nxs*y + (z+1)*nxs*ny + nxs*ny*nz*t;
          int izb = nxs*y + (z-1)*nxs*ny + nxs*ny*nz*t;
          int itf = nxs*y + nxs*ny*z + (t+1)*nxs*ny*nz;
          int itb = nxs*y + nxs*ny*z + (t-1)*nxs*ny*nz;

          __define_preds;

          for(int x = 0; x < nxs; ++x) {
            int ixf = i0+1;
            int ixb = i0-1;
            int iyf = i0+nxs;
            int iyb = i0-nxs;

            if ( y != ny-1 ) {
              const scs_t* __restrict__ in_iyf = in+iyf;
              const glus_t* __restrict__ gy_i0 = gy+i0;
              __prefetch_inp(in_iyf, 0);
              __prefetch_su3(gy_i0, 0);
            }

            //
            // X-forward
            //
            {
#define S(C, S, RI) vecs_t ff_ ## C ## _ ## S ## _ ## RI;
              LOOP_3(LOOP_4, LOOP_2, S);
#undef S
              if ( x != nxs-1 ) {
                const float* __restrict__ in_i0  = ((const float *)(in+i0)) + 1;
                const float* __restrict__ in_ixf = ((const float *)(in+ixf)) + 1-VLENS;
#define S(C, S, RI)                                                     \
                ff_ ## C ## _ ## S ## _ ## RI =                         \
                  for_s(pt,                                             \
                        fload1_s(pt_except_highest, in_i0, dims_scs, C, S, RI), \
                        fload1_s(pt_highest, in_ixf, dims_scs, C, S, RI) \
                        );
                LOOP_3(LOOP_4, LOOP_2, S);
#undef S
              } else {
                const float* __restrict__ in_i0  = ((const float *)(in+i0)) + 1;
#define S(C, S, RI)                                                     \
                ff_ ## C ## _ ## S ## _ ## RI =                         \
                  fload1_s(pt_except_highest, in_i0, dims_scs, C, S, RI);
                LOOP_3(LOOP_4, LOOP_2, S);
#undef S
              }
#define S(C, X, Y)                                                      \
              vecs_t a_ ## C ## _ ## X ## _ ## Y, ua_ ## C ## _ ## X ## _ ## Y;
              LOOP_3(LOOP_2, LOOP_2a, S);
#undef S
              const glus_t* __restrict__ gx_i0  = gx+i0;
              scs_t* __restrict__        out_i0 = out+i0;
              __mult_x_forw_pre_2_vec_(a,ff); // with forward simd-site-shift
              __mult_u_y_vec_(ua,a,(*(gx_i0)));
              __mult_x_forw_pst_vec_((*out_i0),ua);
            }

            if ( y != 0 ) {
              const scs_t* __restrict__  in_iyb = in+iyb;
              const glus_t* __restrict__ gy_iyb = gy+iyb;
              __prefetch_inp(in_iyb, 0);
              __prefetch_su3(gy_iyb, 0);
            }

            //
            // X-backward
            //
            {
              if ( x != 0 ) {
                const float* __restrict__  in_ixb = ((const float *)(in+ixb)) + VLENS-1;
                const float* __restrict__  in_i0  = ((const float *)(in+i0)) + -1;
                const glus_t* __restrict__ gx_i0  = gx+i0;
                const glus_t* __restrict__ gx_ixb = gx+ixb;
                scs_t* __restrict__        out_i0 = out+i0;
#define S(C, S, RI) vecs_t ff_in_ ## C ## _ ## S ## _ ## RI;
                LOOP_3(LOOP_4, LOOP_2, S);
#undef S
#define S(C, S, RI)                                                     \
                ff_in_ ## C ## _ ## S ## _ ## RI =                      \
                  for_s(pt,                                             \
                        fload1_s(pt_lowest, in_ixb, dims_scs, C, S, RI), \
                        fload1_s(pt_except_lowest, in_i0, dims_scs, C, S, RI) \
                        );
                LOOP_3(LOOP_4, LOOP_2, S);
#undef S
#define S(C, X, Y)                                                      \
                vecs_t a_ ## C ## _ ## X ## _ ## Y, ua_ ## C ## _ ## X ## _ ## Y;
                LOOP_3(LOOP_2, LOOP_2a, S);
#undef S
                __mult_x_back_pre_2_vec_(a,ff_in);// with backward simd-site-shift
                __mult_udag_y_2_vec_(ua,a,gx_ixb->c,gx_i0,mid);  // with backward simd-site-shift
                __mult_x_back_pst_vec_((*out_i0),ua);
              } else {
                const float* __restrict__  in_i0  = ((const float *)(in+i0)) + -1;
                const glus_t* __restrict__ gx_i0  = gx+i0;
                scs_t* __restrict__        out_i0 = out+i0;
#define S(C, S, RI) vecs_t ff_in_ ## C ## _ ## S ## _ ## RI;
                LOOP_3(LOOP_4, LOOP_2, S);
#undef S
#define S(C, S, RI)                                                     \
                ff_in_ ## C ## _ ## S ## _ ## RI =                      \
                  fload1_s(pt_except_lowest, in_i0, dims_scs, C, S, RI);
                LOOP_3(LOOP_4, LOOP_2, S);
#undef S
#define S(C, X, Y)                                                      \
                vecs_t a_ ## C ## _ ## X ## _ ## Y, ua_ ## C ## _ ## X ## _ ## Y;
                LOOP_3(LOOP_2, LOOP_2a, S);
#undef S
                __mult_x_back_pre_2_vec_(a,ff_in);// with backward simd-site-shift
                __mult_udag_y_2_vec_(ua,a,,gx_i0->c,edge);  // with backward simd-site-shift
                __mult_x_back_pst_vec_((*out_i0),ua);
              }
            }

            // invalidate temporary values based on i0 and out to suppress register spills
            ASSUME_MODIFIED(i0);
            ASSUME_MODIFIED(out);

            if ( z != nz-1 ) {
              const scs_t* __restrict__  in_izf = in+izf;
              const glus_t* __restrict__ gz_i0  = gz+i0;
              __prefetch_inp(in_izf, 0);
              __prefetch_su3(gz_i0, 0);
            }
            
            //
            // Y-forward
            //
            if ( y != ny-1 ) {
              const glus_t* __restrict__ gy_i0  = gy+i0;
              const scs_t* __restrict__  in_iyf = in + iyf;
              scs_t* __restrict__        out_i0 = out+i0;
#define S(C, X, Y)                                                      \
              vecs_t a_ ## C ## _ ## X ## _ ## Y, ua_ ## C ## _ ## X ## _ ## Y;
              LOOP_3(LOOP_2, LOOP_2a, S);
#undef S
              __mult_y_forw_pre_vec_(a,(*(in_iyf)));
              __mult_u_y_vec_(ua,a,(*(gy_i0)));
              __mult_y_forw_pst_vec_((*out_i0),ua);
            } 

            if ( z != 0 ) {
              const scs_t* __restrict__  in_izb = in+izb;
              const glus_t* __restrict__ gz_izb = gz+izb;
              __prefetch_inp(in_izb, 0);
              __prefetch_su3(gz_izb, 0);
            }

            //
            // Y-backward
            //
            if ( y != 0 ) {
              const glus_t* __restrict__ gy_iyb = gy+iyb;
              const scs_t* __restrict__  in_iyb = in+iyb;
              scs_t* __restrict__        out_i0 = out+i0;
#define S(C, X, Y)                                                      \
              vecs_t a_ ## C ## _ ## X ## _ ## Y, ua_ ## C ## _ ## X ## _ ## Y;
              LOOP_3(LOOP_2, LOOP_2a, S);
#undef S
              __mult_y_back_pre_vec_(a,(*(in_iyb)));
              __mult_udag_y_vec_(ua,a,(*(gy_iyb)));
              __mult_y_back_pst_vec_((*out_i0),ua);
            } 

            if ( t != nt-1 ) {
              const scs_t* __restrict__  in_itf = in+itf;
              const glus_t* __restrict__ gt_i0  = gt+i0;
              __prefetch_inp(in_itf, 0);
              __prefetch_su3(gt_i0, 0);
            }

            //
            // Z-forward
            //
            if ( z != nz-1 ) {
              const glus_t* __restrict__ gz_i0  = gz+i0;
              const scs_t* __restrict__  in_izf = in+izf;
              scs_t* __restrict__        out_i0 = out+i0;
#define S(C, X, Y)                                                      \
              vecs_t a_ ## C ## _ ## X ## _ ## Y, ua_ ## C ## _ ## X ## _ ## Y;
              LOOP_3(LOOP_2, LOOP_2a, S);
#undef S
              __mult_z_forw_pre_vec_(a,(*(in_izf)));
              __mult_u_y_vec_(ua,a,(*(gz_i0)));
              __mult_z_forw_pst_vec_((*out_i0),ua);
            }

            if ( t != 0 ) {
              const scs_t* __restrict__  in_itb = in+itb;
              const glus_t* __restrict__ gt_itb = gt+itb;
              __prefetch_inp(in_itb, 0);
              __prefetch_su3(gt_itb, 0);
            }

            //
            // Z-backward
            //
            if ( z != 0 ) {
              const glus_t* __restrict__ gz_izb = gz+izb;
              const scs_t* __restrict__  in_izb = in+izb;
              scs_t* __restrict__        out_i0 = out+i0;
#define S(C, X, Y)                                                      \
              vecs_t a_ ## C ## _ ## X ## _ ## Y, ua_ ## C ## _ ## X ## _ ## Y;
              LOOP_3(LOOP_2, LOOP_2a, S);
#undef S
              __mult_z_back_pre_vec_(a,(*(in_izb)));
              __mult_udag_y_vec_(ua,a,(*(gz_izb)));
              __mult_z_back_pst_vec_((*out_i0),ua);
            } 

            {
              const clvs_t* __restrict__ cl_i0 = cl+i0;
              __prefetch_clv(cl_i0, 0);
            }

            //
            // T-forward
            //
            if ( t != nt-1 ) { 
              const glus_t* __restrict__ gt_i0  = gt+i0;
              const scs_t* __restrict__  in_itf = in+itf;
              scs_t* __restrict__        out_i0 = out+i0;
#define S(C, X, Y)                                                      \
              vecs_t a_ ## C ## _ ## X ## _ ## Y, ua_ ## C ## _ ## X ## _ ## Y;
              LOOP_3(LOOP_2, LOOP_2a, S);
#undef S
              __mult_t_forw_pre_vec_(a,(*(in_itf)));
              __mult_u_y_vec_(ua,a,(*(gt_i0)));
              __mult_t_forw_pst_vec_((*out_i0),ua);
            } 

            //
            // T-backward
            //
            if ( t != 0 ) { 
              const glus_t* __restrict__ gt_itb = gt+itb;
              const scs_t* __restrict__  in_itb = in+itb;
              scs_t* __restrict__        out_i0 = out+i0;
#define S(C, X, Y)                                                      \
              vecs_t a_ ## C ## _ ## X ## _ ## Y, ua_ ## C ## _ ## X ## _ ## Y;
              LOOP_3(LOOP_2, LOOP_2a, S);
#undef S
              __mult_t_back_pre_vec_(a,(*(in_itb)));
              __mult_udag_y_vec_(ua,a,(*(gt_itb)));
              __mult_t_back_pst_vec_((*out_i0),ua);
            }

            {
              const scs_t* __restrict__  in_ixf = in+ixf+1;
              const glus_t* __restrict__ gx_ixf = gx+ixf+1;
              scs_t* __restrict__        out_i0 = out+i0+1;
              __prefetch_inp(in_ixf, 0);
              __prefetch_su3(gx_ixf, 0);
              __prefetch_out(out_i0, 0);
            }
            ASSUME_MODIFIED(out);
            
            __mult_clvs_vec_m(out[i0].cv, cl[i0].cv);

            {
#define S(C, S, RI)                                                     \
              fstore1_s(pt,                                             \
                        fmadd_s(pt, fload1_s(pt, out_i0, dims_scs, C, S, RI), fdup_s(mkappa), fload1_s(pt, in_i0, dims_scs, C, S, RI)), \
                        out_i0, dims_scs, C, S, RI);
              const scs_t* __restrict__ in_i0 = in+i0;
              scs_t* __restrict__       out_i0 = out+i0;
              LOOP_3(LOOP_4, LOOP_2, S);
#undef S
            }

            i0 ++;
            izf++;
            izb++;
            itf++;
            itb++;
          }
        }
      }
    }
  }

  void ddd_in_s_(scs_t* __restrict__ out, const scs_t* __restrict__ in, const int* DEO) {
#pragma omp parallel
    {
      ddd_in_s_noprl(out, in, DEO);
    }
  }

#ifdef __cplusplus
}
#endif
