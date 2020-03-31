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
#include "qwsintrin.h"
#include "prefetch.h"
#include "util.hh"
#include "addressing.hh"
#include <math.h>

#ifdef __cplusplus
extern "C"{
#endif

  extern int vols;
  extern scs_t *qs;
  extern void ddd_in_s_noprl(scs_t* out, scs_t* in, int *DEO);

  void jinv_ddd_in_s_noprl_(scs_t* x, scs_t* b, int *DEO, int* maxiter){
#ifdef COMPILE_TIME_DIM_SIZE
    const int vols = VOLS;
#else
    // for not to reload after function call
    const int vols = ::vols;
#endif
    scs_t * qs = ::qs;

    pred_t pt = pred_true_all();

    //rvd0 = fload1_s((float)2);
    ddd_in_s_noprl(qs, b, DEO);
#pragma omp for
    for(int i=0; i<vols; i++){
      // for(j=0; j<24; j++){
      //   for(int v=0; v < VLENS; v++) {
      //     x[i].ccs[j].v[v] = 2*b[i].ccs[j].v[v] - qs[i].ccs[j].v[v];
#define S(dummy,idx)                                                    \
      __builtin_prefetch(&(b[i].ccs[(idx)+16]), 0, 3);                    \
      __builtin_prefetch(&(qs[i].ccs[(idx)+16]), 0, 3);                   \
      __builtin_prefetch(&(x[i].ccs[(idx)+16]), 1, 3);                    \
      fstore1_s(pt,                                                     \
                fmsub_s(pt, fdup_s(2), fload1_s(pt, b+i, dims_1d, (idx)*4+0), fload1_s(pt, qs+i, dims_1d, (idx)*4+0)), \
                x+i, dims_1d, (idx)*4+0);                                 \
      fstore1_s(pt,                                                     \
                fmsub_s(pt, fdup_s(2), fload1_s(pt, b+i, dims_1d, (idx)*4+1), fload1_s(pt, qs+i, dims_1d, (idx)*4+1)), \
                x+i, dims_1d, (idx)*4+1);                                 \
      fstore1_s(pt,                                                     \
                fmsub_s(pt, fdup_s(2), fload1_s(pt, b+i, dims_1d, (idx)*4+2), fload1_s(pt, qs+i, dims_1d, (idx)*4+2)), \
                x+i, dims_1d, (idx)*4+2);                                 \
      fstore1_s(pt,                                                     \
                fmsub_s(pt, fdup_s(2), fload1_s(pt, b+i, dims_1d, (idx)*4+3), fload1_s(pt, qs+i, dims_1d, (idx)*4+3)), \
                x+i, dims_1d, (idx)*4+3);
      LOOP_6(S);
#undef S
    }

    for (int iter=1; iter<(*maxiter);iter++){
      // q = Ax
      ddd_in_s_noprl(qs, x, DEO);
      // x = x + b - q
#pragma omp for
      for(int i=0; i<vols; i++){
	// for(j=0; j<24; j++)
	//   for(int v=0; v < VLENS; v++)
	//     x[i].ccs[j].v[v] += b[i].ccs[j].v[v] - qs[i].ccs[j].v[v];
        ASSUME_MODIFIED(x);
        ASSUME_MODIFIED(b);
        ASSUME_MODIFIED(qs);
    
#define S(dummy,idx)                                                    \
        __builtin_prefetch(&(b[i].ccs[(idx)+16]), 0, 3);                  \
        __builtin_prefetch(&(qs[i].ccs[(idx)+16]), 0, 3);                 \
        __builtin_prefetch(&(x[i].ccs[(idx)+16]), 1, 3);                  \
        fstore1_s(pt,                                                   \
                  fadd_s(pt, fload1_s(pt, x+i, dims_1d, (idx)*4+0),       \
                         fsub_s(pt, fload1_s(pt, b+i, dims_1d, (idx)*4+0), fload1_s(pt, qs+i, dims_1d, (idx)*4+0))), \
                  x+i, dims_1d, (idx)*4+0);                               \
        fstore1_s(pt,                                                   \
                  fadd_s(pt, fload1_s(pt, x+i, dims_1d, (idx)*4+1),       \
                         fsub_s(pt, fload1_s(pt, b+i, dims_1d, (idx)*4+1), fload1_s(pt, qs+i, dims_1d, (idx)*4+1))), \
                  x+i, dims_1d, (idx)*4+1);                               \
        fstore1_s(pt,                                                   \
                  fadd_s(pt, fload1_s(pt, x+i, dims_1d, (idx)*4+2),       \
                         fsub_s(pt, fload1_s(pt, b+i, dims_1d, (idx)*4+2), fload1_s(pt, qs+i, dims_1d, (idx)*4+2))), \
                  x+i, dims_1d, (idx)*4+2);                               \
        fstore1_s(pt,                                                   \
                  fadd_s(pt, fload1_s(pt, x+i, dims_1d, (idx)*4+3),       \
                         fsub_s(pt, fload1_s(pt, b+i, dims_1d, (idx)*4+3), fload1_s(pt, qs+i, dims_1d, (idx)*4+3))), \
                  x+i, dims_1d, (idx)*4+3);
        LOOP_6(S);
#undef S
      }
    }//iter
  }//mr
  void jinv_ddd_in_s_(scs_t* x, scs_t* b, int *DEO, int* maxiter){
#pragma omp parallel
    {
      jinv_ddd_in_s_noprl_(x, b, DEO, maxiter);
    }
  }
#ifdef __cplusplus
}
#endif
