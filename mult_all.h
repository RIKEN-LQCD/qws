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
#ifndef _MULT_ALL_H
#define _MULT_ALL_H

#define _ALIGN_SIZE  (4)

#include"util.hh"
#include"qwsintrin.h"

#define __define_preds                                          \
  pred_t pt = pred_true_all();                                  \
  pred_t pt_except_lowest = pred_true_except_lowest();          \
  pred_t pt_lowest = pred_true_lowest();                        \
  pred_t pt_except_highest = pred_true_except_highest();        \
  pred_t pt_highest = pred_true_highest();

#define __store_projscs_vec_body_(dest, src, C, X, Y) \
  fstore1_s(pt, src##_##C##_##X##_##Y, dest, dims_projscs, C, X, Y);
#define __store_projscs_vec_(dest, src) \
  LOOP_3(LOOP_2, LOOP_2a, __store_projscs_vec_body_, dest, src)

#define __load_projscs_vec_body_(pred, dest, src, C, X, Y) \
  dest##_##C##_##X##_##Y = fload1_s(pred, src, dims_projscs, C, X, Y);
#define __load_projscs_vec_pred_(pred, dest, src) \
  LOOP_3(LOOP_2, LOOP_2a, __load_projscs_vec_body_, pred, dest, src)
#define __load_projscs_vec_(dest, src) __load_projscs_vec_pred_(pt, dest, src)

#define __load_scs_vec_body_(dest, src, C, S, R) \
  dest##_##C##_##S##_##R = fload1_s(pt, src, dims_scs, C, S, R);
#define __load_scs_vec_(dest, src) \
  LOOP_3(LOOP_4, LOOP_2, __load_scs_vec_body_, dest, src)

//
// link multiplication on two-spinor
// upy1 = u * y1
// upy2 = u * y2
//
#define __mult_u_y_(upy,py,u) {\
  for (int c = 0; c < 3; ++c) {\
  for (int s = 0; s < 2; ++s) {\
    for (int j = 0; j < VLENS; ++j) {\
      (upy).c[c][s][0][j]  = (u).c[0][c][0][j] * (py).c[0][s][0][j];\
      (upy).c[c][s][1][j]  = (u).c[0][c][0][j] * (py).c[0][s][1][j];\
      (upy).c[c][s][0][j] += (u).c[1][c][0][j] * (py).c[1][s][0][j];\
      (upy).c[c][s][1][j] += (u).c[1][c][0][j] * (py).c[1][s][1][j];\
      (upy).c[c][s][0][j] += (u).c[2][c][0][j] * (py).c[2][s][0][j];\
      (upy).c[c][s][1][j] += (u).c[2][c][0][j] * (py).c[2][s][1][j];\
      (upy).c[c][s][0][j] -= (u).c[0][c][1][j] * (py).c[0][s][1][j];\
      (upy).c[c][s][1][j] += (u).c[0][c][1][j] * (py).c[0][s][0][j];\
      (upy).c[c][s][0][j] -= (u).c[1][c][1][j] * (py).c[1][s][1][j];\
      (upy).c[c][s][1][j] += (u).c[1][c][1][j] * (py).c[1][s][0][j];\
      (upy).c[c][s][0][j] -= (u).c[2][c][1][j] * (py).c[2][s][1][j];\
      (upy).c[c][s][1][j] += (u).c[2][c][1][j] * (py).c[2][s][0][j];\
    }\
  }\
  }\
}

#define __mult_u_y_vec_(upy,py,u) {                             \
    vecs_t tmp0, tmp1, tmp2;                                    \
    tmp0 = fload1_s(pt, (u).c, dims_glus, 0, 0, 0);             \
    tmp1 = fload1_s(pt, (u).c, dims_glus, 0, 1, 0);             \
    tmp2 = fload1_s(pt, (u).c, dims_glus, 0, 2, 0);             \
    upy##_0_0_0 = fmul_s(pt, tmp0, py##_0_0_0);                 \
    upy##_0_0_1 = fmul_s(pt, tmp0, py##_0_0_1);                 \
    upy##_0_1_0 = fmul_s(pt, tmp0, py##_0_1_0);                 \
    upy##_0_1_1 = fmul_s(pt, tmp0, py##_0_1_1);                 \
    upy##_1_0_0 = fmul_s(pt, tmp1, py##_0_0_0);                 \
    upy##_1_0_1 = fmul_s(pt, tmp1, py##_0_0_1);                 \
    upy##_1_1_0 = fmul_s(pt, tmp1, py##_0_1_0);                 \
    upy##_1_1_1 = fmul_s(pt, tmp1, py##_0_1_1);                 \
    upy##_2_0_0 = fmul_s(pt, tmp2, py##_0_0_0);                 \
    upy##_2_0_1 = fmul_s(pt, tmp2, py##_0_0_1);                 \
    upy##_2_1_0 = fmul_s(pt, tmp2, py##_0_1_0);                 \
    upy##_2_1_1 = fmul_s(pt, tmp2, py##_0_1_1);                 \
                                                                \
    tmp0 = fload1_s(pt, (u).c, dims_glus, 1, 0, 0);             \
    tmp1 = fload1_s(pt, (u).c, dims_glus, 1, 1, 0);             \
    tmp2 = fload1_s(pt, (u).c, dims_glus, 1, 2, 0);             \
    upy##_0_0_0 = fmadd_s(pt, tmp0, py##_1_0_0, upy##_0_0_0);   \
    upy##_0_0_1 = fmadd_s(pt, tmp0, py##_1_0_1, upy##_0_0_1);   \
    upy##_0_1_0 = fmadd_s(pt, tmp0, py##_1_1_0, upy##_0_1_0);   \
    upy##_0_1_1 = fmadd_s(pt, tmp0, py##_1_1_1, upy##_0_1_1);   \
    upy##_1_0_0 = fmadd_s(pt, tmp1, py##_1_0_0, upy##_1_0_0);   \
    upy##_1_0_1 = fmadd_s(pt, tmp1, py##_1_0_1, upy##_1_0_1);   \
    upy##_1_1_0 = fmadd_s(pt, tmp1, py##_1_1_0, upy##_1_1_0);   \
    upy##_1_1_1 = fmadd_s(pt, tmp1, py##_1_1_1, upy##_1_1_1);   \
    upy##_2_0_0 = fmadd_s(pt, tmp2, py##_1_0_0, upy##_2_0_0);   \
    upy##_2_0_1 = fmadd_s(pt, tmp2, py##_1_0_1, upy##_2_0_1);   \
    upy##_2_1_0 = fmadd_s(pt, tmp2, py##_1_1_0, upy##_2_1_0);   \
    upy##_2_1_1 = fmadd_s(pt, tmp2, py##_1_1_1, upy##_2_1_1);   \
                                                                \
    tmp0 = fload1_s(pt, (u).c, dims_glus, 2, 0, 0);             \
    tmp1 = fload1_s(pt, (u).c, dims_glus, 2, 1, 0);             \
    tmp2 = fload1_s(pt, (u).c, dims_glus, 2, 2, 0);             \
    upy##_0_0_0 = fmadd_s(pt, tmp0, py##_2_0_0, upy##_0_0_0);   \
    upy##_0_0_1 = fmadd_s(pt, tmp0, py##_2_0_1, upy##_0_0_1);   \
    upy##_0_1_0 = fmadd_s(pt, tmp0, py##_2_1_0, upy##_0_1_0);   \
    upy##_0_1_1 = fmadd_s(pt, tmp0, py##_2_1_1, upy##_0_1_1);   \
    upy##_1_0_0 = fmadd_s(pt, tmp1, py##_2_0_0, upy##_1_0_0);   \
    upy##_1_0_1 = fmadd_s(pt, tmp1, py##_2_0_1, upy##_1_0_1);   \
    upy##_1_1_0 = fmadd_s(pt, tmp1, py##_2_1_0, upy##_1_1_0);   \
    upy##_1_1_1 = fmadd_s(pt, tmp1, py##_2_1_1, upy##_1_1_1);   \
    upy##_2_0_0 = fmadd_s(pt, tmp2, py##_2_0_0, upy##_2_0_0);   \
    upy##_2_0_1 = fmadd_s(pt, tmp2, py##_2_0_1, upy##_2_0_1);   \
    upy##_2_1_0 = fmadd_s(pt, tmp2, py##_2_1_0, upy##_2_1_0);   \
    upy##_2_1_1 = fmadd_s(pt, tmp2, py##_2_1_1, upy##_2_1_1);   \
                                                                \
    tmp0 = fload1_s(pt, (u).c, dims_glus, 0, 0, 1);             \
    tmp1 = fload1_s(pt, (u).c, dims_glus, 0, 1, 1);             \
    tmp2 = fload1_s(pt, (u).c, dims_glus, 0, 2, 1);             \
    upy##_0_0_0 = fnmadd_s(pt, tmp0, py##_0_0_1, upy##_0_0_0);  \
    upy##_0_0_1 = fmadd_s(pt, tmp0, py##_0_0_0, upy##_0_0_1);   \
    upy##_0_1_0 = fnmadd_s(pt, tmp0, py##_0_1_1, upy##_0_1_0);  \
    upy##_0_1_1 = fmadd_s(pt, tmp0, py##_0_1_0, upy##_0_1_1);   \
    upy##_1_0_0 = fnmadd_s(pt, tmp1, py##_0_0_1, upy##_1_0_0);  \
    upy##_1_0_1 = fmadd_s(pt, tmp1, py##_0_0_0, upy##_1_0_1);   \
    upy##_1_1_0 = fnmadd_s(pt, tmp1, py##_0_1_1, upy##_1_1_0);  \
    upy##_1_1_1 = fmadd_s(pt, tmp1, py##_0_1_0, upy##_1_1_1);   \
    upy##_2_0_0 = fnmadd_s(pt, tmp2, py##_0_0_1, upy##_2_0_0);  \
    upy##_2_0_1 = fmadd_s(pt, tmp2, py##_0_0_0, upy##_2_0_1);   \
    upy##_2_1_0 = fnmadd_s(pt, tmp2, py##_0_1_1, upy##_2_1_0);  \
    upy##_2_1_1 = fmadd_s(pt, tmp2, py##_0_1_0, upy##_2_1_1);   \
                                                                \
    tmp0 = fload1_s(pt, (u).c, dims_glus, 1, 0, 1);             \
    tmp1 = fload1_s(pt, (u).c, dims_glus, 1, 1, 1);             \
    tmp2 = fload1_s(pt, (u).c, dims_glus, 1, 2, 1);             \
    upy##_0_0_0 = fnmadd_s(pt, tmp0, py##_1_0_1, upy##_0_0_0);  \
    upy##_0_0_1 = fmadd_s(pt, tmp0, py##_1_0_0, upy##_0_0_1);   \
    upy##_0_1_0 = fnmadd_s(pt, tmp0, py##_1_1_1, upy##_0_1_0);  \
    upy##_0_1_1 = fmadd_s(pt, tmp0, py##_1_1_0, upy##_0_1_1);   \
    upy##_1_0_0 = fnmadd_s(pt, tmp1, py##_1_0_1, upy##_1_0_0);  \
    upy##_1_0_1 = fmadd_s(pt, tmp1, py##_1_0_0, upy##_1_0_1);   \
    upy##_1_1_0 = fnmadd_s(pt, tmp1, py##_1_1_1, upy##_1_1_0);  \
    upy##_1_1_1 = fmadd_s(pt, tmp1, py##_1_1_0, upy##_1_1_1);   \
    upy##_2_0_0 = fnmadd_s(pt, tmp2, py##_1_0_1, upy##_2_0_0);  \
    upy##_2_0_1 = fmadd_s(pt, tmp2, py##_1_0_0, upy##_2_0_1);   \
    upy##_2_1_0 = fnmadd_s(pt, tmp2, py##_1_1_1, upy##_2_1_0);  \
    upy##_2_1_1 = fmadd_s(pt, tmp2, py##_1_1_0, upy##_2_1_1);   \
                                                                \
    tmp0 = fload1_s(pt, (u).c, dims_glus, 2, 0, 1);             \
    tmp1 = fload1_s(pt, (u).c, dims_glus, 2, 1, 1);             \
    tmp2 = fload1_s(pt, (u).c, dims_glus, 2, 2, 1);             \
    upy##_0_0_0 = fnmadd_s(pt, tmp0, py##_2_0_1, upy##_0_0_0);  \
    upy##_0_0_1 = fmadd_s(pt, tmp0, py##_2_0_0, upy##_0_0_1);   \
    upy##_0_1_0 = fnmadd_s(pt, tmp0, py##_2_1_1, upy##_0_1_0);  \
    upy##_0_1_1 = fmadd_s(pt, tmp0, py##_2_1_0, upy##_0_1_1);   \
    upy##_1_0_0 = fnmadd_s(pt, tmp1, py##_2_0_1, upy##_1_0_0);  \
    upy##_1_0_1 = fmadd_s(pt, tmp1, py##_2_0_0, upy##_1_0_1);   \
    upy##_1_1_0 = fnmadd_s(pt, tmp1, py##_2_1_1, upy##_1_1_0);  \
    upy##_1_1_1 = fmadd_s(pt, tmp1, py##_2_1_0, upy##_1_1_1);   \
    upy##_2_0_0 = fnmadd_s(pt, tmp2, py##_2_0_1, upy##_2_0_0);  \
    upy##_2_0_1 = fmadd_s(pt, tmp2, py##_2_0_0, upy##_2_0_1);   \
    upy##_2_1_0 = fnmadd_s(pt, tmp2, py##_2_1_1, upy##_2_1_0);  \
    upy##_2_1_1 = fmadd_s(pt, tmp2, py##_2_1_0, upy##_2_1_1);   \
  }

//
// link multiplication on two-spinor on a site with simd-site last site link
// upy1 = u * y1
// upy2 = u * y2
//
/* #define __mult_u_y_3_(upy,py,u) {                                       \ */
/*     for (int c = 0; c < 3; ++c) {                                       \ */
/*       for (int s = 0; s < 2; ++s) {                                     \ */
/*         (upy).c[c][s][0]  = (u).c[0][c][0][VLENS-1] * (py).c[0][s][0];  \ */
/*         (upy).c[c][s][1]  = (u).c[0][c][0][VLENS-1] * (py).c[0][s][1];  \ */
/*         (upy).c[c][s][0] += (u).c[1][c][0][VLENS-1] * (py).c[1][s][0];  \ */
/*         (upy).c[c][s][1] += (u).c[1][c][0][VLENS-1] * (py).c[1][s][1];  \ */
/*         (upy).c[c][s][0] += (u).c[2][c][0][VLENS-1] * (py).c[2][s][0];  \ */
/*         (upy).c[c][s][1] += (u).c[2][c][0][VLENS-1] * (py).c[2][s][1];  \ */
/*         (upy).c[c][s][0] -= (u).c[0][c][1][VLENS-1] * (py).c[0][s][1];  \ */
/*         (upy).c[c][s][1] += (u).c[0][c][1][VLENS-1] * (py).c[0][s][0];  \ */
/*         (upy).c[c][s][0] -= (u).c[1][c][1][VLENS-1] * (py).c[1][s][1];  \ */
/*         (upy).c[c][s][1] += (u).c[1][c][1][VLENS-1] * (py).c[1][s][0];  \ */
/*         (upy).c[c][s][0] -= (u).c[2][c][1][VLENS-1] * (py).c[2][s][1];  \ */
/*         (upy).c[c][s][1] += (u).c[2][c][1][VLENS-1] * (py).c[2][s][0];  \ */
/*       }                                                                 \ */
/*     }                                                                   \ */
/*   } */

#define __mult_u_y_3_(upy,py,u)                                         \
  {                                                                     \
  (upy).c[0][0][0]  = (u).c[0][0][0][VLENS-1] * (py).c[0][0][0];        \
  (upy).c[0][0][1]  = (u).c[0][0][0][VLENS-1] * (py).c[0][0][1];        \
  (upy).c[0][1][0]  = (u).c[0][0][0][VLENS-1] * (py).c[0][1][0];        \
  (upy).c[0][1][1]  = (u).c[0][0][0][VLENS-1] * (py).c[0][1][1];        \
  (upy).c[1][0][0]  = (u).c[0][1][0][VLENS-1] * (py).c[0][0][0];        \
  (upy).c[1][0][1]  = (u).c[0][1][0][VLENS-1] * (py).c[0][0][1];        \
  (upy).c[1][1][0]  = (u).c[0][1][0][VLENS-1] * (py).c[0][1][0];        \
  (upy).c[1][1][1]  = (u).c[0][1][0][VLENS-1] * (py).c[0][1][1];        \
  (upy).c[2][0][0]  = (u).c[0][2][0][VLENS-1] * (py).c[0][0][0];        \
  (upy).c[2][0][1]  = (u).c[0][2][0][VLENS-1] * (py).c[0][0][1];        \
  (upy).c[2][1][0]  = (u).c[0][2][0][VLENS-1] * (py).c[0][1][0];        \
  (upy).c[2][1][1]  = (u).c[0][2][0][VLENS-1] * (py).c[0][1][1];        \
  (upy).c[0][0][0] += (u).c[1][0][0][VLENS-1] * (py).c[1][0][0];        \
  (upy).c[0][0][1] += (u).c[1][0][0][VLENS-1] * (py).c[1][0][1];        \
  (upy).c[0][1][0] += (u).c[1][0][0][VLENS-1] * (py).c[1][1][0];        \
  (upy).c[0][1][1] += (u).c[1][0][0][VLENS-1] * (py).c[1][1][1];        \
  (upy).c[1][0][0] += (u).c[1][1][0][VLENS-1] * (py).c[1][0][0];        \
  (upy).c[1][0][1] += (u).c[1][1][0][VLENS-1] * (py).c[1][0][1];        \
  (upy).c[1][1][0] += (u).c[1][1][0][VLENS-1] * (py).c[1][1][0];        \
  (upy).c[1][1][1] += (u).c[1][1][0][VLENS-1] * (py).c[1][1][1];        \
  (upy).c[2][0][0] += (u).c[1][2][0][VLENS-1] * (py).c[1][0][0];        \
  (upy).c[2][0][1] += (u).c[1][2][0][VLENS-1] * (py).c[1][0][1];        \
  (upy).c[2][1][0] += (u).c[1][2][0][VLENS-1] * (py).c[1][1][0];        \
  (upy).c[2][1][1] += (u).c[1][2][0][VLENS-1] * (py).c[1][1][1];        \
  (upy).c[0][0][0] += (u).c[2][0][0][VLENS-1] * (py).c[2][0][0];        \
  (upy).c[0][0][1] += (u).c[2][0][0][VLENS-1] * (py).c[2][0][1];        \
  (upy).c[0][1][0] += (u).c[2][0][0][VLENS-1] * (py).c[2][1][0];        \
  (upy).c[0][1][1] += (u).c[2][0][0][VLENS-1] * (py).c[2][1][1];        \
  (upy).c[1][0][0] += (u).c[2][1][0][VLENS-1] * (py).c[2][0][0];        \
  (upy).c[1][0][1] += (u).c[2][1][0][VLENS-1] * (py).c[2][0][1];        \
  (upy).c[1][1][0] += (u).c[2][1][0][VLENS-1] * (py).c[2][1][0];        \
  (upy).c[1][1][1] += (u).c[2][1][0][VLENS-1] * (py).c[2][1][1];        \
  (upy).c[2][0][0] += (u).c[2][2][0][VLENS-1] * (py).c[2][0][0];        \
  (upy).c[2][0][1] += (u).c[2][2][0][VLENS-1] * (py).c[2][0][1];        \
  (upy).c[2][1][0] += (u).c[2][2][0][VLENS-1] * (py).c[2][1][0];        \
  (upy).c[2][1][1] += (u).c[2][2][0][VLENS-1] * (py).c[2][1][1];        \
  (upy).c[0][0][0] -= (u).c[0][0][1][VLENS-1] * (py).c[0][0][1];        \
  (upy).c[0][0][1] += (u).c[0][0][1][VLENS-1] * (py).c[0][0][0];        \
  (upy).c[0][1][0] -= (u).c[0][0][1][VLENS-1] * (py).c[0][1][1];        \
  (upy).c[0][1][1] += (u).c[0][0][1][VLENS-1] * (py).c[0][1][0];        \
  (upy).c[1][0][0] -= (u).c[0][1][1][VLENS-1] * (py).c[0][0][1];        \
  (upy).c[1][0][1] += (u).c[0][1][1][VLENS-1] * (py).c[0][0][0];        \
  (upy).c[1][1][0] -= (u).c[0][1][1][VLENS-1] * (py).c[0][1][1];        \
  (upy).c[1][1][1] += (u).c[0][1][1][VLENS-1] * (py).c[0][1][0];        \
  (upy).c[2][0][0] -= (u).c[0][2][1][VLENS-1] * (py).c[0][0][1];        \
  (upy).c[2][0][1] += (u).c[0][2][1][VLENS-1] * (py).c[0][0][0];        \
  (upy).c[2][1][0] -= (u).c[0][2][1][VLENS-1] * (py).c[0][1][1];        \
  (upy).c[2][1][1] += (u).c[0][2][1][VLENS-1] * (py).c[0][1][0];        \
  (upy).c[0][0][0] -= (u).c[1][0][1][VLENS-1] * (py).c[1][0][1];        \
  (upy).c[0][0][1] += (u).c[1][0][1][VLENS-1] * (py).c[1][0][0];        \
  (upy).c[0][1][0] -= (u).c[1][0][1][VLENS-1] * (py).c[1][1][1];        \
  (upy).c[0][1][1] += (u).c[1][0][1][VLENS-1] * (py).c[1][1][0];        \
  (upy).c[1][0][0] -= (u).c[1][1][1][VLENS-1] * (py).c[1][0][1];        \
  (upy).c[1][0][1] += (u).c[1][1][1][VLENS-1] * (py).c[1][0][0];        \
  (upy).c[1][1][0] -= (u).c[1][1][1][VLENS-1] * (py).c[1][1][1];        \
  (upy).c[1][1][1] += (u).c[1][1][1][VLENS-1] * (py).c[1][1][0];        \
  (upy).c[2][0][0] -= (u).c[1][2][1][VLENS-1] * (py).c[1][0][1];        \
  (upy).c[2][0][1] += (u).c[1][2][1][VLENS-1] * (py).c[1][0][0];        \
  (upy).c[2][1][0] -= (u).c[1][2][1][VLENS-1] * (py).c[1][1][1];        \
  (upy).c[2][1][1] += (u).c[1][2][1][VLENS-1] * (py).c[1][1][0];        \
  (upy).c[0][0][0] -= (u).c[2][0][1][VLENS-1] * (py).c[2][0][1];        \
  (upy).c[0][0][1] += (u).c[2][0][1][VLENS-1] * (py).c[2][0][0];        \
  (upy).c[0][1][0] -= (u).c[2][0][1][VLENS-1] * (py).c[2][1][1];        \
  (upy).c[0][1][1] += (u).c[2][0][1][VLENS-1] * (py).c[2][1][0];        \
  (upy).c[1][0][0] -= (u).c[2][1][1][VLENS-1] * (py).c[2][0][1];        \
  (upy).c[1][0][1] += (u).c[2][1][1][VLENS-1] * (py).c[2][0][0];        \
  (upy).c[1][1][0] -= (u).c[2][1][1][VLENS-1] * (py).c[2][1][1];        \
  (upy).c[1][1][1] += (u).c[2][1][1][VLENS-1] * (py).c[2][1][0];        \
  (upy).c[2][0][0] -= (u).c[2][2][1][VLENS-1] * (py).c[2][0][1];        \
  (upy).c[2][0][1] += (u).c[2][2][1][VLENS-1] * (py).c[2][0][0];        \
  (upy).c[2][1][0] -= (u).c[2][2][1][VLENS-1] * (py).c[2][1][1];        \
  (upy).c[2][1][1] += (u).c[2][2][1][VLENS-1] * (py).c[2][1][0];        \
  }

//
// link multiplication on two-spinor
// upy1 = udag * y1
// upy2 = udag * y2
//
#define __mult_udag_y_(upy,py,u) {\
  for (int c = 0; c < 3; ++c) {\
  for (int s = 0; s < 2; ++s) {\
    for (int j = 0; j < VLENS; ++j) {\
      (upy).c[c][s][0][j]  = (u).c[c][0][0][j] * (py).c[0][s][0][j];\
      (upy).c[c][s][1][j]  = (u).c[c][0][0][j] * (py).c[0][s][1][j];\
      (upy).c[c][s][0][j] += (u).c[c][1][0][j] * (py).c[1][s][0][j];\
      (upy).c[c][s][1][j] += (u).c[c][1][0][j] * (py).c[1][s][1][j];\
      (upy).c[c][s][0][j] += (u).c[c][2][0][j] * (py).c[2][s][0][j];\
      (upy).c[c][s][1][j] += (u).c[c][2][0][j] * (py).c[2][s][1][j];\
      (upy).c[c][s][0][j] += (u).c[c][0][1][j] * (py).c[0][s][1][j];\
      (upy).c[c][s][1][j] -= (u).c[c][0][1][j] * (py).c[0][s][0][j];\
      (upy).c[c][s][0][j] += (u).c[c][1][1][j] * (py).c[1][s][1][j];\
      (upy).c[c][s][1][j] -= (u).c[c][1][1][j] * (py).c[1][s][0][j];\
      (upy).c[c][s][0][j] += (u).c[c][2][1][j] * (py).c[2][s][1][j];\
      (upy).c[c][s][1][j] -= (u).c[c][2][1][j] * (py).c[2][s][0][j];\
    }\
  }\
  }\
}

#define __mult_udag_y_vec_(upy,py,u) {                          \
    vecs_t tmp0, tmp1, tmp2;                                    \
    tmp0 = fload1_s(pt, (u).c, dims_glus, 0, 0, 0);             \
    tmp1 = fload1_s(pt, (u).c, dims_glus, 1, 0, 0);             \
    tmp2 = fload1_s(pt, (u).c, dims_glus, 2, 0, 0);             \
    upy##_0_0_0  = fmul_s(pt, tmp0, py##_0_0_0);                \
    upy##_0_0_1  = fmul_s(pt, tmp0, py##_0_0_1);                \
    upy##_0_1_0  = fmul_s(pt, tmp0, py##_0_1_0);                \
    upy##_0_1_1  = fmul_s(pt, tmp0, py##_0_1_1);                \
    upy##_1_0_0  = fmul_s(pt, tmp1, py##_0_0_0);                \
    upy##_1_0_1  = fmul_s(pt, tmp1, py##_0_0_1);                \
    upy##_1_1_0  = fmul_s(pt, tmp1, py##_0_1_0);                \
    upy##_1_1_1  = fmul_s(pt, tmp1, py##_0_1_1);                \
    upy##_2_0_0  = fmul_s(pt, tmp2, py##_0_0_0);                \
    upy##_2_0_1  = fmul_s(pt, tmp2, py##_0_0_1);                \
    upy##_2_1_0  = fmul_s(pt, tmp2, py##_0_1_0);                \
    upy##_2_1_1  = fmul_s(pt, tmp2, py##_0_1_1);                \
                                                                \
    tmp0 = fload1_s(pt, (u).c, dims_glus, 0, 1, 0);             \
    tmp1 = fload1_s(pt, (u).c, dims_glus, 1, 1, 0);             \
    tmp2 = fload1_s(pt, (u).c, dims_glus, 2, 1, 0);             \
    upy##_0_0_0 = fmadd_s(pt, tmp0, py##_1_0_0, upy##_0_0_0);   \
    upy##_0_0_1 = fmadd_s(pt, tmp0, py##_1_0_1, upy##_0_0_1);   \
    upy##_0_1_0 = fmadd_s(pt, tmp0, py##_1_1_0, upy##_0_1_0);   \
    upy##_0_1_1 = fmadd_s(pt, tmp0, py##_1_1_1, upy##_0_1_1);   \
    upy##_1_0_0 = fmadd_s(pt, tmp1, py##_1_0_0, upy##_1_0_0);   \
    upy##_1_0_1 = fmadd_s(pt, tmp1, py##_1_0_1, upy##_1_0_1);   \
    upy##_1_1_0 = fmadd_s(pt, tmp1, py##_1_1_0, upy##_1_1_0);   \
    upy##_1_1_1 = fmadd_s(pt, tmp1, py##_1_1_1, upy##_1_1_1);   \
    upy##_2_0_0 = fmadd_s(pt, tmp2, py##_1_0_0, upy##_2_0_0);   \
    upy##_2_0_1 = fmadd_s(pt, tmp2, py##_1_0_1, upy##_2_0_1);   \
    upy##_2_1_0 = fmadd_s(pt, tmp2, py##_1_1_0, upy##_2_1_0);   \
    upy##_2_1_1 = fmadd_s(pt, tmp2, py##_1_1_1, upy##_2_1_1);   \
                                                                \
    tmp0 = fload1_s(pt, (u).c, dims_glus, 0, 2, 0);             \
    tmp1 = fload1_s(pt, (u).c, dims_glus, 1, 2, 0);             \
    tmp2 = fload1_s(pt, (u).c, dims_glus, 2, 2, 0);             \
    upy##_0_0_0 = fmadd_s(pt, tmp0, py##_2_0_0, upy##_0_0_0);   \
    upy##_0_0_1 = fmadd_s(pt, tmp0, py##_2_0_1, upy##_0_0_1);   \
    upy##_0_1_0 = fmadd_s(pt, tmp0, py##_2_1_0, upy##_0_1_0);   \
    upy##_0_1_1 = fmadd_s(pt, tmp0, py##_2_1_1, upy##_0_1_1);   \
    upy##_1_0_0 = fmadd_s(pt, tmp1, py##_2_0_0, upy##_1_0_0);   \
    upy##_1_0_1 = fmadd_s(pt, tmp1, py##_2_0_1, upy##_1_0_1);   \
    upy##_1_1_0 = fmadd_s(pt, tmp1, py##_2_1_0, upy##_1_1_0);   \
    upy##_1_1_1 = fmadd_s(pt, tmp1, py##_2_1_1, upy##_1_1_1);   \
    upy##_2_0_0 = fmadd_s(pt, tmp2, py##_2_0_0, upy##_2_0_0);   \
    upy##_2_0_1 = fmadd_s(pt, tmp2, py##_2_0_1, upy##_2_0_1);   \
    upy##_2_1_0 = fmadd_s(pt, tmp2, py##_2_1_0, upy##_2_1_0);   \
    upy##_2_1_1 = fmadd_s(pt, tmp2, py##_2_1_1, upy##_2_1_1);   \
                                                                \
    tmp0 = fload1_s(pt, (u).c, dims_glus, 0, 0, 1);             \
    tmp1 = fload1_s(pt, (u).c, dims_glus, 1, 0, 1);             \
    tmp2 = fload1_s(pt, (u).c, dims_glus, 2, 0, 1);             \
    upy##_0_0_0 = fmadd_s(pt, tmp0, py##_0_0_1, upy##_0_0_0);   \
    upy##_0_0_1 = fnmadd_s(pt, tmp0, py##_0_0_0, upy##_0_0_1);  \
    upy##_0_1_0 = fmadd_s(pt, tmp0, py##_0_1_1, upy##_0_1_0);   \
    upy##_0_1_1 = fnmadd_s(pt, tmp0, py##_0_1_0, upy##_0_1_1);  \
    upy##_1_0_0 = fmadd_s(pt, tmp1, py##_0_0_1, upy##_1_0_0);   \
    upy##_1_0_1 = fnmadd_s(pt, tmp1, py##_0_0_0, upy##_1_0_1);  \
    upy##_1_1_0 = fmadd_s(pt, tmp1, py##_0_1_1, upy##_1_1_0);   \
    upy##_1_1_1 = fnmadd_s(pt, tmp1, py##_0_1_0, upy##_1_1_1);  \
    upy##_2_0_0 = fmadd_s(pt, tmp2, py##_0_0_1, upy##_2_0_0);   \
    upy##_2_0_1 = fnmadd_s(pt, tmp2, py##_0_0_0, upy##_2_0_1);  \
    upy##_2_1_0 = fmadd_s(pt, tmp2, py##_0_1_1, upy##_2_1_0);   \
    upy##_2_1_1 = fnmadd_s(pt, tmp2, py##_0_1_0, upy##_2_1_1);  \
                                                                \
    tmp0 = fload1_s(pt, (u).c, dims_glus, 0, 1, 1);             \
    tmp1 = fload1_s(pt, (u).c, dims_glus, 1, 1, 1);             \
    tmp2 = fload1_s(pt, (u).c, dims_glus, 2, 1, 1);             \
    upy##_0_0_0 = fmadd_s(pt, tmp0, py##_1_0_1, upy##_0_0_0);   \
    upy##_0_0_1 = fnmadd_s(pt, tmp0, py##_1_0_0, upy##_0_0_1);  \
    upy##_0_1_0 = fmadd_s(pt, tmp0, py##_1_1_1, upy##_0_1_0);   \
    upy##_0_1_1 = fnmadd_s(pt, tmp0, py##_1_1_0, upy##_0_1_1);  \
    upy##_1_0_0 = fmadd_s(pt, tmp1, py##_1_0_1, upy##_1_0_0);   \
    upy##_1_0_1 = fnmadd_s(pt, tmp1, py##_1_0_0, upy##_1_0_1);  \
    upy##_1_1_0 = fmadd_s(pt, tmp1, py##_1_1_1, upy##_1_1_0);   \
    upy##_1_1_1 = fnmadd_s(pt, tmp1, py##_1_1_0, upy##_1_1_1);  \
    upy##_2_0_0 = fmadd_s(pt, tmp2, py##_1_0_1, upy##_2_0_0);   \
    upy##_2_0_1 = fnmadd_s(pt, tmp2, py##_1_0_0, upy##_2_0_1);  \
    upy##_2_1_0 = fmadd_s(pt, tmp2, py##_1_1_1, upy##_2_1_0);   \
    upy##_2_1_1 = fnmadd_s(pt, tmp2, py##_1_1_0, upy##_2_1_1);  \
                                                                \
    tmp0 = fload1_s(pt, (u).c, dims_glus, 0, 2, 1);             \
    tmp1 = fload1_s(pt, (u).c, dims_glus, 1, 2, 1);             \
    tmp2 = fload1_s(pt, (u).c, dims_glus, 2, 2, 1);             \
    upy##_0_0_0 = fmadd_s(pt, tmp0, py##_2_0_1, upy##_0_0_0);   \
    upy##_0_0_1 = fnmadd_s(pt, tmp0, py##_2_0_0, upy##_0_0_1);  \
    upy##_0_1_0 = fmadd_s(pt, tmp0, py##_2_1_1, upy##_0_1_0);   \
    upy##_0_1_1 = fnmadd_s(pt, tmp0, py##_2_1_0, upy##_0_1_1);  \
    upy##_1_0_0 = fmadd_s(pt, tmp1, py##_2_0_1, upy##_1_0_0);   \
    upy##_1_0_1 = fnmadd_s(pt, tmp1, py##_2_0_0, upy##_1_0_1);  \
    upy##_1_1_0 = fmadd_s(pt, tmp1, py##_2_1_1, upy##_1_1_0);   \
    upy##_1_1_1 = fnmadd_s(pt, tmp1, py##_2_1_0, upy##_1_1_1);  \
    upy##_2_0_0 = fmadd_s(pt, tmp2, py##_2_0_1, upy##_2_0_0);   \
    upy##_2_0_1 = fnmadd_s(pt, tmp2, py##_2_0_0, upy##_2_0_1);  \
    upy##_2_1_0 = fmadd_s(pt, tmp2, py##_2_1_1, upy##_2_1_0);   \
    upy##_2_1_1 = fnmadd_s(pt, tmp2, py##_2_1_0, upy##_2_1_1);  \
  }

//
// link multiplication on two-spinor with simd-site-shifted gauge link
// upy1(j) = udag(j+VLEN-1) * y1(j)
// upy2(j) = udag(j+VLEN-1) * y2(j)
//
#define __mult_udag_y_2_(upy,py,u) {\
  for (int c = 0; c < 3; ++c) {\
  for (int s = 0; s < 2; ++s) {\
    for (int j = 0; j < VLENS; ++j) {\
      (upy).c[c][s][0][j]  = (u)[c][0][0][j+VLENS-1] * (py).c[0][s][0][j];\
      (upy).c[c][s][1][j]  = (u)[c][0][0][j+VLENS-1] * (py).c[0][s][1][j];\
      (upy).c[c][s][0][j] += (u)[c][1][0][j+VLENS-1] * (py).c[1][s][0][j];\
      (upy).c[c][s][1][j] += (u)[c][1][0][j+VLENS-1] * (py).c[1][s][1][j];\
      (upy).c[c][s][0][j] += (u)[c][2][0][j+VLENS-1] * (py).c[2][s][0][j];\
      (upy).c[c][s][1][j] += (u)[c][2][0][j+VLENS-1] * (py).c[2][s][1][j];\
      (upy).c[c][s][0][j] += (u)[c][0][1][j+VLENS-1] * (py).c[0][s][1][j];\
      (upy).c[c][s][1][j] -= (u)[c][0][1][j+VLENS-1] * (py).c[0][s][0][j];\
      (upy).c[c][s][0][j] += (u)[c][1][1][j+VLENS-1] * (py).c[1][s][1][j];\
      (upy).c[c][s][1][j] -= (u)[c][1][1][j+VLENS-1] * (py).c[1][s][0][j];\
      (upy).c[c][s][0][j] += (u)[c][2][1][j+VLENS-1] * (py).c[2][s][1][j];\
      (upy).c[c][s][1][j] -= (u)[c][2][1][j+VLENS-1] * (py).c[2][s][0][j];\
    }\
  }\
  }\
}

#define __load_backward_mid(b, f, ...)                                  \
  for_s(pt,                                                             \
        fload1_s(pt_lowest, ((float*)b)+VLENS-1, dims_glus, __VA_ARGS__), \
        fload1_s(pt_except_lowest, ((float*)f)-1, dims_glus, __VA_ARGS__) \
        )
#define __load_backward_edge(b, f, ...)                                 \
  fload1_s(pt_except_lowest, ((float*)f)-1, dims_glus, __VA_ARGS__)
#define __mult_udag_y_2_vec_(upy,py,u_b,u_f,pos) {              \
    vecs_t tmp0, tmp1, tmp2;                                    \
    tmp0 = __load_backward_##pos((u_b), (u_f), 0, 0, 0);        \
    tmp1 = __load_backward_##pos((u_b), (u_f), 1, 0, 0);        \
    tmp2 = __load_backward_##pos((u_b), (u_f), 2, 0, 0);        \
    upy##_0_0_0  = fmul_s(pt, tmp0, py##_0_0_0);                \
    upy##_0_0_1  = fmul_s(pt, tmp0, py##_0_0_1);                \
    upy##_0_1_0  = fmul_s(pt, tmp0, py##_0_1_0);                \
    upy##_0_1_1  = fmul_s(pt, tmp0, py##_0_1_1);                \
    upy##_1_0_0  = fmul_s(pt, tmp1, py##_0_0_0);                \
    upy##_1_0_1  = fmul_s(pt, tmp1, py##_0_0_1);                \
    upy##_1_1_0  = fmul_s(pt, tmp1, py##_0_1_0);                \
    upy##_1_1_1  = fmul_s(pt, tmp1, py##_0_1_1);                \
    upy##_2_0_0  = fmul_s(pt, tmp2, py##_0_0_0);                \
    upy##_2_0_1  = fmul_s(pt, tmp2, py##_0_0_1);                \
    upy##_2_1_0  = fmul_s(pt, tmp2, py##_0_1_0);                \
    upy##_2_1_1  = fmul_s(pt, tmp2, py##_0_1_1);                \
                                                                \
    tmp0 = __load_backward_##pos((u_b), (u_f), 0, 1, 0);        \
    tmp1 = __load_backward_##pos((u_b), (u_f), 1, 1, 0);        \
    tmp2 = __load_backward_##pos((u_b), (u_f), 2, 1, 0);        \
    upy##_0_0_0 = fmadd_s(pt, tmp0, py##_1_0_0, upy##_0_0_0);   \
    upy##_0_0_1 = fmadd_s(pt, tmp0, py##_1_0_1, upy##_0_0_1);   \
    upy##_0_1_0 = fmadd_s(pt, tmp0, py##_1_1_0, upy##_0_1_0);   \
    upy##_0_1_1 = fmadd_s(pt, tmp0, py##_1_1_1, upy##_0_1_1);   \
    upy##_1_0_0 = fmadd_s(pt, tmp1, py##_1_0_0, upy##_1_0_0);   \
    upy##_1_0_1 = fmadd_s(pt, tmp1, py##_1_0_1, upy##_1_0_1);   \
    upy##_1_1_0 = fmadd_s(pt, tmp1, py##_1_1_0, upy##_1_1_0);   \
    upy##_1_1_1 = fmadd_s(pt, tmp1, py##_1_1_1, upy##_1_1_1);   \
    upy##_2_0_0 = fmadd_s(pt, tmp2, py##_1_0_0, upy##_2_0_0);   \
    upy##_2_0_1 = fmadd_s(pt, tmp2, py##_1_0_1, upy##_2_0_1);   \
    upy##_2_1_0 = fmadd_s(pt, tmp2, py##_1_1_0, upy##_2_1_0);   \
    upy##_2_1_1 = fmadd_s(pt, tmp2, py##_1_1_1, upy##_2_1_1);   \
                                                                \
    tmp0 = __load_backward_##pos((u_b), (u_f), 0, 2, 0);        \
    tmp1 = __load_backward_##pos((u_b), (u_f), 1, 2, 0);        \
    tmp2 = __load_backward_##pos((u_b), (u_f), 2, 2, 0);        \
    upy##_0_0_0 = fmadd_s(pt, tmp0, py##_2_0_0, upy##_0_0_0);   \
    upy##_0_0_1 = fmadd_s(pt, tmp0, py##_2_0_1, upy##_0_0_1);   \
    upy##_0_1_0 = fmadd_s(pt, tmp0, py##_2_1_0, upy##_0_1_0);   \
    upy##_0_1_1 = fmadd_s(pt, tmp0, py##_2_1_1, upy##_0_1_1);   \
    upy##_1_0_0 = fmadd_s(pt, tmp1, py##_2_0_0, upy##_1_0_0);   \
    upy##_1_0_1 = fmadd_s(pt, tmp1, py##_2_0_1, upy##_1_0_1);   \
    upy##_1_1_0 = fmadd_s(pt, tmp1, py##_2_1_0, upy##_1_1_0);   \
    upy##_1_1_1 = fmadd_s(pt, tmp1, py##_2_1_1, upy##_1_1_1);   \
    upy##_2_0_0 = fmadd_s(pt, tmp2, py##_2_0_0, upy##_2_0_0);   \
    upy##_2_0_1 = fmadd_s(pt, tmp2, py##_2_0_1, upy##_2_0_1);   \
    upy##_2_1_0 = fmadd_s(pt, tmp2, py##_2_1_0, upy##_2_1_0);   \
    upy##_2_1_1 = fmadd_s(pt, tmp2, py##_2_1_1, upy##_2_1_1);   \
                                                                \
    tmp0 = __load_backward_##pos((u_b), (u_f), 0, 0, 1);        \
    tmp1 = __load_backward_##pos((u_b), (u_f), 1, 0, 1);        \
    tmp2 = __load_backward_##pos((u_b), (u_f), 2, 0, 1);        \
    upy##_0_0_0 = fmadd_s(pt, tmp0, py##_0_0_1, upy##_0_0_0);   \
    upy##_0_0_1 = fnmadd_s(pt, tmp0, py##_0_0_0, upy##_0_0_1);  \
    upy##_0_1_0 = fmadd_s(pt, tmp0, py##_0_1_1, upy##_0_1_0);   \
    upy##_0_1_1 = fnmadd_s(pt, tmp0, py##_0_1_0, upy##_0_1_1);  \
    upy##_1_0_0 = fmadd_s(pt, tmp1, py##_0_0_1, upy##_1_0_0);   \
    upy##_1_0_1 = fnmadd_s(pt, tmp1, py##_0_0_0, upy##_1_0_1);  \
    upy##_1_1_0 = fmadd_s(pt, tmp1, py##_0_1_1, upy##_1_1_0);   \
    upy##_1_1_1 = fnmadd_s(pt, tmp1, py##_0_1_0, upy##_1_1_1);  \
    upy##_2_0_0 = fmadd_s(pt, tmp2, py##_0_0_1, upy##_2_0_0);   \
    upy##_2_0_1 = fnmadd_s(pt, tmp2, py##_0_0_0, upy##_2_0_1);  \
    upy##_2_1_0 = fmadd_s(pt, tmp2, py##_0_1_1, upy##_2_1_0);   \
    upy##_2_1_1 = fnmadd_s(pt, tmp2, py##_0_1_0, upy##_2_1_1);  \
                                                                \
    tmp0 = __load_backward_##pos((u_b), (u_f), 0, 1, 1);        \
    tmp1 = __load_backward_##pos((u_b), (u_f), 1, 1, 1);        \
    tmp2 = __load_backward_##pos((u_b), (u_f), 2, 1, 1);        \
    upy##_0_0_0 = fmadd_s(pt, tmp0, py##_1_0_1, upy##_0_0_0);   \
    upy##_0_0_1 = fnmadd_s(pt, tmp0, py##_1_0_0, upy##_0_0_1);  \
    upy##_0_1_0 = fmadd_s(pt, tmp0, py##_1_1_1, upy##_0_1_0);   \
    upy##_0_1_1 = fnmadd_s(pt, tmp0, py##_1_1_0, upy##_0_1_1);  \
    upy##_1_0_0 = fmadd_s(pt, tmp1, py##_1_0_1, upy##_1_0_0);   \
    upy##_1_0_1 = fnmadd_s(pt, tmp1, py##_1_0_0, upy##_1_0_1);  \
    upy##_1_1_0 = fmadd_s(pt, tmp1, py##_1_1_1, upy##_1_1_0);   \
    upy##_1_1_1 = fnmadd_s(pt, tmp1, py##_1_1_0, upy##_1_1_1);  \
    upy##_2_0_0 = fmadd_s(pt, tmp2, py##_1_0_1, upy##_2_0_0);   \
    upy##_2_0_1 = fnmadd_s(pt, tmp2, py##_1_0_0, upy##_2_0_1);  \
    upy##_2_1_0 = fmadd_s(pt, tmp2, py##_1_1_1, upy##_2_1_0);   \
    upy##_2_1_1 = fnmadd_s(pt, tmp2, py##_1_1_0, upy##_2_1_1);  \
                                                                \
    tmp0 = __load_backward_##pos((u_b), (u_f), 0, 2, 1);        \
    tmp1 = __load_backward_##pos((u_b), (u_f), 1, 2, 1);        \
    tmp2 = __load_backward_##pos((u_b), (u_f), 2, 2, 1);        \
    upy##_0_0_0 = fmadd_s(pt, tmp0, py##_2_0_1, upy##_0_0_0);   \
    upy##_0_0_1 = fnmadd_s(pt, tmp0, py##_2_0_0, upy##_0_0_1);  \
    upy##_0_1_0 = fmadd_s(pt, tmp0, py##_2_1_1, upy##_0_1_0);   \
    upy##_0_1_1 = fnmadd_s(pt, tmp0, py##_2_1_0, upy##_0_1_1);  \
    upy##_1_0_0 = fmadd_s(pt, tmp1, py##_2_0_1, upy##_1_0_0);   \
    upy##_1_0_1 = fnmadd_s(pt, tmp1, py##_2_0_0, upy##_1_0_1);  \
    upy##_1_1_0 = fmadd_s(pt, tmp1, py##_2_1_1, upy##_1_1_0);   \
    upy##_1_1_1 = fnmadd_s(pt, tmp1, py##_2_1_0, upy##_1_1_1);  \
    upy##_2_0_0 = fmadd_s(pt, tmp2, py##_2_0_1, upy##_2_0_0);   \
    upy##_2_0_1 = fnmadd_s(pt, tmp2, py##_2_0_0, upy##_2_0_1);  \
    upy##_2_1_0 = fmadd_s(pt, tmp2, py##_2_1_1, upy##_2_1_0);   \
    upy##_2_1_1 = fnmadd_s(pt, tmp2, py##_2_1_0, upy##_2_1_1);  \
  }

//
// link multiplication on two-spinor on a site with the last simd-site gauge
// upy1 = udag(VLEN-1) * y1
// upy2 = udag(VLEN-1) * y2
//
#define __mult_udag_y_3_(upy,py,u) {                                    \
    for (int c = 0; c < 3; ++c) {                                       \
      for (int s = 0; s < 2; ++s) {                                     \
        (upy).c[c][s][0]  = (u).c[c][0][0][VLENS-1] * (py).c[0][s][0];  \
        (upy).c[c][s][1]  = (u).c[c][0][0][VLENS-1] * (py).c[0][s][1];  \
        (upy).c[c][s][0] += (u).c[c][1][0][VLENS-1] * (py).c[1][s][0];  \
        (upy).c[c][s][1] += (u).c[c][1][0][VLENS-1] * (py).c[1][s][1];  \
        (upy).c[c][s][0] += (u).c[c][2][0][VLENS-1] * (py).c[2][s][0];  \
        (upy).c[c][s][1] += (u).c[c][2][0][VLENS-1] * (py).c[2][s][1];  \
        (upy).c[c][s][0] += (u).c[c][0][1][VLENS-1] * (py).c[0][s][1];  \
        (upy).c[c][s][1] -= (u).c[c][0][1][VLENS-1] * (py).c[0][s][0];  \
        (upy).c[c][s][0] += (u).c[c][1][1][VLENS-1] * (py).c[1][s][1];  \
        (upy).c[c][s][1] -= (u).c[c][1][1][VLENS-1] * (py).c[1][s][0];  \
        (upy).c[c][s][0] += (u).c[c][2][1][VLENS-1] * (py).c[2][s][1];  \
        (upy).c[c][s][1] -= (u).c[c][2][1][VLENS-1] * (py).c[2][s][0];  \
      }                                                                 \
    }                                                                   \
  }

//
// X-forward spin pre projection
// py1 = y1 - i y4
// py2 = y2 - i y3
//
#define __mult_x_forw_pre_(py,y) {                                      \
    for (int c = 0; c < 3; ++c) {                                       \
      for (int j = 0; j < VLENS; ++j) {                                 \
        (py).c[c][0][0][j] = (y).c[c][0][0][j] + (y).c[c][3][1][j];     \
        (py).c[c][0][1][j] = (y).c[c][0][1][j] - (y).c[c][3][0][j];     \
        (py).c[c][1][0][j] = (y).c[c][1][0][j] + (y).c[c][2][1][j];     \
        (py).c[c][1][1][j] = (y).c[c][1][1][j] - (y).c[c][2][0][j];     \
      }                                                                 \
    }                                                                   \
  }

//
// X-forward spin pre projection with forward simd-site-shifting
// py1(j) = y1(j+1) - i y4(j+1)
// py2(j) = y2(j+1) - i y3(j+1)
//
#define __mult_x_forw_pre_2_(py,y) {\
  for (int c = 0; c < 3; ++c) {\
    for (int j = 0; j < VLENS; ++j) {\
      (py).c[c][0][0][j] = (y)[c][0][0][j+1] + (y)[c][3][1][j+1];\
      (py).c[c][0][1][j] = (y)[c][0][1][j+1] - (y)[c][3][0][j+1];\
      (py).c[c][1][0][j] = (y)[c][1][0][j+1] + (y)[c][2][1][j+1];\
      (py).c[c][1][1][j] = (y)[c][1][1][j+1] - (y)[c][2][0][j+1];\
    }\
  }\
}

#define S__mult_x_forw_pre_2_(py, y, C)                         \
  py##_##C##_0_0 = fadd_s(pt, y##_##C##_0_0, y##_##C##_3_1);    \
  py##_##C##_0_1 = fsub_s(pt, y##_##C##_0_1, y##_##C##_3_0);    \
  py##_##C##_1_0 = fadd_s(pt, y##_##C##_1_0, y##_##C##_2_1);    \
  py##_##C##_1_1 = fsub_s(pt, y##_##C##_1_1, y##_##C##_2_0);
#define __mult_x_forw_pre_2_vec_(py,y) LOOP_3(S__mult_x_forw_pre_2_, py, y)

//
// X-forward spin pre projection with first simd-site spionr
// py1 = y1(0) - i y4(0)
// py2 = y2(0) - i y3(0)
//
#define __mult_x_forw_pre_3_(py,y) {                            \
    for (int c = 0; c < 3; ++c) {                               \
      (py).c[c][0][0] = (y).c[c][0][0][0] + (y).c[c][3][1][0];  \
      (py).c[c][0][1] = (y).c[c][0][1][0] - (y).c[c][3][0][0];  \
      (py).c[c][1][0] = (y).c[c][1][0][0] + (y).c[c][2][1][0];  \
      (py).c[c][1][1] = (y).c[c][1][1][0] - (y).c[c][2][0][0];  \
    }                                                           \
  }

//
// X-forward spin post reconstruction
// my1 = my1 +   upy1
// my2 = my2 +   upy2
// my3 = my3 + i upy2
// my4 = my4 + i upy1
//
#define __mult_x_forw_pst_(my,upy) {\
  _Pragma("loop norecurrence") \
  for (int c = 0; c < 3; ++c) {\
    for (int j = 0; j < VLENS; ++j) {\
      (my).c[c][0][0][j] = (upy).c[c][0][0][j];\
      (my).c[c][0][1][j] = (upy).c[c][0][1][j];\
      (my).c[c][1][0][j] = (upy).c[c][1][0][j];\
      (my).c[c][1][1][j] = (upy).c[c][1][1][j];\
      (my).c[c][2][0][j] = -(upy).c[c][1][1][j];\
      (my).c[c][2][1][j] = (upy).c[c][1][0][j];\
      (my).c[c][3][0][j] = -(upy).c[c][0][1][j];\
      (my).c[c][3][1][j] = (upy).c[c][0][0][j];\
    }\
  }\
}

#define S__mult_x_forw_pst_(my, upy, C)                                 \
  fstore1_s(pt, upy##_##C##_0_0        , (my).c, dims_scs, C, 0, 0);    \
  fstore1_s(pt, upy##_##C##_0_1        , (my).c, dims_scs, C, 0, 1);    \
  fstore1_s(pt, upy##_##C##_1_0        , (my).c, dims_scs, C, 1, 0);    \
  fstore1_s(pt, upy##_##C##_1_1        , (my).c, dims_scs, C, 1, 1);    \
  fstore1_s(pt, fneg_s(pt, upy##_##C##_1_1), (my).c, dims_scs, C, 2, 0); \
  fstore1_s(pt, upy##_##C##_1_0        , (my).c, dims_scs, C, 2, 1);    \
  fstore1_s(pt, fneg_s(pt, upy##_##C##_0_1), (my).c, dims_scs, C, 3, 0); \
  fstore1_s(pt, upy##_##C##_0_0        , (my).c, dims_scs, C, 3, 1);  
#define __mult_x_forw_pst_vec_(my,upy) LOOP_3(S__mult_x_forw_pst_, my, upy)

//
// X-forward spin post reconstruction at simd-site-shift last site
// my1(VLEN-1) = my1(VLEN-1) +   upy1
// my2(VLEN-1) = my2(VLEN-1) +   upy2
// my3(VLEN-1) = my3(VLEN-1) + i upy2
// my4(VLEN-1) = my4(VLEN-1) + i upy1
//
#define __mult_x_forw_pst_3_(my,upy) {                  \
    for (int c = 0; c < 3; ++c) {                       \
      (my).c[c][0][0][VLENS-1] += (upy).c[c][0][0];     \
      (my).c[c][0][1][VLENS-1] += (upy).c[c][0][1];     \
      (my).c[c][1][0][VLENS-1] += (upy).c[c][1][0];     \
      (my).c[c][1][1][VLENS-1] += (upy).c[c][1][1];     \
      (my).c[c][2][0][VLENS-1] -= (upy).c[c][1][1];     \
      (my).c[c][2][1][VLENS-1] += (upy).c[c][1][0];     \
      (my).c[c][3][0][VLENS-1] -= (upy).c[c][0][1];     \
      (my).c[c][3][1][VLENS-1] += (upy).c[c][0][0];     \
    }                                                   \
  }


//
// X-backward spin pre projection
// py1 = y1 + i y4
// py2 = y2 + i y3
//
#define __mult_x_back_pre_(py,y) {                                      \
    for (int c = 0; c < 3; ++c) {                                       \
      for (int j = 0; j < VLENS; ++j) {                                 \
        (py).c[c][0][0][j] = (y).c[c][0][0][j] - (y).c[c][3][1][j];     \
        (py).c[c][0][1][j] = (y).c[c][0][1][j] + (y).c[c][3][0][j];     \
        (py).c[c][1][0][j] = (y).c[c][1][0][j] - (y).c[c][2][1][j];     \
        (py).c[c][1][1][j] = (y).c[c][1][1][j] + (y).c[c][2][0][j];     \
      }                                                                 \
    }                                                                   \
  }

//
// X-backward spin pre projection with backward simd-site-shifti(ng
// py1(j) = y1(j+VLEN-1) + i y4(j+VLEN-1)
// py2(j) = y2(j+VLEN-1) + i y3(j+VLEN-1)
//
#define __mult_x_back_pre_2_(py,y) {\
  for (int c = 0; c < 3; ++c) {\
    for (int j = 0; j < VLENS; ++j) {\
      (py).c[c][0][0][j] = (y)[c][0][0][j+VLENS-1] - (y)[c][3][1][j+VLENS-1];\
      (py).c[c][0][1][j] = (y)[c][0][1][j+VLENS-1] + (y)[c][3][0][j+VLENS-1];\
      (py).c[c][1][0][j] = (y)[c][1][0][j+VLENS-1] - (y)[c][2][1][j+VLENS-1];\
      (py).c[c][1][1][j] = (y)[c][1][1][j+VLENS-1] + (y)[c][2][0][j+VLENS-1];\
    }\
  }\
}

#define S__mult_x_back_pre_2_(py, y, C)                         \
  py##_##C##_0_0 = fsub_s(pt, y##_##C##_0_0, y##_##C##_3_1);    \
  py##_##C##_0_1 = fadd_s(pt, y##_##C##_0_1, y##_##C##_3_0);    \
  py##_##C##_1_0 = fsub_s(pt, y##_##C##_1_0, y##_##C##_2_1);    \
  py##_##C##_1_1 = fadd_s(pt, y##_##C##_1_1, y##_##C##_2_0);
#define __mult_x_back_pre_2_vec_(py,y) LOOP_3(S__mult_x_back_pre_2_, py, y)

//
// X-backward spin pre projection with last simd-site spinor
// py1 = y1(VLEN-1) + i y4(VLEN-1)
// py2 = y2(VLEN-1) + i y3(VLEN-1)
//
#define __mult_x_back_pre_3_(py,y) {                                    \
    for (int c = 0; c < 3; ++c) {                                       \
      (py).c[c][0][0] = (y).c[c][0][0][VLENS-1] - (y).c[c][3][1][VLENS-1]; \
      (py).c[c][0][1] = (y).c[c][0][1][VLENS-1] + (y).c[c][3][0][VLENS-1]; \
      (py).c[c][1][0] = (y).c[c][1][0][VLENS-1] - (y).c[c][2][1][VLENS-1]; \
      (py).c[c][1][1] = (y).c[c][1][1][VLENS-1] + (y).c[c][2][0][VLENS-1]; \
    }                                                                   \
  }

//
// X-backward spin post reconstruction
// my1 = my1 +   upy1
// my2 = my2 +   upy2
// my3 = my3 - i upy2
// my4 = my4 - i upy1
//
#define __mult_x_back_pst_(my,upy) {\
  _Pragma("loop norecurrence") \
  for (int c = 0; c < 3; ++c) {\
    for (int j = 0; j < VLENS; ++j) {\
      (my).c[c][0][0][j] += (upy).c[c][0][0][j];\
      (my).c[c][0][1][j] += (upy).c[c][0][1][j];\
      (my).c[c][1][0][j] += (upy).c[c][1][0][j];\
      (my).c[c][1][1][j] += (upy).c[c][1][1][j];\
      (my).c[c][2][0][j] += (upy).c[c][1][1][j];\
      (my).c[c][2][1][j] -= (upy).c[c][1][0][j];\
      (my).c[c][3][0][j] += (upy).c[c][0][1][j];\
      (my).c[c][3][1][j] -= (upy).c[c][0][0][j];\
    }\
  }\
}

#define S__mult_x_back_pst_(my, upy, C)                                 \
  fstore1_s(pt, fadd_s(pt, fload1_s(pt, (my).c, dims_scs, C, 0, 0), upy##_##C##_0_0), (my).c, dims_scs, C, 0, 0); \
  fstore1_s(pt, fadd_s(pt, fload1_s(pt, (my).c, dims_scs, C, 0, 1), upy##_##C##_0_1), (my).c, dims_scs, C, 0, 1); \
  fstore1_s(pt, fadd_s(pt, fload1_s(pt, (my).c, dims_scs, C, 1, 0), upy##_##C##_1_0), (my).c, dims_scs, C, 1, 0); \
  fstore1_s(pt, fadd_s(pt, fload1_s(pt, (my).c, dims_scs, C, 1, 1), upy##_##C##_1_1), (my).c, dims_scs, C, 1, 1); \
  fstore1_s(pt, fadd_s(pt, fload1_s(pt, (my).c, dims_scs, C, 2, 0), upy##_##C##_1_1), (my).c, dims_scs, C, 2, 0); \
  fstore1_s(pt, fsub_s(pt, fload1_s(pt, (my).c, dims_scs, C, 2, 1), upy##_##C##_1_0), (my).c, dims_scs, C, 2, 1); \
  fstore1_s(pt, fadd_s(pt, fload1_s(pt, (my).c, dims_scs, C, 3, 0), upy##_##C##_0_1), (my).c, dims_scs, C, 3, 0); \
  fstore1_s(pt, fsub_s(pt, fload1_s(pt, (my).c, dims_scs, C, 3, 1), upy##_##C##_0_0), (my).c, dims_scs, C, 3, 1);  
#define __mult_x_back_pst_vec_(my,upy) LOOP_3(S__mult_x_back_pst_, my, upy)

//
// X-backward spin post reconstruction at simd-site-sift first site
// my1(0) = my1(0) +   upy1
// my2(0) = my2(0) +   upy2
// my3(0) = my3(0) - i upy2
// my4(0) = my4(0) - i upy1
//
#define __mult_x_back_pst_3_(my,upy) {          \
    for (int c = 0; c < 3; ++c) {               \
      (my).c[c][0][0][0] += (upy).c[c][0][0];   \
      (my).c[c][0][1][0] += (upy).c[c][0][1];   \
      (my).c[c][1][0][0] += (upy).c[c][1][0];   \
      (my).c[c][1][1][0] += (upy).c[c][1][1];   \
      (my).c[c][2][0][0] += (upy).c[c][1][1];   \
      (my).c[c][2][1][0] -= (upy).c[c][1][0];   \
      (my).c[c][3][0][0] += (upy).c[c][0][1];   \
      (my).c[c][3][1][0] -= (upy).c[c][0][0];   \
    }                                           \
  }

#define S__mult_x_back_pst_3_vec_(my, upy, C)   \
      my##_##C##_0_0 = fadd_s_m(pt_lowest, my##_##C##_0_0, fdup_s((upy).c[C][0][0]));   \
      my##_##C##_0_1 = fadd_s_m(pt_lowest, my##_##C##_0_1, fdup_s((upy).c[C][0][1]));   \
      my##_##C##_1_0 = fadd_s_m(pt_lowest, my##_##C##_1_0, fdup_s((upy).c[C][1][0]));   \
      my##_##C##_1_1 = fadd_s_m(pt_lowest, my##_##C##_1_1, fdup_s((upy).c[C][1][1]));   \
      my##_##C##_2_0 = fadd_s_m(pt_lowest, my##_##C##_2_0, fdup_s((upy).c[C][1][1]));   \
      my##_##C##_2_1 = fsub_s_m(pt_lowest, my##_##C##_2_1, fdup_s((upy).c[C][1][0]));   \
      my##_##C##_3_0 = fadd_s_m(pt_lowest, my##_##C##_3_0, fdup_s((upy).c[C][0][1]));   \
      my##_##C##_3_1 = fsub_s_m(pt_lowest, my##_##C##_3_1, fdup_s((upy).c[C][0][0]));
#define __mult_x_back_pst_3_vec_(my,upy) LOOP_3(S__mult_x_back_pst_3_vec_, my, upy)

//
// Y-forward spin pre projection
// py1 = y1 - y4
// py2 = y2 + y3
//
#define __mult_y_forw_pre_(py,y) {\
  for (int c = 0; c < 3; ++c) {\
    for (int j = 0; j < VLENS; ++j) {\
      (py).c[c][0][0][j] = (y).c[c][0][0][j] - (y).c[c][3][0][j];\
      (py).c[c][0][1][j] = (y).c[c][0][1][j] - (y).c[c][3][1][j];\
      (py).c[c][1][0][j] = (y).c[c][1][0][j] + (y).c[c][2][0][j];\
      (py).c[c][1][1][j] = (y).c[c][1][1][j] + (y).c[c][2][1][j];\
    }\
  }\
}


#define S__mult_y_forw_pre_(py, y, C)                                   \
  py##_##C##_0_0 = fsub_s(pt, fload1_s(pt, (y).c, dims_scs, C, 0, 0), fload1_s(pt, (y).c, dims_scs, C, 3, 0)); \
  py##_##C##_0_1 = fsub_s(pt, fload1_s(pt, (y).c, dims_scs, C, 0, 1), fload1_s(pt, (y).c, dims_scs, C, 3, 1)); \
  py##_##C##_1_0 = fadd_s(pt, fload1_s(pt, (y).c, dims_scs, C, 1, 0), fload1_s(pt, (y).c, dims_scs, C, 2, 0)); \
  py##_##C##_1_1 = fadd_s(pt, fload1_s(pt, (y).c, dims_scs, C, 1, 1), fload1_s(pt, (y).c, dims_scs, C, 2, 1));
#define __mult_y_forw_pre_vec_(py,y) LOOP_3(S__mult_y_forw_pre_, py, y)

//
// Y-forward spin post reconstruction
// my1 = my1 +   upy1
// my2 = my2 +   upy2
// my3 = my3 +   upy2
// my4 = my4 -   upy1
//
#define __mult_y_forw_pst_(my,upy) {\
  _Pragma("loop norecurrence") \
  for (int c = 0; c < 3; ++c) {\
    for (int j = 0; j < VLENS; ++j) {\
      (my).c[c][0][0][j] += (upy).c[c][0][0][j];\
      (my).c[c][0][1][j] += (upy).c[c][0][1][j];\
      (my).c[c][1][0][j] += (upy).c[c][1][0][j];\
      (my).c[c][1][1][j] += (upy).c[c][1][1][j];\
      (my).c[c][2][0][j] += (upy).c[c][1][0][j];\
      (my).c[c][2][1][j] += (upy).c[c][1][1][j];\
      (my).c[c][3][0][j] -= (upy).c[c][0][0][j];\
      (my).c[c][3][1][j] -= (upy).c[c][0][1][j];\
    }\
  }\
}

#define S__mult_y_forw_pst_(my, upy, C)                                 \
  fstore1_s(pt, fadd_s(pt, fload1_s(pt, (my).c, dims_scs, C, 0, 0), upy##_##C##_0_0), (my).c, dims_scs, C, 0, 0); \
  fstore1_s(pt, fadd_s(pt, fload1_s(pt, (my).c, dims_scs, C, 0, 1), upy##_##C##_0_1), (my).c, dims_scs, C, 0, 1); \
  fstore1_s(pt, fadd_s(pt, fload1_s(pt, (my).c, dims_scs, C, 1, 0), upy##_##C##_1_0), (my).c, dims_scs, C, 1, 0); \
  fstore1_s(pt, fadd_s(pt, fload1_s(pt, (my).c, dims_scs, C, 1, 1), upy##_##C##_1_1), (my).c, dims_scs, C, 1, 1); \
  fstore1_s(pt, fadd_s(pt, fload1_s(pt, (my).c, dims_scs, C, 2, 0), upy##_##C##_1_0), (my).c, dims_scs, C, 2, 0); \
  fstore1_s(pt, fadd_s(pt, fload1_s(pt, (my).c, dims_scs, C, 2, 1), upy##_##C##_1_1), (my).c, dims_scs, C, 2, 1); \
  fstore1_s(pt, fsub_s(pt, fload1_s(pt, (my).c, dims_scs, C, 3, 0), upy##_##C##_0_0), (my).c, dims_scs, C, 3, 0); \
  fstore1_s(pt, fsub_s(pt, fload1_s(pt, (my).c, dims_scs, C, 3, 1), upy##_##C##_0_1), (my).c, dims_scs, C, 3, 1);  
#define __mult_y_forw_pst_vec_(my,upy) LOOP_3(S__mult_y_forw_pst_, my, upy)

//
// Y-backward spin pre projection
// py1 = y1 + y4
// py2 = y2 - y3
//
#define __mult_y_back_pre_(py,y) {\
  for (int c = 0; c < 3; ++c) {\
    for (int j = 0; j < VLENS; ++j) {\
      (py).c[c][0][0][j] = (y).c[c][0][0][j] + (y).c[c][3][0][j];\
      (py).c[c][0][1][j] = (y).c[c][0][1][j] + (y).c[c][3][1][j];\
      (py).c[c][1][0][j] = (y).c[c][1][0][j] - (y).c[c][2][0][j];\
      (py).c[c][1][1][j] = (y).c[c][1][1][j] - (y).c[c][2][1][j];\
    }\
  }\
}

#define S__mult_y_back_pre_(py, y, C)                                   \
  py##_##C##_0_0 = fadd_s(pt, fload1_s(pt, (y).c, dims_scs, C, 0, 0), fload1_s(pt, (y).c, dims_scs, C, 3, 0)); \
  py##_##C##_0_1 = fadd_s(pt, fload1_s(pt, (y).c, dims_scs, C, 0, 1), fload1_s(pt, (y).c, dims_scs, C, 3, 1)); \
  py##_##C##_1_0 = fsub_s(pt, fload1_s(pt, (y).c, dims_scs, C, 1, 0), fload1_s(pt, (y).c, dims_scs, C, 2, 0)); \
  py##_##C##_1_1 = fsub_s(pt, fload1_s(pt, (y).c, dims_scs, C, 1, 1), fload1_s(pt, (y).c, dims_scs, C, 2, 1));
#define __mult_y_back_pre_vec_(py,y) LOOP_3(S__mult_y_back_pre_, py, y)

//
// Y-backward spin post reconstruction
// my1 = my1 +   upy1
// my2 = my2 +   upy2
// my3 = my3 -   upy2
// my4 = my4 +   upy1
//
#define __mult_y_back_pst_(my,upy) {\
  _Pragma("loop norecurrence") \
  for (int c = 0; c < 3; ++c) {\
    for (int j = 0; j < VLENS; ++j) {\
      (my).c[c][0][0][j] += (upy).c[c][0][0][j];\
      (my).c[c][0][1][j] += (upy).c[c][0][1][j];\
      (my).c[c][1][0][j] += (upy).c[c][1][0][j];\
      (my).c[c][1][1][j] += (upy).c[c][1][1][j];\
      (my).c[c][2][0][j] -= (upy).c[c][1][0][j];\
      (my).c[c][2][1][j] -= (upy).c[c][1][1][j];\
      (my).c[c][3][0][j] += (upy).c[c][0][0][j];\
      (my).c[c][3][1][j] += (upy).c[c][0][1][j];\
    }\
  }\
}

#define S__mult_y_back_pst_(my, upy, C)                                 \
  fstore1_s(pt, fadd_s(pt, fload1_s(pt, (my).c, dims_scs, C, 0, 0), upy##_##C##_0_0), (my).c, dims_scs, C, 0, 0); \
  fstore1_s(pt, fadd_s(pt, fload1_s(pt, (my).c, dims_scs, C, 0, 1), upy##_##C##_0_1), (my).c, dims_scs, C, 0, 1); \
  fstore1_s(pt, fadd_s(pt, fload1_s(pt, (my).c, dims_scs, C, 1, 0), upy##_##C##_1_0), (my).c, dims_scs, C, 1, 0); \
  fstore1_s(pt, fadd_s(pt, fload1_s(pt, (my).c, dims_scs, C, 1, 1), upy##_##C##_1_1), (my).c, dims_scs, C, 1, 1); \
  fstore1_s(pt, fsub_s(pt, fload1_s(pt, (my).c, dims_scs, C, 2, 0), upy##_##C##_1_0), (my).c, dims_scs, C, 2, 0); \
  fstore1_s(pt, fsub_s(pt, fload1_s(pt, (my).c, dims_scs, C, 2, 1), upy##_##C##_1_1), (my).c, dims_scs, C, 2, 1); \
  fstore1_s(pt, fadd_s(pt, fload1_s(pt, (my).c, dims_scs, C, 3, 0), upy##_##C##_0_0), (my).c, dims_scs, C, 3, 0); \
  fstore1_s(pt, fadd_s(pt, fload1_s(pt, (my).c, dims_scs, C, 3, 1), upy##_##C##_0_1), (my).c, dims_scs, C, 3, 1);  
#define __mult_y_back_pst_vec_(my,upy) LOOP_3(S__mult_y_back_pst_, my, upy)

#define S__mult_y_back_pst_reg_(my, upy, C)                                 \
  my##_##C##_0_0 = fadd_s(pt, my##_##C##_0_0, fload1_s(pt, (upy).c, dims_projscs, C, 0, 0)); \
  my##_##C##_0_1 = fadd_s(pt, my##_##C##_0_1, fload1_s(pt, (upy).c, dims_projscs, C, 0, 1)); \
  my##_##C##_1_0 = fadd_s(pt, my##_##C##_1_0, fload1_s(pt, (upy).c, dims_projscs, C, 1, 0)); \
  my##_##C##_1_1 = fadd_s(pt, my##_##C##_1_1, fload1_s(pt, (upy).c, dims_projscs, C, 1, 1)); \
  my##_##C##_2_0 = fsub_s(pt, my##_##C##_2_0, fload1_s(pt, (upy).c, dims_projscs, C, 1, 0)); \
  my##_##C##_2_1 = fsub_s(pt, my##_##C##_2_1, fload1_s(pt, (upy).c, dims_projscs, C, 1, 1)); \
  my##_##C##_3_0 = fadd_s(pt, my##_##C##_3_0, fload1_s(pt, (upy).c, dims_projscs, C, 0, 0)); \
  my##_##C##_3_1 = fadd_s(pt, my##_##C##_3_1, fload1_s(pt, (upy).c, dims_projscs, C, 0, 1));  
#define __mult_y_back_pst_vec_reg_(my,upy) LOOP_3(S__mult_y_back_pst_reg_, my, upy)

//
// Z-forward spin pre projection
// py1 = y1 - i y3
// py2 = y2 + i y4
//
#define __mult_z_forw_pre_(py,y) {\
  for (int c = 0; c < 3; ++c) {\
    for (int j = 0; j < VLENS; ++j) {\
      (py).c[c][0][0][j] = (y).c[c][0][0][j] + (y).c[c][2][1][j];\
      (py).c[c][0][1][j] = (y).c[c][0][1][j] - (y).c[c][2][0][j];\
      (py).c[c][1][0][j] = (y).c[c][1][0][j] - (y).c[c][3][1][j];\
      (py).c[c][1][1][j] = (y).c[c][1][1][j] + (y).c[c][3][0][j];\
    }\
  }\
}

#define S__mult_z_forw_pre_(py, y, C)                                   \
  py##_##C##_0_0 = fadd_s(pt, fload1_s(pt, (y).c, dims_scs, C, 0, 0), fload1_s(pt, (y).c, dims_scs, C, 2, 1)); \
  py##_##C##_0_1 = fsub_s(pt, fload1_s(pt, (y).c, dims_scs, C, 0, 1), fload1_s(pt, (y).c, dims_scs, C, 2, 0)); \
  py##_##C##_1_0 = fsub_s(pt, fload1_s(pt, (y).c, dims_scs, C, 1, 0), fload1_s(pt, (y).c, dims_scs, C, 3, 1)); \
  py##_##C##_1_1 = fadd_s(pt, fload1_s(pt, (y).c, dims_scs, C, 1, 1), fload1_s(pt, (y).c, dims_scs, C, 3, 0));
#define __mult_z_forw_pre_vec_(py,y) LOOP_3(S__mult_z_forw_pre_, py, y)

//
// Z-forward spin post reconstruction
// my1 = my1 +   upy1
// my2 = my2 +   upy2
// my3 = my3 + i upy1
// my4 = my4 - i upy2
//
#define __mult_z_forw_pst_(my,upy) {\
  _Pragma("loop norecurrence") \
  for (int c = 0; c < 3; ++c) {\
    for (int j = 0; j < VLENS; ++j) {\
      (my).c[c][0][0][j] += (upy).c[c][0][0][j];\
      (my).c[c][0][1][j] += (upy).c[c][0][1][j];\
      (my).c[c][1][0][j] += (upy).c[c][1][0][j];\
      (my).c[c][1][1][j] += (upy).c[c][1][1][j];\
      (my).c[c][2][0][j] -= (upy).c[c][0][1][j];\
      (my).c[c][2][1][j] += (upy).c[c][0][0][j];\
      (my).c[c][3][0][j] += (upy).c[c][1][1][j];\
      (my).c[c][3][1][j] -= (upy).c[c][1][0][j];\
    }\
  }\
}

#define S__mult_z_forw_pst_(my, upy, C)                                 \
  fstore1_s(pt, fadd_s(pt, fload1_s(pt, (my).c, dims_scs, C, 0, 0), upy##_##C##_0_0), (my).c, dims_scs, C, 0, 0); \
  fstore1_s(pt, fadd_s(pt, fload1_s(pt, (my).c, dims_scs, C, 0, 1), upy##_##C##_0_1), (my).c, dims_scs, C, 0, 1); \
  fstore1_s(pt, fadd_s(pt, fload1_s(pt, (my).c, dims_scs, C, 1, 0), upy##_##C##_1_0), (my).c, dims_scs, C, 1, 0); \
  fstore1_s(pt, fadd_s(pt, fload1_s(pt, (my).c, dims_scs, C, 1, 1), upy##_##C##_1_1), (my).c, dims_scs, C, 1, 1); \
  fstore1_s(pt, fsub_s(pt, fload1_s(pt, (my).c, dims_scs, C, 2, 0), upy##_##C##_0_1), (my).c, dims_scs, C, 2, 0); \
  fstore1_s(pt, fadd_s(pt, fload1_s(pt, (my).c, dims_scs, C, 2, 1), upy##_##C##_0_0), (my).c, dims_scs, C, 2, 1); \
  fstore1_s(pt, fadd_s(pt, fload1_s(pt, (my).c, dims_scs, C, 3, 0), upy##_##C##_1_1), (my).c, dims_scs, C, 3, 0); \
  fstore1_s(pt, fsub_s(pt, fload1_s(pt, (my).c, dims_scs, C, 3, 1), upy##_##C##_1_0), (my).c, dims_scs, C, 3, 1);  
#define __mult_z_forw_pst_vec_(my,upy) LOOP_3(S__mult_z_forw_pst_, my, upy)

//
// Z-backward spin pre projection
// py1 = y1 + i y3
// py2 = y2 - i y4
//
#define __mult_z_back_pre_(py,y) {\
  for (int c = 0; c < 3; ++c) {\
    for (int j = 0; j < VLENS; ++j) {\
      (py).c[c][0][0][j] = (y).c[c][0][0][j] - (y).c[c][2][1][j];\
      (py).c[c][0][1][j] = (y).c[c][0][1][j] + (y).c[c][2][0][j];\
      (py).c[c][1][0][j] = (y).c[c][1][0][j] + (y).c[c][3][1][j];\
      (py).c[c][1][1][j] = (y).c[c][1][1][j] - (y).c[c][3][0][j];\
    }\
  }\
}

#define S__mult_z_back_pre_(py, y, C)                                   \
  py##_##C##_0_0 = fsub_s(pt, fload1_s(pt, (y).c, dims_scs, C, 0, 0), fload1_s(pt, (y).c, dims_scs, C, 2, 1)); \
  py##_##C##_0_1 = fadd_s(pt, fload1_s(pt, (y).c, dims_scs, C, 0, 1), fload1_s(pt, (y).c, dims_scs, C, 2, 0)); \
  py##_##C##_1_0 = fadd_s(pt, fload1_s(pt, (y).c, dims_scs, C, 1, 0), fload1_s(pt, (y).c, dims_scs, C, 3, 1)); \
  py##_##C##_1_1 = fsub_s(pt, fload1_s(pt, (y).c, dims_scs, C, 1, 1), fload1_s(pt, (y).c, dims_scs, C, 3, 0));
#define __mult_z_back_pre_vec_(py,y) LOOP_3(S__mult_z_back_pre_, py, y)

//
// Z-backward spin post reconstruction
// my1 = my1 +   upy1
// my2 = my2 +   upy2
// my3 = my3 - i upy1
// my4 = my4 + i upy2
//
#define __mult_z_back_pst_(my,upy) {\
  _Pragma("loop norecurrence") \
  for (int c = 0; c < 3; ++c) {\
    for (int j = 0; j < VLENS; ++j) {\
      (my).c[c][0][0][j] += (upy).c[c][0][0][j];\
      (my).c[c][0][1][j] += (upy).c[c][0][1][j];\
      (my).c[c][1][0][j] += (upy).c[c][1][0][j];\
      (my).c[c][1][1][j] += (upy).c[c][1][1][j];\
      (my).c[c][2][0][j] += (upy).c[c][0][1][j];\
      (my).c[c][2][1][j] -= (upy).c[c][0][0][j];\
      (my).c[c][3][0][j] -= (upy).c[c][1][1][j];\
      (my).c[c][3][1][j] += (upy).c[c][1][0][j];\
    }\
  }\
}

#define S__mult_z_back_pst_(my, upy, C)                                 \
  fstore1_s(pt, fadd_s(pt, fload1_s(pt, (my).c, dims_scs, C, 0, 0), upy##_##C##_0_0), (my).c, dims_scs, C, 0, 0); \
  fstore1_s(pt, fadd_s(pt, fload1_s(pt, (my).c, dims_scs, C, 0, 1), upy##_##C##_0_1), (my).c, dims_scs, C, 0, 1); \
  fstore1_s(pt, fadd_s(pt, fload1_s(pt, (my).c, dims_scs, C, 1, 0), upy##_##C##_1_0), (my).c, dims_scs, C, 1, 0); \
  fstore1_s(pt, fadd_s(pt, fload1_s(pt, (my).c, dims_scs, C, 1, 1), upy##_##C##_1_1), (my).c, dims_scs, C, 1, 1); \
  fstore1_s(pt, fadd_s(pt, fload1_s(pt, (my).c, dims_scs, C, 2, 0), upy##_##C##_0_1), (my).c, dims_scs, C, 2, 0); \
  fstore1_s(pt, fsub_s(pt, fload1_s(pt, (my).c, dims_scs, C, 2, 1), upy##_##C##_0_0), (my).c, dims_scs, C, 2, 1); \
  fstore1_s(pt, fsub_s(pt, fload1_s(pt, (my).c, dims_scs, C, 3, 0), upy##_##C##_1_1), (my).c, dims_scs, C, 3, 0); \
  fstore1_s(pt, fadd_s(pt, fload1_s(pt, (my).c, dims_scs, C, 3, 1), upy##_##C##_1_0), (my).c, dims_scs, C, 3, 1);  
#define __mult_z_back_pst_vec_(my,upy) LOOP_3(S__mult_z_back_pst_, my, upy)

#define S__mult_z_back_pst_reg_(my, upy, C)                                 \
  my##_##C##_0_0 = fadd_s(pt, my##_##C##_0_0, fload1_s(pt, (upy).c, dims_projscs, C, 0, 0)); \
  my##_##C##_0_1 = fadd_s(pt, my##_##C##_0_1, fload1_s(pt, (upy).c, dims_projscs, C, 0, 1)); \
  my##_##C##_1_0 = fadd_s(pt, my##_##C##_1_0, fload1_s(pt, (upy).c, dims_projscs, C, 1, 0)); \
  my##_##C##_1_1 = fadd_s(pt, my##_##C##_1_1, fload1_s(pt, (upy).c, dims_projscs, C, 1, 1)); \
  my##_##C##_2_0 = fadd_s(pt, my##_##C##_2_0, fload1_s(pt, (upy).c, dims_projscs, C, 0, 1)); \
  my##_##C##_2_1 = fsub_s(pt, my##_##C##_2_1, fload1_s(pt, (upy).c, dims_projscs, C, 0, 0)); \
  my##_##C##_3_0 = fsub_s(pt, my##_##C##_3_0, fload1_s(pt, (upy).c, dims_projscs, C, 1, 1)); \
  my##_##C##_3_1 = fadd_s(pt, my##_##C##_3_1, fload1_s(pt, (upy).c, dims_projscs, C, 1, 0));  
#define __mult_z_back_pst_vec_reg_(my,upy) LOOP_3(S__mult_z_back_pst_reg_, my, upy)

//
// T-forward spin pre projection
// py1 = 2*y3
// py2 = 2*y4
//
#define __mult_t_forw_pre_(py,y) {\
  for (int c = 0; c < 3; ++c) {\
    for (int j = 0; j < VLENS; ++j) {\
      (py).c[c][0][0][j] = 2*(y).c[c][2][0][j];\
      (py).c[c][0][1][j] = 2*(y).c[c][2][1][j];\
      (py).c[c][1][0][j] = 2*(y).c[c][3][0][j];\
      (py).c[c][1][1][j] = 2*(y).c[c][3][1][j];\
    }\
  }\
}

#define S__mult_t_forw_pre_(py, y, C)                                   \
  py##_##C##_0_0 = fadd_s(pt, fload1_s(pt, (y).c, dims_scs, C, 2, 0), fload1_s(pt, (y).c, dims_scs, C, 2, 0)); \
  py##_##C##_0_1 = fadd_s(pt, fload1_s(pt, (y).c, dims_scs, C, 2, 1), fload1_s(pt, (y).c, dims_scs, C, 2, 1)); \
  py##_##C##_1_0 = fadd_s(pt, fload1_s(pt, (y).c, dims_scs, C, 3, 0), fload1_s(pt, (y).c, dims_scs, C, 3, 0)); \
  py##_##C##_1_1 = fadd_s(pt, fload1_s(pt, (y).c, dims_scs, C, 3, 1), fload1_s(pt, (y).c, dims_scs, C, 3, 1));
#define __mult_t_forw_pre_vec_(py,y) LOOP_3(S__mult_t_forw_pre_, py, y)

//
// T-forward spin pre projection with boundary condition
// py1 = fbc *y3
// py2 = fbc *y4
//
// fbc = 2 * phase
//
#define __mult_t_forw_pre_bc_(py,y,fbc) {               \
    for (int c = 0; c < 3; ++c) {                       \
      for (int j = 0; j < VLENS; ++j) {                 \
        (py).c[c][0][0][j] = (fbc)*(y).c[c][2][0][j];   \
        (py).c[c][0][1][j] = (fbc)*(y).c[c][2][1][j];   \
        (py).c[c][1][0][j] = (fbc)*(y).c[c][3][0][j];   \
        (py).c[c][1][1][j] = (fbc)*(y).c[c][3][1][j];   \
      }                                                 \
    }                                                   \
  }

//
// Y-forward spin post reconstruction
// my1 = my1
// my2 = my2
// my3 = my3 +   upy1
// my4 = my4 +   upy2
//
#define __mult_t_forw_pst_(my,upy) {\
  _Pragma("loop norecurrence") \
  for (int c = 0; c < 3; ++c) {\
    for (int j = 0; j < VLENS; ++j) {\
      (my).c[c][2][0][j] += (upy).c[c][0][0][j];\
      (my).c[c][2][1][j] += (upy).c[c][0][1][j];\
      (my).c[c][3][0][j] += (upy).c[c][1][0][j];\
      (my).c[c][3][1][j] += (upy).c[c][1][1][j];\
    }\
  }\
}

#define S__mult_t_forw_pst_(my, upy, C)                                 \
  fstore1_s(pt, fadd_s(pt, fload1_s(pt, (my).c, dims_scs, C, 2, 0), upy##_##C##_0_0), (my).c, dims_scs, C, 2, 0); \
  fstore1_s(pt, fadd_s(pt, fload1_s(pt, (my).c, dims_scs, C, 2, 1), upy##_##C##_0_1), (my).c, dims_scs, C, 2, 1); \
  fstore1_s(pt, fadd_s(pt, fload1_s(pt, (my).c, dims_scs, C, 3, 0), upy##_##C##_1_0), (my).c, dims_scs, C, 3, 0); \
  fstore1_s(pt, fadd_s(pt, fload1_s(pt, (my).c, dims_scs, C, 3, 1), upy##_##C##_1_1), (my).c, dims_scs, C, 3, 1);  
#define __mult_t_forw_pst_vec_(my,upy) LOOP_3(S__mult_t_forw_pst_, my, upy)

//
// Y-forward spin post reconstruction with bounday-condition
// my1 = my1
// my2 = my2
// my3 = my3 +   upy1 * fbc
// my4 = my4 +   upy2 * fbc
//
#define __mult_t_forw_pst_bc_(my,upy,fbc) {                     \
    for (int c = 0; c < 3; ++c) {                               \
      for (int j = 0; j < VLENS; ++j) {                         \
        (my).c[c][2][0][j] += (upy).c[c][0][0][j] * (fbc);      \
        (my).c[c][2][1][j] += (upy).c[c][0][1][j] * (fbc);      \
        (my).c[c][3][0][j] += (upy).c[c][1][0][j] * (fbc);      \
        (my).c[c][3][1][j] += (upy).c[c][1][1][j] * (fbc);      \
      }                                                         \
    }                                                           \
  }

#define S__mult_t_forw_pst_bc_(my, upy, fbc, C)                                 \
  fstore1_s(pt, fmadd_s(pt, fdup_s(fbc), upy##_##C##_0_0, fload1_s(pt, (my).c, dims_scs, C, 2, 0)), (my).c, dims_scs, C, 2, 0); \
  fstore1_s(pt, fmadd_s(pt, fdup_s(fbc), upy##_##C##_0_1, fload1_s(pt, (my).c, dims_scs, C, 2, 1)), (my).c, dims_scs, C, 2, 1); \
  fstore1_s(pt, fmadd_s(pt, fdup_s(fbc), upy##_##C##_1_0, fload1_s(pt, (my).c, dims_scs, C, 3, 0)), (my).c, dims_scs, C, 3, 0); \
  fstore1_s(pt, fmadd_s(pt, fdup_s(fbc), upy##_##C##_1_1, fload1_s(pt, (my).c, dims_scs, C, 3, 1)), (my).c, dims_scs, C, 3, 1);  
#define __mult_t_forw_pst_bc_vec_(my,upy,fbc) LOOP_3(S__mult_t_forw_pst_bc_, my, upy, fbc)

//
// T-backward spin pre projection
// py1 = 2*y1
// py2 = 2*y2
//
#define __mult_t_back_pre_(py,y) {\
  for (int c = 0; c < 3; ++c) {\
    for (int j = 0; j < VLENS; ++j) {\
      (py).c[c][0][0][j] = 2*(y).c[c][0][0][j];\
      (py).c[c][0][1][j] = 2*(y).c[c][0][1][j];\
      (py).c[c][1][0][j] = 2*(y).c[c][1][0][j];\
      (py).c[c][1][1][j] = 2*(y).c[c][1][1][j];\
    }\
  }\
}

#define S__mult_t_back_pre_(py, y, C)                                   \
  py##_##C##_0_0 = fadd_s(pt, fload1_s(pt, (y).c, dims_scs, C, 0, 0), fload1_s(pt, (y).c, dims_scs, C, 0, 0)); \
  py##_##C##_0_1 = fadd_s(pt, fload1_s(pt, (y).c, dims_scs, C, 0, 1), fload1_s(pt, (y).c, dims_scs, C, 0, 1)); \
  py##_##C##_1_0 = fadd_s(pt, fload1_s(pt, (y).c, dims_scs, C, 1, 0), fload1_s(pt, (y).c, dims_scs, C, 1, 0)); \
  py##_##C##_1_1 = fadd_s(pt, fload1_s(pt, (y).c, dims_scs, C, 1, 1), fload1_s(pt, (y).c, dims_scs, C, 1, 1));
#define __mult_t_back_pre_vec_(py,y) LOOP_3(S__mult_t_back_pre_, py, y)

//
// T-backward spin pre projection with boundary condition
// py1 = fbc*y1
// py2 = fbc*y2
//
// fbc = 2 * phase
//
#define __mult_t_back_pre_bc_(py,y,fbc) {               \
    for (int c = 0; c < 3; ++c) {                       \
      for (int j = 0; j < VLENS; ++j) {                 \
        (py).c[c][0][0][j] = (fbc)*(y).c[c][0][0][j];   \
        (py).c[c][0][1][j] = (fbc)*(y).c[c][0][1][j];   \
        (py).c[c][1][0][j] = (fbc)*(y).c[c][1][0][j];   \
        (py).c[c][1][1][j] = (fbc)*(y).c[c][1][1][j];   \
      }                                                 \
    }                                                   \
  }

//
// Y-backward spin post reconstruction
// my1 = my1 +   upy1
// my2 = my2 +   upy2
// my3 = my3
// my4 = my4
//
#define __mult_t_back_pst_(my,upy) {\
  _Pragma("loop norecurrence") \
  for (int c = 0; c < 3; ++c) {\
    for (int j = 0; j < VLENS; ++j) {\
      (my).c[c][0][0][j] += (upy).c[c][0][0][j];\
      (my).c[c][0][1][j] += (upy).c[c][0][1][j];\
      (my).c[c][1][0][j] += (upy).c[c][1][0][j];\
      (my).c[c][1][1][j] += (upy).c[c][1][1][j];\
    }\
  }\
}

#define S__mult_t_back_pst_(my, upy, C)                                 \
  fstore1_s(pt, fadd_s(pt, fload1_s(pt, (my).c, dims_scs, C, 0, 0), upy##_##C##_0_0), (my).c, dims_scs, C, 0, 0); \
  fstore1_s(pt, fadd_s(pt, fload1_s(pt, (my).c, dims_scs, C, 0, 1), upy##_##C##_0_1), (my).c, dims_scs, C, 0, 1); \
  fstore1_s(pt, fadd_s(pt, fload1_s(pt, (my).c, dims_scs, C, 1, 0), upy##_##C##_1_0), (my).c, dims_scs, C, 1, 0); \
  fstore1_s(pt, fadd_s(pt, fload1_s(pt, (my).c, dims_scs, C, 1, 1), upy##_##C##_1_1), (my).c, dims_scs, C, 1, 1);  
#define __mult_t_back_pst_vec_(my,upy) LOOP_3(S__mult_t_back_pst_, my, upy)

//
// Y-backward spin post reconstruction with boundary-condition
// my1 = my1 +   upy1 * fbc
// my2 = my2 +   upy2 * fbc
// my3 = my3
// my4 = my4
//
#define __mult_t_back_pst_bc_(my,upy,fbc) {                     \
    for (int c = 0; c < 3; ++c) {                               \
      for (int j = 0; j < VLENS; ++j) {                         \
        (my).c[c][0][0][j] += (upy).c[c][0][0][j] * (fbc);      \
        (my).c[c][0][1][j] += (upy).c[c][0][1][j] * (fbc);      \
        (my).c[c][1][0][j] += (upy).c[c][1][0][j] * (fbc);      \
        (my).c[c][1][1][j] += (upy).c[c][1][1][j] * (fbc);      \
      }                                                         \
    }                                                           \
  }

#define S__mult_t_back_pst_bc_reg_(my, upy, fbc, C) \
        my##_##C##_0_0 = fmadd_s(pt, fload1_s(pt, (upy).c, dims_projscs, C, 0, 0), fdup_s(fbc), my##_##C##_0_0);      \
        my##_##C##_0_1 = fmadd_s(pt, fload1_s(pt, (upy).c, dims_projscs, C, 0, 1), fdup_s(fbc), my##_##C##_0_1);      \
        my##_##C##_1_0 = fmadd_s(pt, fload1_s(pt, (upy).c, dims_projscs, C, 1, 0), fdup_s(fbc), my##_##C##_1_0);      \
        my##_##C##_1_1 = fmadd_s(pt, fload1_s(pt, (upy).c, dims_projscs, C, 1, 1), fdup_s(fbc), my##_##C##_1_1);
#define __mult_t_back_pst_bc_vec_reg_(my,upy,fbc) LOOP_3(S__mult_t_back_pst_bc_reg_, my, upy, fbc)

#endif
