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
#ifndef CLOVER_S_H
#define CLOVER_S_H

#include "util.hh"
#include "addressing.hh"

void __mult_clvs(rvecs_t (* __restrict__ sc)[4][2], const rvecs_t (* __restrict__ a)[36]);

#define __mult_clvs_vec_m(sc, a) {                                      \
    rvecs_t *x;                                                         \
    alloca_aligned(x, sizeof(rvecs_t)*2*6*2, CLS);               \
                                                                        \
    vecs_t ya_00_l;                                                     \
    vecs_t ya_00_r;                                                     \
    vecs_t ya_01_l;                                                     \
    vecs_t ya_01_r;                                                     \
    vecs_t ya_10_l;                                                     \
    vecs_t ya_10_r;                                                     \
    vecs_t ya_11_l;                                                     \
    vecs_t ya_11_r;                                                     \
    vecs_t yb_00_l;                                                     \
    vecs_t yb_00_r;                                                     \
    vecs_t yb_01_l;                                                     \
    vecs_t yb_01_r;                                                     \
    vecs_t yb_10_l;                                                     \
    vecs_t yb_10_r;                                                     \
    vecs_t yb_11_l;                                                     \
    vecs_t yb_11_r;                                                     \
                                                                        \
    pred_t pt = pred_true_all();                                        \
                                                                        \
    BLOCK_START(0);                                                     \
    fstore1_s(pt, fadd_s(pt, fload1_s(pt, sc, dims_scs, 0, 0, 0), fload1_s(pt, sc, dims_scs, 0, 2, 0)), x, dims_x, 0, 0+0, 0); \
    fstore1_s(pt, fadd_s(pt, fload1_s(pt, sc, dims_scs, 0, 1, 0), fload1_s(pt, sc, dims_scs, 0, 3, 0)), x, dims_x, 0, 3+0, 0); \
    fstore1_s(pt, fsub_s(pt, fload1_s(pt, sc, dims_scs, 0, 0, 0), fload1_s(pt, sc, dims_scs, 0, 2, 0)), x, dims_x, 1, 0+0, 0); \
    fstore1_s(pt, fsub_s(pt, fload1_s(pt, sc, dims_scs, 0, 1, 0), fload1_s(pt, sc, dims_scs, 0, 3, 0)), x, dims_x, 1, 3+0, 0); \
    fstore1_s(pt, fadd_s(pt, fload1_s(pt, sc, dims_scs, 0, 0, 1), fload1_s(pt, sc, dims_scs, 0, 2, 1)), x, dims_x, 0, 0+0, 1); \
    fstore1_s(pt, fadd_s(pt, fload1_s(pt, sc, dims_scs, 0, 1, 1), fload1_s(pt, sc, dims_scs, 0, 3, 1)), x, dims_x, 0, 3+0, 1); \
    fstore1_s(pt, fsub_s(pt, fload1_s(pt, sc, dims_scs, 0, 0, 1), fload1_s(pt, sc, dims_scs, 0, 2, 1)), x, dims_x, 1, 0+0, 1); \
    fstore1_s(pt, fsub_s(pt, fload1_s(pt, sc, dims_scs, 0, 1, 1), fload1_s(pt, sc, dims_scs, 0, 3, 1)), x, dims_x, 1, 3+0, 1); \
    fstore1_s(pt, fadd_s(pt, fload1_s(pt, sc, dims_scs, 1, 0, 0), fload1_s(pt, sc, dims_scs, 1, 2, 0)), x, dims_x, 0, 0+1, 0); \
    fstore1_s(pt, fadd_s(pt, fload1_s(pt, sc, dims_scs, 1, 1, 0), fload1_s(pt, sc, dims_scs, 1, 3, 0)), x, dims_x, 0, 3+1, 0); \
    fstore1_s(pt, fsub_s(pt, fload1_s(pt, sc, dims_scs, 1, 0, 0), fload1_s(pt, sc, dims_scs, 1, 2, 0)), x, dims_x, 1, 0+1, 0); \
    fstore1_s(pt, fsub_s(pt, fload1_s(pt, sc, dims_scs, 1, 1, 0), fload1_s(pt, sc, dims_scs, 1, 3, 0)), x, dims_x, 1, 3+1, 0); \
    fstore1_s(pt, fadd_s(pt, fload1_s(pt, sc, dims_scs, 1, 0, 1), fload1_s(pt, sc, dims_scs, 1, 2, 1)), x, dims_x, 0, 0+1, 1); \
    fstore1_s(pt, fadd_s(pt, fload1_s(pt, sc, dims_scs, 1, 1, 1), fload1_s(pt, sc, dims_scs, 1, 3, 1)), x, dims_x, 0, 3+1, 1); \
    fstore1_s(pt, fsub_s(pt, fload1_s(pt, sc, dims_scs, 1, 0, 1), fload1_s(pt, sc, dims_scs, 1, 2, 1)), x, dims_x, 1, 0+1, 1); \
    fstore1_s(pt, fsub_s(pt, fload1_s(pt, sc, dims_scs, 1, 1, 1), fload1_s(pt, sc, dims_scs, 1, 3, 1)), x, dims_x, 1, 3+1, 1); \
    fstore1_s(pt, fadd_s(pt, fload1_s(pt, sc, dims_scs, 2, 0, 0), fload1_s(pt, sc, dims_scs, 2, 2, 0)), x, dims_x, 0, 0+2, 0); \
    fstore1_s(pt, fadd_s(pt, fload1_s(pt, sc, dims_scs, 2, 1, 0), fload1_s(pt, sc, dims_scs, 2, 3, 0)), x, dims_x, 0, 3+2, 0); \
    fstore1_s(pt, fsub_s(pt, fload1_s(pt, sc, dims_scs, 2, 0, 0), fload1_s(pt, sc, dims_scs, 2, 2, 0)), x, dims_x, 1, 0+2, 0); \
    fstore1_s(pt, fsub_s(pt, fload1_s(pt, sc, dims_scs, 2, 1, 0), fload1_s(pt, sc, dims_scs, 2, 3, 0)), x, dims_x, 1, 3+2, 0); \
    fstore1_s(pt, fadd_s(pt, fload1_s(pt, sc, dims_scs, 2, 0, 1), fload1_s(pt, sc, dims_scs, 2, 2, 1)), x, dims_x, 0, 0+2, 1); \
    fstore1_s(pt, fadd_s(pt, fload1_s(pt, sc, dims_scs, 2, 1, 1), fload1_s(pt, sc, dims_scs, 2, 3, 1)), x, dims_x, 0, 3+2, 1); \
    fstore1_s(pt, fsub_s(pt, fload1_s(pt, sc, dims_scs, 2, 0, 1), fload1_s(pt, sc, dims_scs, 2, 2, 1)), x, dims_x, 1, 0+2, 1); \
    fstore1_s(pt, fsub_s(pt, fload1_s(pt, sc, dims_scs, 2, 1, 1), fload1_s(pt, sc, dims_scs, 2, 3, 1)), x, dims_x, 1, 3+2, 1); \
    BLOCK_END(0);                                                       \
                                                                        \
    BLOCK_START(1);                                                     \
    ya_00_l = fzero_s();                                                \
    ya_00_r = fzero_s();                                                \
    ya_01_l = fzero_s();                                                \
    ya_01_r = fzero_s();                                                \
    ya_10_l = fzero_s();                                                \
    ya_10_r = fzero_s();                                                \
    ya_11_l = fzero_s();                                                \
    ya_11_r = fzero_s();                                                \
    yb_00_l = fzero_s();                                                \
    yb_00_r = fzero_s();                                                \
    yb_01_l = fzero_s();                                                \
    yb_01_r = fzero_s();                                                \
    yb_10_l = fzero_s();                                                \
    yb_10_r = fzero_s();                                                \
    yb_11_l = fzero_s();                                                \
    yb_11_r = fzero_s();                                                \
                                                                        \
    ya_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 0), fload1_s(pt, x, dims_x, 0, 0, 0), ya_00_l); \
    ya_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 0), fload1_s(pt, x, dims_x, 0, 0, 1), ya_01_l); \
    ya_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 0), fload1_s(pt, x, dims_x, 1, 0, 0), ya_10_l); \
    ya_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 0), fload1_s(pt, x, dims_x, 1, 0, 1), ya_11_l); \
    yb_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 1), fload1_s(pt, x, dims_x, 0, 1, 0), yb_00_l); \
    yb_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 1), fload1_s(pt, x, dims_x, 0, 1, 1), yb_01_l); \
    yb_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 1), fload1_s(pt, x, dims_x, 1, 1, 0), yb_10_l); \
    yb_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 1), fload1_s(pt, x, dims_x, 1, 1, 1), yb_11_l); \
                                                                        \
    ya_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 6), fload1_s(pt, x, dims_x, 0, 1, 0), ya_00_l); \
    ya_00_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 7), fload1_s(pt, x, dims_x, 0, 1, 1), ya_00_r); \
    ya_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 6), fload1_s(pt, x, dims_x, 0, 1, 1), ya_01_l); \
    ya_01_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 7), fload1_s(pt, x, dims_x, 0, 1, 0), ya_01_r); \
    ya_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 6), fload1_s(pt, x, dims_x, 1, 1, 0), ya_10_l); \
    ya_10_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 7), fload1_s(pt, x, dims_x, 1, 1, 1), ya_10_r); \
    ya_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 6), fload1_s(pt, x, dims_x, 1, 1, 1), ya_11_l); \
    ya_11_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 7), fload1_s(pt, x, dims_x, 1, 1, 0), ya_11_r); \
    yb_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 6), fload1_s(pt, x, dims_x, 0, 0, 0), yb_00_l); \
    yb_00_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 7), fload1_s(pt, x, dims_x, 0, 0, 1), yb_00_r); \
    yb_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 6), fload1_s(pt, x, dims_x, 0, 0, 1), yb_01_l); \
    yb_01_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 7), fload1_s(pt, x, dims_x, 0, 0, 0), yb_01_r); \
    yb_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 6), fload1_s(pt, x, dims_x, 1, 0, 0), yb_10_l); \
    yb_10_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 7), fload1_s(pt, x, dims_x, 1, 0, 1), yb_10_r); \
    yb_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 6), fload1_s(pt, x, dims_x, 1, 0, 1), yb_11_l); \
    yb_11_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 7), fload1_s(pt, x, dims_x, 1, 0, 0), yb_11_r); \
                                                                        \
    ya_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 8), fload1_s(pt, x, dims_x, 0, 2, 0), ya_00_l); \
    ya_00_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 9), fload1_s(pt, x, dims_x, 0, 2, 1), ya_00_r); \
    ya_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 8), fload1_s(pt, x, dims_x, 0, 2, 1), ya_01_l); \
    ya_01_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 9), fload1_s(pt, x, dims_x, 0, 2, 0), ya_01_r); \
    ya_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 8), fload1_s(pt, x, dims_x, 1, 2, 0), ya_10_l); \
    ya_10_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 9), fload1_s(pt, x, dims_x, 1, 2, 1), ya_10_r); \
    ya_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 8), fload1_s(pt, x, dims_x, 1, 2, 1), ya_11_l); \
    ya_11_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 9), fload1_s(pt, x, dims_x, 1, 2, 0), ya_11_r); \
    yb_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 16), fload1_s(pt, x, dims_x, 0, 2, 0), yb_00_l); \
    yb_00_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 17), fload1_s(pt, x, dims_x, 0, 2, 1), yb_00_r); \
    yb_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 16), fload1_s(pt, x, dims_x, 0, 2, 1), yb_01_l); \
    yb_01_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 17), fload1_s(pt, x, dims_x, 0, 2, 0), yb_01_r); \
    yb_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 16), fload1_s(pt, x, dims_x, 1, 2, 0), yb_10_l); \
    yb_10_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 17), fload1_s(pt, x, dims_x, 1, 2, 1), yb_10_r); \
    yb_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 16), fload1_s(pt, x, dims_x, 1, 2, 1), yb_11_l); \
    yb_11_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 17), fload1_s(pt, x, dims_x, 1, 2, 0), yb_11_r); \
                                                                        \
    ya_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 10), fload1_s(pt, x, dims_x, 0, 3, 0), ya_00_l); \
    ya_00_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 11), fload1_s(pt, x, dims_x, 0, 3, 1), ya_00_r); \
    ya_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 10), fload1_s(pt, x, dims_x, 0, 3, 1), ya_01_l); \
    ya_01_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 11), fload1_s(pt, x, dims_x, 0, 3, 0), ya_01_r); \
    ya_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 10), fload1_s(pt, x, dims_x, 1, 3, 0), ya_10_l); \
    ya_10_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 11), fload1_s(pt, x, dims_x, 1, 3, 1), ya_10_r); \
    ya_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 10), fload1_s(pt, x, dims_x, 1, 3, 1), ya_11_l); \
    ya_11_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 11), fload1_s(pt, x, dims_x, 1, 3, 0), ya_11_r); \
    yb_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 18), fload1_s(pt, x, dims_x, 0, 3, 0), yb_00_l); \
    yb_00_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 19), fload1_s(pt, x, dims_x, 0, 3, 1), yb_00_r); \
    yb_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 18), fload1_s(pt, x, dims_x, 0, 3, 1), yb_01_l); \
    yb_01_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 19), fload1_s(pt, x, dims_x, 0, 3, 0), yb_01_r); \
    yb_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 18), fload1_s(pt, x, dims_x, 1, 3, 0), yb_10_l); \
    yb_10_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 19), fload1_s(pt, x, dims_x, 1, 3, 1), yb_10_r); \
    yb_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 18), fload1_s(pt, x, dims_x, 1, 3, 1), yb_11_l); \
    yb_11_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 19), fload1_s(pt, x, dims_x, 1, 3, 0), yb_11_r); \
                                                                        \
    ya_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 12), fload1_s(pt, x, dims_x, 0, 4, 0), ya_00_l); \
    ya_00_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 13), fload1_s(pt, x, dims_x, 0, 4, 1), ya_00_r); \
    ya_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 12), fload1_s(pt, x, dims_x, 0, 4, 1), ya_01_l); \
    ya_01_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 13), fload1_s(pt, x, dims_x, 0, 4, 0), ya_01_r); \
    ya_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 12), fload1_s(pt, x, dims_x, 1, 4, 0), ya_10_l); \
    ya_10_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 13), fload1_s(pt, x, dims_x, 1, 4, 1), ya_10_r); \
    ya_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 12), fload1_s(pt, x, dims_x, 1, 4, 1), ya_11_l); \
    ya_11_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 13), fload1_s(pt, x, dims_x, 1, 4, 0), ya_11_r); \
    yb_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 20), fload1_s(pt, x, dims_x, 0, 4, 0), yb_00_l); \
    yb_00_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 21), fload1_s(pt, x, dims_x, 0, 4, 1), yb_00_r); \
    yb_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 20), fload1_s(pt, x, dims_x, 0, 4, 1), yb_01_l); \
    yb_01_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 21), fload1_s(pt, x, dims_x, 0, 4, 0), yb_01_r); \
    yb_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 20), fload1_s(pt, x, dims_x, 1, 4, 0), yb_10_l); \
    yb_10_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 21), fload1_s(pt, x, dims_x, 1, 4, 1), yb_10_r); \
    yb_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 20), fload1_s(pt, x, dims_x, 1, 4, 1), yb_11_l); \
    yb_11_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 21), fload1_s(pt, x, dims_x, 1, 4, 0), yb_11_r); \
                                                                        \
    ya_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 14), fload1_s(pt, x, dims_x, 0, 5, 0), ya_00_l); \
    ya_00_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 15), fload1_s(pt, x, dims_x, 0, 5, 1), ya_00_r); \
    ya_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 14), fload1_s(pt, x, dims_x, 0, 5, 1), ya_01_l); \
    ya_01_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 15), fload1_s(pt, x, dims_x, 0, 5, 0), ya_01_r); \
    ya_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 14), fload1_s(pt, x, dims_x, 1, 5, 0), ya_10_l); \
    ya_10_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 15), fload1_s(pt, x, dims_x, 1, 5, 1), ya_10_r); \
    ya_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 14), fload1_s(pt, x, dims_x, 1, 5, 1), ya_11_l); \
    ya_11_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 15), fload1_s(pt, x, dims_x, 1, 5, 0), ya_11_r); \
    yb_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 22), fload1_s(pt, x, dims_x, 0, 5, 0), yb_00_l); \
    yb_00_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 23), fload1_s(pt, x, dims_x, 0, 5, 1), yb_00_r); \
    yb_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 22), fload1_s(pt, x, dims_x, 0, 5, 1), yb_01_l); \
    yb_01_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 23), fload1_s(pt, x, dims_x, 0, 5, 0), yb_01_r); \
    yb_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 22), fload1_s(pt, x, dims_x, 1, 5, 0), yb_10_l); \
    yb_10_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 23), fload1_s(pt, x, dims_x, 1, 5, 1), yb_10_r); \
    yb_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 22), fload1_s(pt, x, dims_x, 1, 5, 1), yb_11_l); \
    yb_11_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 23), fload1_s(pt, x, dims_x, 1, 5, 0), yb_11_r); \
                                                                        \
    ya_00_l = fadd_s(pt, ya_00_l, ya_00_r);                             \
    ya_01_l = fadd_s(pt, ya_01_l, ya_01_r);                             \
    ya_10_l = fadd_s(pt, ya_10_l, ya_10_r);                             \
    ya_11_l = fadd_s(pt, ya_11_l, ya_11_r);                             \
    yb_00_l = fadd_s(pt, yb_00_l, yb_00_r);                             \
    yb_01_l = fadd_s(pt, yb_01_l, yb_01_r);                             \
    yb_10_l = fadd_s(pt, yb_10_l, yb_10_r);                             \
    yb_11_l = fadd_s(pt, yb_11_l, yb_11_r);                             \
                                                                        \
    fstore1_s(pt, fadd_s(pt, ya_00_l, ya_10_l), sc, dims_scs, 0, 0, 0); \
    fstore1_s(pt, fadd_s(pt, ya_01_l, ya_11_l), sc, dims_scs, 0, 0, 1); \
    fstore1_s(pt, fsub_s(pt, ya_00_l, ya_10_l), sc, dims_scs, 0, 2, 0); \
    fstore1_s(pt, fsub_s(pt, ya_01_l, ya_11_l), sc, dims_scs, 0, 2, 1); \
    fstore1_s(pt, fadd_s(pt, yb_00_l, yb_10_l), sc, dims_scs, 1, 0, 0); \
    fstore1_s(pt, fadd_s(pt, yb_01_l, yb_11_l), sc, dims_scs, 1, 0, 1); \
    fstore1_s(pt, fsub_s(pt, yb_00_l, yb_10_l), sc, dims_scs, 1, 2, 0); \
    fstore1_s(pt, fsub_s(pt, yb_01_l, yb_11_l), sc, dims_scs, 1, 2, 1); \
    BLOCK_END(1);                                                       \
                                                                        \
    BLOCK_START(2);                                                     \
    ya_00_l = fzero_s();                                                \
    ya_00_r = fzero_s();                                                \
    ya_01_l = fzero_s();                                                \
    ya_01_r = fzero_s();                                                \
    ya_10_l = fzero_s();                                                \
    ya_10_r = fzero_s();                                                \
    ya_11_l = fzero_s();                                                \
    ya_11_r = fzero_s();                                                \
    yb_00_l = fzero_s();                                                \
    yb_00_r = fzero_s();                                                \
    yb_01_l = fzero_s();                                                \
    yb_01_r = fzero_s();                                                \
    yb_10_l = fzero_s();                                                \
    yb_10_r = fzero_s();                                                \
    yb_11_l = fzero_s();                                                \
    yb_11_r = fzero_s();                                                \
                                                                        \
    ya_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 8), fload1_s(pt, x, dims_x, 0, 0, 0), ya_00_l); \
    ya_00_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 9), fload1_s(pt, x, dims_x, 0, 0, 1), ya_00_r); \
    ya_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 8), fload1_s(pt, x, dims_x, 0, 0, 1), ya_01_l); \
    ya_01_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 9), fload1_s(pt, x, dims_x, 0, 0, 0), ya_01_r); \
    ya_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 8), fload1_s(pt, x, dims_x, 1, 0, 0), ya_10_l); \
    ya_10_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 9), fload1_s(pt, x, dims_x, 1, 0, 1), ya_10_r); \
    ya_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 8), fload1_s(pt, x, dims_x, 1, 0, 1), ya_11_l); \
    ya_11_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 9), fload1_s(pt, x, dims_x, 1, 0, 0), ya_11_r); \
    yb_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 10), fload1_s(pt, x, dims_x, 0, 0, 0), yb_00_l); \
    yb_00_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 11), fload1_s(pt, x, dims_x, 0, 0, 1), yb_00_r); \
    yb_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 10), fload1_s(pt, x, dims_x, 0, 0, 1), yb_01_l); \
    yb_01_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 11), fload1_s(pt, x, dims_x, 0, 0, 0), yb_01_r); \
    yb_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 10), fload1_s(pt, x, dims_x, 1, 0, 0), yb_10_l); \
    yb_10_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 11), fload1_s(pt, x, dims_x, 1, 0, 1), yb_10_r); \
    yb_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 10), fload1_s(pt, x, dims_x, 1, 0, 1), yb_11_l); \
    yb_11_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 11), fload1_s(pt, x, dims_x, 1, 0, 0), yb_11_r); \
                                                                        \
    ya_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 16), fload1_s(pt, x, dims_x, 0, 1, 0), ya_00_l); \
    ya_00_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 17), fload1_s(pt, x, dims_x, 0, 1, 1), ya_00_r); \
    ya_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 16), fload1_s(pt, x, dims_x, 0, 1, 1), ya_01_l); \
    ya_01_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 17), fload1_s(pt, x, dims_x, 0, 1, 0), ya_01_r); \
    ya_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 16), fload1_s(pt, x, dims_x, 1, 1, 0), ya_10_l); \
    ya_10_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 17), fload1_s(pt, x, dims_x, 1, 1, 1), ya_10_r); \
    ya_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 16), fload1_s(pt, x, dims_x, 1, 1, 1), ya_11_l); \
    ya_11_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 17), fload1_s(pt, x, dims_x, 1, 1, 0), ya_11_r); \
    yb_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 18), fload1_s(pt, x, dims_x, 0, 1, 0), yb_00_l); \
    yb_00_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 19), fload1_s(pt, x, dims_x, 0, 1, 1), yb_00_r); \
    yb_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 18), fload1_s(pt, x, dims_x, 0, 1, 1), yb_01_l); \
    yb_01_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 19), fload1_s(pt, x, dims_x, 0, 1, 0), yb_01_r); \
    yb_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 18), fload1_s(pt, x, dims_x, 1, 1, 0), yb_10_l); \
    yb_10_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 19), fload1_s(pt, x, dims_x, 1, 1, 1), yb_10_r); \
    yb_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 18), fload1_s(pt, x, dims_x, 1, 1, 1), yb_11_l); \
    yb_11_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 19), fload1_s(pt, x, dims_x, 1, 1, 0), yb_11_r); \
                                                                        \
    ya_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 2), fload1_s(pt, x, dims_x, 0, 2, 0), ya_00_l); \
    ya_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 2), fload1_s(pt, x, dims_x, 0, 2, 1), ya_01_l); \
    ya_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 2), fload1_s(pt, x, dims_x, 1, 2, 0), ya_10_l); \
    ya_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 2), fload1_s(pt, x, dims_x, 1, 2, 1), ya_11_l); \
    yb_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 3), fload1_s(pt, x, dims_x, 0, 3, 0), yb_00_l); \
    yb_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 3), fload1_s(pt, x, dims_x, 0, 3, 1), yb_01_l); \
    yb_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 3), fload1_s(pt, x, dims_x, 1, 3, 0), yb_10_l); \
    yb_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 3), fload1_s(pt, x, dims_x, 1, 3, 1), yb_11_l); \
                                                                        \
    ya_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 24), fload1_s(pt, x, dims_x, 0, 3, 0), ya_00_l); \
    ya_00_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 25), fload1_s(pt, x, dims_x, 0, 3, 1), ya_00_r); \
    ya_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 24), fload1_s(pt, x, dims_x, 0, 3, 1), ya_01_l); \
    ya_01_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 25), fload1_s(pt, x, dims_x, 0, 3, 0), ya_01_r); \
    ya_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 24), fload1_s(pt, x, dims_x, 1, 3, 0), ya_10_l); \
    ya_10_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 25), fload1_s(pt, x, dims_x, 1, 3, 1), ya_10_r); \
    ya_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 24), fload1_s(pt, x, dims_x, 1, 3, 1), ya_11_l); \
    ya_11_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 25), fload1_s(pt, x, dims_x, 1, 3, 0), ya_11_r); \
    yb_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 24), fload1_s(pt, x, dims_x, 0, 2, 0), yb_00_l); \
    yb_00_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 25), fload1_s(pt, x, dims_x, 0, 2, 1), yb_00_r); \
    yb_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 24), fload1_s(pt, x, dims_x, 0, 2, 1), yb_01_l); \
    yb_01_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 25), fload1_s(pt, x, dims_x, 0, 2, 0), yb_01_r); \
    yb_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 24), fload1_s(pt, x, dims_x, 1, 2, 0), yb_10_l); \
    yb_10_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 25), fload1_s(pt, x, dims_x, 1, 2, 1), yb_10_r); \
    yb_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 24), fload1_s(pt, x, dims_x, 1, 2, 1), yb_11_l); \
    yb_11_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 25), fload1_s(pt, x, dims_x, 1, 2, 0), yb_11_r); \
                                                                        \
    ya_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 26), fload1_s(pt, x, dims_x, 0, 4, 0), ya_00_l); \
    ya_00_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 27), fload1_s(pt, x, dims_x, 0, 4, 1), ya_00_r); \
    ya_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 26), fload1_s(pt, x, dims_x, 0, 4, 1), ya_01_l); \
    ya_01_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 27), fload1_s(pt, x, dims_x, 0, 4, 0), ya_01_r); \
    ya_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 26), fload1_s(pt, x, dims_x, 1, 4, 0), ya_10_l); \
    ya_10_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 27), fload1_s(pt, x, dims_x, 1, 4, 1), ya_10_r); \
    ya_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 26), fload1_s(pt, x, dims_x, 1, 4, 1), ya_11_l); \
    ya_11_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 27), fload1_s(pt, x, dims_x, 1, 4, 0), ya_11_r); \
    yb_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 30), fload1_s(pt, x, dims_x, 0, 4, 0), yb_00_l); \
    yb_00_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 31), fload1_s(pt, x, dims_x, 0, 4, 1), yb_00_r); \
    yb_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 30), fload1_s(pt, x, dims_x, 0, 4, 1), yb_01_l); \
    yb_01_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 31), fload1_s(pt, x, dims_x, 0, 4, 0), yb_01_r); \
    yb_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 30), fload1_s(pt, x, dims_x, 1, 4, 0), yb_10_l); \
    yb_10_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 31), fload1_s(pt, x, dims_x, 1, 4, 1), yb_10_r); \
    yb_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 30), fload1_s(pt, x, dims_x, 1, 4, 1), yb_11_l); \
    yb_11_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 31), fload1_s(pt, x, dims_x, 1, 4, 0), yb_11_r); \
                                                                        \
    ya_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 28), fload1_s(pt, x, dims_x, 0, 5, 0), ya_00_l); \
    ya_00_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 29), fload1_s(pt, x, dims_x, 0, 5, 1), ya_00_r); \
    ya_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 28), fload1_s(pt, x, dims_x, 0, 5, 1), ya_01_l); \
    ya_01_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 29), fload1_s(pt, x, dims_x, 0, 5, 0), ya_01_r); \
    ya_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 28), fload1_s(pt, x, dims_x, 1, 5, 0), ya_10_l); \
    ya_10_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 29), fload1_s(pt, x, dims_x, 1, 5, 1), ya_10_r); \
    ya_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 28), fload1_s(pt, x, dims_x, 1, 5, 1), ya_11_l); \
    ya_11_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 29), fload1_s(pt, x, dims_x, 1, 5, 0), ya_11_r); \
    yb_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 32), fload1_s(pt, x, dims_x, 0, 5, 0), yb_00_l); \
    yb_00_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 33), fload1_s(pt, x, dims_x, 0, 5, 1), yb_00_r); \
    yb_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 32), fload1_s(pt, x, dims_x, 0, 5, 1), yb_01_l); \
    yb_01_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 33), fload1_s(pt, x, dims_x, 0, 5, 0), yb_01_r); \
    yb_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 32), fload1_s(pt, x, dims_x, 1, 5, 0), yb_10_l); \
    yb_10_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 33), fload1_s(pt, x, dims_x, 1, 5, 1), yb_10_r); \
    yb_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 32), fload1_s(pt, x, dims_x, 1, 5, 1), yb_11_l); \
    yb_11_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 33), fload1_s(pt, x, dims_x, 1, 5, 0), yb_11_r); \
                                                                        \
    ya_00_l = fadd_s(pt, ya_00_l, ya_00_r);                             \
    ya_01_l = fadd_s(pt, ya_01_l, ya_01_r);                             \
    ya_10_l = fadd_s(pt, ya_10_l, ya_10_r);                             \
    ya_11_l = fadd_s(pt, ya_11_l, ya_11_r);                             \
    yb_00_l = fadd_s(pt, yb_00_l, yb_00_r);                             \
    yb_01_l = fadd_s(pt, yb_01_l, yb_01_r);                             \
    yb_10_l = fadd_s(pt, yb_10_l, yb_10_r);                             \
    yb_11_l = fadd_s(pt, yb_11_l, yb_11_r);                             \
                                                                        \
    fstore1_s(pt, fadd_s(pt, ya_00_l, ya_10_l), sc, dims_scs, 2, 0, 0); \
    fstore1_s(pt, fadd_s(pt, ya_01_l, ya_11_l), sc, dims_scs, 2, 0, 1); \
    fstore1_s(pt, fsub_s(pt, ya_00_l, ya_10_l), sc, dims_scs, 2, 2, 0); \
    fstore1_s(pt, fsub_s(pt, ya_01_l, ya_11_l), sc, dims_scs, 2, 2, 1); \
    fstore1_s(pt, fadd_s(pt, yb_00_l, yb_10_l), sc, dims_scs, 0, 1, 0); \
    fstore1_s(pt, fadd_s(pt, yb_01_l, yb_11_l), sc, dims_scs, 0, 1, 1); \
    fstore1_s(pt, fsub_s(pt, yb_00_l, yb_10_l), sc, dims_scs, 0, 3, 0); \
    fstore1_s(pt, fsub_s(pt, yb_01_l, yb_11_l), sc, dims_scs, 0, 3, 1); \
    BLOCK_END(2);                                                       \
                                                                        \
    BLOCK_START(3);                                                     \
    ya_00_l = fzero_s();                                                \
    ya_00_r = fzero_s();                                                \
    ya_01_l = fzero_s();                                                \
    ya_01_r = fzero_s();                                                \
    ya_10_l = fzero_s();                                                \
    ya_10_r = fzero_s();                                                \
    ya_11_l = fzero_s();                                                \
    ya_11_r = fzero_s();                                                \
    yb_00_l = fzero_s();                                                \
    yb_00_r = fzero_s();                                                \
    yb_01_l = fzero_s();                                                \
    yb_01_r = fzero_s();                                                \
    yb_10_l = fzero_s();                                                \
    yb_10_r = fzero_s();                                                \
    yb_11_l = fzero_s();                                                \
    yb_11_r = fzero_s();                                                \
                                                                        \
    ya_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 12), fload1_s(pt, x, dims_x, 0, 0, 0), ya_00_l); \
    ya_00_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 13), fload1_s(pt, x, dims_x, 0, 0, 1), ya_00_r); \
    ya_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 12), fload1_s(pt, x, dims_x, 0, 0, 1), ya_01_l); \
    ya_01_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 13), fload1_s(pt, x, dims_x, 0, 0, 0), ya_01_r); \
    ya_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 12), fload1_s(pt, x, dims_x, 1, 0, 0), ya_10_l); \
    ya_10_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 13), fload1_s(pt, x, dims_x, 1, 0, 1), ya_10_r); \
    ya_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 12), fload1_s(pt, x, dims_x, 1, 0, 1), ya_11_l); \
    ya_11_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 13), fload1_s(pt, x, dims_x, 1, 0, 0), ya_11_r); \
    yb_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 14), fload1_s(pt, x, dims_x, 0, 0, 0), yb_00_l); \
    yb_00_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 15), fload1_s(pt, x, dims_x, 0, 0, 1), yb_00_r); \
    yb_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 14), fload1_s(pt, x, dims_x, 0, 0, 1), yb_01_l); \
    yb_01_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 15), fload1_s(pt, x, dims_x, 0, 0, 0), yb_01_r); \
    yb_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 14), fload1_s(pt, x, dims_x, 1, 0, 0), yb_10_l); \
    yb_10_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 15), fload1_s(pt, x, dims_x, 1, 0, 1), yb_10_r); \
    yb_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 14), fload1_s(pt, x, dims_x, 1, 0, 1), yb_11_l); \
    yb_11_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 15), fload1_s(pt, x, dims_x, 1, 0, 0), yb_11_r); \
                                                                        \
    ya_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 20), fload1_s(pt, x, dims_x, 0, 1, 0), ya_00_l); \
    ya_00_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 21), fload1_s(pt, x, dims_x, 0, 1, 1), ya_00_r); \
    ya_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 20), fload1_s(pt, x, dims_x, 0, 1, 1), ya_01_l); \
    ya_01_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 21), fload1_s(pt, x, dims_x, 0, 1, 0), ya_01_r); \
    ya_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 20), fload1_s(pt, x, dims_x, 1, 1, 0), ya_10_l); \
    ya_10_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 21), fload1_s(pt, x, dims_x, 1, 1, 1), ya_10_r); \
    ya_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 20), fload1_s(pt, x, dims_x, 1, 1, 1), ya_11_l); \
    ya_11_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 21), fload1_s(pt, x, dims_x, 1, 1, 0), ya_11_r); \
    yb_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 22), fload1_s(pt, x, dims_x, 0, 1, 0), yb_00_l); \
    yb_00_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 23), fload1_s(pt, x, dims_x, 0, 1, 1), yb_00_r); \
    yb_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 22), fload1_s(pt, x, dims_x, 0, 1, 1), yb_01_l); \
    yb_01_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 23), fload1_s(pt, x, dims_x, 0, 1, 0), yb_01_r); \
    yb_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 22), fload1_s(pt, x, dims_x, 1, 1, 0), yb_10_l); \
    yb_10_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 23), fload1_s(pt, x, dims_x, 1, 1, 1), yb_10_r); \
    yb_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 22), fload1_s(pt, x, dims_x, 1, 1, 1), yb_11_l); \
    yb_11_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 23), fload1_s(pt, x, dims_x, 1, 1, 0), yb_11_r); \
                                                                        \
    ya_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 26), fload1_s(pt, x, dims_x, 0, 2, 0), ya_00_l); \
    ya_00_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 27), fload1_s(pt, x, dims_x, 0, 2, 1), ya_00_r); \
    ya_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 26), fload1_s(pt, x, dims_x, 0, 2, 1), ya_01_l); \
    ya_01_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 27), fload1_s(pt, x, dims_x, 0, 2, 0), ya_01_r); \
    ya_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 26), fload1_s(pt, x, dims_x, 1, 2, 0), ya_10_l); \
    ya_10_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 27), fload1_s(pt, x, dims_x, 1, 2, 1), ya_10_r); \
    ya_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 26), fload1_s(pt, x, dims_x, 1, 2, 1), ya_11_l); \
    ya_11_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 27), fload1_s(pt, x, dims_x, 1, 2, 0), ya_11_r); \
    yb_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 28), fload1_s(pt, x, dims_x, 0, 2, 0), yb_00_l); \
    yb_00_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 29), fload1_s(pt, x, dims_x, 0, 2, 1), yb_00_r); \
    yb_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 28), fload1_s(pt, x, dims_x, 0, 2, 1), yb_01_l); \
    yb_01_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 29), fload1_s(pt, x, dims_x, 0, 2, 0), yb_01_r); \
    yb_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 28), fload1_s(pt, x, dims_x, 1, 2, 0), yb_10_l); \
    yb_10_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 29), fload1_s(pt, x, dims_x, 1, 2, 1), yb_10_r); \
    yb_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 28), fload1_s(pt, x, dims_x, 1, 2, 1), yb_11_l); \
    yb_11_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 29), fload1_s(pt, x, dims_x, 1, 2, 0), yb_11_r); \
                                                                        \
    ya_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 30), fload1_s(pt, x, dims_x, 0, 3, 0), ya_00_l); \
    ya_00_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 31), fload1_s(pt, x, dims_x, 0, 3, 1), ya_00_r); \
    ya_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 30), fload1_s(pt, x, dims_x, 0, 3, 1), ya_01_l); \
    ya_01_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 31), fload1_s(pt, x, dims_x, 0, 3, 0), ya_01_r); \
    ya_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 30), fload1_s(pt, x, dims_x, 1, 3, 0), ya_10_l); \
    ya_10_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 31), fload1_s(pt, x, dims_x, 1, 3, 1), ya_10_r); \
    ya_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 30), fload1_s(pt, x, dims_x, 1, 3, 1), ya_11_l); \
    ya_11_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 31), fload1_s(pt, x, dims_x, 1, 3, 0), ya_11_r); \
    yb_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 32), fload1_s(pt, x, dims_x, 0, 3, 0), yb_00_l); \
    yb_00_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 33), fload1_s(pt, x, dims_x, 0, 3, 1), yb_00_r); \
    yb_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 32), fload1_s(pt, x, dims_x, 0, 3, 1), yb_01_l); \
    yb_01_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 33), fload1_s(pt, x, dims_x, 0, 3, 0), yb_01_r); \
    yb_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 32), fload1_s(pt, x, dims_x, 1, 3, 0), yb_10_l); \
    yb_10_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 33), fload1_s(pt, x, dims_x, 1, 3, 1), yb_10_r); \
    yb_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 32), fload1_s(pt, x, dims_x, 1, 3, 1), yb_11_l); \
    yb_11_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 33), fload1_s(pt, x, dims_x, 1, 3, 0), yb_11_r); \
                                                                        \
    ya_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 4), fload1_s(pt, x, dims_x, 0, 4, 0), ya_00_l); \
    ya_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 4), fload1_s(pt, x, dims_x, 0, 4, 1), ya_01_l); \
    ya_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 4), fload1_s(pt, x, dims_x, 1, 4, 0), ya_10_l); \
    ya_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 4), fload1_s(pt, x, dims_x, 1, 4, 1), ya_11_l); \
    yb_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 5), fload1_s(pt, x, dims_x, 0, 5, 0), yb_00_l); \
    yb_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 5), fload1_s(pt, x, dims_x, 0, 5, 1), yb_01_l); \
    yb_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 5), fload1_s(pt, x, dims_x, 1, 5, 0), yb_10_l); \
    yb_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 5), fload1_s(pt, x, dims_x, 1, 5, 1), yb_11_l); \
                                                                        \
    ya_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 34), fload1_s(pt, x, dims_x, 0, 5, 0), ya_00_l); \
    ya_00_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 35), fload1_s(pt, x, dims_x, 0, 5, 1), ya_00_r); \
    ya_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 34), fload1_s(pt, x, dims_x, 0, 5, 1), ya_01_l); \
    ya_01_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 35), fload1_s(pt, x, dims_x, 0, 5, 0), ya_01_r); \
    ya_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 34), fload1_s(pt, x, dims_x, 1, 5, 0), ya_10_l); \
    ya_10_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 35), fload1_s(pt, x, dims_x, 1, 5, 1), ya_10_r); \
    ya_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 34), fload1_s(pt, x, dims_x, 1, 5, 1), ya_11_l); \
    ya_11_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 35), fload1_s(pt, x, dims_x, 1, 5, 0), ya_11_r); \
    yb_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 34), fload1_s(pt, x, dims_x, 0, 4, 0), yb_00_l); \
    yb_00_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 35), fload1_s(pt, x, dims_x, 0, 4, 1), yb_00_r); \
    yb_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 34), fload1_s(pt, x, dims_x, 0, 4, 1), yb_01_l); \
    yb_01_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 35), fload1_s(pt, x, dims_x, 0, 4, 0), yb_01_r); \
    yb_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 34), fload1_s(pt, x, dims_x, 1, 4, 0), yb_10_l); \
    yb_10_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 35), fload1_s(pt, x, dims_x, 1, 4, 1), yb_10_r); \
    yb_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 34), fload1_s(pt, x, dims_x, 1, 4, 1), yb_11_l); \
    yb_11_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 35), fload1_s(pt, x, dims_x, 1, 4, 0), yb_11_r); \
                                                                        \
    ya_00_l = fadd_s(pt, ya_00_l, ya_00_r);                             \
    ya_01_l = fadd_s(pt, ya_01_l, ya_01_r);                             \
    ya_10_l = fadd_s(pt, ya_10_l, ya_10_r);                             \
    ya_11_l = fadd_s(pt, ya_11_l, ya_11_r);                             \
    yb_00_l = fadd_s(pt, yb_00_l, yb_00_r);                             \
    yb_01_l = fadd_s(pt, yb_01_l, yb_01_r);                             \
    yb_10_l = fadd_s(pt, yb_10_l, yb_10_r);                             \
    yb_11_l = fadd_s(pt, yb_11_l, yb_11_r);                             \
                                                                        \
    fstore1_s(pt, fadd_s(pt, ya_00_l, ya_10_l), sc, dims_scs, 1, 1, 0); \
    fstore1_s(pt, fadd_s(pt, ya_01_l, ya_11_l), sc, dims_scs, 1, 1, 1); \
    fstore1_s(pt, fsub_s(pt, ya_00_l, ya_10_l), sc, dims_scs, 1, 3, 0); \
    fstore1_s(pt, fsub_s(pt, ya_01_l, ya_11_l), sc, dims_scs, 1, 3, 1); \
    fstore1_s(pt, fadd_s(pt, yb_00_l, yb_10_l), sc, dims_scs, 2, 1, 0); \
    fstore1_s(pt, fadd_s(pt, yb_01_l, yb_11_l), sc, dims_scs, 2, 1, 1); \
    fstore1_s(pt, fsub_s(pt, yb_00_l, yb_10_l), sc, dims_scs, 2, 3, 0); \
    fstore1_s(pt, fsub_s(pt, yb_01_l, yb_11_l), sc, dims_scs, 2, 3, 1); \
    BLOCK_END(3);                                                       \
  }                                                                     \

#define __mult_clvs_vec_m_reg(sc, sc_reg_in, a) {                       \
    rvecs_t *x;                                                         \
    alloca_aligned(x, sizeof(rvecs_t)*2*6*2, CLS);               \
                                                                        \
    vecs_t ya_00_l;                                                     \
    vecs_t ya_00_r;                                                     \
    vecs_t ya_01_l;                                                     \
    vecs_t ya_01_r;                                                     \
    vecs_t ya_10_l;                                                     \
    vecs_t ya_10_r;                                                     \
    vecs_t ya_11_l;                                                     \
    vecs_t ya_11_r;                                                     \
    vecs_t yb_00_l;                                                     \
    vecs_t yb_00_r;                                                     \
    vecs_t yb_01_l;                                                     \
    vecs_t yb_01_r;                                                     \
    vecs_t yb_10_l;                                                     \
    vecs_t yb_10_r;                                                     \
    vecs_t yb_11_l;                                                     \
    vecs_t yb_11_r;                                                     \
                                                                        \
    pred_t pt = pred_true_all();                                        \
                                                                        \
    BLOCK_START(0);                                                     \
    fstore1_s(pt, fadd_s(pt, sc_reg_in##_0_0_0, sc_reg_in##_0_2_0), x, dims_x, 0, 0+0, 0); \
    fstore1_s(pt, fadd_s(pt, sc_reg_in##_0_1_0, sc_reg_in##_0_3_0), x, dims_x, 0, 3+0, 0); \
    fstore1_s(pt, fsub_s(pt, sc_reg_in##_0_0_0, sc_reg_in##_0_2_0), x, dims_x, 1, 0+0, 0); \
    fstore1_s(pt, fsub_s(pt, sc_reg_in##_0_1_0, sc_reg_in##_0_3_0), x, dims_x, 1, 3+0, 0); \
    fstore1_s(pt, fadd_s(pt, sc_reg_in##_0_0_1, sc_reg_in##_0_2_1), x, dims_x, 0, 0+0, 1); \
    fstore1_s(pt, fadd_s(pt, sc_reg_in##_0_1_1, sc_reg_in##_0_3_1), x, dims_x, 0, 3+0, 1); \
    fstore1_s(pt, fsub_s(pt, sc_reg_in##_0_0_1, sc_reg_in##_0_2_1), x, dims_x, 1, 0+0, 1); \
    fstore1_s(pt, fsub_s(pt, sc_reg_in##_0_1_1, sc_reg_in##_0_3_1), x, dims_x, 1, 3+0, 1); \
    fstore1_s(pt, fadd_s(pt, sc_reg_in##_1_0_0, sc_reg_in##_1_2_0), x, dims_x, 0, 0+1, 0); \
    fstore1_s(pt, fadd_s(pt, sc_reg_in##_1_1_0, sc_reg_in##_1_3_0), x, dims_x, 0, 3+1, 0); \
    fstore1_s(pt, fsub_s(pt, sc_reg_in##_1_0_0, sc_reg_in##_1_2_0), x, dims_x, 1, 0+1, 0); \
    fstore1_s(pt, fsub_s(pt, sc_reg_in##_1_1_0, sc_reg_in##_1_3_0), x, dims_x, 1, 3+1, 0); \
    fstore1_s(pt, fadd_s(pt, sc_reg_in##_1_0_1, sc_reg_in##_1_2_1), x, dims_x, 0, 0+1, 1); \
    fstore1_s(pt, fadd_s(pt, sc_reg_in##_1_1_1, sc_reg_in##_1_3_1), x, dims_x, 0, 3+1, 1); \
    fstore1_s(pt, fsub_s(pt, sc_reg_in##_1_0_1, sc_reg_in##_1_2_1), x, dims_x, 1, 0+1, 1); \
    fstore1_s(pt, fsub_s(pt, sc_reg_in##_1_1_1, sc_reg_in##_1_3_1), x, dims_x, 1, 3+1, 1); \
    fstore1_s(pt, fadd_s(pt, sc_reg_in##_2_0_0, sc_reg_in##_2_2_0), x, dims_x, 0, 0+2, 0); \
    fstore1_s(pt, fadd_s(pt, sc_reg_in##_2_1_0, sc_reg_in##_2_3_0), x, dims_x, 0, 3+2, 0); \
    fstore1_s(pt, fsub_s(pt, sc_reg_in##_2_0_0, sc_reg_in##_2_2_0), x, dims_x, 1, 0+2, 0); \
    fstore1_s(pt, fsub_s(pt, sc_reg_in##_2_1_0, sc_reg_in##_2_3_0), x, dims_x, 1, 3+2, 0); \
    fstore1_s(pt, fadd_s(pt, sc_reg_in##_2_0_1, sc_reg_in##_2_2_1), x, dims_x, 0, 0+2, 1); \
    fstore1_s(pt, fadd_s(pt, sc_reg_in##_2_1_1, sc_reg_in##_2_3_1), x, dims_x, 0, 3+2, 1); \
    fstore1_s(pt, fsub_s(pt, sc_reg_in##_2_0_1, sc_reg_in##_2_2_1), x, dims_x, 1, 0+2, 1); \
    fstore1_s(pt, fsub_s(pt, sc_reg_in##_2_1_1, sc_reg_in##_2_3_1), x, dims_x, 1, 3+2, 1); \
    BLOCK_END(0);                                                       \
                                                                        \
    BLOCK_START(1);                                                     \
    ya_00_l = fzero_s();                                                \
    ya_00_r = fzero_s();                                                \
    ya_01_l = fzero_s();                                                \
    ya_01_r = fzero_s();                                                \
    ya_10_l = fzero_s();                                                \
    ya_10_r = fzero_s();                                                \
    ya_11_l = fzero_s();                                                \
    ya_11_r = fzero_s();                                                \
    yb_00_l = fzero_s();                                                \
    yb_00_r = fzero_s();                                                \
    yb_01_l = fzero_s();                                                \
    yb_01_r = fzero_s();                                                \
    yb_10_l = fzero_s();                                                \
    yb_10_r = fzero_s();                                                \
    yb_11_l = fzero_s();                                                \
    yb_11_r = fzero_s();                                                \
                                                                        \
    ya_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 0), fload1_s(pt, x, dims_x, 0, 0, 0), ya_00_l); \
    ya_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 0), fload1_s(pt, x, dims_x, 0, 0, 1), ya_01_l); \
    ya_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 0), fload1_s(pt, x, dims_x, 1, 0, 0), ya_10_l); \
    ya_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 0), fload1_s(pt, x, dims_x, 1, 0, 1), ya_11_l); \
    yb_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 1), fload1_s(pt, x, dims_x, 0, 1, 0), yb_00_l); \
    yb_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 1), fload1_s(pt, x, dims_x, 0, 1, 1), yb_01_l); \
    yb_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 1), fload1_s(pt, x, dims_x, 1, 1, 0), yb_10_l); \
    yb_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 1), fload1_s(pt, x, dims_x, 1, 1, 1), yb_11_l); \
                                                                        \
    ya_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 6), fload1_s(pt, x, dims_x, 0, 1, 0), ya_00_l); \
    ya_00_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 7), fload1_s(pt, x, dims_x, 0, 1, 1), ya_00_r); \
    ya_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 6), fload1_s(pt, x, dims_x, 0, 1, 1), ya_01_l); \
    ya_01_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 7), fload1_s(pt, x, dims_x, 0, 1, 0), ya_01_r); \
    ya_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 6), fload1_s(pt, x, dims_x, 1, 1, 0), ya_10_l); \
    ya_10_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 7), fload1_s(pt, x, dims_x, 1, 1, 1), ya_10_r); \
    ya_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 6), fload1_s(pt, x, dims_x, 1, 1, 1), ya_11_l); \
    ya_11_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 7), fload1_s(pt, x, dims_x, 1, 1, 0), ya_11_r); \
    yb_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 6), fload1_s(pt, x, dims_x, 0, 0, 0), yb_00_l); \
    yb_00_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 7), fload1_s(pt, x, dims_x, 0, 0, 1), yb_00_r); \
    yb_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 6), fload1_s(pt, x, dims_x, 0, 0, 1), yb_01_l); \
    yb_01_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 7), fload1_s(pt, x, dims_x, 0, 0, 0), yb_01_r); \
    yb_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 6), fload1_s(pt, x, dims_x, 1, 0, 0), yb_10_l); \
    yb_10_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 7), fload1_s(pt, x, dims_x, 1, 0, 1), yb_10_r); \
    yb_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 6), fload1_s(pt, x, dims_x, 1, 0, 1), yb_11_l); \
    yb_11_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 7), fload1_s(pt, x, dims_x, 1, 0, 0), yb_11_r); \
                                                                        \
    ya_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 8), fload1_s(pt, x, dims_x, 0, 2, 0), ya_00_l); \
    ya_00_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 9), fload1_s(pt, x, dims_x, 0, 2, 1), ya_00_r); \
    ya_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 8), fload1_s(pt, x, dims_x, 0, 2, 1), ya_01_l); \
    ya_01_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 9), fload1_s(pt, x, dims_x, 0, 2, 0), ya_01_r); \
    ya_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 8), fload1_s(pt, x, dims_x, 1, 2, 0), ya_10_l); \
    ya_10_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 9), fload1_s(pt, x, dims_x, 1, 2, 1), ya_10_r); \
    ya_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 8), fload1_s(pt, x, dims_x, 1, 2, 1), ya_11_l); \
    ya_11_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 9), fload1_s(pt, x, dims_x, 1, 2, 0), ya_11_r); \
    yb_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 16), fload1_s(pt, x, dims_x, 0, 2, 0), yb_00_l); \
    yb_00_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 17), fload1_s(pt, x, dims_x, 0, 2, 1), yb_00_r); \
    yb_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 16), fload1_s(pt, x, dims_x, 0, 2, 1), yb_01_l); \
    yb_01_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 17), fload1_s(pt, x, dims_x, 0, 2, 0), yb_01_r); \
    yb_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 16), fload1_s(pt, x, dims_x, 1, 2, 0), yb_10_l); \
    yb_10_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 17), fload1_s(pt, x, dims_x, 1, 2, 1), yb_10_r); \
    yb_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 16), fload1_s(pt, x, dims_x, 1, 2, 1), yb_11_l); \
    yb_11_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 17), fload1_s(pt, x, dims_x, 1, 2, 0), yb_11_r); \
                                                                        \
    ya_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 10), fload1_s(pt, x, dims_x, 0, 3, 0), ya_00_l); \
    ya_00_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 11), fload1_s(pt, x, dims_x, 0, 3, 1), ya_00_r); \
    ya_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 10), fload1_s(pt, x, dims_x, 0, 3, 1), ya_01_l); \
    ya_01_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 11), fload1_s(pt, x, dims_x, 0, 3, 0), ya_01_r); \
    ya_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 10), fload1_s(pt, x, dims_x, 1, 3, 0), ya_10_l); \
    ya_10_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 11), fload1_s(pt, x, dims_x, 1, 3, 1), ya_10_r); \
    ya_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 10), fload1_s(pt, x, dims_x, 1, 3, 1), ya_11_l); \
    ya_11_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 11), fload1_s(pt, x, dims_x, 1, 3, 0), ya_11_r); \
    yb_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 18), fload1_s(pt, x, dims_x, 0, 3, 0), yb_00_l); \
    yb_00_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 19), fload1_s(pt, x, dims_x, 0, 3, 1), yb_00_r); \
    yb_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 18), fload1_s(pt, x, dims_x, 0, 3, 1), yb_01_l); \
    yb_01_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 19), fload1_s(pt, x, dims_x, 0, 3, 0), yb_01_r); \
    yb_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 18), fload1_s(pt, x, dims_x, 1, 3, 0), yb_10_l); \
    yb_10_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 19), fload1_s(pt, x, dims_x, 1, 3, 1), yb_10_r); \
    yb_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 18), fload1_s(pt, x, dims_x, 1, 3, 1), yb_11_l); \
    yb_11_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 19), fload1_s(pt, x, dims_x, 1, 3, 0), yb_11_r); \
                                                                        \
    ya_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 12), fload1_s(pt, x, dims_x, 0, 4, 0), ya_00_l); \
    ya_00_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 13), fload1_s(pt, x, dims_x, 0, 4, 1), ya_00_r); \
    ya_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 12), fload1_s(pt, x, dims_x, 0, 4, 1), ya_01_l); \
    ya_01_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 13), fload1_s(pt, x, dims_x, 0, 4, 0), ya_01_r); \
    ya_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 12), fload1_s(pt, x, dims_x, 1, 4, 0), ya_10_l); \
    ya_10_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 13), fload1_s(pt, x, dims_x, 1, 4, 1), ya_10_r); \
    ya_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 12), fload1_s(pt, x, dims_x, 1, 4, 1), ya_11_l); \
    ya_11_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 13), fload1_s(pt, x, dims_x, 1, 4, 0), ya_11_r); \
    yb_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 20), fload1_s(pt, x, dims_x, 0, 4, 0), yb_00_l); \
    yb_00_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 21), fload1_s(pt, x, dims_x, 0, 4, 1), yb_00_r); \
    yb_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 20), fload1_s(pt, x, dims_x, 0, 4, 1), yb_01_l); \
    yb_01_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 21), fload1_s(pt, x, dims_x, 0, 4, 0), yb_01_r); \
    yb_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 20), fload1_s(pt, x, dims_x, 1, 4, 0), yb_10_l); \
    yb_10_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 21), fload1_s(pt, x, dims_x, 1, 4, 1), yb_10_r); \
    yb_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 20), fload1_s(pt, x, dims_x, 1, 4, 1), yb_11_l); \
    yb_11_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 21), fload1_s(pt, x, dims_x, 1, 4, 0), yb_11_r); \
                                                                        \
    ya_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 14), fload1_s(pt, x, dims_x, 0, 5, 0), ya_00_l); \
    ya_00_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 15), fload1_s(pt, x, dims_x, 0, 5, 1), ya_00_r); \
    ya_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 14), fload1_s(pt, x, dims_x, 0, 5, 1), ya_01_l); \
    ya_01_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 15), fload1_s(pt, x, dims_x, 0, 5, 0), ya_01_r); \
    ya_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 14), fload1_s(pt, x, dims_x, 1, 5, 0), ya_10_l); \
    ya_10_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 15), fload1_s(pt, x, dims_x, 1, 5, 1), ya_10_r); \
    ya_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 14), fload1_s(pt, x, dims_x, 1, 5, 1), ya_11_l); \
    ya_11_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 15), fload1_s(pt, x, dims_x, 1, 5, 0), ya_11_r); \
    yb_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 22), fload1_s(pt, x, dims_x, 0, 5, 0), yb_00_l); \
    yb_00_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 23), fload1_s(pt, x, dims_x, 0, 5, 1), yb_00_r); \
    yb_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 22), fload1_s(pt, x, dims_x, 0, 5, 1), yb_01_l); \
    yb_01_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 23), fload1_s(pt, x, dims_x, 0, 5, 0), yb_01_r); \
    yb_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 22), fload1_s(pt, x, dims_x, 1, 5, 0), yb_10_l); \
    yb_10_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 23), fload1_s(pt, x, dims_x, 1, 5, 1), yb_10_r); \
    yb_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 22), fload1_s(pt, x, dims_x, 1, 5, 1), yb_11_l); \
    yb_11_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 23), fload1_s(pt, x, dims_x, 1, 5, 0), yb_11_r); \
                                                                        \
    ya_00_l = fadd_s(pt, ya_00_l, ya_00_r);                             \
    ya_01_l = fadd_s(pt, ya_01_l, ya_01_r);                             \
    ya_10_l = fadd_s(pt, ya_10_l, ya_10_r);                             \
    ya_11_l = fadd_s(pt, ya_11_l, ya_11_r);                             \
    yb_00_l = fadd_s(pt, yb_00_l, yb_00_r);                             \
    yb_01_l = fadd_s(pt, yb_01_l, yb_01_r);                             \
    yb_10_l = fadd_s(pt, yb_10_l, yb_10_r);                             \
    yb_11_l = fadd_s(pt, yb_11_l, yb_11_r);                             \
                                                                        \
    fstore1_s(pt, fadd_s(pt, ya_00_l, ya_10_l), sc, dims_scs, 0, 0, 0); \
    fstore1_s(pt, fadd_s(pt, ya_01_l, ya_11_l), sc, dims_scs, 0, 0, 1); \
    fstore1_s(pt, fsub_s(pt, ya_00_l, ya_10_l), sc, dims_scs, 0, 2, 0); \
    fstore1_s(pt, fsub_s(pt, ya_01_l, ya_11_l), sc, dims_scs, 0, 2, 1); \
    fstore1_s(pt, fadd_s(pt, yb_00_l, yb_10_l), sc, dims_scs, 1, 0, 0); \
    fstore1_s(pt, fadd_s(pt, yb_01_l, yb_11_l), sc, dims_scs, 1, 0, 1); \
    fstore1_s(pt, fsub_s(pt, yb_00_l, yb_10_l), sc, dims_scs, 1, 2, 0); \
    fstore1_s(pt, fsub_s(pt, yb_01_l, yb_11_l), sc, dims_scs, 1, 2, 1); \
    BLOCK_END(1);                                                       \
                                                                        \
    BLOCK_START(2);                                                     \
    ya_00_l = fzero_s();                                                \
    ya_00_r = fzero_s();                                                \
    ya_01_l = fzero_s();                                                \
    ya_01_r = fzero_s();                                                \
    ya_10_l = fzero_s();                                                \
    ya_10_r = fzero_s();                                                \
    ya_11_l = fzero_s();                                                \
    ya_11_r = fzero_s();                                                \
    yb_00_l = fzero_s();                                                \
    yb_00_r = fzero_s();                                                \
    yb_01_l = fzero_s();                                                \
    yb_01_r = fzero_s();                                                \
    yb_10_l = fzero_s();                                                \
    yb_10_r = fzero_s();                                                \
    yb_11_l = fzero_s();                                                \
    yb_11_r = fzero_s();                                                \
                                                                        \
    ya_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 8), fload1_s(pt, x, dims_x, 0, 0, 0), ya_00_l); \
    ya_00_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 9), fload1_s(pt, x, dims_x, 0, 0, 1), ya_00_r); \
    ya_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 8), fload1_s(pt, x, dims_x, 0, 0, 1), ya_01_l); \
    ya_01_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 9), fload1_s(pt, x, dims_x, 0, 0, 0), ya_01_r); \
    ya_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 8), fload1_s(pt, x, dims_x, 1, 0, 0), ya_10_l); \
    ya_10_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 9), fload1_s(pt, x, dims_x, 1, 0, 1), ya_10_r); \
    ya_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 8), fload1_s(pt, x, dims_x, 1, 0, 1), ya_11_l); \
    ya_11_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 9), fload1_s(pt, x, dims_x, 1, 0, 0), ya_11_r); \
    yb_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 10), fload1_s(pt, x, dims_x, 0, 0, 0), yb_00_l); \
    yb_00_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 11), fload1_s(pt, x, dims_x, 0, 0, 1), yb_00_r); \
    yb_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 10), fload1_s(pt, x, dims_x, 0, 0, 1), yb_01_l); \
    yb_01_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 11), fload1_s(pt, x, dims_x, 0, 0, 0), yb_01_r); \
    yb_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 10), fload1_s(pt, x, dims_x, 1, 0, 0), yb_10_l); \
    yb_10_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 11), fload1_s(pt, x, dims_x, 1, 0, 1), yb_10_r); \
    yb_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 10), fload1_s(pt, x, dims_x, 1, 0, 1), yb_11_l); \
    yb_11_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 11), fload1_s(pt, x, dims_x, 1, 0, 0), yb_11_r); \
                                                                        \
    ya_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 16), fload1_s(pt, x, dims_x, 0, 1, 0), ya_00_l); \
    ya_00_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 17), fload1_s(pt, x, dims_x, 0, 1, 1), ya_00_r); \
    ya_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 16), fload1_s(pt, x, dims_x, 0, 1, 1), ya_01_l); \
    ya_01_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 17), fload1_s(pt, x, dims_x, 0, 1, 0), ya_01_r); \
    ya_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 16), fload1_s(pt, x, dims_x, 1, 1, 0), ya_10_l); \
    ya_10_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 17), fload1_s(pt, x, dims_x, 1, 1, 1), ya_10_r); \
    ya_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 16), fload1_s(pt, x, dims_x, 1, 1, 1), ya_11_l); \
    ya_11_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 17), fload1_s(pt, x, dims_x, 1, 1, 0), ya_11_r); \
    yb_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 18), fload1_s(pt, x, dims_x, 0, 1, 0), yb_00_l); \
    yb_00_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 19), fload1_s(pt, x, dims_x, 0, 1, 1), yb_00_r); \
    yb_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 18), fload1_s(pt, x, dims_x, 0, 1, 1), yb_01_l); \
    yb_01_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 19), fload1_s(pt, x, dims_x, 0, 1, 0), yb_01_r); \
    yb_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 18), fload1_s(pt, x, dims_x, 1, 1, 0), yb_10_l); \
    yb_10_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 19), fload1_s(pt, x, dims_x, 1, 1, 1), yb_10_r); \
    yb_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 18), fload1_s(pt, x, dims_x, 1, 1, 1), yb_11_l); \
    yb_11_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 19), fload1_s(pt, x, dims_x, 1, 1, 0), yb_11_r); \
                                                                        \
    ya_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 2), fload1_s(pt, x, dims_x, 0, 2, 0), ya_00_l); \
    ya_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 2), fload1_s(pt, x, dims_x, 0, 2, 1), ya_01_l); \
    ya_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 2), fload1_s(pt, x, dims_x, 1, 2, 0), ya_10_l); \
    ya_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 2), fload1_s(pt, x, dims_x, 1, 2, 1), ya_11_l); \
    yb_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 3), fload1_s(pt, x, dims_x, 0, 3, 0), yb_00_l); \
    yb_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 3), fload1_s(pt, x, dims_x, 0, 3, 1), yb_01_l); \
    yb_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 3), fload1_s(pt, x, dims_x, 1, 3, 0), yb_10_l); \
    yb_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 3), fload1_s(pt, x, dims_x, 1, 3, 1), yb_11_l); \
                                                                        \
    ya_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 24), fload1_s(pt, x, dims_x, 0, 3, 0), ya_00_l); \
    ya_00_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 25), fload1_s(pt, x, dims_x, 0, 3, 1), ya_00_r); \
    ya_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 24), fload1_s(pt, x, dims_x, 0, 3, 1), ya_01_l); \
    ya_01_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 25), fload1_s(pt, x, dims_x, 0, 3, 0), ya_01_r); \
    ya_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 24), fload1_s(pt, x, dims_x, 1, 3, 0), ya_10_l); \
    ya_10_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 25), fload1_s(pt, x, dims_x, 1, 3, 1), ya_10_r); \
    ya_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 24), fload1_s(pt, x, dims_x, 1, 3, 1), ya_11_l); \
    ya_11_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 25), fload1_s(pt, x, dims_x, 1, 3, 0), ya_11_r); \
    yb_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 24), fload1_s(pt, x, dims_x, 0, 2, 0), yb_00_l); \
    yb_00_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 25), fload1_s(pt, x, dims_x, 0, 2, 1), yb_00_r); \
    yb_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 24), fload1_s(pt, x, dims_x, 0, 2, 1), yb_01_l); \
    yb_01_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 25), fload1_s(pt, x, dims_x, 0, 2, 0), yb_01_r); \
    yb_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 24), fload1_s(pt, x, dims_x, 1, 2, 0), yb_10_l); \
    yb_10_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 25), fload1_s(pt, x, dims_x, 1, 2, 1), yb_10_r); \
    yb_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 24), fload1_s(pt, x, dims_x, 1, 2, 1), yb_11_l); \
    yb_11_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 25), fload1_s(pt, x, dims_x, 1, 2, 0), yb_11_r); \
                                                                        \
    ya_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 26), fload1_s(pt, x, dims_x, 0, 4, 0), ya_00_l); \
    ya_00_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 27), fload1_s(pt, x, dims_x, 0, 4, 1), ya_00_r); \
    ya_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 26), fload1_s(pt, x, dims_x, 0, 4, 1), ya_01_l); \
    ya_01_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 27), fload1_s(pt, x, dims_x, 0, 4, 0), ya_01_r); \
    ya_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 26), fload1_s(pt, x, dims_x, 1, 4, 0), ya_10_l); \
    ya_10_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 27), fload1_s(pt, x, dims_x, 1, 4, 1), ya_10_r); \
    ya_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 26), fload1_s(pt, x, dims_x, 1, 4, 1), ya_11_l); \
    ya_11_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 27), fload1_s(pt, x, dims_x, 1, 4, 0), ya_11_r); \
    yb_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 30), fload1_s(pt, x, dims_x, 0, 4, 0), yb_00_l); \
    yb_00_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 31), fload1_s(pt, x, dims_x, 0, 4, 1), yb_00_r); \
    yb_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 30), fload1_s(pt, x, dims_x, 0, 4, 1), yb_01_l); \
    yb_01_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 31), fload1_s(pt, x, dims_x, 0, 4, 0), yb_01_r); \
    yb_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 30), fload1_s(pt, x, dims_x, 1, 4, 0), yb_10_l); \
    yb_10_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 31), fload1_s(pt, x, dims_x, 1, 4, 1), yb_10_r); \
    yb_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 30), fload1_s(pt, x, dims_x, 1, 4, 1), yb_11_l); \
    yb_11_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 31), fload1_s(pt, x, dims_x, 1, 4, 0), yb_11_r); \
                                                                        \
    ya_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 28), fload1_s(pt, x, dims_x, 0, 5, 0), ya_00_l); \
    ya_00_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 29), fload1_s(pt, x, dims_x, 0, 5, 1), ya_00_r); \
    ya_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 28), fload1_s(pt, x, dims_x, 0, 5, 1), ya_01_l); \
    ya_01_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 29), fload1_s(pt, x, dims_x, 0, 5, 0), ya_01_r); \
    ya_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 28), fload1_s(pt, x, dims_x, 1, 5, 0), ya_10_l); \
    ya_10_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 29), fload1_s(pt, x, dims_x, 1, 5, 1), ya_10_r); \
    ya_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 28), fload1_s(pt, x, dims_x, 1, 5, 1), ya_11_l); \
    ya_11_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 29), fload1_s(pt, x, dims_x, 1, 5, 0), ya_11_r); \
    yb_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 32), fload1_s(pt, x, dims_x, 0, 5, 0), yb_00_l); \
    yb_00_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 33), fload1_s(pt, x, dims_x, 0, 5, 1), yb_00_r); \
    yb_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 32), fload1_s(pt, x, dims_x, 0, 5, 1), yb_01_l); \
    yb_01_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 33), fload1_s(pt, x, dims_x, 0, 5, 0), yb_01_r); \
    yb_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 32), fload1_s(pt, x, dims_x, 1, 5, 0), yb_10_l); \
    yb_10_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 33), fload1_s(pt, x, dims_x, 1, 5, 1), yb_10_r); \
    yb_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 32), fload1_s(pt, x, dims_x, 1, 5, 1), yb_11_l); \
    yb_11_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 33), fload1_s(pt, x, dims_x, 1, 5, 0), yb_11_r); \
                                                                        \
    ya_00_l = fadd_s(pt, ya_00_l, ya_00_r);                             \
    ya_01_l = fadd_s(pt, ya_01_l, ya_01_r);                             \
    ya_10_l = fadd_s(pt, ya_10_l, ya_10_r);                             \
    ya_11_l = fadd_s(pt, ya_11_l, ya_11_r);                             \
    yb_00_l = fadd_s(pt, yb_00_l, yb_00_r);                             \
    yb_01_l = fadd_s(pt, yb_01_l, yb_01_r);                             \
    yb_10_l = fadd_s(pt, yb_10_l, yb_10_r);                             \
    yb_11_l = fadd_s(pt, yb_11_l, yb_11_r);                             \
                                                                        \
    fstore1_s(pt, fadd_s(pt, ya_00_l, ya_10_l), sc, dims_scs, 2, 0, 0); \
    fstore1_s(pt, fadd_s(pt, ya_01_l, ya_11_l), sc, dims_scs, 2, 0, 1); \
    fstore1_s(pt, fsub_s(pt, ya_00_l, ya_10_l), sc, dims_scs, 2, 2, 0); \
    fstore1_s(pt, fsub_s(pt, ya_01_l, ya_11_l), sc, dims_scs, 2, 2, 1); \
    fstore1_s(pt, fadd_s(pt, yb_00_l, yb_10_l), sc, dims_scs, 0, 1, 0); \
    fstore1_s(pt, fadd_s(pt, yb_01_l, yb_11_l), sc, dims_scs, 0, 1, 1); \
    fstore1_s(pt, fsub_s(pt, yb_00_l, yb_10_l), sc, dims_scs, 0, 3, 0); \
    fstore1_s(pt, fsub_s(pt, yb_01_l, yb_11_l), sc, dims_scs, 0, 3, 1); \
    BLOCK_END(2);                                                       \
                                                                        \
    BLOCK_START(3);                                                     \
    ya_00_l = fzero_s();                                                \
    ya_00_r = fzero_s();                                                \
    ya_01_l = fzero_s();                                                \
    ya_01_r = fzero_s();                                                \
    ya_10_l = fzero_s();                                                \
    ya_10_r = fzero_s();                                                \
    ya_11_l = fzero_s();                                                \
    ya_11_r = fzero_s();                                                \
    yb_00_l = fzero_s();                                                \
    yb_00_r = fzero_s();                                                \
    yb_01_l = fzero_s();                                                \
    yb_01_r = fzero_s();                                                \
    yb_10_l = fzero_s();                                                \
    yb_10_r = fzero_s();                                                \
    yb_11_l = fzero_s();                                                \
    yb_11_r = fzero_s();                                                \
                                                                        \
    ya_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 12), fload1_s(pt, x, dims_x, 0, 0, 0), ya_00_l); \
    ya_00_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 13), fload1_s(pt, x, dims_x, 0, 0, 1), ya_00_r); \
    ya_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 12), fload1_s(pt, x, dims_x, 0, 0, 1), ya_01_l); \
    ya_01_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 13), fload1_s(pt, x, dims_x, 0, 0, 0), ya_01_r); \
    ya_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 12), fload1_s(pt, x, dims_x, 1, 0, 0), ya_10_l); \
    ya_10_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 13), fload1_s(pt, x, dims_x, 1, 0, 1), ya_10_r); \
    ya_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 12), fload1_s(pt, x, dims_x, 1, 0, 1), ya_11_l); \
    ya_11_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 13), fload1_s(pt, x, dims_x, 1, 0, 0), ya_11_r); \
    yb_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 14), fload1_s(pt, x, dims_x, 0, 0, 0), yb_00_l); \
    yb_00_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 15), fload1_s(pt, x, dims_x, 0, 0, 1), yb_00_r); \
    yb_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 14), fload1_s(pt, x, dims_x, 0, 0, 1), yb_01_l); \
    yb_01_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 15), fload1_s(pt, x, dims_x, 0, 0, 0), yb_01_r); \
    yb_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 14), fload1_s(pt, x, dims_x, 1, 0, 0), yb_10_l); \
    yb_10_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 15), fload1_s(pt, x, dims_x, 1, 0, 1), yb_10_r); \
    yb_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 14), fload1_s(pt, x, dims_x, 1, 0, 1), yb_11_l); \
    yb_11_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 15), fload1_s(pt, x, dims_x, 1, 0, 0), yb_11_r); \
                                                                        \
    ya_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 20), fload1_s(pt, x, dims_x, 0, 1, 0), ya_00_l); \
    ya_00_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 21), fload1_s(pt, x, dims_x, 0, 1, 1), ya_00_r); \
    ya_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 20), fload1_s(pt, x, dims_x, 0, 1, 1), ya_01_l); \
    ya_01_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 21), fload1_s(pt, x, dims_x, 0, 1, 0), ya_01_r); \
    ya_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 20), fload1_s(pt, x, dims_x, 1, 1, 0), ya_10_l); \
    ya_10_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 21), fload1_s(pt, x, dims_x, 1, 1, 1), ya_10_r); \
    ya_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 20), fload1_s(pt, x, dims_x, 1, 1, 1), ya_11_l); \
    ya_11_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 21), fload1_s(pt, x, dims_x, 1, 1, 0), ya_11_r); \
    yb_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 22), fload1_s(pt, x, dims_x, 0, 1, 0), yb_00_l); \
    yb_00_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 23), fload1_s(pt, x, dims_x, 0, 1, 1), yb_00_r); \
    yb_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 22), fload1_s(pt, x, dims_x, 0, 1, 1), yb_01_l); \
    yb_01_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 23), fload1_s(pt, x, dims_x, 0, 1, 0), yb_01_r); \
    yb_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 22), fload1_s(pt, x, dims_x, 1, 1, 0), yb_10_l); \
    yb_10_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 23), fload1_s(pt, x, dims_x, 1, 1, 1), yb_10_r); \
    yb_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 22), fload1_s(pt, x, dims_x, 1, 1, 1), yb_11_l); \
    yb_11_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 23), fload1_s(pt, x, dims_x, 1, 1, 0), yb_11_r); \
                                                                        \
    ya_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 26), fload1_s(pt, x, dims_x, 0, 2, 0), ya_00_l); \
    ya_00_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 27), fload1_s(pt, x, dims_x, 0, 2, 1), ya_00_r); \
    ya_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 26), fload1_s(pt, x, dims_x, 0, 2, 1), ya_01_l); \
    ya_01_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 27), fload1_s(pt, x, dims_x, 0, 2, 0), ya_01_r); \
    ya_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 26), fload1_s(pt, x, dims_x, 1, 2, 0), ya_10_l); \
    ya_10_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 27), fload1_s(pt, x, dims_x, 1, 2, 1), ya_10_r); \
    ya_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 26), fload1_s(pt, x, dims_x, 1, 2, 1), ya_11_l); \
    ya_11_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 27), fload1_s(pt, x, dims_x, 1, 2, 0), ya_11_r); \
    yb_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 28), fload1_s(pt, x, dims_x, 0, 2, 0), yb_00_l); \
    yb_00_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 29), fload1_s(pt, x, dims_x, 0, 2, 1), yb_00_r); \
    yb_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 28), fload1_s(pt, x, dims_x, 0, 2, 1), yb_01_l); \
    yb_01_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 29), fload1_s(pt, x, dims_x, 0, 2, 0), yb_01_r); \
    yb_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 28), fload1_s(pt, x, dims_x, 1, 2, 0), yb_10_l); \
    yb_10_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 29), fload1_s(pt, x, dims_x, 1, 2, 1), yb_10_r); \
    yb_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 28), fload1_s(pt, x, dims_x, 1, 2, 1), yb_11_l); \
    yb_11_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 29), fload1_s(pt, x, dims_x, 1, 2, 0), yb_11_r); \
                                                                        \
    ya_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 30), fload1_s(pt, x, dims_x, 0, 3, 0), ya_00_l); \
    ya_00_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 31), fload1_s(pt, x, dims_x, 0, 3, 1), ya_00_r); \
    ya_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 30), fload1_s(pt, x, dims_x, 0, 3, 1), ya_01_l); \
    ya_01_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 31), fload1_s(pt, x, dims_x, 0, 3, 0), ya_01_r); \
    ya_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 30), fload1_s(pt, x, dims_x, 1, 3, 0), ya_10_l); \
    ya_10_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 31), fload1_s(pt, x, dims_x, 1, 3, 1), ya_10_r); \
    ya_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 30), fload1_s(pt, x, dims_x, 1, 3, 1), ya_11_l); \
    ya_11_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 31), fload1_s(pt, x, dims_x, 1, 3, 0), ya_11_r); \
    yb_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 32), fload1_s(pt, x, dims_x, 0, 3, 0), yb_00_l); \
    yb_00_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 33), fload1_s(pt, x, dims_x, 0, 3, 1), yb_00_r); \
    yb_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 32), fload1_s(pt, x, dims_x, 0, 3, 1), yb_01_l); \
    yb_01_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 33), fload1_s(pt, x, dims_x, 0, 3, 0), yb_01_r); \
    yb_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 32), fload1_s(pt, x, dims_x, 1, 3, 0), yb_10_l); \
    yb_10_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 33), fload1_s(pt, x, dims_x, 1, 3, 1), yb_10_r); \
    yb_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 32), fload1_s(pt, x, dims_x, 1, 3, 1), yb_11_l); \
    yb_11_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 33), fload1_s(pt, x, dims_x, 1, 3, 0), yb_11_r); \
                                                                        \
    ya_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 4), fload1_s(pt, x, dims_x, 0, 4, 0), ya_00_l); \
    ya_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 4), fload1_s(pt, x, dims_x, 0, 4, 1), ya_01_l); \
    ya_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 4), fload1_s(pt, x, dims_x, 1, 4, 0), ya_10_l); \
    ya_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 4), fload1_s(pt, x, dims_x, 1, 4, 1), ya_11_l); \
    yb_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 5), fload1_s(pt, x, dims_x, 0, 5, 0), yb_00_l); \
    yb_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 5), fload1_s(pt, x, dims_x, 0, 5, 1), yb_01_l); \
    yb_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 5), fload1_s(pt, x, dims_x, 1, 5, 0), yb_10_l); \
    yb_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 5), fload1_s(pt, x, dims_x, 1, 5, 1), yb_11_l); \
                                                                        \
    ya_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 34), fload1_s(pt, x, dims_x, 0, 5, 0), ya_00_l); \
    ya_00_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 35), fload1_s(pt, x, dims_x, 0, 5, 1), ya_00_r); \
    ya_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 34), fload1_s(pt, x, dims_x, 0, 5, 1), ya_01_l); \
    ya_01_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 35), fload1_s(pt, x, dims_x, 0, 5, 0), ya_01_r); \
    ya_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 34), fload1_s(pt, x, dims_x, 1, 5, 0), ya_10_l); \
    ya_10_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 35), fload1_s(pt, x, dims_x, 1, 5, 1), ya_10_r); \
    ya_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 34), fload1_s(pt, x, dims_x, 1, 5, 1), ya_11_l); \
    ya_11_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 35), fload1_s(pt, x, dims_x, 1, 5, 0), ya_11_r); \
    yb_00_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 34), fload1_s(pt, x, dims_x, 0, 4, 0), yb_00_l); \
    yb_00_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 35), fload1_s(pt, x, dims_x, 0, 4, 1), yb_00_r); \
    yb_01_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 34), fload1_s(pt, x, dims_x, 0, 4, 1), yb_01_l); \
    yb_01_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 0, 35), fload1_s(pt, x, dims_x, 0, 4, 0), yb_01_r); \
    yb_10_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 34), fload1_s(pt, x, dims_x, 1, 4, 0), yb_10_l); \
    yb_10_r = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 35), fload1_s(pt, x, dims_x, 1, 4, 1), yb_10_r); \
    yb_11_l = fmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 34), fload1_s(pt, x, dims_x, 1, 4, 1), yb_11_l); \
    yb_11_r = fnmadd_s(pt, fload1_s(pt, a, dims_clv, 1, 35), fload1_s(pt, x, dims_x, 1, 4, 0), yb_11_r); \
                                                                        \
    ya_00_l = fadd_s(pt, ya_00_l, ya_00_r);                             \
    ya_01_l = fadd_s(pt, ya_01_l, ya_01_r);                             \
    ya_10_l = fadd_s(pt, ya_10_l, ya_10_r);                             \
    ya_11_l = fadd_s(pt, ya_11_l, ya_11_r);                             \
    yb_00_l = fadd_s(pt, yb_00_l, yb_00_r);                             \
    yb_01_l = fadd_s(pt, yb_01_l, yb_01_r);                             \
    yb_10_l = fadd_s(pt, yb_10_l, yb_10_r);                             \
    yb_11_l = fadd_s(pt, yb_11_l, yb_11_r);                             \
                                                                        \
    fstore1_s(pt, fadd_s(pt, ya_00_l, ya_10_l), sc, dims_scs, 1, 1, 0); \
    fstore1_s(pt, fadd_s(pt, ya_01_l, ya_11_l), sc, dims_scs, 1, 1, 1); \
    fstore1_s(pt, fsub_s(pt, ya_00_l, ya_10_l), sc, dims_scs, 1, 3, 0); \
    fstore1_s(pt, fsub_s(pt, ya_01_l, ya_11_l), sc, dims_scs, 1, 3, 1); \
    fstore1_s(pt, fadd_s(pt, yb_00_l, yb_10_l), sc, dims_scs, 2, 1, 0); \
    fstore1_s(pt, fadd_s(pt, yb_01_l, yb_11_l), sc, dims_scs, 2, 1, 1); \
    fstore1_s(pt, fsub_s(pt, yb_00_l, yb_10_l), sc, dims_scs, 2, 3, 0); \
    fstore1_s(pt, fsub_s(pt, yb_01_l, yb_11_l), sc, dims_scs, 2, 3, 1); \
    BLOCK_END(3);                                                       \
  }                                                                     \

#endif
