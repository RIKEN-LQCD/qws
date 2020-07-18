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
// FUJITSU CONFIDENTIAL info for target apps for priority issues(Detailed Design(4))

#pragma once

#include <arm_sve.h>
#include"qws.h"
#include"addressing.hh"

typedef svfloat32_t vecs_t;
typedef svbool_t pred_t;

#define fload1_s(pred, head, dims, ...) \
  svld1_vnum_f32(pred, \
		 (((const float*)(head)) + addressing<-8, 7, dims, __VA_ARGS__>::base*VLENS), \
		 (addressing<-8, 7, dims, __VA_ARGS__>::index))

#define fstore1_s(pred, data, head, dims, ...) \
  svst1_vnum_f32(pred, \
		 (((float*)(head)) + addressing<-8, 7, dims, __VA_ARGS__>::base*VLENS), \
		 (addressing<-8, 7, dims, __VA_ARGS__>::index),		\
                 data)

#define fmul_s(pred, x, y) \
  svmul_f32_x(pred, x, y)

#define fadd_s(pred, x, y) \
  svadd_f32_x(pred, x, y)

#define fsub_s(pred, x, y) \
  svsub_f32_x(pred, x, y)

#define fadd_s_m(pred, x, y) \
  svadd_f32_m(pred, x, y)

#define fsub_s_m(pred, x, y) \
  svsub_f32_m(pred, x, y)

// a*b+c
#define fmadd_s(pred, a, b, c) \
  svmla_f32_x(pred, c, a, b)

// -a*b+c
#define fnmadd_s(pred, a, b, c) \
  svmls_f32_x(pred, c, a, b)

// a*b-c;
#define fmsub_s(pred, a, b, c) \
  svnmsb_f32_x(pred, a, b, c)

#define fsum_s(pred, a) \
  svaddv_f32(pred, a)

#define fzero_s() \
  svdup_n_f32(0)

#define fneg_s(pred, x) \
  svneg_f32_x(pred, x)

#define pred_true_all() \
  svptrue_b32()

#define pred_true_except_lowest() \
  svnot_z(svptrue_b32(), svptrue_pat_b32(SV_VL1))

#define pred_true_lowest() \
  svptrue_pat_b32(SV_VL1)

#define pred_true_except_highest() \
  svwhilelt_b32(1, VLENS)

#define pred_true_highest() \
  svnot_z(svptrue_b32(), pred_true_except_highest())

#define for_s(pred, x, y) \
  svreinterpret_f32_s32(svorr_s32_x(pred, svreinterpret_s32_f32(x), svreinterpret_s32_f32(y)))

#define fdup_s(x) \
  svdup_n_f32(x)
