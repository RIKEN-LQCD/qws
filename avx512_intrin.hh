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

#pragma once

#include <immintrin.h>
#include"qws.h"
#include"addressing.hh"

typedef __m512 vecs_t;
typedef __mmask16 pred_t;

#define pred_true_all() \
  0xffff

#define pred_true_except_lowest() \
  0xfffe

#define pred_true_lowest() \
  0x0001

#define pred_true_except_highest() \
  0x7fff

#define pred_true_highest() \
  0x8000

#define fzero_s() \
  _mm512_setzero_ps()

#define fload1_s(pred, head, dims, ...) \
  _mm512_mask_load_ps(_mm512_setzero_ps(), pred, ((const float*)(head)) + addressing<0, 0, dims, __VA_ARGS__>::base*VLENS)

#define fstore1_s(pred, data, head, dims, ...) \
  _mm512_mask_store_ps(((float*)(head)) + addressing<0, 0, dims, __VA_ARGS__>::base*VLENS, pred, data)

#define fadd_s(pred, x, y) \
  _mm512_mask_add_ps(x, pred, x, y)

#define fsub_s(pred, x, y) \
  _mm512_mask_sub_ps(x, pred, x, y)

#define fadd_s_m(pred, x, y) \
  _mm512_mask_add_ps(x, pred, x, y)

#define fsub_s_m(pred, x, y) \
  _mm512_mask_sub_ps(x, pred, x, y)

#define fmul_s(pred, x, y) \
  _mm512_mask_mul_ps(x, pred, x, y)

// a*b+c
#define fmadd_s(pred, a, b, c) \
  _mm512_mask_fmadd_ps(a, pred, b, c)

// -a*b+c
#define fnmadd_s(pred, a, b, c) \
  _mm512_mask_fnmadd_ps(a, pred, b, c)

// a*b-c;
#define fmsub_s(pred, a, b, c) \
  _mm512_mask_fmsub_ps(a, pred, b, c)

#define fsum_s(pred, a) \
  _mm512_mask_reduce_add_ps(pred, a)

#define fneg_s(pred, x) \
  _mm512_mask_mul_ps(x, pred, x, _mm512_set1_ps(-1.f))

#define for_s(pred, x, y) \
  _mm512_castsi512_ps(_mm512_mask_or_epi32(_mm512_castps_si512(x),pred, _mm512_castps_si512(x), _mm512_castps_si512(y)))

#define fdup_s(x) \
  _mm512_set1_ps(x)
