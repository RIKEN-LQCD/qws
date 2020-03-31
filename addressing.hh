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

#include<stddef.h>
#include<inttypes.h>

template<size_t...>
class prod;

// product of first n terms of sequence
// example:
//   prod<3, 2, 1>::p(2) == 1*2   == 2
//   prod<3, 2, 1>::p(3) == 1*2*3 == 6
template<size_t last, size_t... rest>
class prod<last, rest...> {
public:
  static constexpr size_t p(size_t n) {
    return n <= sizeof...(rest) ? prod<rest...>::p(n) : prod<rest...>::p(n) * last;
  }
};

template<>
class prod<> {
public:
  static constexpr size_t p(size_t n) {
    return 1;
  }
};

template<class dims, int64_t...>
class serialize;

// serialize index of multi array
// example:
//   serialized indexes of [0][0], [0][1], [0][2] and [1][0] of x defined as x[2][3] are
//     serialize<prod<2, 3>, 0, 0>::s() == 0,
//     serialize<prod<2, 3>, 0, 1>::s() == 1,
//     serialize<prod<2, 3>, 0, 2>::s() == 2 and
//     serialize<prod<2, 3>, 1, 0>::s() == 3.
template<class dims, int64_t last, int64_t... rest>
class serialize<dims, last, rest...> {
public:
  static constexpr int64_t s() {
    return dims::p(sizeof...(rest)) * last + serialize<dims, rest...>::s();
  }
};

template<class dims>
class serialize<dims> {
public:
  static constexpr int64_t s() {
    return 0;
  }
};

// divide index into base and immediate index
// parameters:
//   imm_min/imm_max: min/max number of the immediate index (inclusive)
//   dims: a type of an instance of prod in which arguments represent the lengths of the multi array
//   p: indexes of the multi array
// example:
//   Let a range of immediate index [-16,15]. then the address x+100 can be accessed as (x+16+2*32) + -12.
//   16+2*32 and -12 are given by addresing<-16, 15, prod<1>, 100>::base and index.
//   When an array is multi dimensional, the arguments are like addressing<-16, 15, prod<2, 3, 4>, 0, 1, 2>.
template <int64_t imm_min, int64_t imm_max, class dims, int64_t... p>
class addressing {
public:
  static constexpr int64_t len_range = imm_max - imm_min + 1;
  static constexpr int64_t spos = serialize<dims, p...>::s();
  static constexpr int64_t base = (spos/len_range)*len_range - imm_min;
  static constexpr int64_t index = spos - base;
};

// preset sizes for qws (The values of leftmost dimensions are meaningless.)
typedef prod<3, 4, 2> dims_scs;
typedef prod<3, 3, 2> dims_glus;
typedef prod<2, 6, 2> dims_x;
typedef prod<2, 36> dims_clv;
typedef prod<3, 2, 2> dims_projscs;
typedef prod<1> dims_1d;
typedef prod<12, 2> dims_cs;
