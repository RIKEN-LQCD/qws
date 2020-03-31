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

// example:
//   -------------------------------------------------
//   #define S(x, y, z) int i_##x##y##z;
//   LOOP_2(LOOP2a, LOOP3, S) // right is inner
//                            // same macro cannot be used twice
//   -------------------------------------------------
//
//   is expanded as
//
//   -------------------------------------------------
//   int i_0_0_0; int i_0_0_1; int i_0_0_2; \
//   int i_0_1_0; int i_0_1_1; int i_0_1_2; \
//   int i_1_0_0; int i_1_0_1; int i_1_0_2; \
//   int i_1_1_0; int i_1_1_1; int i_1_1_2;
//   -------------------------------------------------
#define  LOOP_2(S, ...) S(__VA_ARGS__, 0) S(__VA_ARGS__, 1)
#define LOOP_2a(S, ...) S(__VA_ARGS__, 0) S(__VA_ARGS__, 1)
#define  LOOP_3(S, ...) S(__VA_ARGS__, 0) S(__VA_ARGS__, 1) S(__VA_ARGS__, 2)
#define  LOOP_4(S, ...) S(__VA_ARGS__, 0) S(__VA_ARGS__, 1) S(__VA_ARGS__, 2) S(__VA_ARGS__, 3)
#define  LOOP_6(S, ...) \
  S(__VA_ARGS__, 0) S(__VA_ARGS__, 1) S(__VA_ARGS__, 2) S(__VA_ARGS__, 3) \
  S(__VA_ARGS__, 4) S(__VA_ARGS__, 5)
#define  LOOP_8(S, ...) \
  S(__VA_ARGS__, 0) S(__VA_ARGS__, 1) S(__VA_ARGS__, 2) S(__VA_ARGS__, 3) \
  S(__VA_ARGS__, 4) S(__VA_ARGS__, 5) S(__VA_ARGS__, 6) S(__VA_ARGS__, 7)


#define alloca_aligned(p, size, align)                                  \
  char p##_body[size+align]; p = (decltype(p))((uintptr_t)(p##_body+align) - ((uintptr_t)p##_body)%align);

#ifdef INLINE_ASM_UNAVAILABLE
#ifndef UTIL_C
extern long int true0, true1, true2, true3;
extern long int zero;
#else
long int true0=1, true1=1, true2=1, true3=1;
long int zero=0;
#endif
#define BLOCK_START(x) if(true ## x) {
#define BLOCK_END(x) }
#define ASSUME_MODIFIED(x) \
  x = decltype(x)(((uintptr_t)x) >> zero);
#else
#define BLOCK_START(x) __asm__ volatile ("":::"memory")
#define BLOCK_END(x) 
#define ASSUME_MODIFIED(x) __asm__ volatile ("":"+r"(x))
#endif

#include<tuple>

void dump_asm_data(const int *p, int n, const char *symbol);
void dummy();
std::tuple<long, long> static_sched(long n, long nthrd, long tid);
