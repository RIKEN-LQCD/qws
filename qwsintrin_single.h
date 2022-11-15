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
#ifndef QWS_INTRINSICS_SINGLE
#define QWS_INTRINSICS_SINGLE

#include"qws.h"
#include"addressing.hh"

#ifdef ARCH_POSTK
#include "sve_intrin.hh"
#elif defined(ARCH_AVX512)
#include "avx512_intrin.hh"
#else

typedef rvecs_t vecs_t;
struct pred_t {
  char v[VLENS];
};

//=true
inline pred_t pred_true_all(){
  int i;
  pred_t tmp;
#pragma loop nounroll
  for (i=0;i<VLENS;i++) {
    tmp.v[i] = 1;
  }
  return tmp;
}

//={false, true, true, ..., true}
inline pred_t pred_true_except_lowest(){
  int i;
  pred_t tmp = {};
#pragma loop nounroll
  for (i=1;i<VLENS;i++) {
    tmp.v[i] = 1;
  }
  return tmp;
}

//={true, false, false, ..., false}
inline pred_t pred_true_lowest(){
  pred_t tmp = {};
  tmp.v[0] = 1;
  return tmp;
}

//={true, true, ..., true, false}
inline pred_t pred_true_except_highest(){
  int i;
  pred_t tmp = {};
#pragma loop nounroll
  for (i=0;i<VLENS-1;i++) {
    tmp.v[i] = 1;
  }
  return tmp;
}

//={false, false, ..., false, true}
inline pred_t pred_true_highest(){
  pred_t tmp = {};
  tmp.v[VLENS-1] = 1;
  return tmp;
}

//=0;
inline vecs_t fzero_s(void){
  int i;
  vecs_t tmp;
#pragma loop nounroll
  for (i=0;i<VLENS;i++){tmp.v[i] = 0;}
  return tmp;
}

//=a;
inline vecs_t fcopy_s(const vecs_t &a){
  int i;
  vecs_t tmp;
#pragma loop nounroll
  for (i=0;i<VLENS;i++){tmp.v[i] = a.v[i];}
  return tmp;
}
//=-a;
inline vecs_t mfcopy_s(const vecs_t &a){
  int i;
  vecs_t tmp;
#pragma loop nounroll
  for (i=0;i<VLENS;i++){tmp.v[i] =-a.v[i];}
  return tmp;
}

#define fload1_s(pred, head, dims, ...) \
  fload1_s_(pred, ((const float*)(head)) + addressing<0, 0, dims, __VA_ARGS__>::base*VLENS)

//a[v]=a;
inline vecs_t fload1_s_(const pred_t &p, const float *a){
  int i;
  vecs_t tmp;
#pragma loop nounroll
  for (i=0;i<VLENS;i++) {
    if(p.v[i])
      tmp.v[i] = a[i];
    else
      tmp.v[i] = 0;
  }
  return tmp;
}

#define fstore1_s(pred, data, head, dims, ...) \
  fstore1_s_(pred, data, ((float*)(head)) + addressing<0, 0, dims, __VA_ARGS__>::base*VLENS)

//dest[v]=src[v]
inline void fstore1_s_(const pred_t &p, const vecs_t &src, float *dest){
  int i;
#pragma loop nounroll
  for (i=0;i<VLENS;i++) {
    if(p.v[i])
      dest[i] = src.v[i];
  }
}

// a+b;
inline vecs_t fadd_s(const pred_t &p, const vecs_t &a, const vecs_t &b){
  int i;
  vecs_t tmp;
#pragma loop nounroll
  for (i=0;i<VLENS;i++) {
    if(p.v[i])
      tmp.v[i] = a.v[i]+b.v[i];
    else
      tmp.v[i] = 0;
  }
  return tmp;
}

// a-b;
inline vecs_t fsub_s(const pred_t &p, const vecs_t &a, const vecs_t &b){
  int i;
  vecs_t tmp;
#pragma loop nounroll
  for (i=0;i<VLENS;i++) {
    if(p.v[i])
      tmp.v[i] = a.v[i]-b.v[i];
    else
      tmp.v[i] = 0;
  }
  return tmp;
}

// p ? a+b : a;
inline vecs_t fadd_s_m(const pred_t &p, const vecs_t &a, const vecs_t &b){
  int i;
  vecs_t tmp;
#pragma loop nounroll
  for (i=0;i<VLENS;i++) {
    if(p.v[i])
      tmp.v[i] = a.v[i]+b.v[i];
    else
      tmp.v[i] = a.v[i];
  }
  return tmp;
}

// p ? a-b : a;
inline vecs_t fsub_s_m(const pred_t &p, const vecs_t &a, const vecs_t &b){
  int i;
  vecs_t tmp;
#pragma loop nounroll
  for (i=0;i<VLENS;i++) {
    if(p.v[i])
      tmp.v[i] = a.v[i]-b.v[i];
    else
      tmp.v[i] = a.v[i];
  }
  return tmp;
}

// a*b;
inline vecs_t fmul_s(const pred_t &p, const vecs_t &a, const vecs_t &b){
  int i;
  vecs_t tmp;
#pragma loop nounroll
  for (i=0;i<VLENS;i++) {
    if(p.v[i])
      tmp.v[i] = a.v[i]*b.v[i];
    else
      tmp.v[i] = 0;
  }
  return tmp;
}

inline vecs_t fmul_s(const vecs_t &a, const vecs_t &b) {
  return fmul_s(pred_true_all(), a, b);
}

// a*b+c;
inline vecs_t fmadd_s(const pred_t &p, const vecs_t &a, const vecs_t &b,const vecs_t &c){
  int i;
  vecs_t tmp;
#pragma loop nounroll
  for (i=0;i<VLENS;i++) {
    if(p.v[i])
      tmp.v[i] = a.v[i]*b.v[i]+c.v[i];
    else
      tmp.v[i] = 0;
  }
  return tmp;
}

//-a*b+c;
inline vecs_t fnmadd_s(const pred_t &p, const vecs_t &a, const vecs_t &b,const vecs_t &c){
  int i;
  vecs_t tmp;
#pragma loop nounroll
  for (i=0;i<VLENS;i++) {
    if(p.v[i])
      tmp.v[i] =-a.v[i]*b.v[i]+c.v[i];
    else
      tmp.v[i] = 0;
  }
  return tmp;
}

// a*b-c;
inline vecs_t fmsub_s(const pred_t &p, const vecs_t &a, const vecs_t &b,const vecs_t &c){
  int i;
  vecs_t tmp;
#pragma loop nounroll
  for (i=0;i<VLENS;i++) {
    if(p.v[i])
      tmp.v[i] = a.v[i]*b.v[i]-c.v[i];
    else
      tmp.v[i] = 0;
  }
  return tmp;
}

//-a*b-c;
inline vecs_t fnmsub_s(const pred_t &p, const vecs_t &a, const vecs_t &b,const vecs_t &c){
  int i;
  vecs_t tmp;
#pragma loop nounroll
  for (i=0;i<VLENS;i++) {
    if(p.v[i])
      tmp.v[i] =-a.v[i]*b.v[i]-c.v[i];
    else
      tmp.v[i] = 0;
  }
  return tmp;
}

//=sum(a);
inline float fsum_s(const pred_t &p, const vecs_t &a){
  int i;
  float tmp=0;
#pragma loop nounroll
  for (i=0;i<VLENS;i++) {
    if(p.v[i])
      tmp += a.v[i];
  }
  return tmp;
}

//=-a
inline vecs_t fneg_s(const pred_t &p, const vecs_t &a){
  int i;
  vecs_t tmp;
#pragma loop nounroll
  for (i=0;i<VLENS;i++) {
    if(p.v[i])
      tmp.v[i] = -1 * a.v[i];
    else
      tmp.v[i] = 0;
  }
  return tmp;
}

//=a|b;
inline vecs_t for_s(const pred_t &p, const vecs_t &a, const vecs_t &b){
  int i;
  vecs_t tmp;
#pragma loop nounroll
  for (i=0;i<VLENS;i++) {
    if(p.v[i]) {
      int32_t tmp2 = *((int32_t*)&(a.v[i])) | *((int32_t*)&(b.v[i]));
      tmp.v[i] = *((float*)&tmp2);
    } else {
      tmp.v[i] = 0;
    }
  }
  return tmp;
}

//={x, x, ..., x};
inline vecs_t fdup_s(float x){
  int i;
  vecs_t tmp;
#pragma loop nounroll
  for (i=0;i<VLENS;i++) {
    tmp.v[i] = x;
  }
  return tmp;
}

#endif

#endif // QWS_INTRINSICS_SINGLE
