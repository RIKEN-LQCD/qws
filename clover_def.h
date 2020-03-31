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
#ifndef CLOVER_DEF_H
#define CLOVER_DEF_H

#include "qws.h"
#include "qwsintrin.h"

# define FADD(a, b, r) {			\
    for(int __i = 0; __i < VLEND; __i++) {	\
      r.v[__i] = a.v[__i] + b.v[__i];		\
    }						\
  }

# define FSUB(a, b, r) {			\
    for(int __i = 0; __i < VLEND; __i++) {	\
      r.v[__i] = a.v[__i] - b.v[__i];		\
    }						\
  }

# define FMUL(a, b, r) {			\
    for(int __i = 0; __i < VLEND; __i++) {	\
      r.v[__i] = a.v[__i] * b.v[__i];		\
    }						\
  }

# define FMADD(a, b, c, r) {			\
    for(int __i = 0; __i < VLEND; __i++) {	\
      r.v[__i] = a.v[__i] * b.v[__i] + c.v[__i];	\
    }						\
  }

# define FNMADD(a, b, c, r) {			\
    for(int __i = 0; __i < VLEND; __i++) {	\
      r.v[__i] = - a.v[__i] * b.v[__i] + c.v[__i];	\
    }						\
  }



static inline void __mult_clvd_def(rvecd_t sc[3][4][2], rvecd_t a[2][36]) {
  int c, ri, i;
  rvecd_t x[2][6][2];
  rvecd_t y[2][2];

  for (c=0;c<3;c++){
    for (ri=0;ri<2;ri++){
      FADD(sc[c][0][ri], sc[c][2][ri], x[0][0+c][ri]);
      FADD(sc[c][1][ri], sc[c][3][ri], x[0][3+c][ri]);
      FSUB(sc[c][0][ri], sc[c][2][ri], x[1][0+c][ri]);
      FSUB(sc[c][1][ri], sc[c][3][ri], x[1][3+c][ri]);
    }
  }

  for (i=0;i<2;i++){
    FMUL(a[i][0], x[i][0][0], y[i][0]);
    FMUL(a[i][0], x[i][0][1], y[i][1]);

    FMADD( a[i][ 6], x[i][1][0], y[i][0], y[i][0]);
    FMADD( a[i][ 8], x[i][2][0], y[i][0], y[i][0]);
    FMADD( a[i][10], x[i][3][0], y[i][0], y[i][0]);
    FMADD( a[i][12], x[i][4][0], y[i][0], y[i][0]);
    FMADD( a[i][14], x[i][5][0], y[i][0], y[i][0]);
    FNMADD(a[i][ 7], x[i][1][1], y[i][0], y[i][0]);
    FNMADD(a[i][ 9], x[i][2][1], y[i][0], y[i][0]);
    FNMADD(a[i][11], x[i][3][1], y[i][0], y[i][0]);
    FNMADD(a[i][13], x[i][4][1], y[i][0], y[i][0]);
    FNMADD(a[i][15], x[i][5][1], y[i][0], y[i][0]);

    FMADD( a[i][ 6], x[i][1][1], y[i][1], y[i][1]);
    FMADD( a[i][ 8], x[i][2][1], y[i][1], y[i][1]);
    FMADD( a[i][10], x[i][3][1], y[i][1], y[i][1]);
    FMADD( a[i][12], x[i][4][1], y[i][1], y[i][1]);
    FMADD( a[i][14], x[i][5][1], y[i][1], y[i][1]);
    FMADD( a[i][ 7], x[i][1][0], y[i][1], y[i][1]);
    FMADD( a[i][ 9], x[i][2][0], y[i][1], y[i][1]);
    FMADD( a[i][11], x[i][3][0], y[i][1], y[i][1]);
    FMADD( a[i][13], x[i][4][0], y[i][1], y[i][1]);
    FMADD( a[i][15], x[i][5][0], y[i][1], y[i][1]);
  }

  FADD(y[0][0], y[1][0], sc[0][0][0]);
  FADD(y[0][1], y[1][1], sc[0][0][1]);
  FSUB(y[0][0], y[1][0], sc[0][2][0]);
  FSUB(y[0][1], y[1][1], sc[0][2][1]);

  for (i=0;i<2;i++){
    FMUL(a[i][1], x[i][1][0], y[i][0]);
    FMUL(a[i][1], x[i][1][1], y[i][1]);
    FMADD( a[i][ 6], x[i][0][0], y[i][0], y[i][0]);
    FMADD( a[i][16], x[i][2][0], y[i][0], y[i][0]);
    FMADD( a[i][18], x[i][3][0], y[i][0], y[i][0]);
    FMADD( a[i][20], x[i][4][0], y[i][0], y[i][0]);
    FMADD( a[i][22], x[i][5][0], y[i][0], y[i][0]);
    FMADD( a[i][ 7], x[i][0][1], y[i][0], y[i][0]);
    FNMADD(a[i][17], x[i][2][1], y[i][0], y[i][0]);
    FNMADD(a[i][19], x[i][3][1], y[i][0], y[i][0]);
    FNMADD(a[i][21], x[i][4][1], y[i][0], y[i][0]);
    FNMADD(a[i][23], x[i][5][1], y[i][0], y[i][0]);
    FMADD( a[i][ 6], x[i][0][1], y[i][1], y[i][1]);
    FMADD( a[i][16], x[i][2][1], y[i][1], y[i][1]);
    FMADD( a[i][18], x[i][3][1], y[i][1], y[i][1]);
    FMADD( a[i][20], x[i][4][1], y[i][1], y[i][1]);
    FMADD( a[i][22], x[i][5][1], y[i][1], y[i][1]);
    FNMADD(a[i][ 7], x[i][0][0], y[i][1], y[i][1]);
    FMADD( a[i][17], x[i][2][0], y[i][1], y[i][1]);
    FMADD( a[i][19], x[i][3][0], y[i][1], y[i][1]);
    FMADD( a[i][21], x[i][4][0], y[i][1], y[i][1]);
    FMADD( a[i][23], x[i][5][0], y[i][1], y[i][1]);
  }

  FADD(y[0][0], y[1][0], sc[1][0][0]);
  FADD(y[0][1], y[1][1], sc[1][0][1]);
  FSUB(y[0][0], y[1][0], sc[1][2][0]);
  FSUB(y[0][1], y[1][1], sc[1][2][1]);

  for (i=0;i<2;i++){
    FMUL(a[i][2], x[i][2][0], y[i][0]);
    FMUL(a[i][2], x[i][2][1], y[i][1]);
    FMADD( a[i][ 8], x[i][0][0], y[i][0], y[i][0]);
    FMADD( a[i][16], x[i][1][0], y[i][0], y[i][0]);
    FMADD( a[i][24], x[i][3][0], y[i][0], y[i][0]);
    FMADD( a[i][26], x[i][4][0], y[i][0], y[i][0]);
    FMADD( a[i][28], x[i][5][0], y[i][0], y[i][0]);
    FMADD( a[i][ 9], x[i][0][1], y[i][0], y[i][0]);
    FMADD( a[i][17], x[i][1][1], y[i][0], y[i][0]);
    FNMADD(a[i][25], x[i][3][1], y[i][0], y[i][0]);
    FNMADD(a[i][27], x[i][4][1], y[i][0], y[i][0]);
    FNMADD(a[i][29], x[i][5][1], y[i][0], y[i][0]);
    FMADD( a[i][ 8], x[i][0][1], y[i][1], y[i][1]);
    FMADD( a[i][16], x[i][1][1], y[i][1], y[i][1]);
    FMADD( a[i][24], x[i][3][1], y[i][1], y[i][1]);
    FMADD( a[i][26], x[i][4][1], y[i][1], y[i][1]);
    FMADD( a[i][28], x[i][5][1], y[i][1], y[i][1]);
    FNMADD(a[i][ 9], x[i][0][0], y[i][1], y[i][1]);
    FNMADD(a[i][17], x[i][1][0], y[i][1], y[i][1]);
    FMADD( a[i][25], x[i][3][0], y[i][1], y[i][1]);
    FMADD( a[i][27], x[i][4][0], y[i][1], y[i][1]);
    FMADD( a[i][29], x[i][5][0], y[i][1], y[i][1]);
  }

  FADD(y[0][0], y[1][0], sc[2][0][0]);
  FADD(y[0][1], y[1][1], sc[2][0][1]);
  FSUB(y[0][0], y[1][0], sc[2][2][0]);
  FSUB(y[0][1], y[1][1], sc[2][2][1]);

  for (i=0;i<2;i++){
    FMUL(a[i][3], x[i][3][0], y[i][0]);
    FMUL(a[i][3], x[i][3][1], y[i][1]);
    FMADD( a[i][10], x[i][0][0], y[i][0], y[i][0]);
    FMADD( a[i][18], x[i][1][0], y[i][0], y[i][0]);
    FMADD( a[i][24], x[i][2][0], y[i][0], y[i][0]);
    FMADD( a[i][30], x[i][4][0], y[i][0], y[i][0]);
    FMADD( a[i][32], x[i][5][0], y[i][0], y[i][0]);
    FMADD( a[i][11], x[i][0][1], y[i][0], y[i][0]);
    FMADD( a[i][19], x[i][1][1], y[i][0], y[i][0]);
    FMADD( a[i][25], x[i][2][1], y[i][0], y[i][0]);
    FNMADD(a[i][31], x[i][4][1], y[i][0], y[i][0]);
    FNMADD(a[i][33], x[i][5][1], y[i][0], y[i][0]);
    FMADD( a[i][10], x[i][0][1], y[i][1], y[i][1]);
    FMADD( a[i][18], x[i][1][1], y[i][1], y[i][1]);
    FMADD( a[i][24], x[i][2][1], y[i][1], y[i][1]);
    FMADD( a[i][30], x[i][4][1], y[i][1], y[i][1]);
    FMADD( a[i][32], x[i][5][1], y[i][1], y[i][1]);
    FNMADD(a[i][11], x[i][0][0], y[i][1], y[i][1]);
    FNMADD(a[i][19], x[i][1][0], y[i][1], y[i][1]);
    FNMADD(a[i][25], x[i][2][0], y[i][1], y[i][1]);
    FMADD( a[i][31], x[i][4][0], y[i][1], y[i][1]);
    FMADD( a[i][33], x[i][5][0], y[i][1], y[i][1]);
  }

  FADD(y[0][0], y[1][0], sc[0][1][0]);
  FADD(y[0][1], y[1][1], sc[0][1][1]);
  FSUB(y[0][0], y[1][0], sc[0][3][0]);
  FSUB(y[0][1], y[1][1], sc[0][3][1]);

  for (i=0;i<2;i++){
    FMUL(a[i][4], x[i][4][0], y[i][0]);
    FMUL(a[i][4], x[i][4][1], y[i][1]);
    FMADD( a[i][12], x[i][0][0], y[i][0], y[i][0]);
    FMADD( a[i][20], x[i][1][0], y[i][0], y[i][0]);
    FMADD( a[i][26], x[i][2][0], y[i][0], y[i][0]);
    FMADD( a[i][30], x[i][3][0], y[i][0], y[i][0]);
    FMADD( a[i][34], x[i][5][0], y[i][0], y[i][0]);
    FMADD( a[i][13], x[i][0][1], y[i][0], y[i][0]);
    FMADD( a[i][21], x[i][1][1], y[i][0], y[i][0]);
    FMADD( a[i][27], x[i][2][1], y[i][0], y[i][0]);
    FMADD( a[i][31], x[i][3][1], y[i][0], y[i][0]);
    FNMADD(a[i][35], x[i][5][1], y[i][0], y[i][0]);
    FMADD( a[i][12], x[i][0][1], y[i][1], y[i][1]);
    FMADD( a[i][20], x[i][1][1], y[i][1], y[i][1]);
    FMADD( a[i][26], x[i][2][1], y[i][1], y[i][1]);
    FMADD( a[i][30], x[i][3][1], y[i][1], y[i][1]);
    FMADD( a[i][34], x[i][5][1], y[i][1], y[i][1]);
    FNMADD(a[i][13], x[i][0][0], y[i][1], y[i][1]);
    FNMADD(a[i][21], x[i][1][0], y[i][1], y[i][1]);
    FNMADD(a[i][27], x[i][2][0], y[i][1], y[i][1]);
    FNMADD(a[i][31], x[i][3][0], y[i][1], y[i][1]);
    FMADD( a[i][35], x[i][5][0], y[i][1], y[i][1]);
  }

  FADD(y[0][0], y[1][0], sc[1][1][0]);
  FADD(y[0][1], y[1][1], sc[1][1][1]);
  FSUB(y[0][0], y[1][0], sc[1][3][0]);
  FSUB(y[0][1], y[1][1], sc[1][3][1]);

  for (i=0;i<2;i++){
    FMUL(a[i][5], x[i][5][0], y[i][0]);
    FMUL(a[i][5], x[i][5][1], y[i][1]);
    FMADD( a[i][14], x[i][0][0], y[i][0], y[i][0]);
    FMADD( a[i][22], x[i][1][0], y[i][0], y[i][0]);
    FMADD( a[i][28], x[i][2][0], y[i][0], y[i][0]);
    FMADD( a[i][32], x[i][3][0], y[i][0], y[i][0]);
    FMADD( a[i][34], x[i][4][0], y[i][0], y[i][0]);
    FMADD( a[i][15], x[i][0][1], y[i][0], y[i][0]);
    FMADD( a[i][23], x[i][1][1], y[i][0], y[i][0]);
    FMADD( a[i][29], x[i][2][1], y[i][0], y[i][0]);
    FMADD( a[i][33], x[i][3][1], y[i][0], y[i][0]);
    FMADD( a[i][35], x[i][4][1], y[i][0], y[i][0]);
    FMADD( a[i][14], x[i][0][1], y[i][1], y[i][1]);
    FMADD( a[i][22], x[i][1][1], y[i][1], y[i][1]);
    FMADD( a[i][28], x[i][2][1], y[i][1], y[i][1]);
    FMADD( a[i][32], x[i][3][1], y[i][1], y[i][1]);
    FMADD( a[i][34], x[i][4][1], y[i][1], y[i][1]);
    FNMADD(a[i][15], x[i][0][0], y[i][1], y[i][1]);
    FNMADD(a[i][23], x[i][1][0], y[i][1], y[i][1]);
    FNMADD(a[i][29], x[i][2][0], y[i][1], y[i][1]);
    FNMADD(a[i][33], x[i][3][0], y[i][1], y[i][1]);
    FNMADD(a[i][35], x[i][4][0], y[i][1], y[i][1]);
  }

  FADD(y[0][0], y[1][0], sc[2][1][0]);
  FADD(y[0][1], y[1][1], sc[2][1][1]);
  FSUB(y[0][0], y[1][0], sc[2][3][0]);
  FSUB(y[0][1], y[1][1], sc[2][3][1]);

}
#endif
