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
#ifndef WILSON_DEF_H
#define WILSON_DEF_H
#include "qws.h"

#define __reg_ri(a)				\
  __attribute__((aligned(64))) rvecd_t a ## 0;	\
  __attribute__((aligned(64))) rvecd_t a ## 1;	\

#define __reg_2ri(a)				\
  __reg_ri(a ## 0)				\
  __reg_ri(a ## 1)				\

#define __reg_3ri(a)				\
  __reg_ri(a ## 0)				\
  __reg_ri(a ## 1)				\
  __reg_ri(a ## 2)				\

#define __reg_4ri(a)				\
  __reg_ri(a ## 0)				\
  __reg_ri(a ## 1)				\
  __reg_ri(a ## 2)				\
  __reg_ri(a ## 3)				\
  
#define __reg_32ri(a)				\
  __reg_2ri(a ## 0)				\
  __reg_2ri(a ## 1)				\
  __reg_2ri(a ## 2)				\

#define __reg_33ri(a)				\
  __reg_3ri(a ## 0)				\
  __reg_3ri(a ## 1)				\
  __reg_3ri(a ## 2)				\

#define __reg_34ri(a)				\
  __reg_4ri(a ## 0)				\
  __reg_4ri(a ## 1)				\
  __reg_4ri(a ## 2)				\

#define __load_cs2u(a, s0, s1, in){			\
    a ## 0 ## s0 ## 0 = fcopy( (in).cv[0][0][0] );	\
    a ## 0 ## s0 ## 1 = fcopy( (in).cv[0][0][1] );	\
    a ## 0 ## s1 ## 0 = fcopy( (in).cv[0][1][0] );	\
    a ## 0 ## s1 ## 1 = fcopy( (in).cv[0][1][1] );	\
    a ## 1 ## s0 ## 0 = fcopy( (in).cv[1][0][0] );	\
    a ## 1 ## s0 ## 1 = fcopy( (in).cv[1][0][1] );	\
    a ## 1 ## s1 ## 0 = fcopy( (in).cv[1][1][0] );	\
    a ## 1 ## s1 ## 1 = fcopy( (in).cv[1][1][1] );	\
    a ## 2 ## s0 ## 0 = fcopy( (in).cv[2][0][0] );	\
    a ## 2 ## s0 ## 1 = fcopy( (in).cv[2][0][1] );	\
    a ## 2 ## s1 ## 0 = fcopy( (in).cv[2][1][0] );	\
    a ## 2 ## s1 ## 1 = fcopy( (in).cv[2][1][1] );	\
  }

#define __load_cs2d(a, s0, s1, in){			\
    a ## 0 ## s0 ## 0 = fcopy( (in).cv[0][2][0] );	\
    a ## 0 ## s0 ## 1 = fcopy( (in).cv[0][2][1] );	\
    a ## 0 ## s1 ## 0 = fcopy( (in).cv[0][3][0] );	\
    a ## 0 ## s1 ## 1 = fcopy( (in).cv[0][3][1] );	\
    a ## 1 ## s0 ## 0 = fcopy( (in).cv[1][2][0] );	\
    a ## 1 ## s0 ## 1 = fcopy( (in).cv[1][2][1] );	\
    a ## 1 ## s1 ## 0 = fcopy( (in).cv[1][3][0] );	\
    a ## 1 ## s1 ## 1 = fcopy( (in).cv[1][3][1] );	\
    a ## 2 ## s0 ## 0 = fcopy( (in).cv[2][2][0] );	\
    a ## 2 ## s0 ## 1 = fcopy( (in).cv[2][2][1] );	\
    a ## 2 ## s1 ## 0 = fcopy( (in).cv[2][3][0] );	\
    a ## 2 ## s1 ## 1 = fcopy( (in).cv[2][3][1] );	\
  }

#define __store_cs2u(out, a){			\
    (out).cv[0][0][0] = fcopy( a ## 000 );	\
    (out).cv[0][0][1] = fcopy( a ## 001 );	\
    (out).cv[0][1][0] = fcopy( a ## 010 );	\
    (out).cv[0][1][1] = fcopy( a ## 011 );	\
    (out).cv[1][0][0] = fcopy( a ## 100 );	\
    (out).cv[1][0][1] = fcopy( a ## 101 );	\
    (out).cv[1][1][0] = fcopy( a ## 110 );	\
    (out).cv[1][1][1] = fcopy( a ## 111 );	\
    (out).cv[2][0][0] = fcopy( a ## 200 );	\
    (out).cv[2][0][1] = fcopy( a ## 201 );	\
    (out).cv[2][1][0] = fcopy( a ## 210 );	\
    (out).cv[2][1][1] = fcopy( a ## 211 );	\
  }

#define __store_cs2d(out, a){			\
    (out).cv[0][2][0] = fcopy( a ## 020 );	\
    (out).cv[0][2][1] = fcopy( a ## 021 );	\
    (out).cv[0][3][0] = fcopy( a ## 030 );	\
    (out).cv[0][3][1] = fcopy( a ## 031 );	\
    (out).cv[1][2][0] = fcopy( a ## 120 );	\
    (out).cv[1][2][1] = fcopy( a ## 121 );	\
    (out).cv[1][3][0] = fcopy( a ## 130 );	\
    (out).cv[1][3][1] = fcopy( a ## 131 );	\
    (out).cv[2][2][0] = fcopy( a ## 220 );	\
    (out).cv[2][2][1] = fcopy( a ## 221 );	\
    (out).cv[2][3][0] = fcopy( a ## 230 );	\
    (out).cv[2][3][1] = fcopy( a ## 231 );	\
  }

#define __load_g(g, gp){			\
    g ## 000 = fcopy( (gp).cv[0][0][0] );	\
    g ## 001 = fcopy( (gp).cv[0][0][1] );	\
    g ## 010 = fcopy( (gp).cv[0][1][0] );	\
    g ## 011 = fcopy( (gp).cv[0][1][1] );	\
    g ## 020 = fcopy( (gp).cv[0][2][0] );	\
    g ## 021 = fcopy( (gp).cv[0][2][1] );	\
    g ## 100 = fcopy( (gp).cv[1][0][0] );	\
    g ## 101 = fcopy( (gp).cv[1][0][1] );	\
    g ## 110 = fcopy( (gp).cv[1][1][0] );	\
    g ## 111 = fcopy( (gp).cv[1][1][1] );	\
    g ## 120 = fcopy( (gp).cv[1][2][0] );	\
    g ## 121 = fcopy( (gp).cv[1][2][1] );	\
    g ## 200 = fcopy( (gp).cv[2][0][0] );	\
    g ## 201 = fcopy( (gp).cv[2][0][1] );	\
    g ## 210 = fcopy( (gp).cv[2][1][0] );	\
    g ## 211 = fcopy( (gp).cv[2][1][1] );	\
    g ## 220 = fcopy( (gp).cv[2][2][0] );	\
    g ## 221 = fcopy( (gp).cv[2][2][1] );	\
  }

#define __load_g_dag(g, gp){			\
    g ## 000 = fcopy( (gp).cv[0][0][0] );	\
    g ## 001 =mfcopy( (gp).cv[0][0][1] );	\
    g ## 100 = fcopy( (gp).cv[0][1][0] );	\
    g ## 101 =mfcopy( (gp).cv[0][1][1] );	\
    g ## 200 = fcopy( (gp).cv[0][2][0] );	\
    g ## 201 =mfcopy( (gp).cv[0][2][1] );	\
    g ## 010 = fcopy( (gp).cv[1][0][0] );	\
    g ## 011 =mfcopy( (gp).cv[1][0][1] );	\
    g ## 110 = fcopy( (gp).cv[1][1][0] );	\
    g ## 111 =mfcopy( (gp).cv[1][1][1] );	\
    g ## 210 = fcopy( (gp).cv[1][2][0] );	\
    g ## 211 =mfcopy( (gp).cv[1][2][1] );	\
    g ## 020 = fcopy( (gp).cv[2][0][0] );	\
    g ## 021 =mfcopy( (gp).cv[2][0][1] );	\
    g ## 120 = fcopy( (gp).cv[2][1][0] );	\
    g ## 121 =mfcopy( (gp).cv[2][1][1] );	\
    g ## 220 = fcopy( (gp).cv[2][2][0] );	\
    g ## 221 =mfcopy( (gp).cv[2][2][1] );	\
  }


#define __mult31(b, a, g, i, j)						\
  b ## i ## j ## 0 =   fmul(a ## 0 ## j ## 0, g ## 0 ## i ## 0);	\
  b ## i ## j ## 0 = fnmadd(a ## 0 ## j ## 1, g ## 0 ## i ## 1, b ## i ## j ## 0); \
  b ## i ## j ## 0 =  fmadd(a ## 1 ## j ## 0, g ## 1 ## i ## 0, b ## i ## j ## 0); \
  b ## i ## j ## 0 = fnmadd(a ## 1 ## j ## 1, g ## 1 ## i ## 1, b ## i ## j ## 0); \
  b ## i ## j ## 0 =  fmadd(a ## 2 ## j ## 0, g ## 2 ## i ## 0, b ## i ## j ## 0); \
  b ## i ## j ## 0 = fnmadd(a ## 2 ## j ## 1, g ## 2 ## i ## 1, b ## i ## j ## 0); \
  b ## i ## j ## 1 =   fmul(a ## 0 ## j ## 0, g ## 0 ## i ## 1);	\
  b ## i ## j ## 1 =  fmadd(a ## 0 ## j ## 1, g ## 0 ## i ## 0, b ## i ## j ## 1); \
  b ## i ## j ## 1 =  fmadd(a ## 1 ## j ## 0, g ## 1 ## i ## 1, b ## i ## j ## 1); \
  b ## i ## j ## 1 =  fmadd(a ## 1 ## j ## 1, g ## 1 ## i ## 0, b ## i ## j ## 1); \
  b ## i ## j ## 1 =  fmadd(a ## 2 ## j ## 0, g ## 2 ## i ## 1, b ## i ## j ## 1); \
  b ## i ## j ## 1 =  fmadd(a ## 2 ## j ## 1, g ## 2 ## i ## 0, b ## i ## j ## 1); \

#define __mult_fg(b, a, g)			\
  __mult31(b, a, g, 0, 0)			\
  __mult31(b, a, g, 1, 0)			\
  __mult31(b, a, g, 2, 0)			\
  __mult31(b, a, g, 0, 1)			\
  __mult31(b, a, g, 1, 1)			\
  __mult31(b, a, g, 2, 1)			\



//----------------------------------------------------------------------------------------X forward

//----------------------------------------------------------------------------------------X backward
#define __pxbwd_1col(a, b, col){					\
    a ## col ## 00 = fsub(a ## col ## 00, b ## col ## 11);		\
    a ## col ## 10 = fsub(a ## col ## 10, b ## col ## 01);		\
    a ## col ## 01 = fadd(a ## col ## 01, b ## col ## 10);		\
    a ## col ## 11 = fadd(a ## col ## 11, b ## col ## 00);		\
  }

#define __pxbwd(a, b){			\
    __pxbwd_1col(a, b, 0);		\
    __pxbwd_1col(a, b, 1);		\
    __pxbwd_1col(a, b, 2);		\
  }

#define __xbwd_1col(tmp, col, b){					\
    tmp ## col ## 00 =  fadd(tmp ## col ## 00, b ## col ## 00);		\
    tmp ## col ## 10 =  fadd(tmp ## col ## 10, b ## col ## 10);		\
    tmp ## col ## 20 =  fadd(tmp ## col ## 20, b ## col ## 11);		\
    tmp ## col ## 30 =  fadd(tmp ## col ## 30, b ## col ## 01);		\
    tmp ## col ## 01 =  fadd(tmp ## col ## 01, b ## col ## 01);		\
    tmp ## col ## 11 =  fadd(tmp ## col ## 11, b ## col ## 11);		\
    tmp ## col ## 21 =  fsub(tmp ## col ## 21, b ## col ## 10);		\
    tmp ## col ## 31 =  fsub(tmp ## col ## 31, b ## col ## 00);		\
  }

#define __xbwd(tmp, b){				\
    __xbwd_1col(tmp, 0, b);			\
    __xbwd_1col(tmp, 1, b);			\
    __xbwd_1col(tmp, 2, b);			\
  }

//----------------------------------------------------------------------------------------Y forward
#define __pyfwd_1col(a, b, col){					\
    a ## col ## 00 = fsub(a ## col ## 00, b ## col ## 10);		\
    a ## col ## 10 = fadd(a ## col ## 10, b ## col ## 00);		\
    a ## col ## 01 = fsub(a ## col ## 01, b ## col ## 11);		\
    a ## col ## 11 = fadd(a ## col ## 11, b ## col ## 01);		\
  }

#define __pyfwd(a, b){			\
    __pyfwd_1col(a, b, 0);		\
    __pyfwd_1col(a, b, 1);		\
    __pyfwd_1col(a, b, 2);		\
  }

#define __yfwd_1col(tmp, col, b){					\
    tmp ## col ## 00 =  fadd(tmp ## col ## 00, b ## col ## 00);		\
    tmp ## col ## 10 =  fadd(tmp ## col ## 10, b ## col ## 10);		\
    tmp ## col ## 20 =  fadd(tmp ## col ## 20, b ## col ## 10);		\
    tmp ## col ## 30 =  fsub(tmp ## col ## 30, b ## col ## 00);		\
    tmp ## col ## 01 =  fadd(tmp ## col ## 01, b ## col ## 01);		\
    tmp ## col ## 11 =  fadd(tmp ## col ## 11, b ## col ## 11);		\
    tmp ## col ## 21 =  fadd(tmp ## col ## 21, b ## col ## 11);		\
    tmp ## col ## 31 =  fsub(tmp ## col ## 31, b ## col ## 01);		\
  }

#define __yfwd(tmp, b){				\
    __yfwd_1col(tmp, 0, b);			\
    __yfwd_1col(tmp, 1, b);			\
    __yfwd_1col(tmp, 2, b);			\
  }

//----------------------------------------------------------------------------------------Y backward
#define __pybwd_1col(a, b, col){					\
    a ## col ## 00 = fadd(a ## col ## 00, b ## col ## 10);		\
    a ## col ## 10 = fsub(a ## col ## 10, b ## col ## 00);		\
    a ## col ## 01 = fadd(a ## col ## 01, b ## col ## 11);		\
    a ## col ## 11 = fsub(a ## col ## 11, b ## col ## 01);		\
  }

#define __pybwd(a, b){			\
    __pybwd_1col(a, b, 0);		\
    __pybwd_1col(a, b, 1);		\
    __pybwd_1col(a, b, 2);		\
  }

#define __ybwd_1col(tmp, col, b){					\
    tmp ## col ## 00 =  fadd(tmp ## col ## 00, b ## col ## 00);		\
    tmp ## col ## 10 =  fadd(tmp ## col ## 10, b ## col ## 10);		\
    tmp ## col ## 20 =  fsub(tmp ## col ## 20, b ## col ## 10);		\
    tmp ## col ## 30 =  fadd(tmp ## col ## 30, b ## col ## 00);		\
    tmp ## col ## 01 =  fadd(tmp ## col ## 01, b ## col ## 01);		\
    tmp ## col ## 11 =  fadd(tmp ## col ## 11, b ## col ## 11);		\
    tmp ## col ## 21 =  fsub(tmp ## col ## 21, b ## col ## 11);		\
    tmp ## col ## 31 =  fadd(tmp ## col ## 31, b ## col ## 01);		\
  }

#define __ybwd(tmp, b){				\
    __ybwd_1col(tmp, 0, b);			\
    __ybwd_1col(tmp, 1, b);			\
    __ybwd_1col(tmp, 2, b);			\
  }


//----------------------------------------------------------------------------------------Z forward
#define __pzfwd_1col(a, b, col){		\
    a ## col ## 00 = fadd(a ## col ## 00, b ## col ## 01);			\
    a ## col ## 10 = fsub(a ## col ## 10, b ## col ## 11);			\
    a ## col ## 01 = fsub(a ## col ## 01, b ## col ## 00);			\
    a ## col ## 11 = fadd(a ## col ## 11, b ## col ## 10);			\
  }

#define __pzfwd(a, b){			\
    __pzfwd_1col(a, b, 0);		\
    __pzfwd_1col(a, b, 1);		\
    __pzfwd_1col(a, b, 2);		\
  }

#define __zfwd_1col(tmp, col, b){					\
    tmp ## col ## 00 =  fadd(tmp ## col ## 00, b ## col ## 00);		\
    tmp ## col ## 10 =  fadd(tmp ## col ## 10, b ## col ## 10);		\
    tmp ## col ## 20 =  fsub(tmp ## col ## 20, b ## col ## 01);		\
    tmp ## col ## 30 =  fadd(tmp ## col ## 30, b ## col ## 11);		\
    tmp ## col ## 01 =  fadd(tmp ## col ## 01, b ## col ## 01);		\
    tmp ## col ## 11 =  fadd(tmp ## col ## 11, b ## col ## 11);		\
    tmp ## col ## 21 =  fadd(tmp ## col ## 21, b ## col ## 00);		\
    tmp ## col ## 31 =  fsub(tmp ## col ## 31, b ## col ## 10);		\
  }

#define __zfwd(tmp, b){				\
    __zfwd_1col(tmp, 0, b);			\
    __zfwd_1col(tmp, 1, b);			\
    __zfwd_1col(tmp, 2, b);			\
  }

//----------------------------------------------------------------------------------------Z backward
#define __pzbwd_1col(a, b, col){		\
    a ## col ## 00 = fsub(a ## col ## 00, b ## col ## 01);			\
    a ## col ## 10 = fadd(a ## col ## 10, b ## col ## 11);			\
    a ## col ## 01 = fadd(a ## col ## 01, b ## col ## 00);			\
    a ## col ## 11 = fsub(a ## col ## 11, b ## col ## 10);			\
  }

#define __pzbwd(a, b){			\
    __pzbwd_1col(a, b, 0);		\
    __pzbwd_1col(a, b, 1);		\
    __pzbwd_1col(a, b, 2);		\
  }

#define __zbwd_1col(tmp, col, b){					\
    tmp ## col ## 00 =  fadd(tmp ## col ## 00, b ## col ## 00);		\
    tmp ## col ## 10 =  fadd(tmp ## col ## 10, b ## col ## 10);		\
    tmp ## col ## 20 =  fadd(tmp ## col ## 20, b ## col ## 01);		\
    tmp ## col ## 30 =  fsub(tmp ## col ## 30, b ## col ## 11);		\
    tmp ## col ## 01 =  fadd(tmp ## col ## 01, b ## col ## 01);		\
    tmp ## col ## 11 =  fadd(tmp ## col ## 11, b ## col ## 11);		\
    tmp ## col ## 21 =  fsub(tmp ## col ## 21, b ## col ## 00);		\
    tmp ## col ## 31 =  fadd(tmp ## col ## 31, b ## col ## 10);		\
  }

#define __zbwd(tmp, b){				\
    __zbwd_1col(tmp, 0, b);			\
    __zbwd_1col(tmp, 1, b);			\
    __zbwd_1col(tmp, 2, b);			\
  }


//----------------------------------------------------------------------------------------T forward
#define __tfwd_1col(tmp, col, b){					\
    tmp ## col ## 20 = fmadd(fac, b ## col ## 00, tmp ## col ## 20);	\
    tmp ## col ## 21 = fmadd(fac, b ## col ## 01, tmp ## col ## 21);	\
    tmp ## col ## 30 = fmadd(fac, b ## col ## 10, tmp ## col ## 30);	\
    tmp ## col ## 31 = fmadd(fac, b ## col ## 11, tmp ## col ## 31);	\
  }

#define __tfwd(tmp, b){				\
    __tfwd_1col(tmp, 0, b);			\
    __tfwd_1col(tmp, 1, b);			\
    __tfwd_1col(tmp, 2, b);			\
  }

//----------------------------------------------------------------------------------------T backward
#define __tbwd_1col(tmp, col, b){					\
    tmp ## col ## 00 = fmadd(fac, b ## col ## 00, tmp ## col ## 00);	\
    tmp ## col ## 01 = fmadd(fac, b ## col ## 01, tmp ## col ## 01);	\
    tmp ## col ## 10 = fmadd(fac, b ## col ## 10, tmp ## col ## 10);	\
    tmp ## col ## 11 = fmadd(fac, b ## col ## 11, tmp ## col ## 11);	\
  }

#define __tbwd(tmp, b){				\
    __tbwd_1col(tmp, 0, b);			\
    __tbwd_1col(tmp, 1, b);			\
    __tbwd_1col(tmp, 2, b);			\
  }



#endif
