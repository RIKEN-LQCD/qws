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
#ifndef QWS_H
#define QWS_H

#ifndef VLEND
#define VLEND 8
#endif
#ifndef VLENS
#define VLENS 16
#endif

#define SU3_RECONSTRUCT_D 18
#define SU3_RECONSTRUCT_S 18

#define NEO 2
#define NDIM 4
#define NCOL 3

#ifndef CLS
#define CLS 256
#endif

#ifdef _MPI_
#include <mpi.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef NX
#define NX 32
#endif

#define NXS (NX/2/VLENS)

#ifndef NY
#define NY 4
#endif

#ifndef NZ
#define NZ 4
#endif

#ifndef NT
#define NT 4
#endif

#define VOLS (NXS*NY*NZ*NT)

#ifdef _RDMA_UTOFU
#ifndef _RDMA_
#define _RDMA_
#endif
#endif

#ifdef _RDMA_
#ifndef _MPI_
#error _RDMA_ is defined but _MPI_ is not. 
#endif
#endif

// double precision
typedef struct{
  double v[VLEND];
} rvecd_t;

typedef union {
  double   c[3][3][2][VLEND];
  rvecd_t cv[3][3][2];
} g33d_t, *pg33d_t;

typedef union {
#if SU3_RECONSTRUCT_D == 18
  double   c[3][3][2][VLEND];
#elif SU3_RECONSTRUCT_D == 12
  double   c[2][3][2][VLEND];
#endif
} glud_t, *pglud_t;

typedef union {
  double   c[2][36][VLEND];
  rvecd_t cv[2][36];
} clvd_t, *pclvd_t;

typedef union {
  double   c[3][4][2][VLEND];
  rvecd_t cv[3][4][2];
  rvecd_t cs[12][2];
  rvecd_t ccs[24];
} scd_t;

typedef union {
  double   c[3][2][2][VLEND];
  rvecd_t cv[3][2][2];
} projscd_t;

typedef union {
  double   c[3][2][2];
} projscd1_t;

// single precision
typedef struct{
  float  v[VLENS];
} rvecs_t;

typedef union {
  float   c[3][3][2][VLENS];
  rvecs_t cv[3][3][2];
} g33s_t, *pg33s_t;

typedef union {
#if SU3_RECONSTRUCT_S == 18
  float   c[3][3][2][VLENS];
  float   c_prefetch[18*VLENS];
#elif SU3_RECONSTRUCT_S == 12
  float   c[2][3][2][VLENS];
  float   c_prefetch[12*VLENS];
#endif
} glus_t, *pglus_t;

typedef union {
  float   c[2][36][VLENS];
  float   c_prefetch[2*36*VLENS];
  rvecs_t cv[2][36];
} clvs_t, *pclvs_t;

typedef union {
  float   c[3][4][2][VLENS];
  float   c_prefetch[24*VLENS];
  rvecs_t cv[3][4][2];
  rvecs_t cs[12][2];
  rvecs_t ccs[24];
} scs_t;

typedef union {
  float   c[3][2][2][VLENS];
  float   c_prefetch[12*VLENS];
  rvecs_t cv[3][2][2];
} projscs_t;

typedef union {
  float   c[3][2][2];
} projscs1_t;

// omp block
typedef struct{
  int sy, sz, st;
  int ey, ez, et;
} block_map_t;


#include "qws_xbound.h"

#endif
