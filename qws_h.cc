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
//#define libhalf
#include "qws.h"

#include "mult_all.h"

#ifdef libhalf
#include "half.hpp"
using namespace half_float;
#else
#include <math.h>
typedef __fp16 half;
#endif

#ifdef __cplusplus
extern "C"{
#endif

// Half precision
typedef struct{
  half  v[VLENS];
} rvech_t;

typedef union {
  half   c[3][3][2][VLENS];
  rvech_t  cv[3][3][2];
} g33h_t, *pg33h_t;

typedef union {
  half   c[3][3][2][VLENS];
  half   c_prefetch[18*VLENS];
} gluh_t, *pgluh_t;

typedef union {
  half   c[2][36][VLENS];
  half   c_prefetch[2*36*VLENS];
  rvech_t  cv[2][36];
} clvh_t, *pclvh_t;

union sch_t {
  half   c[3][4][2][VLENS];
  half   c_prefetch[24*VLENS];
  rvech_t  cv[3][4][2];
  rvech_t  cs[12][2];
  rvech_t ccs[24];
  sch_t() {};
};

union projsch_t {
  half   c[3][2][2][VLENS];
  rvech_t  cv[3][2][2];
  projsch_t() {};
};

union projsch1_t {
  half   c[3][2][2];
  projsch1_t() {};
};


#include "clover_h.h"

//#include "prefetch.h"
#include "timing.h"

///////////////////////////////////////////////
// QWS half precision initialization status
///////////////////////////////////////////////
static int m_is_initialised_qws_hf = 0;

////////////////////////////////////////////
// QWS half precision gauge/clover fields
////////////////////////////////////////////
__attribute__((aligned(64))) pgluh_t gluh = nullptr;  // gauge  link (pointer)
__attribute__((aligned(64))) pclvh_t clvh = nullptr;  // clover term (pointer)

/////////////////////////////////
// QWS global variables
/////////////////////////////////
extern double kappa, mkappa;
extern int nx, ny, nz, nt, nxh, nxs;
extern int vols;
extern int domain_e, domain_o;
extern int rank, size, px, py, pz, pt;// process coordinate
extern int pxf, pyf, pzf, ptf;
extern int pxb, pyb, pzb, ptb;
extern int npe[4];
extern double fbc[4][2];

/////////////////////////
// private prototypes
/////////////////////////
void ddd_in_h_(sch_t * __restrict__ out, const sch_t * __restrict__ in, const int* DEO);
void ddd_out_pre_h_(const sch_t * __restrict__  in, const int *idomain);
void ddd_out_pos_h_(sch_t * __restrict__ out, const sch_t * __restrict__ in, const int *idomain, float factor);
void jinv_ddd_in_h_(sch_t * __restrict__   x, const sch_t * __restrict__  b, const int *DEO, const int *maxiter);
void prec_ddd_h_(sch_t *out, const sch_t *in, const int *nsap, const int *nm);
void xbound_h_init_();
void assign_u_s2h(gluh_t *out, const glus_t *in);
void assign_clv_s2h(clvh_t *out, const clvs_t *in);

void qws_h_init_(const glus_t *glus, const clvs_t *clvs)
{

  if ( 0 == m_is_initialised_qws_hf ) {

    if ( nullptr == gluh ) gluh = (pgluh_t)malloc( sizeof(gluh_t)*vols*NDIM*NEO);
    if ( nullptr == clvh ) clvh = (pclvh_t)malloc( sizeof(clvh_t)*vols*NEO);

    xbound_h_init_();

    m_is_initialised_qws_hf = 1;
  }

  assign_u_s2h(   gluh, glus);
  assign_clv_s2h( clvh, clvs);

}

static inline void qws_allreduce(float *array, const int count)
{
#ifdef _MPI_
  MPI_Allreduce(MPI_IN_PLACE,(void *)array,count,MPI_REAL,MPI_SUM,MPI_COMM_WORLD);
#endif
}

void assign_q_s2h(sch_t *out, const scs_t *in, float *norm)
//
// type conversion and data assignment for quark field with normalized norm
//
// out := half(in/norm)
//
//  norm = |in|
//
//   out : half type
//    in : float type
//  norm : norm of in
//
{

  float rtmp = 0.0f;
#pragma omp parallel for reduction(+:rtmp)
  for(int i=0; i<vols*2; i++)
    for(int j=0; j<24; j++)
      for(int v=0; v<VLENS; v++)
        rtmp += in[i].ccs[j].v[v] * in[i].ccs[j].v[v];

  qws_allreduce(&rtmp,1);
//rtmp = 1.0f;
  *norm = sqrtf(rtmp);
  rtmp = 1.0f/(*norm);

#pragma omp parallel for
  for(int i=0; i<vols*2; i++)
    for(int j=0; j<24; j++)
      for(int v=0; v<VLENS; v++)
        out[i].ccs[j].v[v] = half(in[i].ccs[j].v[v]*rtmp);

}


void assign_q_h2s(scs_t *out, const sch_t *in, const float *norm)
//
// type conversion and data assignment for quark field and multiplying normaliztion factor
//
// out := float(in) * norm
//
//  *out : float type
//   *in : half type
//  norm : normalization factor
//
{

#pragma omp parallel for
  for(int i=0; i<vols*2; i++)
    for(int j=0; j<24; j++)
      for(int v=0; v<VLENS; v++)
        out[i].ccs[j].v[v] = float(in[i].ccs[j].v[v]) * (*norm);

}

void assign_u_s2h(gluh_t *out, const glus_t *in)
//
// type conversion and data assignment for link field
//
// out := float(in)
//
//  *out : half type
//   *in : float type
//
{

#pragma omp parallel
  for(int i = 0; i < vols*2*NDIM; ++i) 
    for (int jc = 0; jc < 3; jc++) 
    for (int ic = 0; ic < 3; ic++) 
      for (int ri = 0; ri < 2; ri++) 
        for (int j = 0; j<VLENS; j++) 
          out[i].c[jc][ic][ri][j] = half(in[i].c[jc][ic][ri][j]);

}

void assign_clv_s2h(clvh_t *out, const clvs_t *in)
//
// type conversion and data assignment for clover term
//
// out := float(in)
//
//  *out : float type
//   *in : half type
//
{

#pragma omp parallel
  for(int i = 0; i < vols*2; ++i) 
    for (int ri = 0; ri < 2; ri++)
      for (int ics = 0; ics < 36; ics++)
        for (int j = 0; j<VLENS; j++)
          out[i].c[ri][ics][j] = half(in[i].c[ri][ics][j]);

}


#include "xbound_h.h"
#include "ddd_in_h.h"
#include "ddd_out_h.h"
#include "jinv_ddd_in_h.h"
#include "prec_ddd_s_h.h"

#ifdef __cplusplus
}
#endif
