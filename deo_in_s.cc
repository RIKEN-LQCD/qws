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
#include "qws.h"
#include "wilson_s.h"
#include "clover_s.h"

#ifdef __cplusplus
extern "C"{
#endif


  extern __attribute__((aligned(64))) pglus_t glus;
  extern __attribute__((aligned(64))) pclvs_t clvs;
  //  nh=nx/2, nd=nx/2/8, nxs=nx/2/16
  extern int nxs, ny, nz, nt, px, py, pz, pt;
  extern int vols;
  //---------------------------------------------------------------------------------------- mult D in bulk
  void deo_in_s_(int* pe, int* po, scs_t* out, scs_t* in){
    __attribute__((aligned(64)))  scs_t tmp;
    int x, y, z, t;
    __attribute__((aligned(64))) g33s_t g;
    __attribute__((aligned(64))) projscs_t a, b;
    __attribute__((aligned(64))) float v[4][2][VLENS];

    pglus_t gxf = &glus[vols*0 + NDIM*vols*(*pe)];
    pglus_t gxb = &glus[vols*0 + NDIM*vols*(*po)];
    pglus_t gyf = &glus[vols*1 + NDIM*vols*(*pe)];
    pglus_t gyb = &glus[vols*1 + NDIM*vols*(*po)];
    pglus_t gzf = &glus[vols*2 + NDIM*vols*(*pe)];
    pglus_t gzb = &glus[vols*2 + NDIM*vols*(*po)];
    pglus_t gtf = &glus[vols*3 + NDIM*vols*(*pe)];
    pglus_t gtb = &glus[vols*3 + NDIM*vols*(*po)];
    int jjj= ((*pe) + ny*py + nz*pz + nt*pt)%2;
#pragma omp for private(t,z,y,x,tmp, g,a,b,v) collapse(3) schedule(static)
    for(t=0; t<nt; t++){
      for(z=0; z<nz; z++){
	for(y=0; y<ny; y++){
	  int ick= ( jjj + y     + z + t)%2;
	  for(x=0; x<nxs; x++){
	    __zero_sc_s(tmp.c);
	    int i0  =      x        + nxs*y + nxs*ny*z + nxs*ny*nz*t;
	    int ixf = (    x+1)%nxs + nxs*y + nxs*ny*z + nxs*ny*nz*t;
	    int ixb = (nxs+x-1)%nxs + nxs*y + nxs*ny*z + nxs*ny*nz*t;
	    int iyf = x + (y+1)*nxs + nxs*ny*z + nxs*ny*nz*t;
	    int iyb = x + (y-1)*nxs + nxs*ny*z + nxs*ny*nz*t;
	    int izf = x + nxs*y + (z+1)*nxs*ny + nxs*ny*nz*t;
	    int izb = x + nxs*y + (z-1)*nxs*ny + nxs*ny*nz*t;
	    int itf = x + nxs*y + nxs*ny*z + (t+1)*nxs*ny*nz;
	    int itb = x + nxs*y + nxs*ny*z + (t-1)*nxs*ny*nz;
	    if (ick == 0) {
              __mult_xfwd0_s(tmp.c, gxf, in, i0);
              if(x != 0 ) {
                __mult_xbwd0_s(tmp.c, gxb, in, i0, ixb);
              } else {
#if SU3_RECONSTRUCT_S == 18
                for(int c2=0;c2<3;c2++){
                  for(int c1=0;c1<3;c1++){
                    g.c[c1][c2][0][0] =0;
                    g.c[c1][c2][1][0] =0;
                    for(int j=1;j<VLENS;j++){
                      g.c[c1][c2][0][j] = (*(gxb+i0)).c[c2][c1][0][j-1];
                      g.c[c1][c2][1][j] =-(*(gxb+i0)).c[c2][c1][1][j-1];
                    }
                  }
                }
#elif SU3_RECONSTRUCT_S == 12
                for(int c2=0;c2<2;c2++){
                  for(int c1=0;c1<3;c1++){
                    g.c[c1][c2][0][0] =0;
                    g.c[c1][c2][1][0] =0;
                    for(int j=1;j<VLENS;j++){
                      g.c[c1][c2][0][j] = (*(gxb+i0)).c[c2][c1][0][j-1];
                      g.c[c1][c2][1][j] =-(*(gxb+i0)).c[c2][c1][1][j-1];
                    }
                  }
                }
                for(int j=0;j<VLENS;j++){
                  g.c[0][2][0][j] = g.c[1][0][0][j]*g.c[2][1][0][j]-g.c[1][0][1][j]*g.c[2][1][1][j]-g.c[2][0][0][j]*g.c[1][1][0][j]+g.c[2][0][1][j]*g.c[1][1][1][j];
                  g.c[1][2][0][j] = g.c[2][0][0][j]*g.c[0][1][0][j]-g.c[2][0][1][j]*g.c[0][1][1][j]-g.c[0][0][0][j]*g.c[2][1][0][j]+g.c[0][0][1][j]*g.c[2][1][1][j];
                  g.c[2][2][0][j] = g.c[0][0][0][j]*g.c[1][1][0][j]-g.c[0][0][1][j]*g.c[1][1][1][j]-g.c[1][0][0][j]*g.c[0][1][0][j]+g.c[1][0][1][j]*g.c[0][1][1][j];
                  g.c[0][2][1][j] =-g.c[1][0][0][j]*g.c[2][1][1][j]-g.c[1][0][1][j]*g.c[2][1][0][j]+g.c[2][0][0][j]*g.c[1][1][1][j]+g.c[2][0][1][j]*g.c[1][1][0][j];
                  g.c[1][2][1][j] =-g.c[2][0][0][j]*g.c[0][1][1][j]-g.c[2][0][1][j]*g.c[0][1][0][j]+g.c[0][0][0][j]*g.c[2][1][1][j]+g.c[0][0][1][j]*g.c[2][1][0][j];
                  g.c[2][2][1][j] =-g.c[0][0][0][j]*g.c[1][1][1][j]-g.c[0][0][1][j]*g.c[1][1][0][j]+g.c[1][0][0][j]*g.c[0][1][1][j]+g.c[1][0][1][j]*g.c[0][1][0][j];
                }
#endif
                for(int c=0;c<3;c++){
                  for (int s=0;s<4;s++) {
                    for (int ri=0;ri<2;ri++) {
                      v[s][ri][0]=0; 
                      for(int j=1;j<VLENS;j++){
                        v[s][ri][j]=in[i0].c[c][s][ri][j-1]; 
                      }
                    }
                  }
                  for(int j=0;j<VLENS;j++){
                    a.c[c][0][0][j] = v[0][0][j] - v[3][1][j];
                    a.c[c][1][0][j] = v[1][0][j] - v[2][1][j];
                    a.c[c][0][1][j] = v[0][1][j] + v[3][0][j];
                    a.c[c][1][1][j] = v[1][1][j] + v[2][0][j];
                  }
                }
                __m_fg_s(b.c, a.c, g.c);
                __s_xb_s(tmp.c, b.c);
              }
            } else {
              if(x != nxs-1 ) {
                __mult_xfwd1_s(tmp.c, gxf, in, i0, ixf);
              } else {
                for(int c=0;c<3;c++){
                  for (int s=0;s<4;s++) {
                    for (int ri=0;ri<2;ri++) {
                      for(int j=0;j<VLENS-1;j++){
                        v[s][ri][j]=in[i0].c[c][s][ri][j+1]; 
                      }
                      v[s][ri][VLENS-1]=0; 
                    }
                  }
                  for(int j=0;j<VLENS;j++){
                    a.c[c][0][0][j] = v[0][0][j] + v[3][1][j];
                    a.c[c][1][0][j] = v[1][0][j] + v[2][1][j];
                    a.c[c][0][1][j] = v[0][1][j] - v[3][0][j];
                    a.c[c][1][1][j] = v[1][1][j] - v[2][0][j];
                  }
                }
                __load_su3_s(g.c, gxf, i0);
                __m_fg_s(b.c, a.c, g.c);
                __s_xf_s(tmp.c, b.c);
              }
              __mult_xbwd1_s(tmp.c, gxb, in, i0);
            }
	    if(y != ny-1 ) __mult_yfwd_s(tmp.c, gyf, in, i0, iyf);
	    if(y != 0    ) __mult_ybwd_s(tmp.c, gyb, in,     iyb);
	    if(z != nz-1 ) __mult_zfwd_s(tmp.c, gzf, in, i0, izf);
	    if(z != 0    ) __mult_zbwd_s(tmp.c, gzb, in,     izb);
	    if(t != nt-1 ) __mult_tfwd_s(tmp.c, gtf, in, i0, itf);
	    if(t != 0    ) __mult_tbwd_s(tmp.c, gtb, in,     itb);
	    __store_sc_s(out[i0].c, tmp.c);
	  }
	}
      }
    }
  }
  //---------------------------------------------------------------------------------------- mult D in bulk
  void dee_deo_in_s_(int* pe, int* po, scs_t* out, scs_t* in){
    __attribute__((aligned(64)))  scs_t tmp;
    int x, y, z, t;
    __attribute__((aligned(64))) g33s_t g;
    __attribute__((aligned(64))) projscs_t a, b;
    __attribute__((aligned(64))) float v[4][2][VLENS];

    pglus_t gxf = &glus[vols*0 + NDIM*vols*(*pe)];
    pglus_t gxb = &glus[vols*0 + NDIM*vols*(*po)];
    pglus_t gyf = &glus[vols*1 + NDIM*vols*(*pe)];
    pglus_t gyb = &glus[vols*1 + NDIM*vols*(*po)];
    pglus_t gzf = &glus[vols*2 + NDIM*vols*(*pe)];
    pglus_t gzb = &glus[vols*2 + NDIM*vols*(*po)];
    pglus_t gtf = &glus[vols*3 + NDIM*vols*(*pe)];
    pglus_t gtb = &glus[vols*3 + NDIM*vols*(*po)];
    int jjj= ((*pe) + ny*py + nz*pz + nt*pt)%2;
#pragma omp for private(t,z,y,x,tmp, g,a,b,v) collapse(3) schedule(static)
    for(t=0; t<nt; t++){
      for(z=0; z<nz; z++){
	for(y=0; y<ny; y++){
	  int ick= ( jjj + y     + z + t)%2;
	  for(x=0; x<nxs; x++){
	    __zero_sc_s(tmp.c);
	    int i0  =      x        + nxs*y + nxs*ny*z + nxs*ny*nz*t;
	    int ixf = (    x+1)%nxs + nxs*y + nxs*ny*z + nxs*ny*nz*t;
	    int ixb = (nxs+x-1)%nxs + nxs*y + nxs*ny*z + nxs*ny*nz*t;
	    int iyf = x + (y+1)*nxs + nxs*ny*z + nxs*ny*nz*t;
	    int iyb = x + (y-1)*nxs + nxs*ny*z + nxs*ny*nz*t;
	    int izf = x + nxs*y + (z+1)*nxs*ny + nxs*ny*nz*t;
	    int izb = x + nxs*y + (z-1)*nxs*ny + nxs*ny*nz*t;
	    int itf = x + nxs*y + nxs*ny*z + (t+1)*nxs*ny*nz;
	    int itb = x + nxs*y + nxs*ny*z + (t-1)*nxs*ny*nz;
            if (ick == 0) {
              __mult_xfwd0_s(tmp.c, gxf, in, i0);
              if(x != 0 ) {
                __mult_xbwd0_s(tmp.c, gxb, in, i0, ixb);
              } else {
#if SU3_RECONSTRUCT_S == 18
                for(int c2=0;c2<3;c2++){
                  for(int c1=0;c1<3;c1++){
                    g.c[c1][c2][0][0] =0;
                    g.c[c1][c2][1][0] =0;
                    for(int j=1;j<VLENS;j++){
                      g.c[c1][c2][0][j] = (*(gxb+i0)).c[c2][c1][0][j-1];
                      g.c[c1][c2][1][j] =-(*(gxb+i0)).c[c2][c1][1][j-1];
                    }
                  }
                }
#elif SU3_RECONSTRUCT_S == 12
                for(int c2=0;c2<2;c2++){
                  for(int c1=0;c1<3;c1++){
                    g.c[c1][c2][0][0] =0;
                    g.c[c1][c2][1][0] =0;
                    for(int j=1;j<VLENS;j++){
                      g.c[c1][c2][0][j] = (*(gxb+i0)).c[c2][c1][0][j-1];
                      g.c[c1][c2][1][j] =-(*(gxb+i0)).c[c2][c1][1][j-1];
                    }
                  }
                }
                for(int j=0;j<VLENS;j++){
                  g.c[0][2][0][j] = g.c[1][0][0][j]*g.c[2][1][0][j]-g.c[1][0][1][j]*g.c[2][1][1][j]-g.c[2][0][0][j]*g.c[1][1][0][j]+g.c[2][0][1][j]*g.c[1][1][1][j];
                  g.c[1][2][0][j] = g.c[2][0][0][j]*g.c[0][1][0][j]-g.c[2][0][1][j]*g.c[0][1][1][j]-g.c[0][0][0][j]*g.c[2][1][0][j]+g.c[0][0][1][j]*g.c[2][1][1][j];
                  g.c[2][2][0][j] = g.c[0][0][0][j]*g.c[1][1][0][j]-g.c[0][0][1][j]*g.c[1][1][1][j]-g.c[1][0][0][j]*g.c[0][1][0][j]+g.c[1][0][1][j]*g.c[0][1][1][j];
                  g.c[0][2][1][j] =-g.c[1][0][0][j]*g.c[2][1][1][j]-g.c[1][0][1][j]*g.c[2][1][0][j]+g.c[2][0][0][j]*g.c[1][1][1][j]+g.c[2][0][1][j]*g.c[1][1][0][j];
                  g.c[1][2][1][j] =-g.c[2][0][0][j]*g.c[0][1][1][j]-g.c[2][0][1][j]*g.c[0][1][0][j]+g.c[0][0][0][j]*g.c[2][1][1][j]+g.c[0][0][1][j]*g.c[2][1][0][j];
                  g.c[2][2][1][j] =-g.c[0][0][0][j]*g.c[1][1][1][j]-g.c[0][0][1][j]*g.c[1][1][0][j]+g.c[1][0][0][j]*g.c[0][1][1][j]+g.c[1][0][1][j]*g.c[0][1][0][j];
                }
#endif
                for(int c=0;c<3;c++){
                  for (int s=0;s<4;s++) {
                    for (int ri=0;ri<2;ri++) {
                      v[s][ri][0]=0; 
                      for(int j=1;j<VLENS;j++){
                        v[s][ri][j]=in[i0].c[c][s][ri][j-1]; 
                      }
                    }
                  }
                  for(int j=0;j<VLENS;j++){
                    a.c[c][0][0][j] = v[0][0][j] - v[3][1][j];
                    a.c[c][1][0][j] = v[1][0][j] - v[2][1][j];
                    a.c[c][0][1][j] = v[0][1][j] + v[3][0][j];
                    a.c[c][1][1][j] = v[1][1][j] + v[2][0][j];
                  }
                }
                __m_fg_s(b.c, a.c, g.c);
                __s_xb_s(tmp.c, b.c);
              }
            } else {
              if(x != nxs-1 ) {
                __mult_xfwd1_s(tmp.c, gxf, in, i0, ixf);
              } else {
                for(int c=0;c<3;c++){
                  for (int s=0;s<4;s++) {
                    for (int ri=0;ri<2;ri++) {
                      for(int j=0;j<VLENS-1;j++){
                        v[s][ri][j]=in[i0].c[c][s][ri][j+1]; 
                      }
                      v[s][ri][VLENS-1]=0; 
                    }
                  }
                  for(int j=0;j<VLENS;j++){
                    a.c[c][0][0][j] = v[0][0][j] + v[3][1][j];
                    a.c[c][1][0][j] = v[1][0][j] + v[2][1][j];
                    a.c[c][0][1][j] = v[0][1][j] - v[3][0][j];
                    a.c[c][1][1][j] = v[1][1][j] - v[2][0][j];
                  }
                }
                __load_su3_s(g.c, gxf, i0);
                __m_fg_s(b.c, a.c, g.c);
                __s_xf_s(tmp.c, b.c);
              }
              __mult_xbwd1_s(tmp.c, gxb, in, i0);
            }
	    if(y != ny-1 ) __mult_yfwd_s(tmp.c, gyf, in, i0, iyf);
	    if(y != 0    ) __mult_ybwd_s(tmp.c, gyb, in,     iyb);
	    if(z != nz-1 ) __mult_zfwd_s(tmp.c, gzf, in, i0, izf);
	    if(z != 0    ) __mult_zbwd_s(tmp.c, gzb, in,     izb);
	    if(t != nt-1 ) __mult_tfwd_s(tmp.c, gtf, in, i0, itf);
	    if(t != 0    ) __mult_tbwd_s(tmp.c, gtb, in,     itb);
	    __mult_clvs( tmp.cv, clvs[i0 + vols*(*pe)].cv);
	    __store_sc_s(out[i0].c, tmp.c);
	  }
	}
      }
    }
  }
  void one_minus_dee_deo_in_s_(int* pe, int* po, scs_t* out, scs_t* in, scs_t* in0, float factor){
    __attribute__((aligned(64)))  scs_t tmp;
    int x, y, z, t;
    __attribute__((aligned(64))) g33s_t g;
    __attribute__((aligned(64))) projscs_t a, b;
    __attribute__((aligned(64))) float v[4][2][VLENS];

    pglus_t gxf = &glus[vols*0 + NDIM*vols*(*pe)];
    pglus_t gxb = &glus[vols*0 + NDIM*vols*(*po)];
    pglus_t gyf = &glus[vols*1 + NDIM*vols*(*pe)];
    pglus_t gyb = &glus[vols*1 + NDIM*vols*(*po)];
    pglus_t gzf = &glus[vols*2 + NDIM*vols*(*pe)];
    pglus_t gzb = &glus[vols*2 + NDIM*vols*(*po)];
    pglus_t gtf = &glus[vols*3 + NDIM*vols*(*pe)];
    pglus_t gtb = &glus[vols*3 + NDIM*vols*(*po)];
    int jjj= ((*pe) + ny*py + nz*pz + nt*pt)%2;
#pragma omp for private(t,z,y,x,tmp, g,a,b,v) collapse(3) schedule(static)
    for(t=0; t<nt; t++){
      for(z=0; z<nz; z++){
	for(y=0; y<ny; y++){
	  int ick= ( jjj + y     + z + t)%2;
	  for(x=0; x<nxs; x++){
	    __zero_sc_s(tmp.c);
	    int i0  =      x        + nxs*y + nxs*ny*z + nxs*ny*nz*t;
	    int ixf = (    x+1)%nxs + nxs*y + nxs*ny*z + nxs*ny*nz*t;
	    int ixb = (nxs+x-1)%nxs + nxs*y + nxs*ny*z + nxs*ny*nz*t;
	    int iyf = x + (y+1)*nxs + nxs*ny*z + nxs*ny*nz*t;
	    int iyb = x + (y-1)*nxs + nxs*ny*z + nxs*ny*nz*t;
	    int izf = x + nxs*y + (z+1)*nxs*ny + nxs*ny*nz*t;
	    int izb = x + nxs*y + (z-1)*nxs*ny + nxs*ny*nz*t;
	    int itf = x + nxs*y + nxs*ny*z + (t+1)*nxs*ny*nz;
	    int itb = x + nxs*y + nxs*ny*z + (t-1)*nxs*ny*nz;
            if (ick == 0) {
              __mult_xfwd0_s(tmp.c, gxf, in, i0);
              if(x != 0 ) {
                __mult_xbwd0_s(tmp.c, gxb, in, i0, ixb);
              } else {
#if SU3_RECONSTRUCT_S == 18
                for(int c2=0;c2<3;c2++){
                  for(int c1=0;c1<3;c1++){
                    g.c[c1][c2][0][0] =0;
                    g.c[c1][c2][1][0] =0;
                    for(int j=1;j<VLENS;j++){
                      g.c[c1][c2][0][j] = (*(gxb+i0)).c[c2][c1][0][j-1];
                      g.c[c1][c2][1][j] =-(*(gxb+i0)).c[c2][c1][1][j-1];
                    }
                  }
                }
#elif SU3_RECONSTRUCT_S == 12
                for(int c2=0;c2<2;c2++){
                  for(int c1=0;c1<3;c1++){
                    g.c[c1][c2][0][0] =0;
                    g.c[c1][c2][1][0] =0;
                    for(int j=1;j<VLENS;j++){
                      g.c[c1][c2][0][j] = (*(gxb+i0)).c[c2][c1][0][j-1];
                      g.c[c1][c2][1][j] =-(*(gxb+i0)).c[c2][c1][1][j-1];
                    }
                  }
                }
                for(int j=0;j<VLENS;j++){
                  g.c[0][2][0][j] = g.c[1][0][0][j]*g.c[2][1][0][j]-g.c[1][0][1][j]*g.c[2][1][1][j]-g.c[2][0][0][j]*g.c[1][1][0][j]+g.c[2][0][1][j]*g.c[1][1][1][j];
                  g.c[1][2][0][j] = g.c[2][0][0][j]*g.c[0][1][0][j]-g.c[2][0][1][j]*g.c[0][1][1][j]-g.c[0][0][0][j]*g.c[2][1][0][j]+g.c[0][0][1][j]*g.c[2][1][1][j];
                  g.c[2][2][0][j] = g.c[0][0][0][j]*g.c[1][1][0][j]-g.c[0][0][1][j]*g.c[1][1][1][j]-g.c[1][0][0][j]*g.c[0][1][0][j]+g.c[1][0][1][j]*g.c[0][1][1][j];
                  g.c[0][2][1][j] =-g.c[1][0][0][j]*g.c[2][1][1][j]-g.c[1][0][1][j]*g.c[2][1][0][j]+g.c[2][0][0][j]*g.c[1][1][1][j]+g.c[2][0][1][j]*g.c[1][1][0][j];
                  g.c[1][2][1][j] =-g.c[2][0][0][j]*g.c[0][1][1][j]-g.c[2][0][1][j]*g.c[0][1][0][j]+g.c[0][0][0][j]*g.c[2][1][1][j]+g.c[0][0][1][j]*g.c[2][1][0][j];
                  g.c[2][2][1][j] =-g.c[0][0][0][j]*g.c[1][1][1][j]-g.c[0][0][1][j]*g.c[1][1][0][j]+g.c[1][0][0][j]*g.c[0][1][1][j]+g.c[1][0][1][j]*g.c[0][1][0][j];
                }
#endif
                for(int c=0;c<3;c++){
                  for (int s=0;s<4;s++) {
                    for (int ri=0;ri<2;ri++) {
                      v[s][ri][0]=0; 
                      for(int j=1;j<VLENS;j++){
                        v[s][ri][j]=in[i0].c[c][s][ri][j-1]; 
                      }
                    }
                  }
                  for(int j=0;j<VLENS;j++){
                    a.c[c][0][0][j] = v[0][0][j] - v[3][1][j];
                    a.c[c][1][0][j] = v[1][0][j] - v[2][1][j];
                    a.c[c][0][1][j] = v[0][1][j] + v[3][0][j];
                    a.c[c][1][1][j] = v[1][1][j] + v[2][0][j];
                  }
                }
                __m_fg_s(b.c, a.c, g.c);
                __s_xb_s(tmp.c, b.c);
              }
            } else {
              if(x != nxs-1 ) {
                __mult_xfwd1_s(tmp.c, gxf, in, i0, ixf);
              } else {
                for(int c=0;c<3;c++){
                  for (int s=0;s<4;s++) {
                    for (int ri=0;ri<2;ri++) {
                      for(int j=0;j<VLENS-1;j++){
                        v[s][ri][j]=in[i0].c[c][s][ri][j+1]; 
                      }
                      v[s][ri][VLENS-1]=0; 
                    }
                  }
                  for(int j=0;j<VLENS;j++){
                    a.c[c][0][0][j] = v[0][0][j] + v[3][1][j];
                    a.c[c][1][0][j] = v[1][0][j] + v[2][1][j];
                    a.c[c][0][1][j] = v[0][1][j] - v[3][0][j];
                    a.c[c][1][1][j] = v[1][1][j] - v[2][0][j];
                  }
                }
                __load_su3_s(g.c, gxf, i0);
                __m_fg_s(b.c, a.c, g.c);
                __s_xf_s(tmp.c, b.c);
              }
              __mult_xbwd1_s(tmp.c, gxb, in, i0);
            }
	    if(y != ny-1 ) __mult_yfwd_s(tmp.c, gyf, in, i0, iyf);
	    if(y != 0    ) __mult_ybwd_s(tmp.c, gyb, in,     iyb);
	    if(z != nz-1 ) __mult_zfwd_s(tmp.c, gzf, in, i0, izf);
	    if(z != 0    ) __mult_zbwd_s(tmp.c, gzb, in,     izb);
	    if(t != nt-1 ) __mult_tfwd_s(tmp.c, gtf, in, i0, itf);
	    if(t != 0    ) __mult_tbwd_s(tmp.c, gtb, in,     itb);
	    __mult_clvs( tmp.cv, clvs[i0 + vols*(*pe)].cv);
	    __load_sc_s(out[i0].c, in0[i0].c);
	    __store_add_sc_s(out[i0].c, tmp.c, factor);
	  }
	}
      }
    }
  }
#ifdef __cplusplus
}
#endif
