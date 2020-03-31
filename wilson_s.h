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
#ifndef WILSON_S_H
#define WILSON_S_H
#include "qws.h"

static inline void __zero_sc_s(float out[3][4][2][VLENS]) {
  int c, s, ri, j;
  for(c=0;c<3;c++){
    for(s=0;s<4;s++){
      for(ri=0;ri<2;ri++){
#pragma loop nounroll
	for(j=0;j<VLENS;j++){
	  out[c][s][ri][j] = 0;
	}
      }
    }
  }
}
static inline void __load_sc_s(float out[3][4][2][VLENS], float in[3][4][2][VLENS]) {
  int c, s, ri, j;
  for(c=0;c<3;c++){
    for(s=0;s<4;s++){
      for(ri=0;ri<2;ri++){
#pragma loop nounroll
	for(j=0;j<VLENS;j++){
	  out[c][s][ri][j] = in[c][s][ri][j];
	}
      }
    }
  }
}

static inline void __store_sc_s(float out[3][4][2][VLENS], float sc[3][4][2][VLENS]) {
  int c, s, j;
  for(c=0;c<3;c++){
    for(s=0;s<4;s++){
      //#pragma vector nontemporal(out)
      //#pragma vector always
#pragma loop nounroll
      for(j=0;j<VLENS;j++){
	out[c][s][0][j]=sc[c][s][0][j];
	out[c][s][1][j]=sc[c][s][1][j];
      }
    }
  }
}

static inline void __store_add_sc_s(float out[3][4][2][VLENS], float sc[3][4][2][VLENS], float a) {
  int c, s, j;
  for(c=0;c<3;c++){
    for(s=0;s<4;s++){
#pragma loop nounroll
#pragma loop norecurrence
      for(j=0;j<VLENS;j++){
	out[c][s][0][j]+=sc[c][s][0][j] * a;
	out[c][s][1][j]+=sc[c][s][1][j] * a;
      }
    }
  }
}

static inline void __m_fg_s(float b[3][2][2][VLENS], float a[3][2][2][VLENS], float u[3][3][2][VLENS]) {
  int c1, c2, ud, ri, j;
  for(c1=0;c1<3;c1++){
    for(ud=0;ud<2;ud++){
      for(ri=0;ri<2;ri++){
#pragma loop nounroll
	for(j=0;j<VLENS;j++){
	  b[c1][ud][ri][j] = 0;
	}
      }
    }
  }
  //#pragma unroll(3)
  for(c2=0;c2<3;c2++){
    for(c1=0;c1<3;c1++){
      for(ud=0;ud<2;ud++){
#pragma loop nounroll
	for(j=0;j<VLENS;j++){
	  b[c1][ud][0][j] += a[c2][ud][0][j] * u[c2][c1][0][j] - a[c2][ud][1][j] * u[c2][c1][1][j];
	  b[c1][ud][1][j] += a[c2][ud][0][j] * u[c2][c1][1][j] + a[c2][ud][1][j] * u[c2][c1][0][j];
	}
      }
    }
  }

}

#if SU3_RECONSTRUCT_S == 18
static inline void __load_su3_s(float g[3][3][2][VLENS], pglus_t gp, int i) {
  int c1, c2, j;
  for(c2=0;c2<3;c2++){
    for(c1=0;c1<3;c1++){
      //#pragma simd
#pragma loop nounroll
      for(j=0;j<VLENS;j++){
        g[c2][c1][0][j] = (*(gp+i)).c[c2][c1][0][j];
        g[c2][c1][1][j] = (*(gp+i)).c[c2][c1][1][j];
      }
    }
  }
}
#elif SU3_RECONSTRUCT_S == 12
static inline void __load_su3_s(float g[3][3][2][VLENS], pglus_t gp, int i) {
  int c1, c2, j;
  for(c2=0;c2<2;c2++){
    for(c1=0;c1<3;c1++){
      #pragma simd
#pragma loop nounroll
      for(j=0;j<VLENS;j++){
	g[c2][c1][0][j] = (*(gp+i)).c[c2][c1][0][j];
        g[c2][c1][1][j] = (*(gp+i)).c[c2][c1][1][j];
      }
    }
  }
  for(j=0;j<VLENS;j++){
    g[2][0][0][j] = g[0][1][0][j]*g[1][2][0][j]-g[0][1][1][j]*g[1][2][1][j]-g[0][2][0][j]*g[1][1][0][j]+g[0][2][1][j]*g[1][1][1][j];
    g[2][1][0][j] = g[0][2][0][j]*g[1][0][0][j]-g[0][2][1][j]*g[1][0][1][j]-g[0][0][0][j]*g[1][2][0][j]+g[0][0][1][j]*g[1][2][1][j];
    g[2][2][0][j] = g[0][0][0][j]*g[1][1][0][j]-g[0][0][1][j]*g[1][1][1][j]-g[0][1][0][j]*g[1][0][0][j]+g[0][1][1][j]*g[1][0][1][j];
    g[2][0][1][j] =-g[0][1][0][j]*g[1][2][1][j]-g[0][1][1][j]*g[1][2][0][j]+g[0][2][0][j]*g[1][1][1][j]+g[0][2][1][j]*g[1][1][0][j];
    g[2][1][1][j] =-g[0][2][0][j]*g[1][0][1][j]-g[0][2][1][j]*g[1][0][0][j]+g[0][0][0][j]*g[1][2][1][j]+g[0][0][1][j]*g[1][2][0][j];
    g[2][2][1][j] =-g[0][0][0][j]*g[1][1][1][j]-g[0][0][1][j]*g[1][1][0][j]+g[0][1][0][j]*g[1][0][1][j]+g[0][1][1][j]*g[1][0][0][j];
  }
}
#endif

#if SU3_RECONSTRUCT_S == 18
static inline void __load_su3d0_s(float g[3][3][2][VLENS], pglus_t gp, int i) {
  int c1, c2, j;
  for(c2=0;c2<3;c2++){
    for(c1=0;c1<3;c1++){
#pragma loop nounroll
      for(j=0;j<VLENS;j++){
        g[c1][c2][0][j] = (*(gp+i)).c[c2][c1][0][j];
        g[c1][c2][1][j] =-(*(gp+i)).c[c2][c1][1][j];
      }
    }
  }
}
#elif SU3_RECONSTRUCT_S == 12
static inline void __load_su3d0_s(float g[3][3][2][VLENS], pglus_t gp, int i) {
  int c1, c2, j;
  for(c2=0;c2<2;c2++){
    for(c1=0;c1<3;c1++){
#pragma loop nounroll
      for(j=0;j<VLENS;j++){
        g[c1][c2][0][j] = (*(gp+i)).c[c2][c1][0][j];
        g[c1][c2][1][j] =-(*(gp+i)).c[c2][c1][1][j];
      }
    }
  }
  for(j=0;j<VLENS;j++){
    g[0][2][0][j] = g[1][0][0][j]*g[2][1][0][j]-g[1][0][1][j]*g[2][1][1][j]-g[2][0][0][j]*g[1][1][0][j]+g[2][0][1][j]*g[1][1][1][j];
    g[1][2][0][j] = g[2][0][0][j]*g[0][1][0][j]-g[2][0][1][j]*g[0][1][1][j]-g[0][0][0][j]*g[2][1][0][j]+g[0][0][1][j]*g[2][1][1][j];
    g[2][2][0][j] = g[0][0][0][j]*g[1][1][0][j]-g[0][0][1][j]*g[1][1][1][j]-g[1][0][0][j]*g[0][1][0][j]+g[1][0][1][j]*g[0][1][1][j];
    g[0][2][1][j] =-g[1][0][0][j]*g[2][1][1][j]-g[1][0][1][j]*g[2][1][0][j]+g[2][0][0][j]*g[1][1][1][j]+g[2][0][1][j]*g[1][1][0][j];
    g[1][2][1][j] =-g[2][0][0][j]*g[0][1][1][j]-g[2][0][1][j]*g[0][1][0][j]+g[0][0][0][j]*g[2][1][1][j]+g[0][0][1][j]*g[2][1][0][j];
    g[2][2][1][j] =-g[0][0][0][j]*g[1][1][1][j]-g[0][0][1][j]*g[1][1][0][j]+g[1][0][0][j]*g[0][1][1][j]+g[1][0][1][j]*g[0][1][0][j];
  }
}
#endif

#if SU3_RECONSTRUCT_S == 18
static inline void __load_su3d1_s(float g[3][3][2][VLENS], pglus_t gp, int i0, int ib) {
  int c1, c2, j;
  for(c2=0;c2<3;c2++){
    for(c1=0;c1<3;c1++){
      //#pragma loop nounroll
      g[c1][c2][0][0] = (*(gp+ib)).c[c2][c1][0][VLENS-1];
      g[c1][c2][1][0] =-(*(gp+ib)).c[c2][c1][1][VLENS-1];
      for(j=1;j<VLENS;j++){
        g[c1][c2][0][j] = (*(gp+i0)).c[c2][c1][0][j-1];
        g[c1][c2][1][j] =-(*(gp+i0)).c[c2][c1][1][j-1];
      }
    }
  }
}
#elif SU3_RECONSTRUCT_S == 12
static inline void __load_su3d1_s(float g[3][3][2][VLENS], pglus_t gp, int i0, int ib) {
  int c1, c2, j;
  for(c2=0;c2<2;c2++){
    for(c1=0;c1<3;c1++){
#pragma loop nounroll
      g[c1][c2][0][0] = (*(gp+ib)).c[c2][c1][0][VLENS-1];
      g[c1][c2][1][0] =-(*(gp+ib)).c[c2][c1][1][VLENS-1];
      for(j=1;j<VLENS;j++){
        g[c1][c2][0][j] = (*(gp+i0)).c[c2][c1][0][j-1];
        g[c1][c2][1][j] =-(*(gp+i0)).c[c2][c1][1][j-1];
      }
    }
  }
  for(j=0;j<VLENS;j++){
    g[0][2][0][j] = g[1][0][0][j]*g[2][1][0][j]-g[1][0][1][j]*g[2][1][1][j]-g[2][0][0][j]*g[1][1][0][j]+g[2][0][1][j]*g[1][1][1][j];
    g[1][2][0][j] = g[2][0][0][j]*g[0][1][0][j]-g[2][0][1][j]*g[0][1][1][j]-g[0][0][0][j]*g[2][1][0][j]+g[0][0][1][j]*g[2][1][1][j];
    g[2][2][0][j] = g[0][0][0][j]*g[1][1][0][j]-g[0][0][1][j]*g[1][1][1][j]-g[1][0][0][j]*g[0][1][0][j]+g[1][0][1][j]*g[0][1][1][j];
    g[0][2][1][j] =-g[1][0][0][j]*g[2][1][1][j]-g[1][0][1][j]*g[2][1][0][j]+g[2][0][0][j]*g[1][1][1][j]+g[2][0][1][j]*g[1][1][0][j];
    g[1][2][1][j] =-g[2][0][0][j]*g[0][1][1][j]-g[2][0][1][j]*g[0][1][0][j]+g[0][0][0][j]*g[2][1][1][j]+g[0][0][1][j]*g[2][1][0][j];
    g[2][2][1][j] =-g[0][0][0][j]*g[1][1][1][j]-g[0][0][1][j]*g[1][1][0][j]+g[1][0][0][j]*g[0][1][1][j]+g[1][0][1][j]*g[0][1][0][j];
  }
}
#endif

#define __load_sc0_s(ri, s, c, i0){\
  for(int __load_sc0_i = 0; __load_sc0_i < VLENS; __load_sc0_i++) {\
    v[s][ri][__load_sc0_i]=in[i0].c[c][s][ri][__load_sc0_i];\
  }\
  }

#if VLENS == 16
#define __load_scf_s(ri, s, c, i0, i1){\
  v[s][ri][ 0]=in[i0].c[c][s][ri][ 1];\
  v[s][ri][ 1]=in[i0].c[c][s][ri][ 2];\
  v[s][ri][ 2]=in[i0].c[c][s][ri][ 3];\
  v[s][ri][ 3]=in[i0].c[c][s][ri][ 4];\
  v[s][ri][ 4]=in[i0].c[c][s][ri][ 5];\
  v[s][ri][ 5]=in[i0].c[c][s][ri][ 6];\
  v[s][ri][ 6]=in[i0].c[c][s][ri][ 7];\
  v[s][ri][ 7]=in[i0].c[c][s][ri][ 8];\
  v[s][ri][ 8]=in[i0].c[c][s][ri][ 9];\
  v[s][ri][ 9]=in[i0].c[c][s][ri][10];\
  v[s][ri][10]=in[i0].c[c][s][ri][11];\
  v[s][ri][11]=in[i0].c[c][s][ri][12];\
  v[s][ri][12]=in[i0].c[c][s][ri][13];\
  v[s][ri][13]=in[i0].c[c][s][ri][14];\
  v[s][ri][14]=in[i0].c[c][s][ri][15];\
  v[s][ri][15]=in[i1].c[c][s][ri][ 0];}
#define __load_scb_s(ri, s, c, i0, i1){\
  v[s][ri][ 0]=in[i1].c[c][s][ri][15];\
  v[s][ri][ 1]=in[i0].c[c][s][ri][ 0];\
  v[s][ri][ 2]=in[i0].c[c][s][ri][ 1];\
  v[s][ri][ 3]=in[i0].c[c][s][ri][ 2];\
  v[s][ri][ 4]=in[i0].c[c][s][ri][ 3];\
  v[s][ri][ 5]=in[i0].c[c][s][ri][ 4];\
  v[s][ri][ 6]=in[i0].c[c][s][ri][ 5];\
  v[s][ri][ 7]=in[i0].c[c][s][ri][ 6];\
  v[s][ri][ 8]=in[i0].c[c][s][ri][ 7];\
  v[s][ri][ 9]=in[i0].c[c][s][ri][ 8];\
  v[s][ri][10]=in[i0].c[c][s][ri][ 9];\
  v[s][ri][11]=in[i0].c[c][s][ri][10];\
  v[s][ri][12]=in[i0].c[c][s][ri][11];\
  v[s][ri][13]=in[i0].c[c][s][ri][12];\
  v[s][ri][14]=in[i0].c[c][s][ri][13];\
  v[s][ri][15]=in[i0].c[c][s][ri][14];}
#elif VLENS == 8
#define __load_scf_s(ri, s, c, i0, i1){\
  v[s][ri][0]=in[i0].c[c][s][ri][1];\
  v[s][ri][1]=in[i0].c[c][s][ri][2];\
  v[s][ri][2]=in[i0].c[c][s][ri][3];\
  v[s][ri][3]=in[i0].c[c][s][ri][4];\
  v[s][ri][4]=in[i0].c[c][s][ri][5];\
  v[s][ri][5]=in[i0].c[c][s][ri][6];\
  v[s][ri][6]=in[i0].c[c][s][ri][7];\
  v[s][ri][7]=in[i1].c[c][s][ri][0];}
#define __load_scb_s(ri, s, c, i0, i1){\
  v[s][ri][0]=in[i1].c[c][s][ri][7];\
  v[s][ri][1]=in[i0].c[c][s][ri][0];\
  v[s][ri][2]=in[i0].c[c][s][ri][1];\
  v[s][ri][3]=in[i0].c[c][s][ri][2];\
  v[s][ri][4]=in[i0].c[c][s][ri][3];\
  v[s][ri][5]=in[i0].c[c][s][ri][4];\
  v[s][ri][6]=in[i0].c[c][s][ri][5];\
  v[s][ri][7]=in[i0].c[c][s][ri][6];}
#elif VLENS == 4
#define __load_scf_s(ri, s, c, i0, i1){\
  v[s][ri][0]=in[i0].c[c][s][ri][1];\
  v[s][ri][1]=in[i0].c[c][s][ri][2];\
  v[s][ri][2]=in[i0].c[c][s][ri][3];\
  v[s][ri][3]=in[i1].c[c][s][ri][0];}
#define __load_scb_s(ri, s, c, i0, i1){\
  v[s][ri][0]=in[i1].c[c][s][ri][3];\
  v[s][ri][1]=in[i0].c[c][s][ri][0];\
  v[s][ri][2]=in[i0].c[c][s][ri][1];\
  v[s][ri][3]=in[i0].c[c][s][ri][2];}
#elif VLENS == 2
#define __load_scf_s(ri, s, c, i0, i1){\
  v[s][ri][0]=in[i0].c[c][s][ri][1];\
  v[s][ri][1]=in[i1].c[c][s][ri][0];}
#define __load_scb_s(ri, s, c, i0, i1){\
  v[s][ri][0]=in[i1].c[c][s][ri][1];\
  v[s][ri][1]=in[i0].c[c][s][ri][0];}
#endif

#define __load_sc0_sri_s(ri, s, c, i0){\
  __load_sc0_s(0, 0, c, i0);	     \
  __load_sc0_s(1, 0, c, i0);	     \
  __load_sc0_s(0, 1, c, i0);	     \
  __load_sc0_s(1, 1, c, i0);	     \
  __load_sc0_s(0, 2, c, i0);	     \
  __load_sc0_s(1, 2, c, i0);	     \
  __load_sc0_s(0, 3, c, i0);	     \
  __load_sc0_s(1, 3, c, i0);}

#define __load_scf_sri_s(ri, s, c, i0, ifwd){\
  __load_scf_s(0, 0, c, i0, ifwd);	   \
  __load_scf_s(1, 0, c, i0, ifwd);	   \
  __load_scf_s(0, 1, c, i0, ifwd);	   \
  __load_scf_s(1, 1, c, i0, ifwd);	   \
  __load_scf_s(0, 2, c, i0, ifwd);	   \
  __load_scf_s(1, 2, c, i0, ifwd);	   \
  __load_scf_s(0, 3, c, i0, ifwd);	   \
  __load_scf_s(1, 3, c, i0, ifwd);}

#define __load_scb_sri_s(ri, s, c, i0, ibwd){\
  __load_scb_s(0, 0, c, i0, ibwd);	   \
  __load_scb_s(1, 0, c, i0, ibwd);	   \
  __load_scb_s(0, 1, c, i0, ibwd);	   \
  __load_scb_s(1, 1, c, i0, ibwd);	   \
  __load_scb_s(0, 2, c, i0, ibwd);	   \
  __load_scb_s(1, 2, c, i0, ibwd);	   \
  __load_scb_s(0, 3, c, i0, ibwd);	   \
  __load_scb_s(1, 3, c, i0, ibwd);}


static inline void __s_xf_s(float out[3][4][2][VLENS], float b[3][2][2][VLENS]) {
  int c, j;
  for(c=0;c<3;c++){
#pragma loop nounroll
    for(j=0;j<VLENS;j++){
      out[c][0][0][j] +=  + b[c][0][0][j];
      out[c][1][0][j] +=  + b[c][1][0][j];
      out[c][2][0][j] +=  - b[c][1][1][j];
      out[c][3][0][j] +=  - b[c][0][1][j];
      out[c][0][1][j] +=  + b[c][0][1][j];
      out[c][1][1][j] +=  + b[c][1][1][j];
      out[c][2][1][j] +=  + b[c][1][0][j];
      out[c][3][1][j] +=  + b[c][0][0][j];
    }
  }
}

static inline void __s_xb1_s(float out[3][4][2][VLENS], float b[3][2][2]) {
  int c;
  for(c=0;c<3;c++){
    out[c][0][0][0] +=  + b[c][0][0];
    out[c][1][0][0] +=  + b[c][1][0];
    out[c][2][0][0] +=  + b[c][1][1];
    out[c][3][0][0] +=  + b[c][0][1];
    out[c][0][1][0] +=  + b[c][0][1];
    out[c][1][1][0] +=  + b[c][1][1];
    out[c][2][1][0] +=  - b[c][1][0];
    out[c][3][1][0] +=  - b[c][0][0];
  }
}

static inline void __s_xb_s(float out[3][4][2][VLENS], float b[3][2][2][VLENS]) {
  int c, j;
  for(c=0;c<3;c++){
#pragma loop nounroll
    for(j=0;j<VLENS;j++){
      out[c][0][0][j] +=  + b[c][0][0][j];
      out[c][1][0][j] +=  + b[c][1][0][j];
      out[c][2][0][j] +=  + b[c][1][1][j];
      out[c][3][0][j] +=  + b[c][0][1][j];
      out[c][0][1][j] +=  + b[c][0][1][j];
      out[c][1][1][j] +=  + b[c][1][1][j];
      out[c][2][1][j] +=  - b[c][1][0][j];
      out[c][3][1][j] +=  - b[c][0][0][j];
    }
  }
}


static inline void __s_yf_s(float out[3][4][2][VLENS], float b[3][2][2][VLENS]) {
  int c, j;
  for(c=0;c<3;c++){
#pragma loop nounroll
    for(j=0;j<VLENS;j++){
      out[c][0][0][j] +=  + b[c][0][0][j];
      out[c][1][0][j] +=  + b[c][1][0][j];
      out[c][2][0][j] +=  + b[c][1][0][j];
      out[c][3][0][j] +=  - b[c][0][0][j];
      out[c][0][1][j] +=  + b[c][0][1][j];
      out[c][1][1][j] +=  + b[c][1][1][j];
      out[c][2][1][j] +=  + b[c][1][1][j];
      out[c][3][1][j] +=  - b[c][0][1][j];
    }
  }
}
static inline void __s_yb_s(float out[3][4][2][VLENS], float b[3][2][2][VLENS]) {
  int c, j;
  for(c=0;c<3;c++){
#pragma loop nounroll
    for(j=0;j<VLENS;j++){
      out[c][0][0][j] +=  + b[c][0][0][j];
      out[c][1][0][j] +=  + b[c][1][0][j];
      out[c][2][0][j] +=  - b[c][1][0][j];
      out[c][3][0][j] +=  + b[c][0][0][j];
      out[c][0][1][j] +=  + b[c][0][1][j];
      out[c][1][1][j] +=  + b[c][1][1][j];
      out[c][2][1][j] +=  - b[c][1][1][j];
      out[c][3][1][j] +=  + b[c][0][1][j];
    }
  }
}
static inline void __s_zf_s(float out[3][4][2][VLENS], float b[3][2][2][VLENS]) {
  int c, j;
  for(c=0;c<3;c++){
#pragma loop nounroll
    for(j=0;j<VLENS;j++){
      out[c][0][0][j] +=  + b[c][0][0][j];
      out[c][1][0][j] +=  + b[c][1][0][j];
      out[c][2][0][j] +=  - b[c][0][1][j];
      out[c][3][0][j] +=  + b[c][1][1][j];
      out[c][0][1][j] +=  + b[c][0][1][j];
      out[c][1][1][j] +=  + b[c][1][1][j];
      out[c][2][1][j] +=  + b[c][0][0][j];
      out[c][3][1][j] +=  - b[c][1][0][j];
    }
  }
}
static inline void __s_zb_s(float out[3][4][2][VLENS], float b[3][2][2][VLENS]) {
  int c, j;
  for(c=0;c<3;c++){
#pragma loop nounroll
    for(j=0;j<VLENS;j++){
      out[c][0][0][j] +=  + b[c][0][0][j];
      out[c][1][0][j] +=  + b[c][1][0][j];
      out[c][2][0][j] +=  + b[c][0][1][j];
      out[c][3][0][j] +=  - b[c][1][1][j];
      out[c][0][1][j] +=  + b[c][0][1][j];
      out[c][1][1][j] +=  + b[c][1][1][j];
      out[c][2][1][j] +=  - b[c][0][0][j];
      out[c][3][1][j] +=  + b[c][1][0][j];
    }
  }
}
static inline void __s_tf_s(float out[3][4][2][VLENS], float b[3][2][2][VLENS]) {
  int c, j;
  for(c=0;c<3;c++){
#pragma loop nounroll
    for(j=0;j<VLENS;j++){
      out[c][2][0][j] +=  2 * b[c][0][0][j];
      out[c][3][0][j] +=  2 * b[c][1][0][j];
      out[c][2][1][j] +=  2 * b[c][0][1][j];
      out[c][3][1][j] +=  2 * b[c][1][1][j];
    }
  }
}
static inline void __s_tb_s(float out[3][4][2][VLENS], float b[3][2][2][VLENS]) {
  int c, j;
  for(c=0;c<3;c++){
#pragma loop nounroll
    for(j=0;j<VLENS;j++){
      out[c][0][0][j] +=  2 * b[c][0][0][j];
      out[c][1][0][j] +=  2 * b[c][1][0][j];
      out[c][0][1][j] +=  2 * b[c][0][1][j];
      out[c][1][1][j] +=  2 * b[c][1][1][j];
    }
  }
}
// bc is dummy so far, this will be changed in future
static inline void __s_tf_bc_s(float out[3][4][2][VLENS], float b[3][2][2][VLENS], float bc) {
  int c, j;
  for(c=0;c<3;c++){
#pragma loop nounroll
    for(j=0;j<VLENS;j++){
      out[c][2][0][j] +=  bc * b[c][0][0][j];
      out[c][3][0][j] +=  bc * b[c][1][0][j];
      out[c][2][1][j] +=  bc * b[c][0][1][j];
      out[c][3][1][j] +=  bc * b[c][1][1][j];
    }
  }
}
static inline void __s_tb_bc_s(float out[3][4][2][VLENS], float b[3][2][2][VLENS], float bc) {
  int c, j;
  for(c=0;c<3;c++){
#pragma loop nounroll
    for(j=0;j<VLENS;j++){
      out[c][0][0][j] +=  bc * b[c][0][0][j];
      out[c][1][0][j] +=  bc * b[c][1][0][j];
      out[c][0][1][j] +=  bc * b[c][0][1][j];
      out[c][1][1][j] +=  bc * b[c][1][1][j];
    }
  }
}

static inline void __l_xf_s(float a[3][2][2], scs_t* in, int ifwd) {
  float v[4][2];
  for(int c=0;c<3;c++){
    for (int s=0;s<4;s++) {
      for (int ri=0;ri<2;ri++) {
        v[s][ri]=in[ifwd].c[c][s][ri][0];
      }
    }
    a[c][0][0] = v[0][0] + v[3][1];
    a[c][1][0] = v[1][0] + v[2][1];
    a[c][0][1] = v[0][1] - v[3][0];
    a[c][1][1] = v[1][1] - v[2][0];
  }
}

static inline void __l_yf_s(float a[3][2][2][VLENS], scs_t* in, int ifwd) {
  int c, j;
#pragma loop nounroll
  for(c=0;c<3;c++){
    float v[4][2][VLENS];
    __load_sc0_sri_s(ri, s, c, ifwd);
#pragma loop norecurrence
    for(j=0;j<VLENS;j++){
      a[c][0][0][j] = v[0][0][j] - v[3][0][j];
      a[c][1][0][j] = v[1][0][j] + v[2][0][j];
      a[c][0][1][j] = v[0][1][j] - v[3][1][j];
      a[c][1][1][j] = v[1][1][j] + v[2][1][j];
    }
  }
}
static inline void __l_yb_s(float a[3][2][2][VLENS], scs_t* in, int ibwd) {
  int c, j;
#pragma loop nounroll
  for(c=0;c<3;c++){
    float v[4][2][VLENS];
    __load_sc0_sri_s(ri, s, c, ibwd);
#pragma loop norecurrence
    for(j=0;j<VLENS;j++){
      a[c][0][0][j] = v[0][0][j] + v[3][0][j];
      a[c][1][0][j] = v[1][0][j] - v[2][0][j];
      a[c][0][1][j] = v[0][1][j] + v[3][1][j];
      a[c][1][1][j] = v[1][1][j] - v[2][1][j];
    }
  }
}
static inline void __l_zf_s(float a[3][2][2][VLENS], scs_t* in, int ifwd) {
  int c, j;
#pragma loop nounroll
  for(c=0;c<3;c++){
    float v[4][2][VLENS];
    __load_sc0_sri_s(ri, s, c, ifwd);
#pragma loop norecurrence
    for(j=0;j<VLENS;j++){
      a[c][0][0][j] = v[0][0][j] + v[2][1][j];
      a[c][1][0][j] = v[1][0][j] - v[3][1][j];
      a[c][0][1][j] = v[0][1][j] - v[2][0][j];
      a[c][1][1][j] = v[1][1][j] + v[3][0][j];
    }
  }
}
static inline void __l_zb_s(float a[3][2][2][VLENS], scs_t* in, int ibwd) {
  int c, j;
#pragma loop nounroll
  for(c=0;c<3;c++){
    float v[4][2][VLENS];
    __load_sc0_sri_s(ri, s, c, ibwd);
#pragma loop norecurrence
    for(j=0;j<VLENS;j++){
      a[c][0][0][j] = v[0][0][j] - v[2][1][j];
      a[c][1][0][j] = v[1][0][j] + v[3][1][j];
      a[c][0][1][j] = v[0][1][j] + v[2][0][j];
      a[c][1][1][j] = v[1][1][j] - v[3][0][j];
    }
  }
}
static inline void __l_tf_s(float a[3][2][2][VLENS], scs_t* in, int ifwd) {
  int c, j;
#pragma loop nounroll
  for(c=0;c<3;c++){
    float v[4][2][VLENS];
    __load_sc0_sri_s(ri, s, c, ifwd);
#pragma loop norecurrence
    for(j=0;j<VLENS;j++){
      a[c][0][0][j] = v[2][0][j];
      a[c][1][0][j] = v[3][0][j];
      a[c][0][1][j] = v[2][1][j];
      a[c][1][1][j] = v[3][1][j];
    }
  }
}
static inline void __l_tb_s(float a[3][2][2][VLENS], scs_t* in, int ibwd) {
  int c, j;
#pragma loop nounroll
  for(c=0;c<3;c++){
    float v[4][2][VLENS];
    __load_sc0_sri_s(ri, s, c, ibwd);
#pragma loop norecurrence
    for(j=0;j<VLENS;j++){
      a[c][0][0][j] = v[0][0][j];
      a[c][1][0][j] = v[1][0][j];
      a[c][0][1][j] = v[0][1][j];
      a[c][1][1][j] = v[1][1][j];
    }
  }
}

//----------------------------------------------------------------------------------------X forward
static inline void __mult_xfwd0_s(float out[3][4][2][VLENS], pglus_t gp, scs_t* in, int i0) {
  int c, j;
  __attribute__((aligned(64))) g33s_t g;
  __attribute__((aligned(64))) projscs_t a, b;

  __load_su3_s(g.c, gp, i0);
#pragma loop nounroll
  for(c=0;c<3;c++){
    __attribute__((aligned(64))) float v[4][2][VLENS];
    __load_sc0_sri_s(ri, s, c, i0);
    for(j=0;j<VLENS;j++){
      a.c[c][0][0][j] = v[0][0][j] + v[3][1][j];
      a.c[c][1][0][j] = v[1][0][j] + v[2][1][j];
      a.c[c][0][1][j] = v[0][1][j] - v[3][0][j];
      a.c[c][1][1][j] = v[1][1][j] - v[2][0][j];
    }
  }
  __m_fg_s(b.c, a.c, g.c);
  __s_xf_s(out, b.c);
}

static inline void __mult_xfwd1_s(float out[3][4][2][VLENS], pglus_t gp, scs_t* in, int i0, int ifwd) {
  int c, j;
  __attribute__((aligned(64))) g33s_t g;
  __attribute__((aligned(64))) projscs_t a, b;

  __load_su3_s(g.c, gp, i0);
#pragma loop nounroll
  for(c=0;c<3;c++){
    __attribute__((aligned(64))) float v[4][2][VLENS];
    __load_scf_sri_s(ri, s, c, i0, ifwd);
    for(j=0;j<VLENS;j++){
      a.c[c][0][0][j] = v[0][0][j] + v[3][1][j];
      a.c[c][1][0][j] = v[1][0][j] + v[2][1][j];
      a.c[c][0][1][j] = v[0][1][j] - v[3][0][j];
      a.c[c][1][1][j] = v[1][1][j] - v[2][0][j];
    }
  }
  __m_fg_s(b.c, a.c, g.c);
  __s_xf_s(out, b.c);
}
//----------------------------------------------------------------------------------------X backward
static inline void __mult_xbwd0_s(float out[3][4][2][VLENS], pglus_t gp, scs_t* in, int i0, int ibwd) {
  int c, j;
  __attribute__((aligned(64))) g33s_t g;
  __attribute__((aligned(64))) projscs_t a, b;

  __load_su3d1_s(g.c, gp, i0, ibwd);
#pragma loop nounroll
  for(c=0;c<3;c++){
     __attribute__((aligned(64))) float v[4][2][VLENS];
     __load_scb_sri_s(ri, s, c, i0, ibwd);
    for(j=0;j<VLENS;j++){
      a.c[c][0][0][j] = v[0][0][j] - v[3][1][j];
      a.c[c][1][0][j] = v[1][0][j] - v[2][1][j];
      a.c[c][0][1][j] = v[0][1][j] + v[3][0][j];
      a.c[c][1][1][j] = v[1][1][j] + v[2][0][j];
    }
  }
  __m_fg_s(b.c, a.c, g.c);
  __s_xb_s(out, b.c);
}

static inline void __mult_xbwd1_s(float out[3][4][2][VLENS], pglus_t gp, scs_t* in, int i0) {
  int c, j;
  __attribute__((aligned(64))) g33s_t g;
  __attribute__((aligned(64))) projscs_t a, b;

  __load_su3d0_s(g.c, gp, i0);
#pragma loop nounroll
  for(c=0;c<3;c++){
     __attribute__((aligned(64))) float v[4][2][VLENS];
     __load_sc0_sri_s(ri, s, c, i0);
    for(j=0;j<VLENS;j++){
      a.c[c][0][0][j] = v[0][0][j] - v[3][1][j];
      a.c[c][1][0][j] = v[1][0][j] - v[2][1][j];
      a.c[c][0][1][j] = v[0][1][j] + v[3][0][j];
      a.c[c][1][1][j] = v[1][1][j] + v[2][0][j];
    }
  }
  __m_fg_s(b.c, a.c, g.c);
  __s_xb_s(out, b.c);
}
//----------------------------------------------------------------------------------------X forward dagger
static inline void __mult_dag_xfwd0_s(float out[3][4][2][VLENS], pglus_t gp, scs_t* in, int i0) {
  int c, j;
  __attribute__((aligned(64))) g33s_t g;
  __attribute__((aligned(64))) projscs_t a, b;

  __load_su3_s(g.c, gp, i0);
#pragma loop nounroll
  for(c=0;c<3;c++){
    __attribute__((aligned(64))) float v[4][2][VLENS];
    __load_sc0_sri_s(ri, s, c, i0);
    for(j=0;j<VLENS;j++){
      a.c[c][0][0][j] = v[0][0][j] - v[3][1][j];
      a.c[c][1][0][j] = v[1][0][j] - v[2][1][j];
      a.c[c][0][1][j] = v[0][1][j] + v[3][0][j];
      a.c[c][1][1][j] = v[1][1][j] + v[2][0][j];
    }
  }
  __m_fg_s(b.c, a.c, g.c);
  __s_xb_s(out, b.c);
}

static inline void __mult_dag_xfwd1_s(float out[3][4][2][VLENS], pglus_t gp, scs_t* in, int i0, int ifwd) {
  int c, j;
  __attribute__((aligned(64))) g33s_t g;
  __attribute__((aligned(64))) projscs_t a, b;

  __load_su3_s(g.c, gp, i0);
#pragma loop nounroll
  for(c=0;c<3;c++){
    __attribute__((aligned(64))) float v[4][2][VLENS];
    __load_scf_sri_s(ri, s, c, i0, ifwd);
    for(j=0;j<VLENS;j++){
      a.c[c][0][0][j] = v[0][0][j] - v[3][1][j];
      a.c[c][1][0][j] = v[1][0][j] - v[2][1][j];
      a.c[c][0][1][j] = v[0][1][j] + v[3][0][j];
      a.c[c][1][1][j] = v[1][1][j] + v[2][0][j];
    }
  }
  __m_fg_s(b.c, a.c, g.c);
  __s_xb_s(out, b.c);
}
//----------------------------------------------------------------------------------------X backward dagger
static inline void __mult_dag_xbwd0_s(float out[3][4][2][VLENS], pglus_t gp, scs_t* in, int i0, int ibwd) {
  int c, j;
  __attribute__((aligned(64))) g33s_t g;
  __attribute__((aligned(64))) projscs_t a, b;

  __load_su3d1_s(g.c, gp, i0, ibwd);
#pragma loop nounroll
  for(c=0;c<3;c++){
    __attribute__((aligned(64))) float v[4][2][VLENS];
    __load_scb_sri_s(ri, s, c, i0, ibwd);
    for(j=0;j<VLENS;j++){
      a.c[c][0][0][j] = v[0][0][j] + v[3][1][j];
      a.c[c][1][0][j] = v[1][0][j] + v[2][1][j];
      a.c[c][0][1][j] = v[0][1][j] - v[3][0][j];
      a.c[c][1][1][j] = v[1][1][j] - v[2][0][j];
    }
  }
  __m_fg_s(b.c, a.c, g.c);
  __s_xf_s(out, b.c);
}

static inline void __mult_dag_xbwd1_s(float out[3][4][2][VLENS], pglus_t gp, scs_t* in, int i0) {
  int c, j;
  __attribute__((aligned(64))) g33s_t g;
  __attribute__((aligned(64))) projscs_t a, b;

  __load_su3d0_s(g.c, gp, i0);
#pragma loop nounroll
  for(c=0;c<3;c++){
    __attribute__((aligned(64))) float v[4][2][VLENS];
    __load_sc0_sri_s(ri, s, c, i0);
    for(j=0;j<VLENS;j++){
      a.c[c][0][0][j] = v[0][0][j] + v[3][1][j];
      a.c[c][1][0][j] = v[1][0][j] + v[2][1][j];
      a.c[c][0][1][j] = v[0][1][j] - v[3][0][j];
      a.c[c][1][1][j] = v[1][1][j] - v[2][0][j];
    }
  }
  __m_fg_s(b.c, a.c, g.c);
  __s_xf_s(out, b.c);
}
//-------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------Y forward
static inline void __mult_yfwd_s(float out[3][4][2][VLENS], pglus_t gp, scs_t* in, int i0, int ifwd) {
  __attribute__((aligned(64))) g33s_t g;
  __attribute__((aligned(64))) projscs_t a, b;
  __load_su3_s(g.c, gp, i0);
  __l_yf_s(a.c, in, ifwd);
  __m_fg_s(b.c, a.c, g.c);
  __s_yf_s(out, b.c);
}

//----------------------------------------------------------------------------------------Y backward
static inline void __mult_ybwd_s(float out[3][4][2][VLENS], pglus_t gp, scs_t* in, int ibwd) {
  __attribute__((aligned(64))) g33s_t g;
  __attribute__((aligned(64))) projscs_t a, b;
  __load_su3d0_s(g.c, gp, ibwd);
  __l_yb_s(a.c, in, ibwd);
  __m_fg_s(b.c, a.c, g.c);
  __s_yb_s(out, b.c);
}

//----------------------------------------------------------------------------------------Y forward dagger
static inline void __mult_dag_yfwd_s(float out[3][4][2][VLENS], pglus_t gp, scs_t* in, int i0, int ifwd) {
  __attribute__((aligned(64))) g33s_t g;
  __attribute__((aligned(64))) projscs_t a, b;
  __load_su3_s(g.c, gp, i0);
  __l_yb_s(a.c, in, ifwd);
  __m_fg_s(b.c, a.c, g.c);
  __s_yb_s(out, b.c);
}

//----------------------------------------------------------------------------------------Y backward dagger
static inline void __mult_dag_ybwd_s(float out[3][4][2][VLENS], pglus_t gp, scs_t* in, int ibwd) {
  __attribute__((aligned(64))) g33s_t g;
  __attribute__((aligned(64))) projscs_t a, b;
  __load_su3d0_s(g.c, gp, ibwd);
  __l_yf_s(a.c, in, ibwd);
  __m_fg_s(b.c, a.c, g.c);
  __s_yf_s(out, b.c);
}

//-------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------Z forward
static inline void __mult_zfwd_s(float out[3][4][2][VLENS], pglus_t gp, scs_t* in, int i0, int ifwd) {
  __attribute__((aligned(64))) g33s_t g;
  __attribute__((aligned(64))) projscs_t a, b;
  __load_su3_s(g.c, gp, i0);
  __l_zf_s(a.c, in, ifwd);
  __m_fg_s(b.c, a.c, g.c);
  __s_zf_s(out, b.c);
}

//----------------------------------------------------------------------------------------Z backward
static inline void __mult_zbwd_s(float out[3][4][2][VLENS], pglus_t gp, scs_t* in, int ibwd) {
  __attribute__((aligned(64))) g33s_t g;
  __attribute__((aligned(64))) projscs_t a, b;
  __load_su3d0_s(g.c, gp, ibwd);
  __l_zb_s(a.c, in, ibwd);
  __m_fg_s(b.c, a.c, g.c);
  __s_zb_s(out, b.c);
}
//----------------------------------------------------------------------------------------Z forward dagger
static inline void __mult_dag_zfwd_s(float out[3][4][2][VLENS], pglus_t gp, scs_t* in, int i0, int ifwd) {
  __attribute__((aligned(64))) g33s_t g;
  __attribute__((aligned(64))) projscs_t a, b;
  __load_su3_s(g.c, gp, i0);
  __l_zb_s(a.c, in, ifwd);
  __m_fg_s(b.c, a.c, g.c);
  __s_zb_s(out, b.c);
}

//----------------------------------------------------------------------------------------Z backward dagger
static inline void __mult_dag_zbwd_s(float out[3][4][2][VLENS], pglus_t gp, scs_t* in, int ibwd) {
  __attribute__((aligned(64))) g33s_t g;
  __attribute__((aligned(64))) projscs_t a, b;
  __load_su3d0_s(g.c, gp, ibwd);
  __l_zf_s(a.c, in, ibwd);
  __m_fg_s(b.c, a.c, g.c);
  __s_zf_s(out, b.c);
}

//-------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------T forward
static inline void __mult_tfwd_s(float out[3][4][2][VLENS], pglus_t gp, scs_t* in, int i0, int ifwd) {
  __attribute__((aligned(64))) g33s_t g;
  __attribute__((aligned(64))) projscs_t a, b;
  __load_su3_s(g.c, gp, i0);
  __l_tf_s(a.c, in, ifwd);
  __m_fg_s(b.c, a.c, g.c);
  __s_tf_s(out, b.c);
}

//----------------------------------------------------------------------------------------T backward
static inline void __mult_tbwd_s(float out[3][4][2][VLENS], pglus_t gp, scs_t* in, int ibwd) {
  __attribute__((aligned(64))) g33s_t g;
  __attribute__((aligned(64))) projscs_t a, b;
  __load_su3d0_s(g.c, gp, ibwd);
  __l_tb_s(a.c, in, ibwd);
  __m_fg_s(b.c, a.c, g.c);
  __s_tb_s(out, b.c);
}
//----------------------------------------------------------------------------------------T forward dagger
static inline void __mult_dag_tfwd_s(float out[3][4][2][VLENS], pglus_t gp, scs_t* in, int i0, int ifwd) {
  __attribute__((aligned(64))) g33s_t g;
  __attribute__((aligned(64))) projscs_t a, b;
  __load_su3_s(g.c, gp, i0);
  __l_tb_s(a.c, in, ifwd);
  __m_fg_s(b.c, a.c, g.c);
  __s_tb_s(out, b.c);
}

//----------------------------------------------------------------------------------------T backward dagger
static inline void __mult_dag_tbwd_s(float out[3][4][2][VLENS], pglus_t gp, scs_t* in, int ibwd) {
  __attribute__((aligned(64))) g33s_t g;
  __attribute__((aligned(64))) projscs_t a, b;
  __load_su3d0_s(g.c, gp, ibwd);
  __l_tf_s(a.c, in, ibwd);
  __m_fg_s(b.c, a.c, g.c);
  __s_tf_s(out, b.c);
}

//-------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------X forward pos
static inline void __pos_xfwd_s(float out[3][4][2][VLENS], pglus_t gp, projscs1_t a, int i0) {
  __attribute__((aligned(64))) float g[3][3][2];
  __attribute__((aligned(64))) float b[3][2][2];
  int c1, c2, ud, ri;
#if SU3_RECONSTRUCT_S == 18
  for(c2=0;c2<3;c2++){
    for(c1=0;c1<3;c1++){
      g[c2][c1][0] = (*(gp+i0)).c[c2][c1][0][VLENS-1];
      g[c2][c1][1] = (*(gp+i0)).c[c2][c1][1][VLENS-1];
    }
  }
#elif SU3_RECONSTRUCT_S == 12
  for(c2=0;c2<2;c2++){
    for(c1=0;c1<3;c1++){
      g[c2][c1][0] = (*(gp+i0)).c[c2][c1][0][VLENS-1];
      g[c2][c1][1] = (*(gp+i0)).c[c2][c1][1][VLENS-1];
    }
  }
  g[2][0][0] = g[0][1][0]*g[1][2][0]-g[0][1][1]*g[1][2][1]-g[0][2][0]*g[1][1][0]+g[0][2][1]*g[1][1][1];
  g[2][1][0] = g[0][2][0]*g[1][0][0]-g[0][2][1]*g[1][0][1]-g[0][0][0]*g[1][2][0]+g[0][0][1]*g[1][2][1];
  g[2][2][0] = g[0][0][0]*g[1][1][0]-g[0][0][1]*g[1][1][1]-g[0][1][0]*g[1][0][0]+g[0][1][1]*g[1][0][1];
  g[2][0][1] =-g[0][1][0]*g[1][2][1]-g[0][1][1]*g[1][2][0]+g[0][2][0]*g[1][1][1]+g[0][2][1]*g[1][1][0];
  g[2][1][1] =-g[0][2][0]*g[1][0][1]-g[0][2][1]*g[1][0][0]+g[0][0][0]*g[1][2][1]+g[0][0][1]*g[1][2][0];
  g[2][2][1] =-g[0][0][0]*g[1][1][1]-g[0][0][1]*g[1][1][0]+g[0][1][0]*g[1][0][1]+g[0][1][1]*g[1][0][0];
#endif
  for(c1=0;c1<3;c1++){
    for(ud=0;ud<2;ud++){
      for(ri=0;ri<2;ri++){
        b[c1][ud][ri] = 0;
      }
    }
  }
  //#pragma unroll(3)
  for(c2=0;c2<3;c2++){
    for(c1=0;c1<3;c1++){
      for(ud=0;ud<2;ud++){
        b[c1][ud][0] += a.c[c2][ud][0] * g[c2][c1][0] - a.c[c2][ud][1] * g[c2][c1][1];
        b[c1][ud][1] += a.c[c2][ud][0] * g[c2][c1][1] + a.c[c2][ud][1] * g[c2][c1][0];
      }
    }
  }
  for(int c=0;c<3;c++){
    out[c][0][0][VLENS-1] +=  + b[c][0][0];
    out[c][1][0][VLENS-1] +=  + b[c][1][0];
    out[c][2][0][VLENS-1] +=  - b[c][1][1];
    out[c][3][0][VLENS-1] +=  - b[c][0][1];
    out[c][0][1][VLENS-1] +=  + b[c][0][1];
    out[c][1][1][VLENS-1] +=  + b[c][1][1];
    out[c][2][1][VLENS-1] +=  + b[c][1][0];
    out[c][3][1][VLENS-1] +=  + b[c][0][0];
  }
}
//----------------------------------------------------------------------------------------X backward pre
static inline void __pre_xbwd_s(float out[3][2][2], pglus_t gp, scs_t* in, int ibwd) {
  __attribute__((aligned(64))) float g[3][3][2];
  __attribute__((aligned(64))) float a[3][2][2];
  __attribute__((aligned(64))) float v[4][2];

#if SU3_RECONSTRUCT_S == 18
  for(int c2=0;c2<3;c2++){
    for(int c1=0;c1<3;c1++){
      g[c1][c2][0] = (*(gp+ibwd)).c[c2][c1][0][VLENS-1];
      g[c1][c2][1] =-(*(gp+ibwd)).c[c2][c1][1][VLENS-1];
    }
  }
#elif SU3_RECONSTRUCT_S == 12
  for(int c2=0;c2<2;c2++){
    for(int c1=0;c1<3;c1++){
      g[c1][c2][0] = (*(gp+ibwd)).c[c2][c1][0][VLENS-1];
      g[c1][c2][1] =-(*(gp+ibwd)).c[c2][c1][1][VLENS-1];
    }
  }
  g[0][2][0] = g[1][0][0]*g[2][1][0]-g[1][0][1]*g[2][1][1]-g[2][0][0]*g[1][1][0]+g[2][0][1]*g[1][1][1];
  g[1][2][0] = g[2][0][0]*g[0][1][0]-g[2][0][1]*g[0][1][1]-g[0][0][0]*g[2][1][0]+g[0][0][1]*g[2][1][1];
  g[2][2][0] = g[0][0][0]*g[1][1][0]-g[0][0][1]*g[1][1][1]-g[1][0][0]*g[0][1][0]+g[1][0][1]*g[0][1][1];
  g[0][2][1] =-g[1][0][0]*g[2][1][1]-g[1][0][1]*g[2][1][0]+g[2][0][0]*g[1][1][1]+g[2][0][1]*g[1][1][0];
  g[1][2][1] =-g[2][0][0]*g[0][1][1]-g[2][0][1]*g[0][1][0]+g[0][0][0]*g[2][1][1]+g[0][0][1]*g[2][1][0];
  g[2][2][1] =-g[0][0][0]*g[1][1][1]-g[0][0][1]*g[1][1][0]+g[1][0][0]*g[0][1][1]+g[1][0][1]*g[0][1][0];
#endif
  for(int c=0;c<3;c++){
    for (int s=0;s<4;s++) {
      for (int ri=0;ri<2;ri++) {
        v[s][ri]=in[ibwd].c[c][s][ri][VLENS-1];
      }
    }
    a[c][0][0] = v[0][0] - v[3][1];
    a[c][1][0] = v[1][0] - v[2][1];
    a[c][0][1] = v[0][1] + v[3][0];
    a[c][1][1] = v[1][1] + v[2][0];
  }
  for(int c1=0;c1<3;c1++){
    for(int ud=0;ud<2;ud++){
      for(int ri=0;ri<2;ri++){
        out[c1][ud][ri] = 0;
      }
    }
  }
  for(int c2=0;c2<3;c2++){
    for(int c1=0;c1<3;c1++){
      for(int ud=0;ud<2;ud++){
        out[c1][ud][0] += a[c2][ud][0] * g[c2][c1][0] - a[c2][ud][1] * g[c2][c1][1];
        out[c1][ud][1] += a[c2][ud][0] * g[c2][c1][1] + a[c2][ud][1] * g[c2][c1][0];
      }
    }
  }
}
//----------------------------------------------------------------------------------------Y forward pos
static inline void __pos_yfwd_s(float out[3][4][2][VLENS], pglus_t gp, projscs_t a, int i0) {
  __attribute__((aligned(64))) g33s_t g;
  __attribute__((aligned(64))) projscs_t b;
  __load_su3_s(g.c, gp, i0);
  __m_fg_s(b.c, a.c, g.c);
  __s_yf_s(out, b.c);
}
//----------------------------------------------------------------------------------------Y backward pre
static inline void __pre_ybwd_s(float out[3][2][2][VLENS], pglus_t gp, scs_t* in, int ibwd) {
  __attribute__((aligned(64))) g33s_t g;
  __attribute__((aligned(64))) projscs_t a;
  __load_su3d0_s(g.c, gp, ibwd);
  __l_yb_s(a.c, in, ibwd);
  __m_fg_s(out, a.c, g.c);
}
//----------------------------------------------------------------------------------------Z forward pos
static inline void __pos_zfwd_s(float out[3][4][2][VLENS], pglus_t gp, projscs_t a, int i0) {
  __attribute__((aligned(64))) g33s_t g;
  __attribute__((aligned(64))) projscs_t b;
  __load_su3_s(g.c, gp, i0);
  __m_fg_s(b.c, a.c, g.c);
  __s_zf_s(out, b.c);
}
//----------------------------------------------------------------------------------------Z backward pre
static inline void __pre_zbwd_s(float out[3][2][2][VLENS], pglus_t gp, scs_t* in, int ibwd) {
  __attribute__((aligned(64))) g33s_t g;
  __attribute__((aligned(64))) projscs_t a;
  __load_su3d0_s(g.c, gp, ibwd);
  __l_zb_s(a.c, in, ibwd);
  __m_fg_s(out, a.c, g.c);
}

//----------------------------------------------------------------------------------------T forward pos
static inline void __pos_tfwd_s(float out[3][4][2][VLENS], pglus_t gp, projscs_t a, int i0) {
  __attribute__((aligned(64))) g33s_t g;
  __attribute__((aligned(64))) projscs_t b;
  __load_su3_s(g.c, gp, i0);
  __m_fg_s(b.c, a.c, g.c);
  __s_tf_s(out, b.c);
}
//----------------------------------------------------------------------------------------T forward pos
static inline void __pos_tfwd_bc_s(float out[3][4][2][VLENS], pglus_t gp, projscs_t a, int i0, float bc) {
  __attribute__((aligned(64))) g33s_t g;
  __attribute__((aligned(64))) projscs_t b;
  __load_su3_s(g.c, gp, i0);
  __m_fg_s(b.c, a.c, g.c);
  __s_tf_bc_s(out, b.c, bc);
}
//----------------------------------------------------------------------------------------T backward pre
static inline void __pre_tbwd_s(float out[3][2][2][VLENS], pglus_t gp, scs_t* in, int ibwd) {
  __attribute__((aligned(64))) g33s_t g;
  __attribute__((aligned(64))) projscs_t a;
  __load_su3d0_s(g.c, gp, ibwd);
  __l_tb_s(a.c, in, ibwd);
  __m_fg_s(out, a.c, g.c);
}

//------------------------------------------------------------------------------------------------- 
//----------------------------------------------------------------------------------------Y forward pos dagger
static inline void __pos_dag_yfwd_s(float out[3][4][2][VLENS], pglus_t gp, projscs_t a, int i0) {
  __attribute__((aligned(64))) g33s_t g;
  __attribute__((aligned(64))) projscs_t b;
  __load_su3_s(g.c, gp, i0);
  __m_fg_s(b.c, a.c, g.c);
  __s_yb_s(out, b.c);
}
//----------------------------------------------------------------------------------------Y backward pre dagger
static inline void __pre_dag_ybwd_s(float out[3][2][2][VLENS], pglus_t gp, scs_t* in, int ibwd) {
  __attribute__((aligned(64))) g33s_t g;
  __attribute__((aligned(64))) projscs_t a;
  __load_su3d0_s(g.c, gp, ibwd);
  __l_yf_s(a.c, in, ibwd);
  __m_fg_s(out, a.c, g.c);
}
//----------------------------------------------------------------------------------------Z forward pos dagger
static inline void __pos_dag_zfwd_s(float out[3][4][2][VLENS], pglus_t gp, projscs_t a, int i0) {
  __attribute__((aligned(64))) g33s_t g;
  __attribute__((aligned(64))) projscs_t b;
  __load_su3_s(g.c, gp, i0);
  __m_fg_s(b.c, a.c, g.c);
  __s_zb_s(out, b.c);
}
//----------------------------------------------------------------------------------------Z backward pre dagger
static inline void __pre_dag_zbwd_s(float out[3][2][2][VLENS], pglus_t gp, scs_t* in, int ibwd) {
  __attribute__((aligned(64))) g33s_t g;
  __attribute__((aligned(64))) projscs_t a;
  __load_su3d0_s(g.c, gp, ibwd);
  __l_zf_s(a.c, in, ibwd);
  __m_fg_s(out, a.c, g.c);
}

//----------------------------------------------------------------------------------------T forward pos dagger
static inline void __pos_dag_tfwd_s(float out[3][4][2][VLENS], pglus_t gp, projscs_t a, int i0) {
  __attribute__((aligned(64))) g33s_t g;
  __attribute__((aligned(64))) projscs_t b;
  __load_su3_s(g.c, gp, i0);
  __m_fg_s(b.c, a.c, g.c);
  __s_tb_s(out, b.c);
}
//----------------------------------------------------------------------------------------T forward pos dagger
static inline void __pos_dag_tfwd_bc_s(float out[3][4][2][VLENS], pglus_t gp, projscs_t a, int i0, float bc) {
  __attribute__((aligned(64))) g33s_t g;
  __attribute__((aligned(64))) projscs_t b;
  __load_su3_s(g.c, gp, i0);
  __m_fg_s(b.c, a.c, g.c);
  __s_tb_bc_s(out, b.c, bc);
}
//----------------------------------------------------------------------------------------T backward pre dagger
static inline void __pre_dag_tbwd_s(float out[3][2][2][VLENS], pglus_t gp, scs_t* in, int ibwd) {
  __attribute__((aligned(64))) g33s_t g;
  __attribute__((aligned(64))) projscs_t a;
  __load_su3d0_s(g.c, gp, ibwd);
  __l_tf_s(a.c, in, ibwd);
  __m_fg_s(out, a.c, g.c);
}

#define __pre_xbwd_s_def(o,a,g,j){					\
    o.c[0][0][0]=+(a.c[0][0][0][j]-a.c[0][3][1][j])*g.c[0][0][0][j]-(a.c[0][0][1][j]+a.c[0][3][0][j])*-g.c[0][0][1][j]+(a.c[1][0][0][j]-a.c[1][3][1][j])*g.c[0][1][0][j]-(a.c[1][0][1][j]+a.c[1][3][0][j])*-g.c[0][1][1][j]+(a.c[2][0][0][j]-a.c[2][3][1][j])*g.c[0][2][0][j]-(a.c[2][0][1][j]+a.c[2][3][0][j])*-g.c[0][2][1][j]; \
    o.c[0][0][1]=+(a.c[0][0][0][j]-a.c[0][3][1][j])*-g.c[0][0][1][j]+(a.c[0][0][1][j]+a.c[0][3][0][j])*g.c[0][0][0][j]+(a.c[1][0][0][j]-a.c[1][3][1][j])*-g.c[0][1][1][j]+(a.c[1][0][1][j]+a.c[1][3][0][j])*g.c[0][1][0][j]+(a.c[2][0][0][j]-a.c[2][3][1][j])*-g.c[0][2][1][j]+(a.c[2][0][1][j]+a.c[2][3][0][j])*g.c[0][2][0][j]; \
    o.c[0][1][0]=+(a.c[0][1][0][j]-a.c[0][2][1][j])*g.c[0][0][0][j]-(a.c[0][1][1][j]+a.c[0][2][0][j])*-g.c[0][0][1][j]+(a.c[1][1][0][j]-a.c[1][2][1][j])*g.c[0][1][0][j]-(a.c[1][1][1][j]+a.c[1][2][0][j])*-g.c[0][1][1][j]+(a.c[2][1][0][j]-a.c[2][2][1][j])*g.c[0][2][0][j]-(a.c[2][1][1][j]+a.c[2][2][0][j])*-g.c[0][2][1][j]; \
    o.c[0][1][1]=+(a.c[0][1][0][j]-a.c[0][2][1][j])*-g.c[0][0][1][j]+(a.c[0][1][1][j]+a.c[0][2][0][j])*g.c[0][0][0][j]+(a.c[1][1][0][j]-a.c[1][2][1][j])*-g.c[0][1][1][j]+(a.c[1][1][1][j]+a.c[1][2][0][j])*g.c[0][1][0][j]+(a.c[2][1][0][j]-a.c[2][2][1][j])*-g.c[0][2][1][j]+(a.c[2][1][1][j]+a.c[2][2][0][j])*g.c[0][2][0][j]; \
    o.c[1][0][0]=+(a.c[0][0][0][j]-a.c[0][3][1][j])*g.c[1][0][0][j]-(a.c[0][0][1][j]+a.c[0][3][0][j])*-g.c[1][0][1][j]+(a.c[1][0][0][j]-a.c[1][3][1][j])*g.c[1][1][0][j]-(a.c[1][0][1][j]+a.c[1][3][0][j])*-g.c[1][1][1][j]+(a.c[2][0][0][j]-a.c[2][3][1][j])*g.c[1][2][0][j]-(a.c[2][0][1][j]+a.c[2][3][0][j])*-g.c[1][2][1][j]; \
    o.c[1][0][1]=+(a.c[0][0][0][j]-a.c[0][3][1][j])*-g.c[1][0][1][j]+(a.c[0][0][1][j]+a.c[0][3][0][j])*g.c[1][0][0][j]+(a.c[1][0][0][j]-a.c[1][3][1][j])*-g.c[1][1][1][j]+(a.c[1][0][1][j]+a.c[1][3][0][j])*g.c[1][1][0][j]+(a.c[2][0][0][j]-a.c[2][3][1][j])*-g.c[1][2][1][j]+(a.c[2][0][1][j]+a.c[2][3][0][j])*g.c[1][2][0][j]; \
    o.c[1][1][0]=+(a.c[0][1][0][j]-a.c[0][2][1][j])*g.c[1][0][0][j]-(a.c[0][1][1][j]+a.c[0][2][0][j])*-g.c[1][0][1][j]+(a.c[1][1][0][j]-a.c[1][2][1][j])*g.c[1][1][0][j]-(a.c[1][1][1][j]+a.c[1][2][0][j])*-g.c[1][1][1][j]+(a.c[2][1][0][j]-a.c[2][2][1][j])*g.c[1][2][0][j]-(a.c[2][1][1][j]+a.c[2][2][0][j])*-g.c[1][2][1][j]; \
    o.c[1][1][1]=+(a.c[0][1][0][j]-a.c[0][2][1][j])*-g.c[1][0][1][j]+(a.c[0][1][1][j]+a.c[0][2][0][j])*g.c[1][0][0][j]+(a.c[1][1][0][j]-a.c[1][2][1][j])*-g.c[1][1][1][j]+(a.c[1][1][1][j]+a.c[1][2][0][j])*g.c[1][1][0][j]+(a.c[2][1][0][j]-a.c[2][2][1][j])*-g.c[1][2][1][j]+(a.c[2][1][1][j]+a.c[2][2][0][j])*g.c[1][2][0][j]; \
    o.c[2][0][0]=+(a.c[0][0][0][j]-a.c[0][3][1][j])*g.c[2][0][0][j]-(a.c[0][0][1][j]+a.c[0][3][0][j])*-g.c[2][0][1][j]+(a.c[1][0][0][j]-a.c[1][3][1][j])*g.c[2][1][0][j]-(a.c[1][0][1][j]+a.c[1][3][0][j])*-g.c[2][1][1][j]+(a.c[2][0][0][j]-a.c[2][3][1][j])*g.c[2][2][0][j]-(a.c[2][0][1][j]+a.c[2][3][0][j])*-g.c[2][2][1][j]; \
    o.c[2][0][1]=+(a.c[0][0][0][j]-a.c[0][3][1][j])*-g.c[2][0][1][j]+(a.c[0][0][1][j]+a.c[0][3][0][j])*g.c[2][0][0][j]+(a.c[1][0][0][j]-a.c[1][3][1][j])*-g.c[2][1][1][j]+(a.c[1][0][1][j]+a.c[1][3][0][j])*g.c[2][1][0][j]+(a.c[2][0][0][j]-a.c[2][3][1][j])*-g.c[2][2][1][j]+(a.c[2][0][1][j]+a.c[2][3][0][j])*g.c[2][2][0][j]; \
    o.c[2][1][0]=+(a.c[0][1][0][j]-a.c[0][2][1][j])*g.c[2][0][0][j]-(a.c[0][1][1][j]+a.c[0][2][0][j])*-g.c[2][0][1][j]+(a.c[1][1][0][j]-a.c[1][2][1][j])*g.c[2][1][0][j]-(a.c[1][1][1][j]+a.c[1][2][0][j])*-g.c[2][1][1][j]+(a.c[2][1][0][j]-a.c[2][2][1][j])*g.c[2][2][0][j]-(a.c[2][1][1][j]+a.c[2][2][0][j])*-g.c[2][2][1][j]; \
    o.c[2][1][1]=+(a.c[0][1][0][j]-a.c[0][2][1][j])*-g.c[2][0][1][j]+(a.c[0][1][1][j]+a.c[0][2][0][j])*g.c[2][0][0][j]+(a.c[1][1][0][j]-a.c[1][2][1][j])*-g.c[2][1][1][j]+(a.c[1][1][1][j]+a.c[1][2][0][j])*g.c[2][1][0][j]+(a.c[2][1][0][j]-a.c[2][2][1][j])*-g.c[2][2][1][j]+(a.c[2][1][1][j]+a.c[2][2][0][j])*g.c[2][2][0][j]; \
  }

#define __pre_ybwd_s_def(o,a,g){					\
    for(int j=0; j < VLENS; j++) {					\
      o.c[0][0][0][j]=+(a.c[0][0][0][j]+a.c[0][3][0][j])*g.c[0][0][0][j]-(a.c[0][0][1][j]+a.c[0][3][1][j])*-g.c[0][0][1][j]+(a.c[1][0][0][j]+a.c[1][3][0][j])*g.c[0][1][0][j]-(a.c[1][0][1][j]+a.c[1][3][1][j])*-g.c[0][1][1][j]+(a.c[2][0][0][j]+a.c[2][3][0][j])*g.c[0][2][0][j]-(a.c[2][0][1][j]+a.c[2][3][1][j])*-g.c[0][2][1][j]; \
      o.c[0][0][1][j]=+(a.c[0][0][0][j]+a.c[0][3][0][j])*-g.c[0][0][1][j]+(a.c[0][0][1][j]+a.c[0][3][1][j])*g.c[0][0][0][j]+(a.c[1][0][0][j]+a.c[1][3][0][j])*-g.c[0][1][1][j]+(a.c[1][0][1][j]+a.c[1][3][1][j])*g.c[0][1][0][j]+(a.c[2][0][0][j]+a.c[2][3][0][j])*-g.c[0][2][1][j]+(a.c[2][0][1][j]+a.c[2][3][1][j])*g.c[0][2][0][j]; \
      o.c[0][1][0][j]=+(a.c[0][1][0][j]-a.c[0][2][0][j])*g.c[0][0][0][j]-(a.c[0][1][1][j]-a.c[0][2][1][j])*-g.c[0][0][1][j]+(a.c[1][1][0][j]-a.c[1][2][0][j])*g.c[0][1][0][j]-(a.c[1][1][1][j]-a.c[1][2][1][j])*-g.c[0][1][1][j]+(a.c[2][1][0][j]-a.c[2][2][0][j])*g.c[0][2][0][j]-(a.c[2][1][1][j]-a.c[2][2][1][j])*-g.c[0][2][1][j]; \
      o.c[0][1][1][j]=+(a.c[0][1][0][j]-a.c[0][2][0][j])*-g.c[0][0][1][j]+(a.c[0][1][1][j]-a.c[0][2][1][j])*g.c[0][0][0][j]+(a.c[1][1][0][j]-a.c[1][2][0][j])*-g.c[0][1][1][j]+(a.c[1][1][1][j]-a.c[1][2][1][j])*g.c[0][1][0][j]+(a.c[2][1][0][j]-a.c[2][2][0][j])*-g.c[0][2][1][j]+(a.c[2][1][1][j]-a.c[2][2][1][j])*g.c[0][2][0][j]; \
      o.c[1][0][0][j]=+(a.c[0][0][0][j]+a.c[0][3][0][j])*g.c[1][0][0][j]-(a.c[0][0][1][j]+a.c[0][3][1][j])*-g.c[1][0][1][j]+(a.c[1][0][0][j]+a.c[1][3][0][j])*g.c[1][1][0][j]-(a.c[1][0][1][j]+a.c[1][3][1][j])*-g.c[1][1][1][j]+(a.c[2][0][0][j]+a.c[2][3][0][j])*g.c[1][2][0][j]-(a.c[2][0][1][j]+a.c[2][3][1][j])*-g.c[1][2][1][j]; \
      o.c[1][0][1][j]=+(a.c[0][0][0][j]+a.c[0][3][0][j])*-g.c[1][0][1][j]+(a.c[0][0][1][j]+a.c[0][3][1][j])*g.c[1][0][0][j]+(a.c[1][0][0][j]+a.c[1][3][0][j])*-g.c[1][1][1][j]+(a.c[1][0][1][j]+a.c[1][3][1][j])*g.c[1][1][0][j]+(a.c[2][0][0][j]+a.c[2][3][0][j])*-g.c[1][2][1][j]+(a.c[2][0][1][j]+a.c[2][3][1][j])*g.c[1][2][0][j]; \
      o.c[1][1][0][j]=+(a.c[0][1][0][j]-a.c[0][2][0][j])*g.c[1][0][0][j]-(a.c[0][1][1][j]-a.c[0][2][1][j])*-g.c[1][0][1][j]+(a.c[1][1][0][j]-a.c[1][2][0][j])*g.c[1][1][0][j]-(a.c[1][1][1][j]-a.c[1][2][1][j])*-g.c[1][1][1][j]+(a.c[2][1][0][j]-a.c[2][2][0][j])*g.c[1][2][0][j]-(a.c[2][1][1][j]-a.c[2][2][1][j])*-g.c[1][2][1][j]; \
      o.c[1][1][1][j]=+(a.c[0][1][0][j]-a.c[0][2][0][j])*-g.c[1][0][1][j]+(a.c[0][1][1][j]-a.c[0][2][1][j])*g.c[1][0][0][j]+(a.c[1][1][0][j]-a.c[1][2][0][j])*-g.c[1][1][1][j]+(a.c[1][1][1][j]-a.c[1][2][1][j])*g.c[1][1][0][j]+(a.c[2][1][0][j]-a.c[2][2][0][j])*-g.c[1][2][1][j]+(a.c[2][1][1][j]-a.c[2][2][1][j])*g.c[1][2][0][j]; \
      o.c[2][0][0][j]=+(a.c[0][0][0][j]+a.c[0][3][0][j])*g.c[2][0][0][j]-(a.c[0][0][1][j]+a.c[0][3][1][j])*-g.c[2][0][1][j]+(a.c[1][0][0][j]+a.c[1][3][0][j])*g.c[2][1][0][j]-(a.c[1][0][1][j]+a.c[1][3][1][j])*-g.c[2][1][1][j]+(a.c[2][0][0][j]+a.c[2][3][0][j])*g.c[2][2][0][j]-(a.c[2][0][1][j]+a.c[2][3][1][j])*-g.c[2][2][1][j]; \
      o.c[2][0][1][j]=+(a.c[0][0][0][j]+a.c[0][3][0][j])*-g.c[2][0][1][j]+(a.c[0][0][1][j]+a.c[0][3][1][j])*g.c[2][0][0][j]+(a.c[1][0][0][j]+a.c[1][3][0][j])*-g.c[2][1][1][j]+(a.c[1][0][1][j]+a.c[1][3][1][j])*g.c[2][1][0][j]+(a.c[2][0][0][j]+a.c[2][3][0][j])*-g.c[2][2][1][j]+(a.c[2][0][1][j]+a.c[2][3][1][j])*g.c[2][2][0][j]; \
      o.c[2][1][0][j]=+(a.c[0][1][0][j]-a.c[0][2][0][j])*g.c[2][0][0][j]-(a.c[0][1][1][j]-a.c[0][2][1][j])*-g.c[2][0][1][j]+(a.c[1][1][0][j]-a.c[1][2][0][j])*g.c[2][1][0][j]-(a.c[1][1][1][j]-a.c[1][2][1][j])*-g.c[2][1][1][j]+(a.c[2][1][0][j]-a.c[2][2][0][j])*g.c[2][2][0][j]-(a.c[2][1][1][j]-a.c[2][2][1][j])*-g.c[2][2][1][j]; \
      o.c[2][1][1][j]=+(a.c[0][1][0][j]-a.c[0][2][0][j])*-g.c[2][0][1][j]+(a.c[0][1][1][j]-a.c[0][2][1][j])*g.c[2][0][0][j]+(a.c[1][1][0][j]-a.c[1][2][0][j])*-g.c[2][1][1][j]+(a.c[1][1][1][j]-a.c[1][2][1][j])*g.c[2][1][0][j]+(a.c[2][1][0][j]-a.c[2][2][0][j])*-g.c[2][2][1][j]+(a.c[2][1][1][j]-a.c[2][2][1][j])*g.c[2][2][0][j]; \
    }									\
  }

#define __pre_zbwd_s_def(o,a,g){					\
    for(int j=0; j < VLENS; j++) {					\
      o.c[0][0][0][j]=							\
	+(a.c[0][0][0][j]-a.c[0][2][1][j])*g.c[0][0][0][j]-(a.c[0][0][1][j]+a.c[0][2][0][j])*-g.c[0][0][1][j] \
	+(a.c[1][0][0][j]-a.c[1][2][1][j])*g.c[0][1][0][j]-(a.c[1][0][1][j]+a.c[1][2][0][j])*-g.c[0][1][1][j] \
	+(a.c[2][0][0][j]-a.c[2][2][1][j])*g.c[0][2][0][j]-(a.c[2][0][1][j]+a.c[2][2][0][j])*-g.c[0][2][1][j]; \
      o.c[0][0][1][j]=							\
	+(a.c[0][0][0][j]-a.c[0][2][1][j])*-g.c[0][0][1][j]+(a.c[0][0][1][j]+a.c[0][2][0][j])*g.c[0][0][0][j] \
	+(a.c[1][0][0][j]-a.c[1][2][1][j])*-g.c[0][1][1][j]+(a.c[1][0][1][j]+a.c[1][2][0][j])*g.c[0][1][0][j] \
	+(a.c[2][0][0][j]-a.c[2][2][1][j])*-g.c[0][2][1][j]+(a.c[2][0][1][j]+a.c[2][2][0][j])*g.c[0][2][0][j]; \
      o.c[0][1][0][j]=							\
	+(a.c[0][1][0][j]+a.c[0][3][1][j])*g.c[0][0][0][j]-(a.c[0][1][1][j]-a.c[0][3][0][j])*-g.c[0][0][1][j] \
	+(a.c[1][1][0][j]+a.c[1][3][1][j])*g.c[0][1][0][j]-(a.c[1][1][1][j]-a.c[1][3][0][j])*-g.c[0][1][1][j] \
	+(a.c[2][1][0][j]+a.c[2][3][1][j])*g.c[0][2][0][j]-(a.c[2][1][1][j]-a.c[2][3][0][j])*-g.c[0][2][1][j]; \
      o.c[0][1][1][j]=							\
	+(a.c[0][1][0][j]+a.c[0][3][1][j])*-g.c[0][0][1][j]+(a.c[0][1][1][j]-a.c[0][3][0][j])*g.c[0][0][0][j] \
	+(a.c[1][1][0][j]+a.c[1][3][1][j])*-g.c[0][1][1][j]+(a.c[1][1][1][j]-a.c[1][3][0][j])*g.c[0][1][0][j] \
	+(a.c[2][1][0][j]+a.c[2][3][1][j])*-g.c[0][2][1][j]+(a.c[2][1][1][j]-a.c[2][3][0][j])*g.c[0][2][0][j]; \
      o.c[1][0][0][j]=							\
	+(a.c[0][0][0][j]-a.c[0][2][1][j])*g.c[1][0][0][j]-(a.c[0][0][1][j]+a.c[0][2][0][j])*-g.c[1][0][1][j] \
	+(a.c[1][0][0][j]-a.c[1][2][1][j])*g.c[1][1][0][j]-(a.c[1][0][1][j]+a.c[1][2][0][j])*-g.c[1][1][1][j] \
	+(a.c[2][0][0][j]-a.c[2][2][1][j])*g.c[1][2][0][j]-(a.c[2][0][1][j]+a.c[2][2][0][j])*-g.c[1][2][1][j]; \
      o.c[1][0][1][j]=							\
	+(a.c[0][0][0][j]-a.c[0][2][1][j])*-g.c[1][0][1][j]+(a.c[0][0][1][j]+a.c[0][2][0][j])*g.c[1][0][0][j] \
	+(a.c[1][0][0][j]-a.c[1][2][1][j])*-g.c[1][1][1][j]+(a.c[1][0][1][j]+a.c[1][2][0][j])*g.c[1][1][0][j] \
	+(a.c[2][0][0][j]-a.c[2][2][1][j])*-g.c[1][2][1][j]+(a.c[2][0][1][j]+a.c[2][2][0][j])*g.c[1][2][0][j]; \
      o.c[1][1][0][j]=							\
	+(a.c[0][1][0][j]+a.c[0][3][1][j])*g.c[1][0][0][j]-(a.c[0][1][1][j]-a.c[0][3][0][j])*-g.c[1][0][1][j] \
	+(a.c[1][1][0][j]+a.c[1][3][1][j])*g.c[1][1][0][j]-(a.c[1][1][1][j]-a.c[1][3][0][j])*-g.c[1][1][1][j] \
	+(a.c[2][1][0][j]+a.c[2][3][1][j])*g.c[1][2][0][j]-(a.c[2][1][1][j]-a.c[2][3][0][j])*-g.c[1][2][1][j]; \
      o.c[1][1][1][j]=							\
	+(a.c[0][1][0][j]+a.c[0][3][1][j])*-g.c[1][0][1][j]+(a.c[0][1][1][j]-a.c[0][3][0][j])*g.c[1][0][0][j] \
	+(a.c[1][1][0][j]+a.c[1][3][1][j])*-g.c[1][1][1][j]+(a.c[1][1][1][j]-a.c[1][3][0][j])*g.c[1][1][0][j] \
	+(a.c[2][1][0][j]+a.c[2][3][1][j])*-g.c[1][2][1][j]+(a.c[2][1][1][j]-a.c[2][3][0][j])*g.c[1][2][0][j]; \
      o.c[2][0][0][j]=							\
	+(a.c[0][0][0][j]-a.c[0][2][1][j])*g.c[2][0][0][j]-(a.c[0][0][1][j]+a.c[0][2][0][j])*-g.c[2][0][1][j] \
	+(a.c[1][0][0][j]-a.c[1][2][1][j])*g.c[2][1][0][j]-(a.c[1][0][1][j]+a.c[1][2][0][j])*-g.c[2][1][1][j] \
	+(a.c[2][0][0][j]-a.c[2][2][1][j])*g.c[2][2][0][j]-(a.c[2][0][1][j]+a.c[2][2][0][j])*-g.c[2][2][1][j]; \
      o.c[2][0][1][j]=							\
	+(a.c[0][0][0][j]-a.c[0][2][1][j])*-g.c[2][0][1][j]+(a.c[0][0][1][j]+a.c[0][2][0][j])*g.c[2][0][0][j] \
	+(a.c[1][0][0][j]-a.c[1][2][1][j])*-g.c[2][1][1][j]+(a.c[1][0][1][j]+a.c[1][2][0][j])*g.c[2][1][0][j] \
	+(a.c[2][0][0][j]-a.c[2][2][1][j])*-g.c[2][2][1][j]+(a.c[2][0][1][j]+a.c[2][2][0][j])*g.c[2][2][0][j]; \
      o.c[2][1][0][j]=							\
	+(a.c[0][1][0][j]+a.c[0][3][1][j])*g.c[2][0][0][j]-(a.c[0][1][1][j]-a.c[0][3][0][j])*-g.c[2][0][1][j] \
	+(a.c[1][1][0][j]+a.c[1][3][1][j])*g.c[2][1][0][j]-(a.c[1][1][1][j]-a.c[1][3][0][j])*-g.c[2][1][1][j] \
	+(a.c[2][1][0][j]+a.c[2][3][1][j])*g.c[2][2][0][j]-(a.c[2][1][1][j]-a.c[2][3][0][j])*-g.c[2][2][1][j]; \
      o.c[2][1][0][j]=							\
	+(a.c[0][1][0][j]+a.c[0][3][1][j])*g.c[2][0][0][j]-(a.c[0][1][1][j]-a.c[0][3][0][j])*-g.c[2][0][1][j] \
	+(a.c[1][1][0][j]+a.c[1][3][1][j])*g.c[2][1][0][j]-(a.c[1][1][1][j]-a.c[1][3][0][j])*-g.c[2][1][1][j] \
	+(a.c[2][1][0][j]+a.c[2][3][1][j])*g.c[2][2][0][j]-(a.c[2][1][1][j]-a.c[2][3][0][j])*-g.c[2][2][1][j]; \
      o.c[2][1][1][j]=							\
	+(a.c[0][1][0][j]+a.c[0][3][1][j])*-g.c[2][0][1][j]+(a.c[0][1][1][j]-a.c[0][3][0][j])*g.c[2][0][0][j] \
	+(a.c[1][1][0][j]+a.c[1][3][1][j])*-g.c[2][1][1][j]+(a.c[1][1][1][j]-a.c[1][3][0][j])*g.c[2][1][0][j] \
	+(a.c[2][1][0][j]+a.c[2][3][1][j])*-g.c[2][2][1][j]+(a.c[2][1][1][j]-a.c[2][3][0][j])*g.c[2][2][0][j]; \
    }									\
    }

#define __pre_tbwd_s_def(o,a,g){		\
    for(int j=0; j<VLENS; j++) {					\
      o.c[0][0][0][j]=+a.c[0][0][0][j]*g.c[0][0][0][j]-a.c[0][0][1][j]*-g.c[0][0][1][j]+a.c[1][0][0][j]*g.c[0][1][0][j]-a.c[1][0][1][j]*-g.c[0][1][1][j]+a.c[2][0][0][j]*g.c[0][2][0][j]-a.c[2][0][1][j]*-g.c[0][2][1][j]; \
      o.c[0][1][0][j]=+a.c[0][1][0][j]*g.c[0][0][0][j]-a.c[0][1][1][j]*-g.c[0][0][1][j]+a.c[1][1][0][j]*g.c[0][1][0][j]-a.c[1][1][1][j]*-g.c[0][1][1][j]+a.c[2][1][0][j]*g.c[0][2][0][j]-a.c[2][1][1][j]*-g.c[0][2][1][j]; \
      o.c[0][0][1][j]=+a.c[0][0][0][j]*-g.c[0][0][1][j]+a.c[0][0][1][j]*g.c[0][0][0][j]+a.c[1][0][0][j]*-g.c[0][1][1][j]+a.c[1][0][1][j]*g.c[0][1][0][j]+a.c[2][0][0][j]*-g.c[0][2][1][j]+a.c[2][0][1][j]*g.c[0][2][0][j]; \
      o.c[0][1][1][j]=+a.c[0][1][0][j]*-g.c[0][0][1][j]+a.c[0][1][1][j]*g.c[0][0][0][j]+a.c[1][1][0][j]*-g.c[0][1][1][j]+a.c[1][1][1][j]*g.c[0][1][0][j]+a.c[2][1][0][j]*-g.c[0][2][1][j]+a.c[2][1][1][j]*g.c[0][2][0][j]; \
      o.c[1][0][0][j]=+a.c[0][0][0][j]*g.c[1][0][0][j]-a.c[0][0][1][j]*-g.c[1][0][1][j]+a.c[1][0][0][j]*g.c[1][1][0][j]-a.c[1][0][1][j]*-g.c[1][1][1][j]+a.c[2][0][0][j]*g.c[1][2][0][j]-a.c[2][0][1][j]*-g.c[1][2][1][j]; \
      o.c[1][1][0][j]=+a.c[0][1][0][j]*g.c[1][0][0][j]-a.c[0][1][1][j]*-g.c[1][0][1][j]+a.c[1][1][0][j]*g.c[1][1][0][j]-a.c[1][1][1][j]*-g.c[1][1][1][j]+a.c[2][1][0][j]*g.c[1][2][0][j]-a.c[2][1][1][j]*-g.c[1][2][1][j]; \
      o.c[1][0][1][j]=+a.c[0][0][0][j]*-g.c[1][0][1][j]+a.c[0][0][1][j]*g.c[1][0][0][j]+a.c[1][0][0][j]*-g.c[1][1][1][j]+a.c[1][0][1][j]*g.c[1][1][0][j]+a.c[2][0][0][j]*-g.c[1][2][1][j]+a.c[2][0][1][j]*g.c[1][2][0][j]; \
      o.c[1][1][1][j]=+a.c[0][1][0][j]*-g.c[1][0][1][j]+a.c[0][1][1][j]*g.c[1][0][0][j]+a.c[1][1][0][j]*-g.c[1][1][1][j]+a.c[1][1][1][j]*g.c[1][1][0][j]+a.c[2][1][0][j]*-g.c[1][2][1][j]+a.c[2][1][1][j]*g.c[1][2][0][j]; \
      o.c[2][0][0][j]=+a.c[0][0][0][j]*g.c[2][0][0][j]-a.c[0][0][1][j]*-g.c[2][0][1][j]+a.c[1][0][0][j]*g.c[2][1][0][j]-a.c[1][0][1][j]*-g.c[2][1][1][j]+a.c[2][0][0][j]*g.c[2][2][0][j]-a.c[2][0][1][j]*-g.c[2][2][1][j]; \
      o.c[2][1][0][j]=+a.c[0][1][0][j]*g.c[2][0][0][j]-a.c[0][1][1][j]*-g.c[2][0][1][j]+a.c[1][1][0][j]*g.c[2][1][0][j]-a.c[1][1][1][j]*-g.c[2][1][1][j]+a.c[2][1][0][j]*g.c[2][2][0][j]-a.c[2][1][1][j]*-g.c[2][2][1][j]; \
      o.c[2][0][1][j]=+a.c[0][0][0][j]*-g.c[2][0][1][j]+a.c[0][0][1][j]*g.c[2][0][0][j]+a.c[1][0][0][j]*-g.c[2][1][1][j]+a.c[1][0][1][j]*g.c[2][1][0][j]+a.c[2][0][0][j]*-g.c[2][2][1][j]+a.c[2][0][1][j]*g.c[2][2][0][j]; \
      o.c[2][1][1][j]=+a.c[0][1][0][j]*-g.c[2][0][1][j]+a.c[0][1][1][j]*g.c[2][0][0][j]+a.c[1][1][0][j]*-g.c[2][1][1][j]+a.c[1][1][1][j]*g.c[2][1][0][j]+a.c[2][1][0][j]*-g.c[2][2][1][j]+a.c[2][1][1][j]*g.c[2][2][0][j]; \
    }									\
    }

#define __pos_xfwd_s_def(o,a,g,j){					\
    o.c[0][0][0][j] += + (a.c[0][0][0]*g.c[0][0][0][j]-a.c[0][0][1]*g.c[0][0][1][j]+a.c[1][0][0]*g.c[1][0][0][j]-a.c[1][0][1]*g.c[1][0][1][j]+a.c[2][0][0]*g.c[2][0][0][j]-a.c[2][0][1]*g.c[2][0][1][j]); \
    o.c[0][1][0][j] += + (a.c[0][1][0]*g.c[0][0][0][j]-a.c[0][1][1]*g.c[0][0][1][j]+a.c[1][1][0]*g.c[1][0][0][j]-a.c[1][1][1]*g.c[1][0][1][j]+a.c[2][1][0]*g.c[2][0][0][j]-a.c[2][1][1]*g.c[2][0][1][j]); \
    o.c[0][2][0][j] += - (a.c[0][1][0]*g.c[0][0][1][j]+a.c[0][1][1]*g.c[0][0][0][j]+a.c[1][1][0]*g.c[1][0][1][j]+a.c[1][1][1]*g.c[1][0][0][j]+a.c[2][1][0]*g.c[2][0][1][j]+a.c[2][1][1]*g.c[2][0][0][j]); \
    o.c[0][3][0][j] += - (a.c[0][0][0]*g.c[0][0][1][j]+a.c[0][0][1]*g.c[0][0][0][j]+a.c[1][0][0]*g.c[1][0][1][j]+a.c[1][0][1]*g.c[1][0][0][j]+a.c[2][0][0]*g.c[2][0][1][j]+a.c[2][0][1]*g.c[2][0][0][j]); \
    o.c[0][0][1][j] += + (a.c[0][0][0]*g.c[0][0][1][j]+a.c[0][0][1]*g.c[0][0][0][j]+a.c[1][0][0]*g.c[1][0][1][j]+a.c[1][0][1]*g.c[1][0][0][j]+a.c[2][0][0]*g.c[2][0][1][j]+a.c[2][0][1]*g.c[2][0][0][j]); \
    o.c[0][1][1][j] += + (a.c[0][1][0]*g.c[0][0][1][j]+a.c[0][1][1]*g.c[0][0][0][j]+a.c[1][1][0]*g.c[1][0][1][j]+a.c[1][1][1]*g.c[1][0][0][j]+a.c[2][1][0]*g.c[2][0][1][j]+a.c[2][1][1]*g.c[2][0][0][j]); \
    o.c[0][2][1][j] += + (a.c[0][1][0]*g.c[0][0][0][j]-a.c[0][1][1]*g.c[0][0][1][j]+a.c[1][1][0]*g.c[1][0][0][j]-a.c[1][1][1]*g.c[1][0][1][j]+a.c[2][1][0]*g.c[2][0][0][j]-a.c[2][1][1]*g.c[2][0][1][j]); \
    o.c[0][3][1][j] += + (a.c[0][0][0]*g.c[0][0][0][j]-a.c[0][0][1]*g.c[0][0][1][j]+a.c[1][0][0]*g.c[1][0][0][j]-a.c[1][0][1]*g.c[1][0][1][j]+a.c[2][0][0]*g.c[2][0][0][j]-a.c[2][0][1]*g.c[2][0][1][j]); \
    o.c[1][0][0][j] += + (a.c[0][0][0]*g.c[0][1][0][j]-a.c[0][0][1]*g.c[0][1][1][j]+a.c[1][0][0]*g.c[1][1][0][j]-a.c[1][0][1]*g.c[1][1][1][j]+a.c[2][0][0]*g.c[2][1][0][j]-a.c[2][0][1]*g.c[2][1][1][j]); \
    o.c[1][1][0][j] += + (a.c[0][1][0]*g.c[0][1][0][j]-a.c[0][1][1]*g.c[0][1][1][j]+a.c[1][1][0]*g.c[1][1][0][j]-a.c[1][1][1]*g.c[1][1][1][j]+a.c[2][1][0]*g.c[2][1][0][j]-a.c[2][1][1]*g.c[2][1][1][j]); \
    o.c[1][2][0][j] += - (a.c[0][1][0]*g.c[0][1][1][j]+a.c[0][1][1]*g.c[0][1][0][j]+a.c[1][1][0]*g.c[1][1][1][j]+a.c[1][1][1]*g.c[1][1][0][j]+a.c[2][1][0]*g.c[2][1][1][j]+a.c[2][1][1]*g.c[2][1][0][j]); \
    o.c[1][3][0][j] += - (a.c[0][0][0]*g.c[0][1][1][j]+a.c[0][0][1]*g.c[0][1][0][j]+a.c[1][0][0]*g.c[1][1][1][j]+a.c[1][0][1]*g.c[1][1][0][j]+a.c[2][0][0]*g.c[2][1][1][j]+a.c[2][0][1]*g.c[2][1][0][j]); \
    o.c[1][0][1][j] += + (a.c[0][0][0]*g.c[0][1][1][j]+a.c[0][0][1]*g.c[0][1][0][j]+a.c[1][0][0]*g.c[1][1][1][j]+a.c[1][0][1]*g.c[1][1][0][j]+a.c[2][0][0]*g.c[2][1][1][j]+a.c[2][0][1]*g.c[2][1][0][j]); \
    o.c[1][1][1][j] += + (a.c[0][1][0]*g.c[0][1][1][j]+a.c[0][1][1]*g.c[0][1][0][j]+a.c[1][1][0]*g.c[1][1][1][j]+a.c[1][1][1]*g.c[1][1][0][j]+a.c[2][1][0]*g.c[2][1][1][j]+a.c[2][1][1]*g.c[2][1][0][j]); \
    o.c[1][2][1][j] += + (a.c[0][1][0]*g.c[0][1][0][j]-a.c[0][1][1]*g.c[0][1][1][j]+a.c[1][1][0]*g.c[1][1][0][j]-a.c[1][1][1]*g.c[1][1][1][j]+a.c[2][1][0]*g.c[2][1][0][j]-a.c[2][1][1]*g.c[2][1][1][j]); \
    o.c[1][3][1][j] += + (a.c[0][0][0]*g.c[0][1][0][j]-a.c[0][0][1]*g.c[0][1][1][j]+a.c[1][0][0]*g.c[1][1][0][j]-a.c[1][0][1]*g.c[1][1][1][j]+a.c[2][0][0]*g.c[2][1][0][j]-a.c[2][0][1]*g.c[2][1][1][j]); \
    o.c[2][0][0][j] += + (a.c[0][0][0]*g.c[0][2][0][j]-a.c[0][0][1]*g.c[0][2][1][j]+a.c[1][0][0]*g.c[1][2][0][j]-a.c[1][0][1]*g.c[1][2][1][j]+a.c[2][0][0]*g.c[2][2][0][j]-a.c[2][0][1]*g.c[2][2][1][j]); \
    o.c[2][1][0][j] += + (a.c[0][1][0]*g.c[0][2][0][j]-a.c[0][1][1]*g.c[0][2][1][j]+a.c[1][1][0]*g.c[1][2][0][j]-a.c[1][1][1]*g.c[1][2][1][j]+a.c[2][1][0]*g.c[2][2][0][j]-a.c[2][1][1]*g.c[2][2][1][j]); \
    o.c[2][2][0][j] += - (a.c[0][1][0]*g.c[0][2][1][j]+a.c[0][1][1]*g.c[0][2][0][j]+a.c[1][1][0]*g.c[1][2][1][j]+a.c[1][1][1]*g.c[1][2][0][j]+a.c[2][1][0]*g.c[2][2][1][j]+a.c[2][1][1]*g.c[2][2][0][j]); \
    o.c[2][3][0][j] += - (a.c[0][0][0]*g.c[0][2][1][j]+a.c[0][0][1]*g.c[0][2][0][j]+a.c[1][0][0]*g.c[1][2][1][j]+a.c[1][0][1]*g.c[1][2][0][j]+a.c[2][0][0]*g.c[2][2][1][j]+a.c[2][0][1]*g.c[2][2][0][j]); \
    o.c[2][0][1][j] += + (a.c[0][0][0]*g.c[0][2][1][j]+a.c[0][0][1]*g.c[0][2][0][j]+a.c[1][0][0]*g.c[1][2][1][j]+a.c[1][0][1]*g.c[1][2][0][j]+a.c[2][0][0]*g.c[2][2][1][j]+a.c[2][0][1]*g.c[2][2][0][j]); \
    o.c[2][1][1][j] += + (a.c[0][1][0]*g.c[0][2][1][j]+a.c[0][1][1]*g.c[0][2][0][j]+a.c[1][1][0]*g.c[1][2][1][j]+a.c[1][1][1]*g.c[1][2][0][j]+a.c[2][1][0]*g.c[2][2][1][j]+a.c[2][1][1]*g.c[2][2][0][j]); \
    o.c[2][2][1][j] += + (a.c[0][1][0]*g.c[0][2][0][j]-a.c[0][1][1]*g.c[0][2][1][j]+a.c[1][1][0]*g.c[1][2][0][j]-a.c[1][1][1]*g.c[1][2][1][j]+a.c[2][1][0]*g.c[2][2][0][j]-a.c[2][1][1]*g.c[2][2][1][j]); \
    o.c[2][3][1][j] += + (a.c[0][0][0]*g.c[0][2][0][j]-a.c[0][0][1]*g.c[0][2][1][j]+a.c[1][0][0]*g.c[1][2][0][j]-a.c[1][0][1]*g.c[1][2][1][j]+a.c[2][0][0]*g.c[2][2][0][j]-a.c[2][0][1]*g.c[2][2][1][j]); \
    }

#define __pos_yfwd_s_def(o,a,g){\
    for(int j=0; j<VLENS; j++){						\
      o.c[0][0][0][j] += + (a.c[0][0][0][j]*g.c[0][0][0][j]-a.c[0][0][1][j]*g.c[0][0][1][j]+a.c[1][0][0][j]*g.c[1][0][0][j]-a.c[1][0][1][j]*g.c[1][0][1][j]+a.c[2][0][0][j]*g.c[2][0][0][j]-a.c[2][0][1][j]*g.c[2][0][1][j]); \
      o.c[0][1][0][j] += + (a.c[0][1][0][j]*g.c[0][0][0][j]-a.c[0][1][1][j]*g.c[0][0][1][j]+a.c[1][1][0][j]*g.c[1][0][0][j]-a.c[1][1][1][j]*g.c[1][0][1][j]+a.c[2][1][0][j]*g.c[2][0][0][j]-a.c[2][1][1][j]*g.c[2][0][1][j]); \
      o.c[0][2][0][j] += + (a.c[0][1][0][j]*g.c[0][0][0][j]-a.c[0][1][1][j]*g.c[0][0][1][j]+a.c[1][1][0][j]*g.c[1][0][0][j]-a.c[1][1][1][j]*g.c[1][0][1][j]+a.c[2][1][0][j]*g.c[2][0][0][j]-a.c[2][1][1][j]*g.c[2][0][1][j]); \
      o.c[0][3][0][j] += - (a.c[0][0][0][j]*g.c[0][0][0][j]-a.c[0][0][1][j]*g.c[0][0][1][j]+a.c[1][0][0][j]*g.c[1][0][0][j]-a.c[1][0][1][j]*g.c[1][0][1][j]+a.c[2][0][0][j]*g.c[2][0][0][j]-a.c[2][0][1][j]*g.c[2][0][1][j]); \
      o.c[0][0][1][j] += + (a.c[0][0][0][j]*g.c[0][0][1][j]+a.c[0][0][1][j]*g.c[0][0][0][j]+a.c[1][0][0][j]*g.c[1][0][1][j]+a.c[1][0][1][j]*g.c[1][0][0][j]+a.c[2][0][0][j]*g.c[2][0][1][j]+a.c[2][0][1][j]*g.c[2][0][0][j]); \
      o.c[0][1][1][j] += + (a.c[0][1][0][j]*g.c[0][0][1][j]+a.c[0][1][1][j]*g.c[0][0][0][j]+a.c[1][1][0][j]*g.c[1][0][1][j]+a.c[1][1][1][j]*g.c[1][0][0][j]+a.c[2][1][0][j]*g.c[2][0][1][j]+a.c[2][1][1][j]*g.c[2][0][0][j]); \
      o.c[0][2][1][j] += + (a.c[0][1][0][j]*g.c[0][0][1][j]+a.c[0][1][1][j]*g.c[0][0][0][j]+a.c[1][1][0][j]*g.c[1][0][1][j]+a.c[1][1][1][j]*g.c[1][0][0][j]+a.c[2][1][0][j]*g.c[2][0][1][j]+a.c[2][1][1][j]*g.c[2][0][0][j]); \
      o.c[0][3][1][j] += - (a.c[0][0][0][j]*g.c[0][0][1][j]+a.c[0][0][1][j]*g.c[0][0][0][j]+a.c[1][0][0][j]*g.c[1][0][1][j]+a.c[1][0][1][j]*g.c[1][0][0][j]+a.c[2][0][0][j]*g.c[2][0][1][j]+a.c[2][0][1][j]*g.c[2][0][0][j]); \
      o.c[1][0][0][j] += + (a.c[0][0][0][j]*g.c[0][1][0][j]-a.c[0][0][1][j]*g.c[0][1][1][j]+a.c[1][0][0][j]*g.c[1][1][0][j]-a.c[1][0][1][j]*g.c[1][1][1][j]+a.c[2][0][0][j]*g.c[2][1][0][j]-a.c[2][0][1][j]*g.c[2][1][1][j]); \
      o.c[1][1][0][j] += + (a.c[0][1][0][j]*g.c[0][1][0][j]-a.c[0][1][1][j]*g.c[0][1][1][j]+a.c[1][1][0][j]*g.c[1][1][0][j]-a.c[1][1][1][j]*g.c[1][1][1][j]+a.c[2][1][0][j]*g.c[2][1][0][j]-a.c[2][1][1][j]*g.c[2][1][1][j]); \
      o.c[1][2][0][j] += + (a.c[0][1][0][j]*g.c[0][1][0][j]-a.c[0][1][1][j]*g.c[0][1][1][j]+a.c[1][1][0][j]*g.c[1][1][0][j]-a.c[1][1][1][j]*g.c[1][1][1][j]+a.c[2][1][0][j]*g.c[2][1][0][j]-a.c[2][1][1][j]*g.c[2][1][1][j]); \
      o.c[1][3][0][j] += - (a.c[0][0][0][j]*g.c[0][1][0][j]-a.c[0][0][1][j]*g.c[0][1][1][j]+a.c[1][0][0][j]*g.c[1][1][0][j]-a.c[1][0][1][j]*g.c[1][1][1][j]+a.c[2][0][0][j]*g.c[2][1][0][j]-a.c[2][0][1][j]*g.c[2][1][1][j]); \
      o.c[1][0][1][j] += + (a.c[0][0][0][j]*g.c[0][1][1][j]+a.c[0][0][1][j]*g.c[0][1][0][j]+a.c[1][0][0][j]*g.c[1][1][1][j]+a.c[1][0][1][j]*g.c[1][1][0][j]+a.c[2][0][0][j]*g.c[2][1][1][j]+a.c[2][0][1][j]*g.c[2][1][0][j]); \
      o.c[1][1][1][j] += + (a.c[0][1][0][j]*g.c[0][1][1][j]+a.c[0][1][1][j]*g.c[0][1][0][j]+a.c[1][1][0][j]*g.c[1][1][1][j]+a.c[1][1][1][j]*g.c[1][1][0][j]+a.c[2][1][0][j]*g.c[2][1][1][j]+a.c[2][1][1][j]*g.c[2][1][0][j]); \
      o.c[1][2][1][j] += + (a.c[0][1][0][j]*g.c[0][1][1][j]+a.c[0][1][1][j]*g.c[0][1][0][j]+a.c[1][1][0][j]*g.c[1][1][1][j]+a.c[1][1][1][j]*g.c[1][1][0][j]+a.c[2][1][0][j]*g.c[2][1][1][j]+a.c[2][1][1][j]*g.c[2][1][0][j]); \
      o.c[1][3][1][j] += - (a.c[0][0][0][j]*g.c[0][1][1][j]+a.c[0][0][1][j]*g.c[0][1][0][j]+a.c[1][0][0][j]*g.c[1][1][1][j]+a.c[1][0][1][j]*g.c[1][1][0][j]+a.c[2][0][0][j]*g.c[2][1][1][j]+a.c[2][0][1][j]*g.c[2][1][0][j]); \
      o.c[2][0][0][j] += + (a.c[0][0][0][j]*g.c[0][2][0][j]-a.c[0][0][1][j]*g.c[0][2][1][j]+a.c[1][0][0][j]*g.c[1][2][0][j]-a.c[1][0][1][j]*g.c[1][2][1][j]+a.c[2][0][0][j]*g.c[2][2][0][j]-a.c[2][0][1][j]*g.c[2][2][1][j]); \
      o.c[2][1][0][j] += + (a.c[0][1][0][j]*g.c[0][2][0][j]-a.c[0][1][1][j]*g.c[0][2][1][j]+a.c[1][1][0][j]*g.c[1][2][0][j]-a.c[1][1][1][j]*g.c[1][2][1][j]+a.c[2][1][0][j]*g.c[2][2][0][j]-a.c[2][1][1][j]*g.c[2][2][1][j]); \
      o.c[2][2][0][j] += + (a.c[0][1][0][j]*g.c[0][2][0][j]-a.c[0][1][1][j]*g.c[0][2][1][j]+a.c[1][1][0][j]*g.c[1][2][0][j]-a.c[1][1][1][j]*g.c[1][2][1][j]+a.c[2][1][0][j]*g.c[2][2][0][j]-a.c[2][1][1][j]*g.c[2][2][1][j]); \
      o.c[2][3][0][j] += - (a.c[0][0][0][j]*g.c[0][2][0][j]-a.c[0][0][1][j]*g.c[0][2][1][j]+a.c[1][0][0][j]*g.c[1][2][0][j]-a.c[1][0][1][j]*g.c[1][2][1][j]+a.c[2][0][0][j]*g.c[2][2][0][j]-a.c[2][0][1][j]*g.c[2][2][1][j]); \
      o.c[2][0][1][j] += + (a.c[0][0][0][j]*g.c[0][2][1][j]+a.c[0][0][1][j]*g.c[0][2][0][j]+a.c[1][0][0][j]*g.c[1][2][1][j]+a.c[1][0][1][j]*g.c[1][2][0][j]+a.c[2][0][0][j]*g.c[2][2][1][j]+a.c[2][0][1][j]*g.c[2][2][0][j]); \
      o.c[2][1][1][j] += + (a.c[0][1][0][j]*g.c[0][2][1][j]+a.c[0][1][1][j]*g.c[0][2][0][j]+a.c[1][1][0][j]*g.c[1][2][1][j]+a.c[1][1][1][j]*g.c[1][2][0][j]+a.c[2][1][0][j]*g.c[2][2][1][j]+a.c[2][1][1][j]*g.c[2][2][0][j]); \
      o.c[2][2][1][j] += + (a.c[0][1][0][j]*g.c[0][2][1][j]+a.c[0][1][1][j]*g.c[0][2][0][j]+a.c[1][1][0][j]*g.c[1][2][1][j]+a.c[1][1][1][j]*g.c[1][2][0][j]+a.c[2][1][0][j]*g.c[2][2][1][j]+a.c[2][1][1][j]*g.c[2][2][0][j]); \
      o.c[2][3][1][j] += - (a.c[0][0][0][j]*g.c[0][2][1][j]+a.c[0][0][1][j]*g.c[0][2][0][j]+a.c[1][0][0][j]*g.c[1][2][1][j]+a.c[1][0][1][j]*g.c[1][2][0][j]+a.c[2][0][0][j]*g.c[2][2][1][j]+a.c[2][0][1][j]*g.c[2][2][0][j]); \
    }									\
    }

#define __pos_zfwd_s_def(o,a,g){					\
    for(int j=0; j<VLENS; j++){						\
      o.c[0][0][0][j] += + (a.c[0][0][0][j]*g.c[0][0][0][j]-a.c[0][0][1][j]*g.c[0][0][1][j]+a.c[1][0][0][j]*g.c[1][0][0][j]-a.c[1][0][1][j]*g.c[1][0][1][j]+a.c[2][0][0][j]*g.c[2][0][0][j]-a.c[2][0][1][j]*g.c[2][0][1][j]); \
      o.c[0][1][0][j] += + (a.c[0][1][0][j]*g.c[0][0][0][j]-a.c[0][1][1][j]*g.c[0][0][1][j]+a.c[1][1][0][j]*g.c[1][0][0][j]-a.c[1][1][1][j]*g.c[1][0][1][j]+a.c[2][1][0][j]*g.c[2][0][0][j]-a.c[2][1][1][j]*g.c[2][0][1][j]); \
      o.c[0][2][0][j] += - (a.c[0][0][0][j]*g.c[0][0][1][j]+a.c[0][0][1][j]*g.c[0][0][0][j]+a.c[1][0][0][j]*g.c[1][0][1][j]+a.c[1][0][1][j]*g.c[1][0][0][j]+a.c[2][0][0][j]*g.c[2][0][1][j]+a.c[2][0][1][j]*g.c[2][0][0][j]); \
      o.c[0][3][0][j] += + (a.c[0][1][0][j]*g.c[0][0][1][j]+a.c[0][1][1][j]*g.c[0][0][0][j]+a.c[1][1][0][j]*g.c[1][0][1][j]+a.c[1][1][1][j]*g.c[1][0][0][j]+a.c[2][1][0][j]*g.c[2][0][1][j]+a.c[2][1][1][j]*g.c[2][0][0][j]); \
      o.c[0][0][1][j] += + (a.c[0][0][0][j]*g.c[0][0][1][j]+a.c[0][0][1][j]*g.c[0][0][0][j]+a.c[1][0][0][j]*g.c[1][0][1][j]+a.c[1][0][1][j]*g.c[1][0][0][j]+a.c[2][0][0][j]*g.c[2][0][1][j]+a.c[2][0][1][j]*g.c[2][0][0][j]); \
      o.c[0][1][1][j] += + (a.c[0][1][0][j]*g.c[0][0][1][j]+a.c[0][1][1][j]*g.c[0][0][0][j]+a.c[1][1][0][j]*g.c[1][0][1][j]+a.c[1][1][1][j]*g.c[1][0][0][j]+a.c[2][1][0][j]*g.c[2][0][1][j]+a.c[2][1][1][j]*g.c[2][0][0][j]); \
      o.c[0][2][1][j] += + (a.c[0][0][0][j]*g.c[0][0][0][j]-a.c[0][0][1][j]*g.c[0][0][1][j]+a.c[1][0][0][j]*g.c[1][0][0][j]-a.c[1][0][1][j]*g.c[1][0][1][j]+a.c[2][0][0][j]*g.c[2][0][0][j]-a.c[2][0][1][j]*g.c[2][0][1][j]); \
      o.c[0][3][1][j] += - (a.c[0][1][0][j]*g.c[0][0][0][j]-a.c[0][1][1][j]*g.c[0][0][1][j]+a.c[1][1][0][j]*g.c[1][0][0][j]-a.c[1][1][1][j]*g.c[1][0][1][j]+a.c[2][1][0][j]*g.c[2][0][0][j]-a.c[2][1][1][j]*g.c[2][0][1][j]); \
      o.c[1][0][0][j] += + (a.c[0][0][0][j]*g.c[0][1][0][j]-a.c[0][0][1][j]*g.c[0][1][1][j]+a.c[1][0][0][j]*g.c[1][1][0][j]-a.c[1][0][1][j]*g.c[1][1][1][j]+a.c[2][0][0][j]*g.c[2][1][0][j]-a.c[2][0][1][j]*g.c[2][1][1][j]); \
      o.c[1][1][0][j] += + (a.c[0][1][0][j]*g.c[0][1][0][j]-a.c[0][1][1][j]*g.c[0][1][1][j]+a.c[1][1][0][j]*g.c[1][1][0][j]-a.c[1][1][1][j]*g.c[1][1][1][j]+a.c[2][1][0][j]*g.c[2][1][0][j]-a.c[2][1][1][j]*g.c[2][1][1][j]); \
      o.c[1][2][0][j] += - (a.c[0][0][0][j]*g.c[0][1][1][j]+a.c[0][0][1][j]*g.c[0][1][0][j]+a.c[1][0][0][j]*g.c[1][1][1][j]+a.c[1][0][1][j]*g.c[1][1][0][j]+a.c[2][0][0][j]*g.c[2][1][1][j]+a.c[2][0][1][j]*g.c[2][1][0][j]); \
      o.c[1][3][0][j] += + (a.c[0][1][0][j]*g.c[0][1][1][j]+a.c[0][1][1][j]*g.c[0][1][0][j]+a.c[1][1][0][j]*g.c[1][1][1][j]+a.c[1][1][1][j]*g.c[1][1][0][j]+a.c[2][1][0][j]*g.c[2][1][1][j]+a.c[2][1][1][j]*g.c[2][1][0][j]); \
      o.c[1][0][1][j] += + (a.c[0][0][0][j]*g.c[0][1][1][j]+a.c[0][0][1][j]*g.c[0][1][0][j]+a.c[1][0][0][j]*g.c[1][1][1][j]+a.c[1][0][1][j]*g.c[1][1][0][j]+a.c[2][0][0][j]*g.c[2][1][1][j]+a.c[2][0][1][j]*g.c[2][1][0][j]); \
      o.c[1][1][1][j] += + (a.c[0][1][0][j]*g.c[0][1][1][j]+a.c[0][1][1][j]*g.c[0][1][0][j]+a.c[1][1][0][j]*g.c[1][1][1][j]+a.c[1][1][1][j]*g.c[1][1][0][j]+a.c[2][1][0][j]*g.c[2][1][1][j]+a.c[2][1][1][j]*g.c[2][1][0][j]); \
      o.c[1][2][1][j] += + (a.c[0][0][0][j]*g.c[0][1][0][j]-a.c[0][0][1][j]*g.c[0][1][1][j]+a.c[1][0][0][j]*g.c[1][1][0][j]-a.c[1][0][1][j]*g.c[1][1][1][j]+a.c[2][0][0][j]*g.c[2][1][0][j]-a.c[2][0][1][j]*g.c[2][1][1][j]); \
      o.c[1][3][1][j] += - (a.c[0][1][0][j]*g.c[0][1][0][j]-a.c[0][1][1][j]*g.c[0][1][1][j]+a.c[1][1][0][j]*g.c[1][1][0][j]-a.c[1][1][1][j]*g.c[1][1][1][j]+a.c[2][1][0][j]*g.c[2][1][0][j]-a.c[2][1][1][j]*g.c[2][1][1][j]); \
      o.c[2][0][0][j] += + (a.c[0][0][0][j]*g.c[0][2][0][j]-a.c[0][0][1][j]*g.c[0][2][1][j]+a.c[1][0][0][j]*g.c[1][2][0][j]-a.c[1][0][1][j]*g.c[1][2][1][j]+a.c[2][0][0][j]*g.c[2][2][0][j]-a.c[2][0][1][j]*g.c[2][2][1][j]); \
      o.c[2][1][0][j] += + (a.c[0][1][0][j]*g.c[0][2][0][j]-a.c[0][1][1][j]*g.c[0][2][1][j]+a.c[1][1][0][j]*g.c[1][2][0][j]-a.c[1][1][1][j]*g.c[1][2][1][j]+a.c[2][1][0][j]*g.c[2][2][0][j]-a.c[2][1][1][j]*g.c[2][2][1][j]); \
      o.c[2][2][0][j] += - (a.c[0][0][0][j]*g.c[0][2][1][j]+a.c[0][0][1][j]*g.c[0][2][0][j]+a.c[1][0][0][j]*g.c[1][2][1][j]+a.c[1][0][1][j]*g.c[1][2][0][j]+a.c[2][0][0][j]*g.c[2][2][1][j]+a.c[2][0][1][j]*g.c[2][2][0][j]); \
      o.c[2][3][0][j] += + (a.c[0][1][0][j]*g.c[0][2][1][j]+a.c[0][1][1][j]*g.c[0][2][0][j]+a.c[1][1][0][j]*g.c[1][2][1][j]+a.c[1][1][1][j]*g.c[1][2][0][j]+a.c[2][1][0][j]*g.c[2][2][1][j]+a.c[2][1][1][j]*g.c[2][2][0][j]); \
      o.c[2][0][1][j] += + (a.c[0][0][0][j]*g.c[0][2][1][j]+a.c[0][0][1][j]*g.c[0][2][0][j]+a.c[1][0][0][j]*g.c[1][2][1][j]+a.c[1][0][1][j]*g.c[1][2][0][j]+a.c[2][0][0][j]*g.c[2][2][1][j]+a.c[2][0][1][j]*g.c[2][2][0][j]); \
      o.c[2][1][1][j] += + (a.c[0][1][0][j]*g.c[0][2][1][j]+a.c[0][1][1][j]*g.c[0][2][0][j]+a.c[1][1][0][j]*g.c[1][2][1][j]+a.c[1][1][1][j]*g.c[1][2][0][j]+a.c[2][1][0][j]*g.c[2][2][1][j]+a.c[2][1][1][j]*g.c[2][2][0][j]); \
      o.c[2][2][1][j] += + (a.c[0][0][0][j]*g.c[0][2][0][j]-a.c[0][0][1][j]*g.c[0][2][1][j]+a.c[1][0][0][j]*g.c[1][2][0][j]-a.c[1][0][1][j]*g.c[1][2][1][j]+a.c[2][0][0][j]*g.c[2][2][0][j]-a.c[2][0][1][j]*g.c[2][2][1][j]); \
      o.c[2][3][1][j] += - (a.c[0][1][0][j]*g.c[0][2][0][j]-a.c[0][1][1][j]*g.c[0][2][1][j]+a.c[1][1][0][j]*g.c[1][2][0][j]-a.c[1][1][1][j]*g.c[1][2][1][j]+a.c[2][1][0][j]*g.c[2][2][0][j]-a.c[2][1][1][j]*g.c[2][2][1][j]); \
    }									\
    }

#define __pos_tfwd_s_def(o,a,g,s){		\
    for(int j=0; j<VLENS; j++) {					\
      o.c[0][2][0][j] += 2 * (a.c[0][0+s][0][j]*g.c[0][0][0][j]-a.c[0][0+s][1][j]*g.c[0][0][1][j]+a.c[1][0+s][0][j]*g.c[1][0][0][j]-a.c[1][0+s][1][j]*g.c[1][0][1][j]+a.c[2][0+s][0][j]*g.c[2][0][0][j]-a.c[2][0+s][1][j]*g.c[2][0][1][j]); \
      o.c[0][3][0][j] += 2 * (a.c[0][1+s][0][j]*g.c[0][0][0][j]-a.c[0][1+s][1][j]*g.c[0][0][1][j]+a.c[1][1+s][0][j]*g.c[1][0][0][j]-a.c[1][1+s][1][j]*g.c[1][0][1][j]+a.c[2][1+s][0][j]*g.c[2][0][0][j]-a.c[2][1+s][1][j]*g.c[2][0][1][j]); \
      o.c[0][2][1][j] += 2 * (a.c[0][0+s][0][j]*g.c[0][0][1][j]+a.c[0][0+s][1][j]*g.c[0][0][0][j]+a.c[1][0+s][0][j]*g.c[1][0][1][j]+a.c[1][0+s][1][j]*g.c[1][0][0][j]+a.c[2][0+s][0][j]*g.c[2][0][1][j]+a.c[2][0+s][1][j]*g.c[2][0][0][j]); \
      o.c[0][3][1][j] += 2 * (a.c[0][1+s][0][j]*g.c[0][0][1][j]+a.c[0][1+s][1][j]*g.c[0][0][0][j]+a.c[1][1+s][0][j]*g.c[1][0][1][j]+a.c[1][1+s][1][j]*g.c[1][0][0][j]+a.c[2][1+s][0][j]*g.c[2][0][1][j]+a.c[2][1+s][1][j]*g.c[2][0][0][j]); \
      o.c[1][2][0][j] += 2 * (a.c[0][0+s][0][j]*g.c[0][1][0][j]-a.c[0][0+s][1][j]*g.c[0][1][1][j]+a.c[1][0+s][0][j]*g.c[1][1][0][j]-a.c[1][0+s][1][j]*g.c[1][1][1][j]+a.c[2][0+s][0][j]*g.c[2][1][0][j]-a.c[2][0+s][1][j]*g.c[2][1][1][j]); \
      o.c[1][3][0][j] += 2 * (a.c[0][1+s][0][j]*g.c[0][1][0][j]-a.c[0][1+s][1][j]*g.c[0][1][1][j]+a.c[1][1+s][0][j]*g.c[1][1][0][j]-a.c[1][1+s][1][j]*g.c[1][1][1][j]+a.c[2][1+s][0][j]*g.c[2][1][0][j]-a.c[2][1+s][1][j]*g.c[2][1][1][j]); \
      o.c[1][2][1][j] += 2 * (a.c[0][0+s][0][j]*g.c[0][1][1][j]+a.c[0][0+s][1][j]*g.c[0][1][0][j]+a.c[1][0+s][0][j]*g.c[1][1][1][j]+a.c[1][0+s][1][j]*g.c[1][1][0][j]+a.c[2][0+s][0][j]*g.c[2][1][1][j]+a.c[2][0+s][1][j]*g.c[2][1][0][j]); \
      o.c[1][3][1][j] += 2 * (a.c[0][1+s][0][j]*g.c[0][1][1][j]+a.c[0][1+s][1][j]*g.c[0][1][0][j]+a.c[1][1+s][0][j]*g.c[1][1][1][j]+a.c[1][1+s][1][j]*g.c[1][1][0][j]+a.c[2][1+s][0][j]*g.c[2][1][1][j]+a.c[2][1+s][1][j]*g.c[2][1][0][j]); \
      o.c[2][2][0][j] += 2 * (a.c[0][0+s][0][j]*g.c[0][2][0][j]-a.c[0][0+s][1][j]*g.c[0][2][1][j]+a.c[1][0+s][0][j]*g.c[1][2][0][j]-a.c[1][0+s][1][j]*g.c[1][2][1][j]+a.c[2][0+s][0][j]*g.c[2][2][0][j]-a.c[2][0+s][1][j]*g.c[2][2][1][j]); \
      o.c[2][3][0][j] += 2 * (a.c[0][1+s][0][j]*g.c[0][2][0][j]-a.c[0][1+s][1][j]*g.c[0][2][1][j]+a.c[1][1+s][0][j]*g.c[1][2][0][j]-a.c[1][1+s][1][j]*g.c[1][2][1][j]+a.c[2][1+s][0][j]*g.c[2][2][0][j]-a.c[2][1+s][1][j]*g.c[2][2][1][j]); \
      o.c[2][2][1][j] += 2 * (a.c[0][0+s][0][j]*g.c[0][2][1][j]+a.c[0][0+s][1][j]*g.c[0][2][0][j]+a.c[1][0+s][0][j]*g.c[1][2][1][j]+a.c[1][0+s][1][j]*g.c[1][2][0][j]+a.c[2][0+s][0][j]*g.c[2][2][1][j]+a.c[2][0+s][1][j]*g.c[2][2][0][j]); \
      o.c[2][3][1][j] += 2 * (a.c[0][1+s][0][j]*g.c[0][2][1][j]+a.c[0][1+s][1][j]*g.c[0][2][0][j]+a.c[1][1+s][0][j]*g.c[1][2][1][j]+a.c[1][1+s][1][j]*g.c[1][2][0][j]+a.c[2][1+s][0][j]*g.c[2][2][1][j]+a.c[2][1+s][1][j]*g.c[2][2][0][j]); \
    }									\
    }

#define __pos_tfwd_bc_s_def(o,a,g,bc){					\
    for(int j=0; j<VLENS; j++) {					\
      o.c[0][2][0][j] += bc * (a.c[0][0][0][j]*g.c[0][0][0][j]-a.c[0][0][1][j]*g.c[0][0][1][j]+a.c[1][0][0][j]*g.c[1][0][0][j]-a.c[1][0][1][j]*g.c[1][0][1][j]+a.c[2][0][0][j]*g.c[2][0][0][j]-a.c[2][0][1][j]*g.c[2][0][1][j]); \
      o.c[0][3][0][j] += bc * (a.c[0][1][0][j]*g.c[0][0][0][j]-a.c[0][1][1][j]*g.c[0][0][1][j]+a.c[1][1][0][j]*g.c[1][0][0][j]-a.c[1][1][1][j]*g.c[1][0][1][j]+a.c[2][1][0][j]*g.c[2][0][0][j]-a.c[2][1][1][j]*g.c[2][0][1][j]); \
      o.c[0][2][1][j] += bc * (a.c[0][0][0][j]*g.c[0][0][1][j]+a.c[0][0][1][j]*g.c[0][0][0][j]+a.c[1][0][0][j]*g.c[1][0][1][j]+a.c[1][0][1][j]*g.c[1][0][0][j]+a.c[2][0][0][j]*g.c[2][0][1][j]+a.c[2][0][1][j]*g.c[2][0][0][j]); \
      o.c[0][3][1][j] += bc * (a.c[0][1][0][j]*g.c[0][0][1][j]+a.c[0][1][1][j]*g.c[0][0][0][j]+a.c[1][1][0][j]*g.c[1][0][1][j]+a.c[1][1][1][j]*g.c[1][0][0][j]+a.c[2][1][0][j]*g.c[2][0][1][j]+a.c[2][1][1][j]*g.c[2][0][0][j]); \
      o.c[1][2][0][j] += bc * (a.c[0][0][0][j]*g.c[0][1][0][j]-a.c[0][0][1][j]*g.c[0][1][1][j]+a.c[1][0][0][j]*g.c[1][1][0][j]-a.c[1][0][1][j]*g.c[1][1][1][j]+a.c[2][0][0][j]*g.c[2][1][0][j]-a.c[2][0][1][j]*g.c[2][1][1][j]); \
      o.c[1][3][0][j] += bc * (a.c[0][1][0][j]*g.c[0][1][0][j]-a.c[0][1][1][j]*g.c[0][1][1][j]+a.c[1][1][0][j]*g.c[1][1][0][j]-a.c[1][1][1][j]*g.c[1][1][1][j]+a.c[2][1][0][j]*g.c[2][1][0][j]-a.c[2][1][1][j]*g.c[2][1][1][j]); \
      o.c[1][2][1][j] += bc * (a.c[0][0][0][j]*g.c[0][1][1][j]+a.c[0][0][1][j]*g.c[0][1][0][j]+a.c[1][0][0][j]*g.c[1][1][1][j]+a.c[1][0][1][j]*g.c[1][1][0][j]+a.c[2][0][0][j]*g.c[2][1][1][j]+a.c[2][0][1][j]*g.c[2][1][0][j]); \
      o.c[1][3][1][j] += bc * (a.c[0][1][0][j]*g.c[0][1][1][j]+a.c[0][1][1][j]*g.c[0][1][0][j]+a.c[1][1][0][j]*g.c[1][1][1][j]+a.c[1][1][1][j]*g.c[1][1][0][j]+a.c[2][1][0][j]*g.c[2][1][1][j]+a.c[2][1][1][j]*g.c[2][1][0][j]); \
      o.c[2][2][0][j] += bc * (a.c[0][0][0][j]*g.c[0][2][0][j]-a.c[0][0][1][j]*g.c[0][2][1][j]+a.c[1][0][0][j]*g.c[1][2][0][j]-a.c[1][0][1][j]*g.c[1][2][1][j]+a.c[2][0][0][j]*g.c[2][2][0][j]-a.c[2][0][1][j]*g.c[2][2][1][j]); \
      o.c[2][3][0][j] += bc * (a.c[0][1][0][j]*g.c[0][2][0][j]-a.c[0][1][1][j]*g.c[0][2][1][j]+a.c[1][1][0][j]*g.c[1][2][0][j]-a.c[1][1][1][j]*g.c[1][2][1][j]+a.c[2][1][0][j]*g.c[2][2][0][j]-a.c[2][1][1][j]*g.c[2][2][1][j]); \
      o.c[2][2][1][j] += bc * (a.c[0][0][0][j]*g.c[0][2][1][j]+a.c[0][0][1][j]*g.c[0][2][0][j]+a.c[1][0][0][j]*g.c[1][2][1][j]+a.c[1][0][1][j]*g.c[1][2][0][j]+a.c[2][0][0][j]*g.c[2][2][1][j]+a.c[2][0][1][j]*g.c[2][2][0][j]); \
      o.c[2][3][1][j] += bc * (a.c[0][1][0][j]*g.c[0][2][1][j]+a.c[0][1][1][j]*g.c[0][2][0][j]+a.c[1][1][0][j]*g.c[1][2][1][j]+a.c[1][1][1][j]*g.c[1][2][0][j]+a.c[2][1][0][j]*g.c[2][2][1][j]+a.c[2][1][1][j]*g.c[2][2][0][j]); \
    }									\
    }

#define __mult_xfwd_s_def(o,a,g){			\
    for(int j=0;j<VLENS;j++){						\
      o.c[0][0][0][j] += + (a.c[0][0][0][j]*g.c[0][0][0][j]-a.c[0][0][1][j]*g.c[0][0][1][j]+a.c[1][0][0][j]*g.c[1][0][0][j]-a.c[1][0][1][j]*g.c[1][0][1][j]+a.c[2][0][0][j]*g.c[2][0][0][j]-a.c[2][0][1][j]*g.c[2][0][1][j]); \
      o.c[0][1][0][j] += + (a.c[0][1][0][j]*g.c[0][0][0][j]-a.c[0][1][1][j]*g.c[0][0][1][j]+a.c[1][1][0][j]*g.c[1][0][0][j]-a.c[1][1][1][j]*g.c[1][0][1][j]+a.c[2][1][0][j]*g.c[2][0][0][j]-a.c[2][1][1][j]*g.c[2][0][1][j]); \
      o.c[0][2][0][j] += - (a.c[0][1][0][j]*g.c[0][0][1][j]+a.c[0][1][1][j]*g.c[0][0][0][j]+a.c[1][1][0][j]*g.c[1][0][1][j]+a.c[1][1][1][j]*g.c[1][0][0][j]+a.c[2][1][0][j]*g.c[2][0][1][j]+a.c[2][1][1][j]*g.c[2][0][0][j]); \
      o.c[0][3][0][j] += - (a.c[0][0][0][j]*g.c[0][0][1][j]+a.c[0][0][1][j]*g.c[0][0][0][j]+a.c[1][0][0][j]*g.c[1][0][1][j]+a.c[1][0][1][j]*g.c[1][0][0][j]+a.c[2][0][0][j]*g.c[2][0][1][j]+a.c[2][0][1][j]*g.c[2][0][0][j]); \
      o.c[0][0][1][j] += + (a.c[0][0][0][j]*g.c[0][0][1][j]+a.c[0][0][1][j]*g.c[0][0][0][j]+a.c[1][0][0][j]*g.c[1][0][1][j]+a.c[1][0][1][j]*g.c[1][0][0][j]+a.c[2][0][0][j]*g.c[2][0][1][j]+a.c[2][0][1][j]*g.c[2][0][0][j]); \
      o.c[0][1][1][j] += + (a.c[0][1][0][j]*g.c[0][0][1][j]+a.c[0][1][1][j]*g.c[0][0][0][j]+a.c[1][1][0][j]*g.c[1][0][1][j]+a.c[1][1][1][j]*g.c[1][0][0][j]+a.c[2][1][0][j]*g.c[2][0][1][j]+a.c[2][1][1][j]*g.c[2][0][0][j]); \
      o.c[0][2][1][j] += + (a.c[0][1][0][j]*g.c[0][0][0][j]-a.c[0][1][1][j]*g.c[0][0][1][j]+a.c[1][1][0][j]*g.c[1][0][0][j]-a.c[1][1][1][j]*g.c[1][0][1][j]+a.c[2][1][0][j]*g.c[2][0][0][j]-a.c[2][1][1][j]*g.c[2][0][1][j]); \
      o.c[0][3][1][j] += + (a.c[0][0][0][j]*g.c[0][0][0][j]-a.c[0][0][1][j]*g.c[0][0][1][j]+a.c[1][0][0][j]*g.c[1][0][0][j]-a.c[1][0][1][j]*g.c[1][0][1][j]+a.c[2][0][0][j]*g.c[2][0][0][j]-a.c[2][0][1][j]*g.c[2][0][1][j]); \
      o.c[1][0][0][j] += + (a.c[0][0][0][j]*g.c[0][1][0][j]-a.c[0][0][1][j]*g.c[0][1][1][j]+a.c[1][0][0][j]*g.c[1][1][0][j]-a.c[1][0][1][j]*g.c[1][1][1][j]+a.c[2][0][0][j]*g.c[2][1][0][j]-a.c[2][0][1][j]*g.c[2][1][1][j]); \
      o.c[1][1][0][j] += + (a.c[0][1][0][j]*g.c[0][1][0][j]-a.c[0][1][1][j]*g.c[0][1][1][j]+a.c[1][1][0][j]*g.c[1][1][0][j]-a.c[1][1][1][j]*g.c[1][1][1][j]+a.c[2][1][0][j]*g.c[2][1][0][j]-a.c[2][1][1][j]*g.c[2][1][1][j]); \
      o.c[1][2][0][j] += - (a.c[0][1][0][j]*g.c[0][1][1][j]+a.c[0][1][1][j]*g.c[0][1][0][j]+a.c[1][1][0][j]*g.c[1][1][1][j]+a.c[1][1][1][j]*g.c[1][1][0][j]+a.c[2][1][0][j]*g.c[2][1][1][j]+a.c[2][1][1][j]*g.c[2][1][0][j]); \
      o.c[1][3][0][j] += - (a.c[0][0][0][j]*g.c[0][1][1][j]+a.c[0][0][1][j]*g.c[0][1][0][j]+a.c[1][0][0][j]*g.c[1][1][1][j]+a.c[1][0][1][j]*g.c[1][1][0][j]+a.c[2][0][0][j]*g.c[2][1][1][j]+a.c[2][0][1][j]*g.c[2][1][0][j]); \
      o.c[1][0][1][j] += + (a.c[0][0][0][j]*g.c[0][1][1][j]+a.c[0][0][1][j]*g.c[0][1][0][j]+a.c[1][0][0][j]*g.c[1][1][1][j]+a.c[1][0][1][j]*g.c[1][1][0][j]+a.c[2][0][0][j]*g.c[2][1][1][j]+a.c[2][0][1][j]*g.c[2][1][0][j]); \
      o.c[1][1][1][j] += + (a.c[0][1][0][j]*g.c[0][1][1][j]+a.c[0][1][1][j]*g.c[0][1][0][j]+a.c[1][1][0][j]*g.c[1][1][1][j]+a.c[1][1][1][j]*g.c[1][1][0][j]+a.c[2][1][0][j]*g.c[2][1][1][j]+a.c[2][1][1][j]*g.c[2][1][0][j]); \
      o.c[1][2][1][j] += + (a.c[0][1][0][j]*g.c[0][1][0][j]-a.c[0][1][1][j]*g.c[0][1][1][j]+a.c[1][1][0][j]*g.c[1][1][0][j]-a.c[1][1][1][j]*g.c[1][1][1][j]+a.c[2][1][0][j]*g.c[2][1][0][j]-a.c[2][1][1][j]*g.c[2][1][1][j]); \
      o.c[1][3][1][j] += + (a.c[0][0][0][j]*g.c[0][1][0][j]-a.c[0][0][1][j]*g.c[0][1][1][j]+a.c[1][0][0][j]*g.c[1][1][0][j]-a.c[1][0][1][j]*g.c[1][1][1][j]+a.c[2][0][0][j]*g.c[2][1][0][j]-a.c[2][0][1][j]*g.c[2][1][1][j]); \
      o.c[2][0][0][j] += + (a.c[0][0][0][j]*g.c[0][2][0][j]-a.c[0][0][1][j]*g.c[0][2][1][j]+a.c[1][0][0][j]*g.c[1][2][0][j]-a.c[1][0][1][j]*g.c[1][2][1][j]+a.c[2][0][0][j]*g.c[2][2][0][j]-a.c[2][0][1][j]*g.c[2][2][1][j]); \
      o.c[2][1][0][j] += + (a.c[0][1][0][j]*g.c[0][2][0][j]-a.c[0][1][1][j]*g.c[0][2][1][j]+a.c[1][1][0][j]*g.c[1][2][0][j]-a.c[1][1][1][j]*g.c[1][2][1][j]+a.c[2][1][0][j]*g.c[2][2][0][j]-a.c[2][1][1][j]*g.c[2][2][1][j]); \
      o.c[2][2][0][j] += - (a.c[0][1][0][j]*g.c[0][2][1][j]+a.c[0][1][1][j]*g.c[0][2][0][j]+a.c[1][1][0][j]*g.c[1][2][1][j]+a.c[1][1][1][j]*g.c[1][2][0][j]+a.c[2][1][0][j]*g.c[2][2][1][j]+a.c[2][1][1][j]*g.c[2][2][0][j]); \
      o.c[2][3][0][j] += - (a.c[0][0][0][j]*g.c[0][2][1][j]+a.c[0][0][1][j]*g.c[0][2][0][j]+a.c[1][0][0][j]*g.c[1][2][1][j]+a.c[1][0][1][j]*g.c[1][2][0][j]+a.c[2][0][0][j]*g.c[2][2][1][j]+a.c[2][0][1][j]*g.c[2][2][0][j]); \
      o.c[2][0][1][j] += + (a.c[0][0][0][j]*g.c[0][2][1][j]+a.c[0][0][1][j]*g.c[0][2][0][j]+a.c[1][0][0][j]*g.c[1][2][1][j]+a.c[1][0][1][j]*g.c[1][2][0][j]+a.c[2][0][0][j]*g.c[2][2][1][j]+a.c[2][0][1][j]*g.c[2][2][0][j]); \
      o.c[2][1][1][j] += + (a.c[0][1][0][j]*g.c[0][2][1][j]+a.c[0][1][1][j]*g.c[0][2][0][j]+a.c[1][1][0][j]*g.c[1][2][1][j]+a.c[1][1][1][j]*g.c[1][2][0][j]+a.c[2][1][0][j]*g.c[2][2][1][j]+a.c[2][1][1][j]*g.c[2][2][0][j]); \
      o.c[2][2][1][j] += + (a.c[0][1][0][j]*g.c[0][2][0][j]-a.c[0][1][1][j]*g.c[0][2][1][j]+a.c[1][1][0][j]*g.c[1][2][0][j]-a.c[1][1][1][j]*g.c[1][2][1][j]+a.c[2][1][0][j]*g.c[2][2][0][j]-a.c[2][1][1][j]*g.c[2][2][1][j]); \
      o.c[2][3][1][j] += + (a.c[0][0][0][j]*g.c[0][2][0][j]-a.c[0][0][1][j]*g.c[0][2][1][j]+a.c[1][0][0][j]*g.c[1][2][0][j]-a.c[1][0][1][j]*g.c[1][2][1][j]+a.c[2][0][0][j]*g.c[2][2][0][j]-a.c[2][0][1][j]*g.c[2][2][1][j]); \
    }									\
    }

#define __mult_xbwd_s_def(o,a,g){			\
    for(int j=0;j<VLENS;j++){						\
      o.c[0][0][0][j] += + (a.c[0][0][0][j]*g.c[0][0][0][j]-a.c[0][0][1][j]*g.c[0][0][1][j]+a.c[1][0][0][j]*g.c[1][0][0][j]-a.c[1][0][1][j]*g.c[1][0][1][j]+a.c[2][0][0][j]*g.c[2][0][0][j]-a.c[2][0][1][j]*g.c[2][0][1][j]); \
      o.c[0][1][0][j] += + (a.c[0][1][0][j]*g.c[0][0][0][j]-a.c[0][1][1][j]*g.c[0][0][1][j]+a.c[1][1][0][j]*g.c[1][0][0][j]-a.c[1][1][1][j]*g.c[1][0][1][j]+a.c[2][1][0][j]*g.c[2][0][0][j]-a.c[2][1][1][j]*g.c[2][0][1][j]); \
      o.c[0][2][0][j] += + (a.c[0][1][0][j]*g.c[0][0][1][j]+a.c[0][1][1][j]*g.c[0][0][0][j]+a.c[1][1][0][j]*g.c[1][0][1][j]+a.c[1][1][1][j]*g.c[1][0][0][j]+a.c[2][1][0][j]*g.c[2][0][1][j]+a.c[2][1][1][j]*g.c[2][0][0][j]); \
      o.c[0][3][0][j] += + (a.c[0][0][0][j]*g.c[0][0][1][j]+a.c[0][0][1][j]*g.c[0][0][0][j]+a.c[1][0][0][j]*g.c[1][0][1][j]+a.c[1][0][1][j]*g.c[1][0][0][j]+a.c[2][0][0][j]*g.c[2][0][1][j]+a.c[2][0][1][j]*g.c[2][0][0][j]); \
      o.c[0][0][1][j] += + (a.c[0][0][0][j]*g.c[0][0][1][j]+a.c[0][0][1][j]*g.c[0][0][0][j]+a.c[1][0][0][j]*g.c[1][0][1][j]+a.c[1][0][1][j]*g.c[1][0][0][j]+a.c[2][0][0][j]*g.c[2][0][1][j]+a.c[2][0][1][j]*g.c[2][0][0][j]); \
      o.c[0][1][1][j] += + (a.c[0][1][0][j]*g.c[0][0][1][j]+a.c[0][1][1][j]*g.c[0][0][0][j]+a.c[1][1][0][j]*g.c[1][0][1][j]+a.c[1][1][1][j]*g.c[1][0][0][j]+a.c[2][1][0][j]*g.c[2][0][1][j]+a.c[2][1][1][j]*g.c[2][0][0][j]); \
      o.c[0][2][1][j] += - (a.c[0][1][0][j]*g.c[0][0][0][j]-a.c[0][1][1][j]*g.c[0][0][1][j]+a.c[1][1][0][j]*g.c[1][0][0][j]-a.c[1][1][1][j]*g.c[1][0][1][j]+a.c[2][1][0][j]*g.c[2][0][0][j]-a.c[2][1][1][j]*g.c[2][0][1][j]); \
      o.c[0][3][1][j] += - (a.c[0][0][0][j]*g.c[0][0][0][j]-a.c[0][0][1][j]*g.c[0][0][1][j]+a.c[1][0][0][j]*g.c[1][0][0][j]-a.c[1][0][1][j]*g.c[1][0][1][j]+a.c[2][0][0][j]*g.c[2][0][0][j]-a.c[2][0][1][j]*g.c[2][0][1][j]); \
      o.c[1][0][0][j] += + (a.c[0][0][0][j]*g.c[0][1][0][j]-a.c[0][0][1][j]*g.c[0][1][1][j]+a.c[1][0][0][j]*g.c[1][1][0][j]-a.c[1][0][1][j]*g.c[1][1][1][j]+a.c[2][0][0][j]*g.c[2][1][0][j]-a.c[2][0][1][j]*g.c[2][1][1][j]); \
      o.c[1][1][0][j] += + (a.c[0][1][0][j]*g.c[0][1][0][j]-a.c[0][1][1][j]*g.c[0][1][1][j]+a.c[1][1][0][j]*g.c[1][1][0][j]-a.c[1][1][1][j]*g.c[1][1][1][j]+a.c[2][1][0][j]*g.c[2][1][0][j]-a.c[2][1][1][j]*g.c[2][1][1][j]); \
      o.c[1][2][0][j] += + (a.c[0][1][0][j]*g.c[0][1][1][j]+a.c[0][1][1][j]*g.c[0][1][0][j]+a.c[1][1][0][j]*g.c[1][1][1][j]+a.c[1][1][1][j]*g.c[1][1][0][j]+a.c[2][1][0][j]*g.c[2][1][1][j]+a.c[2][1][1][j]*g.c[2][1][0][j]); \
      o.c[1][3][0][j] += + (a.c[0][0][0][j]*g.c[0][1][1][j]+a.c[0][0][1][j]*g.c[0][1][0][j]+a.c[1][0][0][j]*g.c[1][1][1][j]+a.c[1][0][1][j]*g.c[1][1][0][j]+a.c[2][0][0][j]*g.c[2][1][1][j]+a.c[2][0][1][j]*g.c[2][1][0][j]); \
      o.c[1][0][1][j] += + (a.c[0][0][0][j]*g.c[0][1][1][j]+a.c[0][0][1][j]*g.c[0][1][0][j]+a.c[1][0][0][j]*g.c[1][1][1][j]+a.c[1][0][1][j]*g.c[1][1][0][j]+a.c[2][0][0][j]*g.c[2][1][1][j]+a.c[2][0][1][j]*g.c[2][1][0][j]); \
      o.c[1][1][1][j] += + (a.c[0][1][0][j]*g.c[0][1][1][j]+a.c[0][1][1][j]*g.c[0][1][0][j]+a.c[1][1][0][j]*g.c[1][1][1][j]+a.c[1][1][1][j]*g.c[1][1][0][j]+a.c[2][1][0][j]*g.c[2][1][1][j]+a.c[2][1][1][j]*g.c[2][1][0][j]); \
      o.c[1][2][1][j] += - (a.c[0][1][0][j]*g.c[0][1][0][j]-a.c[0][1][1][j]*g.c[0][1][1][j]+a.c[1][1][0][j]*g.c[1][1][0][j]-a.c[1][1][1][j]*g.c[1][1][1][j]+a.c[2][1][0][j]*g.c[2][1][0][j]-a.c[2][1][1][j]*g.c[2][1][1][j]); \
      o.c[1][3][1][j] += - (a.c[0][0][0][j]*g.c[0][1][0][j]-a.c[0][0][1][j]*g.c[0][1][1][j]+a.c[1][0][0][j]*g.c[1][1][0][j]-a.c[1][0][1][j]*g.c[1][1][1][j]+a.c[2][0][0][j]*g.c[2][1][0][j]-a.c[2][0][1][j]*g.c[2][1][1][j]); \
      o.c[2][0][0][j] += + (a.c[0][0][0][j]*g.c[0][2][0][j]-a.c[0][0][1][j]*g.c[0][2][1][j]+a.c[1][0][0][j]*g.c[1][2][0][j]-a.c[1][0][1][j]*g.c[1][2][1][j]+a.c[2][0][0][j]*g.c[2][2][0][j]-a.c[2][0][1][j]*g.c[2][2][1][j]); \
      o.c[2][1][0][j] += + (a.c[0][1][0][j]*g.c[0][2][0][j]-a.c[0][1][1][j]*g.c[0][2][1][j]+a.c[1][1][0][j]*g.c[1][2][0][j]-a.c[1][1][1][j]*g.c[1][2][1][j]+a.c[2][1][0][j]*g.c[2][2][0][j]-a.c[2][1][1][j]*g.c[2][2][1][j]); \
      o.c[2][2][0][j] += + (a.c[0][1][0][j]*g.c[0][2][1][j]+a.c[0][1][1][j]*g.c[0][2][0][j]+a.c[1][1][0][j]*g.c[1][2][1][j]+a.c[1][1][1][j]*g.c[1][2][0][j]+a.c[2][1][0][j]*g.c[2][2][1][j]+a.c[2][1][1][j]*g.c[2][2][0][j]); \
      o.c[2][3][0][j] += + (a.c[0][0][0][j]*g.c[0][2][1][j]+a.c[0][0][1][j]*g.c[0][2][0][j]+a.c[1][0][0][j]*g.c[1][2][1][j]+a.c[1][0][1][j]*g.c[1][2][0][j]+a.c[2][0][0][j]*g.c[2][2][1][j]+a.c[2][0][1][j]*g.c[2][2][0][j]); \
      o.c[2][0][1][j] += + (a.c[0][0][0][j]*g.c[0][2][1][j]+a.c[0][0][1][j]*g.c[0][2][0][j]+a.c[1][0][0][j]*g.c[1][2][1][j]+a.c[1][0][1][j]*g.c[1][2][0][j]+a.c[2][0][0][j]*g.c[2][2][1][j]+a.c[2][0][1][j]*g.c[2][2][0][j]); \
      o.c[2][1][1][j] += + (a.c[0][1][0][j]*g.c[0][2][1][j]+a.c[0][1][1][j]*g.c[0][2][0][j]+a.c[1][1][0][j]*g.c[1][2][1][j]+a.c[1][1][1][j]*g.c[1][2][0][j]+a.c[2][1][0][j]*g.c[2][2][1][j]+a.c[2][1][1][j]*g.c[2][2][0][j]); \
      o.c[2][2][1][j] += - (a.c[0][1][0][j]*g.c[0][2][0][j]-a.c[0][1][1][j]*g.c[0][2][1][j]+a.c[1][1][0][j]*g.c[1][2][0][j]-a.c[1][1][1][j]*g.c[1][2][1][j]+a.c[2][1][0][j]*g.c[2][2][0][j]-a.c[2][1][1][j]*g.c[2][2][1][j]); \
      o.c[2][3][1][j] += - (a.c[0][0][0][j]*g.c[0][2][0][j]-a.c[0][0][1][j]*g.c[0][2][1][j]+a.c[1][0][0][j]*g.c[1][2][0][j]-a.c[1][0][1][j]*g.c[1][2][1][j]+a.c[2][0][0][j]*g.c[2][2][0][j]-a.c[2][0][1][j]*g.c[2][2][1][j]); \
    }									\
    }

#define __mult_yfwd_s_def(o,a,g){			\
    for(int j=0;j<VLENS;j++){						\
      o.c[0][0][0][j] += + ((a.c[0][0][0][j]-a.c[0][3][0][j])*g.c[0][0][0][j]-(a.c[0][0][1][j]-a.c[0][3][1][j])*g.c[0][0][1][j]+(a.c[1][0][0][j]-a.c[1][3][0][j])*g.c[1][0][0][j]-(a.c[1][0][1][j]-a.c[1][3][1][j])*g.c[1][0][1][j]+(a.c[2][0][0][j]-a.c[2][3][0][j])*g.c[2][0][0][j]-(a.c[2][0][1][j]-a.c[2][3][1][j])*g.c[2][0][1][j]); \
      o.c[0][1][0][j] += + ((a.c[0][1][0][j]+a.c[0][2][0][j])*g.c[0][0][0][j]-(a.c[0][1][1][j]+a.c[0][2][1][j])*g.c[0][0][1][j]+(a.c[1][1][0][j]+a.c[1][2][0][j])*g.c[1][0][0][j]-(a.c[1][1][1][j]+a.c[1][2][1][j])*g.c[1][0][1][j]+(a.c[2][1][0][j]+a.c[2][2][0][j])*g.c[2][0][0][j]-(a.c[2][1][1][j]+a.c[2][2][1][j])*g.c[2][0][1][j]); \
      o.c[0][2][0][j] += + ((a.c[0][1][0][j]+a.c[0][2][0][j])*g.c[0][0][0][j]-(a.c[0][1][1][j]+a.c[0][2][1][j])*g.c[0][0][1][j]+(a.c[1][1][0][j]+a.c[1][2][0][j])*g.c[1][0][0][j]-(a.c[1][1][1][j]+a.c[1][2][1][j])*g.c[1][0][1][j]+(a.c[2][1][0][j]+a.c[2][2][0][j])*g.c[2][0][0][j]-(a.c[2][1][1][j]+a.c[2][2][1][j])*g.c[2][0][1][j]); \
      o.c[0][3][0][j] += - ((a.c[0][0][0][j]-a.c[0][3][0][j])*g.c[0][0][0][j]-(a.c[0][0][1][j]-a.c[0][3][1][j])*g.c[0][0][1][j]+(a.c[1][0][0][j]-a.c[1][3][0][j])*g.c[1][0][0][j]-(a.c[1][0][1][j]-a.c[1][3][1][j])*g.c[1][0][1][j]+(a.c[2][0][0][j]-a.c[2][3][0][j])*g.c[2][0][0][j]-(a.c[2][0][1][j]-a.c[2][3][1][j])*g.c[2][0][1][j]); \
      o.c[0][0][1][j] += + ((a.c[0][0][0][j]-a.c[0][3][0][j])*g.c[0][0][1][j]+(a.c[0][0][1][j]-a.c[0][3][1][j])*g.c[0][0][0][j]+(a.c[1][0][0][j]-a.c[1][3][0][j])*g.c[1][0][1][j]+(a.c[1][0][1][j]-a.c[1][3][1][j])*g.c[1][0][0][j]+(a.c[2][0][0][j]-a.c[2][3][0][j])*g.c[2][0][1][j]+(a.c[2][0][1][j]-a.c[2][3][1][j])*g.c[2][0][0][j]); \
      o.c[0][1][1][j] += + ((a.c[0][1][0][j]+a.c[0][2][0][j])*g.c[0][0][1][j]+(a.c[0][1][1][j]+a.c[0][2][1][j])*g.c[0][0][0][j]+(a.c[1][1][0][j]+a.c[1][2][0][j])*g.c[1][0][1][j]+(a.c[1][1][1][j]+a.c[1][2][1][j])*g.c[1][0][0][j]+(a.c[2][1][0][j]+a.c[2][2][0][j])*g.c[2][0][1][j]+(a.c[2][1][1][j]+a.c[2][2][1][j])*g.c[2][0][0][j]); \
      o.c[0][2][1][j] += + ((a.c[0][1][0][j]+a.c[0][2][0][j])*g.c[0][0][1][j]+(a.c[0][1][1][j]+a.c[0][2][1][j])*g.c[0][0][0][j]+(a.c[1][1][0][j]+a.c[1][2][0][j])*g.c[1][0][1][j]+(a.c[1][1][1][j]+a.c[1][2][1][j])*g.c[1][0][0][j]+(a.c[2][1][0][j]+a.c[2][2][0][j])*g.c[2][0][1][j]+(a.c[2][1][1][j]+a.c[2][2][1][j])*g.c[2][0][0][j]); \
      o.c[0][3][1][j] += - ((a.c[0][0][0][j]-a.c[0][3][0][j])*g.c[0][0][1][j]+(a.c[0][0][1][j]-a.c[0][3][1][j])*g.c[0][0][0][j]+(a.c[1][0][0][j]-a.c[1][3][0][j])*g.c[1][0][1][j]+(a.c[1][0][1][j]-a.c[1][3][1][j])*g.c[1][0][0][j]+(a.c[2][0][0][j]-a.c[2][3][0][j])*g.c[2][0][1][j]+(a.c[2][0][1][j]-a.c[2][3][1][j])*g.c[2][0][0][j]); \
      o.c[1][0][0][j] += + ((a.c[0][0][0][j]-a.c[0][3][0][j])*g.c[0][1][0][j]-(a.c[0][0][1][j]-a.c[0][3][1][j])*g.c[0][1][1][j]+(a.c[1][0][0][j]-a.c[1][3][0][j])*g.c[1][1][0][j]-(a.c[1][0][1][j]-a.c[1][3][1][j])*g.c[1][1][1][j]+(a.c[2][0][0][j]-a.c[2][3][0][j])*g.c[2][1][0][j]-(a.c[2][0][1][j]-a.c[2][3][1][j])*g.c[2][1][1][j]); \
      o.c[1][1][0][j] += + ((a.c[0][1][0][j]+a.c[0][2][0][j])*g.c[0][1][0][j]-(a.c[0][1][1][j]+a.c[0][2][1][j])*g.c[0][1][1][j]+(a.c[1][1][0][j]+a.c[1][2][0][j])*g.c[1][1][0][j]-(a.c[1][1][1][j]+a.c[1][2][1][j])*g.c[1][1][1][j]+(a.c[2][1][0][j]+a.c[2][2][0][j])*g.c[2][1][0][j]-(a.c[2][1][1][j]+a.c[2][2][1][j])*g.c[2][1][1][j]); \
      o.c[1][2][0][j] += + ((a.c[0][1][0][j]+a.c[0][2][0][j])*g.c[0][1][0][j]-(a.c[0][1][1][j]+a.c[0][2][1][j])*g.c[0][1][1][j]+(a.c[1][1][0][j]+a.c[1][2][0][j])*g.c[1][1][0][j]-(a.c[1][1][1][j]+a.c[1][2][1][j])*g.c[1][1][1][j]+(a.c[2][1][0][j]+a.c[2][2][0][j])*g.c[2][1][0][j]-(a.c[2][1][1][j]+a.c[2][2][1][j])*g.c[2][1][1][j]); \
      o.c[1][3][0][j] += - ((a.c[0][0][0][j]-a.c[0][3][0][j])*g.c[0][1][0][j]-(a.c[0][0][1][j]-a.c[0][3][1][j])*g.c[0][1][1][j]+(a.c[1][0][0][j]-a.c[1][3][0][j])*g.c[1][1][0][j]-(a.c[1][0][1][j]-a.c[1][3][1][j])*g.c[1][1][1][j]+(a.c[2][0][0][j]-a.c[2][3][0][j])*g.c[2][1][0][j]-(a.c[2][0][1][j]-a.c[2][3][1][j])*g.c[2][1][1][j]); \
      o.c[1][0][1][j] += + ((a.c[0][0][0][j]-a.c[0][3][0][j])*g.c[0][1][1][j]+(a.c[0][0][1][j]-a.c[0][3][1][j])*g.c[0][1][0][j]+(a.c[1][0][0][j]-a.c[1][3][0][j])*g.c[1][1][1][j]+(a.c[1][0][1][j]-a.c[1][3][1][j])*g.c[1][1][0][j]+(a.c[2][0][0][j]-a.c[2][3][0][j])*g.c[2][1][1][j]+(a.c[2][0][1][j]-a.c[2][3][1][j])*g.c[2][1][0][j]); \
      o.c[1][1][1][j] += + ((a.c[0][1][0][j]+a.c[0][2][0][j])*g.c[0][1][1][j]+(a.c[0][1][1][j]+a.c[0][2][1][j])*g.c[0][1][0][j]+(a.c[1][1][0][j]+a.c[1][2][0][j])*g.c[1][1][1][j]+(a.c[1][1][1][j]+a.c[1][2][1][j])*g.c[1][1][0][j]+(a.c[2][1][0][j]+a.c[2][2][0][j])*g.c[2][1][1][j]+(a.c[2][1][1][j]+a.c[2][2][1][j])*g.c[2][1][0][j]); \
      o.c[1][2][1][j] += + ((a.c[0][1][0][j]+a.c[0][2][0][j])*g.c[0][1][1][j]+(a.c[0][1][1][j]+a.c[0][2][1][j])*g.c[0][1][0][j]+(a.c[1][1][0][j]+a.c[1][2][0][j])*g.c[1][1][1][j]+(a.c[1][1][1][j]+a.c[1][2][1][j])*g.c[1][1][0][j]+(a.c[2][1][0][j]+a.c[2][2][0][j])*g.c[2][1][1][j]+(a.c[2][1][1][j]+a.c[2][2][1][j])*g.c[2][1][0][j]); \
      o.c[1][3][1][j] += - ((a.c[0][0][0][j]-a.c[0][3][0][j])*g.c[0][1][1][j]+(a.c[0][0][1][j]-a.c[0][3][1][j])*g.c[0][1][0][j]+(a.c[1][0][0][j]-a.c[1][3][0][j])*g.c[1][1][1][j]+(a.c[1][0][1][j]-a.c[1][3][1][j])*g.c[1][1][0][j]+(a.c[2][0][0][j]-a.c[2][3][0][j])*g.c[2][1][1][j]+(a.c[2][0][1][j]-a.c[2][3][1][j])*g.c[2][1][0][j]); \
      o.c[2][0][0][j] += + ((a.c[0][0][0][j]-a.c[0][3][0][j])*g.c[0][2][0][j]-(a.c[0][0][1][j]-a.c[0][3][1][j])*g.c[0][2][1][j]+(a.c[1][0][0][j]-a.c[1][3][0][j])*g.c[1][2][0][j]-(a.c[1][0][1][j]-a.c[1][3][1][j])*g.c[1][2][1][j]+(a.c[2][0][0][j]-a.c[2][3][0][j])*g.c[2][2][0][j]-(a.c[2][0][1][j]-a.c[2][3][1][j])*g.c[2][2][1][j]); \
      o.c[2][1][0][j] += + ((a.c[0][1][0][j]+a.c[0][2][0][j])*g.c[0][2][0][j]-(a.c[0][1][1][j]+a.c[0][2][1][j])*g.c[0][2][1][j]+(a.c[1][1][0][j]+a.c[1][2][0][j])*g.c[1][2][0][j]-(a.c[1][1][1][j]+a.c[1][2][1][j])*g.c[1][2][1][j]+(a.c[2][1][0][j]+a.c[2][2][0][j])*g.c[2][2][0][j]-(a.c[2][1][1][j]+a.c[2][2][1][j])*g.c[2][2][1][j]); \
      o.c[2][2][0][j] += + ((a.c[0][1][0][j]+a.c[0][2][0][j])*g.c[0][2][0][j]-(a.c[0][1][1][j]+a.c[0][2][1][j])*g.c[0][2][1][j]+(a.c[1][1][0][j]+a.c[1][2][0][j])*g.c[1][2][0][j]-(a.c[1][1][1][j]+a.c[1][2][1][j])*g.c[1][2][1][j]+(a.c[2][1][0][j]+a.c[2][2][0][j])*g.c[2][2][0][j]-(a.c[2][1][1][j]+a.c[2][2][1][j])*g.c[2][2][1][j]); \
      o.c[2][3][0][j] += - ((a.c[0][0][0][j]-a.c[0][3][0][j])*g.c[0][2][0][j]-(a.c[0][0][1][j]-a.c[0][3][1][j])*g.c[0][2][1][j]+(a.c[1][0][0][j]-a.c[1][3][0][j])*g.c[1][2][0][j]-(a.c[1][0][1][j]-a.c[1][3][1][j])*g.c[1][2][1][j]+(a.c[2][0][0][j]-a.c[2][3][0][j])*g.c[2][2][0][j]-(a.c[2][0][1][j]-a.c[2][3][1][j])*g.c[2][2][1][j]); \
      o.c[2][0][1][j] += + ((a.c[0][0][0][j]-a.c[0][3][0][j])*g.c[0][2][1][j]+(a.c[0][0][1][j]-a.c[0][3][1][j])*g.c[0][2][0][j]+(a.c[1][0][0][j]-a.c[1][3][0][j])*g.c[1][2][1][j]+(a.c[1][0][1][j]-a.c[1][3][1][j])*g.c[1][2][0][j]+(a.c[2][0][0][j]-a.c[2][3][0][j])*g.c[2][2][1][j]+(a.c[2][0][1][j]-a.c[2][3][1][j])*g.c[2][2][0][j]); \
      o.c[2][1][1][j] += + ((a.c[0][1][0][j]+a.c[0][2][0][j])*g.c[0][2][1][j]+(a.c[0][1][1][j]+a.c[0][2][1][j])*g.c[0][2][0][j]+(a.c[1][1][0][j]+a.c[1][2][0][j])*g.c[1][2][1][j]+(a.c[1][1][1][j]+a.c[1][2][1][j])*g.c[1][2][0][j]+(a.c[2][1][0][j]+a.c[2][2][0][j])*g.c[2][2][1][j]+(a.c[2][1][1][j]+a.c[2][2][1][j])*g.c[2][2][0][j]); \
      o.c[2][2][1][j] += + ((a.c[0][1][0][j]+a.c[0][2][0][j])*g.c[0][2][1][j]+(a.c[0][1][1][j]+a.c[0][2][1][j])*g.c[0][2][0][j]+(a.c[1][1][0][j]+a.c[1][2][0][j])*g.c[1][2][1][j]+(a.c[1][1][1][j]+a.c[1][2][1][j])*g.c[1][2][0][j]+(a.c[2][1][0][j]+a.c[2][2][0][j])*g.c[2][2][1][j]+(a.c[2][1][1][j]+a.c[2][2][1][j])*g.c[2][2][0][j]); \
      o.c[2][3][1][j] += - ((a.c[0][0][0][j]-a.c[0][3][0][j])*g.c[0][2][1][j]+(a.c[0][0][1][j]-a.c[0][3][1][j])*g.c[0][2][0][j]+(a.c[1][0][0][j]-a.c[1][3][0][j])*g.c[1][2][1][j]+(a.c[1][0][1][j]-a.c[1][3][1][j])*g.c[1][2][0][j]+(a.c[2][0][0][j]-a.c[2][3][0][j])*g.c[2][2][1][j]+(a.c[2][0][1][j]-a.c[2][3][1][j])*g.c[2][2][0][j]); \
    }									\
    }

#define __mult_ybwd_s_def(o,a,g){			\
    for(int j=0;j<VLENS;j++){						\
      o.c[0][0][0][j] += + ((a.c[0][0][0][j]+a.c[0][3][0][j])*g.c[0][0][0][j]-(a.c[0][0][1][j]+a.c[0][3][1][j])*-g.c[0][0][1][j]+(a.c[1][0][0][j]+a.c[1][3][0][j])*g.c[0][1][0][j]-(a.c[1][0][1][j]+a.c[1][3][1][j])*-g.c[0][1][1][j]+(a.c[2][0][0][j]+a.c[2][3][0][j])*g.c[0][2][0][j]-(a.c[2][0][1][j]+a.c[2][3][1][j])*-g.c[0][2][1][j]); \
      o.c[0][1][0][j] += + ((a.c[0][1][0][j]-a.c[0][2][0][j])*g.c[0][0][0][j]-(a.c[0][1][1][j]-a.c[0][2][1][j])*-g.c[0][0][1][j]+(a.c[1][1][0][j]-a.c[1][2][0][j])*g.c[0][1][0][j]-(a.c[1][1][1][j]-a.c[1][2][1][j])*-g.c[0][1][1][j]+(a.c[2][1][0][j]-a.c[2][2][0][j])*g.c[0][2][0][j]-(a.c[2][1][1][j]-a.c[2][2][1][j])*-g.c[0][2][1][j]); \
      o.c[0][2][0][j] += - ((a.c[0][1][0][j]-a.c[0][2][0][j])*g.c[0][0][0][j]-(a.c[0][1][1][j]-a.c[0][2][1][j])*-g.c[0][0][1][j]+(a.c[1][1][0][j]-a.c[1][2][0][j])*g.c[0][1][0][j]-(a.c[1][1][1][j]-a.c[1][2][1][j])*-g.c[0][1][1][j]+(a.c[2][1][0][j]-a.c[2][2][0][j])*g.c[0][2][0][j]-(a.c[2][1][1][j]-a.c[2][2][1][j])*-g.c[0][2][1][j]); \
      o.c[0][3][0][j] += + ((a.c[0][0][0][j]+a.c[0][3][0][j])*g.c[0][0][0][j]-(a.c[0][0][1][j]+a.c[0][3][1][j])*-g.c[0][0][1][j]+(a.c[1][0][0][j]+a.c[1][3][0][j])*g.c[0][1][0][j]-(a.c[1][0][1][j]+a.c[1][3][1][j])*-g.c[0][1][1][j]+(a.c[2][0][0][j]+a.c[2][3][0][j])*g.c[0][2][0][j]-(a.c[2][0][1][j]+a.c[2][3][1][j])*-g.c[0][2][1][j]); \
      o.c[0][0][1][j] += + ((a.c[0][0][0][j]+a.c[0][3][0][j])*-g.c[0][0][1][j]+(a.c[0][0][1][j]+a.c[0][3][1][j])*g.c[0][0][0][j]+(a.c[1][0][0][j]+a.c[1][3][0][j])*-g.c[0][1][1][j]+(a.c[1][0][1][j]+a.c[1][3][1][j])*g.c[0][1][0][j]+(a.c[2][0][0][j]+a.c[2][3][0][j])*-g.c[0][2][1][j]+(a.c[2][0][1][j]+a.c[2][3][1][j])*g.c[0][2][0][j]); \
      o.c[0][1][1][j] += + ((a.c[0][1][0][j]-a.c[0][2][0][j])*-g.c[0][0][1][j]+(a.c[0][1][1][j]-a.c[0][2][1][j])*g.c[0][0][0][j]+(a.c[1][1][0][j]-a.c[1][2][0][j])*-g.c[0][1][1][j]+(a.c[1][1][1][j]-a.c[1][2][1][j])*g.c[0][1][0][j]+(a.c[2][1][0][j]-a.c[2][2][0][j])*-g.c[0][2][1][j]+(a.c[2][1][1][j]-a.c[2][2][1][j])*g.c[0][2][0][j]); \
      o.c[0][2][1][j] += - ((a.c[0][1][0][j]-a.c[0][2][0][j])*-g.c[0][0][1][j]+(a.c[0][1][1][j]-a.c[0][2][1][j])*g.c[0][0][0][j]+(a.c[1][1][0][j]-a.c[1][2][0][j])*-g.c[0][1][1][j]+(a.c[1][1][1][j]-a.c[1][2][1][j])*g.c[0][1][0][j]+(a.c[2][1][0][j]-a.c[2][2][0][j])*-g.c[0][2][1][j]+(a.c[2][1][1][j]-a.c[2][2][1][j])*g.c[0][2][0][j]); \
      o.c[0][3][1][j] += + ((a.c[0][0][0][j]+a.c[0][3][0][j])*-g.c[0][0][1][j]+(a.c[0][0][1][j]+a.c[0][3][1][j])*g.c[0][0][0][j]+(a.c[1][0][0][j]+a.c[1][3][0][j])*-g.c[0][1][1][j]+(a.c[1][0][1][j]+a.c[1][3][1][j])*g.c[0][1][0][j]+(a.c[2][0][0][j]+a.c[2][3][0][j])*-g.c[0][2][1][j]+(a.c[2][0][1][j]+a.c[2][3][1][j])*g.c[0][2][0][j]); \
      o.c[1][0][0][j] += + ((a.c[0][0][0][j]+a.c[0][3][0][j])*g.c[1][0][0][j]-(a.c[0][0][1][j]+a.c[0][3][1][j])*-g.c[1][0][1][j]+(a.c[1][0][0][j]+a.c[1][3][0][j])*g.c[1][1][0][j]-(a.c[1][0][1][j]+a.c[1][3][1][j])*-g.c[1][1][1][j]+(a.c[2][0][0][j]+a.c[2][3][0][j])*g.c[1][2][0][j]-(a.c[2][0][1][j]+a.c[2][3][1][j])*-g.c[1][2][1][j]); \
      o.c[1][1][0][j] += + ((a.c[0][1][0][j]-a.c[0][2][0][j])*g.c[1][0][0][j]-(a.c[0][1][1][j]-a.c[0][2][1][j])*-g.c[1][0][1][j]+(a.c[1][1][0][j]-a.c[1][2][0][j])*g.c[1][1][0][j]-(a.c[1][1][1][j]-a.c[1][2][1][j])*-g.c[1][1][1][j]+(a.c[2][1][0][j]-a.c[2][2][0][j])*g.c[1][2][0][j]-(a.c[2][1][1][j]-a.c[2][2][1][j])*-g.c[1][2][1][j]); \
      o.c[1][2][0][j] += - ((a.c[0][1][0][j]-a.c[0][2][0][j])*g.c[1][0][0][j]-(a.c[0][1][1][j]-a.c[0][2][1][j])*-g.c[1][0][1][j]+(a.c[1][1][0][j]-a.c[1][2][0][j])*g.c[1][1][0][j]-(a.c[1][1][1][j]-a.c[1][2][1][j])*-g.c[1][1][1][j]+(a.c[2][1][0][j]-a.c[2][2][0][j])*g.c[1][2][0][j]-(a.c[2][1][1][j]-a.c[2][2][1][j])*-g.c[1][2][1][j]); \
      o.c[1][3][0][j] += + ((a.c[0][0][0][j]+a.c[0][3][0][j])*g.c[1][0][0][j]-(a.c[0][0][1][j]+a.c[0][3][1][j])*-g.c[1][0][1][j]+(a.c[1][0][0][j]+a.c[1][3][0][j])*g.c[1][1][0][j]-(a.c[1][0][1][j]+a.c[1][3][1][j])*-g.c[1][1][1][j]+(a.c[2][0][0][j]+a.c[2][3][0][j])*g.c[1][2][0][j]-(a.c[2][0][1][j]+a.c[2][3][1][j])*-g.c[1][2][1][j]); \
      o.c[1][0][1][j] += + ((a.c[0][0][0][j]+a.c[0][3][0][j])*-g.c[1][0][1][j]+(a.c[0][0][1][j]+a.c[0][3][1][j])*g.c[1][0][0][j]+(a.c[1][0][0][j]+a.c[1][3][0][j])*-g.c[1][1][1][j]+(a.c[1][0][1][j]+a.c[1][3][1][j])*g.c[1][1][0][j]+(a.c[2][0][0][j]+a.c[2][3][0][j])*-g.c[1][2][1][j]+(a.c[2][0][1][j]+a.c[2][3][1][j])*g.c[1][2][0][j]); \
      o.c[1][1][1][j] += + ((a.c[0][1][0][j]-a.c[0][2][0][j])*-g.c[1][0][1][j]+(a.c[0][1][1][j]-a.c[0][2][1][j])*g.c[1][0][0][j]+(a.c[1][1][0][j]-a.c[1][2][0][j])*-g.c[1][1][1][j]+(a.c[1][1][1][j]-a.c[1][2][1][j])*g.c[1][1][0][j]+(a.c[2][1][0][j]-a.c[2][2][0][j])*-g.c[1][2][1][j]+(a.c[2][1][1][j]-a.c[2][2][1][j])*g.c[1][2][0][j]); \
      o.c[1][2][1][j] += - ((a.c[0][1][0][j]-a.c[0][2][0][j])*-g.c[1][0][1][j]+(a.c[0][1][1][j]-a.c[0][2][1][j])*g.c[1][0][0][j]+(a.c[1][1][0][j]-a.c[1][2][0][j])*-g.c[1][1][1][j]+(a.c[1][1][1][j]-a.c[1][2][1][j])*g.c[1][1][0][j]+(a.c[2][1][0][j]-a.c[2][2][0][j])*-g.c[1][2][1][j]+(a.c[2][1][1][j]-a.c[2][2][1][j])*g.c[1][2][0][j]); \
      o.c[1][3][1][j] += + ((a.c[0][0][0][j]+a.c[0][3][0][j])*-g.c[1][0][1][j]+(a.c[0][0][1][j]+a.c[0][3][1][j])*g.c[1][0][0][j]+(a.c[1][0][0][j]+a.c[1][3][0][j])*-g.c[1][1][1][j]+(a.c[1][0][1][j]+a.c[1][3][1][j])*g.c[1][1][0][j]+(a.c[2][0][0][j]+a.c[2][3][0][j])*-g.c[1][2][1][j]+(a.c[2][0][1][j]+a.c[2][3][1][j])*g.c[1][2][0][j]); \
      o.c[2][0][0][j] += + ((a.c[0][0][0][j]+a.c[0][3][0][j])*g.c[2][0][0][j]-(a.c[0][0][1][j]+a.c[0][3][1][j])*-g.c[2][0][1][j]+(a.c[1][0][0][j]+a.c[1][3][0][j])*g.c[2][1][0][j]-(a.c[1][0][1][j]+a.c[1][3][1][j])*-g.c[2][1][1][j]+(a.c[2][0][0][j]+a.c[2][3][0][j])*g.c[2][2][0][j]-(a.c[2][0][1][j]+a.c[2][3][1][j])*-g.c[2][2][1][j]); \
      o.c[2][1][0][j] += + ((a.c[0][1][0][j]-a.c[0][2][0][j])*g.c[2][0][0][j]-(a.c[0][1][1][j]-a.c[0][2][1][j])*-g.c[2][0][1][j]+(a.c[1][1][0][j]-a.c[1][2][0][j])*g.c[2][1][0][j]-(a.c[1][1][1][j]-a.c[1][2][1][j])*-g.c[2][1][1][j]+(a.c[2][1][0][j]-a.c[2][2][0][j])*g.c[2][2][0][j]-(a.c[2][1][1][j]-a.c[2][2][1][j])*-g.c[2][2][1][j]); \
      o.c[2][2][0][j] += - ((a.c[0][1][0][j]-a.c[0][2][0][j])*g.c[2][0][0][j]-(a.c[0][1][1][j]-a.c[0][2][1][j])*-g.c[2][0][1][j]+(a.c[1][1][0][j]-a.c[1][2][0][j])*g.c[2][1][0][j]-(a.c[1][1][1][j]-a.c[1][2][1][j])*-g.c[2][1][1][j]+(a.c[2][1][0][j]-a.c[2][2][0][j])*g.c[2][2][0][j]-(a.c[2][1][1][j]-a.c[2][2][1][j])*-g.c[2][2][1][j]); \
      o.c[2][3][0][j] += + ((a.c[0][0][0][j]+a.c[0][3][0][j])*g.c[2][0][0][j]-(a.c[0][0][1][j]+a.c[0][3][1][j])*-g.c[2][0][1][j]+(a.c[1][0][0][j]+a.c[1][3][0][j])*g.c[2][1][0][j]-(a.c[1][0][1][j]+a.c[1][3][1][j])*-g.c[2][1][1][j]+(a.c[2][0][0][j]+a.c[2][3][0][j])*g.c[2][2][0][j]-(a.c[2][0][1][j]+a.c[2][3][1][j])*-g.c[2][2][1][j]); \
      o.c[2][0][1][j] += + ((a.c[0][0][0][j]+a.c[0][3][0][j])*-g.c[2][0][1][j]+(a.c[0][0][1][j]+a.c[0][3][1][j])*g.c[2][0][0][j]+(a.c[1][0][0][j]+a.c[1][3][0][j])*-g.c[2][1][1][j]+(a.c[1][0][1][j]+a.c[1][3][1][j])*g.c[2][1][0][j]+(a.c[2][0][0][j]+a.c[2][3][0][j])*-g.c[2][2][1][j]+(a.c[2][0][1][j]+a.c[2][3][1][j])*g.c[2][2][0][j]); \
      o.c[2][1][1][j] += + ((a.c[0][1][0][j]-a.c[0][2][0][j])*-g.c[2][0][1][j]+(a.c[0][1][1][j]-a.c[0][2][1][j])*g.c[2][0][0][j]+(a.c[1][1][0][j]-a.c[1][2][0][j])*-g.c[2][1][1][j]+(a.c[1][1][1][j]-a.c[1][2][1][j])*g.c[2][1][0][j]+(a.c[2][1][0][j]-a.c[2][2][0][j])*-g.c[2][2][1][j]+(a.c[2][1][1][j]-a.c[2][2][1][j])*g.c[2][2][0][j]); \
      o.c[2][2][1][j] += - ((a.c[0][1][0][j]-a.c[0][2][0][j])*-g.c[2][0][1][j]+(a.c[0][1][1][j]-a.c[0][2][1][j])*g.c[2][0][0][j]+(a.c[1][1][0][j]-a.c[1][2][0][j])*-g.c[2][1][1][j]+(a.c[1][1][1][j]-a.c[1][2][1][j])*g.c[2][1][0][j]+(a.c[2][1][0][j]-a.c[2][2][0][j])*-g.c[2][2][1][j]+(a.c[2][1][1][j]-a.c[2][2][1][j])*g.c[2][2][0][j]); \
      o.c[2][3][1][j] += + ((a.c[0][0][0][j]+a.c[0][3][0][j])*-g.c[2][0][1][j]+(a.c[0][0][1][j]+a.c[0][3][1][j])*g.c[2][0][0][j]+(a.c[1][0][0][j]+a.c[1][3][0][j])*-g.c[2][1][1][j]+(a.c[1][0][1][j]+a.c[1][3][1][j])*g.c[2][1][0][j]+(a.c[2][0][0][j]+a.c[2][3][0][j])*-g.c[2][2][1][j]+(a.c[2][0][1][j]+a.c[2][3][1][j])*g.c[2][2][0][j]); \
    }									\
    }

#define __mult_zfwd_s_def(o,a,g){			\
    for(int j=0;j<VLENS;j++){						\
      o.c[0][0][0][j] += + ((a.c[0][0][0][j]+a.c[0][2][1][j])*g.c[0][0][0][j]-(a.c[0][0][1][j]-a.c[0][2][0][j])*g.c[0][0][1][j]+(a.c[1][0][0][j]+a.c[1][2][1][j])*g.c[1][0][0][j]-(a.c[1][0][1][j]-a.c[1][2][0][j])*g.c[1][0][1][j]+(a.c[2][0][0][j]+a.c[2][2][1][j])*g.c[2][0][0][j]-(a.c[2][0][1][j]-a.c[2][2][0][j])*g.c[2][0][1][j]); \
      o.c[0][1][0][j] += + ((a.c[0][1][0][j]-a.c[0][3][1][j])*g.c[0][0][0][j]-(a.c[0][1][1][j]+a.c[0][3][0][j])*g.c[0][0][1][j]+(a.c[1][1][0][j]-a.c[1][3][1][j])*g.c[1][0][0][j]-(a.c[1][1][1][j]+a.c[1][3][0][j])*g.c[1][0][1][j]+(a.c[2][1][0][j]-a.c[2][3][1][j])*g.c[2][0][0][j]-(a.c[2][1][1][j]+a.c[2][3][0][j])*g.c[2][0][1][j]); \
      o.c[0][2][0][j] += - ((a.c[0][0][0][j]+a.c[0][2][1][j])*g.c[0][0][1][j]+(a.c[0][0][1][j]-a.c[0][2][0][j])*g.c[0][0][0][j]+(a.c[1][0][0][j]+a.c[1][2][1][j])*g.c[1][0][1][j]+(a.c[1][0][1][j]-a.c[1][2][0][j])*g.c[1][0][0][j]+(a.c[2][0][0][j]+a.c[2][2][1][j])*g.c[2][0][1][j]+(a.c[2][0][1][j]-a.c[2][2][0][j])*g.c[2][0][0][j]); \
      o.c[0][3][0][j] += + ((a.c[0][1][0][j]-a.c[0][3][1][j])*g.c[0][0][1][j]+(a.c[0][1][1][j]+a.c[0][3][0][j])*g.c[0][0][0][j]+(a.c[1][1][0][j]-a.c[1][3][1][j])*g.c[1][0][1][j]+(a.c[1][1][1][j]+a.c[1][3][0][j])*g.c[1][0][0][j]+(a.c[2][1][0][j]-a.c[2][3][1][j])*g.c[2][0][1][j]+(a.c[2][1][1][j]+a.c[2][3][0][j])*g.c[2][0][0][j]); \
      o.c[0][0][1][j] += + ((a.c[0][0][0][j]+a.c[0][2][1][j])*g.c[0][0][1][j]+(a.c[0][0][1][j]-a.c[0][2][0][j])*g.c[0][0][0][j]+(a.c[1][0][0][j]+a.c[1][2][1][j])*g.c[1][0][1][j]+(a.c[1][0][1][j]-a.c[1][2][0][j])*g.c[1][0][0][j]+(a.c[2][0][0][j]+a.c[2][2][1][j])*g.c[2][0][1][j]+(a.c[2][0][1][j]-a.c[2][2][0][j])*g.c[2][0][0][j]); \
      o.c[0][1][1][j] += + ((a.c[0][1][0][j]-a.c[0][3][1][j])*g.c[0][0][1][j]+(a.c[0][1][1][j]+a.c[0][3][0][j])*g.c[0][0][0][j]+(a.c[1][1][0][j]-a.c[1][3][1][j])*g.c[1][0][1][j]+(a.c[1][1][1][j]+a.c[1][3][0][j])*g.c[1][0][0][j]+(a.c[2][1][0][j]-a.c[2][3][1][j])*g.c[2][0][1][j]+(a.c[2][1][1][j]+a.c[2][3][0][j])*g.c[2][0][0][j]); \
      o.c[0][2][1][j] += + ((a.c[0][0][0][j]+a.c[0][2][1][j])*g.c[0][0][0][j]-(a.c[0][0][1][j]-a.c[0][2][0][j])*g.c[0][0][1][j]+(a.c[1][0][0][j]+a.c[1][2][1][j])*g.c[1][0][0][j]-(a.c[1][0][1][j]-a.c[1][2][0][j])*g.c[1][0][1][j]+(a.c[2][0][0][j]+a.c[2][2][1][j])*g.c[2][0][0][j]-(a.c[2][0][1][j]-a.c[2][2][0][j])*g.c[2][0][1][j]); \
      o.c[0][3][1][j] += - ((a.c[0][1][0][j]-a.c[0][3][1][j])*g.c[0][0][0][j]-(a.c[0][1][1][j]+a.c[0][3][0][j])*g.c[0][0][1][j]+(a.c[1][1][0][j]-a.c[1][3][1][j])*g.c[1][0][0][j]-(a.c[1][1][1][j]+a.c[1][3][0][j])*g.c[1][0][1][j]+(a.c[2][1][0][j]-a.c[2][3][1][j])*g.c[2][0][0][j]-(a.c[2][1][1][j]+a.c[2][3][0][j])*g.c[2][0][1][j]); \
      o.c[1][0][0][j] += + ((a.c[0][0][0][j]+a.c[0][2][1][j])*g.c[0][1][0][j]-(a.c[0][0][1][j]-a.c[0][2][0][j])*g.c[0][1][1][j]+(a.c[1][0][0][j]+a.c[1][2][1][j])*g.c[1][1][0][j]-(a.c[1][0][1][j]-a.c[1][2][0][j])*g.c[1][1][1][j]+(a.c[2][0][0][j]+a.c[2][2][1][j])*g.c[2][1][0][j]-(a.c[2][0][1][j]-a.c[2][2][0][j])*g.c[2][1][1][j]); \
      o.c[1][1][0][j] += + ((a.c[0][1][0][j]-a.c[0][3][1][j])*g.c[0][1][0][j]-(a.c[0][1][1][j]+a.c[0][3][0][j])*g.c[0][1][1][j]+(a.c[1][1][0][j]-a.c[1][3][1][j])*g.c[1][1][0][j]-(a.c[1][1][1][j]+a.c[1][3][0][j])*g.c[1][1][1][j]+(a.c[2][1][0][j]-a.c[2][3][1][j])*g.c[2][1][0][j]-(a.c[2][1][1][j]+a.c[2][3][0][j])*g.c[2][1][1][j]); \
      o.c[1][2][0][j] += - ((a.c[0][0][0][j]+a.c[0][2][1][j])*g.c[0][1][1][j]+(a.c[0][0][1][j]-a.c[0][2][0][j])*g.c[0][1][0][j]+(a.c[1][0][0][j]+a.c[1][2][1][j])*g.c[1][1][1][j]+(a.c[1][0][1][j]-a.c[1][2][0][j])*g.c[1][1][0][j]+(a.c[2][0][0][j]+a.c[2][2][1][j])*g.c[2][1][1][j]+(a.c[2][0][1][j]-a.c[2][2][0][j])*g.c[2][1][0][j]); \
      o.c[1][3][0][j] += + ((a.c[0][1][0][j]-a.c[0][3][1][j])*g.c[0][1][1][j]+(a.c[0][1][1][j]+a.c[0][3][0][j])*g.c[0][1][0][j]+(a.c[1][1][0][j]-a.c[1][3][1][j])*g.c[1][1][1][j]+(a.c[1][1][1][j]+a.c[1][3][0][j])*g.c[1][1][0][j]+(a.c[2][1][0][j]-a.c[2][3][1][j])*g.c[2][1][1][j]+(a.c[2][1][1][j]+a.c[2][3][0][j])*g.c[2][1][0][j]); \
      o.c[1][0][1][j] += + ((a.c[0][0][0][j]+a.c[0][2][1][j])*g.c[0][1][1][j]+(a.c[0][0][1][j]-a.c[0][2][0][j])*g.c[0][1][0][j]+(a.c[1][0][0][j]+a.c[1][2][1][j])*g.c[1][1][1][j]+(a.c[1][0][1][j]-a.c[1][2][0][j])*g.c[1][1][0][j]+(a.c[2][0][0][j]+a.c[2][2][1][j])*g.c[2][1][1][j]+(a.c[2][0][1][j]-a.c[2][2][0][j])*g.c[2][1][0][j]); \
      o.c[1][1][1][j] += + ((a.c[0][1][0][j]-a.c[0][3][1][j])*g.c[0][1][1][j]+(a.c[0][1][1][j]+a.c[0][3][0][j])*g.c[0][1][0][j]+(a.c[1][1][0][j]-a.c[1][3][1][j])*g.c[1][1][1][j]+(a.c[1][1][1][j]+a.c[1][3][0][j])*g.c[1][1][0][j]+(a.c[2][1][0][j]-a.c[2][3][1][j])*g.c[2][1][1][j]+(a.c[2][1][1][j]+a.c[2][3][0][j])*g.c[2][1][0][j]); \
      o.c[1][2][1][j] += + ((a.c[0][0][0][j]+a.c[0][2][1][j])*g.c[0][1][0][j]-(a.c[0][0][1][j]-a.c[0][2][0][j])*g.c[0][1][1][j]+(a.c[1][0][0][j]+a.c[1][2][1][j])*g.c[1][1][0][j]-(a.c[1][0][1][j]-a.c[1][2][0][j])*g.c[1][1][1][j]+(a.c[2][0][0][j]+a.c[2][2][1][j])*g.c[2][1][0][j]-(a.c[2][0][1][j]-a.c[2][2][0][j])*g.c[2][1][1][j]); \
      o.c[1][3][1][j] += - ((a.c[0][1][0][j]-a.c[0][3][1][j])*g.c[0][1][0][j]-(a.c[0][1][1][j]+a.c[0][3][0][j])*g.c[0][1][1][j]+(a.c[1][1][0][j]-a.c[1][3][1][j])*g.c[1][1][0][j]-(a.c[1][1][1][j]+a.c[1][3][0][j])*g.c[1][1][1][j]+(a.c[2][1][0][j]-a.c[2][3][1][j])*g.c[2][1][0][j]-(a.c[2][1][1][j]+a.c[2][3][0][j])*g.c[2][1][1][j]); \
      o.c[2][0][0][j] += + ((a.c[0][0][0][j]+a.c[0][2][1][j])*g.c[0][2][0][j]-(a.c[0][0][1][j]-a.c[0][2][0][j])*g.c[0][2][1][j]+(a.c[1][0][0][j]+a.c[1][2][1][j])*g.c[1][2][0][j]-(a.c[1][0][1][j]-a.c[1][2][0][j])*g.c[1][2][1][j]+(a.c[2][0][0][j]+a.c[2][2][1][j])*g.c[2][2][0][j]-(a.c[2][0][1][j]-a.c[2][2][0][j])*g.c[2][2][1][j]); \
      o.c[2][1][0][j] += + ((a.c[0][1][0][j]-a.c[0][3][1][j])*g.c[0][2][0][j]-(a.c[0][1][1][j]+a.c[0][3][0][j])*g.c[0][2][1][j]+(a.c[1][1][0][j]-a.c[1][3][1][j])*g.c[1][2][0][j]-(a.c[1][1][1][j]+a.c[1][3][0][j])*g.c[1][2][1][j]+(a.c[2][1][0][j]-a.c[2][3][1][j])*g.c[2][2][0][j]-(a.c[2][1][1][j]+a.c[2][3][0][j])*g.c[2][2][1][j]); \
      o.c[2][2][0][j] += - ((a.c[0][0][0][j]+a.c[0][2][1][j])*g.c[0][2][1][j]+(a.c[0][0][1][j]-a.c[0][2][0][j])*g.c[0][2][0][j]+(a.c[1][0][0][j]+a.c[1][2][1][j])*g.c[1][2][1][j]+(a.c[1][0][1][j]-a.c[1][2][0][j])*g.c[1][2][0][j]+(a.c[2][0][0][j]+a.c[2][2][1][j])*g.c[2][2][1][j]+(a.c[2][0][1][j]-a.c[2][2][0][j])*g.c[2][2][0][j]); \
      o.c[2][3][0][j] += + ((a.c[0][1][0][j]-a.c[0][3][1][j])*g.c[0][2][1][j]+(a.c[0][1][1][j]+a.c[0][3][0][j])*g.c[0][2][0][j]+(a.c[1][1][0][j]-a.c[1][3][1][j])*g.c[1][2][1][j]+(a.c[1][1][1][j]+a.c[1][3][0][j])*g.c[1][2][0][j]+(a.c[2][1][0][j]-a.c[2][3][1][j])*g.c[2][2][1][j]+(a.c[2][1][1][j]+a.c[2][3][0][j])*g.c[2][2][0][j]); \
      o.c[2][0][1][j] += + ((a.c[0][0][0][j]+a.c[0][2][1][j])*g.c[0][2][1][j]+(a.c[0][0][1][j]-a.c[0][2][0][j])*g.c[0][2][0][j]+(a.c[1][0][0][j]+a.c[1][2][1][j])*g.c[1][2][1][j]+(a.c[1][0][1][j]-a.c[1][2][0][j])*g.c[1][2][0][j]+(a.c[2][0][0][j]+a.c[2][2][1][j])*g.c[2][2][1][j]+(a.c[2][0][1][j]-a.c[2][2][0][j])*g.c[2][2][0][j]); \
      o.c[2][1][1][j] += + ((a.c[0][1][0][j]-a.c[0][3][1][j])*g.c[0][2][1][j]+(a.c[0][1][1][j]+a.c[0][3][0][j])*g.c[0][2][0][j]+(a.c[1][1][0][j]-a.c[1][3][1][j])*g.c[1][2][1][j]+(a.c[1][1][1][j]+a.c[1][3][0][j])*g.c[1][2][0][j]+(a.c[2][1][0][j]-a.c[2][3][1][j])*g.c[2][2][1][j]+(a.c[2][1][1][j]+a.c[2][3][0][j])*g.c[2][2][0][j]); \
      o.c[2][2][1][j] += + ((a.c[0][0][0][j]+a.c[0][2][1][j])*g.c[0][2][0][j]-(a.c[0][0][1][j]-a.c[0][2][0][j])*g.c[0][2][1][j]+(a.c[1][0][0][j]+a.c[1][2][1][j])*g.c[1][2][0][j]-(a.c[1][0][1][j]-a.c[1][2][0][j])*g.c[1][2][1][j]+(a.c[2][0][0][j]+a.c[2][2][1][j])*g.c[2][2][0][j]-(a.c[2][0][1][j]-a.c[2][2][0][j])*g.c[2][2][1][j]); \
      o.c[2][3][1][j] += - ((a.c[0][1][0][j]-a.c[0][3][1][j])*g.c[0][2][0][j]-(a.c[0][1][1][j]+a.c[0][3][0][j])*g.c[0][2][1][j]+(a.c[1][1][0][j]-a.c[1][3][1][j])*g.c[1][2][0][j]-(a.c[1][1][1][j]+a.c[1][3][0][j])*g.c[1][2][1][j]+(a.c[2][1][0][j]-a.c[2][3][1][j])*g.c[2][2][0][j]-(a.c[2][1][1][j]+a.c[2][3][0][j])*g.c[2][2][1][j]); \
    }									\
    }

#define __mult_zbwd_s_def(o,a,g){			\
  for(int j=0;j<VLENS;j++){						\
    o.c[0][0][0][j] += + ((a.c[0][0][0][j]-a.c[0][2][1][j])*g.c[0][0][0][j]-(a.c[0][0][1][j]+a.c[0][2][0][j])*-g.c[0][0][1][j]+(a.c[1][0][0][j]-a.c[1][2][1][j])*g.c[0][1][0][j]-(a.c[1][0][1][j]+a.c[1][2][0][j])*-g.c[0][1][1][j]+(a.c[2][0][0][j]-a.c[2][2][1][j])*g.c[0][2][0][j]-(a.c[2][0][1][j]+a.c[2][2][0][j])*-g.c[0][2][1][j]); \
    o.c[0][1][0][j] += + ((a.c[0][1][0][j]+a.c[0][3][1][j])*g.c[0][0][0][j]-(a.c[0][1][1][j]-a.c[0][3][0][j])*-g.c[0][0][1][j]+(a.c[1][1][0][j]+a.c[1][3][1][j])*g.c[0][1][0][j]-(a.c[1][1][1][j]-a.c[1][3][0][j])*-g.c[0][1][1][j]+(a.c[2][1][0][j]+a.c[2][3][1][j])*g.c[0][2][0][j]-(a.c[2][1][1][j]-a.c[2][3][0][j])*-g.c[0][2][1][j]); \
    o.c[0][2][0][j] += + ((a.c[0][0][0][j]-a.c[0][2][1][j])*-g.c[0][0][1][j]+(a.c[0][0][1][j]+a.c[0][2][0][j])*g.c[0][0][0][j]+(a.c[1][0][0][j]-a.c[1][2][1][j])*-g.c[0][1][1][j]+(a.c[1][0][1][j]+a.c[1][2][0][j])*g.c[0][1][0][j]+(a.c[2][0][0][j]-a.c[2][2][1][j])*-g.c[0][2][1][j]+(a.c[2][0][1][j]+a.c[2][2][0][j])*g.c[0][2][0][j]); \
    o.c[0][3][0][j] += - ((a.c[0][1][0][j]+a.c[0][3][1][j])*-g.c[0][0][1][j]+(a.c[0][1][1][j]-a.c[0][3][0][j])*g.c[0][0][0][j]+(a.c[1][1][0][j]+a.c[1][3][1][j])*-g.c[0][1][1][j]+(a.c[1][1][1][j]-a.c[1][3][0][j])*g.c[0][1][0][j]+(a.c[2][1][0][j]+a.c[2][3][1][j])*-g.c[0][2][1][j]+(a.c[2][1][1][j]-a.c[2][3][0][j])*g.c[0][2][0][j]); \
    o.c[0][0][1][j] += + ((a.c[0][0][0][j]-a.c[0][2][1][j])*-g.c[0][0][1][j]+(a.c[0][0][1][j]+a.c[0][2][0][j])*g.c[0][0][0][j]+(a.c[1][0][0][j]-a.c[1][2][1][j])*-g.c[0][1][1][j]+(a.c[1][0][1][j]+a.c[1][2][0][j])*g.c[0][1][0][j]+(a.c[2][0][0][j]-a.c[2][2][1][j])*-g.c[0][2][1][j]+(a.c[2][0][1][j]+a.c[2][2][0][j])*g.c[0][2][0][j]); \
    o.c[0][1][1][j] += + ((a.c[0][1][0][j]+a.c[0][3][1][j])*-g.c[0][0][1][j]+(a.c[0][1][1][j]-a.c[0][3][0][j])*g.c[0][0][0][j]+(a.c[1][1][0][j]+a.c[1][3][1][j])*-g.c[0][1][1][j]+(a.c[1][1][1][j]-a.c[1][3][0][j])*g.c[0][1][0][j]+(a.c[2][1][0][j]+a.c[2][3][1][j])*-g.c[0][2][1][j]+(a.c[2][1][1][j]-a.c[2][3][0][j])*g.c[0][2][0][j]); \
    o.c[0][2][1][j] += - ((a.c[0][0][0][j]-a.c[0][2][1][j])*g.c[0][0][0][j]-(a.c[0][0][1][j]+a.c[0][2][0][j])*-g.c[0][0][1][j]+(a.c[1][0][0][j]-a.c[1][2][1][j])*g.c[0][1][0][j]-(a.c[1][0][1][j]+a.c[1][2][0][j])*-g.c[0][1][1][j]+(a.c[2][0][0][j]-a.c[2][2][1][j])*g.c[0][2][0][j]-(a.c[2][0][1][j]+a.c[2][2][0][j])*-g.c[0][2][1][j]); \
    o.c[0][3][1][j] += + ((a.c[0][1][0][j]+a.c[0][3][1][j])*g.c[0][0][0][j]-(a.c[0][1][1][j]-a.c[0][3][0][j])*-g.c[0][0][1][j]+(a.c[1][1][0][j]+a.c[1][3][1][j])*g.c[0][1][0][j]-(a.c[1][1][1][j]-a.c[1][3][0][j])*-g.c[0][1][1][j]+(a.c[2][1][0][j]+a.c[2][3][1][j])*g.c[0][2][0][j]-(a.c[2][1][1][j]-a.c[2][3][0][j])*-g.c[0][2][1][j]); \
    o.c[1][0][0][j] += + ((a.c[0][0][0][j]-a.c[0][2][1][j])*g.c[1][0][0][j]-(a.c[0][0][1][j]+a.c[0][2][0][j])*-g.c[1][0][1][j]+(a.c[1][0][0][j]-a.c[1][2][1][j])*g.c[1][1][0][j]-(a.c[1][0][1][j]+a.c[1][2][0][j])*-g.c[1][1][1][j]+(a.c[2][0][0][j]-a.c[2][2][1][j])*g.c[1][2][0][j]-(a.c[2][0][1][j]+a.c[2][2][0][j])*-g.c[1][2][1][j]); \
    o.c[1][1][0][j] += + ((a.c[0][1][0][j]+a.c[0][3][1][j])*g.c[1][0][0][j]-(a.c[0][1][1][j]-a.c[0][3][0][j])*-g.c[1][0][1][j]+(a.c[1][1][0][j]+a.c[1][3][1][j])*g.c[1][1][0][j]-(a.c[1][1][1][j]-a.c[1][3][0][j])*-g.c[1][1][1][j]+(a.c[2][1][0][j]+a.c[2][3][1][j])*g.c[1][2][0][j]-(a.c[2][1][1][j]-a.c[2][3][0][j])*-g.c[1][2][1][j]); \
    o.c[1][2][0][j] += + ((a.c[0][0][0][j]-a.c[0][2][1][j])*-g.c[1][0][1][j]+(a.c[0][0][1][j]+a.c[0][2][0][j])*g.c[1][0][0][j]+(a.c[1][0][0][j]-a.c[1][2][1][j])*-g.c[1][1][1][j]+(a.c[1][0][1][j]+a.c[1][2][0][j])*g.c[1][1][0][j]+(a.c[2][0][0][j]-a.c[2][2][1][j])*-g.c[1][2][1][j]+(a.c[2][0][1][j]+a.c[2][2][0][j])*g.c[1][2][0][j]); \
    o.c[1][3][0][j] += - ((a.c[0][1][0][j]+a.c[0][3][1][j])*-g.c[1][0][1][j]+(a.c[0][1][1][j]-a.c[0][3][0][j])*g.c[1][0][0][j]+(a.c[1][1][0][j]+a.c[1][3][1][j])*-g.c[1][1][1][j]+(a.c[1][1][1][j]-a.c[1][3][0][j])*g.c[1][1][0][j]+(a.c[2][1][0][j]+a.c[2][3][1][j])*-g.c[1][2][1][j]+(a.c[2][1][1][j]-a.c[2][3][0][j])*g.c[1][2][0][j]); \
    o.c[1][0][1][j] += + ((a.c[0][0][0][j]-a.c[0][2][1][j])*-g.c[1][0][1][j]+(a.c[0][0][1][j]+a.c[0][2][0][j])*g.c[1][0][0][j]+(a.c[1][0][0][j]-a.c[1][2][1][j])*-g.c[1][1][1][j]+(a.c[1][0][1][j]+a.c[1][2][0][j])*g.c[1][1][0][j]+(a.c[2][0][0][j]-a.c[2][2][1][j])*-g.c[1][2][1][j]+(a.c[2][0][1][j]+a.c[2][2][0][j])*g.c[1][2][0][j]); \
    o.c[1][1][1][j] += + ((a.c[0][1][0][j]+a.c[0][3][1][j])*-g.c[1][0][1][j]+(a.c[0][1][1][j]-a.c[0][3][0][j])*g.c[1][0][0][j]+(a.c[1][1][0][j]+a.c[1][3][1][j])*-g.c[1][1][1][j]+(a.c[1][1][1][j]-a.c[1][3][0][j])*g.c[1][1][0][j]+(a.c[2][1][0][j]+a.c[2][3][1][j])*-g.c[1][2][1][j]+(a.c[2][1][1][j]-a.c[2][3][0][j])*g.c[1][2][0][j]); \
    o.c[1][2][1][j] += - ((a.c[0][0][0][j]-a.c[0][2][1][j])*g.c[1][0][0][j]-(a.c[0][0][1][j]+a.c[0][2][0][j])*-g.c[1][0][1][j]+(a.c[1][0][0][j]-a.c[1][2][1][j])*g.c[1][1][0][j]-(a.c[1][0][1][j]+a.c[1][2][0][j])*-g.c[1][1][1][j]+(a.c[2][0][0][j]-a.c[2][2][1][j])*g.c[1][2][0][j]-(a.c[2][0][1][j]+a.c[2][2][0][j])*-g.c[1][2][1][j]); \
    o.c[1][3][1][j] += + ((a.c[0][1][0][j]+a.c[0][3][1][j])*g.c[1][0][0][j]-(a.c[0][1][1][j]-a.c[0][3][0][j])*-g.c[1][0][1][j]+(a.c[1][1][0][j]+a.c[1][3][1][j])*g.c[1][1][0][j]-(a.c[1][1][1][j]-a.c[1][3][0][j])*-g.c[1][1][1][j]+(a.c[2][1][0][j]+a.c[2][3][1][j])*g.c[1][2][0][j]-(a.c[2][1][1][j]-a.c[2][3][0][j])*-g.c[1][2][1][j]); \
    o.c[2][0][0][j] += + ((a.c[0][0][0][j]-a.c[0][2][1][j])*g.c[2][0][0][j]-(a.c[0][0][1][j]+a.c[0][2][0][j])*-g.c[2][0][1][j]+(a.c[1][0][0][j]-a.c[1][2][1][j])*g.c[2][1][0][j]-(a.c[1][0][1][j]+a.c[1][2][0][j])*-g.c[2][1][1][j]+(a.c[2][0][0][j]-a.c[2][2][1][j])*g.c[2][2][0][j]-(a.c[2][0][1][j]+a.c[2][2][0][j])*-g.c[2][2][1][j]); \
    o.c[2][1][0][j] += + ((a.c[0][1][0][j]+a.c[0][3][1][j])*g.c[2][0][0][j]-(a.c[0][1][1][j]-a.c[0][3][0][j])*-g.c[2][0][1][j]+(a.c[1][1][0][j]+a.c[1][3][1][j])*g.c[2][1][0][j]-(a.c[1][1][1][j]-a.c[1][3][0][j])*-g.c[2][1][1][j]+(a.c[2][1][0][j]+a.c[2][3][1][j])*g.c[2][2][0][j]-(a.c[2][1][1][j]-a.c[2][3][0][j])*-g.c[2][2][1][j]); \
    o.c[2][2][0][j] += + ((a.c[0][0][0][j]-a.c[0][2][1][j])*-g.c[2][0][1][j]+(a.c[0][0][1][j]+a.c[0][2][0][j])*g.c[2][0][0][j]+(a.c[1][0][0][j]-a.c[1][2][1][j])*-g.c[2][1][1][j]+(a.c[1][0][1][j]+a.c[1][2][0][j])*g.c[2][1][0][j]+(a.c[2][0][0][j]-a.c[2][2][1][j])*-g.c[2][2][1][j]+(a.c[2][0][1][j]+a.c[2][2][0][j])*g.c[2][2][0][j]); \
    o.c[2][3][0][j] += - ((a.c[0][1][0][j]+a.c[0][3][1][j])*-g.c[2][0][1][j]+(a.c[0][1][1][j]-a.c[0][3][0][j])*g.c[2][0][0][j]+(a.c[1][1][0][j]+a.c[1][3][1][j])*-g.c[2][1][1][j]+(a.c[1][1][1][j]-a.c[1][3][0][j])*g.c[2][1][0][j]+(a.c[2][1][0][j]+a.c[2][3][1][j])*-g.c[2][2][1][j]+(a.c[2][1][1][j]-a.c[2][3][0][j])*g.c[2][2][0][j]); \
    o.c[2][0][1][j] += + ((a.c[0][0][0][j]-a.c[0][2][1][j])*-g.c[2][0][1][j]+(a.c[0][0][1][j]+a.c[0][2][0][j])*g.c[2][0][0][j]+(a.c[1][0][0][j]-a.c[1][2][1][j])*-g.c[2][1][1][j]+(a.c[1][0][1][j]+a.c[1][2][0][j])*g.c[2][1][0][j]+(a.c[2][0][0][j]-a.c[2][2][1][j])*-g.c[2][2][1][j]+(a.c[2][0][1][j]+a.c[2][2][0][j])*g.c[2][2][0][j]); \
    o.c[2][1][1][j] += + ((a.c[0][1][0][j]+a.c[0][3][1][j])*-g.c[2][0][1][j]+(a.c[0][1][1][j]-a.c[0][3][0][j])*g.c[2][0][0][j]+(a.c[1][1][0][j]+a.c[1][3][1][j])*-g.c[2][1][1][j]+(a.c[1][1][1][j]-a.c[1][3][0][j])*g.c[2][1][0][j]+(a.c[2][1][0][j]+a.c[2][3][1][j])*-g.c[2][2][1][j]+(a.c[2][1][1][j]-a.c[2][3][0][j])*g.c[2][2][0][j]); \
    o.c[2][2][1][j] += - ((a.c[0][0][0][j]-a.c[0][2][1][j])*g.c[2][0][0][j]-(a.c[0][0][1][j]+a.c[0][2][0][j])*-g.c[2][0][1][j]+(a.c[1][0][0][j]-a.c[1][2][1][j])*g.c[2][1][0][j]-(a.c[1][0][1][j]+a.c[1][2][0][j])*-g.c[2][1][1][j]+(a.c[2][0][0][j]-a.c[2][2][1][j])*g.c[2][2][0][j]-(a.c[2][0][1][j]+a.c[2][2][0][j])*-g.c[2][2][1][j]); \
    o.c[2][3][1][j] += + ((a.c[0][1][0][j]+a.c[0][3][1][j])*g.c[2][0][0][j]-(a.c[0][1][1][j]-a.c[0][3][0][j])*-g.c[2][0][1][j]+(a.c[1][1][0][j]+a.c[1][3][1][j])*g.c[2][1][0][j]-(a.c[1][1][1][j]-a.c[1][3][0][j])*-g.c[2][1][1][j]+(a.c[2][1][0][j]+a.c[2][3][1][j])*g.c[2][2][0][j]-(a.c[2][1][1][j]-a.c[2][3][0][j])*-g.c[2][2][1][j]); \
  }									\
  }

#define __mult_tbwd_s_def(o,a,g){					\
  for(int j=0; j<VLENS; j++) {					\
    o.c[0][0][0][j] += 2 * (a.c[0][0][0][j]*g.c[0][0][0][j]-a.c[0][0][1][j]*-g.c[0][0][1][j]+a.c[1][0][0][j]*g.c[0][1][0][j]-a.c[1][0][1][j]*-g.c[0][1][1][j]+a.c[2][0][0][j]*g.c[0][2][0][j]-a.c[2][0][1][j]*-g.c[0][2][1][j]); \
    o.c[0][0][1][j] += 2 * (a.c[0][0][0][j]*-g.c[0][0][1][j]+a.c[0][0][1][j]*g.c[0][0][0][j]+a.c[1][0][0][j]*-g.c[0][1][1][j]+a.c[1][0][1][j]*g.c[0][1][0][j]+a.c[2][0][0][j]*-g.c[0][2][1][j]+a.c[2][0][1][j]*g.c[0][2][0][j]); \
    o.c[0][1][0][j] += 2 * (a.c[0][1][0][j]*g.c[0][0][0][j]-a.c[0][1][1][j]*-g.c[0][0][1][j]+a.c[1][1][0][j]*g.c[0][1][0][j]-a.c[1][1][1][j]*-g.c[0][1][1][j]+a.c[2][1][0][j]*g.c[0][2][0][j]-a.c[2][1][1][j]*-g.c[0][2][1][j]); \
    o.c[0][1][1][j] += 2 * (a.c[0][1][0][j]*-g.c[0][0][1][j]+a.c[0][1][1][j]*g.c[0][0][0][j]+a.c[1][1][0][j]*-g.c[0][1][1][j]+a.c[1][1][1][j]*g.c[0][1][0][j]+a.c[2][1][0][j]*-g.c[0][2][1][j]+a.c[2][1][1][j]*g.c[0][2][0][j]); \
    o.c[1][0][0][j] += 2 * (a.c[0][0][0][j]*g.c[1][0][0][j]-a.c[0][0][1][j]*-g.c[1][0][1][j]+a.c[1][0][0][j]*g.c[1][1][0][j]-a.c[1][0][1][j]*-g.c[1][1][1][j]+a.c[2][0][0][j]*g.c[1][2][0][j]-a.c[2][0][1][j]*-g.c[1][2][1][j]); \
    o.c[1][0][1][j] += 2 * (a.c[0][0][0][j]*-g.c[1][0][1][j]+a.c[0][0][1][j]*g.c[1][0][0][j]+a.c[1][0][0][j]*-g.c[1][1][1][j]+a.c[1][0][1][j]*g.c[1][1][0][j]+a.c[2][0][0][j]*-g.c[1][2][1][j]+a.c[2][0][1][j]*g.c[1][2][0][j]); \
    o.c[1][1][0][j] += 2 * (a.c[0][1][0][j]*g.c[1][0][0][j]-a.c[0][1][1][j]*-g.c[1][0][1][j]+a.c[1][1][0][j]*g.c[1][1][0][j]-a.c[1][1][1][j]*-g.c[1][1][1][j]+a.c[2][1][0][j]*g.c[1][2][0][j]-a.c[2][1][1][j]*-g.c[1][2][1][j]); \
    o.c[1][1][1][j] += 2 * (a.c[0][1][0][j]*-g.c[1][0][1][j]+a.c[0][1][1][j]*g.c[1][0][0][j]+a.c[1][1][0][j]*-g.c[1][1][1][j]+a.c[1][1][1][j]*g.c[1][1][0][j]+a.c[2][1][0][j]*-g.c[1][2][1][j]+a.c[2][1][1][j]*g.c[1][2][0][j]); \
    o.c[2][0][0][j] += 2 * (a.c[0][0][0][j]*g.c[2][0][0][j]-a.c[0][0][1][j]*-g.c[2][0][1][j]+a.c[1][0][0][j]*g.c[2][1][0][j]-a.c[1][0][1][j]*-g.c[2][1][1][j]+a.c[2][0][0][j]*g.c[2][2][0][j]-a.c[2][0][1][j]*-g.c[2][2][1][j]); \
    o.c[2][0][1][j] += 2 * (a.c[0][0][0][j]*-g.c[2][0][1][j]+a.c[0][0][1][j]*g.c[2][0][0][j]+a.c[1][0][0][j]*-g.c[2][1][1][j]+a.c[1][0][1][j]*g.c[2][1][0][j]+a.c[2][0][0][j]*-g.c[2][2][1][j]+a.c[2][0][1][j]*g.c[2][2][0][j]); \
    o.c[2][1][0][j] += 2 * (a.c[0][1][0][j]*g.c[2][0][0][j]-a.c[0][1][1][j]*-g.c[2][0][1][j]+a.c[1][1][0][j]*g.c[2][1][0][j]-a.c[1][1][1][j]*-g.c[2][1][1][j]+a.c[2][1][0][j]*g.c[2][2][0][j]-a.c[2][1][1][j]*-g.c[2][2][1][j]); \
    o.c[2][1][1][j] += 2 * (a.c[0][1][0][j]*-g.c[2][0][1][j]+a.c[0][1][1][j]*g.c[2][0][0][j]+a.c[1][1][0][j]*-g.c[2][1][1][j]+a.c[1][1][1][j]*g.c[2][1][0][j]+a.c[2][1][0][j]*-g.c[2][2][1][j]+a.c[2][1][1][j]*g.c[2][2][0][j]); \
  }									\
  }

#define __su3recon_s_def(g33,g23){		\
  for(int c2=0;c2<2;c2++){			\
    for(int c1=0;c1<3;c1++){					\
      for(int j=0;j<VLENS;j++){					\
	g33.c[c2][c1][0][j] = g23.c[c2][c1][0][j];		\
	g33.c[c2][c1][1][j] = g23.c[c2][c1][1][j];		\
      }								\
    }									\
  }									\
  for(int j=0;j<VLENS;j++){						\
    g33.c[2][0][0][j] = g33.c[0][1][0][j]*g33.c[1][2][0][j]-g33.c[0][1][1][j]*g33.c[1][2][1][j]-g33.c[0][2][0][j]*g33.c[1][1][0][j]+g33.c[0][2][1][j]*g33.c[1][1][1][j]; \
    g33.c[2][1][0][j] = g33.c[0][2][0][j]*g33.c[1][0][0][j]-g33.c[0][2][1][j]*g33.c[1][0][1][j]-g33.c[0][0][0][j]*g33.c[1][2][0][j]+g33.c[0][0][1][j]*g33.c[1][2][1][j]; \
    g33.c[2][2][0][j] = g33.c[0][0][0][j]*g33.c[1][1][0][j]-g33.c[0][0][1][j]*g33.c[1][1][1][j]-g33.c[0][1][0][j]*g33.c[1][0][0][j]+g33.c[0][1][1][j]*g33.c[1][0][1][j]; \
    g33.c[2][0][1][j] =-g33.c[0][1][0][j]*g33.c[1][2][1][j]-g33.c[0][1][1][j]*g33.c[1][2][0][j]+g33.c[0][2][0][j]*g33.c[1][1][1][j]+g33.c[0][2][1][j]*g33.c[1][1][0][j]; \
    g33.c[2][1][1][j] =-g33.c[0][2][0][j]*g33.c[1][0][1][j]-g33.c[0][2][1][j]*g33.c[1][0][0][j]+g33.c[0][0][0][j]*g33.c[1][2][1][j]+g33.c[0][0][1][j]*g33.c[1][2][0][j]; \
    g33.c[2][2][1][j] =-g33.c[0][0][0][j]*g33.c[1][1][1][j]-g33.c[0][0][1][j]*g33.c[1][1][0][j]+g33.c[0][1][0][j]*g33.c[1][0][1][j]+g33.c[0][1][1][j]*g33.c[1][0][0][j]; \
  }									\
  }


#endif
