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
#include "wilson_d.h"
#include "clover_d.h"


#ifdef __cplusplus
extern "C"{
#endif


  extern __attribute__((aligned(64))) pglud_t glud;
  extern __attribute__((aligned(64))) pclvd_t clvd;
  extern double kappa, mkappa;
  extern int nxd, ny, nz, nt;
  extern int vold;// = nxd * ny * nz * nt;
  //---------------------------------------------------------------------------------------- mult D in bulk
  void ddd_in_d_(scd_t* out, scd_t* in, int* DEO){
    __attribute__((aligned(64))) scd_t tmp;
    int x, y, z, t;
    __attribute__((aligned(64))) g33d_t g;
    __attribute__((aligned(64))) projscd_t a, b;
    __attribute__((aligned(64))) double v[4][2][VLEND];

    pglud_t gx = &glud[vold*0 + NDIM*vold*(*DEO)];
    pglud_t gy = &glud[vold*1 + NDIM*vold*(*DEO)];
    pglud_t gz = &glud[vold*2 + NDIM*vold*(*DEO)];
    pglud_t gt = &glud[vold*3 + NDIM*vold*(*DEO)];

#pragma omp parallel for private(t,z,y,x,tmp, g,a,b,v) collapse(3) schedule(static)
    for(t=0; t<nt; t++){
      for(z=0; z<nz; z++){
	for(y=0; y<ny; y++){
	  //int ix = y + ny*z + ny*nz*t;
	  for(x=0; x<nxd; x++){
	    __zero_sc(tmp.c);
	    int i0  =      x        + nxd*y + nxd*ny*z + nxd*ny*nz*t;
	    int ixf = (x+1) + nxd*y + nxd*ny*z + nxd*ny*nz*t;
	    int ixb = (x-1) + nxd*y + nxd*ny*z + nxd*ny*nz*t;
	    int iyf = x + (y+1)*nxd + nxd*ny*z + nxd*ny*nz*t;
	    int iyb = x + (y-1)*nxd + nxd*ny*z + nxd*ny*nz*t;
	    int izf = x + nxd*y + (z+1)*nxd*ny + nxd*ny*nz*t;
	    int izb = x + nxd*y + (z-1)*nxd*ny + nxd*ny*nz*t;
	    int itf = x + nxd*y + nxd*ny*z + (t+1)*nxd*ny*nz;
	    int itb = x + nxd*y + nxd*ny*z + (t-1)*nxd*ny*nz;

	    if(x != nxd-1 ) {
	      __mult_xfwd1(tmp.c, gx, in, i0, ixf);
	    } else {
	      for(int c=0;c<3;c++){
	    	for (int s=0;s<4;s++) {
	    	  for (int ri=0;ri<2;ri++) {
	    	    for(int j=0;j<VLEND-1;j++){
	    	      v[s][ri][j]=in[i0].c[c][s][ri][j+1]; 
	    	    }
	    	    v[s][ri][VLEND-1]=0; 
	    	  }
	    	}
	    	for(int j=0;j<VLEND;j++){
	    	  a.c[c][0][0][j] = v[0][0][j] + v[3][1][j];
	    	  a.c[c][1][0][j] = v[1][0][j] + v[2][1][j];
	    	  a.c[c][0][1][j] = v[0][1][j] - v[3][0][j];
	    	  a.c[c][1][1][j] = v[1][1][j] - v[2][0][j];
	    	}
	      }
	      __load_su3(g.c, gx, i0);
	      __m_fg(b.c, a.c, g.c);
	      __s_xf(tmp.c, b.c);
	    }
	    if(x != 0 ) {
	      __mult_xbwd0(tmp.c, gx, in, i0, ixb);
	    } else {
#if SU3_RECONSTRUCT_D == 18
	      for(int c2=0;c2<3;c2++){
		for(int c1=0;c1<3;c1++){
		  g.c[c1][c2][0][0] =0;
		  g.c[c1][c2][1][0] =0;
		  for(int j=1;j<VLEND;j++){
		    g.c[c1][c2][0][j] = (*(gx+i0)).c[c2][c1][0][j-1];
		    g.c[c1][c2][1][j] =-(*(gx+i0)).c[c2][c1][1][j-1];
		  }
		}
	      }
#elif SU3_RECONSTRUCT_D == 12
	      for(int c2=0;c2<2;c2++){
		for(int c1=0;c1<3;c1++){
		  g.c[c1][c2][0][0] =0;
		  g.c[c1][c2][1][0] =0;
		  for(int j=1;j<VLEND;j++){
		    g.c[c1][c2][0][j] = (*(gx+i0)).c[c2][c1][0][j-1];
		    g.c[c1][c2][1][j] =-(*(gx+i0)).c[c2][c1][1][j-1];
		  }
		}
	      }
	      for(int j=0;j<VLEND;j++){
		g.c[0][2][0][j] = g.c[1][0][0][j]*g.c[2][1][0][j]-g.c[1][0][1][j]*g.c[2][1][1][j]-g.c[2][0][0][j]*g.c[1][1][0][j]+g.c[2][0][1][j]*g.c[1][1][1][j];
		g.c[1][2][0][j] = g.c[2][0][0][j]*g.c[0][1][0][j]-g.c[2][0][1][j]*g.c[0][1][1][j]-g.c[0][0][0][j]*g.c[2][1][0][j]+g.c[0][0][1][j]*g.c[2][1][1][j];
		g.c[2][2][0][j] = g.c[0][0][0][j]*g.c[1][1][0][j]-g.c[0][0][1][j]*g.c[1][1][1][j]-g.c[1][0][0][j]*g.c[0][1][0][j]+g.c[1][0][1][j]*g.c[0][1][1][j];
		g.c[0][2][1][j] =-g.c[1][0][0][j]*g.c[2][1][1][j]-g.c[1][0][1][j]*g.c[2][1][0][j]+g.c[2][0][0][j]*g.c[1][1][1][j]+g.c[2][0][1][j]*g.c[1][1][0][j];
		g.c[1][2][1][j] =-g.c[2][0][0][j]*g.c[0][1][1][j]-g.c[2][0][1][j]*g.c[0][1][0][j]+g.c[0][0][0][j]*g.c[2][1][1][j]+g.c[0][0][1][j]*g.c[2][1][0][j];
		g.c[2][2][1][j] =-g.c[0][0][0][j]*g.c[1][1][1][j]-g.c[0][0][1][j]*g.c[1][1][0][j]+g.c[1][0][0][j]*g.c[0][1][1][j]+g.c[1][0][1][j]*g.c[0][1][0][j];
	      }
#endif
	      /*
	      for(int c2=0;c2<3;c2++){
	    	for(int c1=0;c1<3;c1++){
		  g.c[c1][c2][0][0] =0;
		  g.c[c1][c2][1][0] =0;
	    	  for(int j=1;j<VLEND;j++){
	    	    g.c[c1][c2][0][j] = (*(gx+i0)).c[c2][c1][0][j-1];
	    	    g.c[c1][c2][1][j] =-(*(gx+i0)).c[c2][c1][1][j-1];
	    	  }
	    	}
	      }
	      */
	      for(int c=0;c<3;c++){
	    	for (int s=0;s<4;s++) {
	    	  for (int ri=0;ri<2;ri++) {
	    	    v[s][ri][0]=0; 
	    	    for(int j=1;j<VLEND;j++){
	    	      v[s][ri][j]=in[i0].c[c][s][ri][j-1]; 
	    	    }
	    	  }
	    	}
	    	for(int j=0;j<VLEND;j++){
	    	  a.c[c][0][0][j] = v[0][0][j] - v[3][1][j];
	    	  a.c[c][1][0][j] = v[1][0][j] - v[2][1][j];
	    	  a.c[c][0][1][j] = v[0][1][j] + v[3][0][j];
	    	  a.c[c][1][1][j] = v[1][1][j] + v[2][0][j];
	    	}
	      }
	      __m_fg(b.c, a.c, g.c);
	      __s_xb(tmp.c, b.c);
	    }
	    if(y != ny-1 ) __mult_yfwd(tmp.c, gy, in, i0, iyf);
	    if(y != 0    ) __mult_ybwd(tmp.c, gy, in,     iyb);
	    if(z != nz-1 ) __mult_zfwd(tmp.c, gz, in, i0, izf);
	    if(z != 0    ) __mult_zbwd(tmp.c, gz, in,     izb);
	    if(t != nt-1 ) __mult_tfwd(tmp.c, gt, in, i0, itf);
	    if(t != 0    ) __mult_tbwd(tmp.c, gt, in,     itb);
            __mult_clvd( tmp.cv, clvd[i0 + vold*(*DEO)].cv);
            __load_sc(out[i0].c, in[i0].c);
            __store_add_sc(out[i0].c, tmp.c, mkappa);
	  }
	}
      }
    }
  }
#ifdef __cplusplus
}
#endif
