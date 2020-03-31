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

#ifdef __cplusplus
extern "C"{
#endif

  extern int std_xyzt2i_(int*);
  extern int e_o_(int*);
  extern __attribute__((aligned(64))) pglud_t glud;
  extern __attribute__((aligned(64))) pglus_t glus;
  extern __attribute__((aligned(64))) pclvd_t clvd;
  extern __attribute__((aligned(64))) pclvs_t clvs;
  extern double kappa, kappa2, mkappa;
  extern int nt, nz, ny, nx, nxh, nxd, vold, nxs, vols;
  

  void qws_loadg_bqcd_(double *u, int* volh_tot, double* in_kappa){
    int x, y, z, t, j[4], mu, col1, col2;
    kappa = *in_kappa;
    mkappa= -kappa;
    kappa2= kappa * kappa;
    for(t=0; t<nt; t++){
      for(z=0; z<nz; z++){
	for(y=0; y<ny; y++){
	  for(x=0; x<nx; x++){
	    j[0] = x;
	    j[1] = y;
	    j[2] = z;
	    j[3] = t;
	    int i = std_xyzt2i_(j);
	    int eo = e_o_(j);
	    int ihh = (i-1)*18 +eo*(*volh_tot)*18;
	    int iwsd=  x/2 /VLEND + nxd*y + nxd*ny*z + nxd*ny*nz*t + NDIM*vold*eo;
	    int ixxd= (x/2)%VLEND;
	    int iwss=  x/2 /VLENS + nxs*y + nxs*ny*z + nxs*ny*nz*t + NDIM*vols*eo;
	    int ixxs= (x/2)%VLENS;
	    for(mu=0; mu<NDIM; mu++){
#if SU3_RECONSTRUCT_D == 18
	      for(col1=0; col1<NCOL; col1++){
		for(col2=0; col2<NCOL; col2++){
		  glud[iwsd+ vold*mu].c[col2][col1][0][ixxd] =        u[0 + 2*col1 + 6*col2 + ihh +mu*(*volh_tot)*36];
		  glud[iwsd+ vold*mu].c[col2][col1][1][ixxd] =        u[1 + 2*col1 + 6*col2 + ihh +mu*(*volh_tot)*36];
		}
	      }
#elif SU3_RECONSTRUCT_D == 12
	      for(col1=0; col1<NCOL; col1++){
		for(col2=0; col2<NCOL-1; col2++){
		  glud[iwsd+ vold*mu].c[col2][col1][0][ixxd] =        u[0 + 2*col1 + 6*col2 + ihh +mu*(*volh_tot)*36];
		  glud[iwsd+ vold*mu].c[col2][col1][1][ixxd] =        u[1 + 2*col1 + 6*col2 + ihh +mu*(*volh_tot)*36];
		}
	      }
#endif
#if SU3_RECONSTRUCT_S == 18
	      for(col1=0; col1<NCOL; col1++){
		for(col2=0; col2<NCOL; col2++){
		  glus[iwss+ vols*mu].c[col2][col1][0][ixxs] = (float)u[0 + 2*col1 + 6*col2 + ihh +mu*(*volh_tot)*36];
		  glus[iwss+ vols*mu].c[col2][col1][1][ixxs] = (float)u[1 + 2*col1 + 6*col2 + ihh +mu*(*volh_tot)*36];
		}
	      }
#elif SU3_RECONSTRUCT_S == 12
	      for(col1=0; col1<NCOL; col1++){
		for(col2=0; col2<NCOL-1; col2++){
		  glus[iwss+ vols*mu].c[col2][col1][0][ixxs] = (float)u[0 + 2*col1 + 6*col2 + ihh +mu*(*volh_tot)*36];
		  glus[iwss+ vols*mu].c[col2][col1][1][ixxs] = (float)u[1 + 2*col1 + 6*col2 + ihh +mu*(*volh_tot)*36];
		}
	      }
#endif

	    }
	  }
	}
      }
    }
  }

  void qws_loadg_dd_(double *u, int* volh_tot, double* in_kappa){
    int x, y, z, t, j[4], mu, col1, col2;
    kappa = *in_kappa;
    mkappa= -kappa;
    kappa2= kappa * kappa;
    for(t=0; t<nt; t++){
      for(z=0; z<nz; z++){
	for(y=0; y<ny; y++){
	  for(x=0; x<nx; x++){
	    j[0] = x;
	    j[1] = y;
	    j[2] = z;
	    j[3] = t;
	    int i = std_xyzt2i_(j);
	    int eo = e_o_(j);
	    int ihh = (i-1)*18 +eo*(*volh_tot)*18;
	    int iwsd= (x%nxh)/VLEND + nxd*y + nxd*ny*z + nxd*ny*nz*t + NDIM*vold*(x/nxh);
	    int ixxd= (x%nxh)%VLEND;
	    int iwss= (x%nxh)/VLENS + nxs*y + nxs*ny*z + nxs*ny*nz*t + NDIM*vols*(x/nxh);
	    int ixxs= (x%nxh)%VLENS;
	    //printf("%d %d %d %d    %d %d \n",x, y, z, t,iwsd, ixxd );
	    for(mu=0; mu<NDIM; mu++){
#if SU3_RECONSTRUCT_D == 18
	      for(col1=0; col1<NCOL; col1++){
	    	for(col2=0; col2<NCOL; col2++){
	    	  glud[iwsd+ vold*mu].c[col2][col1][0][ixxd] =        u[0 + 2*col1 + 6*col2 + ihh +mu*(*volh_tot)*36];
	    	  glud[iwsd+ vold*mu].c[col2][col1][1][ixxd] =        u[1 + 2*col1 + 6*col2 + ihh +mu*(*volh_tot)*36];
		}
	      }
#elif SU3_RECONSTRUCT_D == 12
	      for(col1=0; col1<NCOL; col1++){
	    	for(col2=0; col2<NCOL-1; col2++){
	    	  glud[iwsd+ vold*mu].c[col2][col1][0][ixxd] =        u[0 + 2*col1 + 6*col2 + ihh +mu*(*volh_tot)*36];
	    	  glud[iwsd+ vold*mu].c[col2][col1][1][ixxd] =        u[1 + 2*col1 + 6*col2 + ihh +mu*(*volh_tot)*36];
		}
	      }
#endif
#if SU3_RECONSTRUCT_S == 18
	      for(col1=0; col1<NCOL; col1++){
	    	for(col2=0; col2<NCOL; col2++){
	    	  glus[iwss+ vols*mu].c[col2][col1][0][ixxs] = (float)u[0 + 2*col1 + 6*col2 + ihh +mu*(*volh_tot)*36];
	    	  glus[iwss+ vols*mu].c[col2][col1][1][ixxs] = (float)u[1 + 2*col1 + 6*col2 + ihh +mu*(*volh_tot)*36];
	    	}
	      }
#elif SU3_RECONSTRUCT_S == 12
	      for(col1=0; col1<NCOL; col1++){
	    	for(col2=0; col2<NCOL-1; col2++){
	    	  glus[iwss+ vols*mu].c[col2][col1][0][ixxs] = (float)u[0 + 2*col1 + 6*col2 + ihh +mu*(*volh_tot)*36];
	    	  glus[iwss+ vols*mu].c[col2][col1][1][ixxs] = (float)u[1 + 2*col1 + 6*col2 + ihh +mu*(*volh_tot)*36];
	    	}
	      }
#endif
	    }
	  }
	}
      }
    }
  }

  void qws_loadc_bqcd_(double *c, int* volh){
    int x, y, z, t, j[4], ud;
    for(t=0; t<nt; t++){
      for(z=0; z<nz; z++){
	for(y=0; y<ny; y++){
	  for(x=0; x<nx; x++){
	    j[0] = x;
	    j[1] = y;
	    j[2] = z;
	    j[3] = t;
	    int i = std_xyzt2i_(j);
	    int eo = e_o_(j);
	    int ihh = (i-1)*72 +eo*(*volh)*72;
	    int iwsd=  x/2 /VLEND + nxd*y + nxd*ny*z + nxd*ny*nz*t + vold*eo;
	    int ixxd= (x/2)%VLEND;
	    int iwss=  x/2 /VLENS + nxs*y + nxs*ny*z + nxs*ny*nz*t + vols*eo;
	    int ixxs= (x/2)%VLENS;
	    for(ud=0;ud<2;ud++){
	      clvd[iwsd].c[ud][ 0][ixxd] = c[ 0 + ud*36 + ihh];
	      clvd[iwsd].c[ud][ 1][ixxd] = c[ 1 + ud*36 + ihh];
	      clvd[iwsd].c[ud][ 2][ixxd] = c[20 + ud*36 + ihh];
	      clvd[iwsd].c[ud][ 3][ixxd] = c[21 + ud*36 + ihh];
	      clvd[iwsd].c[ud][ 4][ixxd] = c[32 + ud*36 + ihh];
	      clvd[iwsd].c[ud][ 5][ixxd] = c[33 + ud*36 + ihh];
	      clvd[iwsd].c[ud][ 6][ixxd] = c[ 2 + ud*36 + ihh];
	      clvd[iwsd].c[ud][ 7][ixxd] = c[ 3 + ud*36 + ihh];
	      clvd[iwsd].c[ud][ 8][ixxd] = c[ 4 + ud*36 + ihh];
	      clvd[iwsd].c[ud][ 9][ixxd] = c[ 5 + ud*36 + ihh];
	      clvd[iwsd].c[ud][10][ixxd] = c[ 6 + ud*36 + ihh];
	      clvd[iwsd].c[ud][11][ixxd] = c[ 7 + ud*36 + ihh];
	      clvd[iwsd].c[ud][12][ixxd] = c[ 8 + ud*36 + ihh];
	      clvd[iwsd].c[ud][13][ixxd] = c[ 9 + ud*36 + ihh];
	      clvd[iwsd].c[ud][14][ixxd] = c[10 + ud*36 + ihh];
	      clvd[iwsd].c[ud][15][ixxd] = c[11 + ud*36 + ihh];
	      clvd[iwsd].c[ud][16][ixxd] = c[12 + ud*36 + ihh];
	      clvd[iwsd].c[ud][17][ixxd] = c[13 + ud*36 + ihh];
	      clvd[iwsd].c[ud][18][ixxd] = c[14 + ud*36 + ihh];
	      clvd[iwsd].c[ud][19][ixxd] = c[15 + ud*36 + ihh];
	      clvd[iwsd].c[ud][20][ixxd] = c[16 + ud*36 + ihh];
	      clvd[iwsd].c[ud][21][ixxd] = c[17 + ud*36 + ihh];
	      clvd[iwsd].c[ud][22][ixxd] = c[18 + ud*36 + ihh];
	      clvd[iwsd].c[ud][23][ixxd] = c[19 + ud*36 + ihh];
	      clvd[iwsd].c[ud][24][ixxd] = c[22 + ud*36 + ihh];
	      clvd[iwsd].c[ud][25][ixxd] = c[23 + ud*36 + ihh];
	      clvd[iwsd].c[ud][26][ixxd] = c[24 + ud*36 + ihh];
	      clvd[iwsd].c[ud][27][ixxd] = c[25 + ud*36 + ihh];
	      clvd[iwsd].c[ud][28][ixxd] = c[26 + ud*36 + ihh];
	      clvd[iwsd].c[ud][29][ixxd] = c[27 + ud*36 + ihh];
	      clvd[iwsd].c[ud][30][ixxd] = c[28 + ud*36 + ihh];
	      clvd[iwsd].c[ud][31][ixxd] = c[29 + ud*36 + ihh];
	      clvd[iwsd].c[ud][32][ixxd] = c[30 + ud*36 + ihh];
	      clvd[iwsd].c[ud][33][ixxd] = c[31 + ud*36 + ihh];
	      clvd[iwsd].c[ud][34][ixxd] = c[34 + ud*36 + ihh];
	      clvd[iwsd].c[ud][35][ixxd] = c[35 + ud*36 + ihh];
	    }
	    // single
	    for(ud=0;ud<2;ud++){
	      clvs[iwss].c[ud][ 0][ixxs] = (float)c[ 0 + ud*36 + ihh];
	      clvs[iwss].c[ud][ 1][ixxs] = (float)c[ 1 + ud*36 + ihh];
	      clvs[iwss].c[ud][ 2][ixxs] = (float)c[20 + ud*36 + ihh];
	      clvs[iwss].c[ud][ 3][ixxs] = (float)c[21 + ud*36 + ihh];
	      clvs[iwss].c[ud][ 4][ixxs] = (float)c[32 + ud*36 + ihh];
	      clvs[iwss].c[ud][ 5][ixxs] = (float)c[33 + ud*36 + ihh];
	      clvs[iwss].c[ud][ 6][ixxs] = (float)c[ 2 + ud*36 + ihh];
	      clvs[iwss].c[ud][ 7][ixxs] = (float)c[ 3 + ud*36 + ihh];
	      clvs[iwss].c[ud][ 8][ixxs] = (float)c[ 4 + ud*36 + ihh];
	      clvs[iwss].c[ud][ 9][ixxs] = (float)c[ 5 + ud*36 + ihh];
	      clvs[iwss].c[ud][10][ixxs] = (float)c[ 6 + ud*36 + ihh];
	      clvs[iwss].c[ud][11][ixxs] = (float)c[ 7 + ud*36 + ihh];
	      clvs[iwss].c[ud][12][ixxs] = (float)c[ 8 + ud*36 + ihh];
	      clvs[iwss].c[ud][13][ixxs] = (float)c[ 9 + ud*36 + ihh];
	      clvs[iwss].c[ud][14][ixxs] = (float)c[10 + ud*36 + ihh];
	      clvs[iwss].c[ud][15][ixxs] = (float)c[11 + ud*36 + ihh];
	      clvs[iwss].c[ud][16][ixxs] = (float)c[12 + ud*36 + ihh];
	      clvs[iwss].c[ud][17][ixxs] = (float)c[13 + ud*36 + ihh];
	      clvs[iwss].c[ud][18][ixxs] = (float)c[14 + ud*36 + ihh];
	      clvs[iwss].c[ud][19][ixxs] = (float)c[15 + ud*36 + ihh];
	      clvs[iwss].c[ud][20][ixxs] = (float)c[16 + ud*36 + ihh];
	      clvs[iwss].c[ud][21][ixxs] = (float)c[17 + ud*36 + ihh];
	      clvs[iwss].c[ud][22][ixxs] = (float)c[18 + ud*36 + ihh];
	      clvs[iwss].c[ud][23][ixxs] = (float)c[19 + ud*36 + ihh];
	      clvs[iwss].c[ud][24][ixxs] = (float)c[22 + ud*36 + ihh];
	      clvs[iwss].c[ud][25][ixxs] = (float)c[23 + ud*36 + ihh];
	      clvs[iwss].c[ud][26][ixxs] = (float)c[24 + ud*36 + ihh];
	      clvs[iwss].c[ud][27][ixxs] = (float)c[25 + ud*36 + ihh];
	      clvs[iwss].c[ud][28][ixxs] = (float)c[26 + ud*36 + ihh];
	      clvs[iwss].c[ud][29][ixxs] = (float)c[27 + ud*36 + ihh];
	      clvs[iwss].c[ud][30][ixxs] = (float)c[28 + ud*36 + ihh];
	      clvs[iwss].c[ud][31][ixxs] = (float)c[29 + ud*36 + ihh];
	      clvs[iwss].c[ud][32][ixxs] = (float)c[30 + ud*36 + ihh];
	      clvs[iwss].c[ud][33][ixxs] = (float)c[31 + ud*36 + ihh];
	      clvs[iwss].c[ud][34][ixxs] = (float)c[34 + ud*36 + ihh];
	      clvs[iwss].c[ud][35][ixxs] = (float)c[35 + ud*36 + ihh];
	    }
	  }
	}
      }
    }
  }

  void qws_loadc_dd_(double *c, int* volh){
    int x, y, z, t, j[4], ud;
    for(t=0; t<nt; t++){
      for(z=0; z<nz; z++){
	for(y=0; y<ny; y++){
	  for(x=0; x<nx; x++){
	    j[0] = x;
	    j[1] = y;
	    j[2] = z;
	    j[3] = t;
	    int i = std_xyzt2i_(j);
	    int eo = e_o_(j);
	    int ihh = (i-1)*72 +eo*(*volh)*72;
	    int iwsd= (x%nxh)/VLEND + nxd*y + nxd*ny*z + nxd*ny*nz*t + vold*(x/nxh);
	    int ixxd= (x%nxh)%VLEND;
	    int iwss= (x%nxh)/VLENS + nxs*y + nxs*ny*z + nxs*ny*nz*t + vols*(x/nxh);
	    int ixxs= (x%nxh)%VLENS;
	    for(ud=0;ud<2;ud++){
	      clvd[iwsd].c[ud][ 0][ixxd] = c[ 0 + ud*36 + ihh];
	      clvd[iwsd].c[ud][ 1][ixxd] = c[ 1 + ud*36 + ihh];
	      clvd[iwsd].c[ud][ 2][ixxd] = c[20 + ud*36 + ihh];
	      clvd[iwsd].c[ud][ 3][ixxd] = c[21 + ud*36 + ihh];
	      clvd[iwsd].c[ud][ 4][ixxd] = c[32 + ud*36 + ihh];
	      clvd[iwsd].c[ud][ 5][ixxd] = c[33 + ud*36 + ihh];
	      clvd[iwsd].c[ud][ 6][ixxd] = c[ 2 + ud*36 + ihh];
	      clvd[iwsd].c[ud][ 7][ixxd] = c[ 3 + ud*36 + ihh];
	      clvd[iwsd].c[ud][ 8][ixxd] = c[ 4 + ud*36 + ihh];
	      clvd[iwsd].c[ud][ 9][ixxd] = c[ 5 + ud*36 + ihh];
	      clvd[iwsd].c[ud][10][ixxd] = c[ 6 + ud*36 + ihh];
	      clvd[iwsd].c[ud][11][ixxd] = c[ 7 + ud*36 + ihh];
	      clvd[iwsd].c[ud][12][ixxd] = c[ 8 + ud*36 + ihh];
	      clvd[iwsd].c[ud][13][ixxd] = c[ 9 + ud*36 + ihh];
	      clvd[iwsd].c[ud][14][ixxd] = c[10 + ud*36 + ihh];
	      clvd[iwsd].c[ud][15][ixxd] = c[11 + ud*36 + ihh];
	      clvd[iwsd].c[ud][16][ixxd] = c[12 + ud*36 + ihh];
	      clvd[iwsd].c[ud][17][ixxd] = c[13 + ud*36 + ihh];
	      clvd[iwsd].c[ud][18][ixxd] = c[14 + ud*36 + ihh];
	      clvd[iwsd].c[ud][19][ixxd] = c[15 + ud*36 + ihh];
	      clvd[iwsd].c[ud][20][ixxd] = c[16 + ud*36 + ihh];
	      clvd[iwsd].c[ud][21][ixxd] = c[17 + ud*36 + ihh];
	      clvd[iwsd].c[ud][22][ixxd] = c[18 + ud*36 + ihh];
	      clvd[iwsd].c[ud][23][ixxd] = c[19 + ud*36 + ihh];
	      clvd[iwsd].c[ud][24][ixxd] = c[22 + ud*36 + ihh];
	      clvd[iwsd].c[ud][25][ixxd] = c[23 + ud*36 + ihh];
	      clvd[iwsd].c[ud][26][ixxd] = c[24 + ud*36 + ihh];
	      clvd[iwsd].c[ud][27][ixxd] = c[25 + ud*36 + ihh];
	      clvd[iwsd].c[ud][28][ixxd] = c[26 + ud*36 + ihh];
	      clvd[iwsd].c[ud][29][ixxd] = c[27 + ud*36 + ihh];
	      clvd[iwsd].c[ud][30][ixxd] = c[28 + ud*36 + ihh];
	      clvd[iwsd].c[ud][31][ixxd] = c[29 + ud*36 + ihh];
	      clvd[iwsd].c[ud][32][ixxd] = c[30 + ud*36 + ihh];
	      clvd[iwsd].c[ud][33][ixxd] = c[31 + ud*36 + ihh];
	      clvd[iwsd].c[ud][34][ixxd] = c[34 + ud*36 + ihh];
	      clvd[iwsd].c[ud][35][ixxd] = c[35 + ud*36 + ihh];
	    }
	    // single
	    for(ud=0;ud<2;ud++){
	      clvs[iwss].c[ud][ 0][ixxs] = (float)c[ 0 + ud*36 + ihh];
	      clvs[iwss].c[ud][ 1][ixxs] = (float)c[ 1 + ud*36 + ihh];
	      clvs[iwss].c[ud][ 2][ixxs] = (float)c[20 + ud*36 + ihh];
	      clvs[iwss].c[ud][ 3][ixxs] = (float)c[21 + ud*36 + ihh];
	      clvs[iwss].c[ud][ 4][ixxs] = (float)c[32 + ud*36 + ihh];
	      clvs[iwss].c[ud][ 5][ixxs] = (float)c[33 + ud*36 + ihh];
	      clvs[iwss].c[ud][ 6][ixxs] = (float)c[ 2 + ud*36 + ihh];
	      clvs[iwss].c[ud][ 7][ixxs] = (float)c[ 3 + ud*36 + ihh];
	      clvs[iwss].c[ud][ 8][ixxs] = (float)c[ 4 + ud*36 + ihh];
	      clvs[iwss].c[ud][ 9][ixxs] = (float)c[ 5 + ud*36 + ihh];
	      clvs[iwss].c[ud][10][ixxs] = (float)c[ 6 + ud*36 + ihh];
	      clvs[iwss].c[ud][11][ixxs] = (float)c[ 7 + ud*36 + ihh];
	      clvs[iwss].c[ud][12][ixxs] = (float)c[ 8 + ud*36 + ihh];
	      clvs[iwss].c[ud][13][ixxs] = (float)c[ 9 + ud*36 + ihh];
	      clvs[iwss].c[ud][14][ixxs] = (float)c[10 + ud*36 + ihh];
	      clvs[iwss].c[ud][15][ixxs] = (float)c[11 + ud*36 + ihh];
	      clvs[iwss].c[ud][16][ixxs] = (float)c[12 + ud*36 + ihh];
	      clvs[iwss].c[ud][17][ixxs] = (float)c[13 + ud*36 + ihh];
	      clvs[iwss].c[ud][18][ixxs] = (float)c[14 + ud*36 + ihh];
	      clvs[iwss].c[ud][19][ixxs] = (float)c[15 + ud*36 + ihh];
	      clvs[iwss].c[ud][20][ixxs] = (float)c[16 + ud*36 + ihh];
	      clvs[iwss].c[ud][21][ixxs] = (float)c[17 + ud*36 + ihh];
	      clvs[iwss].c[ud][22][ixxs] = (float)c[18 + ud*36 + ihh];
	      clvs[iwss].c[ud][23][ixxs] = (float)c[19 + ud*36 + ihh];
	      clvs[iwss].c[ud][24][ixxs] = (float)c[22 + ud*36 + ihh];
	      clvs[iwss].c[ud][25][ixxs] = (float)c[23 + ud*36 + ihh];
	      clvs[iwss].c[ud][26][ixxs] = (float)c[24 + ud*36 + ihh];
	      clvs[iwss].c[ud][27][ixxs] = (float)c[25 + ud*36 + ihh];
	      clvs[iwss].c[ud][28][ixxs] = (float)c[26 + ud*36 + ihh];
	      clvs[iwss].c[ud][29][ixxs] = (float)c[27 + ud*36 + ihh];
	      clvs[iwss].c[ud][30][ixxs] = (float)c[28 + ud*36 + ihh];
	      clvs[iwss].c[ud][31][ixxs] = (float)c[29 + ud*36 + ihh];
	      clvs[iwss].c[ud][32][ixxs] = (float)c[30 + ud*36 + ihh];
	      clvs[iwss].c[ud][33][ixxs] = (float)c[31 + ud*36 + ihh];
	      clvs[iwss].c[ud][34][ixxs] = (float)c[34 + ud*36 + ihh];
	      clvs[iwss].c[ud][35][ixxs] = (float)c[35 + ud*36 + ihh];
	    }
	  }
	}
      }
    }
  }

  void qws_loadfd_bqcd_(scd_t *out, double *in){
    int x, y, z, t, j[4], s, c;
    for(t=0; t<nt; t++){
      for(z=0; z<nz; z++){
	for(y=0; y<ny; y++){
	  for(x=0; x<nx; x++){
	    j[0] = x;
	    j[1] = y;
	    j[2] = z;
	    j[3] = t;
	    int i = std_xyzt2i_(j);
	    int eo = e_o_(j);
	    if (eo == 0) {
	      int ihh = (i-1)*24;
	      int iws =  x/2 /VLEND + nxd*y + nxd*ny*z + nxd*ny*nz*t;
	      int ixx = (x/2)%VLEND;
	      for(c=0; c<NCOL; c++){
		for(s=0; s<NDIM; s++){
		  out[iws].c[c][s][0][ixx] = in[0 + 2*s + 8*c + ihh];
		  out[iws].c[c][s][1][ixx] = in[1 + 2*s + 8*c + ihh];
		}
	      }
	    }
	  }
	}
      }
    }
  }
  
  void qws_loadfs_bqcd_(scs_t *out, double *in){
    int x, y, z, t, j[4], s, c;
    for(t=0; t<nt; t++){
      for(z=0; z<nz; z++){
	for(y=0; y<ny; y++){
	  for(x=0; x<nx; x++){
	    j[0] = x;
	    j[1] = y;
	    j[2] = z;
	    j[3] = t;
	    int i = std_xyzt2i_(j);
	    int eo = e_o_(j);
	    if (eo == 0) {
	      int ihh = (i-1)*24;
	      int iws =  x/2 /VLENS + nxs*y + nxs*ny*z + nxs*ny*nz*t;
	      int ixx = (x/2)%VLENS;
	      for(c=0; c<NCOL; c++){
		for(s=0; s<NDIM; s++){
		  out[iws].c[c][s][0][ixx] = (float)in[0 + 2*s + 8*c + ihh];
		  out[iws].c[c][s][1][ixx] = (float)in[1 + 2*s + 8*c + ihh];
		}
	      }
	    }
	  }
	}
      }
    }
  }
  
  void qws_storefd_bqcd_(double *out, scd_t *in){
    int x, y, z, t, j[4], s, c;
    for(t=0; t<nt; t++){
      for(z=0; z<nz; z++){
	for(y=0; y<ny; y++){
	  for(x=0; x<nx; x++){
	    j[0] = x;
	    j[1] = y;
	    j[2] = z;
	    j[3] = t;
	    int i = std_xyzt2i_(j);
	    int eo = e_o_(j);
	    if (eo == 0) {
	      int ihh = (i-1)*24;
	      int iws =  x/2 /VLEND + nxd*y + nxd*ny*z + nxd*ny*nz*t;
	      int ixx = (x/2)%VLEND;
	      for(c=0; c<NCOL; c++){
		for(s=0; s<NDIM; s++){
		  out[0 + 2*s + 8*c + ihh] = in[iws].c[c][s][0][ixx];
		  out[1 + 2*s + 8*c + ihh] = in[iws].c[c][s][1][ixx];
		}
	      }
	    }
	  }
	}
      }
    }
  }
  
  void qws_storefs_bqcd_(double *out, scs_t *in){
    int x, y, z, t, j[4], s, c;
    for(t=0; t<nt; t++){
      for(z=0; z<nz; z++){
	for(y=0; y<ny; y++){
	  for(x=0; x<nx; x++){
	    j[0] = x;
	    j[1] = y;
	    j[2] = z;
	    j[3] = t;
	    int i = std_xyzt2i_(j);
	    int eo = e_o_(j);
	    if (eo == 0) {
	      int ihh = (i-1)*24;
	      int iws =  x/2 /VLENS + nxs*y + nxs*ny*z + nxs*ny*nz*t;
	      int ixx = (x/2)%VLENS;
	      for(c=0; c<NCOL; c++){
		for(s=0; s<NDIM; s++){
		  out[0 + 2*s + 8*c + ihh] = (double)in[iws].c[c][s][0][ixx];
		  out[1 + 2*s + 8*c + ihh] = (double)in[iws].c[c][s][1][ixx];
		}
	      }
	    }
	  }
	}
      }
    }
  }
  
  void qws_loadfd_dd_(scd_t *out, double *ine, double *ino){
    int x, y, z, t, j[4], s, c;
    for(t=0; t<nt; t++){
      for(z=0; z<nz; z++){
	for(y=0; y<ny; y++){
	  for(x=0; x<nx; x++){
	    j[0] = x;
	    j[1] = y;
	    j[2] = z;
	    j[3] = t;
	    int i = std_xyzt2i_(j);
	    int eo = e_o_(j);
	    int ihh = (i-1)*24;
	    int iws = (x%nxh)/VLEND + nxd*y + nxd*ny*z + nxd*ny*nz*t + vold*(x/nxh);
	    int ixx = (x%nxh)%VLEND;
	    if (eo == 0) {
	      for(c=0; c<NCOL; c++){
		for(s=0; s<NDIM; s++){
		  out[iws].c[c][s][0][ixx] = ine[0 + 2*s + 8*c + ihh];
		  out[iws].c[c][s][1][ixx] = ine[1 + 2*s + 8*c + ihh];
		}
	      }
	    } else {
	      for(c=0; c<NCOL; c++){
		for(s=0; s<NDIM; s++){
		  out[iws].c[c][s][0][ixx] = ino[0 + 2*s + 8*c + ihh];
		  out[iws].c[c][s][1][ixx] = ino[1 + 2*s + 8*c + ihh];
		}
	      }
	    }
	  }
	}
      }
    }
  }

  void qws_storefd_dd_(double *oute, double *outo, scd_t *in){
    int x, y, z, t, j[4], s, c;
    for(t=0; t<nt; t++){
      for(z=0; z<nz; z++){
	for(y=0; y<ny; y++){
	  for(x=0; x<nx; x++){
	    j[0] = x;
	    j[1] = y;
	    j[2] = z;
	    j[3] = t;
	    int i = std_xyzt2i_(j);
	    int eo = e_o_(j);
	    int ihh = (i-1)*24;
	    int iws = (x%nxh)/VLEND + nxd*y + nxd*ny*z + nxd*ny*nz*t + vold*(x/nxh);
	    int ixx = (x%nxh)%VLEND;
	    if (eo == 0) {
	      for(c=0; c<NCOL; c++){
		for(s=0; s<NDIM; s++){
		  oute[0 + 2*s + 8*c + ihh] = in[iws].c[c][s][0][ixx];
		  oute[1 + 2*s + 8*c + ihh] = in[iws].c[c][s][1][ixx];
		}
	      }
	    } else {
	      for(c=0; c<NCOL; c++){
		for(s=0; s<NDIM; s++){
		  outo[0 + 2*s + 8*c + ihh] = in[iws].c[c][s][0][ixx];
		  outo[1 + 2*s + 8*c + ihh] = in[iws].c[c][s][1][ixx];
		}
	      }
	    }
	  }
	}
      }
    }
  }

  void qws_loadfs_dd_(scs_t *out, double *ine, double *ino){
    int x, y, z, t, j[4], s, c;
    for(t=0; t<nt; t++){
      for(z=0; z<nz; z++){
	for(y=0; y<ny; y++){
	  for(x=0; x<nx; x++){
	    j[0] = x;
	    j[1] = y;
	    j[2] = z;
	    j[3] = t;
	    int i = std_xyzt2i_(j);
	    int eo = e_o_(j);
	    int ihh = (i-1)*24;
	    int iws = (x%nxh)/VLENS + nxs*y + nxs*ny*z + nxs*ny*nz*t + vols*(x/nxh);
	    int ixx = (x%nxh)%VLENS;
	    if (eo == 0) {
	      for(c=0; c<NCOL; c++){
		for(s=0; s<NDIM; s++){
		  out[iws].c[c][s][0][ixx] = (float)ine[0 + 2*s + 8*c + ihh];
		  out[iws].c[c][s][1][ixx] = (float)ine[1 + 2*s + 8*c + ihh];
		}
	      }
	    } else {
	      for(c=0; c<NCOL; c++){
		for(s=0; s<NDIM; s++){
		  out[iws].c[c][s][0][ixx] = (float)ino[0 + 2*s + 8*c + ihh];
		  out[iws].c[c][s][1][ixx] = (float)ino[1 + 2*s + 8*c + ihh];
		}
	      }
	    }
	  }
	}
      }
    }
  }

  void qws_storefs_dd_(double *oute, double *outo, scs_t *in){
    int x, y, z, t, j[4], s, c;
    for(t=0; t<nt; t++){
      for(z=0; z<nz; z++){
	for(y=0; y<ny; y++){
	  for(x=0; x<nx; x++){
	    j[0] = x;
	    j[1] = y;
	    j[2] = z;
	    j[3] = t;
	    int i = std_xyzt2i_(j);
	    int eo = e_o_(j);
	    int ihh = (i-1)*24;
	    int iws = (x%nxh)/VLENS + nxs*y + nxs*ny*z + nxs*ny*nz*t + vols*(x/nxh);
	    int ixx = (x%nxh)%VLENS;
	    if (eo == 0) {
	      for(c=0; c<NCOL; c++){
		for(s=0; s<NDIM; s++){
		  oute[0 + 2*s + 8*c + ihh] = (double)in[iws].c[c][s][0][ixx];
		  oute[1 + 2*s + 8*c + ihh] = (double)in[iws].c[c][s][1][ixx];
		}
	      }
	    } else {
	      for(c=0; c<NCOL; c++){
		for(s=0; s<NDIM; s++){
		  outo[0 + 2*s + 8*c + ihh] = (double)in[iws].c[c][s][0][ixx];
		  outo[1 + 2*s + 8*c + ihh] = (double)in[iws].c[c][s][1][ixx];
		}
	      }
	    }
	  }
	}
      }
    }
  }


#ifdef HALF_PREC
#include "qws_h.h"
  void bqcd_qws_h_init_(){qws_h_init_(glus,clvs);}
#endif

#ifdef __cplusplus
}
#endif
