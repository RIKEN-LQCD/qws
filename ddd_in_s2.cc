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
//#define PREFETCH
#define _PFI 1
#include "qws.h"
#include "clover_s.h"
#include "mult_all.h"
#include "prefetch.h"

#ifdef __cplusplus
extern "C"{
#endif


  extern __attribute__((aligned(_ALIGN_SIZE))) pglus_t glus;
  extern __attribute__((aligned(_ALIGN_SIZE))) pclvs_t clvs;
  extern double kappa, mkappa;
  extern int nxs, ny, nz, nt;
  extern int vols;
  //---------------------------------------------------------------------------------------- mult D in bulk
  void ddd_in_s_(scs_t* __restrict__ out, const scs_t* __restrict__ in, const int* DEO) {

    const __restrict__ pglus_t gx = &glus[vols*0 + NDIM*vols*(*DEO)];
    const __restrict__ pglus_t gy = &glus[vols*1 + NDIM*vols*(*DEO)];
    const __restrict__ pglus_t gz = &glus[vols*2 + NDIM*vols*(*DEO)];
    const __restrict__ pglus_t gt = &glus[vols*3 + NDIM*vols*(*DEO)];
    const __restrict__ pclvs_t cl = &clvs[              vols*(*DEO)];

        //
        // X
        //
#pragma omp parallel for collapse(3) schedule(static)
    for(int t = 0; t < nt; ++t) {
    for(int z = 0; z < nz; ++z) {
    for(int y = 0; y < ny; ++y) {

      for(int x = 0; x < nxs; ++x) {

        int i0  =  x    + nxs*y + nxs*ny*z + nxs*ny*nz*t;
        int ixf = (x+1) + nxs*y + nxs*ny*z + nxs*ny*nz*t;
        int ixb = (x-1) + nxs*y + nxs*ny*z + nxs*ny*nz*t;

#ifdef PREFETCH
	//int pfx = (i0+_PFI)/nxs;
#endif

        scs_t tmp0 __attribute__((aligned(_ALIGN_SIZE)));

        for (int c = 0; c < 3; c++) {
        for (int s = 0; s < 4; s++) {
          for (int ri = 0; ri < 2; ri++) {
            for (int j = 0; j<VLENS; j++) {
              tmp0.c[c][s][ri][j] = 0.0f;
            }
          }
        }
        }

        //
        // X-forward
        //
        {
          __attribute__((aligned(_ALIGN_SIZE))) projscs_t a,ua;
          float ff[3][4][2][VLENS*2] __attribute__((aligned(4)));

          if ( x != nxs-1 ) {

            for(int c = 0; c < 3; ++c) {
            for(int s = 0; s < 4; ++s) {
              for(int ri = 0; ri < 2; ++ri) { 

                for(int j = 0; j < VLENS; ++j) { 
                  ff[c][s][ri][j      ] = ((float(*)[4][2][VLENS])((in+i0)->c))[c][s][ri][j]; 
                  ff[c][s][ri][j+VLENS] = ((float(*)[4][2][VLENS])((in+ixf)->c))[c][s][ri][j];
                }
              }
            }
            }

          } else {

            for(int c = 0; c < 3; ++c) {
            for(int s = 0; s < 4; ++s) {
              for(int ri = 0; ri < 2; ++ri) { 

                for(int j = 0; j < VLENS; ++j) { 
                  ff[c][s][ri][j      ] = ((float(*)[4][2][VLENS])((in+i0)->c))[c][s][ri][j]; 
                  ff[c][s][ri][j+VLENS] = 0.0f;
                }
              }
            }
            }

          }

          __mult_x_forw_pre_2_(a,ff); // with forward simd-site-shift
          __mult_u_y_(ua,a,(*(gx + i0)));
          __mult_x_forw_pst_(tmp0,ua);

        }

        //
        // X-backward
        //
        {
          __attribute__((aligned(_ALIGN_SIZE))) projscs_t a,ua;
          float ff_in[3][4][2][VLENS*2] __attribute__((aligned(4)));
          float ff_g [3][3][2][VLENS*2] __attribute__((aligned(4)));

          if ( x != 0 ) {

            for(int c = 0; c < 3; ++c) {
            for(int s = 0; s < 4; ++s) {
              for(int ri = 0; ri < 2; ++ri) { 
                for(int j = 0; j < VLENS; ++j) { 
                  ff_in[c][s][ri][j      ] = ((float(*)[4][2][VLENS])((in+ixb)->c))[c][s][ri][j];
                  ff_in[c][s][ri][j+VLENS] = ((float(*)[4][2][VLENS])((in+i0)->c))[c][s][ri][j];
                }
              }
            }
            }

            //
            // load link U.  WITHOUT taking Hermitian-conjugate
            //
            for (int c2 = 0; c2 < 3; ++c2) {
            for (int c1 = 0; c1 < 3; ++c1) {
              for(int ri = 0; ri < 2; ++ri) { 
                for(int j = 0; j < VLENS; ++j) { 
                  ff_g[c2][c1][ri][j      ] = ((float(*)[3][2][VLENS])((gx+ixb)->c))[c2][c1][ri][j];
                  ff_g[c2][c1][ri][j+VLENS] = ((float(*)[3][2][VLENS])((gx+i0)->c))[c2][c1][ri][j];
                }
              }
            }
            }

          } else {

            for (int c = 0; c < 3; ++c) {
            for (int s = 0; s < 4; ++s) {
              for(int ri = 0; ri < 2; ++ri) { 
                for(int j = 0; j < VLENS; ++j) { 
                  ff_in[c][s][ri][j      ] = 0.0f;
                  ff_in[c][s][ri][j+VLENS] = ((float(*)[4][2][VLENS])((in+i0)->c))[c][s][ri][j];
                }
              }
            }
            }

            //
            // load link U.  WITHOUT taking Hermitian-conjugate
            //
            for (int c2 = 0; c2 < 3; ++c2) {
            for (int c1 = 0; c1 < 3; ++c1) {
              for(int ri = 0; ri < 2; ++ri) { 
                for(int j = 0; j < VLENS; ++j) { 
                  ff_g[c2][c1][ri][j      ] = 0.0f;
                  ff_g[c2][c1][ri][j+VLENS] = ((float(*)[3][2][VLENS])((gx+i0)->c))[c2][c1][ri][j];
                }
              }
            }
            }

          }

          __mult_x_back_pre_2_(a,ff_in);// with backward simd-site-shift
          __mult_udag_y_2_(ua,a,ff_g);  // with backward simd-site-shift
          __mult_x_back_pst_(tmp0,ua);

        }

	__prefetch_out(out, (i0+_PFI));

#pragma loop norecurrence
	for (int c = 0; c < 3; c++) {
          for (int s = 0; s < 4; s++) {
            for (int ri = 0; ri < 2; ri++) {
              for (int j = 0; j < VLENS; j++) {
                ((float(*)[4][2][VLENS])((out + i0)->c))[c][s][ri][j] = ((float(*)[4][2][VLENS])(tmp0.c))[c][s][ri][j] ;
              }
            }
          }
	}

      }

    }
    }
    }



        //
        // Y
        //
#pragma omp parallel for collapse(3) schedule(static)
    for(int t = 0; t < nt; ++t) {
    for(int z = 0; z < nz; ++z) {
    for(int y = 0; y < ny; ++y) {

      for(int x = 0; x < nxs; ++x) {

        int i0  =  x    + nxs*y + nxs*ny*z + nxs*ny*nz*t;
        int iyf = x + (y+1)*nxs + nxs*ny*z + nxs*ny*nz*t;
        int iyb = x + (y-1)*nxs + nxs*ny*z + nxs*ny*nz*t;

#ifdef PREFETCH
	int pfy = (i0+_PFI - pft*nt - pfz*nz*nt)/nxs;
#endif

        scs_t tmp0 __attribute__((aligned(_ALIGN_SIZE)));

        for (int c = 0; c < 3; c++) {
        for (int s = 0; s < 4; s++) {
          for (int ri = 0; ri < 2; ri++) {
            for (int j = 0; j<VLENS; j++) {
              tmp0.c[c][s][ri][j] = 0.0f;
            }
          }
        }
        }

        //
        // Y-forward
        //
#ifdef PREFETCH
	if ( pfy != ny-1) {
	  __prefetch_inp(in, (iyf+_PFI));
	  __prefetch_su3(gy, (i0 +_PFI));
	}
#endif
        if ( y != ny-1 ) {

          __attribute__((aligned(_ALIGN_SIZE))) projscs_t a,ua;

          __mult_y_forw_pre_(a,(*(in + iyf)));
          __mult_u_y_(ua,a,(*(gy + i0)));
          __mult_y_forw_pst_(tmp0,ua);

        } 

        //
        // Y-backward
        //
#ifdef PREFETCH
	if ( pfy != 0) {
	  __prefetch_inp(in, (iyb+_PFI));
	  __prefetch_su3(gy, (iyb+_PFI));
	}
#endif
        if ( y != 0 ) {

          __attribute__((aligned(_ALIGN_SIZE))) projscs_t a,ua;

          __mult_y_back_pre_(a,(*(in + iyb)));
          __mult_udag_y_(ua,a,(*(gy + iyb)));
          __mult_y_back_pst_(tmp0,ua);

        } 

	__prefetch_out(out, (i0+_PFI));
#pragma loop norecurrence
        for (int c = 0; c < 3; c++) {
          for (int s = 0; s < 4; s++) {
            for (int ri = 0; ri < 2; ri++) {
              for (int j = 0; j < VLENS; j++) {
                ((float(*)[4][2][VLENS])((out + i0)->c))[c][s][ri][j] += ((float(*)[4][2][VLENS])(tmp0.c))[c][s][ri][j] ;
              }
            }
          }
        }

      }

    }
    }
    }


        //
        // Z
        //
#pragma omp parallel for collapse(3) schedule(static)
    for(int t = 0; t < nt; ++t) {
    for(int z = 0; z < nz; ++z) {
    for(int y = 0; y < ny; ++y) {

      for(int x = 0; x < nxs; ++x) {

        int i0  =  x    + nxs*y + nxs*ny*z + nxs*ny*nz*t;
        int izf = x + nxs*y + (z+1)*nxs*ny + nxs*ny*nz*t;
        int izb = x + nxs*y + (z-1)*nxs*ny + nxs*ny*nz*t;
#ifdef PREFETCH
	int pfz = (i0+_PFI - pft*nt)/(nxs*ny);
#endif

        scs_t tmp0 __attribute__((aligned(_ALIGN_SIZE)));

        for (int c = 0; c < 3; c++) {
        for (int s = 0; s < 4; s++) {
          for (int ri = 0; ri < 2; ri++) {
            for (int j = 0; j<VLENS; j++) {
              tmp0.c[c][s][ri][j] = 0.0f;
            }
          }
        }
        }

        //
        // Z-forward
        //
#ifdef PREFETCH
	if ( pfz != nz-1) {
	  __prefetch_inp(in, (izf+_PFI));
	  __prefetch_su3(gz, (i0 +_PFI));
	}
#endif
        if ( z != nz-1 ) {

          __attribute__((aligned(_ALIGN_SIZE))) projscs_t a,ua;

          __mult_z_forw_pre_(a,(*(in + izf)));
          __mult_u_y_(ua,a,(*(gz + i0)));
          __mult_z_forw_pst_(tmp0,ua);

        }

        //
        // Z-backward
        //
#ifdef PREFETCH
	if ( pfz != 0) {
	  __prefetch_inp(in, (izb+_PFI));
	  __prefetch_su3(gz, (izb+_PFI));
	}
#endif
        if ( z != 0 ) {

          __attribute__((aligned(_ALIGN_SIZE))) projscs_t a,ua;

          __mult_z_back_pre_(a,(*(in + izb)));
          __mult_udag_y_(ua,a,(*(gz + izb)));
          __mult_z_back_pst_(tmp0,ua);

        } 

        __prefetch_out(out, (i0+_PFI));
#pragma loop norecurrence
        for (int c = 0; c < 3; c++) {
          for (int s = 0; s < 4; s++) {
            for (int ri = 0; ri < 2; ri++) {
              for (int j = 0; j < VLENS; j++) {
                ((float(*)[4][2][VLENS])((out + i0)->c))[c][s][ri][j] += ((float(*)[4][2][VLENS])(tmp0.c))[c][s][ri][j] ;
              }
            }
          }
        }

      }

    }
    }
    }


        //
        // T
        //
#pragma omp parallel for collapse(3) schedule(static)
    for(int t = 0; t < nt; ++t) {
    for(int z = 0; z < nz; ++z) {
    for(int y = 0; y < ny; ++y) {

      for(int x = 0; x < nxs; ++x) {

        int i0  =  x    + nxs*y + nxs*ny*z + nxs*ny*nz*t;
        int itf = x + nxs*y + nxs*ny*z + (t+1)*nxs*ny*nz;
        int itb = x + nxs*y + nxs*ny*z + (t-1)*nxs*ny*nz;

#ifdef PREFETCH
	int pft = (i0+_PFI)/(nxs*ny*nz);
#endif

        scs_t tmp0 __attribute__((aligned(_ALIGN_SIZE)));

        for (int c = 0; c < 3; c++) {
        for (int s = 0; s < 4; s++) {
          for (int ri = 0; ri < 2; ri++) {
            for (int j = 0; j<VLENS; j++) {
              tmp0.c[c][s][ri][j] = 0.0f;
            }
          }
        }
        }

        //
        // T-forward
        //
#ifdef PREFETCH
	if ( pft != nt-1) {
	  __prefetch_inp(in, (itf+_PFI));
	  __prefetch_su3(gt, (i0 +_PFI));
	}
#endif
        if ( t != nt-1 ) { 

          __attribute__((aligned(_ALIGN_SIZE))) projscs_t a,ua;

          __mult_t_forw_pre_(a,(*(in + itf)));
          __mult_u_y_(ua,a,(*(gt + i0)));
          __mult_t_forw_pst_(tmp0,ua);

        } 

        //
        // T-backward
        //
#ifdef PREFETCH
	if ( pft != 0) {
	  __prefetch_inp(in, (itb+_PFI));
	  __prefetch_su3(gt, (itb+_PFI));
	}
#endif
        if ( t != 0 ) { 

          __attribute__((aligned(_ALIGN_SIZE))) projscs_t a,ua;

          __mult_t_back_pre_(a,(*(in + itb)));
          __mult_udag_y_(ua,a,(*(gt + itb)));
          __mult_t_back_pst_(tmp0,ua);

        }

        __prefetch_out(out, (i0+_PFI));
#pragma loop norecurrence
        for (int c = 0; c < 3; c++) {
          for (int s = 0; s < 4; s++) {
            for (int ri = 0; ri < 2; ri++) {
              for (int j = 0; j < VLENS; j++) {
                ((float(*)[4][2][VLENS])((out + i0)->c))[c][s][ri][j] += ((float(*)[4][2][VLENS])(tmp0.c))[c][s][ri][j] ;
              }
            }
          }
        }

      }

    }
    }
    }




        //
        // CLV
        //
#pragma omp parallel for collapse(3) schedule(static)
    for(int t = 0; t < nt; ++t) {
    for(int z = 0; z < nz; ++z) {
    for(int y = 0; y < ny; ++y) {

      for(int x = 0; x < nxs; ++x) {

        int i0  =  x    + nxs*y + nxs*ny*z + nxs*ny*nz*t;

	__prefetch_clv(cl, (i0+_PFI));
	__prefetch_out(out, (i0+_PFI));
        __mult_clvs(out[i0].cv, cl[i0].cv);

#pragma loop norecurrence
        for (int c = 0; c < 3; ++c){
        for (int s = 0; s < 4; ++s){
          for (int ri = 0; ri < 2; ++ri){
            for (int j = 0; j < VLENS; ++j){
              ((float(*)[4][2][VLENS])((out + i0)->c))[c][s][ri][j] = ((float(*)[4][2][VLENS])((in + i0)->c))[c][s][ri][j] + ((float(*)[4][2][VLENS])((out + i0)->c))[c][s][ri][j] * (float)mkappa;
            }
          }
        }
        }

      }

    }
    }
    }


  }
#ifdef __cplusplus
}
#endif
