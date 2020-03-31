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
#include <math.h>

#ifdef __cplusplus
extern "C"{
#endif

  extern int nt, nz, ny, nx, nxh, nxd, vold, nxs, vols;
  
  void fermi_reorder_s2d_dd_(scd_t *out, const scs_t *in){
    int x, y, z, t, s, c, ri;
#pragma omp for private(t,z,y,x,s,c,ri) collapse(4) schedule(static)
    for(t=0; t<nt; t++){
      for(z=0; z<nz; z++){
	for(y=0; y<ny; y++){
	  for(x=0; x<nx; x++){
	    int id = (x%nxh)/VLEND + nxd*y + nxd*ny*z + nxd*ny*nz*t + vold*(x/nxh);
	    int jd = (x%nxh)%VLEND;
	    int is = (x%nxh)/VLENS + nxs*y + nxs*ny*z + nxs*ny*nz*t + vols*(x/nxh);
	    int js = (x%nxh)%VLENS;
	    for(c=0; c<NCOL; c++){
	      for(s=0; s<NDIM; s++){
		for (ri=0; ri<2; ri++){
		  out[id].c[c][s][ri][jd] = (double)in[is].c[c][s][ri][js];
		}
	      }
	    }
	  }
	}
      }
    }
  }

  void fermi_reorder_d2s_dd_(scs_t *out, const scd_t *in){
    int x, y, z, t, s, c, ri;
#pragma omp for private(t,z,y,x,s,c,ri) collapse(4) schedule(static)
    for(t=0; t<nt; t++){
      for(z=0; z<nz; z++){
	for(y=0; y<ny; y++){
	  for(x=0; x<nx; x++){
	    int id = (x%nxh)/VLEND + nxd*y + nxd*ny*z + nxd*ny*nz*t + vold*(x/nxh);
	    int jd = (x%nxh)%VLEND;
	    int is = (x%nxh)/VLENS + nxs*y + nxs*ny*z + nxs*ny*nz*t + vols*(x/nxh);
	    int js = (x%nxh)%VLENS;
	    for(c=0; c<NCOL; c++){
	      for(s=0; s<NDIM; s++){
		for (ri=0; ri<2; ri++){
		  out[is].c[c][s][ri][js] = (float)in[id].c[c][s][ri][jd];
		}
	      }
	    }
	  }
	}
      }
    }
  }

  
  void fermi_reorder_s2d_dd_renorm_(scd_t *out, const scs_t *in, const double *dnorm){
  //
  // convert to double from single with nrom
  //
  // out_DP = in_SP*dnorm
  //
  //   out : quark field in double precision, normalization recovred (output)
  //    in : quark field in single precision (input)
  // dnorm :  norm in double precision (input)
  //
    int x, y, z, t, s, c, ri;
#pragma omp for private(t,z,y,x,s,c,ri) collapse(4) schedule(static)
    for(t=0; t<nt; t++){
      for(z=0; z<nz; z++){
	for(y=0; y<ny; y++){
	  for(x=0; x<nx; x++){
	    int id = (x%nxh)/VLEND + nxd*y + nxd*ny*z + nxd*ny*nz*t + vold*(x/nxh);
	    int jd = (x%nxh)%VLEND;
	    int is = (x%nxh)/VLENS + nxs*y + nxs*ny*z + nxs*ny*nz*t + vols*(x/nxh);
	    int js = (x%nxh)%VLENS;
	    for(c=0; c<NCOL; c++){
	      for(s=0; s<NDIM; s++){
		for (ri=0; ri<2; ri++){
		  out[id].c[c][s][ri][jd] = ((double)in[is].c[c][s][ri][js])*(*dnorm);
		}
	      }
	    }
	  }
	}
      }
    }
  }

  void fermi_reorder_d2s_dd_renorm_(scs_t *out, const scd_t *in, double *dnorm){
  //
  // convert to single from double with normalization
  //
  // out_SP = in_DP/|in_DP|
  // dnorm = |in_DP|
  //
  //   out : quark field in single precision, normazlised (output)
  //    in : quark field in double precision (input)
  // dnorm :  norm in double precision (output)
  //
    int x, y, z, t, s, c, ri;
    double dtmp;

    dtmp = 0.0e0;
#pragma omp parallel for private(t,z,y,x,s,c,ri) collapse(4) reduction(+:dtmp) schedule(static) 
    for(t=0; t<nt; t++){
      for(z=0; z<nz; z++){
	for(y=0; y<ny; y++){
	  for(x=0; x<nx; x++){
	    int id = (x%nxh)/VLEND + nxd*y + nxd*ny*z + nxd*ny*nz*t + vold*(x/nxh);
	    int jd = (x%nxh)%VLEND;

	    for(c=0; c<NCOL; c++){
	      for(s=0; s<NDIM; s++){
		for (ri=0; ri<2; ri++){
                  dtmp += in[id].c[c][s][ri][jd]*in[id].c[c][s][ri][jd];
		}
	      }
	    }

	  }
	}
      }
    }
#ifdef _MPI_
    MPI_Allreduce(MPI_IN_PLACE,(void *)&dtmp,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
    dtmp = sqrt(dtmp);
    *dnorm = dtmp;
    dtmp = 1.0e0/dtmp;

#pragma omp parallel for private(t,z,y,x,s,c,ri) collapse(4) schedule(static)
    for(t=0; t<nt; t++){
      for(z=0; z<nz; z++){
	for(y=0; y<ny; y++){
	  for(x=0; x<nx; x++){
	    int id = (x%nxh)/VLEND + nxd*y + nxd*ny*z + nxd*ny*nz*t + vold*(x/nxh);
	    int jd = (x%nxh)%VLEND;
	    int is = (x%nxh)/VLENS + nxs*y + nxs*ny*z + nxs*ny*nz*t + vols*(x/nxh);
	    int js = (x%nxh)%VLENS;
	    for(c=0; c<NCOL; c++){
	      for(s=0; s<NDIM; s++){
		for (ri=0; ri<2; ri++){
		  out[is].c[c][s][ri][js] = (float)(in[id].c[c][s][ri][jd]*dtmp);
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
