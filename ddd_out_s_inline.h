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

//
// X-forward (top end)
//
#define _mult_forw_x_recv_ {\
\
  int ix = y + ny*z + ny*nz*t;\
\
  __attribute__((aligned(_ALIGN_SIZE))) projscs1_t ua;\
\
  __mult_u_y_3_(ua,(*(xfs_recv + ix)),(*(gx + i0)));\
  __mult_x_forw_pst_3_(tmp,ua);\
\
}

//
// X-backward (bottom end)
//
#define _mult_back_x_recv_ {\
\
  int ix = y + ny*z + ny*nz*t;\
\
  projscs1_t *xbs_recvi = xbs_recv + ix;\
\
  __mult_x_back_pst_3_(tmp,(*xbs_recvi));\
\
}

//
// Y-forward  (top end)
//
#define _mult_forw_y_recv_ {\
\
  int i2 = x +nxs*z +nxs*nz*t;\
\
  __attribute__((aligned(_ALIGN_SIZE))) projscs_t ua;\
\
  __mult_u_y_(ua,(*(yfs_recv + i2)),(*(gy + i0)));\
  __mult_y_forw_pst_(tmp,ua);\
\
}

//
// Y-backward (bottom end)
//
#define _mult_back_y_recv_ {\
\
  int i3 = x +nxs*z +nxs*nz*t;\
  projscs_t *ybs_recvi = ybs_recv + i3;\
\
  __mult_y_back_pst_(tmp,(*ybs_recvi));\
\
}

//
// Z-forward (top end)
//
#define _mult_forw_z_recv_ {\
\
  int i4 = x +nxs*y +nxs*ny*t;\
\
  __attribute__((aligned(_ALIGN_SIZE))) projscs_t ua;\
\
  __mult_u_y_(ua,(*(zfs_recv + i4)),(*(gz + i0)));\
  __mult_z_forw_pst_(tmp,ua);\
\
}

//
// Z-backward (bottom end)
//
#define _mult_back_z_recv_ {\
\
  int i5 = x +nxs*y +nxs*ny*t;\
  projscs_t *zbs_recvi = zbs_recv + i5;\
\
  __mult_z_back_pst_(tmp,(*zbs_recvi));\
\
}

//
// T-forward (top end)
//
#define _mult_forw_t_recv_ {\
\
  int i6 = x +nxs*y +nxs*ny*z;\
\
  __attribute__((aligned(_ALIGN_SIZE))) projscs_t ua;\
  float tbc_fwd = ((float)fbc[3][0])*0.5f;\
\
  __mult_u_y_(ua,(*(tfs_recv + i6)),(*(gt + i0)));\
  __mult_t_forw_pst_bc_(tmp,ua,tbc_fwd);\
\
}

//
// T-backward (bottom end)
//
#define _mult_back_t_recv_ {\
\
  int i7 = x +nxs*y +nxs*ny*z;\
  projscs_t *tbs_recvi = tbs_recv + i7;\
  float tbc_bwd = ((float)fbc[3][1])*0.5f;\
\
  __mult_t_back_pst_bc_(tmp,(*tbs_recvi),tbc_bwd);\
\
}

//
// mult clover and accumulate 
//
#define _mult_clover_and_accumulate_recv_ {\
\
  __mult_clvs( tmp.cv, clvs[i0 + vols*(*idomain)].cv);\
\
  scs_t *outi = out + i0;\
  for (int c = 0; c < 3; c++) {\
  for (int s = 0; s < 4; s++) {\
    for (int j = 0; j < VLENS; j++) {\
      outi->c[c][s][0][j] += tmp.c[c][s][0][j] * factor;\
      outi->c[c][s][1][j] += tmp.c[c][s][1][j] * factor;\
    }\
  }\
  }\
}
