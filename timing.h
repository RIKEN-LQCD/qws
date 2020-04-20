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
#ifndef TIMI_H
#define TIMI_H

#include"profiler.h"

#ifndef PROF_SELECTIVE
#define IF_PROF(A) A
#else
#define IF_PROF(A) if(prof_flag){A}
#ifdef MAIN
int prof_flag;
#else
extern int prof_flag;
#endif
#endif

#ifdef PROF_TARGET
# define _TIC_
# define _TOC_ 
# define _BCG_ITER_TIC_
# define _BCG_ITER_TOC_ 
# define _BCG_DDS_ITER_TIC_
# define _BCG_DDS_ITER_TOC_ 
# define _BCG_PRECDDS_TIC_
# define _BCG_PRECDDS_TOC_ 
# define _BCG_PRECDDS_ITER_TIC_
# define _BCG_PRECDDS_ITER_TOC_ 
# define _BCG_PRECDDS_ITER_REDUC1_TIC_
# define _BCG_PRECDDS_ITER_REDUC1_TOC_ 
# define _BCG_PRECDDS_ITER_REDUC2_TIC_
# define _BCG_PRECDDS_ITER_REDUC2_TOC_ 
# define _BCG_PRECDDS_ITER_REDUC3_TIC_
# define _BCG_PRECDDS_ITER_REDUC3_TOC_ 
# define _BCG_PRECDDS_ITER_BARRIER_BEFORE_REDUC1_TIC_
# define _BCG_PRECDDS_ITER_BARRIER_BEFORE_REDUC1_TOC_ 
# define _BCG_PRECDDS_ITER_BARRIER_BEFORE_REDUC2_TIC_
# define _BCG_PRECDDS_ITER_BARRIER_BEFORE_REDUC2_TOC_ 
# define _BCG_PRECDDS_ITER_BARRIER_BEFORE_REDUC3_1_TIC_
# define _BCG_PRECDDS_ITER_BARRIER_BEFORE_REDUC3_1_TOC_ 
# define _BCG_PRECDDS_ITER_BARRIER_BEFORE_REDUC3_2_TIC_
# define _BCG_PRECDDS_ITER_BARRIER_BEFORE_REDUC3_2_TOC_ 
# define _MCG_ITER_TIC_
# define _MCG_ITER_TOC_ 
# define _CG_ITER_TIC_
# define _CG_ITER_TOC_ 
# define _DEO_IN_TIC_
# define _DEO_IN_TOC_ 
# define _DEE_DEO_IN_TIC_
# define _DEE_DEO_IN_TOC_ 
# define _DEO_OUT_PRE_TIC_
# define _DEO_OUT_PRE_TOC_ 
# define _DEE_DEO_OUT_POST_TIC_
# define _DEE_DEO_OUT_POST_TOC_ 
# define _MTILDE_TIC_
# define _MTILDE_TOC_ 
# define _PREC_S_TIC_
# define _PREC_S_TOC_ 
# define _PREC_DDD_S_TIC_
# define _PREC_DDD_S_TOC_ 
# define _JINV_DDD_IN_S_TIC_
# define _JINV_DDD_IN_S_TOC_
# define _DDD_IN_S_TIC_
# define _DDD_IN_S_TOC_ 
# define _DDD_OUT_PRE_S_TIC_
# define _DDD_OUT_PRE_S_TOC_ 
# define _DDD_OUT_POS_S_TIC_
# define _DDD_OUT_POS_S_TOC_
# define _S_MULT_WD_DEO_OUT_RECV_HPC_CALC_TIC_
# define _S_MULT_WD_DEO_OUT_RECV_HPC_CALC_TOC_
# define _S_MULT_WD_DEO_OUT_SEND_HPC_TIC_
# define _S_MULT_WD_DEO_OUT_SEND_HPC_TOC_
# define _COMM_TIC_
# define _COMM_TOC_
# define _OTHER_CALC_TIC_
# define _OTHER_CALC_TOC_
# define _OVERLAPPED_CALC_TIC_
# define _OVERLAPPED_CALC_TOC_
# define _COMLIB_ISEND_ALL_C_TIC_
# define _COMLIB_ISEND_ALL_C_TOC_
# define _COMLIB_SEND_WAIT_ALL_C_TIC_
# define _COMLIB_SEND_WAIT_ALL_C_TOC_

# define TARGET_ALL        1
# define TARGET_JINV       2
# define TARGET_IN         3
# define TARGET_PRE        4
# define TARGET_POS        5
# define TARGET_OTHER      6
# define TARGET_ALL_CALC   7
# define TARGET_OVERLAPPED 8
# define TARGET_SEND       9
# define TARGET_SEND_POST 10
# define TARGET_RECV      11
# define TARGET_REDUC1    12
# define TARGET_REDUC2    13
# define TARGET_REDUC3    14

# if PROF_TARGET == TARGET_ALL
#  undef _BCG_PRECDDS_ITER_TIC_
#  undef _BCG_PRECDDS_ITER_TOC_
#  define _BCG_PRECDDS_ITER_TIC_                IF_PROF(PROF_START_SRL("all"))
#  define _BCG_PRECDDS_ITER_TOC_                IF_PROF(PROF_STOP_SRL("all"))
# elif PROF_TARGET == TARGET_JINV
#  undef _JINV_DDD_IN_S_TIC_
#  undef _JINV_DDD_IN_S_TOC_
#  define _JINV_DDD_IN_S_TIC_                   IF_PROF(PROF_START_SRL("jinv"))
#  define _JINV_DDD_IN_S_TOC_                   IF_PROF(PROF_STOP_SRL("jinv"))
# elif PROF_TARGET == TARGET_IN
#  undef _DDD_IN_S_TIC_
#  undef _DDD_IN_S_TOC_
#  define _DDD_IN_S_TIC_                        IF_PROF(PROF_START_SRL("in"))
#  define _DDD_IN_S_TOC_                        IF_PROF(PROF_STOP_SRL("in"))
# elif PROF_TARGET == TARGET_PRE
#  undef _S_MULT_WD_DEO_OUT_SEND_HPC_TIC_
#  undef _S_MULT_WD_DEO_OUT_SEND_HPC_TOC_
#  define _S_MULT_WD_DEO_OUT_SEND_HPC_TIC_      IF_PROF(PROF_START_SRL("pre"))
#  define _S_MULT_WD_DEO_OUT_SEND_HPC_TOC_      IF_PROF(PROF_STOP_SRL("pre"))
# elif PROF_TARGET == TARGET_POS
#  undef _S_MULT_WD_DEO_OUT_RECV_HPC_CALC_TIC_
#  undef _S_MULT_WD_DEO_OUT_RECV_HPC_CALC_TOC_
#  define _S_MULT_WD_DEO_OUT_RECV_HPC_CALC_TIC_ IF_PROF(PROF_START_SRL("pos"))
#  define _S_MULT_WD_DEO_OUT_RECV_HPC_CALC_TOC_ IF_PROF(PROF_STOP_SRL("pos"))
# elif PROF_TARGET == TARGET_OTHER
#  undef _BCG_PRECDDS_ITER_TIC_
#  undef _BCG_PRECDDS_ITER_TOC_
#  undef _BCG_PRECDDS_ITER_REDUC1_TIC_
#  undef _BCG_PRECDDS_ITER_REDUC1_TOC_
#  undef _BCG_PRECDDS_ITER_REDUC2_TIC_
#  undef _BCG_PRECDDS_ITER_REDUC2_TOC_
#  undef _BCG_PRECDDS_ITER_REDUC3_TIC_
#  undef _BCG_PRECDDS_ITER_REDUC3_TOC_
#  undef _OTHER_CALC_TIC_
#  undef _OTHER_CALC_TOC_
#  define _BCG_PRECDDS_ITER_TIC_                IF_PROF(PROF_START("other"))
#  define _BCG_PRECDDS_ITER_TOC_                IF_PROF(PROF_STOP("other"))
#  define _BCG_PRECDDS_ITER_REDUC1_TIC_         IF_PROF(PROF_STOP("other"))
#  define _BCG_PRECDDS_ITER_REDUC1_TOC_         IF_PROF(PROF_START("other"))
#  define _BCG_PRECDDS_ITER_REDUC2_TIC_         IF_PROF(PROF_STOP("other"))
#  define _BCG_PRECDDS_ITER_REDUC2_TOC_         IF_PROF(PROF_START("other"))
#  define _BCG_PRECDDS_ITER_REDUC3_TIC_         IF_PROF(PROF_STOP ("other"))
#  define _BCG_PRECDDS_ITER_REDUC3_TOC_         IF_PROF(PROF_START("other"))
#  define _OTHER_CALC_TIC_                      IF_PROF(PROF_START_SRL("other"))
#  define _OTHER_CALC_TOC_                      IF_PROF(PROF_STOP_SRL("other"))
# elif PROF_TARGET == TARGET_ALL_CALC
#  undef _BCG_PRECDDS_ITER_TIC_
#  undef _BCG_PRECDDS_ITER_TOC_
#  undef _BCG_PRECDDS_ITER_REDUC1_TIC_
#  undef _BCG_PRECDDS_ITER_REDUC1_TOC_
#  undef _BCG_PRECDDS_ITER_REDUC2_TIC_
#  undef _BCG_PRECDDS_ITER_REDUC2_TOC_
#  undef _BCG_PRECDDS_ITER_REDUC3_TIC_
#  undef _BCG_PRECDDS_ITER_REDUC3_TOC_
#  undef _S_MULT_WD_DEO_OUT_SEND_HPC_TOC_
#  undef _S_MULT_WD_DEO_OUT_RECV_HPC_CALC_TIC_
#  undef _COMLIB_SEND_WAIT_ALL_C_TIC_
#  undef _COMLIB_SEND_WAIT_ALL_C_TOC_
#  define _BCG_PRECDDS_ITER_TIC_                IF_PROF(PROF_START("all_calc"))
#  define _BCG_PRECDDS_ITER_TOC_                IF_PROF(PROF_STOP("all_calc"))
#  define _BCG_PRECDDS_ITER_REDUC1_TIC_         IF_PROF(PROF_STOP("all_calc"))
#  define _BCG_PRECDDS_ITER_REDUC1_TOC_         IF_PROF(PROF_START("all_calc"))
#  define _BCG_PRECDDS_ITER_REDUC2_TIC_         IF_PROF(PROF_STOP("all_calc"))
#  define _BCG_PRECDDS_ITER_REDUC2_TOC_         IF_PROF(PROF_START("all_calc"))
#  define _BCG_PRECDDS_ITER_REDUC3_TIC_         IF_PROF(PROF_STOP ("all_calc"))
#  define _BCG_PRECDDS_ITER_REDUC3_TOC_         IF_PROF(PROF_START("all_calc"))

// The regions from _S_MULT_WD_DEO_OUT_SEND_HPC_TOC_ to _S_MULT_WD_DEO_OUT_RECV_HPC_CALC_TIC_
// and from _COMLIB_SEND_WAIT_ALL_C_TIC_ to _COMLIB_SEND_WAIT_ALL_C_TOC_
// include communications and caluculations which overlap with the communications,
// so it is not counted as calculation time.
#  define _S_MULT_WD_DEO_OUT_SEND_HPC_TOC_      IF_PROF(PROF_STOP_SRL("all_calc"))
#  define _S_MULT_WD_DEO_OUT_RECV_HPC_CALC_TIC_ IF_PROF(PROF_START_SRL("all_calc"))

#  define _COMLIB_SEND_WAIT_ALL_C_TIC_          IF_PROF(PROF_STOP_SRL("all_calc"))
#  define _COMLIB_SEND_WAIT_ALL_C_TOC_          IF_PROF(PROF_START_SRL("all_calc"))

# elif PROF_TARGET == TARGET_OVERLAPPED
#  undef _OVERLAPPED_CALC_TIC_
#  undef _OVERLAPPED_CALC_TOC_
#  define _OVERLAPPED_CALC_TIC_ IF_PROF(PROF_START_SRL("overlapped"))
#  define _OVERLAPPED_CALC_TOC_ IF_PROF(PROF_STOP_SRL("overlapped"))
# elif PROF_TARGET == TARGET_RECV
#  undef _COMLIB_IRECV_ALL_C_TIC_
#  undef _COMLIB_IRECV_ALL_C_TOC_
#  define _COMLIB_IRECV_ALL_C_TIC_ IF_PROF(PROF_START_SRL("comlib_irecv_all_c"))
#  define _COMLIB_IRECV_ALL_C_TOC_ IF_PROF(PROF_STOP_SRL("comlib_irecv_all_c"))
# elif PROF_TARGET == TARGET_SEND
#  undef _COMLIB_ISEND_ALL_C_TIC_
#  undef _COMLIB_ISEND_ALL_C_TOC_
#  define _COMLIB_ISEND_ALL_C_TIC_ IF_PROF(PROF_START_SRL("comlib_isend_all_c"))
#  define _COMLIB_ISEND_ALL_C_TOC_ IF_PROF(PROF_STOP_SRL("comlib_isend_all_c"))
# elif PROF_TARGET == TARGET_SEND_POST
#  undef _COMLIB_SEND_WAIT_ALL_C_TIC_
#  undef _COMLIB_SEND_WAIT_ALL_C_TOC_
#  define _COMLIB_SEND_WAIT_ALL_C_TIC_ IF_PROF(PROF_START_SRL("comlib_isend_wait_all_c"))
#  define _COMLIB_SEND_WAIT_ALL_C_TOC_ IF_PROF(PROF_STOP_SRL("comlib_isend_wait_all_c"))
# elif PROF_TARGET == TARGET_REDUC1
#  undef _BCG_PRECDDS_ITER_REDUC1_TIC_
#  undef _BCG_PRECDDS_ITER_REDUC1_TOC_
#  define _BCG_PRECDDS_ITER_REDUC1_TIC_ IF_PROF(PROF_START_SRL("bicgstab_precdd_s_iter_reduc1_"))
#  define _BCG_PRECDDS_ITER_REDUC1_TOC_ IF_PROF(PROF_STOP_SRL("bicgstab_precdd_s_iter_reduc1_"))
# elif PROF_TARGET == TARGET_REDUC2
#  undef _BCG_PRECDDS_ITER_REDUC2_TIC_
#  undef _BCG_PRECDDS_ITER_REDUC2_TOC_
#  define _BCG_PRECDDS_ITER_REDUC2_TIC_ IF_PROF(PROF_START_SRL("bicgstab_precdd_s_iter_reduc2_"))
#  define _BCG_PRECDDS_ITER_REDUC2_TOC_ IF_PROF(PROF_STOP_SRL("bicgstab_precdd_s_iter_reduc2_"))
# elif PROF_TARGET == TARGET_REDUC3
#  undef _BCG_PRECDDS_ITER_REDUC3_TIC_
#  undef _BCG_PRECDDS_ITER_REDUC3_TOC_
#  define _BCG_PRECDDS_ITER_REDUC3_TIC_ IF_PROF(PROF_START_SRL("bicgstab_precdd_s_iter_reduc3_"))
#  define _BCG_PRECDDS_ITER_REDUC3_TOC_ IF_PROF(PROF_STOP_SRL("bicgstab_precdd_s_iter_reduc3_"))
# else

#  error "invalid PROF_TARGET"
# endif

#else

# define _TIC_                                  IF_PROF(PROF_START(__func__))
# define _TOC_                                  IF_PROF( PROF_STOP(__func__))
# define _BCG_ITER_TIC_                         IF_PROF(PROF_START("bicgstab_iter_"))
# define _BCG_ITER_TOC_                         IF_PROF( PROF_STOP("bicgstab_iter_"))
# define _BCG_DDS_ITER_TIC_                     IF_PROF(PROF_START("bicgstab_dd_s_iter_"))
# define _BCG_DDS_ITER_TOC_                     IF_PROF( PROF_STOP("bicgstab_dd_s_iter_"))
# define _BCG_PRECDDS_TIC_                      IF_PROF(PROF_START("bicgstab_precdd_s_"))
# define _BCG_PRECDDS_TOC_                      IF_PROF( PROF_STOP("bicgstab_precdd_s_"))
# define _BCG_PRECDDS_ITER_TIC_                 IF_PROF(PROF_START("bicgstab_precdd_s_iter_"))
# define _BCG_PRECDDS_ITER_TOC_                 IF_PROF( PROF_STOP("bicgstab_precdd_s_iter_"))
# define _BCG_PRECDDS_ITER_REDUC1_TIC_          IF_PROF(PROF_START("bicgstab_precdd_s_iter_reduc1_"))
# define _BCG_PRECDDS_ITER_REDUC1_TOC_          IF_PROF( PROF_STOP("bicgstab_precdd_s_iter_reduc1_"))
# define _BCG_PRECDDS_ITER_REDUC2_TIC_          IF_PROF(PROF_START("bicgstab_precdd_s_iter_reduc2_"))
# define _BCG_PRECDDS_ITER_REDUC2_TOC_          IF_PROF( PROF_STOP("bicgstab_precdd_s_iter_reduc2_"))
# define _BCG_PRECDDS_ITER_REDUC3_TIC_          IF_PROF(PROF_START("bicgstab_precdd_s_iter_reduc3_"))
# define _BCG_PRECDDS_ITER_REDUC3_TOC_          IF_PROF( PROF_STOP("bicgstab_precdd_s_iter_reduc3_"))
# define _BCG_PRECDDS_ITER_BARRIER_BEFORE_REDUC1_TIC_   IF_PROF(PROF_START("bicgstab_precdd_s_iter_bar_reduc1_"))
# define _BCG_PRECDDS_ITER_BARRIER_BEFORE_REDUC1_TOC_   IF_PROF( PROF_STOP("bicgstab_precdd_s_iter_bar_reduc1_"))
# define _BCG_PRECDDS_ITER_BARRIER_BEFORE_REDUC2_TIC_   IF_PROF(PROF_START("bicgstab_precdd_s_iter_bar_reduc2_"))
# define _BCG_PRECDDS_ITER_BARRIER_BEFORE_REDUC2_TOC_   IF_PROF( PROF_STOP("bicgstab_precdd_s_iter_bar_reduc2_"))
# define _BCG_PRECDDS_ITER_BARRIER_BEFORE_REDUC3_1_TIC_ IF_PROF(PROF_START("bicgstab_precdd_s_iter_bar_reduc3_1_"))
# define _BCG_PRECDDS_ITER_BARRIER_BEFORE_REDUC3_1_TOC_ IF_PROF( PROF_STOP("bicgstab_precdd_s_iter_bar_reduc3_1_"))
# define _BCG_PRECDDS_ITER_BARRIER_BEFORE_REDUC3_2_TIC_ IF_PROF(PROF_START("bicgstab_precdd_s_iter_bar_reduc3_2_"))
# define _BCG_PRECDDS_ITER_BARRIER_BEFORE_REDUC3_2_TOC_ IF_PROF( PROF_STOP("bicgstab_precdd_s_iter_bar_reduc3_2_"))
# define _MCG_ITER_TIC_                         IF_PROF(PROF_START("mcg_iter_"))
# define _MCG_ITER_TOC_                         IF_PROF( PROF_STOP("mcg_iter_"))
# define _CG_ITER_TIC_                          IF_PROF(PROF_START("cg_iter_"))
# define _CG_ITER_TOC_                          IF_PROF( PROF_STOP("cg_iter_"))
# define _DEO_IN_TIC_                           IF_PROF(PROF_START("deo_in_"))
# define _DEO_IN_TOC_                           IF_PROF( PROF_STOP("deo_in_"))
# define _DEE_DEO_IN_TIC_                       IF_PROF(PROF_START("dee_deo_in_"))
# define _DEE_DEO_IN_TOC_                       IF_PROF( PROF_STOP("dee_deo_in_"))
# define _DEO_OUT_PRE_TIC_                      IF_PROF(PROF_START("deo_out_pre_"))
# define _DEO_OUT_PRE_TOC_                      IF_PROF( PROF_STOP("deo_out_pre_"))
# define _DEE_DEO_OUT_POST_TIC_                 IF_PROF(PROF_START("dee_deo_out_post_"))
# define _DEE_DEO_OUT_POST_TOC_                 IF_PROF( PROF_STOP("dee_deo_out_post_"))
# define _MTILDE_TIC_                           IF_PROF(PROF_START("mtilde"))
# define _MTILDE_TOC_                           IF_PROF( PROF_STOP("mtilde"))
# define _PREC_S_TIC_                           IF_PROF(PROF_START("prec_s_"))
# define _PREC_S_TOC_                           IF_PROF( PROF_STOP("prec_s_"))
# define _PREC_DDD_S_TIC_                       IF_PROF(PROF_START("prec_ddd_s_"))
# define _PREC_DDD_S_TOC_                       IF_PROF( PROF_STOP("prec_ddd_s_"))
# define _DDD_IN_S_TIC_                         IF_PROF(PROF_START("ddd_in_s_"))
# define _DDD_IN_S_TOC_                         IF_PROF( PROF_STOP("ddd_in_s_"))
# define _JINV_DDD_IN_S_TIC_                    IF_PROF(PROF_START_SRL("jinv_ddd_in_s_"))
# define _JINV_DDD_IN_S_TOC_                    IF_PROF( PROF_STOP_SRL("jinv_ddd_in_s_"))
# define _DDD_OUT_PRE_S_TIC_                    IF_PROF(PROF_START("ddd_out_pre_s_"))
# define _DDD_OUT_PRE_S_TOC_                    IF_PROF( PROF_STOP("ddd_out_pre_s_"))
# define _DDD_OUT_POS_S_TIC_                    IF_PROF(PROF_START("ddd_out_pos_s_"))
# define _DDD_OUT_POS_S_TOC_                    IF_PROF( PROF_STOP("ddd_out_pos_s_"))
# define _S_MULT_WD_DEO_OUT_RECV_HPC_CALC_TIC_  IF_PROF(PROF_START_SRL("s_mult_wd_deo_out_recv_hpc_calc_"))
# define _S_MULT_WD_DEO_OUT_RECV_HPC_CALC_TOC_  IF_PROF( PROF_STOP_SRL("s_mult_wd_deo_out_recv_hpc_calc_"))
# define _S_MULT_WD_DEO_OUT_SEND_HPC_TIC_       IF_PROF(PROF_START_SRL("s_mult_wd_deo_out_send_hpc_"))
# define _S_MULT_WD_DEO_OUT_SEND_HPC_TOC_       IF_PROF( PROF_STOP_SRL("s_mult_wd_deo_out_send_hpc_"))
# define _OTHER_CALC_TIC_
# define _OTHER_CALC_TOC_
# define _COMM_TIC_
# define _COMM_TOC_
# define _OVERLAPPED_CALC_TIC_
# define _OVERLAPPED_CALC_TOC_
# define _COMLIB_ISEND_ALL_C_TIC_
# define _COMLIB_ISEND_ALL_C_TOC_
# define _COMLIB_SEND_WAIT_ALL_C_TIC_
# define _COMLIB_SEND_WAIT_ALL_C_TOC_
#endif

#endif
