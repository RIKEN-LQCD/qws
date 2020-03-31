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
#ifndef PROFILER_H
#define PROFILER_H

#ifdef _CHECK_PA
#include"readpic.h"
#define PROF_INIT         set_pic_env(); reset_pic_prl();
#define PROF_FINALIZE     add_pic_i_prl(0); output_prof(1);

// to be used in outside of parallel regions. *_prl() functions are parallelized
#define PROF_START(a)     start_pic_sysusr_prl();
#define PROF_STOP(a)      stop_pic_prl();

// to be used in inside of parallel regions
#define PROF_START_SRL(a) start_pic_sysusr();
#define PROF_STOP_SRL(a)  stop_pic();

#define PROF_START_ALL
#define PROF_STOP_ALL

#elif defined(_FPCOLL)
#include <fjcoll.h>
#define PROF_INIT
#define PROF_FINALIZE
#define PROF_START(a)     start_collection(a);
#define PROF_STOP(a)      stop_collection(a);
#define PROF_START_SRL(a) start_collection(a);
#define PROF_STOP_SRL(a)  stop_collection(a);
#define PROF_START_ALL
#define PROF_STOP_ALL

#elif defined(_FAPP)
#include <fj_tool/fapp.h>
#define PROF_INIT
#define PROF_FINALIZE
#define PROF_START(a)     fapp_start(a,1,0);
#define PROF_STOP(a)      fapp_stop(a,1,0);
#define PROF_START_SRL(a) fapp_start(a,1,0);
#define PROF_STOP_SRL(a)  fapp_stop(a,1,0);
#define PROF_START_ALL
#define PROF_STOP_ALL

#elif defined(_CHECK_TIMING)
#ifdef __cplusplus
extern "C" {
#endif
#ifdef __cplusplus
  void check_timing_ (const char *);
}
#endif
#define PROF_INIT
#define PROF_FINALIZE
#define PROF_START(a)     check_timing_(a);
#define PROF_STOP(a)      check_timing_(a);
#define PROF_START_SRL(a) check_timing_(a);
#define PROF_STOP_SRL(a)  check_timing_(a);
#define PROF_START_ALL
#define PROF_STOP_ALL

#elif defined(_CHECK_TIMING2)

#ifdef __cplusplus
extern "C" {
#endif
void prof_time_start(void);
void prof_time_stop(void);
long int prof_time_get_usec(void);
long int prof_time_get_count(void);
void prof_time_init(void);
void prof_time_fin(void);
void prof_time_print(void);
void prof_time_print_prec(int);

void prof_time_start_(void);
void prof_time_stop_(void);
void prof_time_init_(void);
void prof_time_fin_(void);
void prof_time_print_(void);
void prof_time_get_usec_(long int*);
void prof_time_get_count_(long int*);
void prof_time_print_prec_(int*);
#ifdef __cplusplus
}
#endif

#define PROF_INIT         prof_time_init(); prof_time_print_prec(10000);
#define PROF_FINALIZE     prof_time_print(); prof_time_fin();
#define PROF_START_ALL       
#define PROF_STOP_ALL        
#define PROF_START(a)     prof_time_start();
#define PROF_STOP(a)      prof_time_stop();
#define PROF_START_SRL(a) prof_time_start();
#define PROF_STOP_SRL(a)  prof_time_stop();

#else

#define PROF_INIT
#define PROF_FINALIZE
#define PROF_START(a)
#define PROF_STOP(a)
#define PROF_START_SRL(a)
#define PROF_STOP_SRL(a)
#define PROF_START_ALL
#define PROF_STOP_ALL
#endif

#endif
