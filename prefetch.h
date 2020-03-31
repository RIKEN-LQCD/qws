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
#pragma once

#if VLENS == 16 && !defined(DISABLE_PREFETCH)
// The cache line size is assumed to be 256 bytes.
// a is assumed to be aligned on 256-byte boundary. Then the last elements of su3, out and clv are
// guaranteed to be included in the last prefetched line.
#define __prefetch_su3(a,i){                                            \
    __builtin_prefetch((void*)&((*(a+i)).c_prefetch[64*0]), 0, 3);      \
    __builtin_prefetch((void*)&((*(a+i)).c_prefetch[64*1]), 0, 3);      \
    __builtin_prefetch((void*)&((*(a+i)).c_prefetch[64*2]), 0, 3);      \
    __builtin_prefetch((void*)&((*(a+i)).c_prefetch[64*3]), 0, 3);      \
    __builtin_prefetch((void*)&((*(a+i)).c_prefetch[64*4]), 0, 3);      \
  }
#define __prefetch_inp(a,i){                                    \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*0], 0, 3);       \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*1], 0, 3);       \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*2], 0, 3);       \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*3], 0, 3);       \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*4], 0, 3);       \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*5], 0, 3);       \
  }
#define __prefetch_out(a,i){                                    \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*0], 1, 3);       \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*1], 1, 3);       \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*2], 1, 3);       \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*3], 1, 3);       \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*4], 1, 3);       \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*5], 1, 3);       \
  }
#define __prefetch_clv(a,i){                                    \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*0], 0, 3);       \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*1], 0, 3);       \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*2], 0, 3);       \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*3], 0, 3);       \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*4], 0, 3);       \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*5], 0, 3);       \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*6], 0, 3);       \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*7], 0, 3);       \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*8], 0, 3);       \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*9], 0, 3);       \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*10], 0, 3);      \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*11], 0, 3);      \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*12], 0, 3);      \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*13], 0, 3);      \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*14], 0, 3);      \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*15], 0, 3);      \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*16], 0, 3);      \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*17], 0, 3);      \
  }
#define __prefetch_sendbuf(a,i){                                        \
    __builtin_prefetch((void*)&((*(a+i)).c_prefetch[64*0]), 1, 3);      \
    __builtin_prefetch((void*)&((*(a+i)).c_prefetch[64*1]), 1, 3);      \
    __builtin_prefetch((void*)&((*(a+i)).c_prefetch[64*2]), 1, 3);      \
  }
#define __prefetch_recvbuf(a,i){                                        \
    __builtin_prefetch((void*)&((*(a+i)).c_prefetch[64*0]), 0, 3);      \
    __builtin_prefetch((void*)&((*(a+i)).c_prefetch[64*1]), 0, 3);      \
    __builtin_prefetch((void*)&((*(a+i)).c_prefetch[64*2]), 0, 3);      \
  }
#elif VLENS == 8 && !defined(DISABLE_PREFETCH)
#define __prefetch_su3(a,i){                                            \
    __builtin_prefetch((void*)&((*(a+i)).c_prefetch[64*0]), 0, 3);      \
    __builtin_prefetch((void*)&((*(a+i)).c_prefetch[64*1]), 0, 3);      \
    __builtin_prefetch((void*)&((*(a+i)).c_prefetch[64*2]), 0, 3);      \
  }
#define __prefetch_inp(a,i){                                    \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*0], 0, 3);       \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*1], 0, 3);       \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*2], 0, 3);       \
  }
#define __prefetch_out(a,i){                                    \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*0], 1, 3);       \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*1], 1, 3);       \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*2], 1, 3);       \
  }
#define __prefetch_clv(a,i){                                    \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*0], 0, 3);       \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*1], 0, 3);       \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*2], 0, 3);       \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*3], 0, 3);       \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*4], 0, 3);       \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*5], 0, 3);       \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*6], 0, 3);       \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*7], 0, 3);       \
    __builtin_prefetch(&(*(a+i)).c_prefetch[64*8], 0, 3);       \
  }
#define __prefetch_sendbuf(a,i){                                        \
    __builtin_prefetch((void*)&((*(a+i)).c_prefetch[64*0]), 1, 3);      \
    __builtin_prefetch((void*)&((*(a+i)).c_prefetch[64*1]), 1, 3);      \
    __builtin_prefetch((void*)&((*(a+i)).c_prefetch[64*2]), 1, 3);      \
  }
#define __prefetch_recvbuf(a,i){                                        \
    __builtin_prefetch((void*)&((*(a+i)).c_prefetch[64*0]), 0, 3);      \
    __builtin_prefetch((void*)&((*(a+i)).c_prefetch[64*1]), 0, 3);      \
    __builtin_prefetch((void*)&((*(a+i)).c_prefetch[64*2]), 0, 3);      \
  }
#else
#define __prefetch_su3(a, i)
#define __prefetch_inp(a, i)
#define __prefetch_out(a, i)
#define __prefetch_clv(a, i)
#define __prefetch_sendbuf(a,i)
#define __prefetch_recvbuf(a,i)

#endif
