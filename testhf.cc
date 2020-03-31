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
#include <stdio.h>
#include <stdint.h>
#include <bitset>

typedef  int16_t s16;
typedef  int32_t s32;
typedef unsigned short u16;
typedef unsigned int   u32;

//#include "half_float.h"
//#define libhalf
#ifdef libhalf
#include "half.hpp"

//using namespace HalfFloat;
using namespace half_float;
#else
typedef __fp16 half;
#endif
typedef half float16;

typedef union {
  float f;
  u32   u;
} floatbit;

int main(void)
{
  float   as = 2.9452903e-5f;
  float16 ah = float16(as);

  floatbit asb; asb.f =       as;
  floatbit ahb; ahb.f = float(ah);

  std::bitset<32> as_bit(asb.u);
  std::bitset<32> ah_bit(ahb.u);
  
  printf("as = %18.7e\n",as);
  printf("as = %X\n",    asb.u);
  printf("as = %s\n",    as_bit.to_string().c_str());

  u32 nsign =  ((asb.u) >> 31);
  u32 nexp  = (((asb.u) >> 23) & 0x000000FF) - 127;
  u32 nfrac = ( (asb.u) & 0x007FFFFF);
  printf("sign %d\n",nsign);
  printf("exp  %d\n",nexp);
  printf("frac %s\n",std::bitset<32>(nfrac).to_string().c_str());

  printf("ah = %18.7e\n",float(ah));
  printf("ah = %X\n",    ahb.u);
  printf("ah = %s\n",    ah_bit.to_string().c_str());

}
