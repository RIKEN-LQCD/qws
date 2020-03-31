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
#ifndef JINV_DDD_IN_H_H
#define JINV_DDD_IN_H_H

void jinv_ddd_in_h_(sch_t * __restrict__ x, const sch_t * __restrict__ b, const int *DEO, const int *maxiter)
//
//  Multiply approximate inverse of  Wilson/Clover operator in a domain block
//
//  xe = Aee be   or   xo = Aoo bo
//
//  Aee, Aoo : approximate inverse for (Dee)^-1 and (Doo)^-1, respectively
//
//        x : quark field in a even/odd domain (output)
//        b : quark field in a even/odd domain (input)
//      DEO : even/odd block index (0 for even, 1 for odd)
//  maxiter : Jacobbi tieration count for approximate inverse
//
{
  __attribute__((aligned(64))) static sch_t *q;

  if (q==0) q = (sch_t*)malloc( sizeof(sch_t) * vols);
  //rvech_t rvd0;
  //rvd0 = fload1_s((float)2);

  ///////////////////
  //  q = Ab
  ///////////////////
  ddd_in_h_(q, b, DEO);

  ///////////////////
  // x = 2 b - q
  ///////////////////
#pragma omp parallel for
  for(int i=0; i<vols; i++){

    for(int j=0; j<24; j++){
      for(int v=0; v < VLENS; v++) {
        x[i].ccs[j].v[v] = 2.0f * b[i].ccs[j].v[v] - q[i].ccs[j].v[v];
      }
    }

  }

  //////////////////////////
  // Jacobbi iteration
  //////////////////////////
  for (int iter=1; iter<(*maxiter);iter++){

    ///////////////////
    // q = Ax
    ///////////////////
    ddd_in_h_(q, x, DEO);

    ///////////////////
    // x = x + b - q
    ///////////////////
#pragma omp parallel for
    for(int i=0; i<vols; i++){
      for(int j=0; j<24; j++){
        for(int v=0; v < VLENS; v++) {
          x[i].ccs[j].v[v] += b[i].ccs[j].v[v] - q[i].ccs[j].v[v];
        }
      }
    }

  }//iter

  // free(q);

}

#endif
