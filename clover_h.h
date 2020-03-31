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
#ifndef _CLOVER_H_H
#define _CLOVER_H_H

void __mult_clvh(rvech_t sc[3][4][2], rvech_t a[2][36]) {
  rvech_t x[2][6][2];
  rvech_t y[2][2];

  for (int c=0;c<3;c++)
    for (int ri=0;ri<2;ri++)
      for(int v=0; v < VLENS; v++) {
        x[0][0+c][ri].v[v] = sc[c][0][ri].v[v] + sc[c][2][ri].v[v];
        x[0][3+c][ri].v[v] = sc[c][1][ri].v[v] + sc[c][3][ri].v[v];
        x[1][0+c][ri].v[v] = sc[c][0][ri].v[v] - sc[c][2][ri].v[v];
        x[1][3+c][ri].v[v] = sc[c][1][ri].v[v] - sc[c][3][ri].v[v];
      }

  for (int i=0;i<2;i++){
    for (int v=0; v < VLENS; v++) {
      y[i][0].v[v] = a[i][0].v[v] * x[i][0][0].v[v] + 
        a[i][ 6].v[v] * x[i][1][0].v[v] +
        a[i][ 8].v[v] * x[i][2][0].v[v] +
        a[i][10].v[v] * x[i][3][0].v[v] +
        a[i][12].v[v] * x[i][4][0].v[v] +
        a[i][14].v[v] * x[i][5][0].v[v] -
        a[i][ 7].v[v] * x[i][1][1].v[v] - 
        a[i][ 9].v[v] * x[i][2][1].v[v] -
        a[i][11].v[v] * x[i][3][1].v[v] -
        a[i][13].v[v] * x[i][4][1].v[v] -
        a[i][15].v[v] * x[i][5][1].v[v];

      y[i][1].v[v] = a[i][0].v[v] * x[i][0][1].v[v] +
        a[i][ 6].v[v] * x[i][1][1].v[v] +
        a[i][ 8].v[v] * x[i][2][1].v[v] +
        a[i][10].v[v] * x[i][3][1].v[v] +
        a[i][12].v[v] * x[i][4][1].v[v] +
        a[i][14].v[v] * x[i][5][1].v[v] +
        a[i][ 7].v[v] * x[i][1][0].v[v] +
        a[i][ 9].v[v] * x[i][2][0].v[v] +
        a[i][11].v[v] * x[i][3][0].v[v] +
        a[i][13].v[v] * x[i][4][0].v[v] +
        a[i][15].v[v] * x[i][5][0].v[v];
    }
  }

  for (int v=0; v < VLENS; v++) {
    sc[0][0][0].v[v] = y[0][0].v[v] +  y[1][0].v[v];
    sc[0][0][1].v[v] = y[0][1].v[v] +  y[1][1].v[v];
    sc[0][2][0].v[v] = y[0][0].v[v] -  y[1][0].v[v];
    sc[0][2][1].v[v] = y[0][1].v[v] -  y[1][1].v[v];
  }

  for (int i=0;i<2;i++){
    for (int v=0; v < VLENS; v++) {
      y[i][0].v[v] = a[i][1].v[v] *  x[i][1][0].v[v]
        + a[i][ 6].v[v] *  x[i][0][0].v[v] 
        + a[i][16].v[v] *  x[i][2][0].v[v] 
        + a[i][18].v[v] *  x[i][3][0].v[v] 
        + a[i][20].v[v] *  x[i][4][0].v[v] 
        + a[i][22].v[v] *  x[i][5][0].v[v] 
        + a[i][ 7].v[v] *  x[i][0][1].v[v] 
        - a[i][17].v[v] *  x[i][2][1].v[v]
        - a[i][19].v[v] *  x[i][3][1].v[v]
        - a[i][21].v[v] *  x[i][4][1].v[v]
        - a[i][23].v[v] *  x[i][5][1].v[v];
      y[i][1].v[v] = a[i][1].v[v] *  x[i][1][1].v[v]
        + a[i][ 6].v[v] *  x[i][0][1].v[v] 
        + a[i][16].v[v] *  x[i][2][1].v[v] 
        + a[i][18].v[v] *  x[i][3][1].v[v] 
        + a[i][20].v[v] *  x[i][4][1].v[v] 
        + a[i][22].v[v] *  x[i][5][1].v[v] 
        - a[i][ 7].v[v] *  x[i][0][0].v[v] 
        + a[i][17].v[v] *  x[i][2][0].v[v] 
        + a[i][19].v[v] *  x[i][3][0].v[v] 
        + a[i][21].v[v] *  x[i][4][0].v[v] 
        + a[i][23].v[v] *  x[i][5][0].v[v];
    }
  }

  for (int v=0; v < VLENS; v++) {
    sc[1][0][0].v[v] = y[0][0].v[v] +  y[1][0].v[v];
    sc[1][0][1].v[v] = y[0][1].v[v] +  y[1][1].v[v];
    sc[1][2][0].v[v] = y[0][0].v[v] -  y[1][0].v[v];
    sc[1][2][1].v[v] = y[0][1].v[v] -  y[1][1].v[v];
  }

  for (int i=0;i<2;i++){
    for (int v=0; v < VLENS; v++) {
      y[i][0].v[v] = a[i][2].v[v] *  x[i][2][0].v[v]
        + a[i][ 8].v[v] *  x[i][0][0].v[v]
        + a[i][16].v[v] *  x[i][1][0].v[v]
        + a[i][24].v[v] *  x[i][3][0].v[v]
        + a[i][26].v[v] *  x[i][4][0].v[v]
        + a[i][28].v[v] *  x[i][5][0].v[v]
        + a[i][ 9].v[v] *  x[i][0][1].v[v]
        + a[i][17].v[v] *  x[i][1][1].v[v]
        - a[i][25].v[v] *  x[i][3][1].v[v]
        - a[i][27].v[v] *  x[i][4][1].v[v]
        - a[i][29].v[v] *  x[i][5][1].v[v];
      y[i][1].v[v] = a[i][2].v[v] *  x[i][2][1].v[v]
        + a[i][ 8].v[v] *  x[i][0][1].v[v]
        + a[i][16].v[v] *  x[i][1][1].v[v]
        + a[i][24].v[v] *  x[i][3][1].v[v]
        + a[i][26].v[v] *  x[i][4][1].v[v]
        + a[i][28].v[v] *  x[i][5][1].v[v]
        - a[i][ 9].v[v] *  x[i][0][0].v[v]
        - a[i][17].v[v] *  x[i][1][0].v[v]
        + a[i][25].v[v] *  x[i][3][0].v[v]
        + a[i][27].v[v] *  x[i][4][0].v[v]
        + a[i][29].v[v] *  x[i][5][0].v[v];
    }
  }

  for (int v=0; v < VLENS; v++) {
    sc[2][0][0].v[v] = y[0][0].v[v] +  y[1][0].v[v];
    sc[2][0][1].v[v] = y[0][1].v[v] +  y[1][1].v[v];
    sc[2][2][0].v[v] = y[0][0].v[v] -  y[1][0].v[v];
    sc[2][2][1].v[v] = y[0][1].v[v] -  y[1][1].v[v];
  }

  for (int i=0;i<2;i++){
    for (int v=0; v < VLENS; v++) {
      y[i][0].v[v] = a[i][3].v[v] *  x[i][3][0].v[v]
        + a[i][10].v[v] *  x[i][0][0].v[v]
        + a[i][18].v[v] *  x[i][1][0].v[v]
        + a[i][24].v[v] *  x[i][2][0].v[v]
        + a[i][30].v[v] *  x[i][4][0].v[v]
        + a[i][32].v[v] *  x[i][5][0].v[v]
        + a[i][11].v[v] *  x[i][0][1].v[v]
        + a[i][19].v[v] *  x[i][1][1].v[v]
        + a[i][25].v[v] *  x[i][2][1].v[v]
        - a[i][31].v[v] *  x[i][4][1].v[v]
        - a[i][33].v[v] *  x[i][5][1].v[v];
      y[i][1].v[v] = a[i][3].v[v] *  x[i][3][1].v[v]
        + a[i][10].v[v] *  x[i][0][1].v[v]
        + a[i][18].v[v] *  x[i][1][1].v[v]
        + a[i][24].v[v] *  x[i][2][1].v[v]
        + a[i][30].v[v] *  x[i][4][1].v[v]
        + a[i][32].v[v] *  x[i][5][1].v[v]
        - a[i][11].v[v] *  x[i][0][0].v[v]
        - a[i][19].v[v] *  x[i][1][0].v[v]
        - a[i][25].v[v] *  x[i][2][0].v[v]
        + a[i][31].v[v] *  x[i][4][0].v[v]
        + a[i][33].v[v] *  x[i][5][0].v[v];
    }
  }

  for (int v=0; v < VLENS; v++) {
    sc[0][1][0].v[v] = y[0][0].v[v] +  y[1][0].v[v];
    sc[0][1][1].v[v] = y[0][1].v[v] +  y[1][1].v[v];
    sc[0][3][0].v[v] = y[0][0].v[v] -  y[1][0].v[v];
    sc[0][3][1].v[v] = y[0][1].v[v] -  y[1][1].v[v];
  }

  for (int i=0;i<2;i++){
    for (int v=0; v < VLENS; v++) {
      y[i][0].v[v] = a[i][4].v[v] *  x[i][4][0].v[v]
        + a[i][12].v[v] *  x[i][0][0].v[v]
        + a[i][20].v[v] *  x[i][1][0].v[v]
        + a[i][26].v[v] *  x[i][2][0].v[v]
        + a[i][30].v[v] *  x[i][3][0].v[v]
        + a[i][34].v[v] *  x[i][5][0].v[v]
        + a[i][13].v[v] *  x[i][0][1].v[v]
        + a[i][21].v[v] *  x[i][1][1].v[v]
        + a[i][27].v[v] *  x[i][2][1].v[v]
        + a[i][31].v[v] *  x[i][3][1].v[v]
        - a[i][35].v[v] *  x[i][5][1].v[v];
      y[i][1].v[v] = a[i][4].v[v] *  x[i][4][1].v[v]
        + a[i][12].v[v] *  x[i][0][1].v[v]
        + a[i][20].v[v] *  x[i][1][1].v[v]
        + a[i][26].v[v] *  x[i][2][1].v[v]
        + a[i][30].v[v] *  x[i][3][1].v[v]
        + a[i][34].v[v] *  x[i][5][1].v[v]
        - a[i][13].v[v] *  x[i][0][0].v[v]
        - a[i][21].v[v] *  x[i][1][0].v[v]
        - a[i][27].v[v] *  x[i][2][0].v[v]
        - a[i][31].v[v] *  x[i][3][0].v[v]
        + a[i][35].v[v] *  x[i][5][0].v[v];
    }
  }

  for (int v=0; v < VLENS; v++) {
    sc[1][1][0].v[v] = y[0][0].v[v] +  y[1][0].v[v];
    sc[1][1][1].v[v] = y[0][1].v[v] +  y[1][1].v[v];
    sc[1][3][0].v[v] = y[0][0].v[v] -  y[1][0].v[v];
    sc[1][3][1].v[v] = y[0][1].v[v] -  y[1][1].v[v];
  }

  for (int i=0;i<2;i++){
    for (int v=0; v < VLENS; v++) {
      y[i][0].v[v] = a[i][5].v[v] *  x[i][5][0].v[v]
        + a[i][14].v[v] *  x[i][0][0].v[v]
        + a[i][22].v[v] *  x[i][1][0].v[v]
        + a[i][28].v[v] *  x[i][2][0].v[v]
        + a[i][32].v[v] *  x[i][3][0].v[v]
        + a[i][34].v[v] *  x[i][4][0].v[v]
        + a[i][15].v[v] *  x[i][0][1].v[v]
        + a[i][23].v[v] *  x[i][1][1].v[v]
        + a[i][29].v[v] *  x[i][2][1].v[v]
        + a[i][33].v[v] *  x[i][3][1].v[v]
        + a[i][35].v[v] *  x[i][4][1].v[v];
      y[i][1].v[v] = a[i][5].v[v] *  x[i][5][1].v[v]
        + a[i][14].v[v] *  x[i][0][1].v[v]
        + a[i][22].v[v] *  x[i][1][1].v[v]
        + a[i][28].v[v] *  x[i][2][1].v[v]
        + a[i][32].v[v] *  x[i][3][1].v[v]
        + a[i][34].v[v] *  x[i][4][1].v[v]
        - a[i][15].v[v] *  x[i][0][0].v[v]
        - a[i][23].v[v] *  x[i][1][0].v[v]
        - a[i][29].v[v] *  x[i][2][0].v[v]
        - a[i][33].v[v] *  x[i][3][0].v[v]
        - a[i][35].v[v] *  x[i][4][0].v[v];
    }
  }

  for (int v=0; v < VLENS; v++) {
    sc[2][1][0].v[v] = y[0][0].v[v] +  y[1][0].v[v];
    sc[2][1][1].v[v] = y[0][1].v[v] +  y[1][1].v[v];
    sc[2][3][0].v[v] = y[0][0].v[v] -  y[1][0].v[v];
    sc[2][3][1].v[v] = y[0][1].v[v] -  y[1][1].v[v];
  }
}

#endif
