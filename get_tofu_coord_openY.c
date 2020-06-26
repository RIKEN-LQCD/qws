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
/*
   rank map for 4-dim system
     Tofu X: 24
     Tofu Y: open
     Tofu Z: 24
     Tofu A,B,C: 2x3x2

 */
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include "rankmap_list.h"
#include "get_tofu_coord_common.h"


int get_tofu_coord_and_tni_openY(const int myrank, const uint8_t *my_coords,
                                 const uint8_t *coords_org, const uint8_t *coords_size,
                                 const uint8_t *coords_min, const uint8_t *coords_max,
                                 int *rank_coord, int *rank_size,
                                 uint8_t (*positive_neighbor_coords)[6], int *pos_rank_in_node,
                                 uint8_t (*negative_neighbor_coords)[6], int *neg_rank_in_node,
                                 int *tni_list,
                                 int DirX, int DirY, int DirZ){

  const int DirA=DirA_;
  const int DirB=DirB_;
  const int DirC=DirC_;


  if(myrank==0){
    printf("rankmap for open Y tofu axis (20200622)\n");
    printf("  coords_org: %d %d %d %d %d %d\n", coords_org[0], coords_org[1], coords_org[2], coords_org[3], coords_org[4], coords_org[5]);
    printf("  coords_size: %d %d %d %d %d %d\n", coords_size[0], coords_size[1], coords_size[2], coords_size[3], coords_size[4], coords_size[5]);
    printf("  coords_min: %d %d %d %d %d %d\n", coords_min[0], coords_min[1], coords_min[2], coords_min[3], coords_min[4], coords_min[5]);
    printf("  coords_max: %d %d %d %d %d %d\n", coords_max[0], coords_max[1], coords_max[2], coords_max[3], coords_max[4], coords_max[5]);
    fflush(0);
  }


  // policy
  //   TZ+ tni=0
  //   TZ- tni=1
  //
  /////////////////////////////////////////////
  // QX
  /////////////////////////////////////////////
  int TA[6] ={0,1,1,1,0,0};
  int TZc[6]={0,0,1,2,2,1};
  int size_x = 6;  // node size in x-direction

  // starting from "o"
  // QX+
  //   0     0
  //  +---> +---> +
  //  ^           | 0   ^ TA
  // 1|          \|     |
  //  o <---+ <---+     ---->TZc
  //       1     1
  int TNI_Xp[6]={1,0,0,0,1,1};

  // QX-  (coordinate: clockwise, sending direction: counter-clockwise)
  //       1    1
  //  + <---+ <---+
  // 1|           ^      ^ TA
  // \|           |0     |
  //  o---> +---> +      ---->TZc
  //    0     0
  int TNI_Xm[6]={0,1,1,1,0,0};


  /////////////////////////////////////////////
  // QY
  /////////////////////////////////////////////
  //int TC[6] ={0,1,1,1,0,0};
  //int TZd[6]={0,0,1,2,2,1};
  int size_TZd=coords_size[DirZ]/3;
  int size_y = 2*size_TZd;  // node size in x-direction

  int TC[16] ={-1};
  int TZd[16]={-1};
  int n=0;
  TC[n]=0;
  TZd[n]=0;
  n++;
  for(int z=0; z<size_TZd; z++){
    TC[n]=1;
    TZd[n]=z;
    n++;
  }
  for(int z=size_TZd-1;  z>0; z--){
    TC[n]=0;
    TZd[n]=z;
    n++;
  }
  assert(size_y == n);

  // QY+  : uses peridicity of TZd
  //        the loop back to the same node is not depicted in the figure
  //      0            0
  //  +---> +     + ---> +
  //  ^     |0    ^      |0          ^ TC
  // 0|    \|     |0    \|           |
  //  o     +---> +      +---> ...   ---->TZd
  //          0            0
  //  int TNI_Yp[12]={2,0, 2,0, 2,0, 2,0, 2,0,... 2,0};  // 2 for loop back;  TC, TZd
  int TNI_Yp[32]={0};
  n=0;
  while(n<size_y*2){
    TNI_Yp[n] = 2; n++; // loop back
    TNI_Yp[n] = 0; n++;
  }

  // QY-  : uses peridicity of TZd
  //        the loop back to the same node is not depicted in the figure
  //       1           1
  //  + <---+     + <---+
  // 1|     ^     |1    ^           ^ TC
  // \|     |1   \|     |1          |
  //  o     + <---+     + <--- ...  ---->TZd
  //            1
  //int TNI_Ym[12]={1,3, 1,3,..., 1,3};  // 3 for loop back
  int TNI_Ym[32]={1};
  n=0;
  while(n<size_y*2){
    TNI_Ym[n]=1; n++;
    TNI_Ym[n]=3; n++;  // loop back
  }



  /////////////////////////////////////////////
  // QZ
  /////////////////////////////////////////////
  int TX[24]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
  int size_z=24;  // node size in z-direction

  // QZ+: use torus in TX
  //      the loop back to the same node is not depicted in th figure
  //   2     2         2
  //  o---> +---> ... +---> +
  //  ^                     | 2
  //  |_____________________|
  //
  int TNI_Zp[48];
  n=0;
  while(n<size_z*2){
    TNI_Zp[n]=2; n++;   // 2 for loop back as well; TX
  }

  // QZ-: use torus in TX
  //      the loop back to the same node is not depicted in th figure
  //       3     3        3
  //  o <---+ <--- ... <---+
  // 3|                    ^
  //  |____________________|
  //
  int TNI_Zm[48];
  n=0;
  while(n<size_z*2){
    TNI_Zm[n]=3; n++;  // 3 for loop back as well
  }


  /////////////////////////////////////////////
  // QT
  /////////////////////////////////////////////
  int TB[66]={-1};
  int TY[66]={-1};
  int size_t=3*coords_size[DirY]; // node size in t-direction
  //  Y=1: TB= 0 12;                   TY= 0 00
  //  Y=2: TB= 0 1221 0;               TY= 0 0011 1
  //  Y=3: TB= 0 1221 12 00;           TY= 0 0011 22 21
  //  Y=4: TB= 0 1221 1221 000;        TY= 0 0011 2233 321
  //  Y=5: TB= 0 1221 1221 12 0000;    TY= 0 0011 2233 44 4321
  //  Y=6: TB= 0 1221 1221 1221 00000; TY= 0 0011 2233 4455 54321
  // ...
  ///
  int TB4[4]={1,2,2,1};
  n=0;
  TB[n]=0;
  TY[n]=0;
  n++;
  int ii=0;
  for(int i=0; i<2*coords_size[DirY];i++){
    TB[n]=TB4[ii];
    TY[n]=i/2;
    ii=(ii+1)%4;
    n++;
  }
  for(int i=coords_size[DirY]-1; i>0; i--){
    TB[n]=0;
    TY[n]=i;
    n++;
  }

  // QT+  : uses peridicity of TB
  //
  //   4            4               4
  //  +---> +     +---> +      ...---> +
  //  ^     |4    ^     |4    ^        |4       ^ TB
  // 4|    \|     |4   \|     |       \|        |
  //  +     +---> +     + --->+        +        ---->TY
  //  ^        4            4          |4
  // 4|                               \|
  //  o <---+ <---+ <---+ ...  ... <---+
  //       4     4     4              4
  //
  //  or      +---
  //          ^   |
  //        4 |   |4
  //      ... +   |
  //              |
  //              |
  //  ... <---+<--
  //        4
  int TNI_Tp[66]={4}; // TB, TY
  for(int i=0; i<size_t; i++) {TNI_Tp[i] = 4; }

  // QT-  : uses peridicity of TB
  //
  //     5            5               5
  //  + <---+     + <---+      ... <---+
  //  |5    ^     |5    ^     |5       ^        ^ TB
  // \|     |5   \|     |5   \|        |5       |
  //  +     + <---+     + <---+        +        ---->TY
  //  |5         5           5         ^
  // \|                                |5
  //  o---> +---> +---> + ...  ...---> +
  //   5      5    5               5
  //
  //  or      +<--
  //         5|   |
  //         \|   |
  //      ... +   |
  //              |
  //              |5
  //  ...---> +---
  //      5
  int TNI_Tm[66]={5}; // TB, TY
  for(int i=0; i<size_t; i++) {TNI_Tm[i] = 5; }

  const int Dir[6]={DirA,DirB,DirC,DirX,DirY,DirZ};
  const int *Tmap[7]={TA,TB,TC,TX,TY,TZc,TZd};
  const int Nsize[4]={size_x,size_y,size_z,size_t};
  print_tofu(myrank, Dir, Tmap, Nsize, coords_size);

  // find the logical rank coordinate of this rank and the tofu coordinates of the neighbors
  int mapid=set_neighbors(myrank, my_coords,
                          coords_org, coords_size, coords_min, coords_max,
                          rank_coord, rank_size,
                          positive_neighbor_coords, pos_rank_in_node,
                          negative_neighbor_coords, neg_rank_in_node,
                          Dir, Tmap, Nsize);
  if(mapid <0 ){ // something is wrong with rank coordinate
    return mapid;
  }

  // tni list
  const int *tni_list_full[8]={TNI_Xp, TNI_Xm, TNI_Yp, TNI_Ym, TNI_Zp, TNI_Zm, TNI_Tp, TNI_Tm};
  for(int dir2=0; dir2<8; dir2++){
    int coord=rank_coord[dir2/2];
    tni_list[dir2]=tni_list_full[dir2][coord];
  }
  if(myrank==0){
    printf("tni map (rankid=%d)\n", myrank);
    for(int dir2=0; dir2<8; dir2++){
      printf(" dir=%d:", dir2);
      for(int i=0; i<rank_size[dir2/2]; i++) {
        printf(" %d",tni_list_full[dir2][i]);
      }
      printf("\n");
    }
    fflush(stdout);
  }

  return RANKMAP_OPEN_Y;
}
