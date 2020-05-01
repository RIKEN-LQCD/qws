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

 */
//#ifdef _UTOFU_RANKMAP
#ifdef _UTOFU_RANKMAP

#include <stdint.h>
#include <stdio.h>
#include "rankmap_list.h"

#ifdef _USE_MPI
#include <mpi.h>
#endif

enum {DirX_=0, DirY_, DirZ_, DirA_, DirB_, DirC_};
const char *tofu_char="XYZABC";


/*
 helper functions
*/
void get_relative_coords(int *rcoords, int *size,
                         const uint8_t *my_coords,
                         const uint8_t *coords_org, const uint8_t *coords_size,
                         const uint8_t *coords_min, const uint8_t *coords_max){
  for(int i=0; i<6; i++){
    if(my_coords[i] < coords_org[i]){
      rcoords[i] = (coords_max[i] + 1 + my_coords[i] - coords_min[i]) - coords_org[i];
    } else {
      rcoords[i] = my_coords[i] - coords_org[i];
    }
    size[i]=coords_size[i];
  }
}

void to_absoulte_coords(uint8_t (*neighbor_coords)[6],
                        const uint8_t *coords_org, const uint8_t *coords_size,
                        const uint8_t *coords_min, const uint8_t *coords_max){

  for(int dir=0; dir <4; dir++){
    for(int i=0; i<6; i++){
      neighbor_coords[dir][i]+=coords_org[i];
      if(neighbor_coords[dir][i]>coords_max[i]){
        neighbor_coords[dir][i]-=(coords_max[i]+1);
        neighbor_coords[dir][i]+=coords_min[i];
      }
    }
  }
}


int node_coords_is_ok(int myrank, const int *rcoords, int node_x, int node_y, int node_z, int node_t){
  if(node_x<0){
    fprintf(stderr, "rank %d: cannot happen: no node_x is found: rcoords=%d %d %d %d %d %d\n", myrank,
            rcoords[0],rcoords[1], rcoords[2], rcoords[3], rcoords[4], rcoords[5]);
    return 0;
  }
  if(node_y<0){
    fprintf(stderr, "rank %d: cannot happen: no node_y is found: rcoords=%d %d %d %d %d %d\n", myrank,
            rcoords[0],rcoords[1], rcoords[2], rcoords[3], rcoords[4], rcoords[5]);
    return 0;
  }
  if(node_z<0){
    fprintf(stderr, "rank %d: cannot happen: no node_z is found: rcoords=%d %d %d %d %d %d\n", myrank,
            rcoords[0],rcoords[1], rcoords[2], rcoords[3], rcoords[4], rcoords[5]);
    return 0;
  }
  if(node_t<0){
    fprintf(stderr, "rank %d: cannot happen: no node_t is found: rcoords=%d %d %d %d %d %d\n", myrank,
            rcoords[0],rcoords[1], rcoords[2], rcoords[3], rcoords[4], rcoords[5]);
    return 0;
  }
  return 1;
}


int lookup1(int rcoord1, const int *T1, int size){
  int node=-1;
  for(int i=0; i<size; i++){
    if(rcoord1 == T1[i]){
      node=i;
      break;
    }
  }
  return node;
}

int lookup2(int rcoord1, const int *T1, int rcoord2, const int *T2, int size){
  int node=-1;
  for(int i=0; i<size; i++){
    if(rcoord1 == T1[i]
       && rcoord2 == T2[i] ){
      node=i;
      break;
    }
  }
  return node;
}

#ifdef _DEBUG_TOFU_COORD
void dump_coords(int myrank, const uint8_t *my_coords,
                 const uint8_t (*positive_neighbor_coords)[6], const uint8_t (*negative_neighbor_coords)[6]){
  for(int dir=0; dir <4; dir++){
    const uint8_t *p=positive_neighbor_coords[dir];
    const uint8_t *n=negative_neighbor_coords[dir];
    const uint8_t *m=my_coords;
    fprintf(stderr, "  rank=%d [%d %d %d %d %d %d], dir=%d: positive_ncoords = %d %d %d %d %d %d, nagative_ncoords = %d %d %d %d %d %d\n",
            myrank, m[0], m[1], m[2], m[3], m[4], m[5],
            dir, p[0],p[1],p[2],p[3],p[4],p[5], n[0],n[1],n[2],n[3],n[4],n[5]);
  }
}
#endif


int get_tofu_coord_4d_eval2Jan(const int myrank, const uint8_t *my_coords,
                               const uint8_t *coords_org, const uint8_t *coords_size,
                               const uint8_t *coords_min, const uint8_t *coords_max,
                               int *rank_coord, int *rank_size,
                               uint8_t (*positive_neighbor_coords)[6], int *pos_rank_in_node,
                               uint8_t (*negative_neighbor_coords)[6], int *neg_rank_in_node,
                               int DirX, int DirY, int DirZ){

  const int DirA=DirA_;
  const int DirB=DirB_;
  const int DirC=DirC_;

  // 2x2 in each node
  const int proc_per_node=4;
  const int pny_in_node=2;
  const int pnz_in_node=2;

  if(myrank==0){
    printf("rankmap for Fugaku evaluation environment 2\n");
    printf("  coords_org: %d %d %d %d %d %d\n", coords_org[0], coords_org[1], coords_org[2], coords_org[3], coords_org[4], coords_org[5]);
    printf("  coords_size: %d %d %d %d %d %d\n", coords_size[0], coords_size[1], coords_size[2], coords_size[3], coords_size[4], coords_size[5]);
    printf("  coords_min: %d %d %d %d %d %d\n", coords_min[0], coords_min[1], coords_min[2], coords_min[3], coords_min[4], coords_min[5]);
    printf("  coords_max: %d %d %d %d %d %d\n", coords_max[0], coords_max[1], coords_max[2], coords_max[3], coords_max[4], coords_max[5]);
    fflush(0);
  }

  // for x
  int TA[6] ={0,1,1,1,0,0};
  int TZc[6]={0,0,1,2,2,1};
  int size_x = 6;  // node size in x-direction

  // for y
  //int TC[4]={0,1,1,0};
  //int TZd[4]={0,0,1,1} (  = {0,0,3,3}/3 )
  //int TC[6]={0,1,1,1,0,0};
  //int TZd[6]={0,0,1,2,2,1}( = {0,0,3,6,6,3}/3 )
  //...
  //int TC[14]={0,1,1,1,1,1,1,1,0,0,0,0,0,0};
  //int TZd[14]={0,0,1,2,3,4,5,6,6,5,4,3,2,1} ( = {0,0,3,6,9,12,15,18,18,15,12,9,6,3}/3 )
  // for size=16: TZ is torus
  int TC[16] ={0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0};
  int TZd[16]={0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7};

  int size_TZd=coords_size[DirZ]/3;
  int n=0;
  if(size_TZd==8){
    n=16;
  }else{
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
  }
  int size_y=n;  // node size in y-direction

  // for z
  int TB[3]={0,1,2};
  int size_z=3;  // node size in z-direction

  // for t
  int TX[24];
  int TY[24];
  //  int TX[4]={0,1o,1,0}; // for size[DirX]=2
  //  int TY[4]={0,0,1,1};
  //  int size_t=4;  // node size in t-direction
  n=0;
  TX[n]=0;
  TY[n]=0;
  n++;
  for(int t=0; t<coords_size[DirX]; t++){
    TX[n]=t;
    TY[n]=1;
    n++;
  }
  for(int t=coords_size[DirX]-1; t>0; t--){
    TX[n]=t;
    TY[n]=0;
    n++;
  }
  int size_t=n;

  int rank_in_node=myrank % proc_per_node;
  int py=rank_in_node%pny_in_node;
  int pz=rank_in_node/pny_in_node;

  // rank size in each direction
  rank_size[0] = size_x;
  rank_size[1] = pny_in_node*size_y;
  rank_size[2] = pnz_in_node*size_z;
  rank_size[3] = size_t;

  int rcoords[6]; // relative tofu coordinates
  int size[6];    // tofu size
  get_relative_coords(rcoords, size, my_coords, coords_org, coords_size, coords_min, coords_max);

  if(myrank==0){
    //    printf("rankmap for Fugaku evaluation environment 2\n");
    printf("tofu size (original): [A,B,C,X,Y,Z]=[%d,%d,%d,%d,%d,%d]\n",
            size[DirA_],size[DirB_],size[DirC_],size[DirX_],size[DirY_],size[DirZ_]);
    printf("  rotate: X,Y,Z --> %c,%c,%c\n", tofu_char[DirX], tofu_char[DirY], tofu_char[DirZ]);
    printf("tofu size (rotated): [A,B,C,X,Y,Z]=[%d,%d,%d,%d,%d,%d]\n",
            size[DirA],size[DirB],size[DirC],size[DirX],size[DirY],size[DirZ]);
    printf("----- node coordinates -----\n");
    printf("map for x:  TA  = ");
    for(int i=0; i<size_x; i++){printf("%3d",TA[i]); }
    printf("\n");
    printf("            TZc = ");
    for(int i=0; i<size_x; i++){printf("%3d",TZc[i]); }
    printf("\n");

    printf("map for y:  TC  = ");
    for(int i=0; i<size_y; i++){printf("%3d",TC[i]); }
    printf("\n");
    printf("            TZd = ");
    for(int i=0; i<size_y; i++){printf("%3d",TZd[i]); }
    printf("\n");

    printf("map for z:  TB  = ");
    for(int i=0; i<size_z; i++){printf("%3d",TB[i]); }
    printf("\n");

    printf("map for t:  TX  = ");
    for(int i=0; i<size_t; i++){printf("%3d",TX[i]); }
    printf("\n");
    printf("            TY  = ");
    for(int i=0; i<size_t; i++){printf("%3d",TY[i]); }
    printf("\n");
    printf("inner node size: 1 2 2 1\n");
    printf("----------------------------\n");
    fflush(0);
  }

  // divide Tz into continuous and discrete
  int rcoords_Zc=rcoords[DirZ]%3; // continuous
  int rcoords_Zd=rcoords[DirZ]/3; // discrete

  // look up the node coordinate
  int node_x=lookup2(rcoords[DirA], TA, rcoords_Zc, TZc, size_x);
  int node_y=lookup2(rcoords[DirC], TC, rcoords_Zd, TZd, size_y);
  int node_z=lookup1(rcoords[DirB], TB, size_z);
  int node_t=lookup2(rcoords[DirX], TX, rcoords[DirY], TY, size_t);
  if(!node_coords_is_ok(myrank, rcoords, node_x, node_y, node_z, node_t)){
    return RANKMAP_NODE_NOT_FOUND;
  }

  // logical coordinate
  rank_coord[0]=node_x;
  rank_coord[1]=2 * node_y + py;
  rank_coord[2]=2 * node_z + pz;
  rank_coord[3]=node_t;

  ///////////////////////////////////////////////////////////////////
  //
  // neighbors
  //
  for(int i=0; i<6; i++){
    for(int dir=0; dir<4; dir++){
      positive_neighbor_coords[dir][i]=rcoords[i];
      negative_neighbor_coords[dir][i]=rcoords[i];
    }
  }

  //
  // x-direction
  //
  int dir=0;
  int xf=(rank_coord[dir]+1) % rank_size[dir];
  int xb=(rank_coord[dir]-1+rank_size[dir]) % rank_size[dir];
  positive_neighbor_coords[dir][DirA]=TA[xf];
  positive_neighbor_coords[dir][DirZ]=3*rcoords_Zd + TZc[xf];
  negative_neighbor_coords[dir][DirA]=TA[xb];
  negative_neighbor_coords[dir][DirZ]=3*rcoords_Zd + TZc[xb];
  pos_rank_in_node[dir]=rank_in_node;
  neg_rank_in_node[dir]=rank_in_node;

  //
  // y-direction : has inner node ranks
  //
  dir=1;
  int yf=(rank_coord[dir]+1) % rank_size[dir];
  int yb=(rank_coord[dir]-1+rank_size[dir]) % rank_size[dir];
  int pyf=yf % pny_in_node;
  int pyb=yb % pny_in_node;
  yf/=pny_in_node;
  yb/=pny_in_node;
  positive_neighbor_coords[dir][DirC]=TC[yf];
  positive_neighbor_coords[dir][DirZ]=3*TZd[yf] + rcoords_Zc;
  negative_neighbor_coords[dir][DirC]=TC[yb];
  negative_neighbor_coords[dir][DirZ]=3*TZd[yb] + rcoords_Zc;
  pos_rank_in_node[dir] = pyf + pz*pny_in_node;
  neg_rank_in_node[dir] = pyb + pz*pny_in_node;

  //
  // z-direction : has inner node ranks
  //
  dir=2;
  int zf=(rank_coord[dir]+1) % rank_size[dir];
  int zb=(rank_coord[dir]-1+rank_size[dir]) % rank_size[dir];
  int pzf=zf % pnz_in_node;
  int pzb=zb % pnz_in_node;
  zf/=pnz_in_node;
  zb/=pnz_in_node;
  positive_neighbor_coords[dir][DirB]=TB[zf];
  negative_neighbor_coords[dir][DirB]=TB[zb];
  pos_rank_in_node[dir] = py + pzf*pny_in_node;
  neg_rank_in_node[dir] = py + pzb*pny_in_node;

  //
  // t-direction
  //
  dir=3;
  int tf=(rank_coord[dir]+1) % size_t;
  int tb=(rank_coord[dir]-1+size_t) % size_t;
  positive_neighbor_coords[dir][DirX]=TX[tf];
  positive_neighbor_coords[dir][DirY]=TY[tf];
  negative_neighbor_coords[dir][DirX]=TX[tb];
  negative_neighbor_coords[dir][DirY]=TY[tb];
  pos_rank_in_node[dir]=rank_in_node;
  neg_rank_in_node[dir]=rank_in_node;

  to_absoulte_coords(positive_neighbor_coords, coords_org, coords_size, coords_min, coords_max);
  to_absoulte_coords(negative_neighbor_coords, coords_org, coords_size, coords_min, coords_max);

#ifdef _DEBUG_TOFU_COORD
  dump_coords(myrank, my_coords, positive_neighbor_coords, negative_neighbor_coords);
#endif

  return RANKMAP_EVAL2Jan;
}



int get_tofu_coord_4d_eval2Feb(const int myrank, const uint8_t *my_coords,
                               const uint8_t *coords_org, const uint8_t *coords_size,
                               const uint8_t *coords_min, const uint8_t *coords_max,
                               int *rank_coord, int *rank_size,
                               uint8_t (*positive_neighbor_coords)[6], int *pos_rank_in_node,
                               uint8_t (*negative_neighbor_coords)[6], int *neg_rank_in_node,
                               int DirX, int DirY, int DirZ){

  const int DirA=DirA_;
  const int DirB=DirB_;
  const int DirC=DirC_;

  // 2x2 in each node
  const int proc_per_node=4;
  const int pny_in_node=2;
  const int pnz_in_node=2;

  if(myrank==0){
    printf("rankmap for Fugaku evaluation environment 2 (2nd round)\n");
    printf("  coords_org: %d %d %d %d %d %d\n", coords_org[0], coords_org[1], coords_org[2], coords_org[3], coords_org[4], coords_org[5]);
    printf("  coords_size: %d %d %d %d %d %d\n", coords_size[0], coords_size[1], coords_size[2], coords_size[3], coords_size[4], coords_size[5]);
    printf("  coords_min: %d %d %d %d %d %d\n", coords_min[0], coords_min[1], coords_min[2], coords_min[3], coords_min[4], coords_min[5]);
    printf("  coords_max: %d %d %d %d %d %d\n", coords_max[0], coords_max[1], coords_max[2], coords_max[3], coords_max[4], coords_max[5]);
    fflush(0);
  }

  // for x
  int TA[6] ={0,1,1,1,0,0};
  int TZc[6]={0,0,1,2,2,1};
  int size_x = 6;  // node size in x-direction

  // for y
  //int TC[4]={0,1,1,0};
  //int TZd[4]={0,0,1,1} (  = {0,0,3,3}/3 )
  //int TC[6]={0,1,1,1,0,0};
  //int TZd[6]={0,0,1,2,2,1}( = {0,0,3,6,6,3}/3 )
  //...
  //int TC[14]={0,1,1,1,1,1,1,1,0,0,0,0,0,0};
  //int TZd[14]={0,0,1,2,3,4,5,6,6,5,4,3,2,1} ( = {0,0,3,6,9,12,15,18,18,15,12,9,6,3}/3 )
  // for size=16: TZ is torus
  int TC[16] ={0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0};
  int TZd[16]={0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7};

  int size_TZd=coords_size[DirZ]/3;
  int n=0;
  if(size_TZd==8){
    n=16;
  }else{
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
  }
  int size_y=n;  // node size in y-direction

  // for z
  int TX[24]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
  int size_z=24;  // node size in z-direction

  // for t: to realize 192 lattice sites with local 3 sites, we need 64 nodes
  //        thus at least 22 nodes in TY direction
  int TB[64]={0, 1,2,2,1, 1,2,2,1, 1,2,2,1, 1,2,2,1, 1,2,2,1, 1,2,2,1,
              1,2,2,1, 1,2,2,1, 1,2,2,1, 1,2,2,1, 1,1,0,
              0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0};
  int TY[64]={0,0,1,1, 2,2,3,3, 4,4,5,5, 6,6,7,7, 8,8,9,9, 10,10,11,11,
              12,12,13,13, 14,14,15,15, 16,16,17,17, 18,18,19,19, 20,21,21,
              20,19,18,17, 16,15,14,13, 12,11,10,9, 8,7,6,5, 4,3,2,1};

  int size_t=3*coords_size[DirY]; // node size in t-direction
  if(size_t>64){
    size_t=64;
  } else {
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
  }
  int rank_in_node=myrank % proc_per_node;
  int py=rank_in_node%pny_in_node;
  int pz=rank_in_node/pny_in_node;

  // rank size in each direction
  rank_size[0] = size_x;
  rank_size[1] = pny_in_node*size_y;
  rank_size[2] = pnz_in_node*size_z;
  rank_size[3] = size_t;

  int rcoords[6]; // relative tofu coordinates
  int size[6];    // tofu size
  get_relative_coords(rcoords, size, my_coords, coords_org, coords_size, coords_min, coords_max);

  if(myrank==0){
    printf("tofu size (original): [A,B,C,X,Y,Z]=[%d,%d,%d,%d,%d,%d]\n",
            size[DirA_],size[DirB_],size[DirC_],size[DirX_],size[DirY_],size[DirZ_]);
    printf("  rotate: X,Y,Z --> %c,%c,%c\n", tofu_char[DirX], tofu_char[DirY], tofu_char[DirZ]);
    printf("tofu size (rotated): [A,B,C,X,Y,Z]=[%d,%d,%d,%d,%d,%d]\n",
            size[DirA],size[DirB],size[DirC],size[DirX],size[DirY],size[DirZ]);
    printf("----- node coordinates -----\n");
    printf("map for x:  TA  = ");
    for(int i=0; i<size_x; i++){printf("%3d",TA[i]); }
    printf("\n");
    printf("            TZc = ");
    for(int i=0; i<size_x; i++){printf("%3d",TZc[i]); }
    printf("\n");

    printf("map for y:  TC  = ");
    for(int i=0; i<size_y; i++){printf("%3d",TC[i]); }
    printf("\n");
    printf("            TZd = ");
    for(int i=0; i<size_y; i++){printf("%3d",TZd[i]); }
    printf("\n");

    printf("map for z:  TX  = ");
    for(int i=0; i<size_z; i++){printf("%3d",TX[i]); }
    printf("\n");

    printf("map for t:  TB  = ");
    for(int i=0; i<size_t; i++){printf("%3d",TB[i]); }
    printf("\n");
    printf("            TY  = ");
    for(int i=0; i<size_t; i++){printf("%3d",TY[i]); }
    printf("\n");
    printf("inner node size: 1 2 2 1\n");
    printf("----------------------------\n");
    fflush(0);
  }

  // divide Tz into continuous and discrete
  int rcoords_Zc=rcoords[DirZ]%3; // continuous
  int rcoords_Zd=rcoords[DirZ]/3; // discrete

  // look up the node coordinate
  int node_x=lookup2(rcoords[DirA], TA, rcoords_Zc, TZc, size_x);
  int node_y=lookup2(rcoords[DirC], TC, rcoords_Zd, TZd, size_y);
  int node_z=lookup1(rcoords[DirX], TX, size_z);
  int node_t=lookup2(rcoords[DirB], TB, rcoords[DirY], TY, size_t);
  if(!node_coords_is_ok(myrank, rcoords, node_x, node_y, node_z, node_t)){
    return RANKMAP_NODE_NOT_FOUND;
  }


  // logical coordinate
  rank_coord[0]=node_x;
  rank_coord[1]=pny_in_node * node_y + py;
  rank_coord[2]=pny_in_node * node_z + pz;
  rank_coord[3]=node_t;

  ///////////////////////////////////////////////////////////////////
  //
  // neighbors
  //
  for(int i=0; i<6; i++){
    for(int dir=0; dir<4; dir++){
      positive_neighbor_coords[dir][i]=rcoords[i];
      negative_neighbor_coords[dir][i]=rcoords[i];
    }
  }

  //
  // x-direction
  //
  int dir=0;
  int xf=(rank_coord[dir]+1) % rank_size[dir];
  int xb=(rank_coord[dir]-1+rank_size[dir]) % rank_size[dir];
  positive_neighbor_coords[dir][DirA]=TA[xf];
  positive_neighbor_coords[dir][DirZ]=3*rcoords_Zd + TZc[xf];
  negative_neighbor_coords[dir][DirA]=TA[xb];
  negative_neighbor_coords[dir][DirZ]=3*rcoords_Zd + TZc[xb];
  pos_rank_in_node[dir]=rank_in_node;
  neg_rank_in_node[dir]=rank_in_node;

  //
  // y-direction : has inner node ranks
  //
  dir=1;
  int yf=(rank_coord[dir]+1) % rank_size[dir];
  int yb=(rank_coord[dir]-1+rank_size[dir]) % rank_size[dir];
  int pyf=yf % pny_in_node;
  int pyb=yb % pny_in_node;
  yf/=pny_in_node;
  yb/=pny_in_node;
  pos_rank_in_node[dir] = pyf + pz*pny_in_node;
  neg_rank_in_node[dir] = pyb + pz*pny_in_node;

  positive_neighbor_coords[dir][DirC]=TC[yf];
  positive_neighbor_coords[dir][DirZ]=3*TZd[yf] + rcoords_Zc;
  negative_neighbor_coords[dir][DirC]=TC[yb];
  negative_neighbor_coords[dir][DirZ]=3*TZd[yb] + rcoords_Zc;

  //
  // z-direction : has inner node ranks
  //
  dir=2;
  int zf=(rank_coord[dir]+1) % rank_size[dir];
  int zb=(rank_coord[dir]-1+rank_size[dir]) % rank_size[dir];
  int pzf=zf % pnz_in_node;
  int pzb=zb % pnz_in_node;
  zf/=pnz_in_node;
  zb/=pnz_in_node;
  pos_rank_in_node[dir] = py + pzf*pny_in_node;
  neg_rank_in_node[dir] = py + pzb*pny_in_node;
  positive_neighbor_coords[dir][DirX]=TX[zf];
  negative_neighbor_coords[dir][DirX]=TX[zb];

  //
  // t-direction
  //
  dir=3;
  int tf=(rank_coord[dir]+1) % size_t;
  int tb=(rank_coord[dir]-1+size_t) % size_t;
  positive_neighbor_coords[dir][DirB]=TB[tf];
  positive_neighbor_coords[dir][DirY]=TY[tf];
  negative_neighbor_coords[dir][DirB]=TB[tb];
  negative_neighbor_coords[dir][DirY]=TY[tb];
  pos_rank_in_node[dir]=rank_in_node;
  neg_rank_in_node[dir]=rank_in_node;

  to_absoulte_coords(positive_neighbor_coords, coords_org, coords_size, coords_min, coords_max);
  to_absoulte_coords(negative_neighbor_coords, coords_org, coords_size, coords_min, coords_max);

#ifdef _DEBUG_TOFU_COORD
  dump_coords(myrank, my_coords, positive_neighbor_coords, negative_neighbor_coords);
#endif

  return RANKMAP_EVAL2Feb;
}


int get_tofu_coord_4d_240rack(const int myrank, const uint8_t *my_coords,
			      const uint8_t *coords_org, const uint8_t *coords_size,
			      const uint8_t *coords_min, const uint8_t *coords_max,
			      int *rank_coord, int *rank_size,
			      uint8_t (*positive_neighbor_coords)[6], int *pos_rank_in_node,
			      uint8_t (*negative_neighbor_coords)[6], int *neg_rank_in_node,
			      int DirX, int DirY, int DirZ){

  const int DirA=DirA_;
  const int DirB=DirB_;
  const int DirC=DirC_;

  // 2x2 in each node
  const int proc_per_node=4;
  const int pny_in_node=2;
  const int pnz_in_node=2;

  if(myrank==0){
    printf("rankmap for 240 racks\n");
    printf("  coords_org: %d %d %d %d %d %d\n", coords_org[0], coords_org[1], coords_org[2], coords_org[3], coords_org[4], coords_org[5]);
    printf("  coords_size: %d %d %d %d %d %d\n", coords_size[0], coords_size[1], coords_size[2], coords_size[3], coords_size[4], coords_size[5]);
    printf("  coords_min: %d %d %d %d %d %d\n", coords_min[0], coords_min[1], coords_min[2], coords_min[3], coords_min[4], coords_min[5]);
    printf("  coords_max: %d %d %d %d %d %d\n", coords_max[0], coords_max[1], coords_max[2], coords_max[3], coords_max[4], coords_max[5]);
    fflush(0);
  }

  // for x
  int TA[6] ={0,1,1,1,0,0};
  int TZc[6]={0,0,1,2,2,1};
  int size_x = 6;  // node size in x-direction

  // for y
  //int TC[4]={0,1,1,0};
  //int TZd[4]={0,0,1,1} (  = {0,0,3,3}/3 )
  //int TC[6]={0,1,1,1,0,0};
  //int TZd[6]={0,0,1,2,2,1}( = {0,0,3,6,6,3}/3 )
  //...
  //int TC[14]={0,1,1,1,1,1,1,1,0,0,0,0,0,0};
  //int TZd[14]={0,0,1,2,3,4,5,6,6,5,4,3,2,1} ( = {0,0,3,6,9,12,15,18,18,15,12,9,6,3}/3 )
  // for size=16: TZ is torus
  int TC[16] ={0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0};
  int TZd[16]={0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7};

  int size_TZd=coords_size[DirZ]/3;
  int n=0;
  if(size_TZd==8){
    n=16;
  }else{
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
  }
  int size_y=n;  // node size in y-direction

  // for z
  int TB[3]={0,1,2};
  int size_z=3;  // node size in z-direction

  // for t
  int TX[320]={0};
  int TY[320]={0};

  n=0;
  TX[n]=0;
  TY[n]=0;
  n++;
  for(int y=0; y<16; y+=2){
    for(int x=1; x<20; x++){
      TX[n]=x;
      TY[n]=y;
      n++;
    }
    for(int x=19; x>0; x--){
      TX[n]=x;
      TY[n]=y+1;
      n++;
    }
  }
  for(int y=15; y>0; y--){
    TX[n]=0;
    TY[n]=y;
    n++;
  }
  int size_t=n;

  int rank_in_node=myrank % proc_per_node;
  int py=rank_in_node%pny_in_node;
  int pz=rank_in_node/pny_in_node;

  // rank size in each direction
  rank_size[0] = size_x;
  rank_size[1] = pny_in_node*size_y;
  rank_size[2] = pnz_in_node*size_z;
  rank_size[3] = size_t;

  int rcoords[6]; // relative tofu coordinates
  int size[6];    // tofu size
  get_relative_coords(rcoords, size, my_coords, coords_org, coords_size, coords_min, coords_max);

  if(myrank==0){
    //    printf("rankmap for Fugaku evaluation environment 2\n");
    printf("tofu size (original): [A,B,C,X,Y,Z]=[%d,%d,%d,%d,%d,%d]\n",
            size[DirA_],size[DirB_],size[DirC_],size[DirX_],size[DirY_],size[DirZ_]);
    printf("  rotate: X,Y,Z --> %c,%c,%c\n", tofu_char[DirX], tofu_char[DirY], tofu_char[DirZ]);
    printf("tofu size (rotated): [A,B,C,X,Y,Z]=[%d,%d,%d,%d,%d,%d]\n",
            size[DirA],size[DirB],size[DirC],size[DirX],size[DirY],size[DirZ]);
    printf("----- node coordinates -----\n");
    printf("map for x:  TA  = ");
    for(int i=0; i<size_x; i++){printf("%3d",TA[i]); }
    printf("\n");
    printf("            TZc = ");
    for(int i=0; i<size_x; i++){printf("%3d",TZc[i]); }
    printf("\n");

    printf("map for y:  TC  = ");
    for(int i=0; i<size_y; i++){printf("%3d",TC[i]); }
    printf("\n");
    printf("            TZd = ");
    for(int i=0; i<size_y; i++){printf("%3d",TZd[i]); }
    printf("\n");

    printf("map for z:  TB  = ");
    for(int i=0; i<size_z; i++){printf("%3d",TB[i]); }
    printf("\n");

    printf("map for t:  TX  = ");
    for(int i=0; i<size_t; i++){printf("%3d",TX[i]); }
    printf("\n");
    printf("            TY  = ");
    for(int i=0; i<size_t; i++){printf("%3d",TY[i]); }
    printf("\n");
    printf("inner node size: 1 2 2 1\n");
    printf("----------------------------\n");
    fflush(0);
  }

  // divide Tz into continuous and discrete
  int rcoords_Zc=rcoords[DirZ]%3; // continuous
  int rcoords_Zd=rcoords[DirZ]/3; // discrete

  // look up the node coordinate
  int node_x=lookup2(rcoords[DirA], TA, rcoords_Zc, TZc, size_x);
  int node_y=lookup2(rcoords[DirC], TC, rcoords_Zd, TZd, size_y);
  int node_z=lookup1(rcoords[DirB], TB, size_z);
  int node_t=lookup2(rcoords[DirX], TX, rcoords[DirY], TY, size_t);
  if(!node_coords_is_ok(myrank, rcoords, node_x, node_y, node_z, node_t)){
    return RANKMAP_NODE_NOT_FOUND;
  }

  // logical coordinate
  rank_coord[0]=node_x;
  rank_coord[1]=2 * node_y + py;
  rank_coord[2]=2 * node_z + pz;
  rank_coord[3]=node_t;

  ///////////////////////////////////////////////////////////////////
  //
  // neighbors
  //
  for(int i=0; i<6; i++){
    for(int dir=0; dir<4; dir++){
      positive_neighbor_coords[dir][i]=rcoords[i];
      negative_neighbor_coords[dir][i]=rcoords[i];
    }
  }

  //
  // x-direction
  //
  int dir=0;
  int xf=(rank_coord[dir]+1) % rank_size[dir];
  int xb=(rank_coord[dir]-1+rank_size[dir]) % rank_size[dir];
  positive_neighbor_coords[dir][DirA]=TA[xf];
  positive_neighbor_coords[dir][DirZ]=3*rcoords_Zd + TZc[xf];
  negative_neighbor_coords[dir][DirA]=TA[xb];
  negative_neighbor_coords[dir][DirZ]=3*rcoords_Zd + TZc[xb];
  pos_rank_in_node[dir]=rank_in_node;
  neg_rank_in_node[dir]=rank_in_node;

  //
  // y-direction : has inner node ranks
  //
  dir=1;
  int yf=(rank_coord[dir]+1) % rank_size[dir];
  int yb=(rank_coord[dir]-1+rank_size[dir]) % rank_size[dir];
  int pyf=yf % pny_in_node;
  int pyb=yb % pny_in_node;
  yf/=pny_in_node;
  yb/=pny_in_node;
  positive_neighbor_coords[dir][DirC]=TC[yf];
  positive_neighbor_coords[dir][DirZ]=3*TZd[yf] + rcoords_Zc;
  negative_neighbor_coords[dir][DirC]=TC[yb];
  negative_neighbor_coords[dir][DirZ]=3*TZd[yb] + rcoords_Zc;
  pos_rank_in_node[dir] = pyf + pz*pny_in_node;
  neg_rank_in_node[dir] = pyb + pz*pny_in_node;

  //
  // z-direction : has inner node ranks
  //
  dir=2;
  int zf=(rank_coord[dir]+1) % rank_size[dir];
  int zb=(rank_coord[dir]-1+rank_size[dir]) % rank_size[dir];
  int pzf=zf % pnz_in_node;
  int pzb=zb % pnz_in_node;
  zf/=pnz_in_node;
  zb/=pnz_in_node;
  positive_neighbor_coords[dir][DirB]=TB[zf];
  negative_neighbor_coords[dir][DirB]=TB[zb];
  pos_rank_in_node[dir] = py + pzf*pny_in_node;
  neg_rank_in_node[dir] = py + pzb*pny_in_node;

  //
  // t-direction
  //
  dir=3;
  int tf=(rank_coord[dir]+1) % size_t;
  int tb=(rank_coord[dir]-1+size_t) % size_t;
  positive_neighbor_coords[dir][DirX]=TX[tf];
  positive_neighbor_coords[dir][DirY]=TY[tf];
  negative_neighbor_coords[dir][DirX]=TX[tb];
  negative_neighbor_coords[dir][DirY]=TY[tb];
  pos_rank_in_node[dir]=rank_in_node;
  neg_rank_in_node[dir]=rank_in_node;

  to_absoulte_coords(positive_neighbor_coords, coords_org, coords_size, coords_min, coords_max);
  to_absoulte_coords(negative_neighbor_coords, coords_org, coords_size, coords_min, coords_max);

#ifdef _DEBUG_TOFU_COORD
  dump_coords(myrank, my_coords, positive_neighbor_coords, negative_neighbor_coords);
#endif

  return RANKMAP_240RACK;
}


int get_tofu_coord_4d(const int myrank, const uint8_t *my_coords, const uint8_t *coords_org, const uint8_t *coords_size,
                      const uint8_t *coords_min, const uint8_t *coords_max,
                      int *rank_coord, int *rank_size,
                      uint8_t (*positive_neighbor_coords)[6], int *pos_rank_in_node,
                      uint8_t (*negative_neighbor_coords)[6], int *neg_rank_in_node){

  // check the tofu size
  int size[6];    // tofu size
  int vol=1;
  for(int i=0; i<6; i++){
    size[i]=coords_size[i];
    vol*=size[i];
  }

  if(myrank==0){
    printf("get_tofu_coord_4d: [X,Y,Z,A,B,C]=[%d,%d,%d,%d,%d,%d]\n", size[DirX_], size[DirY_], size[DirZ_], size[DirA_], size[DirB_], size[DirC_]);
    fflush(0);
  }

#ifdef _USE_RANKMAP_240RACK
  { // for 240 rack
    if(size[DirA_] !=2 || size[DirB_] !=3 || size[DirC_] !=2
       || size[DirX_] != 20 || size[DirY_] !=16 || size[DirZ_] !=24 ){
      if(myrank==0){
      fprintf(stderr, "bad mpi size: must be [X,Y,Z,A,B,C]=[20,16,24,2,3,2]  (but [%d,%d,%d,%d,%d,%d])\n", size[DirX_], size[DirY_], size[DirZ_], size[DirA_], size[DirB_], size[DirC_]);
      }
#ifdef _USE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort();
#else
      abort();
#endif
    }

    return get_tofu_coord_4d_240rack(myrank, my_coords, coords_org, coords_size, coords_min, coords_max,
				     rank_coord, rank_size,
				     positive_neighbor_coords, pos_rank_in_node,
				     negative_neighbor_coords, neg_rank_in_node,
				     DirX_, DirY_, DirZ_);
  }
#endif

  if(size[DirA_] !=2 || size[DirB_] !=3 || size[DirC_] !=2 ){
    if(myrank==0){
      fprintf(stderr, "bad mpi size: must be [X,Y,Z,A,B,C]=[*,*,*,2,3,2]  (but [%d,%d,%d,%d,%d,%d])\n", size[DirX_], size[DirY_], size[DirZ_], size[DirA_], size[DirB_], size[DirC_]);
    }
    return -1;
  }


  if(  // (X,Y,Z,A,B,C)=(24,*,24,2,3,2)
     size[DirX_] == 24 && size[DirZ_] == 24
     && size[DirY_]>1) {
    if(myrank==0){
       printf("get_tofu_coord_4d: map as (X,Y,Z,A,B,C)=(24,*,24,2,3,2)\n");
     }
    return get_tofu_coord_4d_eval2Feb(myrank, my_coords, coords_org, coords_size, coords_min, coords_max,
                                      rank_coord, rank_size,
                                      positive_neighbor_coords, pos_rank_in_node,
                                      negative_neighbor_coords, neg_rank_in_node,
                                      DirX_, DirY_, DirZ_);
  }

  if(   // for 1 rack: (X,Y,Z,A,B,C)=(*,2,24,2,3,2)
     size[DirX_] > 1 && size[DirY_] == 2
     && size[DirZ_] == 24 ) {
    if(myrank==0){
       printf("get_tofu_coord_4d: map as (X,Y,Z,A,B,C)=(*,2,24,2,3,2)\n");
     }
    return get_tofu_coord_4d_eval2Jan(myrank, my_coords, coords_org, coords_size, coords_min, coords_max,
                                      rank_coord, rank_size,
                                      positive_neighbor_coords, pos_rank_in_node,
                                      negative_neighbor_coords, neg_rank_in_node,
                                      DirX_, DirY_, DirZ_);
  }

  if(   // for 1 rack: (X,Y,Z,A,B,C)=(24,2,*,2,3,2)
     size[DirX_] ==24  && size[DirY_] == 2
     && size[DirZ_] > 1 ) {
    if(myrank==0){
      printf("get_tofu_coord_4d: map as (X,Y,Z,A,B,C)=(24,2,*,2,3,2)\n");
    }
    return get_tofu_coord_4d_eval2Jan(myrank, my_coords, coords_org, coords_size, coords_min, coords_max,
                                      rank_coord, rank_size,
                                      positive_neighbor_coords, pos_rank_in_node,
                                      negative_neighbor_coords, neg_rank_in_node,
                                      DirZ_, DirY_, DirX_);
  }

  if(   // for 1 rack: (X,Y,Z,A,B,C)=(2,*,3n,2,3,2)
     size[DirX_] == 2 && size[DirY_] > 1
     && size[DirZ_] > 5 && (size[DirZ_] %3 ==0) ) {
    if(myrank==0){
       printf("get_tofu_coord_4d: map as (X,Y,Z,A,B,C)=(2,*,3n,2,3,2)\n");
     }
    return get_tofu_coord_4d_eval2Jan(myrank, my_coords, coords_org, coords_size, coords_min, coords_max,
                                      rank_coord, rank_size,
                                      positive_neighbor_coords, pos_rank_in_node,
                                      negative_neighbor_coords, neg_rank_in_node,
                                      DirY_, DirX_, DirZ_);
  }

  if(   // for 1 rack: (X,Y,Z,A,B,C)=(*,2,3n,2,3,2)
     size[DirX_] > 1 && size[DirY_] == 2
     && size[DirZ_] > 5 && (size[DirZ_] %3 ==0) ) {
    if(myrank==0){
       printf("get_tofu_coord_4d: map as (X,Y,Z,A,B,C)=(*,2,3n,2,3,2)\n");
     }
    return get_tofu_coord_4d_eval2Jan(myrank, my_coords, coords_org, coords_size, coords_min, coords_max,
                                      rank_coord, rank_size,
                                      positive_neighbor_coords, pos_rank_in_node,
                                      negative_neighbor_coords, neg_rank_in_node,
                                      DirX_, DirY_, DirZ_);
  }

  if(   // for 1 rack: (X,Y,Z,A,B,C)=(3n,*,2,2,3,2)
     size[DirZ_] == 2 && size[DirY_] >= 2
     && size[DirX_] > 5 && (size[DirX_] %3 ==0) ) {
    if(myrank==0){
       printf("get_tofu_coord_4d: map as (X,Y,Z,A,B,C)=(3n,*,2,2,3,2)\n");
     }
    return get_tofu_coord_4d_eval2Jan(myrank, my_coords, coords_org, coords_size, coords_min, coords_max,
                                      rank_coord, rank_size,
                                      positive_neighbor_coords, pos_rank_in_node,
                                      negative_neighbor_coords, neg_rank_in_node,
                                      DirY_, DirZ_, DirX_);
  }

  if(   // for 1 rack: (X,Y,Z,A,B,C)=(3n,2,*,2,3,2)
     size[DirZ_] > 1 && size[DirY_] == 2
     && size[DirX_] > 5 && (size[DirX_] %3 ==0) ) {
    if(myrank==0){
       printf("get_tofu_coord_4d: map as (X,Y,Z,A,B,C)=(3n,2,*,2,3,2)\n");
     }
    return get_tofu_coord_4d_eval2Jan(myrank, my_coords, coords_org, coords_size, coords_min, coords_max,
                                      rank_coord, rank_size,
                                      positive_neighbor_coords, pos_rank_in_node,
                                      negative_neighbor_coords, neg_rank_in_node,
                                      DirY_, DirZ_, DirX_);
  }




  if(myrank==0){
    fprintf(stderr, "no rank map is found:  [A,B,C,X,Y,Z]=[%d,%d,%d,%d,%d,%d]\n", size[DirA_], size[DirB_], size[DirC_], size[DirX_], size[DirY_], size[DirZ_]);
  }
  return -1;

}


#endif
