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
#include "rankmap_list.h"
#include "get_tofu_coord_common.h"
#include <stdint.h>
#include <stdio.h>

extern  FILE *para_outputfile;

/*
 helper functions
*/
void get_relative_coords(int *rcoords, int *size,
                         const uint8_t *my_coords,
                         const uint8_t *coords_org, const uint8_t *coords_size,
                         const uint8_t *coords_min, const uint8_t *coords_max){
  for(int i=0; i<6; i++){
    if(my_coords[i] < coords_org[i]){ // using torus
      rcoords[i] = (coords_max[i] + 1 + my_coords[i] - coords_min[i]) - coords_org[i];
    } else {
      rcoords[i] = my_coords[i] - coords_org[i];
    }
    size[i]=coords_size[i];
  }
}

void to_absolute_coords(uint8_t (*neighbor_coords)[6],
                        const uint8_t *coords_org, const uint8_t *coords_size,
                        const uint8_t *coords_min, const uint8_t *coords_max){

  for(int dir=0; dir <4; dir++){
    for(int i=0; i<6; i++){
      neighbor_coords[dir][i]+=coords_org[i];
      if(neighbor_coords[dir][i]>coords_max[i]){ // using the torus
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

void dump_coords(int myrank, const uint8_t *my_coords,
                 const uint8_t (*positive_neighbor_coords)[6], const uint8_t (*negative_neighbor_coords)[6]){
#ifdef _DEBUG_TOFU_COORD
  for(int dir=0; dir <4; dir++){
    const uint8_t *p=positive_neighbor_coords[dir];
    const uint8_t *n=negative_neighbor_coords[dir];
    const uint8_t *m=my_coords;
    fprintf(stderr, "  rank=%d [%d %d %d %d %d %d], dir=%d: positive_ncoords = %d %d %d %d %d %d, nagative_ncoords = %d %d %d %d %d %d\n",
            myrank, m[0], m[1], m[2], m[3], m[4], m[5],
            dir, p[0],p[1],p[2],p[3],p[4],p[5], n[0],n[1],n[2],n[3],n[4],n[5]);
  }
#endif
}


void print_coords_and_tni(int myrank, const uint8_t *my_coords,
                          const int *rank_coord,
                          const uint8_t (*positive_neighbor_coords)[6], const uint8_t (*negative_neighbor_coords)[6],
                          const int *tni_list){
  if(para_outputfile){
    fprintf(para_outputfile, "rank=%d\n", myrank);
    fprintf(para_outputfile, "   Logical coordinate:  %d %d %d %d\n",
            rank_coord[0], rank_coord[1], rank_coord[2], rank_coord[3]);
    fprintf(para_outputfile, "      Tofu coordinate:  %2d %2d %2d %2d %2d %2d\n",
            my_coords[0], my_coords[1], my_coords[2], my_coords[3], my_coords[4], my_coords[5]);

    for(int dir=0; dir <4; dir++){
      const uint8_t *p=positive_neighbor_coords[dir];
      const uint8_t *n=negative_neighbor_coords[dir];
      fprintf(para_outputfile, "dir=%d:  forward, tofu = %2d %2d %2d %2d %2d %2d , tni = %d\n",
              dir, p[0],p[1],p[2],p[3],p[4],p[5], tni_list[2*dir]);
      fprintf(para_outputfile, "dir=%d: backward, tofu = %2d %2d %2d %2d %2d %2d , tni = %d\n",
              dir, n[0],n[1],n[2],n[3],n[4],n[5], tni_list[2*dir+1]);
    }
    fflush(para_outputfile);
  }
}

void print_tni(int myrank, const int *rank_size, const int **tni_list_array){
  if(myrank==0){
    printf("----------------------------\n");
    printf("TNI_list\n");
    printf("  QX+: ");
    for(int i=0; i<rank_size[0]; i++) { printf(" %d", tni_list_array[0][i]); }
    printf("\n");
    printf("  QX-: ");
    for(int i=0; i<rank_size[0]; i++) { printf(" %d", tni_list_array[1][i]); }
    printf("\n");
    printf("  QY+: ");
    for(int i=0; i<rank_size[1]; i++) { printf(" %d", tni_list_array[2][i]); }
    printf("\n");
    printf("  QY-: ");
    for(int i=0; i<rank_size[1]; i++) { printf(" %d", tni_list_array[3][i]); }
    printf("\n");
    printf("  QZ+: ");
    for(int i=0; i<rank_size[2]; i++) { printf(" %d", tni_list_array[4][i]); }
    printf("\n");
    printf("  QZ-: ");
    for(int i=0; i<rank_size[2]; i++) { printf(" %d", tni_list_array[5][i]); }
    printf("\n");
    printf("  QT+: ");
    for(int i=0; i<rank_size[3]; i++) { printf(" %d", tni_list_array[6][i]); }
    printf("\n");
    printf("  QT-: ");
    for(int i=0; i<rank_size[3]; i++) { printf(" %d", tni_list_array[7][i]); }
    printf("\n");
    printf("----------------------------\n");
    fflush(stdout);
  }

}



void print_tofu_openXYZ(const int myrank, const int **Tmap, const int *Nsize, const uint8_t *coords_size,
                        const int DirX, const int DirY, const int DirZ){
  const int *TA=Tmap[0];
  const int *TB=Tmap[1];
  const int *TC=Tmap[2];
  const int *TX=Tmap[3];
  const int *TY=Tmap[4];
  const int *TZc=Tmap[5];
  const int *TZd=Tmap[6];

  const int size_x=Nsize[0];
  const int size_y=Nsize[1];
  const int size_z=Nsize[2];
  const int size_t=Nsize[3];

  const int DirA=DirA_;
  const int DirB=DirB_;
  const int DirC=DirC_;

  if(myrank==0){
    int size[6]={coords_size[0], coords_size[1], coords_size[2], coords_size[3], coords_size[4], coords_size[5]};
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

}



int set_neighbors(const int myrank, const uint8_t *my_coords,
                  const uint8_t *coords_org, const uint8_t *coords_size,
                  const uint8_t *coords_min, const uint8_t *coords_max,
                  int *rank_coord, int *rank_size,
                  uint8_t (*positive_neighbor_coords)[6], int *pos_rank_in_node,
                  uint8_t (*negative_neighbor_coords)[6], int *neg_rank_in_node,
                  const int *Dir,
                  const int **Tmap, const int *Nsize){

  // no need to rotate coords_org, coords_size, coords_min, coords_max
  // but need to rotate (i.e, DirX and DirB might be flipped) positive_neighbor_coords and negative_neighbor_coords

  const int proc_per_node=4;
  const int pny_in_node=2;
  const int pnz_in_node=2;

  const int *TA=Tmap[0];
  const int *TB=Tmap[1];
  const int *TC=Tmap[2];
  const int *TX=Tmap[3];
  const int *TY=Tmap[4];
  const int *TZc=Tmap[5];
  const int *TZd=Tmap[6];

  const int size_x=Nsize[0];
  const int size_y=Nsize[1];
  const int size_z=Nsize[2];
  const int size_t=Nsize[3];

  const int DirA=Dir[0];
  const int DirB=Dir[1];
  const int DirC=Dir[2];
  const int DirX=Dir[3];
  const int DirY=Dir[4];
  const int DirZ=Dir[5];


  ///////////////////////////////////////////////////////////////////
  //
  // find the 4d rank coordiante of this rank
  //
  // rank_in_node:
  // ^ QZ
  // |   2  3
  // |   0  1
  // ---------> QY
  // 2x2 in each node
  int rank_in_node=myrank % proc_per_node;
  int py=rank_in_node%pny_in_node;
  int pz=rank_in_node/pny_in_node;

  // rank size in each direction
  rank_size[0] = size_x;
  rank_size[1] = pny_in_node*size_y;
  rank_size[2] = pnz_in_node*size_z;
  rank_size[3] = size_t;

  ///////////////////////////////////////////////////////////////////
  //
  // move to the relative tofu coordiantes
  //
  int rcoords[6]; // relative tofu coordinates
  int size[6];    // tofu size
  get_relative_coords(rcoords, size, my_coords, coords_org, coords_size, coords_min, coords_max);

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
  positive_neighbor_coords[dir][DirX]=TX[zf];
  negative_neighbor_coords[dir][DirX]=TX[zb];
  pos_rank_in_node[dir] = py + pzf*pny_in_node;
  neg_rank_in_node[dir] = py + pzb*pny_in_node;

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

  ///////////////////////////////////////////////////////////////////
  //
  // back to the absolute tofu coordiantes
  //
  to_absolute_coords(positive_neighbor_coords, coords_org, coords_size, coords_min, coords_max);
  to_absolute_coords(negative_neighbor_coords, coords_org, coords_size, coords_min, coords_max);

  return 1;
}



void print_tofu(const int myrank, const int *Dir, const int **Tmap, const int *Nsize, const uint8_t *coords_size){
  // no need to rotate coords_size

  const int size_x=Nsize[0];
  const int size_y=Nsize[1];
  const int size_z=Nsize[2];
  const int size_t=Nsize[3];

  const int DirA=DirA_;
  const int DirC=DirC_;
  const int DirB=DirB_;

  if(Dir[3] == DirB_){ // shuriked map: TX and TB are flipped
    const int *TA=Tmap[0];
    const int *TX=Tmap[1]; //  not TB
    const int *TC=Tmap[2];
    const int *TB=Tmap[3]; //  not TX
    const int *TY=Tmap[4];
    const int *TZc=Tmap[5];
    const int *TZd=Tmap[6];
    const int DirX=Dir[1]; //  not Dir[3]
    const int DirY=Dir[4];
    const int DirZ=Dir[5];

    if(myrank==0){
      int size[6]={coords_size[0], coords_size[1], coords_size[2], coords_size[3], coords_size[4], coords_size[5]};
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
      fflush(stdout);
    }

  } else {   // original map for full system
    const int *TA=Tmap[0];
    const int *TB=Tmap[1];
    const int *TC=Tmap[2];
    const int *TX=Tmap[3];
    const int *TY=Tmap[4];
    const int *TZc=Tmap[5];
    const int *TZd=Tmap[6];
    const int DirX=Dir[3];
    const int DirY=Dir[4];
    const int DirZ=Dir[5];

    if(myrank==0){
      int size[6]={coords_size[0], coords_size[1], coords_size[2], coords_size[3], coords_size[4], coords_size[5]};
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
      fflush(stdout);
    }

  }


}
