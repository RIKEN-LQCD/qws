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

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <mpi-ext.h>
#include <utofu.h>
#include "rankmap_list.h"
#include "get_tofu_coord_common.h"




#define TOFU_MAX_IN_1AXIS 32
int check_tofu_volume(const uint8_t *my_coords, uint8_t *coords_org, uint8_t *coords_size, uint8_t *coords_min, uint8_t *coords_max, int *np){


  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, np);

  MPI_Allreduce(my_coords, coords_min, 6, MPI_INTEGER1, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(my_coords, coords_max, 6, MPI_INTEGER1, MPI_MAX, MPI_COMM_WORLD);

  int my_occupation[6][TOFU_MAX_IN_1AXIS]={0};
  int occupied[6][TOFU_MAX_IN_1AXIS]={0};
  for(int i=0; i<6; i++){
    my_occupation[i][my_coords[i]]=1;
  }
  MPI_Allreduce(my_occupation, occupied, 6*TOFU_MAX_IN_1AXIS,
                MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  // count the number of nodes in each direction
  for(int i=0; i<6; i++) {
    int size=0;
    for(int n=0; n<TOFU_MAX_IN_1AXIS; n++){
      if(occupied[i][n]>0) { size++; }
    }
    coords_size[i] = size;
  }
  // look for the origin
  for(int i=0; i<6; i++){
    int org=0;
    int prev=0;
    for(int n=TOFU_MAX_IN_1AXIS; n>0; n--){
      if(occupied[i][n-1] < prev) {
        org=n;
        break;
      }
      prev=occupied[i][n-1];
    }
    coords_org[i]=org;
  }

  //
  // periodic condition:
  //   if coord[i] > coord_max[i]
  //     coord[i] = coords_min[i] + coord[i] % coord_max[i]
  //
  //   if coord[i] < coord_org[i]
  //     coord[i] = coord_max[i] + coord[i] - coord_min[i]
  //
  //  if tofu in i-th dir is torus and
  //     sytem         7 8 9 10 11 12 13 14 15
  //     *:using/-not  * * - -  -  -  -  *  *
  //  coords_min[i]=7
  //  coords_max[i]=15
  //  coords_org[i]=14
  //  coords_size[i]=4
  //    (N.B. in the above case, naive min/max of tofu coordinate
  //     becomes 7 and 15, respectively.)

  // check the FJMPI interface as well
  //int rel_coords[6];
  //  FJMPI_Topology_get_coords(MPI_COMM_WORLD, myrank, FJMPI_TOFU_REL, 6, rel_coords);
  //  printf("rank %d: tofu coords sys = %d %d %d %d %d %d; rel = %d %d %d %d %d %d\n",
  //         myrank,
  //         my_coords[0],my_coords[1],my_coords[2],my_coords[3],my_coords[4],my_coords[5],
  //         rel_coords[0],rel_coords[1],rel_coords[2],rel_coords[3],rel_coords[4],rel_coords[5]);

  int tofu_vol=1;
  for(int i=0; i<6; i++){
    tofu_vol*=coords_size[i];
  }
  if(tofu_vol*4 != *np){
    if(myrank==0){
      fprintf(stderr, "warning: allocated size is not (hyper-)rectangluer:\n");
      fprintf(stderr, "       np = %d, tofu_vol * 4 = %d\n", *np, tofu_vol*4);
      fprintf(stderr, "       dir: min  max  size origin\n");
      for(int i=0; i<6; i++){
        fprintf(stderr, "        %d: %3d  %3d  %3d  %3d\n", i, coords_min[i], coords_max[i], coords_size[i], coords_org[i]);
      }
    }
    return RANKMAP_TOPOLOGY_Y;

    // or find a largest hyper rectangluar which fit in the given nodes
    //  np = ...
  }

  if(coords_size[DirA_] !=2 || coords_size[DirB_] !=3 || coords_size[DirC_] !=2 ){
    if(myrank==0){
      fprintf(stderr, "bad mpi size: must be [X,Y,Z,A,B,C]=[*,*,*,2,3,2]  (but [%d,%d,%d,%d,%d,%d])\n",
              coords_size[DirX_], coords_size[DirY_], coords_size[DirZ_], coords_size[DirA_], coords_size[DirB_], coords_size[DirC_]);
    }
    return -1;
  }

  if(myrank==0){
    fprintf(stderr, "       np = %d, tofu_vol * 4 = %d\n",*np, tofu_vol*4);
    fprintf(stderr, "       dir: min  max  size origin\n");
    for(int i=0; i<6; i++){
      fprintf(stderr, "        %d: %3d  %3d  %3d  %3d\n", i, coords_min[i], coords_max[i], coords_size[i], coords_org[i]);
    }
  }

  return 0;
}



//int naive_rank_asignment(int *rank_coord, int *rank_size);

//
// choose and get the rankmap
//
int call_get_tofu_coord_and_tni(const int myrank, const uint8_t *my_coords,
                                const uint8_t *coords_org, const uint8_t *coords_size,
                                const uint8_t *coords_min, const uint8_t *coords_max,
                                int *rank_coord, int *rank_size,
                                uint8_t (*positive_neighbor_coords)[6], int *pos_rank_in_node,
                                uint8_t (*negative_neighbor_coords)[6], int *neg_rank_in_node,
                                int pre_mapid,
                                int *tni_list){
  int size[6];    // tofu size
  for(int i=0; i<6; i++){
    size[i]=coords_size[i];
  }

  if( pre_mapid == RANKMAP_TOPOLOGY_Y ){
    if ( size[DirX_] == 24 && size[DirZ_] == 24
         && size[DirA_] == 2 && size[DirC_] == 2){
      if(myrank==0){
        printf("get_tofu_coord_and_tni: RANKMAP_TOPOLOGY_Y is given, (X,Y,Z,A,B,C)=(24,*,24,2,*,2)\n");
        fprintf(stderr, "map for topology_y is not yet ready\n");
      }
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
      /*
      return get_tofu_coord_and_tni_topology_y(myrank, my_coords, coords_org, coords_size, coords_min, coords_max,
                                               rank_coord, rank_size,
                                               positive_neighbor_coords, pos_rank_in_node,
                                               negative_neighbor_coords, neg_rank_in_node,
                                               tni_list,
                                               DirX_, DirY_, DirZ_);
      */
    }
    if(myrank==0){
      printf("get_tofu_coord_and_tni: RANKMAP_TOPOLOGY_Y requires (X,Y,Z,A,B,C)=(24,*,24,2,*,2), but X and/or Y are not.\n");
    }
    return -1;
  }


  // open Y
  if(  // (X,Y,Z,A,B,C)=(24,*,24,2,3,2)
     (size[DirX_] == 24)
     && (size[DirZ_] == 24) ) {
    if(myrank==0){
       printf("get_tofu_coord_and_tni: map as (X,Y,Z,A,B,C)=(24,*,24,2,3,2)\n");
     }
    return get_tofu_coord_and_tni_openY(myrank, my_coords, coords_org, coords_size, coords_min, coords_max,
                                        rank_coord, rank_size,
                                        positive_neighbor_coords, pos_rank_in_node,
                                        negative_neighbor_coords, neg_rank_in_node,
                                        tni_list,
                                        DirX_, DirY_, DirZ_);
  }

  // open X,Y
  if(  // (X,Y,Z,A,B,C)=(*,*,24,2,3,2) with X x Y = even
     (size[DirX_] * size[DirY_]  % 2 == 0)
     && (size[DirZ_] == 24) ) {
    if(myrank==0){
       printf("get_tofu_coord_and_tni: map as (X,Y,Z,A,B,C)=(*,*,24,2,3,2) with X x Y = even\n");
     }
    return get_tofu_coord_and_tni_openXY(myrank, my_coords, coords_org, coords_size, coords_min, coords_max,
                                         rank_coord, rank_size,
                                         positive_neighbor_coords, pos_rank_in_node,
                                         negative_neighbor_coords, neg_rank_in_node,
                                         tni_list,
                                         DirX_, DirY_, DirZ_);
  }

  if(  // (X,Y,Z,A,B,C)=(24,*,*,2,3,2) with Z x Y = even
     (size[DirZ_] * size[DirY_]  % 2 == 0)
     && (size[DirY_] == 24) ) {
    if(myrank==0){
       printf("get_tofu_coord_and_tni: map as (X,Y,Z,A,B,C)=(24,*,*,2,3,2) with Z x Y = even\n");
     }
    return get_tofu_coord_and_tni_openXY(myrank, my_coords, coords_org, coords_size, coords_min, coords_max,
                                         rank_coord, rank_size,
                                         positive_neighbor_coords, pos_rank_in_node,
                                         negative_neighbor_coords, neg_rank_in_node,
                                         tni_list,
                                         DirZ_, DirY_, DirX_);
  }

  // open X,Y,Z
  if(  // (X,Y,Z,A,B,C)=(*,*,3n,2,3,2) with X x Y = even
     (size[DirX_] * size[DirY_]  % 2 == 0)
     && (size[DirZ_] % 3 == 0) ) {
    if(myrank==0){
       printf("get_tofu_coord_and_tni: map as (X,Y,Z,A,B,C)=(*,*,3n,2,3,2) with X x Y = even\n");
     }
    return get_tofu_coord_and_tni_openXYZ(myrank, my_coords, coords_org, coords_size, coords_min, coords_max,
                                          rank_coord, rank_size,
                                          positive_neighbor_coords, pos_rank_in_node,
                                          negative_neighbor_coords, neg_rank_in_node,
                                          tni_list,
                                          DirX_, DirY_, DirZ_);
  }

  if(    // (X,Y,Z,A,B,C)=(3n*,*,*,2,3,2) with Z x Y = even
     (size[DirZ_] * size[DirY_]  % 2 == 0)
     && size[DirX_] % 3 == 0) {
    if(myrank==0){
       printf("get_tofu_coord_and_tni: map as (X,Y,Z,A,B,C)=(3n,*,*,2,3,2) with Z x Y = even\n");
     }
    return get_tofu_coord_and_tni_openXYZ(myrank, my_coords, coords_org, coords_size, coords_min, coords_max,
                                          rank_coord, rank_size,
                                          positive_neighbor_coords, pos_rank_in_node,
                                          negative_neighbor_coords, neg_rank_in_node,
                                          tni_list,
                                          DirZ_, DirY_, DirX_);
  }




  // not found
  return -1;
}

int get_tofu_coord_and_tni(uint8_t *my_coords, int *rank_coord, int *rank_size,
                           uint8_t (*positive_neighbor_coords)[6], int *pos_rank_in_node,
                           uint8_t (*negative_neighbor_coords)[6], int *neg_rank_in_node,
                           int *tni_list){

  int myrank, np_available;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  int rc;

  // get tofu coordinate
  rc=utofu_query_my_coords(my_coords);
  if(UTOFU_SUCCESS != rc){
    fprintf(stderr, "rank %d: Failed at utofu_querry_my_coords()! rc=%d\n", myrank, rc);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  // check if the assigned network is a hyper-rectangular
  //   it also obtain the assinged tofu coordiate
  //   N.B:  if the XXX axis use the torus nature,
  //         coords_min[XXX] + coords_size[XXX] != coords_max[XXX]

  uint8_t coords_org[6], coords_size[6];
  uint8_t coords_min[6], coords_max[6];
  int pre_mapid=check_tofu_volume(my_coords, coords_org, coords_size,
                                  coords_min, coords_max, &np_available);
  if(pre_mapid<0){
    fprintf(stderr, "rank %d: Failed at check_tofu_volume()! err=%d\n", myrank, pre_mapid);
    //    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    return pre_mapid;
  }

  int mapid=call_get_tofu_coord_and_tni(myrank, my_coords, coords_org, coords_size, coords_min, coords_max,
                                        rank_coord, rank_size,
                                        positive_neighbor_coords, pos_rank_in_node,
                                        negative_neighbor_coords, neg_rank_in_node,
                                        pre_mapid,
                                        tni_list);


  if(mapid<0){
    if(myrank==0){
      fprintf(stderr, "no rank map is found:  [A,B,C,X,Y,Z]=[%d,%d,%d,%d,%d,%d]\n", coords_size[DirA_], coords_size[DirB_], coords_size[DirC_], coords_size[DirX_], coords_size[DirY_], coords_size[DirZ_]);
    }
  } else {
    print_coords_and_tni(myrank, my_coords, rank_coord,
                         positive_neighbor_coords, negative_neighbor_coords,
                         tni_list);

    if(myrank==0){
      int flag=0;
      for(int i=0; i<4; i++){
        if(rank_coord[i] != 0) { flag++; }
      }
      if(flag !=0 ){
        fprintf(stderr, "WARNING! master rank (id=0) is not the logical origin of the rankmap.");
      }
    }
    // do NOT shift the process coordinates, as the shift may break the implict assumption for TNI assignment
    //int r[4];
    //    for(int i=0; i<dim; i++){
    //      r[i]=rank_coord[i];
    //    }
    //    MPI_Bcast(r, sizeof(int)*dim, MPI_BYTE, 0, MPI_COMM_WORLD);
    //    for(int i=0; i<dim; i++){
    //      rank_coord[i]=(rank_coord[i]-r[i]+rank_size[i]) % rank_size[i];
  }

  return mapid;
}

