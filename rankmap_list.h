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
// macros for rankmaps
#ifndef RANKMAP_LIST_H_INCLUDED
#define RANKMAP_LIST_H_INCLUDED

#define RANKMAP_NODE_NOT_FOUND -999
#define RANKMAP_BAD_TOFU_SIZE -9999
#define RANKMAP_EVAL2Jan 1 // January
#define RANKMAP_EVAL2Feb 2 // Februay
#define RANKMAP_240RACK  3 // 20x16x24
#define RANKMAP_TOPOLOGY_Y 4 // use predefined TBxTY

#define RANKMAP_OPEN_XYZ 10
#define RANKMAP_OPEN_XY  11
#define RANKMAP_OPEN_Y   12

enum {DirX_=0, DirY_, DirZ_, DirA_, DirB_, DirC_};

#include <stdint.h>

int get_tofu_coord_and_tni_openXYZ(const int myrank, const uint8_t *my_coords,
                                   const uint8_t *coords_org, const uint8_t *coords_size,
                                   const uint8_t *coords_min, const uint8_t *coords_max,
                                   int *rank_coord, int *rank_size,
                                   uint8_t (*positive_neighbor_coords)[6], int *pos_rank_in_node,
                                   uint8_t (*negative_neighbor_coords)[6], int *neg_rank_in_node,
                                   int *tni_list,
                                   int DirX, int DirY, int DirZ);

int get_tofu_coord_and_tni_openXY(const int myrank, const uint8_t *my_coords,
                                  const uint8_t *coords_org, const uint8_t *coords_size,
                                  const uint8_t *coords_min, const uint8_t *coords_max,
                                  int *rank_coord, int *rank_size,
                                  uint8_t (*positive_neighbor_coords)[6], int *pos_rank_in_node,
                                  uint8_t (*negative_neighbor_coords)[6], int *neg_rank_in_node,
                                  int *tni_list,
                                  int DirX, int DirY, int DirZ);

int get_tofu_coord_and_tni_openY(const int myrank, const uint8_t *my_coords,
                                 const uint8_t *coords_org, const uint8_t *coords_size,
                                 const uint8_t *coords_min, const uint8_t *coords_max,
                                 int *rank_coord, int *rank_size,
                                 uint8_t (*positive_neighbor_coords)[6], int *pos_rank_in_node,
                                 uint8_t (*negative_neighbor_coords)[6], int *neg_rank_in_node,
                                 int *tni_list,
                                 int DirX, int DirY, int DirZ);

int get_tofu_coord_and_tni_openY_topology(const int myrank, const uint8_t *my_coords,
                                          const uint8_t *coords_org, const uint8_t *coords_size,
                                          const uint8_t *coords_min, const uint8_t *coords_max,
                                          int *rank_coord, int *rank_size,
                                          uint8_t (*positive_neighbor_coords)[6], int *pos_rank_in_node,
                                          uint8_t (*negative_neighbor_coords)[6], int *neg_rank_in_node,
                                          int *tni_list,
                                          int DirX, int DirY, int DirZ);

#endif
