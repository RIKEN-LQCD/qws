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
#include <stdint.h>
#include <stdio.h>
#include "rankmap_list.h"


int get_tni_set(int *tni,
                const int myrank,
                const int *rank_coords,
                const int *rank_size,
                const int **TNI_list, const int *TNI_size){

  for(int i=0; i<4; i++){
    int rc0=rank_coords[i];
    int rc1=rank_coords[i];
    if(TNI_size[2*i]==1) {rc0=0; }
    if(TNI_size[2*i+1]==1) {rc1=0; }
    tni[2*i]=TNI_list[2*i][rc0];
    tni[2*i+1]=TNI_list[2*i+1][rc1];
  }

  if(myrank==0){
    printf("tni map\n");
    for(int dir=0; dir<8; dir++){
      printf(" dir=%d:", dir);
      for(int i=0, ii=0, step=(TNI_size[dir]==1)?0:1;
          i<rank_size[dir/2]; i++, ii+=step) {
        printf(" %d",TNI_list[dir][ii]);
      };
      printf("\n");
    }
    fflush(0);
  }
  return 0;
}

int get_tni_default(int *tni, const int myrank,
                    const int *rank_coords,
                    const int *rank_size){
  if(myrank==0){
    printf("using the default tni assignment\n");
  }
  int TNI_Xp[1]={0}; // from +x direction
  int TNI_Xm[1]={0}; // from -x direction
  int TNI_Yp[1]={1}; // from +y direction
  int TNI_Ym[1]={2}; // from -y direction
  int TNI_Zp[1]={3}; // from +z direction
  int TNI_Zm[1]={4}; // from -z direction
  int TNI_Tp[1]={5}; // from +t direction
  int TNI_Tm[1]={5}; // from -t direction
  const int *TNI_list[8]={TNI_Xp, TNI_Xm, TNI_Yp, TNI_Ym,
                          TNI_Zp, TNI_Zm, TNI_Tp, TNI_Tm};
  const int TNI_size[8]={1,1,1,1,1,1,1,1};
  return get_tni_set(tni, myrank, rank_coords, rank_size, TNI_list, TNI_size);
}

int get_tni_eval2Jan(int *tni,
                  const int myrank,
                  const int *rank_coords,
                  const int *rank_size){

  if(myrank==0){
    printf("using tni assignment for fugaku evalutaion 2 (Jan. 20200228)\n");
  }
  if(rank_size[0] != 6
     || rank_size[1] < 8 || rank_size[1] % 2 == 1 || rank_size[1] > 32
     || rank_size[2] != 6
     || rank_size[3] < 4 || rank_size[3] % 2 == 1 ){
    if(myrank==0){
      fprintf(stderr, "bad rank size: %d %d %d %d\n",
              rank_size[0], rank_size[1], rank_size[2], rank_size[3]);
    }
    return -1;
  }


  // from +x direction = to -x direction etc.
  //   tni[0]  from +x, to -x
  //   tni[1]  from -x, to +x
  //     (comlib[0]: tni[0] from +x, tni[0] to -x)
  //     (comlib[1]: tni[1] from -x, tni[1] to +x)
  // after swap_vcq_for_sending:
  //   tni[0]  from +x, to +x
  //   tni[1]  from -x, to -x
  //     (comlib[0]: tni[0] from +x, tni[1] to -x)
  //     (comlib[1]: tni[1] from -x, tni[0] to +x)


  //int TNI_Xp[6]={0,0,0,0,0,0};
  //int TNI_Xm[6]={1,0,1,1,0,1};
  int TNI_Xp[6]={1,0,0,0,1,1};
  int TNI_Xm[6]={0,1,1,1,0,0};
  int TNI_Yp[32]={2,1,2,1, 2,1,2,1, 2,1,2,1, 2,1,2,1,
                  2,1,2,1, 2,1,2,1, 2,1,2,1, 2,1,2,1}; // 2 for loop back
  int TNI_Ym[32]={0,3,0,3, 0,3,0,3, 0,3,0,3, 0,3,0,3,
                  0,3,0,3, 0,3,0,3, 0,3,0,3, 0,3,0,3}; // 3 for loop back
  int TNI_Zp[1]={2};     // 2 for loop back as well
  int TNI_Zm[1]={3};     // 3 for loop back as well
  int TNI_Tp[1]={4};
  int TNI_Tm[1]={5};
  const int *TNI_list[8]={TNI_Xp, TNI_Xm, TNI_Yp, TNI_Ym,
                          TNI_Zp, TNI_Zm, TNI_Tp, TNI_Tm};
  const int TNI_size[8]={6,6,32,32,1,1,1,1};

  return get_tni_set(tni, myrank, rank_coords, rank_size, TNI_list, TNI_size);
}

int get_tni_eval2Feb(int *tni,
                  const int myrank,
                  const int *rank_coords,
                  const int *rank_size){

  if(myrank==0){
    printf("using tni assignment for fugaku evalutaion 2 (Feb. 20200228)\n");
  }
  if(rank_size[0] != 6
     || rank_size[1] < 8 || rank_size[1] % 2 == 1 || rank_size[1] > 32
     || rank_size[2] != 48
     || rank_size[3] > 12 ){
    if(myrank==0){
      fprintf(stderr, "bad rank size: %d %d %d %d\n",
              rank_size[0], rank_size[1], rank_size[2], rank_size[3]);
    }
    return -1;
  }

  // from +x direction = to -x direction etc.
  //   tni[0]  from +x, to -x
  //   tni[1]  from -x, to +x
  //     (comlib[0]: tni[0] from +x, tni[0] to -x)
  //     (comlib[1]: tni[1] from -x, tni[1] to +x)
  // after swap_vcq_for_sending:
  //   tni[0]  from +x, to +x
  //   tni[1]  from -x, to -x
  //     (comlib[0]: tni[0] from +x, tni[1] to -x)
  //     (comlib[1]: tni[1] from -x, tni[0] to +x)

  //  int TNI_Xp[6]={0,0,0,0,0,0};
  //  int TNI_Xm[6]={1,0,1,1,0,1};
  int TNI_Xp[6]={1,0,0,0,1,1};
  int TNI_Xm[6]={0,1,1,1,0,0};
  int TNI_Yp[32]={2,1,2,1, 2,1,2,1, 2,1,2,1, 2,1,2,1,
                  2,1,2,1, 2,1,2,1, 2,1,2,1, 2,1,2,1}; // 2 for loop back
  int TNI_Ym[32]={0,3,0,3, 0,3,0,3, 0,3,0,3, 0,3,0,3,
                  0,3,0,3, 0,3,0,3, 0,3,0,3, 0,3,0,3}; // 3 for loop back
  int TNI_Zp[1]={2};     // 2 for loop back as well
  int TNI_Zm[1]={3};     // 3 for loop back as well
  int TNI_Tp[1]={4};
  int TNI_Tm[1]={5};

  const int *TNI_list[8]={TNI_Xp, TNI_Xm, TNI_Yp, TNI_Ym,
                          TNI_Zp, TNI_Zm, TNI_Tp, TNI_Tm};
  const int TNI_size[8]={6,6,32,32,1,1,1,1};
  return get_tni_set(tni, myrank, rank_coords, rank_size,
                     TNI_list, TNI_size);
}


int get_tni_240rack(int *tni,
                  const int myrank,
                  const int *rank_coords,
                  const int *rank_size){

  if(myrank==0){
    printf("using tni assignment for 240/330 rack\n");
  }
  if(rank_size[0] != 6
     || rank_size[1] != 32
     || rank_size[2] != 6
     //     || rank_size[3] != 320 // 240 rack
     ){
    if(myrank==0){
      fprintf(stderr, "bad rank size: %d %d %d %d\n",
              rank_size[0], rank_size[1], rank_size[2], rank_size[3]);
    }
    abort();
  }


  int TNI_Xp[6]={1,0,0,0,1,1};
  int TNI_Xm[6]={0,1,1,1,0,0};
  int TNI_Yp[32]={2,1,2,1, 2,1,2,1, 2,1,2,1, 2,1,2,1,
		  2,1,2,1, 2,1,2,1, 2,1,2,1, 2,1,2,1}; // 2 for loop back
  int TNI_Ym[32]={0,3,0,3, 0,3,0,3, 0,3,0,3, 0,3,0,3,
                  0,3,0,3, 0,3,0,3, 0,3,0,3, 0,3,0,3}; // 3 for loop back
  int TNI_Zp[1]={2};     // 2 for loop back as well
  int TNI_Zm[1]={3};     // 3 for loop back as well
  int TNI_Tp[1]={4};
  int TNI_Tm[1]={5};

  const int *TNI_list[8]={TNI_Xp, TNI_Xm, TNI_Yp, TNI_Ym,
                          TNI_Zp, TNI_Zm, TNI_Tp, TNI_Tm};
  const int TNI_size[8]={6,6,32,32,1,1,1,1};

  return get_tni_set(tni, myrank, rank_coords, rank_size, TNI_list, TNI_size);


}

int get_tni_topology_y(int *tni,
                       const int myrank,
                       const int *rank_coords,
                       const int *rank_size){

  if(myrank==0){
    printf("using tni assignment for topology_y (20200507)\n");
  }
  // skip ranksize test

  // from +x direction = to -x direction etc.
  //   tni[0]  from +x, to -x
  //   tni[1]  from -x, to +x
  //     (comlib[0]: tni[0] from +x, tni[0] to -x)
  //     (comlib[1]: tni[1] from -x, tni[1] to +x)
  // after swap_vcq_for_sending:
  //   tni[0]  from +x, to +x
  //   tni[1]  from -x, to -x
  //     (comlib[0]: tni[0] from +x, tni[1] to -x)
  //     (comlib[1]: tni[1] from -x, tni[0] to +x)

  //  int TNI_Xp[6]={0,0,0,0,0,0};
  //  int TNI_Xm[6]={1,0,1,1,0,1};
  int TNI_Xp[6]={1,0,0,0,1,1};
  int TNI_Xm[6]={0,1,1,1,0,0};
  int TNI_Yp[32]={2,1,2,1, 2,1,2,1, 2,1,2,1, 2,1,2,1,
                  2,1,2,1, 2,1,2,1, 2,1,2,1, 2,1,2,1}; // 2 for loop back
  int TNI_Ym[32]={0,3,0,3, 0,3,0,3, 0,3,0,3, 0,3,0,3,
                  0,3,0,3, 0,3,0,3, 0,3,0,3, 0,3,0,3}; // 3 for loop back
  int TNI_Zp[1]={2};     // 2 for loop back as well
  int TNI_Zm[1]={3};     // 3 for loop back as well
  int TNI_Tp[1]={4};
  int TNI_Tm[1]={5};

  const int *TNI_list[8]={TNI_Xp, TNI_Xm, TNI_Yp, TNI_Ym,
                          TNI_Zp, TNI_Zm, TNI_Tp, TNI_Tm};
  const int TNI_size[8]={6,6,32,32,1,1,1,1};
  return get_tni_set(tni, myrank, rank_coords, rank_size,
                     TNI_list, TNI_size);
}

int get_tni_list(int *tni,
                 const int myrank,
                 const int *myrank_coords, const int *myrank_size,
                 const int flag){
  if(myrank==0){
    printf("get_tni_list: mapid=%d\n", flag);
    fflush(0);
  }

  if(flag<0){
    return get_tni_default(tni, myrank, myrank_coords, myrank_size);
  } else if (flag==RANKMAP_EVAL2Jan){
      return get_tni_eval2Jan(tni, myrank, myrank_coords, myrank_size);
  } else if (flag==RANKMAP_EVAL2Feb){
      return get_tni_eval2Feb(tni, myrank, myrank_coords, myrank_size);
  } else if (flag==RANKMAP_TOPOLOGY_Y){
      return get_tni_topology_y(tni, myrank, myrank_coords, myrank_size);
  } else if (flag==RANKMAP_240RACK){
      return get_tni_240rack(tni, myrank, myrank_coords, myrank_size);
  } else { // unkown
    if(myrank==0){
      fprintf(stderr, "unknown rankmap flag is given, using default tni list: flag=%d\n", flag);
    }
    return get_tni_default(tni, myrank, myrank_coords, myrank_size);
  }

  return -1; // cannot get here
}
