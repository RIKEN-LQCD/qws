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
#include <assert.h>

// tofu coordinate etc.: up to 4-dim
volatile int m_rankmap_lib_neighbor_ranks[8];
volatile int m_rankmap_lib_process_size[4];
volatile int m_rankmap_lib_process_coord[4];
volatile uint8_t m_rankmap_lib_neighbor_tofu[8][6];
volatile uint8_t m_rankmap_lib_mytofu[6];
volatile uint8_t m_rankmap_lib_mytofu_offset[6];
volatile int m_rankmap_lib_myrank_in_node;
volatile int m_rankmap_lib_dim;
volatile int m_rankmap_lib_neighbor_rank_in_node[8];
volatile int m_rankmap_lib_tni_list[8];


// get_tofu_coord.c
int get_tofu_coord_and_tni(uint8_t *my_coords, int *rank_coord, int *rank_size,
                   uint8_t (*positive_neighbor_coords)[6], int *pos_rank_in_node,
                   uint8_t (*negative_neighbor_coords)[6], int *neg_rank_in_node,
                   int *tni_list);

//int get_tofu_coord(uint8_t *my_coords, uint8_t *coords_offset,
//                   int *rank_coord, int *rank_size, 
//                   uint8_t (*positive_neighbor_coords)[6], int *pos_rank_in_node, 
//                   uint8_t (*negative_neighbor_coords)[6], int *neg_rank_in_node, 
//                   const int dim);

int get_tni_list_default(int *tni,
                         const int myrank,
                         const int *myrank_coords, const int *myrank_size,
                         const int flag);


///////////////////////////////////////////////
//
// helper private functions
//
///////////////////////////////////////////////

//
// sanity check: is cq_id correct?
//
void check_cq_id(utofu_vcq_hdl_t *local_vcq_hdl,
                 uint8_t *coords, int cq_id_guess, int comp_id, int tni_id){
  int rc;
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  // obtain id for local vcq
  utofu_vcq_id_t local_vcq_id0;
  rc=utofu_query_vcq_id(*local_vcq_hdl, &local_vcq_id0);
  if(UTOFU_SUCCESS != rc){
    fprintf(stderr, "rank %d: Failed at utofu_query_vcq_id()! rc=%d\n", myrank, rc);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  // extract information of the current vcq
  uint8_t coords0[6];
  utofu_tni_id_t tni_id0;
  utofu_cq_id_t  cq_id0;
  uint16_t extra_val0;
  rc=utofu_query_vcq_info(local_vcq_id0, coords0, &tni_id0, &cq_id0, &extra_val0);
  if(UTOFU_SUCCESS != rc){
    fprintf(stderr, "rank %d: Failed at utofu_query_vcq_info()! rc=%d\n", myrank, rc);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
#ifdef _RANKMAP_DEBUG
  fprintf(stderr, "rank %d: local_vcq_id=%d, coords=[%d %d %d %d %d %d], tni_id=%d, cq_id=%d, extra_val=%d\n",
          myrank,
          local_vcq_id0, coords0[0],coords0[1],coords0[2],coords0[3],coords0[4],coords0[5],
          tni_id0, cq_id0, extra_val0);
#endif

  // another vcq_id: is this the same as the origianl local_vcq_id0?
  utofu_vcq_id_t local_vcq_id2;
  //    int cq_id2=3*(m_rdma_comlib_myrank % proc_per_node);  // guess
  rc=utofu_construct_vcq_id(coords, tni_id, cq_id_guess, comp_id, &local_vcq_id2);
  if(UTOFU_SUCCESS != rc){
    fprintf(stderr, "rank %d: Failed at utofu_construct_vcq_id()! rc=%d, coords=%d %d %d %d %d %d\n", myrank, rc,
            coords[0],coords[1], coords[2], coords[3], coords[4], coords[5]);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  // check if the caclated vcq_id is ok
  if(local_vcq_id0 != local_vcq_id2){
    fprintf(stderr, "rank %d: imporooper vcq_id: local_vcq_id0=%d, local_vcq_id2=%d; cq_id0=%d, cq_id_guess=%d\n", myrank, 
            local_vcq_id0, local_vcq_id2, cq_id0, cq_id_guess);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

}

int check_neighbors(){
  int err=0;
  int dim=m_rankmap_lib_dim;
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  uint8_t sbuffer[6], rbuffer[6];
  for(int i=0; i<6; i++){
    sbuffer[i]=m_rankmap_lib_mytofu[i];
  }
  const uint8_t *sbuf=sbuffer;
  uint8_t *rbuf=rbuffer;
#ifdef _RANKMAP_DEBUG
  if(dim==2){
    fprintf(stderr, "rank %d: neighbors = %d %d %d %d\n",
            myrank,
            m_rankmap_lib_neighbor_ranks[0],
            m_rankmap_lib_neighbor_ranks[1],
            m_rankmap_lib_neighbor_ranks[2],
            m_rankmap_lib_neighbor_ranks[3]);
  } else if (dim==4){
    fprintf(stderr, "rank %d: neighbors = %d %d %d %d\n",
            myrank,
            m_rankmap_lib_neighbor_ranks[0],
            m_rankmap_lib_neighbor_ranks[1],
            m_rankmap_lib_neighbor_ranks[2],
            m_rankmap_lib_neighbor_ranks[3],
            m_rankmap_lib_neighbor_ranks[4],
            m_rankmap_lib_neighbor_ranks[5],
            m_rankmap_lib_neighbor_ranks[6],
            m_rankmap_lib_neighbor_ranks[7]);
  }
#endif
  MPI_Barrier(MPI_COMM_WORLD);

  for(int dir=0; dir<dim; dir++){

    // recieve from the negative direction
    int srank=m_rankmap_lib_neighbor_ranks[2*dir];
    int rrank=m_rankmap_lib_neighbor_ranks[2*dir+1];
    int tag=2*dir+1;
    //fprintf(stderr, "rank %d, dir=%d: negative: srank=%d, rrank=%d, Sendrecv...\n", myrank, dir, srank, rrank);
    MPI_Sendrecv(sbuf, 6, MPI_UINT8_T, srank, tag,
                 rbuf, 6, MPI_UINT8_T, rrank, tag,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //fprintf(stderr, "rank %d, dir=%d: negative: srank=%d, rrank=%d, Sendrecv, done...\n", myrank, dir, srank, rrank);

    const volatile uint8_t *tmp=&m_rankmap_lib_neighbor_tofu[2*dir+1][0];
    if( tmp[0] !=rbuf[0] || tmp[1] !=rbuf[1] ||
        tmp[2] !=rbuf[2] || tmp[3] !=rbuf[3] ||
        tmp[4] !=rbuf[4] || tmp[5] !=rbuf[5] ){
      fprintf(stderr, "rank %d: bad tofu coordinate for dir=%d: calc= %d %d %d %d %d %d, recv=%d %d %d %d %d %d [from %d]\n",
              myrank, 2*dir+1,
              tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5],
              rbuf[0],rbuf[1],rbuf[2],rbuf[3],rbuf[4],rbuf[5],
              rrank);
      err=1;
    }

    // recieve from the positive direction
    srank=m_rankmap_lib_neighbor_ranks[2*dir+1];
    rrank=m_rankmap_lib_neighbor_ranks[2*dir];
    tag=2*dir;
    //fprintf(stderr, "rank %d, dir=%d: positive: srank=%d, rrank=%d, Sendrecv...\n", myrank, dir, srank, rrank);
    MPI_Sendrecv(sbuf, 6, MPI_UINT8_T, srank, tag,
                 rbuf, 6, MPI_UINT8_T, rrank, tag,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //fprintf(stderr, "rank %d, dir=%d: positive: srank=%d, rrank=%d, Sendrecv, done\n", myrank, dir, srank, rrank);

    tmp=&m_rankmap_lib_neighbor_tofu[2*dir][0];
    if( tmp[0] !=rbuf[0] || tmp[1] !=rbuf[1] ||
        tmp[2] !=rbuf[2] || tmp[3] !=rbuf[3] ||
        tmp[4] !=rbuf[4] || tmp[5] !=rbuf[5] ){
      fprintf(stderr, "rank %d: bad tofu coordinate for dir=%d: calc= %d %d %d %d %d %d, recv=%d %d %d %d %d %d [from %d]\n",
              myrank, 2*dir,
              tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5],
              rbuf[0],rbuf[1],rbuf[2],rbuf[3],rbuf[4],rbuf[5],
              rrank);
    }
  }

  return err;
}

//
// obtain rank id of the logical neighbors
//
void exchange_ranks(uint8_t *my_coords,
                    int *pos_ranks, int *neg_ranks,
                    uint8_t (*pos_coords)[6],  const int *pos_rank_in_node,
                    const int dim, const int proc_per_node){

  assert(dim<=4);
  utofu_stadd_t stadd[4];
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  int rc;
  //
  // crealte local vcq
  //
  int tmp=myrank % proc_per_node;
  int comp_id=tmp;  // num_cmp_ids is 8
  int cq_id=3*tmp;
  int tni_id=tmp;
  struct utofu_onesided_caps *onesided_caps;
  { // obtain available features
    int rc=utofu_query_onesided_caps(tni_id, &onesided_caps);
    if(rc != UTOFU_SUCCESS){
      fprintf(stderr, "rank %d: Failed at utofu_query_onesided_caps()\n", myrank);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
  }
  assert(comp_id < onesided_caps->num_cmp_ids);

  utofu_vcq_hdl_t local_vcq_hdl;
  //unsigned long int vcq_flag=UTOFU_VCQ_FLAG_THREAD_SAFE | UTOFU_VCQ_FLAG_EXCLUSIVE;
  //unsigned long int vcq_flag=UTOFU_VCQ_FLAG_EXCLUSIVE;
  unsigned long int vcq_flag=UTOFU_VCQ_FLAG_THREAD_SAFE;
  //unsigned long int vcq_flag=0UL;
  rc=utofu_create_vcq_with_cmp_id(tni_id, comp_id, vcq_flag, &local_vcq_hdl );
  if(UTOFU_SUCCESS != rc){
    fprintf(stderr, "rank %d: Failed at utofu_create_vcq_with_comp_id()! rc=%d\n", myrank, rc);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  //
  // check if the cq_id reproduces the id for local_vcq_hdl
  //
  check_cq_id(&local_vcq_hdl, my_coords, cq_id, comp_id, tni_id);

  //
  // regsiter the buffer to stadd: the local rank is written to the buffer
  //
  uint32_t *ranks_addr;
  size_t stag_alignement=onesided_caps->stag_address_alignment;
  int asize=stag_alignement/sizeof(uint32_t);
  posix_memalign((void**)&ranks_addr, stag_alignement, dim*stag_alignement);
  for(int dir=0; dir<dim; dir++){
    uint32_t *addr=ranks_addr + dir*asize;
    addr[0]=myrank;  // set the local rank
    unsigned int stag=dir;
    unsigned long flag_stadd=0;

    rc=utofu_reg_mem_with_stag(local_vcq_hdl, addr, 256, stag, flag_stadd, &stadd[dir]);
    if(UTOFU_SUCCESS != rc){
      fprintf(stderr, "rank %d: Failed at utofu_reg_mem_with_stag() for local send buffer! rc=%d\n", myrank, rc);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
#ifdef  _RANKMAP_DEBUG
    fprintf(stderr, "rank %d: dir=%d, stadd=%ul\n", myrank, dir, stadd[dir]);
#endif
  }
  // make sure that the stadd is registered in the remote
  MPI_Barrier(MPI_COMM_WORLD);

  ////////////////////////////////////////////////////////////////////////////
  // armw swap:
  //   the buffer is specivied with stadd
  //
  //   local  op_value --> remote buffer
  //   remote   buffer --> local  notice.rmt_val (the value BEFORE overwritten with op_value)
  //   remote    edata --> local  notice.edata
  //
  // swap with posivie logical direction:
  //   before: op_value, buffer = local rank
  //   after:   buffer         = rank in negative dirction
  //            notice.rmt_val = rank in positive dirction
  //   [edata: direction]


  //
  // prepare arwm communcation
  //
  void *desc;
  posix_memalign((void**)&desc, 8, onesided_caps->max_toq_desc_size*dim); //  64=max_toq_desc_size
  void *this_desc=desc;
  size_t desc_size=0;
  for(int dir=0; dir<dim; dir++){
    int tmp=pos_rank_in_node[dir] % proc_per_node;
    int comp_id=tmp;  // num_cmp_ids is 8
    int cq_id=3*tmp;  // this is a guess
    int tni_id=tmp;
    unsigned int stag=dir;
    utofu_vcq_id_t pos_vcq_id;
    utofu_stadd_t pos_stadd;

    // calculate remote_vcq_id
    rc=utofu_construct_vcq_id(pos_coords[dir], tni_id, cq_id, comp_id, &pos_vcq_id);
    if(UTOFU_SUCCESS != rc){
      const uint8_t *c=pos_coords[dir];
      fprintf(stderr, "rank %d: Failed at utofu_construct_vcq_id()! rc=%d, dir=%d, pos_coords=%d %d %d %d %d %d\n", myrank, rc,
              dir, c[0], c[1], c[2], c[3], c[4], c[5]);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
#ifdef  _RANKMAP_DEBUG
    fprintf(stderr, "rank %d: dir=%d, pos_vcq_id=%ul\n", myrank,
            dir, pos_vcq_id);
#endif

    // calculate remote_stadd
    rc=utofu_query_stadd(pos_vcq_id, stag, &pos_stadd);
    if(UTOFU_SUCCESS != rc){
      fprintf(stderr, "rank %d: Failed at utofu_querry_stadd()! rc=%d\n", myrank, rc);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
#ifdef  _RANKMAP_DEBUG
    fprintf(stderr, "rank %d: dir=%d, pos_stadd=%ul\n", myrank,
            dir, pos_stadd);
#endif

    uint32_t op_value=myrank;
    uint64_t edata=dir;
    const unsigned long int flag_armw4 = UTOFU_ONESIDED_FLAG_LOCAL_MRQ_NOTICE | UTOFU_ONESIDED_FLAG_STRONG_ORDER;

    size_t this_desc_size;
    rc=utofu_prepare_armw4(local_vcq_hdl, pos_vcq_id, UTOFU_ARMW_OP_SWAP, op_value, pos_stadd, edata, flag_armw4, this_desc, &this_desc_size);
    this_desc=(char*)this_desc+this_desc_size;
    desc_size+=this_desc_size;
  }


  //
  // excute arwm communication
  //
  void *cbdata;
  while(1){
    rc = utofu_post_toq(local_vcq_hdl, desc, desc_size, &cbdata);
    if(rc != UTOFU_ERR_BUSY) { break; }
    rc = utofu_poll_tcq(local_vcq_hdl, 0UL, &cbdata);
  }
  if(rc != UTOFU_SUCCESS){
    fprintf(stderr, "rank %d: Failed at utofu_post_toq()! rc=%d\n", myrank, rc);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  //
  // poll the data
  //
  for(int i=0; i<dim; i++){
    struct utofu_mrq_notice notice;
    //    fprintf(stderr, "rank %d: dir=%d, polling mrq...\n", myrank, i);
    do {
      rc = utofu_poll_mrq(local_vcq_hdl, 0UL, &notice);
    } while (rc == UTOFU_ERR_NOT_FOUND);
    if(rc != UTOFU_SUCCESS){
      fprintf(stderr,"rank %d: %s at %d : i=%d, utofu_poll_mrq, ERROR: %d\n",
              myrank, __FILE__, __LINE__, i, rc);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    } else {
      switch(notice.notice_type){
      case UTOFU_MRQ_TYPE_LCL_ARMW:
        break;
        //      case UTOFU_MRQ_TYPE_LCL_PUT:
        //      break;
      default:
        {
          const char *type_str;
          type_str="unkown";
          if(notice.notice_type == UTOFU_MRQ_TYPE_LCL_GET) {type_str="LCL_GET";}
          if(notice.notice_type == UTOFU_MRQ_TYPE_RMT_GET) {type_str="RMT_GET";}
          if(notice.notice_type == UTOFU_MRQ_TYPE_LCL_ARMW) {type_str="LCL_ARMW";}
          if(notice.notice_type == UTOFU_MRQ_TYPE_RMT_ARMW) {type_str="RMT_ARMW";}
          fprintf(stderr,"rank %d: %s at %d : utofu_poll_mrq, notice_type=%s\n",
                  myrank, __FILE__, __LINE__, type_str);
          MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
      }
    }
    //fprintf(stderr, "rank %d: dir=%d, polling mrq, done\n", myrank, i);

    uint64_t this_dir=notice.edata;
    uint32_t rank_pos=notice.rmt_value;
    pos_ranks[this_dir]=rank_pos;
  }
  // make sure all the buffers are updated
  MPI_Barrier(MPI_COMM_WORLD);
  if(myrank==0){
    //  fprintf(stderr, "exchange of rank id with neighbors, done.\n");
  }

  //
  // set the value obtained
  //
  for(int dir=0; dir<dim; dir++){
    neg_ranks[dir]=ranks_addr[dir*asize];
  }

  //
  // release the resources
  //
  for(int dir=0; dir<dim; dir++){
    utofu_dereg_mem(local_vcq_hdl, stadd[dir], 0);
  }
  free(ranks_addr);
  free(desc);
  utofu_free_vcq(local_vcq_hdl);

}


//
// obtain rank id of the logical neighbors
//
void exchange_ranks_with_mpi(int *pos_ranks, int *neg_ranks,
                             const int *rank_coord, const int *rank_size,
                             const int dim){

  assert(dim<=4);
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  int num_ranks;
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

  int *rank_list_send=(int*)malloc(sizeof(int)*num_ranks);
  int *rank_list     =(int*)malloc(sizeof(int)*num_ranks);
  if(!rank_list_send){
    fprintf(stderr, "rank %d: Failed at malloc for rank_list0.\n", myrank);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  if(!rank_list){
    fprintf(stderr, "rank %d: Failed at malloc for rank_list0.\n", myrank);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  memset(rank_list_send, 0, sizeof(int)*num_ranks);

  int myid=0;
  int stride=1;
  for(int i=0; i<dim; i++){
    myid+=stride*rank_coord[i];
    stride*=rank_size[i];
  }
  rank_list_send[myid]=myrank;

  MPI_Allreduce(rank_list_send, rank_list, num_ranks, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  for(int dir=0; dir<dim; dir++){
    int pos_id=0;
    int neg_id=0;
    int stride=1;
    for(int i=0; i<dim; i++){
      int pos_coord=rank_coord[i];
      int neg_coord=rank_coord[i];
      if(i==dir){
        pos_coord=(rank_coord[i]+1) % rank_size[i];
        neg_coord=(rank_coord[i]-1+rank_size[i]) % rank_size[i];
      }
      pos_id+=stride*pos_coord;
      neg_id+=stride*neg_coord;
      stride*=rank_size[i];
    }
    pos_ranks[dir]=rank_list[pos_id];
    neg_ranks[dir]=rank_list[neg_id];
  }

  free(rank_list_send);
  free(rank_list);
}


//
// utilities for rankmap and tni assignemnt
//
void rankmap_lib_get_rankmap(int *myrank_coord, int *neighbors, int *process_size) {
  for(int i=0; i < m_rankmap_lib_dim; i++){
    myrank_coord[i]=m_rankmap_lib_process_coord[i];
    neighbors[2*i]=m_rankmap_lib_neighbor_ranks[2*i];
    neighbors[2*i+1]=m_rankmap_lib_neighbor_ranks[2*i+1];
    process_size[i]=m_rankmap_lib_process_size[i];
  }
}



int rankmap_lib_set_rankmap4d() {

  static const int dim=4;
  static const int proc_per_node=4;
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  //  fprintf(stderr, "rank %d: rankmap_lib_set_rankmap4d is called\n", myrank);

  assert(proc_per_node==4);
  int err=0;

  m_rankmap_lib_myrank_in_node=myrank % proc_per_node;
  m_rankmap_lib_dim=dim;

  //
  // get tofu coordinate of this and logical neghboring nodes
  //
  uint8_t pos_coords[4][6];  // target Tofu coordinate: posivie neighbors
  int pos_rank_in_node[4];   // rank id in the node:    posivie neighbors
  uint8_t neg_coords[4][6];  // target Tofu coordinate: negative neighbors
  int neg_rank_in_node[4];   // rank id in the node:    negative neighbors
  uint8_t my_coords[6];      // Tofu coordinate of this rank
  int rank_coord[4];         // logical rank coordiante
  int rank_size[4];          // logical rank size
  int tni_list[8];

#ifdef _USE_RANKMAP
  int mapid=get_tofu_coord_and_tni(my_coords, rank_coord, rank_size, &pos_coords[0], pos_rank_in_node, &neg_coords[0], neg_rank_in_node, tni_list);
#else
  int mapid=-1;
#endif
  if(mapid<0){
    //fprintf(stderr, "rank %d: Failed at get_tofu_coord()! err=%d\n", myrank, err);
    if(myrank==0){
      fprintf(stderr, "WARNING: no rank map is found.\n");
    }
    get_tni_list_default(tni_list, myrank, rank_coord, rank_size, mapid);
    for(int dir2=0; dir2<8; dir2++){m_rankmap_lib_tni_list[dir2]=tni_list[dir2];}

    return mapid;
    //    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  if(myrank==0){
    printf("rdma_utofu_comlib: 4d rankmap, size = %d x %d x %d x %d\n", rank_size[0], rank_size[1], rank_size[2], rank_size[3]);
    fflush(stdout);
  }

#ifdef _RANKMAP_DEBUG
  fprintf(stderr, "rank %d: tofu coord = (%d %d %d %d %d %d):%d, +x:(%d %d %d %d %d %d):%d, +y:(%d %d %d %d %d %d):%d\n",
          myrank,
          my_coords[0], my_coords[1], my_coords[2], my_coords[3], my_coords[4], my_coords[5], tmp,
          pos_coords[0][0], pos_coords[0][1], pos_coords[0][2], pos_coords[0][3], pos_coords[0][4], pos_coords[0][5], rank_in_node[0],
          pos_coords[1][0], pos_coords[1][1], pos_coords[1][2], pos_coords[1][3], pos_coords[1][4], pos_coords[1][5], rank_in_node[1]);

  fprintf(stderr, "rank %d: tni_id=%d, cq_id=%d, comp_id=%d\n", myrank,
          tni_id, cq_id, comp_id);
#endif


  //
  // get rank id of the logical neghboring process
  //
  int pos_ranks[4]; // positive neighbor
  int neg_ranks[4]; // negative neighbor

  //  exchange_ranks(my_coords, pos_ranks, neg_ranks,
  //                 pos_coords, pos_rank_in_node,
  //                 dim, proc_per_node);

  exchange_ranks_with_mpi(pos_ranks, neg_ranks, rank_coord, rank_size, dim);


  //
  // set the result to the global variables
  //
  for(int i=0; i<6; i++){
    m_rankmap_lib_mytofu[i]=my_coords[i];
  }

  for(int i=0; i<8; i++){
    for(int j=0; j<6; j++) {
      m_rankmap_lib_neighbor_tofu[i][j]=99;
    }}
  for(int dir=0; dir<dim; dir++){
    m_rankmap_lib_process_coord[dir]=rank_coord[dir];
    m_rankmap_lib_process_size[dir]=rank_size[dir];
    m_rankmap_lib_neighbor_ranks[2*dir  ] = pos_ranks[dir];
    m_rankmap_lib_neighbor_ranks[2*dir+1] = neg_ranks[dir];
    m_rankmap_lib_neighbor_rank_in_node[2*dir  ] = pos_ranks[dir] % proc_per_node;
    m_rankmap_lib_neighbor_rank_in_node[2*dir+1] = neg_ranks[dir] % proc_per_node;
    m_rankmap_lib_tni_list[2*dir]=tni_list[2*dir];
    m_rankmap_lib_tni_list[2*dir+1]=tni_list[2*dir+1];
    for(int i=0; i<6; i++){
      m_rankmap_lib_neighbor_tofu[2*dir  ][i] = pos_coords[dir][i];
      m_rankmap_lib_neighbor_tofu[2*dir+1][i] = neg_coords[dir][i];
    }
    // bug? w/o this dumping, m_rankmap_lib_neighbor is not properly updated;
    // or violating strict aliasing rule (throught accesing slice of array)?
    // [ to aviod this bug, m_rankmap_lib_neighbor is changed to volatile ]
    //    uint8_t *p=m_rankmap_lib_neighbor_tofu[2*dir  ];
    //    uint8_t *n=m_rankmap_lib_neighbor_tofu[2*dir+1];
    //    fprintf(stderr, "fuga: I am %d, dir=%d: positive_ncoods = %d %d %d %d %d, negative_ncoords = %d %d %d %d %d %d\n",
    //        myrank, dir, p[0],p[1],p[2],p[3],p[4],p[5], n[0],n[1],n[2],n[3],n[4],n[5]);

  }

#ifdef DUMP_RANKMAP
  for(int dir=0; dir <dim; dir++){
    volatile uint8_t *p=m_rankmap_lib_neighbor_tofu[2*dir  ];
    volatile uint8_t *n=m_rankmap_lib_neighbor_tofu[2*dir+1];
    volatile uint8_t *m=m_rankmap_lib_mytofu;
    fprintf(stderr, "  rankmap: rank=%d [%d %d %d %d %d %d], dir=%d: positive_ncoords = %d %d %d %d %d %d (rank %d), nagative_ncoords = %d %d %d %d %d %d (rank %d)\n",
            myrank, m[0], m[1], m[2], m[3], m[4], m[5],
            dir,
            p[0],p[1],p[2],p[3],p[4],p[5], pos_ranks[dir], 
            n[0],n[1],n[2],n[3],n[4],n[5], neg_ranks[dir]);

  }
#endif

  if(myrank==0){
    fprintf(stderr, "checking neighbor ranks...\n");
  }
  MPI_Barrier(MPI_COMM_WORLD);
  err=check_neighbors();
  if(myrank==0){
    fprintf(stderr, "checking neighbor ranks, done: err=%d\n", err);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  if(err){
    fprintf(stderr, "rank %d: error at check_neighbors()\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  if(myrank==0){
    printf("rdma_utofu_comlib: 4d rankmap is ready: mapid=%d\n", mapid);
    fflush(stdout);
  }

  return mapid;
}

int rankmap_lib_get_tni_list(int *tni_list, const int *flag){
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  if(myrank==0){
    fprintf(stderr, "rankmap_lib: calling get_tni_list for dim=%d, flag=%d\n", m_rankmap_lib_dim, *flag);
  }
  for(int dir2=0; dir2<2*m_rankmap_lib_dim; dir2++){
    tni_list[dir2] = m_rankmap_lib_tni_list[dir2];
  }

  return 0;

}

void rankmap_lib_dump_rankmap(){
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  volatile int *rank_coord=m_rankmap_lib_process_coord;
  volatile int *nranks= m_rankmap_lib_neighbor_ranks;

  if(m_rankmap_lib_dim==2){
    fprintf(stderr, "rank %d: (%d,%d)  +x:%d -x:%d +y:%d -y:%d\n",myrank,
            rank_coord[0],rank_coord[1], nranks[0], nranks[1], nranks[2], nranks[3]);
  } else if(m_rankmap_lib_dim==4){
    fprintf(stderr, "rank %d: (%d,%d,%d,%d)  +x:%d -x:%d +y:%d -y:%d +z:%d -z:%d +t:%d -t:%d\n",myrank,
            rank_coord[0],rank_coord[1],rank_coord[2],rank_coord[3],
            nranks[0],nranks[1], nranks[2],nranks[3], nranks[4],nranks[5], nranks[6],nranks[7]);
  }

}

//
