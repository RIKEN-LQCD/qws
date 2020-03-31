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
  rdma commucation with double buffering:
      wrapper functions for rdma communication
  Issaku Kanamori kanamori@hiroshima-u.ac.jp
  2018 Dec. 17, 2019 Feb. 6

  Copying and distribution of this file, with or without modification,
  are permitted in any medium without royalty provided the copyright
  notice and this notice are preserved.  This file is offered as-is,
  without any warranty.


  Usage:


void rdma_comlib_2buf::comlib_init()
void rdma_comlib_2buf::comlib_finalize()
  initialize/finalize the rdma communication
  NOTE: Re-initialization does not work in this FJMPI RDMA implementation.

member functions:

void init(const int tni_id, const int dst_rank, const int rcv_rank, const size_t size)
  allocate the buffer
  tni_id:  specifies the Tofu Netowork Interface (TNI), [0...5]
  dst_rank: "target" rank for the 1-to-1 communication
  rcv_rank: rank of which "target" is this rank
  size:    buffer size in byte

void finalize()
  deallocate the buffer

void* sbuff()
  retunrs pointer to the send buffer

void* rbuff()
  recieve buffer with the current parity

void irecv()
  asynchronous RDMA get the data from reciev buffer with the current parity

void isendrecv()
   asynchronous RDMA send/receive with the current parity

void isend_check()
   check and wait for RDMA isend finishes

void irecv_check()
   check and wait for RDMA irecieve finishes

void irecv_ok()
   set ok status to the receive buffer with the current parity

void reset_comm()
   clear the has_started flag set by isendrecv(), in order to supress *_check()

void flip_parity()
   flip the parity for the buffering

int get_parity()
   return the parity

static void swap_vcq_for_sending(rdma_comlib_2buf &buf1, rdma_comlib_2buf &buf2)
   swap the interanl vcqs used in sending.
   suppose that before the swap:
     vcq1(= vcq in buf1) sends to, e.g., +A and recieves from -A.
     vcq2(= vcq in buf2) sends to, e.g., -A and recieves from +A.
   after the swap
     vcq1 sends to -A and recieves from -A   (now communicates with only -A)
     vcq2 sends to +A and recieves from +A   (now communicates with only +A)

void excahge_rdma_comlib_targets(rdma_comlib_2buf &buf1, rdma_comlib_2buf &buf2)
   [deprecated] swap the targets to send.  If buf1 and buf2 send data to the opposite directions
    to each other, out-going and in-coming trafic to/from the same process become to share the vcq (and TNI)

// interface to rankmap
static int set_rankmap(int dim)
  prepare rankmap and its internal id.
  (in) dim: dimension of the system, [2,4]

static int get_rank_map(int *rank_coord, int *neighbor_ranks, int *rank_size)
  obtain the prepared the rankmap
  (out) rank_coord: 2 or 4 dim array to specify the logical rank cooridanate
  (out) neighbor_ranks: 4 or 8 dim array to spedcify the logical neghbors.  Directions are [+x, -x, +y, -y,...]
  (out) rank_size: logical size of 2 or 4 dim ranks.

 */
#ifndef RDMA_COMLIB_2BUF_H
#define RDMA_COMLIB_2BUF_H
#include <stdlib.h>
#include <mpi.h>

//extern "C" {
//#include "rdma_utofu_comlib.h"
//}

#define RDMA_2BUF_USE_SHARED_COMMUNICATOR

#ifdef _UTOFU_RDMA
//extern "C" {
//#include "rdma_utofu_comlib.h"
//}
// data type defined in the rdma_utofu_comlib.h
typedef struct{
  volatile int *sbuff;
  volatile int *rbuff;
  //  char dummy[96];  // the size of dummy is environment dependent
  char dummy[128];  // the size of dummy is environment dependent
} rdma_data_t;


// interface to call functions in rdma_utofu_comlib.c
extern "C" {
  void rdma_comlib_init(void);
  void rdma_comlib_finalize(void);

  void rdma_comlib_new_ext(rdma_data_t *id, const int *nic_id, const int *dst_rank, const int *rcv_rank, const size_t *size, const int *has_ext);
  void rdma_comlib_delete_ext(rdma_data_t *id, const int *has_ext);
  void rdma_comlib_new_duplicate(rdma_data_t *id, const rdma_data_t *id_org, const size_t *size, const int *has_external); // new communicataor with the same vcq
  void rdma_comlib_exchange_targets(rdma_data_t *id1, rdma_data_t *id2);
  void rdma_comlib_swap_vcq_for_sending(rdma_data_t *id1, rdma_data_t *id2);

  void rdma_comlib_isendrecv(rdma_data_t *id);
  void rdma_comlib_isend_check(rdma_data_t *id);
  void get_rdma_comlib_isend_finish_tag(int *nicid, int *tag);

  void rdma_comlib_irecv(rdma_data_t *id);
  void rdma_comlib_irecv_check(rdma_data_t *id);
  void rdma_comlib_irecv_ok(rdma_data_t *id);
  void rdma_comlib_clear_mrq(rdma_data_t *id);

  int rdma_comlib_get_ssize(const rdma_data_t *id);

  #include "rankmap_lib.h"
}

namespace {

  void rdma_comlib_new_2buf_wrapper(rdma_data_t *id, const int nic_id, const int dst_rank, const int rcv_rank, const size_t size, const int has_ext){
#ifdef RDMA_2BUF_USE_SHARED_COMMUNICATOR
  rdma_comlib_new_ext(&id[0], &nic_id, &dst_rank, &rcv_rank, &size, &has_ext);
  rdma_comlib_new_duplicate(&id[1], &id[0], &size, &has_ext);
#else
  rdma_comlib_new_ext(&id[0], &nic_id, &dst_rank, &rcv_rank, &size, &has_ext);
  rdma_comlib_new_ext(&id[1], &nic_id, &dst_rank, &rcv_rank, &size, &has_ext);
#endif
  }

}

#else
//extern "C" {
//#include "rdma_comlib.h"
//}
typedef struct{
  volatile int *sbuff;
  volatile int *rbuff;
  char dummy[104];  // the size of dummy is environment dependent
} rdma_data_t;

// interface to call functions in rdma_comlib.c
extern "C" {
  void rdma_comlib_init(void);
  void rdma_comlib_finalize(void);

  void rdma_comlib_new_ext(rdma_data_t *id, const int *nic_id, const int *dst_rank, const size_t *size, const int *has_ext);
  void rdma_comlib_delete_ext(rdma_data_t *id, const int *has_ext);

  void rdma_comlib_isendrecv(rdma_data_t *id);
  void rdma_comlib_isend_check(rdma_data_t *id);
  void get_rdma_comlib_isend_finish_tag(int *nicid, int *tag);

  void rdma_comlib_irecv(rdma_data_t *id);
  void rdma_comlib_irecv_check(rdma_data_t *id);
  void rdma_comlib_irecv_ok(rdma_data_t *id);

  int rdma_comlib_get_ssize(const rdma_data_t *id);
}

namespace {
  void rdma_comlib_clear_mrq(rdma_data_t* id) {}
  void rdma_comlib_new_ext_wrapper(rdma_data_t *id, const int nic_id, const int dst_rank, const int rcv_rank, const size_t size, const int has_ext){
    rdma_comlib_new_ext(id, &nic_id, &dst_rank, &size, &has_ext);
  }

}
#endif

class rdma_comlib_2buf{

 public:

 rdma_comlib_2buf(): m_parity(0), m_is_allocated(false) { }
  ~rdma_comlib_2buf() {
    finalize();  // delete the allocated buffers
  }

  // wrapper to rdma_comlib: init/finalize
  static void comlib_init() {
    rdma_comlib_init();
  }

  static void comlib_finalize(void) {
    rdma_comlib_finalize();
  }

  // todo: separate rank map/tni assignment to indepedent class
  static int set_rankmap(int dim){
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if(dim==2){
      m_rankmap_id=rankmap_lib_set_rankmap2d();
    } else if(dim==4) {
      m_rankmap_id=rankmap_lib_set_rankmap4d();
    } else {
      fprintf(stderr, "err in making rankmap: unkown dim = %d\n",dim);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    if(m_rankmap_id<0){
      if(myrank==0){
        fprintf(stderr, "WARNING: err in making rankmap, using the default map.\n");
      }
    }
    return m_rankmap_id;
  }

  static int get_rank_map(int *rank_coord, int *neighbor_ranks,
                          int *rank_size){
    //fprintf(stderr, "get_neighbor_ranks: m_rankmap_id=%d\n",m_rankmap_id);
    if(m_rankmap_id>=0) {
      rankmap_lib_get_rankmap(rank_coord, neighbor_ranks, rank_size);
    }
    return m_rankmap_id;
  }


  static void get_tni_ids(int *tni_ids, int m_rankmap_id){
    rankmap_lib_get_tni_list(tni_ids, &m_rankmap_id);
  }

  static void swap_vcq_for_sending(rdma_comlib_2buf &buf1, rdma_comlib_2buf &buf2){
    buf1.swap_vcq_for_sending(buf2);
  };

  // wapper to rdma_comlib:  new/delete
  void init(const int tni_id, const int dst_rank, const int rcv_rank, const size_t size);
  void finalize();

  void exchange_targets(rdma_comlib_2buf &buf2){
    rdma_comlib_exchange_targets(&m_rdma_data[0], &(buf2.m_rdma_data[0]));
    rdma_comlib_exchange_targets(&m_rdma_data[1], &(buf2.m_rdma_data[1]));
  }

  void swap_vcq_for_sending(rdma_comlib_2buf &buf2){
    rdma_comlib_swap_vcq_for_sending(&m_rdma_data[0], &(buf2.m_rdma_data[0]));
    rdma_comlib_swap_vcq_for_sending(&m_rdma_data[1], &(buf2.m_rdma_data[1]));
  }

  // accesser to the buffers
  void* sbuff() const
  {
    return (void*) m_rdma_data[m_parity].sbuff;
  }
  void* rbuff() const
  {
    return (void*) m_rdma_data[m_parity].rbuff;
  }

  // wrapper to rdma_comlib, cont'd
  void irecv()
  {
    rdma_comlib_irecv( &m_rdma_data[m_parity] );
  }

  void isendrecv()
  {
    rdma_comlib_isendrecv( &m_rdma_data[m_parity] );
    m_has_started[m_parity] = true;
  }

  void isend_check()
  {
    if(m_has_started[m_parity]){
      rdma_comlib_isend_check( &m_rdma_data[m_parity] );
    }
  }

  void irecv_check()
  {
    // poll and clean mrq for the previous communication
    //   this is needed to avoid flooding the mrq
    rdma_comlib_clear_mrq( &m_rdma_data[1-m_parity] );
    // wait for the recv buffer is ready
    rdma_comlib_irecv_check( &m_rdma_data[m_parity] );
  }

  void irecv_ok()
  {
    rdma_comlib_irecv_ok( &m_rdma_data[m_parity] );
  }

  void reset_comm()
  {
    m_has_started[m_parity]=false;
  }

  // parity for the double buffering
  void flip_parity() {
    m_parity = 1-m_parity;
  }

  int get_parity() const {
    return m_parity;
  }

 private:
  volatile int m_parity;
  volatile bool m_has_started[2];
  bool m_is_allocated;
  rdma_data_t m_rdma_data[2];

 private:
  volatile static int m_rankmap_id;

};


#endif


