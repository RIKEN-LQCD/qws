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
////////////////////////////////////////////////////////////////////////////////
//
//  The original version is by Ken-Ichi Ishikawa <ishikawa@theo.phys.sci.hiroshima-u.ac.jp>
//
//  2018 Dec. 17:
//    modification of the original version
//    allows the user allocate the buffers in the outside
//    Issaku Kanamori <kanamori@hiroshima-u.ac.jp>
//
//  2019 Nov. 28:
//    UTOFU version.  Removed interfaces which are not used anymore.
//
//  2020 Jan. 9:
//    added rdma_comlib_new_duplicate
//
//  Copying and distribution of this file, with or without modification,
//  are permitted in any medium without royalty provided the copyright
//  notice and this notice are preserved.  This file is offered as-is,
//  without any warranty.
//
#ifndef _RDMA_UTOFU_COMLIB_H
#define _RDMA_UTOFU_COMLIB_H

#include <stdio.h>
#include <stdint.h>

typedef struct {
  uintptr_t local_vcq_hdl;       // handle to the local virtual control queue [utofu_vcq_hdl_t]
  uint64_t target_vcq_id;        // id of the hadle to the remote virtual control queue [utofu_vcq_id_t]
  int local_rank;                // local MPI rank
  int target_rank;               // target MPI rank (send to)
  int from_rank;                 // MPI rank from which data comes
  int tni_id;                    // TNI ID
  int remaining_remote_mrq;      // number of mrq to poll (0 or 1)
  int dummy;                     // padding
  int *ref_count;                // reference counter
} rdma_comlib_communicator;

typedef struct {
  volatile int *sbuff;           // pointer to int array (additional +1 component is required to recv check)
  volatile int *rbuff;           // pointer to int array (additional +1 component is required to recv check)
  size_t length;                 // size of array (in byte data part)
  size_t data_size;              // size in byte (data + additilanl commponent)
  size_t data_len_int;           // the same as data size but in sizeof(int)

  rdma_comlib_communicator comm;
  uint64_t sending_stadd;         // RDMA address (stadd) of local send buffer. [utofu_stadd_t]
  uint64_t local_stadd_rbuff;     // RDMA address (stadd) of local reciev buffer. [utofu_stadd_t]
  uint64_t target_stadd_rbuff;    // RDMA address (stadd) of target receive buffer. [utofu_stadd_t]
  uintptr_t sending_vcq_hdl;
  uint64_t local_stadd_sbuff;
  int tag;                        // message tag
  int last_component;             // checking the data arrival
  FILE *fp;
} rdma_comlib_data;



void rdma_comlib_init(void);
void rdma_comlib_finalize(void);
void rdma_comlib_new(rdma_comlib_data *id, const int *tni_id, const int *dst_rank, const int *rcv_rank, const size_t *size);
void rdma_comlib_new_ext(rdma_comlib_data *id, const int *tni_id, const int *dst_rank, const int *rcv_rank, const size_t *size, const int *has_external);
void rdma_comlib_new_duplicate(rdma_comlib_data *id, const rdma_comlib_data *id_org, const size_t *size, const int *has_external); // new communicataor with the same vcq
void rdma_comlib_exchange_targets(rdma_comlib_data *id1, rdma_comlib_data *id2);

void rdma_comlib_delete(rdma_comlib_data *id);
void rdma_comlib_delete_ext(rdma_comlib_data *id, const int *has_external);

void rdma_comlib_isendrecv(rdma_comlib_data *id);
void rdma_comlib_clear_mrq(rdma_comlib_data *id);
void rdma_comlib_isend_check(rdma_comlib_data *id);
//void get_rdma_comlib_isend_finish_tag(int *tni_id, int *tag);

void rdma_comlib_irecv(rdma_comlib_data *id);
void rdma_comlib_irecv_check(rdma_comlib_data *id);
void rdma_comlib_irecv_ok(rdma_comlib_data *id);

//int rdma_comlib_get_ssize(const rdma_comlib_data *id);


#endif
