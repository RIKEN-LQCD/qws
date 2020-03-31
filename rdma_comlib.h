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
//  Copying and distribution of this file, with or without modification,
//  are permitted in any medium without royalty provided the copyright
//  notice and this notice are preserved.  This file is offered as-is,
//  without any warranty.
//
#ifndef _RDMA_COMLIB_H
#define _RDMA_COMLIB_H

#include <stdio.h>

typedef struct {
  volatile int *sbuff;           // pointer to int array (additional +1 component is required to recv check)
  volatile int *rbuff;           // pointer to int array (additional +1 component is required to recv check)
  size_t length;                 // size of array (in byte data part)
  size_t data_size;              // size in byte (data + additilanl commonent)
  size_t data_len_int;
  int memid_sbuff;               // RDMA MEM ID for local send buffer
  int memid_rbuff;               // RDMA MEM ID for local receive buffer
  int local_rank;                // local MPI rank
  int target_rank;               // target MPI rank (send to)
  uint64_t addr_local_sbuff;     // RDMA address of local send buffer.
  uint64_t addr_target_rbuff;    // RDMA address of target receiv buffer.
  int nic_id;                    // nic id
  int put_options;               // RDMA put options (nic id)
  int get_options;               // RDMA get options (nic id)
  int poll_options;              // RDMA put_poll options (nic id)
  int tag;                       // message tag
  int dummy;
  struct FJMPI_Rdma_cq  cq;
  FILE *fp;
} rdma_comlib_data;



void rdma_comlib_init(void);
void rdma_comlib_finalize(void);
void rdma_comlib_new(rdma_comlib_data *id, const int nic_id, const int dst_rank, const size_t *size);
void rdma_comlib_new_ext(rdma_comlib_data *id, const int nic_id, const int dst_rank, const size_t *size, const int has_external);
void rdma_comlib_delete(rdma_comlib_data *id);
void rdma_comlib_delete_ext(rdma_comlib_data *id, const int has_external);

void rdma_comlib_isendrecv(rdma_comlib_data *id);
void rdma_comlib_isend_check(rdma_comlib_data *id);
void get_rdma_comlib_isend_finish_tag(int *nicid, int *tag);

void rdma_comlib_irecv(rdma_comlib_data *id);
void rdma_comlib_irecv_check(rdma_comlib_data *id);
void rdma_comlib_irecv_ok(rdma_comlib_data *id);

int rdma_comlib_get_ssize(const rdma_comlib_data *id);


#endif
