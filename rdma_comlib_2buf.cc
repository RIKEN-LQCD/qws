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
      wrapper for functions rdma communication
  Issaku Kanamori kanamori@hiroshima-u.ac.jp
  2018 Dec. 17, 2019 Feb. 6

  UTOFU version
  Issaku Kanamori kanamori-i@riken.jp
  2019 Nov. 29

  Copying and distribution of this file, with or without modification,
  are permitted in any medium without royalty provided the copyright
  notice and this notice are preserved.  This file is offered as-is,
  without any warranty.

 */


#ifdef _UTOFU_RDMA

#include "qws.h"
#include "rdma_comlib_2buf.h"

volatile int rdma_comlib_2buf::m_rankmap_id=-999;

// initializatoin of the buffer
void rdma_comlib_2buf::init(const int tni_id, const int dst_rank, const int rcv_rank, const size_t size)
{

  if(m_is_allocated){  // free the buffer if allocated
    finalize();
  }

  //size_t length=size;                // data length in byte
  size_t data_size = size+(size_t)4;   // additional last 4-byte for flag
  const size_t alignment=CLS;          // cf. utofu_onesided_caps: stag_address_alignment  = 256
  size_t data_size_aligned=alignment*((data_size-1)/alignment+1);
  size_t length_aligned=data_size_aligned-(size_t)4;

  // send buffer: shared with both parities
  //int *sbuff = (int *)malloc(data_size);
  int *sbuff;
  posix_memalign((void**)&sbuff, alignment, data_size_aligned);
  m_rdma_data[0].sbuff=sbuff;
  m_rdma_data[1].sbuff=sbuff;

  // reciev buffer: each parity has its own buffer
  //  m_rdma_data[0].rbuff = (int *)malloc(data_size);
  //  m_rdma_data[1].rbuff = (int *)malloc(data_size);
  int *rbuff0, *rbuff1;
  posix_memalign((void**)&rbuff0, alignment, data_size_aligned);
  posix_memalign((void**)&rbuff1, alignment, data_size_aligned);
  m_rdma_data[0].rbuff = rbuff0;
  m_rdma_data[1].rbuff = rbuff1;

  // initialize the other data for th rdma comm.
  int use_external=1;
#ifdef _UTOFU_RDMA
  rdma_comlib_new_2buf_wrapper(m_rdma_data, tni_id, dst_rank, rcv_rank, length_aligned, use_external);
#else
  rdma_comlib_new_ext_wrapper(&m_rdma_data[0], tni_id, dst_rank, length_aligned, use_external);
  rdma_comlib_new_ext_wrapper(&m_rdma_data[1], tni_id, dst_rank, length_aligned, use_external);
#endif
  // reset the parity
  m_parity=0;

  // flag to check if the commucation has started
  m_has_started[0]=false;
  m_has_started[1]=false;

// set the flag
  m_is_allocated=true;

#ifdef RDMA_2BUF_USE_SHARED_COMMUNICATOR
  int myrank=0;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  if(myrank==0){
    printf("rdma_comlib_2buf: using shrared communicator for both buffers\n");
  }
#endif


}

// cleaning up the buffer
void rdma_comlib_2buf::finalize()
{
  if(!m_is_allocated){ // check if the buffer is allocated
    return;
  }

  int use_external=1;
  rdma_comlib_delete_ext(&m_rdma_data[0], &use_external);
  rdma_comlib_delete_ext(&m_rdma_data[1], &use_external);

  // send buffer: shared with both parities
  free((void*)m_rdma_data[0].sbuff);
  m_rdma_data[0].sbuff = NULL;
  m_rdma_data[1].sbuff = NULL;

  // reciev buffer: each parity has its own buffer
  free((void*)m_rdma_data[0].rbuff);
  free((void*)m_rdma_data[1].rbuff);
  m_rdma_data[0].rbuff = NULL;
  m_rdma_data[1].rbuff = NULL;

  // reset the parity
  m_parity=0;

  // clear the flag
  m_is_allocated=false;
}


#endif
