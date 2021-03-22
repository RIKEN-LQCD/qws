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

#ifdef _UTOFU_RDMA
////////////////////////////////////////////////////////////////////////////////
// RDMA comlib using UTOFU for FUGAKU
//   it assumes the FUGAKU version of strong ordering in the communication
//
// rdma_comlib_data {
//
//  volatile int *sbuff;      // pointer to int array used for send buffer used by user.
//  volatile int *rbuff;      // pointer to int array used for receive buffer used by user.
//
//  other comoponents should not be modified by user.
//
// }
//
//void rdma_comlib_init(void);
//  this initialize RDMA functionality. This should be called after MPI_Init
//  before using RDMA functions.
//
//void rdma_comlib_finalize(void);
//  this terminate RDMA functionality. This is called before MPI_Finalize.
//
//void rdma_comlib_new_ext(rdma_comlib_data *id, const int tni_id, const int dst_rank, const int rcv_rank, const size_t *size, const int has_external);
//  Create vcq, and create (if has_external==0)/register the buffers for RDMA 1-to-1 communication.
//    - The date in id->sbuff will be send to dst_rank.
//    - The date comes from rcv_rank will be stored to id->rbuff.
//    - The communication goes throough via TNI with tni_id.
//    - The unit of the data size (*size) is in byte (sizeof()).
//  The data buffers for send/receive is located in *id. If has_external==0, they are allocated
//  in this comlib_new_ext.  Otherwize assumes they are already allocated in advace.
//  The data to be sent should be stored in id->sbuff (pointer to int).
//  The data received will be be stored in id->rbuff (pointer to int).
//  All send/receive/check for the data communication are applied to the data structure *id.
//
//  rdma_comlib_data *id : pointer to rdma comlib data structure.
//  const int tni_id     : specify TNI id [0...5] for this communication.
//  const int dst_rank   : specify send destination MPI RANK
//  const int rcv_rank   : specify MPI RANK from which the data comes
//  const size_t *size   : pointer to data size in byte.
//
//void rdma_comlib_new(rdma_comlib_data *id, const int *tni_id, const int *dst_rank, const int *rcv_rank, const size_t *size);
//  the same as rdma_comlib_new_ext with has_external=0 so that it allocate the buffer memory
//
//void rdma_comlib_new_duplicate(rdma_comlib_data *id, const rdma_comlib_data *id_org, const size_t *size, const int *has_external);
//  using vcq information in id_org, create (if has_external==0)/register the buffers for RDMA 1-to-1 communication.
//
//void rdma_comlib_delete_ext(rdma_comlib_data *id, const int *has_external);
//  Delete rdma comlib data structure. If has_external==0, the communication buffer stored
//
//void rdma_comlib_delete(rdma_comlib_data *id);
//  the same as rdma_comlib_delete_ext with has_external=0; it always dellocate the buffer
//
//void rdma_comlib_isendrecv(rdma_comlib_data *id);
//  Start asynchronous send and recieve by rdma put communication. This is non-blocking.
//  The data in id->sbuff will be send to id->rbuff in id->taget_rank.
//  After this function:
//   -the user should not modefy the contents of id->sbuff
//    before rdma_comlib_isend_check.
//   -(in the target rank) the user should not use the contents
//    of id->rbuff before rdma_comlib_irecv_check.
//
//void rdma_comlib_isend_check(rdma_comlib_data *id);
//  Check wheather ashynchorouns send finishs.
//  After calling this function user can touch id->sbuff.
//
//void rdma_comlib_irecv_check(rdma_comlib_data *id);
//  Check wheather ashynchorouns receive finishs.
//  After calling this function user read data from id->rbuff.
//
//void rdma_comlib_irecv_ok(rdma_comlib_data *id);
//  Clear recieved flag for the receive buffer. User should call this function after
//  fnishing reading data form id->rbuff in order to check the new date has arrived
//  in this buffer.
//
//void rdma_comlib_swap_vcq_for_sending(rdma_comlib_data *id1, rdma_comlib_data *id2);
//  swap vcq used in sending
//   suppose that before the swap:
//     vcq1(= vcq in id1) sends to, e.g., +A and recieves from -A.
//     vcq2(= vcq in id2) sends to, e.g., -A and recieves from +A.
//   after the swap
//     vcq1 sends to -A and recieves from -A   (now communicates with only -A)
//     vcq2 sends to +A and recieves from +A   (now communicates with only +A)

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
//  Copying and distribution of this file, with or without modification,
//  are permitted in any medium without royalty provided the copyright
//  notice and this notice are preserved.  This file is offered as-is,
//  without any warranty.
//

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <mpi-ext.h>
#include <utofu.h>

//#define  _RDMA_DEBUG
//#define RDMA_NO_REMOTE_MRQ_POLLING
#define RDMA_USE_CACHE_INJECTION

const int MAXNUM_TNIID = 6;             //  Maximum nuber of TNI ID
const int MAX_RDMA_DATASIZE = 16777212; //  RDMA put/get MAX data size in byte
const int MAXNUM_TAG   = 255;            //  RDMA put/get MAX number of message tag
const int TAG_OFFSET   = 128;            //  offset for the message tag
const int MAX_LAST_COMPONENT  = 256;  //  maximum value of the watchdog data for receive data polling.

static int m_rdma_comlib_is_initialized = 0;
static int m_rdma_comlib_myrank = 0; // MPI LOCAL RANK is stored.
static int m_rdma_comlib_tag = 128;    // counts RDMA message tag. [0..127] + TAG_OFFSET

#include "rdma_utofu_comlib.h"


int rdma_comlib_get_ssize(const rdma_comlib_data *id)
{
   return (*id).length;
}

void rdma_comlib_init(void)
{
//
// Initialize UTOFU communication
//
// assumes MPI is already initialized
  if (0 == m_rdma_comlib_is_initialized) {

    MPI_Comm_rank(MPI_COMM_WORLD,&m_rdma_comlib_myrank);

    if ( 0 == m_rdma_comlib_myrank ) {
      printf("%%%% UTOFU interface for one sided communcation  (RDMA) is used.\n");
      printf("sizeof(rdma_comlib_data) = %d, sizeof(int*) = %d\n", sizeof(rdma_comlib_data), sizeof(int*));
#ifdef RDMA_NO_REMOTE_MRQ_POLLING
      printf("RDMA_NO_REMOTE_MRQ_POLLING is defined.\n");
#endif
      fflush(stdout);
    }
    m_rdma_comlib_is_initialized = 1;

    // check if we can use full TNIs
    size_t num_tnis;          // number of available TNIs
    utofu_tni_id_t *tni_ids;  // array of TNI IDs
    if(UTOFU_SUCCESS != utofu_get_onesided_tnis(&tni_ids, &num_tnis)){
      fprintf(stderr, "rank %d: Failed at utofu_get_onesided_tnis()!\n", m_rdma_comlib_myrank);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    { // obtain available features and dump them
      struct utofu_onesided_caps *onesided_caps;
      int rc=utofu_query_onesided_caps(*tni_ids, &onesided_caps);
      if(rc != UTOFU_SUCCESS){
        fprintf(stderr, "rank %d: Failed at utofu_query_onesided_caps()\n", m_rdma_comlib_myrank);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
      }
#ifdef _RDMA_DEBUG
      // dump caps...
      if(m_rdma_comlib_myrank==0){
        int rank=m_rdma_comlib_myrank;
        fprintf(stderr, "rank %d: utofu_onesided_caps: flags = %ul\n", rank, onesided_caps->flags);
        fprintf(stderr, "rank %d: utofu_onesided_caps: arwm_ops  = %ul\n", rank, onesided_caps->armw_ops);
        fprintf(stderr, "rank %d: utofu_onesided_caps: num_cmp_ids  = %ud\n", rank, onesided_caps->num_cmp_ids);
        fprintf(stderr, "rank %d: utofu_onesided_caps: num_reserved_stags  = %ud\n", rank, onesided_caps->num_reserved_stags);
        fprintf(stderr, "rank %d: utofu_onesided_caps: cache_line_size  = %ul\n", rank, onesided_caps->cache_line_size);
        fprintf(stderr, "rank %d: utofu_onesided_caps: stag_address_alignment  = %ul\n", rank, onesided_caps->stag_address_alignment);
        fprintf(stderr, "rank %d: utofu_onesided_caps: max_toq_desc_size  = %ul\n", rank, onesided_caps->max_toq_desc_size);
        fprintf(stderr, "rank %d: utofu_onesided_caps: max_putget_size  = %ul\n", rank, onesided_caps->max_putget_size);
        fprintf(stderr, "rank %d: utofu_onesided_caps: max_piggyback_size  = %ul\n", rank, onesided_caps->max_piggyback_size);
        fprintf(stderr, "rank %d: utofu_onesided_caps: max_edata_size  = %ul\n", rank, onesided_caps->max_edata_size);
        fprintf(stderr, "rank %d: utofu_onesided_caps: max_mtu  = %ul\n", rank, onesided_caps->max_mtu);
        fprintf(stderr, "rank %d: utofu_onesided_caps: max_gap  = %ul\n", rank, onesided_caps->max_gap);
      }
#endif
    }

    if(num_tnis < MAXNUM_TNIID){
      fprintf(stderr, "rank %d: only %d TNIs are available, aborting...\n", m_rdma_comlib_myrank, num_tnis);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    free(tni_ids);
  }
}

void rdma_comlib_finalize(void) {
  if (1 == m_rdma_comlib_is_initialized) {
    // nothing to do
  }
}


void rdma_comlib_communicator_new(rdma_comlib_communicator *comm, const int *tni_id, const int *dst_rank, const int *rcv_rank){

  int err=0;
  if ( (MAXNUM_TNIID <= *tni_id) || ( *tni_id < 0 ) ) {
    fprintf(stderr,"TNI ID should be greater than -1 and less than MAXNUM_TNIID. 0 <= tni_id = %d < %d = MAXNUM_TNIID\n",tni_id,MAXNUM_TNIID);
    err=1;
  }
  if(err){
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  //
  // allocate reference counter
  //
  comm->ref_count=malloc(sizeof(int));

  //
  // set values
  //
  comm->tni_id=*tni_id;
  comm->target_rank=*dst_rank;
  comm->from_rank=*rcv_rank;
  comm->local_rank=m_rdma_comlib_myrank;
  *comm->ref_count=0;  // not yet reffered from the other instance
  comm->remaining_remote_mrq=0;

  //
  // create virtual control queue (local vcq)
  //
  int rc;
  //  unsigned long int vcq_flag=UTOFU_VCQ_FLAG_THREAD_SAFE | UTOFU_VCQ_FLAG_EXCLUSIVE;
  //unsigned long int vcq_flag=UTOFU_VCQ_FLAG_EXCLUSIVE;
  unsigned long int vcq_flag=UTOFU_VCQ_FLAG_THREAD_SAFE;
  //unsigned long int vcq_flag=0UL;
  rc=utofu_create_vcq(*tni_id, vcq_flag, &comm->local_vcq_hdl );
  if(UTOFU_SUCCESS != rc){
    fprintf(stderr, "rank %d: Failed at utofu_create_vcq()! rc=%d\n", m_rdma_comlib_myrank, rc);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  //
  // obtain id for local vcq
  //
  utofu_vcq_id_t local_vcq_id;
  rc=utofu_query_vcq_id(comm->local_vcq_hdl, &local_vcq_id);
  if(UTOFU_SUCCESS != rc){
    fprintf(stderr, "rank %d: Failed at utofu_query_vcq_id()! rc=%d\n", m_rdma_comlib_myrank, rc);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  //
  // send vcq_id to the remote target
  //
  uint64_t tmp_sbuff;
  uint64_t tmp_rbuff;
  int tag = ((m_rdma_comlib_tag-TAG_OFFSET) % (MAXNUM_TAG-TAG_OFFSET)) + TAG_OFFSET;
  tmp_sbuff=local_vcq_id;
#ifdef _RDMA_DEBUG
  fprintf(stderr, "rank %d: calling MPI_Sendrecv: tmp_sbuff=%lu\n", m_rdma_comlib_myrank, tmp_sbuff);
#endif
  MPI_Sendrecv(&tmp_sbuff, 1, MPI_UINT64_T, comm->from_rank, tag,
               &tmp_rbuff, 1, MPI_UINT64_T, comm->target_rank, tag,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  comm->target_vcq_id=tmp_rbuff;

  //
  // set the default path
  //
  rc=utofu_set_vcq_id_path(&comm->target_vcq_id, NULL);
  if(UTOFU_SUCCESS != rc){
    fprintf(stderr, "rank %d: Failed at utofu_set_vcq_id_path()! rc=%d\n", m_rdma_comlib_myrank, rc);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

}

void rdma_comlib_communicator_delete(rdma_comlib_communicator *comm){
  // release the resource
#ifdef _RDMA_DEBUG
  fprintf(stderr, "rank %d: rdma_comlib_communicator_delete: comm=%p, comm->ref_count=%p\n", m_rdma_comlib_myrank, comm, comm->ref_count);
#endif
  if(comm->ref_count == NULL){
    return;
  }
#ifdef _RDMA_DEBUG
  fprintf(stderr, "rank %d: rdma_comlib_communicator_delete: *comm->ref_count=%d\n", m_rdma_comlib_myrank, *comm->ref_count);
#endif

  if(*comm->ref_count > 0) { (*comm->ref_count)--; }

  if(*comm->ref_count==0){
    utofu_free_vcq(comm->local_vcq_hdl);
    comm->tni_id=0;
    comm->target_rank=0;
    comm->from_rank=0;
    comm->local_rank=0;
    free(comm->ref_count);
    comm->ref_count=NULL;
  }
}


void rdma_comlib_new_impl(rdma_comlib_data *id, const size_t *size, const int *has_external)
{
//
// Initialize comlib ID, set destination rank, local-send/recv buffer, data size(in 4byte)
// this version does not malloc,  the resource must be managed by the user
//
//   *id         : comlib_id
//                 user put/get data in id->sbuff : send buffef, id->rbuff : recv buffef
//   nic_id      : TNI ID [0...5]
//   dst_rank    : destination MPI rank
//   rcv_rank    : MPI rank from which data comes
//   size        : data size to be send/recv'd in byte unit.
//
//   uses own buffer:      has_external  = 0
//   uses extenral buffer: has_external != 0


#ifdef _RDMA_DEBUG
  id->fp = stderr;
  fprintf(id->fp,"rank %d: rdma_comlib_new_impl start.\n",m_rdma_comlib_myrank);
#endif

  // check arguments
  int err=0;
  if(*has_external){  // check if buffers are allocated externally
    if(id->sbuff == NULL){
      fprintf(stderr, "sbuff is not allocated\n");
      err=1;
    };
    if(id->rbuff == NULL){
      fprintf(stderr, "fbuff is not allocated\n");
      err=1;
    };
  }
  if(err){
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  //
  // set communicator
  //
  (*id->comm.ref_count)++;
#ifdef _RDMA_DEBUG
  fprintf(stderr, "rank %d: id->comm.ref_count=%p, *id->comm.ref_count=%d\n", m_rdma_comlib_myrank, id->comm.ref_count, *id->comm.ref_count);
#endif

  //
  // compute send/receive data size including polling buffer.
  //
  id->length = *size; // data length in byte

  id->data_size = (id->length + (size_t)4);  // additional last 4-byte component is used for polling data receive.

  if ( MAX_RDMA_DATASIZE < id->data_size ) {
    fprintf(stderr,"RDMA data size is too large.  data size should be data_size = %d < %d = MAX_RDMA_DATASIZE\n",id->data_size,MAX_RDMA_DATASIZE);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  if ( (id->data_size % ((size_t)4)) != 0 ) {
    fprintf(stderr,"RDMA data size shoule be a multiple of 4 bytes.  data_size = %d\n",id->data_size);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  id->data_len_int = (id->data_size) / ((size_t)4);

  //
  // allocate local send and local receive buffers
  //
  if( !(*has_external) ){  // skip if the resource is managed outside
    id->sbuff = (int *)malloc(id->data_size);
    id->rbuff = (int *)malloc(id->data_size);
  }


  //
  // set RDMA message tag
  //
  id->tag = ((m_rdma_comlib_tag-TAG_OFFSET + MAXNUM_TAG) % (MAXNUM_TAG-TAG_OFFSET)) + TAG_OFFSET;
  m_rdma_comlib_tag++;


  //
  // register local send buffer
  //
  int rc;
  rc=utofu_reg_mem(id->comm.local_vcq_hdl, (void *)(id->sbuff), id->data_size, 0, &id->local_stadd_sbuff);
  if(UTOFU_SUCCESS != rc){
    fprintf(stderr, "rank %d: Failed at utofu_reg_mem() for local send buffer! rc=%d\n", m_rdma_comlib_myrank, rc);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  id->sending_vcq_hdl=id->comm.local_vcq_hdl;
  id->sending_stadd=id->local_stadd_sbuff;

  //
  // register local receive buffer
  //
  rc=utofu_reg_mem(id->comm.local_vcq_hdl, (void *)(id->rbuff), id->data_size, 0, &id->local_stadd_rbuff);
  if(UTOFU_SUCCESS != rc){
    fprintf(stderr, "rank %d: Failed at utofu_reg_mem() for local recieve buffer! rc=%d\n", m_rdma_comlib_myrank, rc);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  //
  // send stadd for rbuffer to the remote target
  //
  uint64_t tmp_sbuff;
  uint64_t tmp_rbuff;
  tmp_sbuff=id->local_stadd_rbuff;
  //  fprintf(stderr, "rank %d: calling MPI_Sendrecv...    from=%d, target=%d, tag=%d\n",
  //          m_rdma_comlib_myrank, id->comm.from_rank, id->comm.target_rank, id->tag);
  MPI_Sendrecv(&tmp_sbuff, 1, MPI_UINT64_T, id->comm.from_rank, id->tag,
               &tmp_rbuff, 1, MPI_UINT64_T, id->comm.target_rank, id->tag,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  //  fprintf(stderr, "rank %d: calling MPI_Sendrecv, done from=%d, target=%d, tag=%d\n",
  //          m_rdma_comlib_myrank, id->comm.from_rank, id->comm.target_rank, id->tag);

  id->target_stadd_rbuff=tmp_rbuff;


  //
  // clear watchdog data in the additional last component
  //
  id->sbuff[id->data_len_int-1] = 0;
  id->rbuff[id->data_len_int-1] = MAX_LAST_COMPONENT-1;
  id->last_component=id->rbuff[id->data_len_int-1];

  MPI_Barrier(MPI_COMM_WORLD);

#ifdef _RDMA_DEBUG
  fprintf(id->fp,"message length=%ld  size=%ld  len_int=%ld\n", id->length, id->data_size, id->data_len_int); fflush(NULL);
  fprintf(id->fp,"local_rank = %3d, target_rank = %3d ",id->comm.local_rank, id->comm.target_rank);
  fprintf(id->fp,"local_stadd_sbuff = %ld, target_stadd_rbuff = %ld\n",
                 id->local_stadd_sbuff,id->target_stadd_rbuff);
  fprintf(id->fp,"rank %d: rdma_comlib_new OK.\n",m_rdma_comlib_myrank); fflush(NULL);
  fflush(NULL);
#endif

}


void rdma_comlib_new(rdma_comlib_data *id, const int *tni_id, const int *dst_rank, const int *rcv_rank, const size_t *size)
{
//
// Initialize comlib ID, set destination rank, local-send/recv buffer, data size(in 4byte)
//
//   *id         : comlib_id
//                 user put/get data in id->sbuff : send buffef, id->rbuff : recv buffef
//   *tni_id     : TNI ID [0...5]
//   *dst_rank   : destination MPI rank
//   *rcv_rank   : MPI rank from which data comes
//   *size       : data size to be send/recv'd in byte unit.
//
  // calls extended one with external buffer
  const int has_external_buffer=0;
  rdma_comlib_new_ext(id, tni_id, dst_rank, rcv_rank, size, &has_external_buffer);
}


void rdma_comlib_new_ext(rdma_comlib_data *id, const int *tni_id, const int *dst_rank, const int *rcv_rank, const size_t *size, const int *has_external)
{
//
// Initialize comlib ID, set destination rank, local-send/recv buffer, data size(in 4byte)
// this version does not malloc,  the resource must be managed by the user
//
//   *id         : comlib_id
//                 user put/get data in id->sbuff : send buffef, id->rbuff : recv buffef
//   *nic_id     : TNI ID [0...5]
//   *dst_rank   : destination MPI rank
//   *rcv_rank   : MPI rank from which data comes
//   *size       : data size to be send/recv'd in byte unit.
//
//   uses own buffer:      *has_external  = 0
//   uses extenral buffer: *has_external != 0


#ifdef _RDMA_DEBUG
  id->fp = stderr;
  fprintf(id->fp,"rank %d: rdma_comlib_new start.\n",m_rdma_comlib_myrank);
#endif

  // check arguments
  int err=0;
  if(*has_external){  // check if buffers are allocated externally
    if(id->sbuff == NULL){
      fprintf(stderr, "sbuff is not allocated\n");
      err=1;
    };
    if(id->rbuff == NULL){
      fprintf(stderr, "rbuff is not allocated\n");
      err=1;
    };
  }
  if(err){
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  //
  // create vcq
  //
  rdma_comlib_communicator_new(&id->comm, tni_id, dst_rank, rcv_rank);

  rdma_comlib_new_impl(id, size, has_external);
#ifdef _RDMA_DEBUG
  fprintf(id->fp,"rank %d: rdma_comlib_new_ext, done.: id->comm=%p\n",m_rdma_comlib_myrank, id->comm);
#endif


}

void rdma_comlib_new_duplicate(rdma_comlib_data *id, const rdma_comlib_data *id_org, const size_t *size, const int *has_external){
//
// Initialize comlib ID, using the same vcq as id_org
// this version does not malloc,  the resource must be managed by the user
//
//   *id         : comlib_id
//                 user put/get data in id->sbuff : send buffef, id->rbuff : recv buffef
//   *id_org     : comlib_id
//                 vcq info in this entry is copied to id
//   *size       : data size to be send/recv'd in byte unit.
//
//   uses own buffer:      *has_external  = 0
//   uses extenral buffer: *has_external != 0


#ifdef _RDMA_DEBUG
  id->fp = stderr;
  fprintf(id->fp,"rank %d: rdma_comlib_new_duplicate start. : id_org->comm=%p\n",m_rdma_comlib_myrank, &id_org->comm);
#endif

  id->comm=id_org->comm; // copy including the pointer to the reference counter
  rdma_comlib_new_impl(id, size, has_external);
}

void rdma_comlib_swap_vcq_for_sending(rdma_comlib_data *id1, rdma_comlib_data *id2){
  int rc1, rc2;
  rc1=utofu_dereg_mem(id1->comm.local_vcq_hdl, id1->local_stadd_sbuff, 0);
  rc2=utofu_dereg_mem(id2->comm.local_vcq_hdl, id2->local_stadd_sbuff, 0);

  if(rc1 != UTOFU_SUCCESS){
    fprintf(stderr, "rank %d: Failed at utofu_dreg_mem() for id1\n", m_rdma_comlib_myrank);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  if(rc2 != UTOFU_SUCCESS){
    fprintf(stderr, "rank %d: Failed at utofu_dreg_mem() for id2\n", m_rdma_comlib_myrank);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  // data in id1 : send through local_vcq and local_stadd in id2
  // data in id2 : send through local_vcq and local_stadd in id1
  rc1=utofu_reg_mem(id1->comm.local_vcq_hdl, (void *)(id2->sbuff), id2->data_size, 0, &id1->local_stadd_sbuff);
  rc2=utofu_reg_mem(id2->comm.local_vcq_hdl, (void *)(id1->sbuff), id1->data_size, 0, &id2->local_stadd_sbuff);
  if(rc1 != UTOFU_SUCCESS){
    fprintf(stderr, "rank %d: Failed at utofu_reg_mem() for id1\n", m_rdma_comlib_myrank);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  if(rc2 != UTOFU_SUCCESS){
    fprintf(stderr, "rank %d: Failed at utofu_reg_mem() for id2\n", m_rdma_comlib_myrank);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  // tell actuall vcq and stadd
  id1->sending_vcq_hdl=id2->comm.local_vcq_hdl;
  id2->sending_vcq_hdl=id1->comm.local_vcq_hdl;
  id1->sending_stadd=id2->local_stadd_sbuff;
  id2->sending_stadd=id1->local_stadd_sbuff;
}

void rdma_comlib_delete(rdma_comlib_data *id)
{
//
// Delete comlib ID
//
  rdma_comlib_delete_ext(id, 0);
}

void rdma_comlib_delete_ext(rdma_comlib_data *id, const int *has_external)
{
//
// Delete comlib ID
//
//   uses own buffer:      *has_external  = 0
//   uses extenral buffer: *has_external != 0

  // Note that this function must be called BEFORE MPI_Finalize
  //  Also it is user's responsiblity to make sure that all the
  //  relevant communications are finished.
  //
  //  MPI_Barrier(MPI_COMM_WORLD);

  // clean all mrq (remote mrq only)
  rdma_comlib_clear_mrq(id);

  // release the resource
  // N.B. local_stadd_sbuffer belongs to the local_vcq_hdl but may be used
  //      in the other comlib instances:
  utofu_dereg_mem(id->comm.local_vcq_hdl, id->local_stadd_sbuff, 0);
  utofu_dereg_mem(id->comm.local_vcq_hdl, id->local_stadd_rbuff, 0);
  rdma_comlib_communicator_delete(&id->comm);
  if( !(*has_external) ){
    free((void *)(id->sbuff));
    free((void *)(id->rbuff));
    id->sbuff = NULL;
    id->rbuff = NULL;
  }

  // cleaning
  id->length = 0;
  id->data_size = 0;
  id->data_len_int = 0;
  id->tag = 0;

  id->sending_stadd = 0;
  id->local_stadd_sbuff = 0;
  id->local_stadd_rbuff = 0;
  id->target_stadd_rbuff = 0;
  id->last_component=0;

  id->sending_vcq_hdl = 0;
}


void rdma_comlib_isendrecv(rdma_comlib_data *id)
{
//
// asynchronous RDMA send/receive start with comlib ID
//

  //
  // reset watchdog data for send.
  //
  id->sbuff[id->data_len_int-1] = (id->sbuff[id->data_len_int-1]+1) % MAX_LAST_COMPONENT;

  //
  // Start data send
  //
#ifdef RDMA_NO_REMOTE_MRQ_POLLING
  const unsigned long int mrq_flag=0UL;
#else
  const unsigned long int mrq_flag
    = UTOFU_ONESIDED_FLAG_REMOTE_MRQ_NOTICE;
  //  UTOFU_ONESIDED_FLAG_LOCAL_MRQ_NOTICE |
#endif
#ifdef RDMA_USE_CACHE_INJECTION
  const unsigned long int cache_injection_flag=UTOFU_ONESIDED_FLAG_CACHE_INJECTION;
#else
  const unsigned long int cache_injection_flag=0UL;
#endif
  const unsigned long int send_flags
    = UTOFU_ONESIDED_FLAG_TCQ_NOTICE | UTOFU_ONESIDED_FLAG_STRONG_ORDER
    | mrq_flag | cache_injection_flag;

  uint64_t edata=0;    // for mrq polling; the value is not used
  uintptr_t cbvalue=0; // for tcq polling; the value is not used

  int rc=0;
  rc=utofu_put(id->sending_vcq_hdl, id->comm.target_vcq_id, id->sending_stadd, id->target_stadd_rbuff, id->data_size,
edata, send_flags, (void *)cbvalue);
  //id->comm.remaining_remote_mrq++;
  if ( rc != UTOFU_SUCCESS){
    fprintf(stderr,"rank %d, %s at %d : rdma_comlib_isendrecv, utofu_put ERROR: %d\n",m_rdma_comlib_myrank ,__FILE__, __LINE__, rc);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

#ifdef _RDMA_DEBUG
  fprintf(id->fp,"rank %d: rdma_comlib_isendrecv OK.\n",m_rdma_comlib_myrank); fflush(NULL);
#endif

}


void rdma_comlib_isend_check(rdma_comlib_data *id)
{
//
// check/wait for RDMA isend finish.
//

#ifdef _RDMA_DEBUG
  fprintf(id->fp,"rank %d: rdma_comlib_isend_check start.\n",m_rdma_comlib_myrank); fflush(NULL);
#endif

  // tcq polling
  void *cbdata;
  int rc;
  do {
    rc = utofu_poll_tcq(id->sending_vcq_hdl, 0, &cbdata);
  } while (rc == UTOFU_ERR_NOT_FOUND);

  // check if the polling is successfully finished
  if(rc != UTOFU_SUCCESS){
    fprintf(stderr,"rank %d: %s at %d : rdma_comlib_isend_check, ERROR: %d, target_rank=%d\n",
            m_rdma_comlib_myrank, __FILE__, __LINE__, rc, id->comm.target_rank);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  // skip checking the call back data
  //   if(cdbdata==...){}

#ifdef _RDMA_DEBUG
  fprintf(id->fp,"rank %d: rdma_comlib_isend_check OK.\n",m_rdma_comlib_myrank); fflush(NULL);
#endif
}


void rdma_comlib_clear_mrq(rdma_comlib_data *id)
{
#ifdef RDMA_NO_REMOTE_MRQ_POLLING
  return;
#else

#ifdef _RDMA_DEBUG
  fprintf(id->fp,"rank %d: rdma_comlib_clear_mrq start, remaining_remote_mrq=%d.\n",m_rdma_comlib_myrank, id->comm.remaining_remote_mrq);
  fflush(NULL);
#endif

  for(int i=0; i < id->comm.remaining_remote_mrq; i++){
    struct utofu_mrq_notice notice;
    int rc=0;

    do {
      rc = utofu_poll_mrq(id->comm.local_vcq_hdl, 0UL, &notice);
    } while (rc == UTOFU_ERR_NOT_FOUND);
    if(rc != UTOFU_SUCCESS){
      fprintf(stderr,"rank %d: %s at %d : utofu_poll_mrq, ERROR: %d\n",
              m_rdma_comlib_myrank, __FILE__, __LINE__, rc);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    } else {
      switch(notice.notice_type){
      case UTOFU_MRQ_TYPE_RMT_PUT:
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
                  m_rdma_comlib_myrank, __FILE__, __LINE__, type_str);
          MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
      }
    }
  }
  id->comm.remaining_remote_mrq=0;

#ifdef _RDMA_DEBUG
  fprintf(id->fp,"rank %d: rdma_comlib_clear_mrq done.\n",m_rdma_comlib_myrank);
  fflush(NULL);
#endif
#endif
}

void rdma_comlib_irecv_check(rdma_comlib_data *id)
{
//
// check/wait for RDMA irecv finish.
//
#ifdef _RDMA_DEBUG
  fprintf(id->fp,"rank %d: rdma_comlib_irecv_check start.\n",m_rdma_comlib_myrank); fflush(NULL);
#endif
  int watch_dog=id->last_component;
  id->comm.remaining_remote_mrq++;
  while( id->rbuff[id->data_len_int-1] == watch_dog ){ };

  id->last_component=id->rbuff[id->data_len_int-1];  // update the watch dog data
#ifdef _RDMA_DEBUG
  fprintf(stderr, "rank: %d, local_vcq_hdl=%p: new id->last_compoient_id=%d\n",
          m_rdma_comlib_myrank, id->comm.local_vcq_hdl, id->last_component, id->rbuff[id->data_len_int-1]);
  fprintf(id->fp,"rank %d: rdma_comlib_irecv_check OK.\n",m_rdma_comlib_myrank); fflush(NULL);
#endif

}

void rdma_comlib_irecv_ok(rdma_comlib_data *id)
{
//
// set ok status for receive buffer: nothing to do
//
#ifdef _RDMA_DEBUG
  fprintf(id->fp,"rank %d: rdma_comlib_irecv_ok start.\n",m_rdma_comlib_myrank); fflush(NULL);
#endif
  //
  // reset watchdog data.
  //   moved to rdma_comlib_irecv_check
  // id->last_component=id->rbuff[id->data_len_int-1];
  //id->rbuff[id->data_len_int-1]=0;
#ifdef _RDMA_DEBUG
  fprintf(id->fp,"rank %d: rdma_comlib_irecv_ok OK.\n",m_rdma_comlib_myrank); fflush(NULL);
#endif
}

#endif
