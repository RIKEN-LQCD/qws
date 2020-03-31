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
#ifdef _HPC_RDMA
////////////////////////////////////////////////////////////////////////////////
// RDMA comlib using FJMPI extention 
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
//  Re-initialization does not work in this FJMPI RDMA implementation.
//
//void rdma_comlib_finalize(void);
//  this terminate RDMA functionality. This is called before MPI_Finalize.
//  Re-initialization does not work in this FJMPI RDMA implementation.
//
//void rdma_comlib_new(rdma_comlib_data *id, const int nic_id, const int dst_rank, const size_t *size);
//  Creaate, register and allocate the communication buffers and data for RDMA 1-to-1 communication.
//  MPI Process called this function can send data to another MPI process specified by "dst_rank"
//  throough virtual NIC with "nic_id". The unit of the data size (*size) is in byte (sizeof()).
//  The data buffers for send/receive is located in *id and allocated by this comlib_new.
//  The data to be sent should be stored in id->sbuff (pointer to int).
//  The data received will be be stored in id->rbuff (pointer to int).
//  All send/receive/check for the data communication are applied to the data structure *id.
//
//  rdma_comlib_data *id : pointer to rdma comlib data structure.
//  const int nic_id     : specify virtual NIC id [0...3] for this communication.
//  const int dst_rank   : specify send destination MPI RANK
//  const size_t *size   : pointer to data size in byte.
//
//  Note:
//    Re-creation for the *id does not work in this FJMPI RDMA implementation.
//
//void rdma_comlib_delete(rdma_comlib_data *id);
//  Delete rdma comlib data structure. The communication buffer stored in *id is deallocated. 
//  Note:
//    Re-creation(rdma_comlib_new) for the *id after rdma_comlib_delete does not work 
//    in this FJMPI RDMA implementation because the FJMPI RDMA MEMID cannot be freed.
//
//void rdma_comlib_isendrecv(rdma_comlib_data *id);
//  Start send/receive communication specified by *id. This will almost non-blocking.
//  In this function, receiver-readiness is checked (synchronous at this point).
//  Sending data is asynchronous.
//  id->sbuff and id->rbuff should not be touched before checking send/reseive finishing.
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
//  Set OK flag for the receive buffer. User should call this function after
//  fnishing reading data form id->rbuff in order to notify the sender can send
//  new data to this buffer.
//
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
// void rdma_comlib_new_ext(rdma_comlib_data *id, const int nic_id, const int dst_rank, const size_t *size, const int has_external)
// void rdma_comlib_delete_ext(rdma_comlib_data *id, const int has_external)
//   they are almost the same as rdma_comlib_new, rdma_comlib_delete,
//   but depeding on has_external:
//     has_exteranl == 0 : allocate/deallocate buffer( = the same as before)
//     has_exteranl != 0 : assumes that the buffers have been already allocated/ will
//                         be deleted by the user
//
//


#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <mpi-ext.h>

//#define  _RDMA_DEBUG


const int MAXNUM_NICID = 4;             //  Maximum nuber of NIC id
const int MAX_RDMA_DATASIZE = 16777212; //  RDMA put/get MAX data size in byte
const int MAXNUM_TAG   = 15;            //  RDMA put/get MAX number of message tag
const int MAXNUM_MEMID = 511;           //  RDMA MEM ID MAX number
const int LAST_COMPONENT  = 0xDEADBEEF;  //  watchdog data for receive data polling.
const int ASK_TARGET_STAT = 0xFFFFFFFF;  //  watchdog data for target buffer status polling.

static int m_rdma_comlib_is_initialized = 0;
static int m_rdma_comlib_myrank = 0; // MPI LOCAL RANK is stored.
static int m_rdma_comlib_memid = 0;  // counts RDMA MEM ID.      [0..511]
static int m_rdma_comlib_tag = 0;    // counts RDMA message tag. [0..15]

#include "rdma_comlib.h"

int rdma_comlib_get_ssize(const rdma_comlib_data *id)
{
   return (*id).length;
}

void rdma_comlib_init(void)
{ 
//
// Initialize RDMA Fujitsu MPI extention.
//
  if (0 == m_rdma_comlib_is_initialized) {

    FJMPI_Rdma_init(); 
    MPI_Comm_rank(MPI_COMM_WORLD,&m_rdma_comlib_myrank);
    m_rdma_comlib_is_initialized = 1;

    if ( 0 == m_rdma_comlib_myrank ) {
      printf("%%%% FJMPI Extension (RDMA) is used.\n");
    }
  }
}

void rdma_comlib_finalize(void) { 
  if (1 == m_rdma_comlib_is_initialized) {
    FJMPI_Rdma_finalize(); 
  }
}

void rdma_comlib_new(rdma_comlib_data *id, const int nic_id, const int dst_rank, const size_t *size)
{
//
// Initialize comlib ID, set destination rank, local-send/recv buffer, data size(in 4byte)
//
//   *id         : comlib_id
//                 user put/get data in id->sbuff : send buffef, id->rbuff : recv buffef
//   nic_id      : virtual NIC ID [0...3]
//   dst_rank    : destination MPI rank
//   size        : data size to be send/recv'd in byte unit.
//
  // calls extended one with external buffer
  rdma_comlib_new_ext(id, nic_id, dst_rank, size, 0);
}


void rdma_comlib_new_ext(rdma_comlib_data *id, const int nic_id, const int dst_rank, const size_t *size, const int has_external)
{
//
// Initialize comlib ID, set destination rank, local-send/recv buffer, data size(in 4byte)
// this version does not malloc,  the resource must be managed by the user
//
//   *id         : comlib_id
//                 user put/get data in id->sbuff : send buffef, id->rbuff : recv buffef
//   nic_id      : virtual NIC ID [0...3]
//   dst_rank    : destination MPI rank
//   size        : data size to be send/recv'd in byte unit.
//
//   uses own buffer:      has_external  = 0
//   uses extenral buffer: has_external != 0

  
#ifdef _RDMA_DEBUG
  id->fp = stderr;
  fprintf(id->fp,"rank %d: rdma_comlib_new start.\n",m_rdma_comlib_myrank);
#endif

  if(has_external){  // check if external buffers are allocated
    if(id->sbuff == NULL){
      fprintf(stderr, "sbuff is not allocated\n");
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    };
    if(id->rbuff == NULL){
      fprintf(stderr, "fbuff is not allocated\n");
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    };
  }

  //
  // set NIC ID options
  //
  if ( (MAXNUM_NICID <= nic_id) || ( nic_id < 0 ) ) {
    fprintf(stderr,"NIC ID should be greater than -1 and less than MAXNUM_NICID. 0 <= nic_id = %d < %d = MAXNUM_NICID\n",nic_id,MAXNUM_NICID);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  id->nic_id = nic_id;

  switch(id->nic_id){
    case 0:
      id->put_options = (FJMPI_RDMA_LOCAL_NIC0 | FJMPI_RDMA_REMOTE_NIC0 | FJMPI_RDMA_STRONG_ORDER);
      id->get_options = (FJMPI_RDMA_LOCAL_NIC0 | FJMPI_RDMA_REMOTE_NIC0 | FJMPI_RDMA_STRONG_ORDER);
      id->poll_options = FJMPI_RDMA_NIC0;
      break;
    case 1:
      id->put_options = (FJMPI_RDMA_LOCAL_NIC1 | FJMPI_RDMA_REMOTE_NIC1 | FJMPI_RDMA_STRONG_ORDER);
      id->get_options = (FJMPI_RDMA_LOCAL_NIC1 | FJMPI_RDMA_REMOTE_NIC1 | FJMPI_RDMA_STRONG_ORDER);
      id->poll_options = FJMPI_RDMA_NIC1;
      break;
    case 2:
      id->put_options = (FJMPI_RDMA_LOCAL_NIC2 | FJMPI_RDMA_REMOTE_NIC2 | FJMPI_RDMA_STRONG_ORDER);
      id->get_options = (FJMPI_RDMA_LOCAL_NIC2 | FJMPI_RDMA_REMOTE_NIC2 | FJMPI_RDMA_STRONG_ORDER);
      id->poll_options = FJMPI_RDMA_NIC2;
      break;
    case 3:
      id->put_options = (FJMPI_RDMA_LOCAL_NIC3 | FJMPI_RDMA_REMOTE_NIC3 | FJMPI_RDMA_STRONG_ORDER);
      id->get_options = (FJMPI_RDMA_LOCAL_NIC3 | FJMPI_RDMA_REMOTE_NIC3 | FJMPI_RDMA_STRONG_ORDER);
      id->poll_options = FJMPI_RDMA_NIC3;
      break;
    default:
      break;
  }

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
  if(!has_external){  // skip if the resource is managed outside
    id->sbuff = (int *)malloc(id->data_size);
    id->rbuff = (int *)malloc(id->data_size);
  }

  //
  // set MPI rank
  //
  id->local_rank = m_rdma_comlib_myrank;
  id->target_rank = dst_rank;

  //
  // set RDMA message tag
  //
  id->tag = (m_rdma_comlib_tag % MAXNUM_TAG);
  m_rdma_comlib_tag++;

  //
  // set RDMA MEM ID for local send buffer
  //
  id->memid_sbuff = (m_rdma_comlib_memid % MAXNUM_MEMID);
  m_rdma_comlib_memid++;

  //
  // set RDMA MEM ID for local receive buffer
  //
  id->memid_rbuff = (m_rdma_comlib_memid % MAXNUM_MEMID);
  m_rdma_comlib_memid++;

  //
  // register and get the RDMA address of the local source buffer
  //
  // Note: RDMA MEM ID cannot be freed. 
  //
  id->addr_local_sbuff = FJMPI_Rdma_reg_mem(id->memid_sbuff, (void *)(id->sbuff), id->data_size);

  //
  // register the local receive buffer (RDMA exposed to other ranks)
  //
  // Note: RDMA MEM ID cannot be freed. 
  //
  FJMPI_Rdma_reg_mem(id->memid_rbuff, (void *)(id->rbuff), id->data_size);

  //
  // get the RDMA address of target receive buffer
  //
  while((id->addr_target_rbuff = FJMPI_Rdma_get_remote_addr(id->target_rank, id->memid_rbuff)) == FJMPI_RDMA_ERROR);

  //
  // clear watchdog data in the additional last component
  //
  id->sbuff[id->data_len_int-1] = 0;
  id->rbuff[id->data_len_int-1] = 0;

  MPI_Barrier(MPI_COMM_WORLD);

#ifdef _RDMA_DEBUG
  fprintf(id->fp,"message length=%ld  size=%ld  len_int=%ld\n", id->length, id->data_size, id->data_len_int); fflush(NULL);
  fprintf(id->fp,"local_rank = %3d, target_rank = %3d ",id->local_rank,id->target_rank);
  fprintf(id->fp,"addr_local_sbuff = %ld, addr_target_rbuff = %ld\n",
                 id->addr_local_sbuff,id->addr_target_rbuff);
  fprintf(id->fp,"rank %d: rdma_comlib_new OK.\n",m_rdma_comlib_myrank); fflush(NULL);
  fflush(NULL);
#endif

}

void rdma_comlib_delete(rdma_comlib_data *id)
{
//
// Delete comlib ID
//
  rdma_comlib_delete_ext(id, 0);
}

void rdma_comlib_delete_ext(rdma_comlib_data *id, const int has_external)
{
//
// Delete comlib ID
//
//   uses own buffer:      has_external  = 0
//   uses extenral buffer: has_external != 0
  //
  // RDMA MEM ID registered cannot be freed. 
  //

  id->nic_id = 0;
  id->length = 0;
  id->data_size = 0;
  id->data_len_int = 0;
  id->local_rank = 0;
  id->target_rank = 0;
  id->tag = 0;
  m_rdma_comlib_tag--;

  id->memid_sbuff = 0;

  id->memid_rbuff = 0;

  id->addr_local_sbuff = 0;
  id->addr_target_rbuff = 0;

  if(!has_external){
    free((void *)(id->sbuff));
    free((void *)(id->rbuff));
    id->sbuff = NULL;
    id->rbuff = NULL;
  }
}

void rdma_comlib_irecv(rdma_comlib_data *id)
{
#ifdef _RDMA_DEBUG
  fprintf(id->fp,"rank %d: rdma_comlib_irecv start.\n",m_rdma_comlib_myrank); fflush(NULL);
#endif

  //
  // check remote buffer status
  //
  int remote_ok = 0;

  id->sbuff[id->data_len_int-1] = ASK_TARGET_STAT;

  while( 1 ) {
  //
  // Wait until target buffer becomes ready.
  //

    //
    // get last component(wachdog) data from remote receive buffer
    //
    if ( FJMPI_Rdma_get(id->target_rank, id->tag,
                        id->addr_target_rbuff + sizeof(int)*(id->data_len_int-1),
                        id->addr_local_sbuff  + sizeof(int)*(id->data_len_int-1),
                        sizeof(int), id->get_options) ){
      fprintf(stderr,"%s at %d : rdma_comlib_irecv, FJMPI_Rdma_get ERROR\n",__FILE__,__LINE__);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    //
    // Check tagert buffer status
    //
    while( 1 ) {

      if (id->sbuff[id->data_len_int-1] == 0) {
        //
        // target buffer is ready for receive.
        //
        remote_ok = 1;
        break;
      }

      if (id->sbuff[id->data_len_int-1] == LAST_COMPONENT) {
        //
        // target buffer is still in use.
        //
        id->sbuff[id->data_len_int-1] = ASK_TARGET_STAT;
        break;
      }

    }

    //
    // Clear RDMA get queue.
    //
    while( FJMPI_Rdma_poll_cq(id->poll_options,NULL) != FJMPI_RDMA_NOTICE );

    if ( 1 == remote_ok ) break;
  }

#ifdef _RDMA_DEBUG
  fprintf(id->fp,"rank %d: rdma_comlib_irecv OK.\n",m_rdma_comlib_myrank); fflush(NULL);
#endif

}

void rdma_comlib_isendrecv(rdma_comlib_data *id)
{
//
// asynchronous RDMA send/receive start with comlib ID
//

  //
  // reset watchdog data for send.
  //
  id->sbuff[id->data_len_int-1] = LAST_COMPONENT;

  //
  // Start data send
  //
  if ( FJMPI_Rdma_put(id->target_rank, id->tag,
                      id->addr_target_rbuff,
                      id->addr_local_sbuff,
                      id->data_size, id->put_options) ) {
    fprintf(stderr,"%s at %d : rdma_comlib_isendrecv, FJMPI_Rdma_put ERROR\n",__FILE__,__LINE__);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

#ifdef _RDMA_DEBUG
  fprintf(id->fp,"rank %d: rdma_comlib_isendrecv OK.\n",m_rdma_comlib_myrank); fflush(NULL);
#endif

}

void get_rdma_comlib_isend_finish_tag(int *nicid, int *tag)
{
  struct FJMPI_Rdma_cq cq;
  while( FJMPI_Rdma_poll_cq(*nicid, &cq) != FJMPI_RDMA_NOTICE );
  *tag = cq.tag;
}

void rdma_comlib_isend_check(rdma_comlib_data *id)
{
//
// check/wait for RDMA isend finish.
//

#ifdef _RDMA_DEBUG
  fprintf(id->fp,"rank %d: rdma_comlib_isend_check start.\n",m_rdma_comlib_myrank); fflush(NULL);
#endif

  while( FJMPI_Rdma_poll_cq(id->poll_options, NULL) != FJMPI_RDMA_NOTICE );

/*
  struct FJMPI_Rdma_cq cq;
  while( FJMPI_Rdma_poll_cq(id->poll_options, &cq) != FJMPI_RDMA_NOTICE );
  if ( cq.pid != id->target_rank) {
    fprintf(stderr,"%s at %d : rdma_comlib_isend_check, CQ RANK ERROR, cq.pid=%d, target_rank=%d\n",
                    __FILE__,__LINE__,(id->cq).pid,(id->target_rank));
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  while( FJMPI_Rdma_poll_cq(id->poll_options, &(id->cq)) != FJMPI_RDMA_NOTICE );
  if ( ((id->cq).pid != (id->target_rank)) || ((id->cq).tag != (id->tag)) ) {
    fprintf(stderr,"%s at %d : rdma_comlib_isend_check, CQ ERROR, cq.pid=%d, target_rank=%d, cq.tag=%d, tag=%d\n",
                    __FILE__,__LINE__,(id->cq).pid,(id->target_rank),(id->cq).tag,(id->tag));
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
*/

#ifdef _RDMA_DEBUG
  fprintf(id->fp,"rank %d: rdma_comlib_isend_check OK.\n",m_rdma_comlib_myrank); fflush(NULL);
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

  while( id->rbuff[id->data_len_int-1] != LAST_COMPONENT );

#ifdef _RDMA_DEBUG
  fprintf(id->fp,"rank %d: rdma_comlib_irecv_check OK.\n",m_rdma_comlib_myrank); fflush(NULL);
#endif
}

void rdma_comlib_irecv_ok(rdma_comlib_data *id)
{
//
// set ok status for receive buffer
//
#ifdef _RDMA_DEBUG
  fprintf(id->fp,"rank %d: rdma_comlib_irecv_ok start.\n",m_rdma_comlib_myrank); fflush(NULL);
#endif
  //
  // reset watchdog data.
  //
  id->rbuff[id->data_len_int-1] = 0;

#ifdef _RDMA_DEBUG
  fprintf(id->fp,"rank %d: rdma_comlib_irecv_ok OK.\n",m_rdma_comlib_myrank); fflush(NULL);
#endif
}

#endif
