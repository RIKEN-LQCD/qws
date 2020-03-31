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
//#include "single_fields.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <sys/time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define MAXCHAR 64

//extern int nodeid_f;
int nodeid_f=0;

typedef
struct timing_node
{
  struct timing_node *next;
  char *id;
  unsigned long call_number;
  double total_time;
  double recent_start_clock;
} timing_list;

timing_list *init_timing_node (void);
double get_clock (void);

void check_timing_ (const char *id)
{
  static timing_list *tl = NULL;
  double clock_now = get_clock ();
  timing_list *cur;

  if (tl == NULL) tl = init_timing_node ();
  cur = tl;
#ifdef _OPENMP
  if (omp_get_thread_num() !=0) return;
#endif
  if (id == NULL)
    {
      printf ("%8s\t%-32s\t%8s\t%12s\t%12s\n",
              "rank","func_id", "calls", "total(s)", "average(s)");
      for (;;)
        {
          timing_list *old = cur;
          if (cur == NULL || cur->id == NULL) break;
          printf ("%8d\t%-32s\t%8ld\t%12.6e\t%12.6e\n",
                  nodeid_f,
                  cur->id,
                  cur->call_number,
                  cur->total_time,
                  cur->total_time / (double) cur->call_number);
          cur = cur->next;
          free (old);
        }

      return;
    }

  for (;;)
    {
      if (cur->id == NULL)          /* it's a new id */
        {
          cur->id = (char*)id;
          cur->recent_start_clock = clock_now;
          break;                /* finished creating the new id */
        }
      else                      /* looking at cur */
        {
          if (strncmp (id, cur->id, MAXCHAR) == 0) /* this id */
            {
              if (cur->recent_start_clock == 0) /* start */
                {
                  cur->recent_start_clock = clock_now;
                }
              else              /* end */
                {
                  cur->call_number++;
                  cur->total_time += clock_now - cur->recent_start_clock;
                  cur->recent_start_clock = 0;
                }
              break;            /* finished this id */
            }
          else                  /* not this one */
            {
              if (cur->next == NULL) /* end */
                {
                  cur = init_timing_node ();
                  cur->next = tl;
                  tl = cur;
                }
              else              /* not end */
                {
                  cur = cur->next;
                }
            }
        }
    }
}

timing_list *
init_timing_node (void)
{
  timing_list *cur = (timing_list *) malloc (sizeof (timing_list));
  cur->next = NULL;
  cur->id = NULL;
  cur->call_number = 0;
  cur->total_time = 0;
  cur->recent_start_clock = 0;
  return cur;
}

double
get_clock (void)
{
  struct timeval tp;
  //clock_gettime is not tested yet on K computer
  //  clock_gettime(CLOCK_REALTIME, &tp);
  //  return ((double) tp.tv_sec + (double) tp.tv_nsec * 1e-9);

  gettimeofday(&tp, NULL);
  return ((double) tp.tv_sec + (double) tp.tv_usec * 1e-6);
}

///////////////////////////////////////////////////////
// print out timing statistics 
// exposed to Fortran with underscore convension
///////////////////////////////////////////////////////
void print_timing_(void) { check_timing_ (NULL); }

