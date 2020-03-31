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
#include <stddef.h>
#include "eml_lib.h"

#define MAX_SIZE 20000000

static unsigned char malloc_area[MAX_SIZE] __attribute__((aligned(256)));
extern unsigned char *malloc_p = malloc_area;
static size_t malloc_len;
  
int eml_posix_memalign( void **out_p, int    w, size_t size)
{
  unsigned char *p = malloc_p;
  size_t pad;
  pad =(size_t)p % w;
  p += (w - pad);
  if (malloc_len+size+w-pad > MAX_SIZE) return 1;
  malloc_len += ( size + (w - pad ) );
  malloc_p += (size + (w - pad ) );
  *out_p = p;
  return 0;
}

void *eml_malloc( size_t size )
{
  unsigned char *p;
  if(eml_posix_memalign((void**)&p, 16, size))
    return NULL;

  return p;
}

void *eml_memcpy(void *s1, const void *s2, size_t n)
{
  char *__restrict__ p1 = (char *)s1;
  const char *__restrict__ p2 = (const char *)s2;

  for(long int i=0; i<(long int)n; i++) {
    p1[i] = p2[i];
  }
  return (s1);
}

void eml_free( void *p )
{
  return;
}

int eml_printf(const char *fmt, ...)
{
  return 0;
}

int eml_exit(int in)
{
  return 0;
}
