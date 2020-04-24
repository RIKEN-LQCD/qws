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

// Note
// currenly available for one region power measurement

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

#include "pwr.h"
PWR_Cntxt cntxt = NULL;
PWR_Obj obj = NULL;
int rc;

unsigned long count = 0;
double total_power= 0.0;
double energy_old = 0.0;
PWR_Time ts_old = 0;

extern FILE *para_outputfile;

void power_api_init(){
  //rc = PWR_CntxtInit(PWR_CNTXT_DEFAULT, PWR_ROLE_APP, "app", &cntxt);
  rc = PWR_CntxtInit(PWR_CNTXT_FX1000, PWR_ROLE_APP, "app", &cntxt);
  if (rc != PWR_RET_SUCCESS) {
    printf("CntxtInit Failed\n");
    exit(1);
  }
  rc = PWR_CntxtGetObjByName(cntxt, "plat.node", &obj);
  if (rc != PWR_RET_SUCCESS) {
    printf("CntxtGetObjByName Failed\n");
    exit(1);
  }
}

void power_api_finalize(){
  fprintf(para_outputfile, "bicgstab_precdd_s_iter_ : ave_power = %lf\n", total_power/(double)(count/2));
  PWR_CntxtDestroy(cntxt);
}

void power_measure_power(){
  double energy;
  PWR_Time ts;

  //rc = PWR_ObjAttrGetValue(obj, PWR_ATTR_ENERGY, &energy, &ts);
  rc = PWR_ObjAttrGetValue(obj, PWR_ATTR_MEASURED_ENERGY, &energy, &ts);
  if (rc != PWR_RET_SUCCESS) {
    printf("ObjAttrGetValue Failed (rc = %d)\n", rc);
    exit(1);
  }

  if (count%2==0){
    energy_old = energy;
    ts_old = ts;
  } else {
    total_power += (energy - energy_old) / ((ts - ts_old) / 1000000000.0);
    //double ave_power =  (energy - energy_old) / ((ts - ts_old) / 1000000000.0);
    //fprintf(para_outputfile, "bicgstab_precdd_s_iter_ : ave_power = %lf\n", ave_power);
  }
  count++;
}
