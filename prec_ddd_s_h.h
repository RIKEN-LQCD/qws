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
#ifndef _PREC_DDD_H_H
#define _PREC_DDD_H_H

void prec_s_(scs_t *out, const scs_t *in, const int *nsap, const int *nm);


void      ddd_in_s_(scs_t *out, const scs_t *in, const int *domain);
void jinv_ddd_in_s_(scs_t   *x, const scs_t  *b, const int *domain, const int *maxiter);
void ddd_out_pre_s_(            const scs_t *in, const int *domain);
void ddd_out_pos_s_(scs_t *out, const scs_t *in, const int *domain, const float factor);

void debug_half(scs_t *xs, const scs_t *bs, const int *nsap, const int *nm);

void assign_mult_mSAP_s_h_(scs_t *xs, const scs_t *bs, const int *nsap, const int *nm)
//
// Multiply Half precision SAP preconditioner
//
// input and output are converted to half/single precision from single/half precision, respectively
//
// xs = conv_SP( MSAP_HP * conv_HP(bs) )
//
// MSAP_HP : half presision SAP preconditioner
//
// MSAP = Ksap(sum_{j=0}^{nsap} (1-DKsap)^j) \sim D^{-1}
//
//    xs : output quark field in single precision
//    bs :  input quark field in single precision
//  nsap : sap fixed iteration count
//    nm : Jacobbi fixed iteration count for the approximate inverse of Dee/Doo in a even/odd domain
//
//
{
  static sch_t *x = nullptr;
  static sch_t *b = nullptr;
  static sch_t *s = nullptr;
  static sch_t *q = nullptr;
  if (nullptr == x) x = (sch_t *)malloc( sizeof(sch_t) * vols*2);
  if (nullptr == b) b = (sch_t *)malloc( sizeof(sch_t) * vols*2);
  if (nullptr == s) s = (sch_t *)malloc( sizeof(sch_t) * vols*2);
  if (nullptr == q) q = (sch_t *)malloc( sizeof(sch_t) * vols*2);

        sch_t * __restrict__ xe = &x[vols*domain_e];
        sch_t * __restrict__ xo = &x[vols*domain_o];
  const sch_t * __restrict__ be = &b[vols*domain_e];
  const sch_t * __restrict__ bo = &b[vols*domain_o];
        sch_t * __restrict__ se = &s[vols*domain_e];
        sch_t * __restrict__ so = &s[vols*domain_o];
        sch_t * __restrict__ qe = &q[vols*domain_e];
        sch_t * __restrict__ qo = &q[vols*domain_o];

  if (*nsap > 10) { printf("%s : %s nsap > 10\n",__FILE__,__func__); exit(1); }
  if (*nm   > 10) { printf("%s : %s nn   > 10\n",__FILE__,__func__); exit(1); }

#if 0
//  prec_s_(xs, bs, nsap, nm);
#else
//  debug_half(xs,bs,nsap,nm);

  //////////////////////////////
  // convert float to half
  // b <= bs
  //////////////////////////////
  float bsnorm;
  assign_q_s2h( b, bs, &bsnorm); 

#pragma omp parallel for
  for(int i=0; i<vols*2; i++)
    for(int j=0; j<24; j++)
      for(int v=0; v<VLENS; v++)
        s[i].ccs[j].v[v] = b[i].ccs[j].v[v];


  if (npe[1]==1 || npe[2]==1 || npe[3]==1){
#pragma omp parallel for
    for(int i=0; i<vols; i++)
      for(int j=0; j<24; j++)
	for(int v=0; v<VLENS; v++)
	  xo[i].ccs[j].v[v] = 0;
  }

  ///////////////////////////////////
  // SAP iteration
  ///////////////////////////////////
  for (int isap=1; isap < *nsap; isap++) {

    ///////////////////////////////////
    // xe = Aee se
    ///////////////////////////////////
    jinv_ddd_in_h_(xe, se, &domain_e, nm);

    ///////////////////////////////////
    // Send for Doe xe
    ///////////////////////////////////
    ddd_out_pre_h_( x, &domain_o);
    
    ///////////////////////////////////
    //  qe = Dee xe
    ///////////////////////////////////
    ddd_in_h_( qe, xe, &domain_e);

    ///////////////////////////////////
    // se = se + be - qe
    ///////////////////////////////////
#pragma omp parallel for
    for(int i=0; i<vols; i++)
      for(int j=0; j<24; j++)
        for(int v=0; v<VLENS; v++)
          se[i].ccs[j].v[v] += be[i].ccs[j].v[v] - qe[i].ccs[j].v[v];

    ///////////////////////////////////
    // so = so - Doe xe  (Recv)
    ///////////////////////////////////
    ddd_out_pos_h_( so, x, &domain_o, (float)kappa);
   
    ///////////////////////////////////
    // xo = Aoo so
    ///////////////////////////////////
    jinv_ddd_in_h_( xo, so, &domain_o, nm);

    ///////////////////////////////////
    // Send for Deo xo
    ///////////////////////////////////
    ddd_out_pre_h_( x, &domain_e);

    ///////////////////////////////////
    // qo = Doo xo
    ///////////////////////////////////
    ddd_in_h_( qo, xo, &domain_o);

    ///////////////////////////////////
    // so = so + bo - qo
    ///////////////////////////////////
#pragma omp parallel for
    for(int i=0; i<vols; i++)
      for(int j=0; j<24; j++)
        for(int v=0; v<VLENS; v++)
          so[i].ccs[j].v[v] += bo[i].ccs[j].v[v] - qo[i].ccs[j].v[v];

    ///////////////////////////////////
    // se = se - Deo xo  (Recv)
    ///////////////////////////////////
    ddd_out_pos_h_( se, x, &domain_e, (float)kappa);

  } // end for isap

  if (npe[1]==1 || npe[2]==1 || npe[3]==1){
#pragma omp parallel for
    for(int i=0; i<vols; i++)
      for(int j=0; j<24; j++)
	for(int v=0; v<VLENS; v++)
	  xo[i].ccs[j].v[v] = 0;
  }

  ///////////////////////////////////
  // xe = Aee se
  ///////////////////////////////////
  jinv_ddd_in_h_( xe, se, &domain_e, nm);

  ///////////////////////////////////
  // Send for Doe xe
  ///////////////////////////////////
  ddd_out_pre_h_( x, &domain_o);

  ///////////////////////////////////
  // so = so - Doe xe
  ///////////////////////////////////
  ddd_out_pos_h_( so, x, &domain_o, (float)kappa);

  ///////////////////////////////////
  // xo = Aoo so
  ///////////////////////////////////
  jinv_ddd_in_h_( xo, so, &domain_o, nm);

  //////////////////////////////
  // convert half to float
  // xs <= x
  //////////////////////////////
  assign_q_h2s( xs, x, &bsnorm); 

#endif

}

void assign_mult_wd_s_(scs_t *x, const scs_t *b)
//
// Multiply Wilson/Clover operator (single precision)
//
// x = D b
//
{
        scs_t * __restrict__ xe = &x[vols*domain_e];
        scs_t * __restrict__ xo = &x[vols*domain_o];
  const scs_t * __restrict__ be = &b[vols*domain_e];
  const scs_t * __restrict__ bo = &b[vols*domain_o];

  ///////////////////////
  // Send for Deo bo
  ///////////////////////
  ddd_out_pre_s_( b, &domain_e);

  ///////////////////////
  // xe = Dee be
  ///////////////////////
  ddd_in_s_( xe, be, &domain_e);

  ///////////////////////
  // xe = xe + Deo bo
  //    = xe - kappa Meo bo
  ///////////////////////
  ddd_out_pos_s_( xe, b, &domain_e, (float)mkappa);

  ///////////////////////
  // Send for Doe be
  ///////////////////////
  ddd_out_pre_s_( b, &domain_o);

  ///////////////////////
  // xo = Dee bo
  ///////////////////////
  ddd_in_s_( xo, bo, &domain_o);

  ///////////////////////
  // xo = xo + Doe be 
  //    = xo - kappa Moe be
  ///////////////////////
  ddd_out_pos_s_( xo, b, &domain_o, (float)mkappa);

}

void debug_half(scs_t *xs, const scs_t *bs, const int *nsap, const int *nm)
//
// debug
//
{
  const float tol = 1.0e-5f;

  static sch_t *x = nullptr;
  static sch_t *b = nullptr;
  if (nullptr == x) x = (sch_t *)malloc( sizeof(sch_t) * vols*2);
  if (nullptr == b) b = (sch_t *)malloc( sizeof(sch_t) * vols*2);

  static scs_t *xxs = nullptr;
  if (nullptr == xxs) xxs = (scs_t *)malloc( sizeof(scs_t) * vols*2);

        sch_t * __restrict__ xe = &x[vols*domain_e];
        sch_t * __restrict__ xo = &x[vols*domain_o];
  const sch_t * __restrict__ be = &b[vols*domain_e];
  const sch_t * __restrict__ bo = &b[vols*domain_o];

  //////////////////////
  // check
  // xe = Dee be
  // xo = Doo bo
  //////////////////////
 
  ddd_in_s_( &xxs[vols*domain_e], &bs[vols*domain_e], &domain_e);
  ddd_in_s_( &xxs[vols*domain_o], &bs[vols*domain_o], &domain_o);

  float bsnorm;
  assign_q_s2h( b, bs, &bsnorm); 
  ddd_in_h_( xe, be, &domain_e);
  ddd_in_h_( xo, bo, &domain_o);
  assign_q_h2s( xs, x, &bsnorm); 

  for(int i=0; i<vols*2; i++)
    for(int j=0; j<24; j++)
      for(int v=0; v<VLENS; v++){
        float rr = xxs[i].ccs[j].v[v] - xs[i].ccs[j].v[v];
        if ( rr*rr > tol*tol  && rank == 0) 
          printf("CHKDEE: %3d %3d %3d %4d (SP)x = %14.6e (HP)x = %14.6e diff = %14.6e\n",
                  rank,v,j,i,xxs[i].ccs[j].v[v],
                              xs[i].ccs[j].v[v],rr);
      }

  //////////////////////
  // check
  // xe = Aee be
  // xo = Aoo bo
  //////////////////////
 
  jinv_ddd_in_s_( &xxs[vols*domain_e], &bs[vols*domain_e], &domain_e, nm);
  jinv_ddd_in_s_( &xxs[vols*domain_o], &bs[vols*domain_o], &domain_o, nm);

  assign_q_s2h( b, bs, &bsnorm); 
  jinv_ddd_in_h_( xe, be, &domain_e, nm);
  jinv_ddd_in_h_( xo, bo, &domain_o, nm);
  assign_q_h2s( xs, x, &bsnorm); 

  for(int i=0; i<vols*2; i++)
    for(int j=0; j<24; j++)
      for(int v=0; v<VLENS; v++){
        float rr = xxs[i].ccs[j].v[v] - xs[i].ccs[j].v[v];
        if ( rr*rr > tol*tol  && rank == 0) 
          printf("CHKAEE: %3d %3d %3d %4d (SP)x = %14.6e (HP)x = %14.6e diff = %14.6e\n",
                  rank,v,j,i,xxs[i].ccs[j].v[v],
                              xs[i].ccs[j].v[v],rr);
      }

  //////////////////////
  // check
  // xe = be - Deo bo
  // xo = bo - Doe be
  //////////////////////

#pragma omp parallel for
  for(int i=0; i<vols*2; i++)
    for(int j=0; j<24; j++)
      for(int v=0; v<VLENS; v++)
        xxs[i].ccs[j].v[v] = bs[i].ccs[j].v[v];

  ddd_out_pre_s_(                      bs, &domain_e);
  ddd_out_pos_s_( &xxs[vols*domain_e], bs, &domain_e, (float)kappa);
  ddd_out_pre_s_(                      bs, &domain_o);
  ddd_out_pos_s_( &xxs[vols*domain_o], bs, &domain_o, (float)kappa);

  assign_q_s2h( b, bs, &bsnorm); 
#pragma omp parallel for
  for(int i=0; i<vols*2; i++)
    for(int j=0; j<24; j++)
      for(int v=0; v<VLENS; v++)
        x[i].ccs[j].v[v] = b[i].ccs[j].v[v];

  ddd_out_pre_h_(     b, &domain_e);
  ddd_out_pos_h_( xe, b, &domain_e, (float)kappa);
  ddd_out_pre_h_(     b, &domain_o);
  ddd_out_pos_h_( xo, b, &domain_o, (float)kappa);
  assign_q_h2s( xs, x, &bsnorm); 

  for(int i=0; i<vols*2; i++)
    for(int j=0; j<24; j++)
      for(int v=0; v<VLENS; v++){
        float rr = xxs[i].ccs[j].v[v] - xs[i].ccs[j].v[v];
        if ( rr*rr > tol*tol  && rank == 0) 
          printf("CHKDEO: %3d %3d %3d %4d (SP)x = %14.6e (HP)x = %14.6e diff = %14.6e\n",
                  rank,v,j,i,xxs[i].ccs[j].v[v],
                              xs[i].ccs[j].v[v],rr);
      }
}

#endif
