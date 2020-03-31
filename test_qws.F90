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
!===============================================================================
!
! test_qws.F90
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 2014,2015 Yoshifumi Nakamura
!
! This file is part of BQCD -- Berlin Quantum ChromoDynamics program
!
! BQCD is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! BQCD is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with BQCD.  If not, see <http://www.gnu.org/licenses/>.
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine test_qws()
  call test_qws_eo()
!  call comm_finalize();  stop

  call test_qws_dd()
  call comm_finalize();  stop
end
!-------------------------------------------------------------------------------
subroutine test_qws_dd()
  use module_vol
  use module_action
  use module_field
  use module_lattice
  use module_input
  use module_function_decl
  implicit none
  !dir$ attributes align: 64:: ine, ino, tmpe, tmpo, rind, routd, rins, routs, ttts, ttts2, tttd
  real(8), dimension (24,volh_tot) :: ine, ino, tmpe, tmpo, oute, outo
  real(8), dimension (24*vol) :: rind, routd, tttd
  real(4), dimension (24*vol) :: rins, routs, ttts, ttts2
  integer :: i, j, cid, mid, iter, maxiter, maxiter_s, nsap, nm, block_size(4)
  real(8) :: kappa, tol, tol_s
  block_size =1

  mid = 1
  cid = action%mtilde(mid)%cid
  kappa = action%mtilde(mid)%kappa

  call qws_init(nx,ny,nz,nt, npe, bc_fermions, se, so, block_size)
  call flip_bc(gauge(2)%u) ! now gauge(2)%u is not flipped BC
  call qws_loadg_dd(gauge(2)%u(1,1,1,0,1), volh_tot, kappa)
  call flip_bc(gauge(2)%u) ! now flip again for fortran d()
  call qws_loadc_dd(clover(cid)%cfield%i(1,1,0), volh)

  call ran_gauss_volh(4 * 3, ine, 0.5_8, 0)!!; ine=0;ine(1,1)=1
  call ran_gauss_volh(4 * 3, ino, 0.5_8, 1)!!; ino=0
  call qws_loadfd_dd(rind, ine, ino)
  call qws_loadfs_dd(rins, ine, ino)

  call d(se, so, oute, ino, gauge(2)%u(1,1,1,0,1))
  call clover_inv_mult(oute, cid, se)
  call sc_xpby(oute, ine, - kappa)
  call d(so, se, outo, ine, gauge(2)%u(1,1,1,0,1))
  call clover_inv_mult(outo, cid, so)
  call sc_xpby(outo, ino, - kappa)


  call ddd_d(routd, rind)
  call qws_storefd_dd(tmpe, tmpo, routd)
  call compare_sc(oute, tmpe, 24*volh)
  call compare_sc(outo, tmpo, 24*volh)

  !!call comm_finalize()
  !!stop

  call ddd_s(routs, rins)
  call qws_storefs_dd(tmpe, tmpo, routs)
  call compare_sc(oute, tmpe, 24*volh)
  call compare_sc(outo, tmpo, 24*volh)

!------------------------------- static solver
  maxiter=100
  i = 0
  call jinv_ddd_in_s(ttts, rins, i, maxiter)
  call ddd_in_s(routs, ttts, i)
  i=1
  call jinv_ddd_in_s(ttts(1+24*volh), rins(1+24*volh), i, maxiter)
  call ddd_in_s(routs(1+24*volh), ttts(1+24*volh), i)

  call qws_storefs_dd(tmpe, tmpo, routs)
  call compare_sc(ine, tmpe, 24*volh)
  call compare_sc(ino, tmpo, 24*volh)

!------------------------------- bicgstab
  if (my_pe()==0)write(0,*)"bicgstab_dd_s check"
  maxiter=100
  ttts = rins
  call bicgstab_dd_s(ttts, rins, i, maxiter)
  call ddd_s(routs, ttts)
  call qws_storefs_dd(tmpe, tmpo, routs)
  call compare_sc(ine, tmpe, 24*volh)
  call compare_sc(ino, tmpo, 24*volh)
  if (my_pe()==0)write(0,*)"bicgstab_dd_s iter",i

!------------------------------- bicgstab
  if (my_pe()==0)write(0,*)"bicgstab_precdd_s check"
  maxiter=100
  ttts = rins
  nsap=1
  nm=1
  tol = 1d-14
  call bicgstab_precdd_s(ttts, rins, tol, i, maxiter, nsap,nm)
  call prec_ddd_s(routs, ttts,nsap,nm)
  call qws_storefs_dd(tmpe, tmpo, routs)
  call compare_sc(ine, tmpe, 24*volh)
  call compare_sc(ino, tmpo, 24*volh)
  if (my_pe()==0)write(0,*)"bicgstab_precdd_s iter",i

!------------------------------- bicgstab
  if (my_pe()==0)write(0,*)"bicgstab_precdd_s 2nd check"
  maxiter=100
  ttts = rins
  nsap=1
  nm=1
  call bicgstab_precdd_s(ttts, rins, tol, i, maxiter, nsap,nm)
  call prec_s(ttts2, ttts,nsap,nm)
  call ddd_s(routs, ttts2)
  call qws_storefs_dd(tmpe, tmpo, routs)
  call compare_sc(ine, tmpe, 24*volh)
  call compare_sc(ino, tmpo, 24*volh)
  if (my_pe()==0)write(0,*)"bicgstab_precdd_s iter",i

!------------------------------- bicgstab_dd_d
  if (my_pe()==0)write(0,*)"bicgstab_dd_d check"
  maxiter=100
  tttd = rind
  call bicgstab_dd_d(tttd, rind, i, maxiter)
  call ddd_d(routd, tttd)
  call qws_storefd_dd(tmpe, tmpo, routd)
  call compare_sc(ine, tmpe, 24*volh)
  call compare_sc(ino, tmpo, 24*volh)
  if (my_pe()==0)write(0,*)"bicgstab_dd_d iter",i

!------------------------------- bicgstab_dd_mix
  if (my_pe()==0)write(0,*)"bicgstab_dd_mix check"
  maxiter  = 10
  maxiter_s= 1000
  tttd = rind
  nsap=4!5
  nm=2
  tol  = 1d-14
  tol_s= 1d-6

  call bicgstab_dd_mix(tttd, rind, tol, i, maxiter, tol_s, maxiter_s, nsap, nm)
  call ddd_d(routd, tttd)
  call qws_storefd_dd(tmpe, tmpo, routd)
  call compare_sc(ine, tmpe, 24*volh)
  call compare_sc(ino, tmpo, 24*volh)
  if (my_pe()==0)write(0,*)"bicgstab_dd_mix iter",i

#ifdef HALF_PREC
!------------------------------- bicgstab_dd_mix2_hf
  if (my_pe()==0)write(0,*)"bicgstab_dd_mix2_hf check"
  maxiter  = 10
  maxiter_s= 1000
  tttd = rind
  nsap=4!5
  nm=2
  tol  = 1d-14
  tol_s= 1d-6

  call bqcd_qws_h_init()
  call bicgstab_dd_mix2_hf(tttd, rind, tol, i, maxiter, tol_s, maxiter_s, nsap, nm)
  call qws_storefd_dd(tmpe, tmpo, routd)
  call compare_sc(ine, tmpe, 24*volh)
  call compare_sc(ino, tmpo, 24*volh)
  if (my_pe()==0)write(0,*)"bicgstab_dd_mix2_hf iter",i
#endif

  call comm_finalize()
  stop

end

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine test_qws_eo()
  use module_vol
  use module_action
  use module_field
  use module_lattice
  use module_input
  use module_function_decl
  implicit none
  !dir$ attributes align: 64:: in, out, tmp, tmp1, rind, routd, rins, routs
  real(8), dimension (24*volh_tot) :: in, out, tmp, tmp1
  real(8), dimension (24*volh) :: rind, routd
  real(4), dimension (24*volh) :: rins, routs
  integer :: i, j, cid, mid, iter, maxiter, block_size(4)
  real(8) :: kappa, rtmp0, rtmp1


  integer, parameter :: nshift = 10
  real(8), dimension (24*volh_tot, nshift) :: xx
  real(8), dimension (24*volh, nshift) :: xxqws
  real(8) :: sigma(nshift)
  block_size = 1

  mid = 1
  cid = action%mtilde(mid)%cid
  kappa = action%mtilde(mid)%kappa

  call qws_init(nx,ny,nz,nt, npe, bc_fermions, se, so, block_size)
  call flip_bc(gauge(2)%u) ! now gauge(2)%u is not flipped BC
  call qws_loadg_bqcd(gauge(2)%u(1,1,1,0,1), volh_tot, kappa)
  call flip_bc(gauge(2)%u) ! now flip again for fortran d()
  call qws_loadc_bqcd(clover(cid)%cfield%i(1,1,0), volh)


  call ran_gauss_volh(4 * 3, in, 0.5_8, 0) !!;in=0;in(1)=1.0_8
  call qws_loadfd_bqcd(rind,in)
  call qws_loadfs_bqcd(rins,in)


  !!------------------------------------------------------- single prec.

  !! check d
  call d(se, so, out, in, gauge(2)%u(1,1,1,0,1))
  call deo_s(se, so, routs, rins)
  call qws_storefs_bqcd(tmp1, routs)
  call compare_sc(out, tmp1, 24*volh)

  !! check clv
  tmp=in
  routs=rins
  call clover_mult_ao(clover(cid)%cfield%i(1,1,se), tmp, volh)
  call clv_s(se, routs)
  call qws_storefs_bqcd(tmp1, routs)
  call compare_sc(tmp1, tmp, 24*volh)

  call d(se, so, out, in, gauge(2)%u(1,1,1,0,1))
  call clover_mult_ao(clover(cid)%cfield%i(1,1,se), out, volh)
  call dee_deo_s(se, so, routs, rins)
  call qws_storefs_bqcd(tmp1, routs)
  call compare_sc(out, tmp1, 24*volh)

  call wmul(out, in, mid)
  call mtilde_s(routs, rins)
  call qws_storefs_bqcd(tmp1, routs)
  call compare_sc(out, tmp1, 24*volh)

  !!call comm_finalize()
  !!stop

  !!------------------------------------------------------- double prec.
  !! check d
  call d(se, so, out, in, gauge(2)%u(1,1,1,0,1))
  call deo_vm(se, so, routd, rind)
  call qws_storefd_bqcd(tmp1, routd)
  call compare_sc(out, tmp1, 24*volh)

  !!call comm_finalize()
  !!stop

  !! check d dag
  call d_dag(se, so, out, in, gauge(2)%u(1,1,1,0,1))
  call deo_dag_vm(se, so, routd, rind)
  call qws_storefd_bqcd(tmp1, routd)
  call compare_sc(out, tmp1, 24*volh)

  !!call comm_finalize()
  !!stop

  !! check clv
  tmp=in
  routd=rind
  call clover_mult_ao(clover(cid)%cfield%i(1,1,se), tmp, volh)
  call clvd_vm(se, routd)
  call qws_storefd_bqcd(tmp1, routd)
  call compare_sc(tmp1, tmp, 24*volh)

  !!call comm_finalize()
  !!stop

  !! check wmul
  call wmul(out, in, mid)
  call mtilde_vm(routd, rind)
  call qws_storefd_bqcd(tmp1, routd)
  call compare_sc(out, tmp1, 24*volh)

  !! check wdag
  call wdag(out, in, mid)
  call mtilde_dag_vm(routd, rind)
  call qws_storefd_bqcd(tmp1, routd)
  call compare_sc(out, tmp1, 24*volh)

  !! check wdagw
  call wdagw(out, in, mid)
  call mtdagmt_vm(routd, rind)
  call qws_storefd_bqcd(tmp1, routd)
  call compare_sc(out, tmp1, 24*volh)

  !! check bicgstab
  out = in
  call qws_bc(out, in, mid)
  maxiter =101
  routd = rind
  call bicgstab_vm(routd, rind, iter, maxiter)
  call qws_storefd_bqcd(tmp1, routd)
  call compare_sc(out, tmp1, 24*volh)

  !! check cg
  out = in
  call qws_cg(out, in, mid)
  maxiter =101
  routd = rind
  call cg_vm(routd, rind, iter, maxiter)
  call qws_storefd_bqcd(tmp1, routd)
!  call compare_sc(out, tmp1, 24*volh)

  if (my_pe()==0)write(0,*)"test qws_mcg"
  !! check mcg
  sigma( 1)=1d-6
  sigma( 2)=1d-5
  sigma( 3)=1d-4
  sigma( 4)=1d-3
  sigma( 5)=1d-2
  sigma( 6)=1d-1
  sigma( 7)=1
  sigma( 8)=10
  sigma( 9)=20
  sigma(10)=50

  call qws_mcg(xx, in, mid, nshift, sigma)
  do j = 1, nshift
      call wdagw(out, xx(1, j), mid)
      do i = 1, 24*volh
         out(i) = out(i) + sigma(j) * xx(i, j) - in(i)
      enddo
      rtmp0 = sc_norm2(out)    ; call global_sum_vec(1, rtmp0)
      rtmp1 = sc_norm2(xx(1,j)); call global_sum_vec(1, rtmp1)
      if (my_pe()==0)write(0,*)"qws_mcg:", j, rtmp0, rtmp1
  enddo

  maxiter=101
  call mcg_vm(xxqws, rind, nshift, sigma, maxiter)
  do j = 1, nshift
      call qws_storefd_bqcd(tmp1, xxqws(1,j))
      call wdagw(out, tmp1, mid)
      do i = 1, 24*volh
         out(i) = out(i) + sigma(j) * tmp1(i) - in(i)
      enddo
      rtmp0 = sc_norm2(out)  ; call global_sum_vec(1, rtmp0)
      rtmp1 = sc_norm2(tmp1) ; call global_sum_vec(1, rtmp1)
      if (my_pe()==0)write(0,*)"mcg_vm:", j, rtmp0, rtmp1
  enddo

end

!-------------------------------------------------------------------------------
subroutine qws_bc(x, b, id)
  use module_vol
  implicit none
  complex(8), dimension (12*volh_tot) :: x, b ,r0, r, p, q, t
  complex(8) :: alpha, beta, omega, rho0, rho, ctmp
  integer :: id, iter, i
  real(8) :: rnorm, bnorm
  complex(8), external :: sc_cdotc
  real(8), external :: sc_norm2

  bnorm = sc_norm2(b) ;  call global_sum_vec(1, bnorm)

  call wmul(q, x, id)

  r = b - q
  r0 = r
  p = r
  rho0 = sc_cdotc(r0,r) ;  call global_sum_vec(2, rho0)

  do iter = 1, 100

     call wmul(q, p, id)
     ctmp = sc_cdotc(r0,q) ; call global_sum_vec(2, ctmp)
     alpha = rho0 / ctmp

     do i = 1, volh*12
        x(i) = x(i) + alpha * p(i)
        r(i) = r(i) - alpha * q(i)
     enddo

     rnorm = sc_norm2(r) ; call global_sum_vec(1, rnorm)
!     write(*,*)iter, rnorm
     if (sqrt(rnorm/bnorm)<1e-14) exit

     call wmul(t, r, id)
     ctmp = sc_cdotc(t,t) ; call global_sum_vec(2, ctmp)
     omega = sc_cdotc(t,r); call global_sum_vec(2, omega)
     omega = omega / ctmp

     do i = 1, volh*12
        x(i) = x(i) + omega * r(i)
        r(i) = r(i) - omega * t(i)
     enddo

     rnorm = sc_norm2(r) ; call global_sum_vec(1, rnorm)
!     write(*,*)iter, rnorm
     if (sqrt(rnorm/bnorm)<1e-14) exit

     rho = sc_cdotc(r0,r) ; call global_sum_vec(2, rho)
     beta = ( alpha / omega ) * ( rho / rho0 )
     rho0 = rho

     do i = 1, volh*12
        p(i) = p(i) - omega * q(i) 
     enddo
     do i = 1, volh*12
        p(i) = r(i) + beta * p(i) 
     enddo

  enddo
end

!-------------------------------------------------------------------------------
subroutine qws_cg(x, b, id)
  use module_vol
  implicit none
  real(8), dimension (24*volh_tot) :: x, b ,r, p, q
  real(8) :: alpha, beta, rho0, rho
  integer :: id, iter, i
  real(8) :: rnorm, bnorm, rtmp
  real(8), external :: sc_dot
  real(8), external :: sc_norm2

  bnorm = sc_norm2(b) ; call global_sum_vec(1, bnorm)
  call wdagw(q, x, id)

  r = b - q
  p = r
  rnorm= sc_norm2(r) ; call global_sum_vec(1, rnorm)
  rho0 = rnorm

  do iter = 1, 100
     call wdagw(q, p, id)
     rtmp = sc_dot(p,q) ; call global_sum_vec(1, rtmp)
     alpha = rho0 / rtmp

     do i = 1, volh*24
        x(i) = x(i) + alpha * p(i)
        r(i) = r(i) - alpha * q(i)
     enddo

     rnorm = sc_norm2(r) ; call global_sum_vec(1, rnorm)
!     write(*,*)iter, rnorm
     if (sqrt(rnorm/bnorm)<1e-14) exit

     rho = rnorm
     beta = rho / rho0
     rho0 = rho

     do i = 1, volh*24
        p(i) = r(i) + beta * p(i) 
     enddo

  enddo
end


!-------------------------------------------------------------------------------
subroutine qws_mcg(xx, b, id, n, sigma) 
  use module_vol
  implicit none
  integer :: n, id, i, j, iter
  real(8) :: sigma(n), sxi0(n), sxi1(n), sxi2(n)
  real(8) :: alpha, alpha0, beta, rho, rho0, bnorm, rnorm, rtmp
  real(8), dimension (24*volh_tot,n) :: xx, pp
  real(8), dimension (24*volh_tot) :: r, p, q, b
  real(8), external :: sc_dot, sc_norm2

  alpha0=1
  beta  =0
  sxi0  =1
  sxi1  =1

  r = b
  p = r
  rho0  = sc_norm2(r) ; call global_sum_vec(1, rho0)
  bnorm = rho0

  do i=1,n
     call sc_copy(pp(1,i),r)
     call sc_zero(xx(1,i))
  enddo

  do iter = 1, 100
     call wdagw(q, p, id)
     rtmp =sc_dot(p, q) ; call global_sum_vec(1, rtmp)
     alpha= - rho0 / rtmp

     do i = 1, volh*24
        r(i) = r(i) + alpha * q(i)
     enddo

     rnorm = sc_norm2(r) ; call global_sum_vec(1, rnorm)

     do j = 1,n
!!        if ( n_con(j) == 0) then
           sxi2(j)=sxi1(j)*sxi0(j)*alpha0/(alpha*beta*(sxi0(j)-sxi1(j))+sxi0(j)*alpha0*(1-sigma(j)*alpha))
           rtmp = alpha*sxi2(j)/sxi1(j)
           do i = 1, 24*volh
              xx(i,j) = xx(i,j) - rtmp * pp(i,j)
           enddo
!!           if (sqrt(sxi2(j)*sxi1(j)*rnorm/bnorm)< 1e-14) niter_shift(j)= niter
!!           if (sqrt(sxi2(j)*sxi1(j)*rnorm/bnorm)< 1e-14) n_con(j)=1
!!        endif
     enddo

     if (sqrt(sxi2(1)*sxi1(1)*rnorm/bnorm)< 1e-14) exit

     rho  = rnorm
     beta = rho / rho0
     rho0 = rho

     do i = 1, volh*24
        p(i) = r(i) + beta * p(i)
     enddo

     do j = 1,n
!!        if ( n_con(j) == 0) then
           rtmp = beta * (sxi2(j)/sxi1(j))**2
           do i = 1, 24*volh
              pp(i, j) = rtmp*pp(i, j) + sxi2(j)*r(i)
           enddo
           sxi0(j) = sxi1(j)
           sxi1(j) = sxi2(j)
!!        endif
     enddo

     alpha0 = alpha
  enddo

end
