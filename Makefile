##****************************************************************************************
##
##  Copyright (c) 2015-2020, Yoshifumi Nakamura <nakamura@riken.jp>
##  Copyright (c) 2015-2020, Yuta Mukai         <mukai.yuta@fujitsu.com>
##  Copyright (c) 2018-2020, Ken-Ichi Ishikawa  <ishikawa@theo.phys.sci.hirosima-u.ac.jp>
##  Copyright (c) 2019-2020, Issaku Kanamori    <kanamori-i@riken.jp>
##
##
##  All rights reserved.
##
##  Redistribution and use in source and binary forms, with or without
##  modification, are permitted provided that the following conditions are
##  met:
##
##  * Redistributions of source code must retain the above copyright
##    notice, this list of conditions and the following disclaimer.
##
##  * Redistributions in binary form must reproduce the above copyright
##    notice, this list of conditions and the following disclaimer listed
##    in this license in the documentation and/or other materials
##    provided with the distribution.
##
##  * Neither the name of the copyright holders nor the names of its
##    contributors may be used to endorse or promote products derived from
##    this software without specific prior written permission.
##
##  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
##  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
##  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
##  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
##  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
##  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
##  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
##  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
##  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
##  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
##  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##
##----------------------------------------------------------------------------------------
##  ACKNOWLEDGMENT
##
##  This software has been developed in a co-design working group for the lattice QCD
##  supported by MEXT's programs for the Development and Improvement for the Next
##  Generation Ultra High-Speed Computer System, under its Subsidies for Operating the
##  Specific Advanced Large Research Facilities, and Priority Issue 9
##  (Elucidation of the Fundamental Laws and Evolution of the Universe) to be tackled by
##  using the Supercomputer Fugaku.
##
##****************************************************************************************
debug     =
mpi       =1
omp       =1
#compiler [openmpi-gnu, gnu, intel, fujitsu_native, fujitsu_cross]
compiler  =fujitsu_cross
#arch [fx100, postk, skylake, ofp, thunderx2, simulator]
arch      =postk
#profiler [timing, fapp, fpcoll, pa] (nondisclousure: timing2)
profiler  =timing
#timing2_path=$(HOME)/opt/timing2.o
prof_selective=
#target [jinv, in, pre, pos, other, all_calc, overlapped, send, send_post]
#target    =all_calc
half_prec =
#path to half precision library required in non clang mode
#libhalf=$(HOME)/opt/half-1.12.0/include
PREFIX    =.
#rdma [,utofu, utofu_threaded, mpi_rankmap]  (todo: fjmpi)
rdma      =
#clangmode : clang mode for fujitsu compiler
clang     =
#Power API : Power API for Fugaku and FX1000
powerapi  =1
#===============================================================================
FPP       = cpp -C -P --sysroot=.
ARFLAGS   = 
#FFLAGS    = -I$(HOME)/src/bqcdi/modules 
CPPFLAGS  = $(cppflags)

PATH_LIBREADPIC=$(HOME)


ifeq ($(compiler),openmpi-gnu)
  ifdef mpi
    CC        = mpicc
    CXX       = mpicxx
    F90       = mpif90
    MYFLAGS   = -D_MPI_
  else
    CC        = gcc
    CXX       = g++
    F90       = gfortran
  endif
  MYFLAGS    += -D_CHECK_TIMING
  CFLAGS      = -g
  CXXFLAGS    = $(CFLAGS) -std=gnu++11
  CXXFLAGS_A  = $(CFLAGS) -std=gnu++11
  LDFLAGS=
  SYSLIBS=
  ifdef omp
    CFLAGS     += -fopenmp
    CXXFLAGS   += -fopenmp
    CXXFLAGS_A += -fopenmp
  endif
else ifeq ($(compiler),fujitsu_native)
  ifdef mpi
    CC        = mpicc
    CXX       = mpicxx
    F90       = mpif90
    MYFLAGS   = -D_MPI_
  else
    CC        = fcc
    CXX       = FCC
    F90       = frt
  endif
  CFLAGS      = -Xg -Kfast,restp=all,optmsg=2,ocl,preex,noeval,noprefetch,noswp -Nline,lst=t
  CXXFLAGS    = $(CFLAGS) -std=gnu++11
  CXXFLAGS_A  = $(CFLAGS) -std=gnu++11
  LDFLAGS     =
  SYSLIBS     =
  ifdef omp
    CFLAGS   += -Kopenmp
    CXXFLAGS += -Kopenmp
    CXXFLAGS_A += -Kopenmp
  endif
else ifeq ($(compiler),intel)
  ifdef mpi
    CC        = mpiicc
    CXX       = mpiicpc
    F90       = mpiifort
    MYFLAGS   = -D_MPI_
    MYFLAGS  += -D_NO_OMP_SINGLE  # error avoiding
  else
    CC        = icc
    CXX       = icpc
    F90       = ifort
  endif
  ifdef omp
    CC       += -qopenmp
    CXX      += -qopenmp
  endif
  MYFLAGS    += -D_CHECK_TIMING
  ifeq ($(arch),skylake)
     CFLAGS   = -O3 -xCORE-AVX512 -fno-alias -qopt-zmm-usage=high -Wno-unknown-pragmas
  else ifeq ($(arch),ofp)
     CFLAGS   = -O3 -xKNL -fno-alias -qopt-zmm-usage=high -Wno-unknown-pragmas
  endif
  CXXFLAGS    = $(CFLAGS) -std=gnu++11
  CXXFLAGS_A  = $(CFLAGS) -std=gnu++11
  CFLAGS     += -std=gnu99
  LDFLAGS     =
  SYSLIBS     =
else ifeq ($(compiler),fujitsu_cross)
  ifdef mpi
    CC        = mpifccpx
    CXX       = mpiFCCpx
    F90       = mpifrtpx
    RDMA_FLAGS=
    ifeq ($(rdma),utofu)
        RDMA_FLAGS= -D_RDMA_ -D_UTOFU_RDMA -D_UTOFU_RANKMAP
    endif
    ifeq ($(rdma),utofu_threaded)
        RDMA_FLAGS= -D_RDMA_ -D_UTOFU_RDMA -D_UTOFU_RANKMAP -D_THREADED_RDMA
    endif
    ifeq ($(rdma),mpi_rankmap)
        RDMA_FLAGS= -D_USE_RANKMAP -D_RDMA_ -D_UTOFU_RANKMAP
    endif
    ifeq ($(rdma),fjmpi)
        RDMA_FLAGS= -D_RDMA_ -D_HPC_RDMA
    endif
    MYFLAGS   = -D_MPI_ $(RDMA_FLAGS)
#    MYFLAGS   = -D_MPI_ -D_RDMA_ -D_HPC_RDMA  # with RDMA
#    MYFLAGS   = -D_MPI_  # without RDMA

    MYFLAGS    += -D_NO_OMP_SINGLE  # error avoiding
  else
    CC        = fccpx
    CXX       = FCCpx
    F90       = frtpx
  endif
  ifndef clang
    CFLAGS      = -Kfast,restp=all,optmsg=2,ocl,preex,noprefetch,noswp -Nline,lst=t
    CFLAGS      = -Kfast,restp=all,optmsg=2,ocl,preex,noprefetch,noswp -Nline,lst=t -Nfjomplib -Kilfunc=loop -Krdconv=2
    CFLAGS      = -Kfast,restp=all,optmsg=2,ocl,preex,noprefetch,noswp -Nline,lst=t -Nnofjprof -Nfjomplib -Kilfunc=loop -Krdconv=2
    CXXFLAGS    = $(CFLAGS) -std=gnu++11
    CXXFLAGS_A  = $(CFLAGS) -std=gnu++11 -Knosch_pre_ra,nosch_post_ra -Knoeval
  else
    ## clang mode
    CFLAGS      = -Nclang -mcpu=a64fx+sve -Ofast  -Koptmsg=2,preex,noprefetch -Nnofjprof  -Kilfunc=loop -ffj-line -ffj-lst=t 
    CXXFLAGS    = $(CFLAGS) -stdlib=libc++
    CXXFLAGS_A  = $(CFLAGS) -stdlib=libc++
  endif
  LDFLAGS     =
  SYSLIBS     = -ltofucom
  ifdef powerapi
    MYFLAGS    += -D_POWER_API_
    CFLAGS     += -I/opt/FJSVtcs/pwrm/aarch64/include
    CFLAGS_A   += -I/opt/FJSVtcs/pwrm/aarch64/include
    CXXFLAGS   += -I/opt/FJSVtcs/pwrm/aarch64/include
    CXXFLAGS_A += -I/opt/FJSVtcs/pwrm/aarch64/include
    SYSLIBS    += -L/opt/FJSVtcs/pwrm/aarch64/lib64 -lpwr
  endif
  ifdef omp
    CC       += -Kopenmp
    CXX      += -Kopenmp
  endif
else ifeq ($(compiler),gnu)
    CC        = gcc
    CXX       = g++
    F90       = gfortran
  ifdef omp
    CC       += -fopenmp
    CXX      += -fopenmp
  endif
  MYFLAGS    += -D_CHECK_TIMING -D_MPI_
  CFLAGS      = -g
  CXXFLAGS    = $(CFLAGS) -std=gnu++11
  CXXFLAGS_A  = $(CFLAGS) -std=gnu++11
  LDFLAGS=
  SYSLIBS=
else
  $(error Unknown compiler: $(compiler))
endif

ifeq ($(arch),fx100)
       vlend = 4
       vlens = 8
       strong_prefetch = 1
else ifeq ($(arch),postk)
       vlend = 8
       vlens = 16
       ifndef clang
       MYFLAGS += -DARCH_POSTK
# -DCOMPILE_TIME_DIM_SIZE -DNX=32 -DNY=6 -DNZ=4 -DNT=3 -DINLINE_ASM_UNAVAILABLE
       else
       MYFLAGS +=              -DCOMPILE_TIME_DIM_SIZE -DNX=32 -DNY=6 -DNZ=4 -DNT=3 -DINLINE_ASM_UNAVAILABLE
       endif
else ifeq ($(arch),thunderx2)
       vlend=2
       vlens=4
else  ifeq ($(arch),skylake)
       vlend=8
       vlens=16
else ifeq ($(arch),ofp)
       vlend=8
       vlens=16
else ifeq ($(arch),simulator)
       vlend=2
       vlens=2
else
  $(error Unknown architecture: $(arch))
endif

ifdef vlens
  MYFLAGS += -DVLENS=$(vlens)
endif
ifdef vlend
  MYFLAGS += -DVLEND=$(vlend)
endif

.SUFFIXES:
.SUFFIXES: .a .o .cc .c .F90 .s .S

.cc.o:
	$(CXX) $(CXXFLAGS) $(MYFLAGS) $(CPPFLAGS) -c -o $@ $< 
.c.o:
	$(CC) $(CFLAGS) $(MYFLAGS) $(CPPFLAGS) -c -o $@ $<

.cc.s:
	$(CXX) $(CXXFLAGS) $(MYFLAGS) $(CPPFLAGS) -S $< 

.F90.o:
	$(FPP) $(CPPFLAGS) $(MYFLAGS) $< > $*.f90
	$(F90) -c $(FFLAGS) $*.f90

.S.o:
	$(CC) $(CPPFLAGS) $< -c

# communications
ifdef mpi
  ifndef rdma
    OBJS_XBOUND = qws_xbound_mpi.o
    OBJS_COMM = $(OBJS_XBOUND)
  endif
  ifeq ($(rdma),utofu)
    OBJS_XBOUND = qws_xbound_rdma.o
    OBJS_COMM = $(OBJS_XBOUND) rdma_utofu_comlib.o rankmap_lib_utofu.o get_tni_4d.o get_tofu_coord.o get_tofu_coord_4d.o rdma_comlib_2buf.o
  endif
  ifeq ($(rdma),utofu_threaded)
    OBJS_XBOUND = qws_xbound_rdma.o
    OBJS_COMM = $(OBJS_XBOUND) rdma_utofu_comlib.o rankmap_lib_utofu.o get_tni_4d.o get_tofu_coord.o get_tofu_coord_4d.o rdma_comlib_2buf.o
  endif
  ifeq ($(rdma),mpi_rankmap)
    OBJS_XBOUND = qws_xbound_mpi.o
    OBJS_COMM = $(OBJS_XBOUND) rdma_utofu_comlib.o rankmap_lib_utofu.o get_tni_4d.o get_tofu_coord.o get_tofu_coord_4d.o rdma_comlib_2buf.o
  endif
  ifeq ($(rdma),fjmpi)
    OBJS_XBOUND = qws_xbound_rdma.o
    OBJS_COMM = $(OBJS_XBOUND) rdma_comlib.o rdma_comlib_2buf.o
  endif
else
  OBJS_XBOUND = qws_xbound_nompi.o
  OBJS_COMM = $(OBJS_XBOUND)
endif


OBJS := \
	timing.o \
	qws.o \
	bicgstab.o \
	deo_in_d.o \
	deo_in_s.o \
	deo_out_d.o \
	deo_out_s.o \
	deo_dag_in.o \
	deo_dag_out.o \
	cg.o \
	mcg.o \
	ddd_in_d.o \
	ddd_in_s.o \
	ddd_out_d.o \
	ddd_out_s.o \
	static_solver.o \
	bicgstab_dd_s.o \
	bicgstab_dd_d.o \
	bicgstab_precdd_s.o \
	bicgstab_dd_mix.o \
	reordering.o \
	clover_s.o \
	util.o \

MAIN = main.o

ifdef kernelize
  OBJS := \
    report.o \
    eml_lib.o \
    util.o \
    data.o \
  MAIN = main_kern.o
endif

OBJS += $(OBJS_COMM)

ifdef strong_prefetch
  objs_with_prefetch = qws.o ddd_in_s.o ddd_out_s.o static_solver.o bicgstab_precdd_s.o
  OBJS := $(filter-out $(objs_with_prefetch),$(OBJS))
  OBJS += $(objs_with_prefetch:.o=_strp.s)
endif

ifeq ($(target),jinv)
  CFLAGS += -DPROF_TARGET=TARGET_JINV
endif
ifeq ($(target),in)
  CFLAGS += -DPROF_TARGET=TARGET_IN
endif
ifeq ($(target),pre)
  CFLAGS += -DPROF_TARGET=TARGET_PRE
endif
ifeq ($(target),pos)
  CFLAGS += -DPROF_TARGET=TARGET_POS
endif
ifeq ($(target),other)
  CFLAGS += -DPROF_TARGET=TARGET_OTHER
endif
ifeq ($(target),all_calc)
  CFLAGS += -DPROF_TARGET=TARGET_ALL_CALC
endif
ifeq ($(target),overlapped)
  CFLAGS += -DPROF_TARGET=TARGET_OVERLAPPED
endif
ifeq ($(target),send)
  CFLAGS += -DPROF_TARGET=TARGET_SEND
endif
ifeq ($(target),send_post)
  CFLAGS += -DPROF_TARGET=TARGET_SEND_POST
endif

ifeq ($(profiler),timing)
  CFLAGS += -D_CHECK_TIMING
endif
ifeq ($(profiler),timing2)
  CFLAGS += -D_CHECK_TIMING2
  OBJS   += $(timing2_path)
  CFLAGS += -I.
endif
ifeq ($(profiler),fpcoll)
  CFLAGS += -D_FPCOLL
endif
ifeq ($(profiler),fapp)
  CFLAGS += -D_FAPP -Nlinkprof
endif
ifeq ($(profiler),pa)
  CFLAGS += -D_CHECK_PA -I$(PATH_LIBREADPIC)/libreadpic
  LDFLAGS += -lreadpic -L$(PATH_LIBREADPIC)/libreadpic/sparc64xifx
endif

ifeq ($(prof_selective),1)
  CFLAGS += -DPROF_SELECTIVE
endif


ifdef debug
  MYFLAGS    += -D_DEBUG_
endif

ifdef half_prec
  OBJS += \
	bicgstab_precdd_s_h.o \
	bicgstab_dd_mix2_hf.o \
	qws_h.o
  MYFLAGS  += -DHALF_PREC
ifdef libhalf
  MYFLAGS  += -I$(libhalf) -Dlibhalf
endif
endif


ifdef FFLAGS
fast:
	make main
	make libqws.a
else
fast:
	make main
endif

DELIVERABLES=main libqws.a

install: $(DELIVERABLES)
	install $(DELIVERABLES) $(PREFIX)

#prepare:
#	sed -i 's/^\#define VLEND.*/#define VLEND $(vlend)/' ./qws.h
#	sed -i 's/^\#define VLENS.*/#define VLENS $(vlens)/' ./qws.h

main:$(OBJS) $(MAIN) $(LDFLAGS)
	$(CXX) -o $@ $(MAIN) $(OBJS) $(SYSLIBS) $(CXXFLAGS) $(LDFLAGS)

ifdef half_prec
libqws.a:$(OBJS) test_qws.o qws_bqcd.o
	$(AR) $(ARFLAGS) rv $@ $(OBJS) test_qws.o qws_bqcd.o
else
libqws.a:$(OBJS) test_qws.o qws_bqcd.o
	$(AR) $(ARFLAGS) rv $@ $(OBJS) test_qws.o qws_bqcd.o
endif

clean:
	rm -f *.s *.o *.f90 *.lst
	rm -f main libqws.a

qws.o: qws.h clover_d.h clover_s.h clover_def.h
qws_xbound_nompi.o: qws.h
qws_xbound_mpi.o: qws.h
qws_xbound_rdma.o: qws.h rdma_utofu_comlib.h rdma_comlib_2buf.h get_tni.h rankmap_lib_utofu.o get_tni_4d.o get_tofu_coord.o get_tofu_coord_4d.o
qws_bqcd.o: qws.h
bicgstab.o: qws.h
bicgstab_dd_s.o: qws.h
bicgstab_dd_d.o: qws.h
bicgstab_dd_mix.o: qws.h
bicgstab_precdd_s.o: qws.h
cg.o:qws.h
mcg.o:qws.h
static_solver.o:qws.h
reordering.o:qws.h

bicgstab_dd_mix2_hf.o : qws.h
bicgstab_precdd_s_h.o : qws.h
qws_h.o: qws_h.cc qws.h clover_h.h  mult_all.h xbound_h.h ddd_in_h.h ddd_out_h.h ddd_out_h_0_inline.h  jinv_ddd_in_h.h prec_ddd_s_h.h

deo_in_d.o:qws.h wilson_def.h wilson_d.h clover_d.h
deo_in_s.o:qws.h wilson_def.h wilson_s.h clover_s.h
deo_out_d.o:qws.h wilson_def.h wilson_d.h clover_d.h
deo_out_s.o:qws.h wilson_def.h wilson_s.h clover_s.h
deo_dag_in.o:qws.h wilson_def.h wilson_d.h clover_d.h
deo_dag_out.o:qws.h wilson_def.h wilson_d.h clover_d.h

ddd_in_d.o:qws.h wilson_def.h wilson_d.h clover_d.h 
ddd_in_s.o:qws.h wilson_def.h wilson_s.h clover_s.h mult_all.h prefetch.h
	$(CXX) $(CXXFLAGS_A) $(MYFLAGS) $(CPPFLAGS) -c -o ddd_in_s.o ddd_in_s.cc
ddd_out_d.o:qws.h wilson_def.h wilson_d.h clover_d.h
ddd_out_s.o:qws.h wilson_def.h wilson_s.h clover_s.h mult_all.h ddd_out_s_inline.h
	$(CXX) $(CXXFLAGS_A) $(MYFLAGS) $(CPPFLAGS) -c -o ddd_out_s.o ddd_out_s.cc
clover_s.o:qws.h wilson_def.h wilson_s.h clover_s.h mult_all.h prefetch.h
	$(CXX) $(CXXFLAGS_A) $(MYFLAGS) $(CPPFLAGS) -c -o clover_s.o clover_s.cc
ddd_in_s.s:qws.h wilson_def.h wilson_s.h clover_s.h mult_all.h prefetch.h
clover_s.s:qws.h wilson_def.h wilson_s.h clover_s.h mult_all.h prefetch.h
clover_s2.s:qws.h wilson_def.h wilson_s.h clover_s.h mult_all.h prefetch.h

ddd_out_s.s:qws.h wilson_def.h wilson_s.h clover_s.h mult_all.h ddd_out_s_inline.h
qws.s: qws.h clover_d.h clover_s.h clover_def.h
static_solver.s:qws.h
bicgstab_precdd_s.s: qws.h
rdma_comlib.o: rdma_comlib.h
rdma_utofu_comlib.o: rdma_utofu_comlib.h
rdma_comlib_2buf.o: rdma_utofu_comlib.h rdma_comlib_2buf.h

%_strp.s: %.s
	sed -e "s/\(prefetch\t.*,\)0$$/\120/g" -e "s/\(prefetch\t.*,\)2$$/\122/g" $< > $@
