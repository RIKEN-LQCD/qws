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
qws0.1.5 tune20160920内容

1.Makefile修正内容
 ・手動プリフェッチをストロングプリフェッチになるよう修正
 ・-Knoprefetch,noswpの追加

2.bicgstab_dd_mix.cc修正内容
 ・配列のアラインメント調整

3.bicgstab_precdd_s.ccの修正内容
 ・配列のアラインメント調整
 ・手動プリフェッチの挿入
 
4.ddd_in_s.ccの修正内容
 ・ループ分割
 ・スレッドインバランス解消
 ・手動プリフェッチ挿入
 ・parallelリージョン拡大のための修正
 
5.ddd_out_s_0.ccの修正内容
 ・区間ddd_out_pre_s_のタスク並列化
 ・手動プリフェッチ挿入
 ・区間ddd_out_pos_s_のスレッドインバランス解消
 ・parallelリージョン拡大のための修正
 
6.mult_all.hの修正内容
 ・最適化促進のためのpragma挿入
 ・作業配列削減のための修正

7.no_timer_copy.shの修正内容
 ・parallelリージョン拡大のための修正
 
8.profiler.hの修正内容
 ・parallelリージョン拡大のための修正
 
9.qws.ccの修正内容
 ・parallelリージョン拡大のための修正
 ・プリフェッチ挿入
 ・配列のアラインメント調整
 
10.static_solver.ccの修正内容
 ・parallelリージョン拡大のための修正
 ・プリフェッチ挿入
 
11.time.hの修正内容
 ・parallelリージョン拡大のための修正