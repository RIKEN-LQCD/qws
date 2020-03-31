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
qws-0.1.8 tune20181009 変更内容
2019/4/23
富士通株式会社

* qws-0.1.8 からの主な変更点
  1. ddd_in_s, ddd_out_pre_s, ddd_out_pos_s のチューニング
     - スレッド間負荷均衡化のためにループ分割したりスレッドの担当座標をループに
       よって変更したりしていたものを、単純な1つの4重ループに戻した。
       (評価対象の条件では、負荷不均衡の不利益より、整数演算処理の単純化とメモリアクセス
       の局所化による利益の方が大きい。)
     - 1つのループにした結果、ループの次のイテレーションに対するプリフェッチでは
       キャッシュから溢れるため、同イテレーション内の次の方向の計算に対してプリフェッチ
       するよう変更。
     - qwsintrin.h に定義されたベクトルイントリンジック向けのインターフェースについて、
       - SVE 向けの実装を追加、
       - マスクの引数を追加、
       - マスク演算、符号反転、論理和、
       - ロード/ストアのアドレスをベースとオフセットに分けて指定できるよう引数を追加。
     - ddd_in_s の X 方向の計算時の不連続アクセスをベクトルイントリンジックを使って
       効率化。
     - 統一化のため、他の処理もベクトルイントリンジックを使うよう変更。(ただしその結果、
       SVE 以外は汎用のループによる実装しかないため性能低下)
     - 整数レジスタ使用数削減のため、
       - 変数が書き換えられたかのように見せかけることでレジスタ割り付けを調整、
       - mult_clvs をインライン化、
       - 問題サイズをコンパイル時定数化。
     - mult_all の計算も mult_clvs のように命令レベルの理想的なスケジューリングに
       合うようソースレベルで並べ替え。
     - ロード/ストアのアドレッシングを効率化するため、命令レベルの理想的な即値を
       イントリンジックの即値に相当する引数にコンパイル時定数として指定。
       (実行時間的には効果が小さい可能性あり)
     - ddd_out_s の後方向の計算時に、 scs_t 型の一時領域を常にレジスタに置く。
       (実行時間的には効果が小さい可能性あり)
  2. その他のチューニング、性能のための修正
     - bicgstab_precdd_s のリダクションのループを、ベクトルイントリンジック関数
       を使ってベクトルのレーン間のリダクションの回数が少なくなるようチューニング。
     - prec_ddd_s_, static_solver の冗長なスレッド間バリアの削除。
     - glus, clvs のキャッシュラインサイズへのアライン漏れを修正。
  3. 測定区間の新しい分割を追加。測定のオーバーヘッドを抑えるため、通信とオーバーラップ
     しない区間(all_calc)とする区間(overlapped)の2つに分割。
  4. カーネル化して bicgstab_precdd_s のみ実行する際、通常の実行時に保存しておいた
     初期値を使用し、計算結果を比較するよう修正。
     -DDUMP_DATA_FOR_KERNELIZE で通常実行することで計算結果とアセンブリ形式の
     初期値が表示される。アセンブルしてカーネル実行時にリンクして使用する。
     いくつかのの条件で得られた初期値を data/data_* に置いた。
