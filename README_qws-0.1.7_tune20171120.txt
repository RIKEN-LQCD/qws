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
qws-0.1.7 tune20171120 変更内容
2017/11/20
富士通株式会社

* qws-0.1.7 からの主な変更内容
  1. シミュレーションのためのカーネル化とイタレーション縮小用コードの追加
     カーネル化、イタレーション縮小は CA9-170718-F-1-progress.pdf 、
     結果確認は CA9-170314-F-2-simulation.pdf の方法で実施。
  2. ポスト京向けチューニング
     mult_clvs() のチューニング(CA9-170523-F-1-progress.pdf)と、
     VLENS=16 用のプリフェッチ(CA9-170314-F-1-kernel_in.tar.gz)を実施。
     その他細かいチューニングを実施。
  3. リファクタリング
     今回は、ほとんど性能影響がないことが明らかかつ簡単な範囲で実施。

* make 方法の変更内容
  1. 変数 kernelize を追加。定義されているとき、カーネル化かつイタレーション縮小(1/500)した
     実行ファイルを生成する。(例: $ make kernelize=1)
  2. 変数 vlens を追加。 VLENS を指定する。デフォルトは 8。
  3. 変数 tune を追加。 tune=1 のとき、ポスト京向けチューニングコードを使用する。
  4. 変数 strong_prefetch を追加。strong_prefetch を定義したときプリフェッチをストロングに置き換える。
     デフォルトは無効。(SPARC のみ有効)
  5. 変数 arch を追加。指定アーキ向けのフラグの設定などを行う。デフォルトは FX100。
     arch=fx100: vlens=8, プリフェッチのストロングへの置き換えを行う。
     arch=postk: vlens=16, tune=1
  例:
  $ make arch=fx100 target=jinv             # アプリ全体推定用の設定
  $ make kernelize=1 arch=postk target=jinv # ポスト京シミュレーション向けの設定
  $ make kernelize=1 tune=1 vlens=16 \      # ポスト京シミュレーション結果と比較するための、
         strong_prefetch=1 target=jinv      # 性能予測ツール用の FX100 性能情報取得向けの設定

* ファイル単位の変更内容
  * 追加
    * main_kern.cc
      カーネル化向けの main() が書かれたファイル。(他は、同じファイルの中で ifdef で切り替え)
    * eml_lib.h, eml_lib.cc
      シミュレーション向けの、使用している標準ライブラリ関数を置き換える簡易なエミュレーション版。
    * report.h, report.c, init.hh, io.h
      カーネル化時の結果確認用の関数。初期値用の乱数生成、二乗和の計算、誤差の表示など。
    * ddd_out_s.cc, ddd_out_s_inline.h
      タイマーあり/なし版を作る際に、今まで no_timer_copy.sh で関数を複製していた処理を止めて、
      テンプレートで異なるインスタンスにするよう修正。
      (なお、あり/なし版で別のインスタンスにするのは、ラインプロファイルなどを取得する際に区別しやすくするため)
  * 削除
    * ddd_out_s_0.cc, ddd_out_s_0_inline.h, no_timer_copy.sh
      上記の通り修正したので削除、ファイル名を分かりやすくするためにリネーム。
  * 変更
    * Makefile
      上記の make 方法の変更のための修正。
      ストロングプリフェッチへの置き換え、ddd_out_s.cc 関連のリファクタリング。
    * bicgstab_precdd_s.cc, bicgstab_dd_mix.cc, main.c, qws.cc, qws.h
      カーネル化向けの修正、リファクタリング。
    * clover_s.cc
      ポスト京向けチューニングの追加。
    * ddd_in_s.cc
      VLENS=16 向けのプリフェッチ追加、プリフェッチの修正。
    * profiler.h, time.h
      カーネル化時の測定のための修正。(測定の前にウォームアップのために一度実行するため)
      PROF_START/STOP のインタフェースに区間名を指定できるよう修正。(他の、区間名指定可能なプロファイラと合わせるため)
    * static_solver.cc
      グローバル変数の再読み込みを防ぐチューニング。

* 実行結果
  * 32x3x4x4 の test_double/single_prec_functions が次の条件で CASE2 と比較してそれぞれ
    12/6 桁目の +-1 の範囲に収まることを確認。
    * FX100 2 proc x 12 thrd, vlens=8
    * FX100 2 proc x 12 thrd, vlens=16
    * FX100 2 proc x 12 thrd, vlens=8, tune=1
    * FX100 2 proc x 12 thrd, vlens=16, tune=1

* 実行時間
  * 32x6x6x2 の各区間の FX100 実行時間が以下の条件で大きく変わらないことを確認
    * qws-0.1.7 オリジナル (CA9-170228-R-2-qws-0.1.7-PA/32x6x6x2/run_32x6x6x2.sh.8503.*/rank0/output_prof_1.csv)
    * tune20171120 通常実行 (arch=fx100)
      この条件では本質的な変更を行っていないため差がないことを確認
    * tune20171120 カーネル化 (arch=fx100, kernelize=1)
      カーネル化によって実行時間が変わっていないことを確認
  * 詳細な結果(1イタレーションあたり, マイクロ秒, ブレにより多少の差はあり)
    条件                     in   jinv  pre  pos  other
    ---------------------------------------------------
    qws-0.1.7_オリジナル     275  625   60   184  96
    tune20171120_通常実行    279  634   63   190  93
    tune20171120_カーネル化  275  633   60   184  90

以上
