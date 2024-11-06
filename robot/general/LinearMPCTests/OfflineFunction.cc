//! @file OfflineFunction.cc
//! @brief ARCS6 オフライン計算用メインコード
//! @date 2024/06/25
//! @author Yokokura, Yuki
//!
//! @par オフライン計算用のメインコード
//! - 「make offline」でコンパイルすると，いつものARCS制御用のコードは走らずに，
//!    このソースコードのみが走るようになる。
//! - ARCSライブラリはもちろんそのままいつも通り使用可能。
//! - 従って，オフラインで何か計算をしたいときに，このソースコードに記述すれば良い。
//!
// Copyright (C) 2011-2024 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

// 基本のインクルードファイル
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <complex>
#include <array>

// 追加のARCSライブラリをここに記述
#include "ArcsMatrix.hh"
#include "CsvManipulator.hh"
#include "LinearMPC.hh"

using namespace ARCS;



//! @brief エントリポイント
//! @return 終了ステータス
int main(void){
	printf("ARCS OFFLINE CALCULATION MODE\n");
	printf("Linear MPC tests \n");

	ArcsMat<2,2> A;
	ArcsMat<2,1> B;

	LinearMPC<2,1,5,3,true,0> mpcNoConstraintsTest(A,B);
	LinearMPC<2,1,5,3,true,2> mpcWithConstraintsTest(A,B);
	LinearMPC<2,1,5,3,true,0> mpcNoConstraintsTestAgain(A,B);





	// ここにオフライン計算のコードを記述
	
	return EXIT_SUCCESS;	// 正常終了
}




