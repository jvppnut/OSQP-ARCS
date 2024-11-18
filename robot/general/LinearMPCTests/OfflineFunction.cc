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
using namespace ArcsMatrix;

//! @brief エントリポイント
//! @return 終了ステータス
int main(void){
	printf("ARCS OFFLINE CALCULATION MODE\n");
	printf("Linear MPC tests \n");

	constexpr size_t p_hor = 5;
	constexpr size_t c_hor = 3;
	constexpr size_t n_states = 2;
	constexpr size_t m_inputs = 1;
	constexpr size_t g_outputs = 1;
	constexpr bool constraint_inputrates = true;
	constexpr bool constraint_outputs = false;

	ArcsMat<n_states,n_states> A;
	ArcsMat<n_states,m_inputs> B;
	ArcsMat<g_outputs,n_states> C;
	int i = 3;
	double w_u = 0.1;
	double w_y = 0.1;
	double w_du = 0.05;
	auto identTest = 0.5*eye<3,3>();

	// disp(identTest);
	// disp(identTest^i);
	A.Set(0.5, 0,
		  0.2, 1);

	B.Set(0.3,
			0);

	C.Set(1, 0);

	ArcsMat<n_states,1> x0 = {0,
							  0};

	ArcsMat<m_inputs,1> u_z1 = {0};

	ArcsMat<g_outputs*p_hor,1> Y_REF = ones<g_outputs*p_hor,1>();

















	// ArcsMat<p_hor+1, 2> AbarTest;	
/*
	disp(AbarTest);
	setsubmatrix(AbarTest, C, 1, 1);
	setsubmatrix(AbarTest, C*A, 2, 1);
	setsubmatrix(AbarTest, C*(A^2), 3, 1);
	setsubmatrix(AbarTest, C*(A^3), 4, 1);
	setsubmatrix(AbarTest, C*(A^4), 5, 1);
	disp(AbarTest);
	*/


	LinearMPC<n_states,m_inputs,g_outputs,p_hor,c_hor,constraint_inputrates,constraint_outputs> mpcNoConstraintsTest(A,B,C, w_u, w_y, w_du, x0, u_z1, Y_REF);
	/*
	LinearMPC<2,1,5,3,true,0> mpcNoConstraintsTest(A,B);
	LinearMPC<2,1,5,3,true,2> mpcWithConstraintsTest(A,B);
	LinearMPC<2,1,5,3,true,0> mpcNoConstraintsTestAgain(A,B);
	*/
	// ここにオフライン計算のコードを記述	
	return EXIT_SUCCESS;	// 正常終了


		
}




