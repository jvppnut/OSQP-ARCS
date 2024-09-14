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
#include "OSQP_Solver.hh"

using namespace ARCS;

//! @brief エントリポイント
//! @return 終了ステータス
int main(void){
	printf("ARCS OFFLINE CALCULATION MODE\n");
	printf("OSQP Test File\n");


	//OSQP_Solver<2,2> testSolver;


	ArcsMat<3,3> Perror = {
		1, 2, 3,
		4, 5, 6,
		7, 8, 9
	};




	ArcsMat<2,2> P = {			// 宣言と同時に値をセットする場合
		4, 1,
		1, 2
	};


	//double a = P(1,1);
	//printf("P(1,1): %f", a);

	ArcsMat<2,1> q = {
		1,
		1
	};

	ArcsMat<3,2> A = {			// 宣言と同時に値をセットする場合
		1, 1,
		1, 0,
		0, 1
	};

//Doesn't return OSQP error of index bound

	ArcsMat<3,1> l = {
		1,
		0,
		0
	};

	ArcsMat<3,1> u = {
		1,
		0.7,
		0.7
	};

	

//Returns OSQP error of index bound
/*
	ArcsMat<3,1> l = {
		1,
		0.7,
		0.7
	};

	ArcsMat<3,1> u = {
		1,
		0,
		0
	};
*/

	//OSQP_Solver<2,3> testSolver(Perror,A,q,l,u); //Returns size error
	OSQP_Solver<2,3> testSolver(P,A,q,l,u);
	std::array<OSQPFloat,2> solution;
	OSQPInt exitflag;
	exitflag = testSolver.solve(solution);
	ArcsMat<2,1> solVector;	
	solVector.LoadArray(solution);
	printf("Exitflag: %lld \n",exitflag);
	disp(solVector);



	// ここにオフライン計算のコードを記述
	
	return EXIT_SUCCESS;	// 正常終了
}

