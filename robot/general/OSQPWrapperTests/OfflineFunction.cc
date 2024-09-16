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

	OSQPInt exitflag;
	ArcsMat<2,1> solVector;	

	/* -------- SOLVE TEST ------------*/
	ArcsMat<2,2> P = {			// 宣言と同時に値をセットする場合
		4, 1,
		1, 2
	};


	ArcsMat<2,1> q = {
		1,
		1
	};

	ArcsMat<3,2> A = {			// 宣言と同時に値をセットする場合
		1, 1,
		1, 0,
		0, 1
	};


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


	OSQP_Solver<2,3> testSolver(P,A,q,l,u);
	std::array<OSQPFloat,2> solution;
	exitflag = testSolver.solve(solution);
	solVector.LoadArray(solution);
	printf("Exitflag: %lld \n",exitflag);
	disp(solVector);


	/* -------- UPDATE Q VECTOR TEST ------------*/
	ArcsMat<2,1> qnew = {
		0.3,
		0.9
	};

	exitflag = testSolver.Update_qVec(qnew);
	exitflag = testSolver.solve(solution);
	solVector.LoadArray(solution);
	printf("Exitflag: %lld \n",exitflag);
	disp(solVector);



	/* -------- UPDATE QLU VECTORS TEST ------------*/

	ArcsMat<2,1> qnew2 = {
		0.5,
		0.5
	};

	ArcsMat<3,1> lnew = {
		-OSQP_INFTY,
		-OSQP_INFTY,
		-OSQP_INFTY
	};

	ArcsMat<3,1> unew = {
		OSQP_INFTY,
		OSQP_INFTY,
		OSQP_INFTY
	};

	exitflag = testSolver.Update_Vecs(qnew2,lnew,unew);
	exitflag = testSolver.solve(solution);
	solVector.LoadArray(solution);
	printf("Exitflag: %lld \n",exitflag);
	disp(solVector);

	/* -------- UPDATE P MATRIX TEST ------------*/

		ArcsMat<2,2> Pnew = {			
		.3, 0.005,
		0.005, .2
	};

	exitflag = testSolver.Update_Vecs(q,l,u);
	exitflag = testSolver.Update_PMatrix(Pnew);
	exitflag = testSolver.solve(solution);
	solVector.LoadArray(solution);
	printf("Exitflag: %lld \n",exitflag);
	disp(solVector);

	/* -------- UPDATE A MATRIX TEST ------------*/

	ArcsMat<3,2> Anew = {			
		1.2, 1.1,
		1.5, 0,
		0, 0.8
	};

	exitflag = testSolver.Update_Vecs(q,l,u);
	exitflag = testSolver.Update_PMatrix(P);
	exitflag = testSolver.Update_AMatrix(Anew);
	exitflag = testSolver.solve(solution);
	solVector.LoadArray(solution);
	printf("Exitflag: %lld \n",exitflag);
	disp(solVector);


	// ここにオフライン計算のコードを記述
	
	return EXIT_SUCCESS;	// 正常終了
}

