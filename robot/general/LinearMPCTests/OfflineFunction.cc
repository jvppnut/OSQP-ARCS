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
#include "ArcsControl.hh"
#include "LinearMPC.hh"
#include "OSQP_Solver.hh"


using namespace ARCS;
using namespace ArcsMatrix;
using namespace ArcsControl;

//! @brief エントリポイント
//! @return 終了ステータス
int main(void){
	printf("ARCS OFFLINE CALCULATION MODE\n");
	printf("Linear MPC tests \n");







	constexpr size_t p_hor = 5;
	constexpr size_t c_hor = 3;
	constexpr size_t n_states = 2;
	// constexpr size_t m_inputs = 1;
	constexpr size_t m_inputs = 2;
	constexpr size_t g_outputs = 1;
	constexpr bool constraint_inputrates = false;
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

	// B.Set(0.3,
	// 		0);

		B.Set(0.3, 0.4,
			0.7, 0.0);

	C.Set(1, 0);


	ArcsMat<n_states,1> x0 = {1,
							  3};

	// ArcsMat<m_inputs,1> u_z1 = {0};

	ArcsMat<m_inputs,1> u_z1 = {2,
								3};


	// ArcsMat<m_inputs,1> u_max = {10};
	 ArcsMat<m_inputs,1> u_max = {10,
	 							  20};
	ArcsMat<m_inputs,1> u_min = -u_max;

		 ArcsMat<m_inputs,1> du_max = {500,
	 							  1000};
	ArcsMat<m_inputs,1> du_min = -du_max;

	ArcsMat<g_outputs,1> y_max = {333};
	ArcsMat<g_outputs,1> y_min = -y_max;

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





	// LinearMPC<n_states,m_inputs,g_outputs,p_hor,c_hor,false,false> mpcNoConstraintsTest(A,B,C, w_u, w_y, w_du, x0, u_z1, Y_REF, u_min, u_max);
	// LinearMPC<n_states,m_inputs,g_outputs,p_hor,c_hor,true, false> mpcRatesTest(A,B,C, w_u, w_y, w_du, x0, u_z1, Y_REF, u_min, u_max, du_min, du_max);
	// LinearMPC<n_states,m_inputs,g_outputs,p_hor,c_hor,false, true> mpcOutputsTest(A,B,C, w_u, w_y, w_du, x0, u_z1, Y_REF, u_min, u_max, y_min, y_max);
	LinearMPC<n_states,m_inputs,g_outputs,p_hor,c_hor,true, true> mpcAllConstraintsTest(A,B,C, w_u, w_y, w_du, x0, u_z1, Y_REF, u_min, u_max, du_min, du_max, y_min, y_max);
	/*
	LinearMPC<2,1,5,3,true,0> mpcNoConstraintsTest(A,B);
	LinearMPC<2,1,5,3,true,2> mpcWithConstraintsTest(A,B);
	LinearMPC<2,1,5,3,true,0> mpcNoConstraintsTestAgain(A,B);
	*/
	// ここにオフライン計算のコードを記述	
	// mpcNoConstraintsTest.checkMatrices();
	// mpcRatesTest.checkMatrices();
	// mpcOutputsTest.checkMatrices();
	// mpcAllConstraintsTest.checkMatrices();

	// mpcAllConstraintsTest.update(2.0*Y_REF, 0.5*x0, 0.5*u_z1);
	// mpcAllConstraintsTest.update(Y_REF, x0, u_z1);
	// mpcAllConstraintsTest.checkMatrices();


	// MatExport MatFile1("save_test.mat");// MATファイルを新規作成
	// MatFile1.Save("Aexp1", A);		// MATファイルに行列データを書き出し




	// --------------- EXAMPLE 1 - 1 mass (inertia) system test -------


	///////  (Verification of Hessian, gradient vector, constraint matrix and vectors) /////
	constexpr double Ts_ex1 = 0.01;
	constexpr double Tsim_ex1 = 20;
	constexpr int imax_ex1 = Tsim_ex1/Ts_ex1; 
	constexpr double Jm = 0.1; //Motor inertia 
	constexpr double Kt = 1.0;  //Torque constant
	constexpr double Dm = 1e-2;
	constexpr double u1c_ex1 = 5.0;  //Input constraint i.e. current limit
	constexpr double du1c_ex1 = 100.0; //Input rate of change constraint
	constexpr double y1c_ex1 = 0.201*M_PI; //State x2 constraint: Angular position.
	constexpr size_t n_ex1= 2; //Dimension of state vector: We only have x1: ang. velocity and x2: ang. position
	constexpr size_t m_ex1 = 1; //Dimension of input vector: we only have input current
	constexpr size_t g_ex1 = 1; //Dimension of output vector
	constexpr double wY_ex1 = 100.0;
	constexpr double wU_ex1 = 0.0;
	constexpr double wDU_ex1 = 1.0;
	constexpr size_t phor_ex1 = 60;
	constexpr size_t chor_ex1 = 4;

	ArcsMat<n_ex1,n_ex1> Ad_ex1;
	ArcsMat<n_ex1,m_ex1> Bd_ex1;
	ArcsMat<n_ex1,n_ex1> A_ex1 = {-Dm/Jm, 0,
								1,	0};
	ArcsMat<n_ex1,m_ex1> B_ex1 = {Kt/Jm, 
									0};
	ArcsMat<g_ex1,n_ex1> C_ex1 = {0, 1};
	Discretize(A_ex1,B_ex1,Ad_ex1,Bd_ex1,Ts_ex1);

	ArcsMat<n_ex1,1> x0_ex1 = {0,
								0};
	ArcsMat<m_ex1,1> uz1_ex1 = {0};
	ArcsMat<g_ex1*phor_ex1,1> YREF_ex1(0);
	// disp(A_ex1);
	// disp(Ad_ex1);
	// disp(Bd_ex1);
	// disp(B_ex1);



	ArcsMat<m_ex1,1> umax_ex1 = {u1c_ex1};
	ArcsMat<m_ex1,1> umin_ex1 = -umax_ex1;
	ArcsMat<m_ex1,1> dumax_ex1 = {du1c_ex1};
	ArcsMat<m_ex1,1> dumin_ex1 = -dumax_ex1;
	ArcsMat<g_ex1,1> ymax_ex1 = {y1c_ex1};
	ArcsMat<g_ex1,1> ymin_ex1 = -ymax_ex1;

	LinearMPC<n_ex1,m_ex1,g_ex1,phor_ex1,chor_ex1,true, true> mpc_ex1(Ad_ex1,Bd_ex1,C_ex1,
	 wU_ex1, wY_ex1, wDU_ex1, x0_ex1, uz1_ex1, YREF_ex1, 
	 umin_ex1, umax_ex1, dumin_ex1, dumax_ex1, ymin_ex1, ymax_ex1);

	mpc_ex1.testOutputToMAT("ex1.mat");

	//Test now with non-equal-to-zero x0, uz1 and y_ref
	x0_ex1.Set(1.2, 
				0.11);
	uz1_ex1.Set(1.6);
	YREF_ex1.FillAll(0.15*M_PI);

	mpc_ex1.update(YREF_ex1, x0_ex1, uz1_ex1);
	mpc_ex1.testOutputToMAT("ex2.mat");

	// --------------- EXAMPLE 2 - 1 mass (inertia) control test -------
	constexpr double Ts_ex2 = 0.01; //1 ms
	constexpr double Tsim_ex2 = 20; //20 secs
	constexpr int imax_ex2 = Tsim_ex2/Ts_ex2; 
	constexpr double u1c_ex2 = 5.0;  //Input constraint i.e. current limit
	constexpr double du1c_ex2 = 10000.0; //Input rate of change constraint
	// constexpr double y1c_ex2 = 0.201*M_PI; //State x2 constraint: Angular position.
	constexpr double y1c_ex2 = 1.05*M_PI; //State x2 constraint: Angular position.
	//NOTE: Solver loses feasibility for y1c_ex2 = 1.05*M_PI, yref = M_PI, P_HOR = 20, C_HOR = 4
	//Gotta check why this occurs (slack variable is supposed to shield us from this!)
	constexpr size_t n_ex2= 2; //Dimension of state vector: We only have x1: ang. velocity and x2: ang. position
	constexpr size_t m_ex2 = 1; //Dimension of input vector: we only have input current
	constexpr size_t g_ex2 = 1; //Dimension of output vector
	constexpr double wY_ex2 = 100.0;
	constexpr double wU_ex2 = 0.0;
	constexpr double wDU_ex2 = 1.0;
	constexpr size_t phor_ex2 = 60;
	constexpr size_t chor_ex2 = 4;




	//---Variables
	//For storage
	ArcsMat<imax_ex2,1> t_ex2;
	ArcsMat<imax_ex2,1> iq_ex2;
	ArcsMat<imax_ex2,1> theta_ref_ex2;
	ArcsMat<imax_ex2,1> theta_m_ex2;
	ArcsMat<imax_ex2,1> omega_m_ex2;
	ArcsMat<imax_ex2,1> slackvar_ex2;
	ArcsMat<imax_ex2,1> solruntime_ex2;

	//For control
	double tsim_ex2 = 0.0;
	double theta_ref = 0.0;
	ArcsMat<n_ex2,1> x0_ex2 = {0,
								0};
	ArcsMat<m_ex2,1> uz1_ex2 = {0};
	ArcsMat<g_ex2*phor_ex2,1> YREF_ex2(0);
	ArcsMat<m_ex2,1> umax_ex2 = {u1c_ex2};
	ArcsMat<m_ex2,1> umin_ex2 = -umax_ex2;
	ArcsMat<m_ex2,1> dumax_ex2 = {du1c_ex2};
	ArcsMat<m_ex2,1> dumin_ex2 = -dumax_ex2;
	ArcsMat<g_ex2,1> ymax_ex2 = {y1c_ex2};
	ArcsMat<g_ex2,1> ymin_ex2 = -ymax_ex2;		
	ArcsMat<n_ex2,1> xnext_ex2;
	ArcsMat<m_ex2,1> u_ex2;
	ArcsMat<phor_ex2*m_ex2,1> U_opt;
	ArcsMat<phor_ex2*n_ex2,1> X_pred;
	double slackvar = 0.0;
	OSQP_Status solver_status; 
	//With input rate and output constraints
	LinearMPC<n_ex2,m_ex2,g_ex2,phor_ex2,chor_ex2,true, true> mpc_ex2(Ad_ex1,Bd_ex1,C_ex1,
	 wU_ex2, wY_ex2, wDU_ex2, x0_ex2, uz1_ex2, YREF_ex2, 
	 umin_ex2, umax_ex2, dumin_ex2, dumax_ex2, ymin_ex2, ymax_ex2);

	//Without output constraints
	// LinearMPC<n_ex2,m_ex2,g_ex2,phor_ex2,chor_ex2,true, false> mpc_ex2(Ad_ex1,Bd_ex1,C_ex1,
	//  wU_ex2, wY_ex2, wDU_ex2, x0_ex2, uz1_ex2, YREF_ex2, 
	//  umin_ex2, umax_ex2, dumin_ex2, dumax_ex2);

	//Without input rate and output constraints
	// LinearMPC<n_ex2,m_ex2,g_ex2,phor_ex2,chor_ex2,false, false> mpc_ex2(Ad_ex1,Bd_ex1,C_ex1,
	//  wU_ex2, wY_ex2, wDU_ex2, x0_ex2, uz1_ex2, YREF_ex2, 
	//  umin_ex2, umax_ex2);






	////// Control loop //////////////
	for(size_t i = 1; i<=imax_ex2; i++)
	{
		if(tsim_ex2 >= 2.0)
		{
			// theta_ref = 0.2*M_PI;
			theta_ref = M_PI;
		}

		YREF_ex2.FillAll(theta_ref);
		// disp(YREF_ex2);
		mpc_ex2.update(YREF_ex2, x0_ex2, uz1_ex2);
		mpc_ex2.solve(U_opt, X_pred, slackvar, solver_status);
	


		u_ex2 = getsubmatrix<m_ex2,1>(U_opt,1,1);	//Extract optimal solution for present time

		printf("t: %f [s], u: %f [A], theta_ref: %f [rad], theta_m: %f [rad], OSQP run time: %f [us] \n",tsim_ex2, u_ex2(1,1), theta_ref, x0_ex2(2,1), 1e6*solver_status.run_time);

		xnext_ex2 = Ad_ex1*x0_ex2 + Bd_ex1*u_ex2;	//Advance simulation by one time step

		x0_ex2 = xnext_ex2;	//Pass past values
		uz1_ex2 = u_ex2;


		//Store data
		t_ex2(i,1) = tsim_ex2;
		iq_ex2(i,1) = u_ex2(1,1);		
		omega_m_ex2(i,1) = xnext_ex2(1,1);
		theta_m_ex2(i,1) = xnext_ex2(2,1);
		theta_ref_ex2(i,1) = theta_ref;
		slackvar_ex2(i,1) = slackvar;
		solruntime_ex2(i,1) = solver_status.run_time;
		
		tsim_ex2 += Ts_ex2;	//Advance time

	}



	MatExport toMATLAB("mpcControlEx2.mat");
	toMATLAB.Save("t_vec_exp", t_ex2);
	toMATLAB.Save("u_vec_exp", iq_ex2);
	toMATLAB.Save("thetam_vec_exp", theta_m_ex2);
	toMATLAB.Save("thetaref_vec_exp", theta_ref_ex2);
	toMATLAB.Save("omegam_vec_exp", omega_m_ex2);
	toMATLAB.Save("slack_vec_exp", slackvar_ex2);
	toMATLAB.Save("solruntime_vec_exp", solruntime_ex2);

	//TODO....
	// --------------- EXAMPLE 3 - MIMO system test -------

	 


	
	

	return EXIT_SUCCESS;	// 正常終了


		
}




