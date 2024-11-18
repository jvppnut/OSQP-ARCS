//! @file LinearMPC.cc
//! @brief Linear MPC implementation class
//!
//! Linear MPC based on OSQP solver using sparse formulation
//!
//! @date 2024/9/21
//! @author Juan Padron
//
// Copyright (C) 2011-2024 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.


#ifndef LINEARMPC
#define LINEARMPC

#include <type_traits>
#include <cstdio>
#include <cmath>
#include <cassert>
#include "ArcsMatrix.hh"
#include "OSQP_Solver.hh"

// ARCS組込み用マクロ
#ifdef ARCS_IN
	// ARCSに組み込まれる場合
	#include "ARCSassert.hh"
#else
	// ARCSに組み込まれない場合
	#define arcs_assert(a) (assert(a))
#endif

namespace ARCS{
using namespace ArcsMatrix;

template <size_t N_STATES, size_t M_INPUTS, size_t G_OUTPUTS, size_t P_HOR, size_t C_HOR, bool CONSTRAINT_INPUTRATES = false, bool CONSTRAINT_OUTPUTS = false>
    class LinearMPC{

		public:




		// //Constructor
		// LinearMPC(ArcsMat<N_STATES,N_STATES> A, ArcsMat<N_STATES,M_INPUTS> B, ArcsMat<G_OUTPUTS,N_STATES> C, double w_u, double w_y, double w_du):
		// P_mat(),
		// A_mat(),
		// q_vec(),
		// l_vec(),
		// u_vec(),
		// q1(),
		// q2(),
		// q3(),
		// mpcSolver()
		// {

		// 	printf("Testing: Entering constructor \n");


		// 	//BUILD P matrix -----------------

		// 	auto WU = w_u*eye<M_INPUTS*P_HOR,M_INPUTS*P_HOR>();
		// 	auto WY = w_y*eye<G_OUTPUTS*(P_HOR+1), G_OUTPUTS*(P_HOR+1)>();
		// 	auto WDU = w_du*eye<M_INPUTS*P_HOR,M_INPUTS*P_HOR>();

		// 	ArcsMat<G_OUTPUTS*(P_HOR+1),N_STATES> ABAR;
		// 	ArcsMat<G_OUTPUTS,N_STATES> tempMat;
		// 	for(size_t i = 0; i <= P_HOR; i++){
		// 		tempMat = C*(A^i);		
		// 		setsubmatrix(ABAR, tempMat, G_OUTPUTS*i+1,1);
		// 	}
		// 	disp(ABAR);

		// 	ArcsMat<G_OUTPUTS*(P_HOR+1),M_INPUTS*P_HOR> BBAR;
		// 	ArcsMat<G_OUTPUTS,M_INPUTS> tempMat2;
		// 	for(size_t i = 1; i <= P_HOR; i++){
		// 		for(size_t j = 0; j <= i-1; j++){
		// 			tempMat2 = C*(A^(i-1-j))*B;
		// 			setsubmatrix(BBAR, tempMat2, G_OUTPUTS*i + 1, M_INPUTS*j + 1);
		// 		}
		// 	}
		// 	disp(BBAR);

		// 	ArcsMat<P_HOR,C_HOR> Tbar;
		// 	setsubmatrix(Tbar, eye<C_HOR,C_HOR>(), 1, 1);
		// 	setsubmatrix(Tbar, ones<P_HOR-C_HOR,1>(), 1+C_HOR, C_HOR);
		// 	disp(Tbar);
		// 	ArcsMat<M_INPUTS*P_HOR,M_INPUTS*C_HOR> TBM;
		// 	Kron(Tbar,eye<M_INPUTS,M_INPUTS>(),TBM);
		// 	disp(TBM);

		// 	ArcsMat<M_INPUTS*P_HOR,M_INPUTS*P_HOR> D_DU;
		// 	auto d_DU = eye<P_HOR,P_HOR>();
		// 	for(size_t i=2; i<=P_HOR; i++){
		// 		d_DU(i,i-1) = -1.0;
		// 	}
		// 	Kron(d_DU, eye<M_INPUTS,M_INPUTS>(), D_DU);			
		// 	disp(d_DU);
		// 	disp(D_DU);

		// 	auto DBAR_DU = D_DU*TBM;
		// 	auto BBBAR = BBAR*TBM;
		// 	disp(DBAR_DU);
		// 	disp(BBBAR);

		// 	auto Q = tp(BBBAR)*WY*BBBAR + tp(DBAR_DU)*WDU*DBAR_DU + tp(TBM)*WU*TBM;
		// 	disp(Q);

		// 	setsubmatrix(P_mat, Q, 1, 1);
		// 	P_mat(M_INPUTS*C_HOR+1,M_INPUTS*C_HOR+1) = SLACK_EPS;
		// 	disp(P_mat);


		// 	//BUILD partial q vector -----------------

		// 	q1 = tp(ABAR)*WY*BBBAR;
		// 	q2 = WY*BBBAR;
		// 	q3 = WDU*DBAR_DU;

		// 	disp(q1);
		// 	disp(q2);
		// 	disp(q3);


		// 	//BUILD A matrix -----------------
		//     setsubmatrix(A_mat, eye<M_INPUTS*C_HOR,M_INPUTS*C_HOR>(),1,1);	//Input constraints
		// 	disp(A_mat);

		// 	// ArcsMat<M_INPUTS*C_HOR,M_INPUTS*C_HOR> D_DU_ineq;
		// 	//D_DU_ineq = getsubmatrix(D_DU,M_INPUTS*C_HOR, M_INPUTS*C_HOR);
		// 	//getsubmatrix(D_DU,D_DU_ineq,M_INPUTS*C_HOR, M_INPUTS*C_HOR);
		// 	// disp(D_DU);
		// 	// disp(D_DU_ineq);
		// 	// disp(( getsubmatrix<M_INPUTS*C_HOR,M_INPUTS*C_HOR>(D_DU,1,1) ));

		// 	if(CONSTRAINT_INPUTRATES){
		// 		auto D_DU_ineq = getsubmatrix<M_INPUTS*C_HOR,M_INPUTS*C_HOR>(D_DU,1,1);
		// 		setsubmatrix(A_mat, D_DU_ineq,M_INPUTS*C_HOR+1,1); //Input rate constraints
		// 		disp(D_DU_ineq);				
		// 		disp(A_mat);
		// 	}



		// 	A_mat(constraintsSize(),M_INPUTS*C_HOR+1) = 1;

		// 	disp(A_mat);
			

			
	

			
		// }



		//Constructor
		LinearMPC(ArcsMat<N_STATES,N_STATES> A, ArcsMat<N_STATES,M_INPUTS> B, ArcsMat<G_OUTPUTS,N_STATES> C, double w_u, double w_y, double w_du, ArcsMat<N_STATES,1> x0, ArcsMat<M_INPUTS,1> u_z1, ArcsMat<G_OUTPUTS*P_HOR,1> Y_REF):
		P_mat(),
		A_mat(),
		q_vec(),
		l_vec(),
		u_vec(),
		q_vec_r1(),
		q_vec_r2(),
		mpcSolver()
		{

			printf("Testing: Entering constructor \n");


			//------------ P matrix -----------------

			auto QU = w_u*eye<M_INPUTS*P_HOR,M_INPUTS*P_HOR>();
			auto QY_pre = w_y*eye<G_OUTPUTS*P_HOR, G_OUTPUTS*P_HOR>();
			auto QDU_pre = w_du*eye<M_INPUTS*P_HOR,M_INPUTS*P_HOR>();
			auto CBAR = Kron(eye<P_HOR,P_HOR>(),C);
			ArcsMat<M_INPUTS*P_HOR,M_INPUTS*P_HOR> D_DU;
			auto d_DU = eye<P_HOR,P_HOR>();
			for(size_t i=2; i<=P_HOR; i++){
				d_DU(i,i-1) = -1.0;
			}
			Kron(d_DU, eye<M_INPUTS,M_INPUTS>(), D_DU);		

			auto QY = tp(CBAR)*QY_pre*CBAR;
			auto QDU = tp(D_DU)*QDU_pre*D_DU;

			//Populate Hessian matrix P 
			setsubmatrix(P_mat, QY, 1, 1);
			setsubmatrix(P_mat, QU+QDU, 1+P_HOR*N_STATES, 1+P_HOR*N_STATES);
			P_mat(P_HOR*(N_STATES+M_INPUTS)+1,P_HOR*(N_STATES+M_INPUTS)+1) = SLACK_EPS;
			disp(P_mat);


			
			//------------ q vector -----------------
			
			ArcsMat<M_INPUTS*P_HOR,1> b0;
			setsubmatrix(b0, u_z1, 1, 1);
			q_vec_r1 = -2*tp(CBAR)*tp(QY_pre);
			q_vec_r2 = -2*tp(D_DU)*tp(QDU);	
			setsubmatrix(q_vec, q_vec_r1*Y_REF, 1, 1);
			setsubmatrix(q_vec, q_vec_r2*b0, 1, 1);

			disp(q_vec);


			//-----------A matrix ---------------
			ArcsMat<P_HOR,P_HOR> dyn_base = zeros<P_HOR,P_HOR>();
			for(size_t i=2; i<=P_HOR; i++){
				dyn_base(i,i-1) = -1.0;
			}
			auto Mx_DYN = Kron(eye<P_HOR,P_HOR>(),eye<N_STATES,N_STATES>()) + Kron(dyn_base,A);
			auto Mu_DYN = Kron(-eye<P_HOR,P_HOR>(),B);
			auto Meps_DYN = zeros<N_STATES*P_HOR,1>();

			ArcsMat<P_HOR*N_STATES,P_HOR*(N_STATES+M_INPUTS)+1> M_DYN;
			setsubmatrix(M_DYN, Mx_DYN, 1, 1);
			setsubmatrix(M_DYN, Mu_DYN, P_HOR*N_STATES+1, 1);
			setsubmatrix(M_DYN, Meps_DYN, P_HOR*(N_STATES+M_INPUTS)+1, 1);
			disp(M_DYN);

			
	

			
		}

		private:





		//Calculates number of constraints for constraint vectors and matrix computation
		static constexpr std::size_t constraintsSize() {		
			std::size_t num_of_constraints = M_INPUTS*C_HOR + 1;	

			if (CONSTRAINT_INPUTRATES==true) {
				num_of_constraints += M_INPUTS*C_HOR; // Add input rate constraints
			}
			if(CONSTRAINT_OUTPUTS==true){
				num_of_constraints += 2*G_OUTPUTS*(P_HOR+1);
			}

			return num_of_constraints;
    	}

		
		//Declaration of matrices for solver and solver itself
		const double SLACK_EPS = 1e5;


		ArcsMat<P_HOR*(N_STATES+M_INPUTS)+1,P_HOR*(N_STATES+M_INPUTS)+1> P_mat;
		ArcsMat<constraintsSize(),M_INPUTS*C_HOR+1> A_mat;
		ArcsMat<P_HOR*(N_STATES+M_INPUTS)+1,1> q_vec;
		ArcsMat<constraintsSize(),1> l_vec;
		ArcsMat<constraintsSize(),1> u_vec;		
		ArcsMat<P_HOR*N_STATES, P_HOR*G_OUTPUTS> q_vec_r1;
		ArcsMat<P_HOR*M_INPUTS, P_HOR*M_INPUTS> q_vec_r2;


		OSQP_Solver<M_INPUTS*C_HOR+1,constraintsSize()> mpcSolver;
    };


}

#endif
