//! @file LinearMPC.cc
//! @brief Linear MPC implementation class
//!
//! Linear MPC based on OSQP solver
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


template <size_t N_STATES, size_t M_INPUTS, size_t P_HOR, size_t C_HOR, bool CONSTRAINT_INPUTRATES = false, size_t CONSTRAINED_STATES = 0>
    class LinearMPC{

		public:


/*
		template<int M = CONSTRAINED_STATES, typename std::enable_if_t<(M == 0), int> = 0>
		LinearMPC(int a):
		P_solver(),
		A_solver(),
		q_solver(),
		l_solver(),
		u_solver(),
		mpcSolver()
		{
			printf("Testing: Constructor without state constraints \n");
		}

		template<int M = CONSTRAINED_STATES, typename std::enable_if_t<(M > 0), int> = 0>
		LinearMPC(int b, int c):
		P_solver(),
		A_solver(),
		q_solver(),
		l_solver(),
		u_solver(),
		mpcSolver()
		{
			printf("Testing: Constructor with state constraints \n");
		}
*/		


		//Constructor for case when we don't have neither input rate or state constraints
		template<size_t M = CONSTRAINED_STATES, bool L = CONSTRAINT_INPUTRATES, std::enable_if_t<(M == 0), size_t> = 0, std::enable_if_t<(L == false), size_t> = 0>
		LinearMPC(ArcsMat<N_STATES,N_STATES> A, ArcsMat<N_STATES,M_INPUTS> B):
		P_solver(),
		A_solver(),
		q_solver(),
		l_solver(),
		u_solver(),
		mpcSolver()
		{
			printf("Testing: Constructor without state constraints \n");
		}


		//Constructor for case when we have input rate constraints but no state constraints
		template<size_t M = CONSTRAINED_STATES, bool L = CONSTRAINT_INPUTRATES, std::enable_if_t<(M == 0), size_t> = 0, std::enable_if_t<(L == true), size_t> = 0>
		LinearMPC(ArcsMat<N_STATES,N_STATES> A, ArcsMat<N_STATES,M_INPUTS> B):
		P_solver(),
		A_solver(),
		q_solver(),
		l_solver(),
		u_solver(),
		mpcSolver()
		{
			printf("Testing: Constructor with state constraints \n");
		}
		

		//Constructor for case when we have state constraints but no input rate constraints
		template<size_t M = CONSTRAINED_STATES, bool L = CONSTRAINT_INPUTRATES, std::enable_if_t<(M > 0), size_t> = 0, std::enable_if_t<(L == false), size_t> = 0>
		LinearMPC(ArcsMat<N_STATES,N_STATES> A, ArcsMat<N_STATES,M_INPUTS> B):
		P_solver(),
		A_solver(),
		q_solver(),
		l_solver(),
		u_solver(),
		mpcSolver()
		{
			printf("Testing: Constructor with state constraints \n");
		}


		//Constructor for case when we have state constraints and input rate constraints
		template<size_t M = CONSTRAINED_STATES, bool L = CONSTRAINT_INPUTRATES, std::enable_if_t<(M > 0), size_t> = 0, std::enable_if_t<(L == true), size_t> = 0>
		LinearMPC(ArcsMat<N_STATES,N_STATES> A, ArcsMat<N_STATES,M_INPUTS> B):
		P_solver(),
		A_solver(),
		q_solver(),
		l_solver(),
		u_solver(),
		mpcSolver()
		{
			printf("Testing: Constructor with state constraints \n");
		}
		


		private:


		//Generates matrix P and vector q for cost function
		void generateCostFunction(const ArcsMat<N_STATES,N_STATES>& A, const ArcsMat<N_STATES,M_INPUTS>& B, double w_y, double w_u, double w_du, double w_eps)
		{
			ArcsMat<P_HOR,P_HOR> WY = eye<P_HOR,P_HOR>()*wy;
		}

		//Generates constraint matrix A and constraint lower bound vector "l" and upper bound vector "u"
		void generateConstraintMatrixAndVectors(const ArcsMat<N_STATES,N_STATES>& A, const ArcsMat<N_STATES,M_INPUTS>& B, const ArcsMat<M_INPUTS,1>& U_MAX, const ArcsMat<M_INPUTS,1>& U_MIN, const ArcsMat<M_INPUTS,1>& DU_MAX, const ArcsMat<M_INPUTS,1>& DU_MIN, const ArcsMat<CONSTRAINED_STATES,1>& X_MAX, const ArcsMat<CONSTRAINED_STATES,1>& X_MIN, const ArcsMat<CONSTRAINED_STATES,N_STATES>& T)
		{

		}


	

		static constexpr std::size_t constraintsSize() {		
			std::size_t num_of_constraints = M_INPUTS*C_HOR + CONSTRAINED_STATES*P_HOR + 1;	

			if (CONSTRAINT_INPUTRATES==true) {
				num_of_constraints += M_INPUTS*C_HOR; // Add input rate constraints
			}

			return num_of_constraints;
    	}

		
		//Declaration of matrices for solver and solver itself
		ArcsMat<M_INPUTS*C_HOR+1,M_INPUTS*C_HOR+1> P_solver;
		ArcsMat<constraintsSize(),M_INPUTS*C_HOR+1> A_solver;
		ArcsMat<M_INPUTS*C_HOR+1,1> q_solver;
		ArcsMat<constraintsSize(),1> l_solver;
		ArcsMat<constraintsSize(),1> u_solver;		
		OSQP_Solver<M_INPUTS*C_HOR+1,constraintsSize()> mpcSolver;





    };


}

#endif
