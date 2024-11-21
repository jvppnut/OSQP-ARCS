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

		//Main constructor, initializes all parts of the optimization problem that are independent of 
		//the constraint selection flags (independent of CONSTRAINT_INPUTRATES and CONSTRAINT_OUTPUTS)
		LinearMPC(const ArcsMat<N_STATES,N_STATES>& A, const ArcsMat<N_STATES,M_INPUTS>& B, const ArcsMat<G_OUTPUTS,N_STATES>& C,
		 double w_u, double w_y, double w_du,
		 const ArcsMat<N_STATES,1>& x0, const ArcsMat<M_INPUTS,1>& u_z1, const ArcsMat<G_OUTPUTS*P_HOR,1>& Y_REF,
		 const ArcsMat<M_INPUTS,1>& u_min, const ArcsMat<M_INPUTS,1>& u_max):
		P_mat(),
		A_mat(),
		q_vec(),
		l_vec(),
		u_vec(),
		q_vec_r1(),
		q_vec_r2(),
		A_stored(),
		du_min_stored(),
		du_max_stored(),
		qpSolver()
		{

			printf("Testing: Entering constructor \n");




			arcs_assert(w_u > 0);
			arcs_assert(w_du > 0);
			arcs_assert(w_y > 0);


			//Some necessary storage
			A_stored = A;


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
			// disp(P_mat);


			
			//------------ q vector -----------------
			
			//TODO: Fix submatrix sizes
			ArcsMat<M_INPUTS*P_HOR,1> b0;
			setsubmatrix(b0, u_z1, 1, 1);
			q_vec_r1 = -2*tp(CBAR)*tp(QY_pre);
			q_vec_r2 = -2*tp(D_DU)*tp(QDU);	
			setsubmatrix(q_vec, q_vec_r1*Y_REF, 1, 1);
			setsubmatrix(q_vec, q_vec_r2*b0, 1+P_HOR*N_STATES, 1);

			// disp(q_vec);


			//-----------A matrix ---------------

			//--- Dynamic constraints
			ArcsMat<P_HOR*N_STATES,P_HOR*(N_STATES+M_INPUTS)+1> M_DYN;
			ArcsMat<P_HOR*N_STATES,1> l_DYN;
			ArcsMat<P_HOR,P_HOR> dyn_base = zeros<P_HOR,P_HOR>();
			for(size_t i=2; i<=P_HOR; i++){
				dyn_base(i,i-1) = -1.0;
			}
			auto Mx_DYN = Kron(eye<P_HOR,P_HOR>(),eye<N_STATES,N_STATES>()) + Kron(dyn_base,A);
			auto Mu_DYN = Kron(-eye<P_HOR,P_HOR>(),B);
			
			setsubmatrix(M_DYN, Mx_DYN, 1, 1);
			setsubmatrix(M_DYN, Mu_DYN, 1, P_HOR*N_STATES+1);			
			setsubmatrix(l_DYN, A*x0, 1, 1);
		

			//--- Input constraints
			ArcsMat<P_HOR*M_INPUTS,P_HOR*(N_STATES+M_INPUTS)+1> M_U;	
			auto Mu_U = eye<P_HOR*M_INPUTS,P_HOR*M_INPUTS>();

			setsubmatrix(M_U, Mu_U, 1, P_HOR*N_STATES+1);
			auto l_U = Kron(ones<P_HOR,1>(),u_min);
			auto u_U = Kron(ones<P_HOR,1>(),u_max);


			//-- Control horizon constraint (in the case that input rate constraints are not enabled)

			/*TODO: Blocking movements should be implemented by direct substitution (condensed formulation)
				 in order to reduce the number of variables to optimize. I was too lazy at the moment so just
				  implemented it as additional constraints
			*/
			if constexpr(P_HOR-C_HOR>0)
			{
				if(CONSTRAINT_INPUTRATES == false)
				{

					ArcsMat<(P_HOR-C_HOR)*M_INPUTS,P_HOR*(N_STATES+M_INPUTS)+1> M_CHOR;	
					ArcsMat<(P_HOR-C_HOR)*M_INPUTS,1> l_CHOR;

					ArcsMat<P_HOR-C_HOR,P_HOR+1-C_HOR> block_base;
					for(size_t i=1; i<=P_HOR-C_HOR; i++)
					{
					block_base(i,i) = -1.0;
					block_base(i,i+1) = 1.0;
					}

					// disp(block_base);
					auto block_mat = Kron(block_base,eye<M_INPUTS,M_INPUTS>());
					setsubmatrix(M_CHOR, block_mat, 1, 1 + P_HOR*N_STATES + (C_HOR-1)*M_INPUTS);
					// disp(block_mat);
					// disp(M_CHOR);

					//Stack on constraint matrix A_mat
					setsubmatrix(A_mat, M_CHOR, 1+P_HOR*(N_STATES+M_INPUTS), 1);
					// disp(A_mat);
				}


			}

			//-----------Stack everything together ---------------

			//Populate constraint matrix A_mat
			setsubmatrix(A_mat, M_DYN, 1, 1);
			setsubmatrix(A_mat, M_U, 1+P_HOR*N_STATES, 1);
			A_mat(constraintsSize(),1+P_HOR*(M_INPUTS+N_STATES)) = 1;

			//Populate lower and upper constraints vectors l_vec and u_vec
			setsubmatrix(u_vec, l_DYN, 1, 1); 
			setsubmatrix(u_vec, u_U, 1+P_HOR*N_STATES, 1);
			u_vec(constraintsSize(),1) = OSQP_INFTY; //Slack variable upper bound (infinity)

			setsubmatrix(l_vec, l_DYN, 1, 1);
			setsubmatrix(l_vec, l_U, 1+P_HOR*N_STATES, 1);
			l_vec(constraintsSize(),1) = 0;		//Slack variable lower bound (zero, must be positive)

			// disp(A_mat);
			// disp(u_vec);
			// disp(l_vec);

			//Initialize solver rightaway if we don't have input rate and output constraints
			//Else, we have to wait until A_mat and u_vec, l_vec are properly populated with
			//their corresponding constraints
			if constexpr(!CONSTRAINT_INPUTRATES && !CONSTRAINT_OUTPUTS)
			{
				initializeSolver();
			}
			
		}


		//Constructor for input rate constraints
		//Using SFINAE here to enable constructor only and only when input rate constraints
		// are enabled through template argument CONSTRAINT_INPUTRATES
		template<bool CRATES = CONSTRAINT_INPUTRATES, bool COUTPUTS = CONSTRAINT_OUTPUTS, std::enable_if_t<(CRATES && !COUTPUTS), size_t> = 0>
		// template<std::enable_if_t<(CONSTRAINT_INPUTRATES), size_t> = 0>
		// template <typename = std::enable_if_t<CONSTRAINT_INPUTRATES && !CONSTRAINT_OUTPUTS>>
		LinearMPC(const ArcsMat<N_STATES,N_STATES>& A, const ArcsMat<N_STATES,M_INPUTS>& B, const ArcsMat<G_OUTPUTS,N_STATES>& C,
		 double w_u, double w_y, double w_du,
		 const ArcsMat<N_STATES,1>& x0, const ArcsMat<M_INPUTS,1>& u_z1, const ArcsMat<G_OUTPUTS*P_HOR,1>& Y_REF,
		 const ArcsMat<M_INPUTS,1>& u_min, const ArcsMat<M_INPUTS,1>& u_max,
		 const ArcsMat<M_INPUTS,1>& du_min, const ArcsMat<M_INPUTS,1>& du_max):
		 LinearMPC(A, B, C, w_u, w_y, w_du, x0, u_z1, Y_REF, u_min, u_max)
		{
			printf("Testing: Entering constructor for when input rates constraints are activated \n");

			//Call basic constructor
			// (This doesn't work as it only creates a local copy that is deleted)
			// Basic constructor needs to be called through constructor delegation as done above
			// LinearMPC(A, B, C, w_u, w_y, w_du, x0, u_z1, Y_REF, u_min, u_max);


			
			//Some necessary storage
			du_min_stored = du_min;
			du_max_stored = du_max;

			//Stack input rate constraints
			ArcsMat<M_INPUTS*P_HOR,M_INPUTS*P_HOR> D_DU;
			ArcsMat<M_INPUTS*P_HOR,1> l_DU;
			ArcsMat<M_INPUTS*P_HOR,1> u_DU;
			ArcsMat<M_INPUTS*P_HOR,1> b0;
			setsubmatrix(b0, u_z1, 1, 1);
			auto d_DU = eye<P_HOR,P_HOR>();
			for(size_t i=2; i<=P_HOR; i++){
				d_DU(i,i-1) = -1.0;
			}
			Kron(d_DU, eye<M_INPUTS,M_INPUTS>(), D_DU);		
			setsubmatrix(A_mat, D_DU, 1+P_HOR*N_STATES + P_HOR*M_INPUTS, P_HOR*N_STATES+1);
			setsubmatrix(l_DU,	Kron(ones<C_HOR,1>(),du_min)	,1,1);
			setsubmatrix(u_DU,	Kron(ones<C_HOR,1>(),du_max)	,1,1);	

			setsubmatrix(l_vec,	l_DU + b0,1+P_HOR*N_STATES + P_HOR*M_INPUTS,1);
			setsubmatrix(u_vec,	u_DU + b0,1+P_HOR*N_STATES + P_HOR*M_INPUTS,1);

			//Initialize solver
			initializeSolver();

		}


		//Constructor for output constraints
		//Using SFINAE here to enable constructor only and only when output constraints
		// are enabled through template argument CONSTRAINT_INPUTRATES
		template<bool CRATES = CONSTRAINT_INPUTRATES, bool COUTPUTS = CONSTRAINT_OUTPUTS, std::enable_if_t<(!CRATES && COUTPUTS), size_t> = 0>
		// template<std::enable_if_t<(CONSTRAINT_INPUTRATES), size_t> = 0>
		// template <typename = std::enable_if_t<CONSTRAINT_INPUTRATES && !CONSTRAINT_OUTPUTS>>
		LinearMPC(const ArcsMat<N_STATES,N_STATES>& A, const ArcsMat<N_STATES,M_INPUTS>& B, const ArcsMat<G_OUTPUTS,N_STATES>& C,
		 double w_u, double w_y, double w_du,
		 const ArcsMat<N_STATES,1>& x0, const ArcsMat<M_INPUTS,1>& u_z1, const ArcsMat<G_OUTPUTS*P_HOR,1>& Y_REF,
		 const ArcsMat<M_INPUTS,1>& u_min, const ArcsMat<M_INPUTS,1>& u_max,
		 const ArcsMat<G_OUTPUTS,1>& y_min, const ArcsMat<G_OUTPUTS,1>& y_max):
		 LinearMPC(A, B, C, w_u, w_y, w_du, x0, u_z1, Y_REF, u_min, u_max)
		{
			printf("Testing: Entering constructor for when output constraints are activated \n");



			//Call basic constructor
			// (This doesn't work as it only creates a local copy that is deleted)
			// Basic constructor needs to be called through constructor delegation as done above
			// LinearMPC(A, B, C, w_u, w_y, w_du, x0, u_z1, Y_REF, u_min, u_max);

			//Stack output constraints
			ArcsMat<P_HOR*G_OUTPUTS,P_HOR*N_STATES> Mx_Y = Kron(eye<P_HOR,P_HOR>(),C);
			ArcsMat<P_HOR*G_OUTPUTS,1> u_Y_top = Kron(ones<P_HOR,1>(),y_max);		
			ArcsMat<P_HOR*G_OUTPUTS,1> l_Y_top(-OSQP_INFTY); 	
			ArcsMat<P_HOR*G_OUTPUTS,1> u_Y_bottom(OSQP_INFTY); 	
			ArcsMat<P_HOR*G_OUTPUTS,1> l_Y_bottom = Kron(ones<P_HOR,1>(),y_min);

			setsubmatrix(A_mat, Mx_Y, 1+P_HOR*N_STATES + P_HOR*M_INPUTS + (P_HOR-C_HOR)*M_INPUTS, 1);
			setsubmatrix(A_mat, Mx_Y, 1+P_HOR*N_STATES + P_HOR*M_INPUTS + (P_HOR-C_HOR)*M_INPUTS + P_HOR*G_OUTPUTS, 1);

			setsubmatrix(u_vec,u_Y_top,1+P_HOR*N_STATES + P_HOR*M_INPUTS + (P_HOR-C_HOR)*M_INPUTS, 1);
			setsubmatrix(u_vec, u_Y_bottom, 1+P_HOR*N_STATES + P_HOR*M_INPUTS + (P_HOR-C_HOR)*M_INPUTS + P_HOR*G_OUTPUTS, 1);
			
			setsubmatrix(l_vec,l_Y_top,1+P_HOR*N_STATES + P_HOR*M_INPUTS + (P_HOR-C_HOR)*M_INPUTS, 1);
			setsubmatrix(l_vec, l_Y_bottom, 1+P_HOR*N_STATES + P_HOR*M_INPUTS + (P_HOR-C_HOR)*M_INPUTS + P_HOR*G_OUTPUTS, 1);
			
			//Initialize solver
			initializeSolver();
		}


		//Constructor for input rate AND output constraints
		//Using SFINAE here to enable constructor only and only when input rate AND output constraints
		// are enabled through template argument CONSTRAINT_INPUTRATES and CONSTRAINT_OUTPUTS
		template<bool CRATES = CONSTRAINT_INPUTRATES, bool COUTPUTS = CONSTRAINT_OUTPUTS, std::enable_if_t<(CRATES && COUTPUTS), size_t> = 0>
		// template<std::enable_if_t<(CONSTRAINT_INPUTRATES), size_t> = 0>
		// template <typename = std::enable_if_t<CONSTRAINT_INPUTRATES && !CONSTRAINT_OUTPUTS>>
		LinearMPC(const ArcsMat<N_STATES,N_STATES>& A, const ArcsMat<N_STATES,M_INPUTS>& B, const ArcsMat<G_OUTPUTS,N_STATES>& C,
		 double w_u, double w_y, double w_du,
		 const ArcsMat<N_STATES,1>& x0, const ArcsMat<M_INPUTS,1>& u_z1, const ArcsMat<G_OUTPUTS*P_HOR,1>& Y_REF,
		 const ArcsMat<M_INPUTS,1>& u_min, const ArcsMat<M_INPUTS,1>& u_max,
		 const ArcsMat<M_INPUTS,1>& du_min, const ArcsMat<M_INPUTS,1>& du_max,
		 const ArcsMat<G_OUTPUTS,1>& y_min, const ArcsMat<G_OUTPUTS,1>& y_max):
		 LinearMPC(A, B, C, w_u, w_y, w_du, x0, u_z1, Y_REF, u_min, u_max)
		{
			printf("Testing: Entering constructor for when input rate and output constraints are activated \n");

			//TODO: I could have skipped the input rate constraints definition and reduce the amount of code by just
			//invoking the constructor specialized for input rate constraints.
			//However, the input rate constraints-only constructor is disabled inside this condition
			// (because both constraint bool flags are activated)
			//Maybe we should just remove SFINAE?

			//----------- Input rate constraints ---------------
			
			//Some necessary storage
			du_min_stored = du_min;
			du_max_stored = du_max;

			//Stack input rate constraints
			ArcsMat<M_INPUTS*P_HOR,M_INPUTS*P_HOR> D_DU;
			ArcsMat<M_INPUTS*P_HOR,1> l_DU;
			ArcsMat<M_INPUTS*P_HOR,1> u_DU;
			ArcsMat<M_INPUTS*P_HOR,1> b0;
			setsubmatrix(b0, u_z1, 1, 1);
			auto d_DU = eye<P_HOR,P_HOR>();
			for(size_t i=2; i<=P_HOR; i++){
				d_DU(i,i-1) = -1.0;
			}
			Kron(d_DU, eye<M_INPUTS,M_INPUTS>(), D_DU);		
			setsubmatrix(A_mat, D_DU, 1+P_HOR*N_STATES + P_HOR*M_INPUTS, P_HOR*N_STATES+1);
			setsubmatrix(l_DU,	Kron(ones<C_HOR,1>(),du_min)	,1,1);
			setsubmatrix(u_DU,	Kron(ones<C_HOR,1>(),du_max)	,1,1);	

			setsubmatrix(l_vec,	l_DU + b0,1+P_HOR*N_STATES + P_HOR*M_INPUTS,1);
			setsubmatrix(u_vec,	u_DU + b0,1+P_HOR*N_STATES + P_HOR*M_INPUTS,1);

			//----------- Output constraints ---------------

			//Stack output constraints
			ArcsMat<P_HOR*G_OUTPUTS,P_HOR*N_STATES> Mx_Y = Kron(eye<P_HOR,P_HOR>(),C);
			ArcsMat<P_HOR*G_OUTPUTS,1> u_Y_top = Kron(ones<P_HOR,1>(),y_max);		
			ArcsMat<P_HOR*G_OUTPUTS,1> l_Y_top(-OSQP_INFTY); 	
			ArcsMat<P_HOR*G_OUTPUTS,1> u_Y_bottom(OSQP_INFTY); 	
			ArcsMat<P_HOR*G_OUTPUTS,1> l_Y_bottom = Kron(ones<P_HOR,1>(),y_min);

			setsubmatrix(A_mat, Mx_Y, 1+P_HOR*N_STATES + 2*P_HOR*M_INPUTS, 1);
			setsubmatrix(A_mat, Mx_Y, 1+P_HOR*N_STATES + 2*P_HOR*M_INPUTS + P_HOR*G_OUTPUTS, 1);

			setsubmatrix(u_vec,u_Y_top,1+P_HOR*N_STATES + 2*P_HOR*M_INPUTS, 1);
			setsubmatrix(u_vec, u_Y_bottom, 1+P_HOR*N_STATES + 2*P_HOR*M_INPUTS + P_HOR*G_OUTPUTS, 1);
			
			setsubmatrix(l_vec,l_Y_top,1+P_HOR*N_STATES + 2*P_HOR*M_INPUTS, 1);
			setsubmatrix(l_vec, l_Y_bottom, 1+P_HOR*N_STATES + 2*P_HOR*M_INPUTS + P_HOR*G_OUTPUTS, 1);	

			//Initialize solver
			initializeSolver();		

		}

		//TODO: Will delete soon
		void checkMatrices(){
			printf("Test to check what happens to matrices after constructor \n");
			disp(P_mat);
			disp(q_vec);
			disp(l_vec);
			disp(u_vec);
			disp(A_mat);
			disp(A_stored);
			disp(du_min_stored);
			disp(du_max_stored);
		}

		void update(const ArcsMat<P_HOR*G_OUTPUTS,1>& Y_REF, const ArcsMat<N_STATES,1>& x0, const ArcsMat<M_INPUTS,1>& u_z1){

			OSQPInt exitflag = 0;

			//------------ update q vector -----------------
			
		
			ArcsMat<M_INPUTS*P_HOR,1> b0;
			setsubmatrix(b0, u_z1, 1, 1);
			setsubmatrix(q_vec, q_vec_r1*Y_REF, 1, 1);
			setsubmatrix(q_vec, q_vec_r2*b0, 1+P_HOR*N_STATES, 1);


			//----------- update l and u vectors ----------

			//Dynamics constraints
			ArcsMat<N_STATES,1> x0_computed = A_stored*x0;
			setsubmatrix(u_vec, x0_computed, 1, 1);
			setsubmatrix(l_vec, x0_computed, 1, 1);

			//Input rate constraints (if any)
			if constexpr(CONSTRAINT_INPUTRATES)
			{
				setsubmatrix(u_vec, du_max_stored + u_z1, 1+P_HOR*N_STATES + P_HOR*M_INPUTS, 1);
				setsubmatrix(l_vec, du_min_stored + u_z1, 1+P_HOR*N_STATES + P_HOR*M_INPUTS, 1);
			}

			//Update OSQP solver with new vectors
			exitflag = qpSolver.Update_Vecs(q_vec,l_vec,u_vec);
			arcs_assert(exitflag==0);

		}


		private:

		//Initialize OSQP solver with computed P, A matrices and q, l, u vectors
		//TODO: Add support for other solvers? (such as qpOASES)
		void initializeSolver()
		{
			qpSolver.initializeSolver(P_mat, A_mat, q_vec, l_vec, u_vec);
		}


		//Start OSQP solver with computed P, A matrices and q, l, u vectors
		void solve(ArcsMat<P_HOR*M_INPUTS,1>& U_opt, ArcsMat<P_HOR*N_STATES,1> X_predicted,
		 double& slack_var, OSQP_status& solver_status)
		{

		}

		//Calculates number of constraints for constraint vectors and matrix computation
		static constexpr std::size_t constraintsSize() {		

			//Parameter asserts
			static_assert(P_HOR>0, "LinearMPC: Prediction horizon must be positive");
			static_assert(C_HOR>0, "LinearMPC: Control horizon must be positive");
			static_assert(P_HOR-C_HOR>0, "LinearMPC: Control horizon must be smaller than prediction horizon");
			static_assert(N_STATES>0, "LinearMPC: Number of states must be positive");
			static_assert(M_INPUTS>0, "LinearMPC: Number of inputs must be positive");
			static_assert(G_OUTPUTS>0, "LinearMPC: Number of outputs must be positive");

			// std::size_t num_of_constraints = M_INPUTS*C_HOR + 1;	
			std::size_t num_of_constraints = P_HOR*(N_STATES + M_INPUTS) + 1;

			if(CONSTRAINT_INPUTRATES==false){
				num_of_constraints += (P_HOR-C_HOR)*M_INPUTS; //Just for blocking movement constraints
			}
			else{
				num_of_constraints += P_HOR*M_INPUTS; // For blocking movements AND input rate constraints
			}


			if(CONSTRAINT_OUTPUTS==true){
				num_of_constraints += 2*P_HOR*G_OUTPUTS;
			}

			return num_of_constraints;
    	}

		
		//Declaration of matrices for solver and solver itself
		const double SLACK_EPS = 1e5;


		ArcsMat<P_HOR*(N_STATES+M_INPUTS)+1,P_HOR*(N_STATES+M_INPUTS)+1> P_mat;
		ArcsMat<constraintsSize(),P_HOR*(N_STATES+M_INPUTS)+1> A_mat;
		ArcsMat<P_HOR*(N_STATES+M_INPUTS)+1,1> q_vec;
		ArcsMat<constraintsSize(),1> l_vec;
		ArcsMat<constraintsSize(),1> u_vec;		
		ArcsMat<P_HOR*N_STATES, P_HOR*G_OUTPUTS> q_vec_r1;
		ArcsMat<P_HOR*M_INPUTS, P_HOR*M_INPUTS> q_vec_r2;

		//TODO: du_min and du_max should be conditionally declared only when input rate constraints are enabled
		//Unfortunately I haven't figured how to do this using SFINAE or some other approach
		ArcsMat<N_STATES,N_STATES> A_stored;
		ArcsMat<M_INPUTS,1> du_min_stored;
		ArcsMat<M_INPUTS,1> du_max_stored;
		OSQP_Solver<P_HOR*(N_STATES+M_INPUTS)+1,constraintsSize()> qpSolver;



    };


}

#endif
