//! @file LinearMPC.cc
//! @brief Linear MPC implementation class
//!
//! Linear MPC based on OSQP solver using sparse formulation 
//! (Try to avoid condensed formulation for faster solving speed)
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
#include <array>
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

//! @brief Class for implementing linear MPC with box constraints for LTI systems
//! @tparam	N_STATES	Number of plant states (Dimension of x vector)
//! @tparam M_INPUTS	Number of plant inputs (Dimension of u vector)
//! @tparam G_OUTPUTS	Number of plant outputs to be controlled (Dimension of y vector)
//! @tparam P_HOR	Prediction horizon
//! @tparam C_HOR	Control horizon, must be smaller than prediction horizon
//! @tparam CONSTRAINT_INPUTRATES	Flag to enable input rate constraints
//! @tparam CONSTRAINT_OUTPUTS	Flag to enable output constraints
template <size_t N_STATES, size_t M_INPUTS, size_t G_OUTPUTS, size_t P_HOR, size_t C_HOR, bool CONSTRAINT_INPUTRATES = false, bool CONSTRAINT_OUTPUTS = false>
    class LinearMPC{

		public:

		//! @brief  Main constructor, called mainly when input rate and output constraints are disabled
		//! @param[in] A Plant model discrete-time state matrix (A matrix)
        //! @param[in] B Plant model discrete-time input matrix (B matrix)
        //! @param[in] C Plant model output matrix (C matrix)
        //! @param[in] w_u Input cost weight
        //! @param[in] w_y Output cost weight
        //! @param[in] w_du Input rate cost weight
		//! @param[in] x0 Initial plant state
        //! @param[in] u_z1 Previous sample input to plant
        //! @param[in] Y_REF Reference vector with references from k=0 to k=P_HOR stacked in a single vector
		//! @param[in] u_min Vector with lower bounds values for each input
        //! @param[in] u_max Vector with upper bounds values for each input
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

			//printf("Testing: Entering constructor \n");




			arcs_assert(w_u >= 0);
			arcs_assert(w_du >= 0);
			arcs_assert(w_y >= 0);


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
			P_mat = 2*P_mat;	//To eliminate the default 1/2 multiplier of OSQP solver (is this needed?)
			// disp(P_mat);


			
			//------------ q vector -----------------
			
			//TODO: Fix submatrix sizes
			ArcsMat<M_INPUTS*P_HOR,1> b0;
			setsubmatrix(b0, u_z1, 1, 1);
			q_vec_r1 = -2*tp(CBAR)*tp(QY_pre);
			q_vec_r2 = -2*tp(D_DU)*tp(QDU_pre);	
			// disp(D_DU);
			// disp(QDU);
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


		//! @brief  Constructor for case when only input rate constraints are enabled
		//! @param[in] A Plant model discrete-time state matrix (A matrix)
        //! @param[in] B Plant model discrete-time input matrix (B matrix)
        //! @param[in] C Plant model output matrix (C matrix)
        //! @param[in] w_u Input cost weight
        //! @param[in] w_y Output cost weight
        //! @param[in] w_du Input rate cost weight
		//! @param[in] x0 Initial plant state
        //! @param[in] u_z1 Previous sample input to plant
        //! @param[in] Y_REF Reference vector with references from k=0 to k=P_HOR stacked in a single vector
		//! @param[in] u_min Vector with lower bounds values for each input
        //! @param[in] u_max Vector with upper bounds values for each input
		//! @param[in] du_min Vector with lower bounds for the change of rate of each input
        //! @param[in] du_max Vector with upper bounds for the change of rate of each input
		template<bool CRATES = CONSTRAINT_INPUTRATES, bool COUTPUTS = CONSTRAINT_OUTPUTS, std::enable_if_t<(CRATES && !COUTPUTS), size_t> = 0>
		LinearMPC(const ArcsMat<N_STATES,N_STATES>& A, const ArcsMat<N_STATES,M_INPUTS>& B, const ArcsMat<G_OUTPUTS,N_STATES>& C,
		 double w_u, double w_y, double w_du,
		 const ArcsMat<N_STATES,1>& x0, const ArcsMat<M_INPUTS,1>& u_z1, const ArcsMat<G_OUTPUTS*P_HOR,1>& Y_REF,
		 const ArcsMat<M_INPUTS,1>& u_min, const ArcsMat<M_INPUTS,1>& u_max,
		 const ArcsMat<M_INPUTS,1>& du_min, const ArcsMat<M_INPUTS,1>& du_max):
		 LinearMPC(A, B, C, w_u, w_y, w_du, x0, u_z1, Y_REF, u_min, u_max)
		{
			//printf("Testing: Entering constructor for when input rates constraints are activated \n");

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


		//! @brief  Constructor for case when only output constraints are enabled
		//! @param[in] A Plant model discrete-time state matrix (A matrix)
        //! @param[in] B Plant model discrete-time input matrix (B matrix)
        //! @param[in] C Plant model output matrix (C matrix)
        //! @param[in] w_u Input cost weight
        //! @param[in] w_y Output cost weight
        //! @param[in] w_du Input rate cost weight
		//! @param[in] x0 Initial plant state
        //! @param[in] u_z1 Previous sample input to plant
        //! @param[in] Y_REF Reference vector with references from k=0 to k=P_HOR stacked in a single vector
		//! @param[in] u_min Vector with lower bounds values for each input
        //! @param[in] u_max Vector with upper bounds values for each input
		//! @param[in] y_min Vector with lower bounds for each output state (Note: soft-constrained)
        //! @param[in] y_max Vector with upper bounds for each output state (Note: soft-constrained)
		template<bool CRATES = CONSTRAINT_INPUTRATES, bool COUTPUTS = CONSTRAINT_OUTPUTS, std::enable_if_t<(!CRATES && COUTPUTS), size_t> = 0>
		LinearMPC(const ArcsMat<N_STATES,N_STATES>& A, const ArcsMat<N_STATES,M_INPUTS>& B, const ArcsMat<G_OUTPUTS,N_STATES>& C,
		 double w_u, double w_y, double w_du,
		 const ArcsMat<N_STATES,1>& x0, const ArcsMat<M_INPUTS,1>& u_z1, const ArcsMat<G_OUTPUTS*P_HOR,1>& Y_REF,
		 const ArcsMat<M_INPUTS,1>& u_min, const ArcsMat<M_INPUTS,1>& u_max,
		 const ArcsMat<G_OUTPUTS,1>& y_min, const ArcsMat<G_OUTPUTS,1>& y_max):
		 LinearMPC(A, B, C, w_u, w_y, w_du, x0, u_z1, Y_REF, u_min, u_max)
		{
			//printf("Testing: Entering constructor for when output constraints are activated \n");



			//Call basic constructor
			// (This doesn't work as it only creates a local copy that is deleted)
			// Basic constructor needs to be called through constructor delegation as done above
			// LinearMPC(A, B, C, w_u, w_y, w_du, x0, u_z1, Y_REF, u_min, u_max);

			//Stack output constraints
			ArcsMat<P_HOR*G_OUTPUTS,P_HOR*N_STATES> Mx_Y = Kron(eye<P_HOR,P_HOR>(),C);
			ArcsMat<P_HOR*G_OUTPUTS,1> Meps_Y(1);
			ArcsMat<P_HOR*G_OUTPUTS,1> u_Y_top = Kron(ones<P_HOR,1>(),y_max);		
			ArcsMat<P_HOR*G_OUTPUTS,1> l_Y_top(-OSQP_INFTY); 	
			ArcsMat<P_HOR*G_OUTPUTS,1> u_Y_bottom(OSQP_INFTY); 	
			ArcsMat<P_HOR*G_OUTPUTS,1> l_Y_bottom = Kron(ones<P_HOR,1>(),y_min);

			setsubmatrix(A_mat, Mx_Y, 1+P_HOR*N_STATES + P_HOR*M_INPUTS + (P_HOR-C_HOR)*M_INPUTS, 1);
			setsubmatrix(A_mat, Meps_Y, 1+P_HOR*N_STATES + P_HOR*M_INPUTS + (P_HOR-C_HOR)*M_INPUTS, 1+P_HOR*(M_INPUTS+N_STATES));
			setsubmatrix(A_mat, Mx_Y, 1+P_HOR*N_STATES + P_HOR*M_INPUTS + (P_HOR-C_HOR)*M_INPUTS + P_HOR*G_OUTPUTS, 1);
			setsubmatrix(A_mat, Meps_Y, 1+P_HOR*N_STATES + P_HOR*M_INPUTS + (P_HOR-C_HOR)*M_INPUTS + P_HOR*G_OUTPUTS, 1+P_HOR*(M_INPUTS+N_STATES));

			setsubmatrix(u_vec,u_Y_top,1+P_HOR*N_STATES + P_HOR*M_INPUTS + (P_HOR-C_HOR)*M_INPUTS, 1);
			setsubmatrix(u_vec, u_Y_bottom, 1+P_HOR*N_STATES + P_HOR*M_INPUTS + (P_HOR-C_HOR)*M_INPUTS + P_HOR*G_OUTPUTS, 1);
			
			setsubmatrix(l_vec,l_Y_top,1+P_HOR*N_STATES + P_HOR*M_INPUTS + (P_HOR-C_HOR)*M_INPUTS, 1);
			setsubmatrix(l_vec, l_Y_bottom, 1+P_HOR*N_STATES + P_HOR*M_INPUTS + (P_HOR-C_HOR)*M_INPUTS + P_HOR*G_OUTPUTS, 1);
			
			//Initialize solver
			initializeSolver();
		}


		//! @brief  Constructor for case when only output constraints are enabled
		//! @param[in] A Plant model discrete-time state matrix (A matrix)
        //! @param[in] B Plant model discrete-time input matrix (B matrix)
        //! @param[in] C Plant model output matrix (C matrix)
        //! @param[in] w_u Input cost weight
        //! @param[in] w_y Output cost weight
        //! @param[in] w_du Input rate cost weight
		//! @param[in] x0 Initial plant state
        //! @param[in] u_z1 Previous sample input to plant
        //! @param[in] Y_REF Reference vector with references from k=0 to k=P_HOR stacked in a single vector
		//! @param[in] u_min Vector with lower bounds values for each input
        //! @param[in] u_max Vector with upper bounds values for each input
		//! @param[in] du_min Vector with lower bounds for the change of rate of each input
        //! @param[in] du_max Vector with upper bounds for the change of rate of each input
		//! @param[in] y_min Vector with lower bounds for each output state (Note: soft-constrained)
        //! @param[in] y_max Vector with upper bounds for each output state (Note: soft-constrained)
		template<bool CRATES = CONSTRAINT_INPUTRATES, bool COUTPUTS = CONSTRAINT_OUTPUTS, std::enable_if_t<(CRATES && COUTPUTS), size_t> = 0>
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
			ArcsMat<P_HOR*G_OUTPUTS,1> Meps_Y(1);
			ArcsMat<P_HOR*G_OUTPUTS,1> u_Y_top = Kron(ones<P_HOR,1>(),y_max);		
			ArcsMat<P_HOR*G_OUTPUTS,1> l_Y_top(-OSQP_INFTY); 	
			ArcsMat<P_HOR*G_OUTPUTS,1> u_Y_bottom(OSQP_INFTY); 	
			ArcsMat<P_HOR*G_OUTPUTS,1> l_Y_bottom = Kron(ones<P_HOR,1>(),y_min);

			setsubmatrix(A_mat, Mx_Y, 1+P_HOR*N_STATES + 2*P_HOR*M_INPUTS, 1);
			setsubmatrix(A_mat, Meps_Y, 1+P_HOR*N_STATES + 2*P_HOR*M_INPUTS, 1+P_HOR*(M_INPUTS+N_STATES));
			setsubmatrix(A_mat, Mx_Y, 1+P_HOR*N_STATES + 2*P_HOR*M_INPUTS + P_HOR*G_OUTPUTS, 1);
			setsubmatrix(A_mat, Meps_Y, 1+P_HOR*N_STATES + 2*P_HOR*M_INPUTS + P_HOR*G_OUTPUTS, 1+P_HOR*(M_INPUTS+N_STATES));

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

		//TODO: Will delete soon
		void testOutputToMAT(const std::string fileName)
		{
			std::string outputString = "LinearMPC: Printing matrices to specified mat file: " + fileName;
			// printf("LinearMPC: Printing matrices to specified mat file \n");
			printf("%s \n",outputString.c_str());
			MatExport MatFile1(fileName);
			MatFile1.Save("Amat_exp", A_mat);
			MatFile1.Save("Pmat_exp", P_mat);
			MatFile1.Save("qvec_exp", q_vec);
			MatFile1.Save("lvec_exp", l_vec);
			MatFile1.Save("uvec_exp", u_vec);
		}

		//! @brief Updates MPC internal quadratic program with new reference vector, new initial state and previous input
		//! @param[in]	Y_REF	New stacked reference vector
		//! @param[in]	x0	New initial state
		//! @param[in]	u_z1	New previous sample input
		void update(const ArcsMat<P_HOR*G_OUTPUTS,1>& Y_REF, const ArcsMat<N_STATES,1>& x0, const ArcsMat<M_INPUTS,1>& u_z1){

			OSQPInt exitflag = 0;

			//------------ update q vector -----------------
			
		
			ArcsMat<M_INPUTS*P_HOR,1> b0;
			setsubmatrix(b0, u_z1, 1, 1);
			setsubmatrix(q_vec, q_vec_r1*Y_REF, 1, 1);
			//disp(q_vec_r2*b0);
			// disp(q_vec_r2);
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

		//! @brief Solves quadratic program (Returns by reference)
		//! @param[in]	U_opt	Vector to store optimal solution vector 
		//! @param[in]	X_predicted Vector to store predicted states
		//! @param[in]	slack_var Reference to store slack amount
		//! @param[in]	solver_status OSQP_Status object reference to store solver status
		void solve(ArcsMat<P_HOR*M_INPUTS,1>& U_opt, ArcsMat<P_HOR*N_STATES,1> X_predicted,
		 double& slack_var, OSQP_Status& solver_status)
		{
			std::array<OSQPFloat,P_HOR*(N_STATES+M_INPUTS)+1> solution_array;
			std::array<OSQPFloat,P_HOR*M_INPUTS> U_array;
			std::array<OSQPFloat,P_HOR*N_STATES> X_array;			
			OSQPInt exitflag = 0;

			//Solve quadratic program and verify that a feasible solution was obtained
			exitflag = qpSolver.solve(solution_array);
			arcs_assert(exitflag==0);

			//Populate U, X arrays and slack variable
			for(size_t i=0; i<P_HOR*N_STATES; i++)
			{
				X_array[i] = solution_array[i];
			}
			
			for(size_t i=P_HOR*N_STATES; i<P_HOR*(N_STATES+M_INPUTS); i++)
			{
				U_array[i-P_HOR*N_STATES] = solution_array[i];
			}

			slack_var = solution_array[P_HOR*(N_STATES+M_INPUTS)];

			//Load arrays into respective matrices and get also solver status
			U_opt.LoadArray(U_array);
			X_predicted.LoadArray(X_array);
			solver_status = qpSolver.getSolverStatus();		

		}

		//! @brief Enables/disables internal solver warm start
		//! @param[in]	warmStartFlag Flag to enable/disable solver warm start
		void setWarmStart(const bool warmStartFlag)
		{
			OSQPInt exitflag = 0;

			//Enable/disable warm start and verify that no error occurred
			exitflag = qpSolver.setWarmStart(warmStartFlag);
			arcs_assert(exitflag==0);
		}

		//! @brief Sets maximum number of solver iterations to solve quadratic program
		//! @param[in]	maxIterations Number of max. iterations
		void setMaxIterationsSolver(const OSQPInt maxIterations)
		{
			OSQPInt exitflag = 0;

			//Set number of max iterations for solver and verify that no error occurred
			exitflag = qpSolver.setMaxIterations(maxIterations);
			arcs_assert(exitflag==0);
		}


		private:

		//! @brief Initializes internal solver (OSQP solver)
		void initializeSolver()
		{
			qpSolver.initializeSolver(P_mat, A_mat, q_vec, l_vec, u_vec);
		}



		//! @brief Meta function for calculating size of constraint vectors
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
		const double SLACK_EPS = 1e5;	//I think this should be made variable?


		ArcsMat<P_HOR*(N_STATES+M_INPUTS)+1,P_HOR*(N_STATES+M_INPUTS)+1> P_mat; //Hessian matrix
		ArcsMat<constraintsSize(),P_HOR*(N_STATES+M_INPUTS)+1> A_mat;	//Constraint matrix
		ArcsMat<P_HOR*(N_STATES+M_INPUTS)+1,1> q_vec;	//Gradient vector
		ArcsMat<constraintsSize(),1> l_vec;	//Lower bound vector
		ArcsMat<constraintsSize(),1> u_vec;	//Upper bound vector		
		ArcsMat<P_HOR*N_STATES, P_HOR*G_OUTPUTS> q_vec_r1;
		ArcsMat<P_HOR*M_INPUTS, P_HOR*M_INPUTS> q_vec_r2;

		//TODO: du_min and du_max should be conditionally declared only when input rate constraints are enabled
		//Unfortunately I haven't figured how to do this using SFINAE or some other approach
		//So we declare them as class members independently of whether input rate constraints are enabled or not
		ArcsMat<N_STATES,N_STATES> A_stored;
		ArcsMat<M_INPUTS,1> du_min_stored;
		ArcsMat<M_INPUTS,1> du_max_stored;
		OSQP_Solver<P_HOR*(N_STATES+M_INPUTS)+1,constraintsSize()> qpSolver;

    };


}

#endif
