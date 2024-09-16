//! @file OSQP_Solver.hh
//! @brief C++ ARCS wrapper for OSQP solver (Template class)
//!
//! C++ wrapper for using OSQP solver on ARCS
//! Solves the quadratic program
//!     minimize    0.5*x'*P*x + q'*x
//!     subject to      l <= Ax <= u
//!                 
//! @date 2024/9/9
//! @author Juan Padron
//
// Copyright (C) 2011-2024 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.


//TODO Add positive definiteness (and thus convexity) check for P matrix
//TODO Add lower and upper bound checking (Upper bounds must be larger than lower bounds)
//TODO Add error code return for solve() function
//TODO Add functions for changing specific (not all) settings. Which ones are most important? Gotta check documentation
//TODO Add function for observing solver status

#ifndef OSQPSOLVER_ARCS
#define OSQPSOLVER_ARCS

#include <cstdio>
#include <memory>
#include <osqp.h>
#include <functional>
#include <array>
#include <vector>
#include <cmath>
#include <cassert>
#include "ArcsMatrix.hh"

// ARCS組込み用マクロ
#ifdef ARCS_IN
	// ARCSに組み込まれる場合
	#include "ARCSassert.hh"
#else
	// ARCSに組み込まれない場合
	#define arcs_assert(a) (assert(a))
#endif

namespace ARCS{

//! @brief OSQP solver wrapper class
//! @tparam	N_VARS          Number of variables to optimize
//! @tparam M_CONSTRAINTS	Number of constraints
template <size_t N_VARS, size_t M_CONSTRAINTS>
    class OSQP_Solver{


        public:



        //! @brief OSQP_Solver class constructor
        //! @param[in] Pmat	   Hessian matrix P (quadratic cost) - Must be positive semidefinite symmetric
        //! @param[in] Amat    Constraint matrix A
        //! @param[in] qVec    Linear cost vector q
        //! @param[in] lVec    Constraint lower bound vector l
        //! @param[in] uVec    Constraint upper bound vector u - Its values must be greater than lVec
        OSQP_Solver(const ArcsMat<N_VARS,N_VARS>& Pmat, const ArcsMat<M_CONSTRAINTS,N_VARS>& Amat, const ArcsMat<N_VARS,1>& qVec, const ArcsMat<M_CONSTRAINTS,1>& lVec, const ArcsMat<M_CONSTRAINTS,1>& uVec):
            solver(nullptr, OSQPSolver_deleter),
            settings(nullptr, OSQPSettings_deleter),
            P(nullptr, OSQPCscMatrix_deleter),
            A(nullptr, OSQPCscMatrix_deleter),
            q(),
            l(),
            u(),
            P_elements(),
            P_rows(),
            P_cols(),
            P_nnz(0),
            A_elements(),
            A_rows(),
            A_cols(),
            A_nnz(0)
        {


            OSQPInt exitflag;

            //Checks if P is symmetric and positive definite
            arcs_assert(CheckSymmetryPositiveDef(Pmat));
            

            //Convert [q,l,u] ArcsMat vectors to arrays, [P,A] ArcsMat matrices to CSC format
            qVec.StoreArray(q);
            lVec.StoreArray(l);
            uVec.StoreArray(u);
            convertPtoCsc(Pmat);
            convertAtoCsc(Amat);

            OSQPSettings* tmp_settings = static_cast<OSQPSettings*>(malloc(sizeof(OSQPSettings)));
            OSQPCscMatrix* tmp_P = static_cast<OSQPCscMatrix*>(malloc(sizeof(OSQPCscMatrix)));
            OSQPCscMatrix* tmp_A = static_cast<OSQPCscMatrix*>(malloc(sizeof(OSQPCscMatrix)));
            OSQPSolver* tmp_solver = static_cast<OSQPSolver*>(malloc(sizeof(OSQPSolver)));


            P.reset(tmp_P);
            A.reset(tmp_A);
            settings.reset(tmp_settings);   

            //Set OSQP default settings but enable polishing and disable verbose
            osqp_set_default_settings(settings.get());
            settings->polishing = 1;
            settings->verbose = 0;
            //settings->eps_abs = 1e-5;
            //settings->eps_rel = 1e-5;
            //settings->max_iter = 10000;

            
            //Load matrix data on the OSQPCscMatrix structs P and A using OSQP csc_set_data
            
            csc_set_data(P.get(), N_VARS, N_VARS, P_nnz, P_elements.data(), P_rows.data(), P_cols.data());
            csc_set_data(A.get(), M_CONSTRAINTS, N_VARS, A_nnz, A_elements.data(), A_rows.data(), A_cols.data());

            //Setup OSQP solver and store on temp. pointer 
            exitflag = osqp_setup(&tmp_solver, P.get(), q.data(), A.get(), l.data(), u.data(), M_CONSTRAINTS, N_VARS, settings.get());
            
            //Check if there is 
            arcs_assert(exitflag == 0);
            solver.reset(tmp_solver);
            
      
            /*
            disp(Pmat);
            convertPtoCsc(Pmat);
            testPrintCscVectors(P_elements,P_rows,P_cols,P_nnz);
            disp(Amat);
            convertAtoCsc(Amat);
            testPrintCscVectors(A_elements,A_rows,A_cols,A_nnz);
            */
            
            printf("OSQP_Solver instantiation created \n");
        }

        




        //! @brief OSQP_Solver class destructor
        ~OSQP_Solver(void)
        {
            printf("OSQP_Solver instantiation destroyed \n");
        }

        //! @brief  Minimizes 0.5*x'*P*x + q'*x subject to constraint l <= Ax <= u
        //! @param[out] solution_array	Optimal solution x* as an std::array
        //! @return OSQP exitflag
        OSQPInt solve(std::array<OSQPFloat,N_VARS>& solution_array){
            OSQPInt exitflag;

            exitflag = osqp_solve(solver.get());
            OSQPFloat* sol_vector;

            if(exitflag == 0){
                sol_vector = solver->solution->x;
                for(size_t i = 0; i<N_VARS; ++i) solution_array[i] = sol_vector[i];
            }

            return exitflag;
        }


        //! @brief  Updates linear cost vector q
        //! @param[in] qVec New linear cost vector
        //! @return OSQP exitflag
        OSQPInt Update_qVec(const ArcsMat<N_VARS,1>& qVec){
            OSQPInt exitflag;

            //Check if solver has been initialized (if its pointer is not null)
            //arcs_assert(static_cast<bool>(solver));     
            qVec.StoreArray(q);
            exitflag = osqp_update_data_vec(solver.get(), q.data(), OSQP_NULL, OSQP_NULL);

            return exitflag;            
        }

        //! @brief  Updates constraint lower bound vector l
        //! @param[in] lVec New constraint lower bound vector
        //! @return OSQP exitflag
        OSQPInt Update_lVec(const ArcsMat<M_CONSTRAINTS,1>& lVec)
        {
            OSQPInt exitflag;

            //Check if solver has been initialized (if its pointer is not null)
            //arcs_assert(static_cast<bool>(solver));     
            lVec.StoreArray(l);
            exitflag = osqp_update_data_vec(solver.get(), OSQP_NULL, l.data(), OSQP_NULL);

            return exitflag;            
        }

        //! @brief  Updates constraint upper bound vector u
        //! @param[in] uVec New constraint upper bound vector
        //! @return OSQP exitflag
        OSQPInt Update_uVec(const ArcsMat<M_CONSTRAINTS,1>& uVec)
        {
            OSQPInt exitflag;

            //Check if solver has been initialized (if its pointer is not null)
            //arcs_assert(static_cast<bool>(solver));    
            uVec.StoreArray(u);
            exitflag = osqp_update_data_vec(solver.get(), OSQP_NULL, OSQP_NULL, u.data());

            return exitflag;            
        }

        //! @brief  Updates linear cost vector 
        //! @param[in] qVec New linear cost q vector
        //! @param[in] lVec New constraint lower bound l vector
        //! @param[in] uVec New constraint upper bound u vector
        //! @return OSQP exitflag
        OSQPInt Update_Vecs(const ArcsMat<N_VARS,1>& qVec, const ArcsMat<M_CONSTRAINTS,1>& lVec, const ArcsMat<M_CONSTRAINTS,1>& uVec)
        {
            OSQPInt exitflag;

            //Check if solver has been initialized (if its pointer is not null)
            //arcs_assert(static_cast<bool>(solver));    
            qVec.StoreArray(q);
            lVec.StoreArray(l);
            uVec.StoreArray(u);
            exitflag = osqp_update_data_vec(solver.get(), q.data(), l.data(), u.data());

            return exitflag;    
        }

        //! @brief  Updates hessian matrix (quadratic cost)
        //! @param[in] P_new New Hessian P matrix
        //! @return OSQP exitflag
        OSQPInt Update_PMatrix(const ArcsMat<N_VARS,N_VARS>& P_new)
        {
            OSQPInt exitflag;
            //Check first if matrix is symmetric and positive definite
            arcs_assert(CheckSymmetryPositiveDef(P_new));

            //Next, convert matrix to CSC format                                      


            //Temporary storage for our CSC vectors
            int cols = N_VARS;
            std::vector<OSQPFloat> P_elements_tmp;
            std::vector<OSQPInt> P_rows_tmp;
            std::vector<OSQPInt> P_cols_tmp;
            OSQPInt P_nnz_tmp = 0;



            // Initialize the P_cols array
            P_cols_tmp.resize(cols + 1, 0);

            //CSC conversion
            
            for (int col = 0; col < cols; ++col) {
                P_cols_tmp[col] = P_elements_tmp.size();  // Mark the start index of the column

                for (int row = 0; row < col+1; ++row) {
                    OSQPFloat val = P_new(row+1,col+1); //ArcsMat works with 1-base indexing                  


                    if (std::abs(val)>zero_eps) {                        
                        P_elements_tmp.push_back(val);            // Store non-zero value
                        P_rows_tmp.push_back(row);       // Store the corresponding row index
                        P_nnz_tmp++; //Count number of non-zero elements
                    }
                    
                }
            }
            P_cols_tmp[cols]=P_elements_tmp.size(); //Needed for filling last indexing element  

            //We compare now if the sparsity pattern of P_new is the same as P by checking row and col indexes
            arcs_assert(P_rows_tmp==P_rows);
            arcs_assert(P_cols_tmp==P_cols);
            arcs_assert(P_nnz_tmp==P_nnz);

            //As sparsity pattern is the same, copy only vector with new elements
            P_elements = P_elements_tmp;

            //Finally, update the P matrix data with our new CSC vectors
            exitflag = osqp_update_data_mat(solver.get(),
                                        P_elements.data(), OSQP_NULL, P_nnz,
                                        OSQP_NULL, OSQP_NULL, 0);

            return exitflag;            

        }

        //! @brief  Updates constraint matrix A 
        //! @param[in] A_new New constraint matrix
        //! @return OSQP exitflag
        OSQPInt Update_AMatrix(const ArcsMat<M_CONSTRAINTS,N_VARS>& A_new)
        {
            OSQPInt exitflag;

            //Convert new A matrix to CSC format                                      


            int rows = M_CONSTRAINTS;
            int cols = N_VARS;
            //Temporary storage for our CSC vectors
            std::vector<OSQPFloat> A_elements_tmp;
            std::vector<OSQPInt> A_rows_tmp;
            std::vector<OSQPInt> A_cols_tmp;
            OSQPInt A_nnz_tmp = 0;


            //Clear all elements of our CSC arrays
            //We do this because the A matrix may change during runtime and thus we should empty the vector structs beforehand

            // Initialize the A_cols array
            A_cols_tmp.resize(cols + 1, 0);

            // Iterate over each column
            
            for (int col = 0; col < cols; ++col) {
                A_cols_tmp[col] = A_elements_tmp.size();  // Mark the start index of the column

                for (int row = 0; row < rows; ++row) {
                    OSQPFloat val = A_new(row+1,col+1); //ArcsMat works with 1-base indexing                  

                    if (std::abs(val)>zero_eps) {                        
                        A_elements_tmp.push_back(val);            // Store non-zero value
                        A_rows_tmp.push_back(row);       // Store the corresponding row index
                        A_nnz_tmp++; //Count number of non-zero elements
                    }
                    
                }
            }
            A_cols_tmp[cols]=A_elements_tmp.size(); //Needed for filling last indexing element 

            //We compare now if the sparsity pattern of P_new is the same as P by checking row and col indexes
            arcs_assert(A_rows_tmp==A_rows);
            arcs_assert(A_cols_tmp==A_cols);
            arcs_assert(A_nnz_tmp==A_nnz);

            //As sparsity pattern is the same, copy only vector with new elements
            A_elements = A_elements_tmp;

            //Finally, update the P matrix data with our new CSC vectors
            exitflag = osqp_update_data_mat(solver.get(),
                                        OSQP_NULL, OSQP_NULL, 0,
                                        A_elements.data(), OSQP_NULL, A_nnz);

            return exitflag;            

        }




        private:

        //! @brief  Checks if matrix is symmetric and positive definite
        //! @param[in] Pcheck <NVARSxNVARS> square matrix
        //! @return true if symmetric and pos. def., elsewhere false
        bool CheckSymmetryPositiveDef(const ArcsMat<N_VARS,N_VARS>& Pcheck)
        {
            
            //Check first if matrix is symmetric
            //Iterate on upper triangular part (Pcheck is indexed from 1 as per ArcsMat def)
            for (size_t i = 0; i < N_VARS; ++i) {
                for (size_t j = i + 1; j < N_VARS; ++j) { 
                    if (Pcheck(i+1,j+1) != Pcheck(j+1,i+1)){
                        return false;
                    }
                }
            }


            //If symmetry check is passed, then check if its positive definite
            //TODO: Implement positive definiteness check (Too lazy to do now...)
            return true;
        }

        //! @brief  Converts Hessian matrix to CSC (Compressed Sparse Column) format
        //! @param[in] Pmat <N_VARSxN_VARS> square matrix in ArcsMat type
        void convertPtoCsc(const ArcsMat<N_VARS,N_VARS>& Pmat)
        {                       
            int cols = N_VARS;


            //Clear all elements of our CSC arrays
            //We do this because the P matrix may change during runtime and thus we should empty the vector structs beforehand
            P_elements.clear();
            P_rows.clear();
            P_cols.clear();
            P_nnz=0;



            // Initialize the P_cols array
            P_cols.resize(cols + 1, 0);

            // Iterate over each row whose index is greater or equal than the column index
            // This is in order to get only the upper triangular elements as specified by OSQP C library
            
            for (int col = 0; col < cols; ++col) {
                P_cols[col] = P_elements.size();  // Mark the start index of the column

                for (int row = 0; row < col+1; ++row) {
                    OSQPFloat val = Pmat(row+1,col+1); //ArcsMat works with 1-base indexing                  

                    if (std::abs(val)>zero_eps) {                        
                        P_elements.push_back(val);            // Store non-zero value
                        P_rows.push_back(row);       // Store the corresponding row index
                        P_nnz++; //Count number of non-zero elements
                    }
                    
                }
            }
            P_cols[cols]=P_elements.size(); //Needed for filling last indexing element               
              
        }

        //! @brief  Converts constraint matrix A to CSC (Compressed Sparse Column) format
        //! @param[in] Amat Constraint matrix in ArcsMat type
        void convertAtoCsc(const ArcsMat<M_CONSTRAINTS,N_VARS>& Amat)
        {
            int rows = M_CONSTRAINTS;
            int cols = N_VARS;


            //Clear all elements of our CSC arrays
            //We do this because the A matrix may change during runtime and thus we should empty the vector structs beforehand
            A_elements.clear();
            A_rows.clear();
            A_cols.clear();
            A_nnz=0;

            // Initialize the P_cols array
            A_cols.resize(cols + 1, 0);

            // Iterate over each column
            
            for (int col = 0; col < cols; ++col) {
                A_cols[col] = A_elements.size();  // Mark the start index of the column

                for (int row = 0; row < rows; ++row) {
                    OSQPFloat val = Amat(row+1,col+1); //ArcsMat works with 1-base indexing                  

                    if (std::abs(val)>zero_eps) {                        
                        A_elements.push_back(val);            // Store non-zero value
                        A_rows.push_back(row);       // Store the corresponding row index
                        A_nnz++; //Count number of non-zero elements
                    }
                    
                }
            }
            A_cols[cols]=A_elements.size(); //Needed for filling last indexing element             

            
        }

        //! @brief  Only for internal testing, prints CSC vectors associated to a matrix
        //! @param[in] values std::vector object containing elements of matrix
        //! @param[in] rows std::vector object containing row indexes
        //! @param[in] cols std::vector object containing column indexes
        //! @param[in] nnz  Number of non-zero elements
        static void testPrintCscVectors(std::vector<OSQPFloat> values, std::vector<OSQPInt> rows, std::vector<OSQPInt> cols, OSQPInt nnz)
        {
            printf("Values: \n");
            printf("[ ");
            for (size_t i = 0; i < values.size(); ++i) {
                printf(" %f ,", values[i]);  // Print the rest of the elements with spaces
            }
            printf(" ] \n");

            printf("Rows: \n");
            printf("[ ");
            for (size_t i = 0; i < rows.size(); ++i) {
                printf(" %lld ,", rows[i]);  // Print the rest of the elements with spaces
            }
            printf(" ] \n");

            printf("Columns: \n");
            printf("[ ");
            for (size_t i = 0; i < cols.size(); ++i) {
                printf(" %lld ,", cols[i]);  // Print the rest of the elements with spaces
            }
            printf(" ] \n");

            printf("Number of non-zero elements: %lld \n",nnz);

        }

        //! @brief  Deleter function for OSQPCscMatrix* type unique_ptr
        static void OSQPCscMatrix_deleter(OSQPCscMatrix* mat)
        {
            free(mat);
            printf("CsCMatrix ptr destroyed \n");
        }

        //! @brief  Deleter function for OSQPSolver* type unique_ptr
        static void OSQPSolver_deleter(OSQPSolver* slv)
        {
            osqp_cleanup(slv);
            printf("Solver ptr destroyed \n");
        }

        //! @brief  Deleter function for OSQPSettings* type unique_ptr
        static void OSQPSettings_deleter(OSQPSettings* stt)
        {
            free(stt);
            printf("Settings ptr destroyed \n");
        }


        //Member constants and variables
        const double zero_eps = 1e-16;  //For zeroing

        
        std::unique_ptr<OSQPSolver, std::function<void(OSQPSolver*)>> solver;
        std::unique_ptr<OSQPSettings, std::function<void(OSQPSettings*)>> settings;
        std::unique_ptr<OSQPCscMatrix, std::function<void(OSQPCscMatrix*)>> P;
        std::unique_ptr<OSQPCscMatrix, std::function<void(OSQPCscMatrix*)>> A;
        std::array<OSQPFloat,N_VARS> q;
        std::array<OSQPFloat,M_CONSTRAINTS> l;
        std::array<OSQPFloat,M_CONSTRAINTS> u;   
        std::vector<OSQPFloat> P_elements;
        std::vector<OSQPInt> P_rows;
        std::vector<OSQPInt> P_cols;
        OSQPInt P_nnz;
        std::vector<OSQPFloat> A_elements;
        std::vector<OSQPInt> A_rows;
        std::vector<OSQPInt> A_cols;
        OSQPInt A_nnz;

          
    


    };

}

#endif