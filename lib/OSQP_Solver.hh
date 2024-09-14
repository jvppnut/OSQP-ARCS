//! @file OSQP_Solver.hh
//! @brief C++ ARCS wrapper for OSQP solver (Template class)
//!
//! C++ wrapper for using OSQP solver on ARCS
//!
//! @date 2024/9/9
//! @author Juan Padron
//
// Copyright (C) 2011-2024 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.


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


        OSQP_Solver(void):
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
            printf("OSQP_Solver instantiation created \n");
        }


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

        





        ~OSQP_Solver(void)
        {
            printf("OSQP_Solver instantiation destroyed \n");
        }


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





        private:





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

                    if (abs(val)>zero_eps) {                        
                        P_elements.push_back(val);            // Store non-zero value
                        P_rows.push_back(row);       // Store the corresponding row index
                        P_nnz++; //Count number of non-zero elements
                    }
                    
                }
            }
            P_cols[cols]=P_elements.size(); //Needed for filling last indexing element               
              
        }

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

                    if (abs(val)>zero_eps) {                        
                        A_elements.push_back(val);            // Store non-zero value
                        A_rows.push_back(row);       // Store the corresponding row index
                        A_nnz++; //Count number of non-zero elements
                    }
                    
                }
            }
            A_cols[cols]=A_elements.size(); //Needed for filling last indexing element             

            
        }

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


        static void OSQPCscMatrix_deleter(OSQPCscMatrix* mat)
        {
            free(mat);
            printf("CsCMatrix ptr destroyed \n");
        }

        static void OSQPSolver_deleter(OSQPSolver* slv)
        {
            osqp_cleanup(slv);
            printf("Solver ptr destroyed \n");
        }

        static void OSQPSettings_deleter(OSQPSettings* stt)
        {
            free(stt);
            printf("Settings ptr destroyed \n");
        }



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