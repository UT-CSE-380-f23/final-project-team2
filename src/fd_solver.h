#pragma once
#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_splinalg.h>
#include <gsl/gsl_spblas.h>

/*

FDM Solver Abstract Class 

*/
class FDSolver{
  protected:
    gsl_spmatrix *A{};
    gsl_vector *f{}, *u{}; //Class variables for solving Au = f
    const size_t num_nodes;
    const size_t num_nodes_no_bndry;
    const int dim;
    const int order;
    const bool solver_method; /* Variables used for the Gauss-Seidel iteration */
    const bool verify;
    const bool debug;
    const double tol;                   /* solution relative tolerance */
    const size_t max_iter;                     /* maximum iterations */
    const int nnz;
    const size_t matrix_length;
    typedef const double (FDSolver::*fn_element)(const gsl_vector*, const int&); // This is pointer to function for  defining element
    fn_element solver_function;
    void iterative_solve();
    virtual void construct_matrix() = 0; // pure virtual function defined in derived class
    inline const double jacobi_element(const gsl_vector *u_prev, const int &j);
    inline const double gauss_sidel_element(const gsl_vector *u_prev, const int &j);
  public:
    FDSolver();
    FDSolver(const size_t& num_nodes, const int& dim, const bool& solver_method, const int& order, const size_t& nnz);
    FDSolver(const size_t& num_nodes, const int& dim, const bool& solver_method, const int& order, const bool& verify, const bool& debug, const size_t& nnz);
    FDSolver(const size_t& num_nodes, const int& dim, const bool& solver_method, const int& order, const bool& verify, const bool& debug, const size_t& nnz, const double& tol, const size_t& max_iter);
    ~FDSolver();
    //void (*solver_method)(const size_t& N, gsl_spmatrix& A, gsl_vector& f, const gsl_vector& u)
    void system_solve();
    //adding order as a parameter to the construct matrix function
    //virtual void construct_matrix(const int& order)=0; // pure virtual function defined in derived class
    //virtual void construct_matrix()=0; // pure virtual function defined in derived class
    //const double jacobi_element(const gsl_vector* u_prev, const int& j);
    //const double gauss_sidel_element(const gsl_vector* u_prev, const int& j);
    const std::string solver_method_to_string();
    void output_L2_norm();

};
