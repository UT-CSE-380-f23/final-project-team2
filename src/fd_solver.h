#pragma once
#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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
    gsl_vector *f{}, *u{};
    const size_t num_nodes;
    const size_t num_nodes_no_bndry;
    const size_t dim;
    const size_t order;
    const bool solver_method; /* Variables used for the Gauss-Seidel iteration */
    const double tol;                   /* solution relative tolerance */
    const int max_iter;                     /* maximum iterations */
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
    FDSolver(const size_t& num_nodes, const size_t& dim, const bool& solver_method, const size_t& order, const int& nnz);
    FDSolver(const size_t& num_nodes, const size_t& dim, const bool& solver_method, const size_t& order, const int& nnz, const double& tol, 
    const int& max_iter);
    ~FDSolver();
    //void (*solver_method)(const size_t& N, gsl_spmatrix& A, gsl_vector& f, const gsl_vector& u)
    void system_solve();
};
