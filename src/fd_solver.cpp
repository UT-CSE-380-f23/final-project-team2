#include "fd_solver.h"

/*
** Default constructor
*/
FDSolver::FDSolver():FDSolver(0,0,true, 2, 0, 1.0e-12, 1000000){
}; 

/*
** Parametrized constructor
*/
FDSolver::FDSolver(const size_t& num_nodes, const size_t& dim, const bool& solver_method, const size_t& order, const int& nnz):FDSolver(num_nodes, dim, solver_method, order, nnz, 1.0e-12, 1000000){
};

FDSolver::FDSolver(const size_t& num_nodes, const size_t& dim, const bool& solver_method, const size_t& order, const int& nnz, const double& tol, 
    const int& max_iter) : num_nodes(num_nodes),dim(dim), solver_method(solver_method), order(order), tol(tol), max_iter(max_iter), num_nodes_no_bndry(num_nodes-2), matrix_length(std::pow(this->num_nodes_no_bndry, this->dim)), nnz(nnz){
        //Allocate memory for matrix and vector
    // testing!!
    std::cout << "Order of the method that we are using: " << order << std::endl;
    this->A = gsl_spmatrix_alloc_nzmax(this->matrix_length ,this->matrix_length, this->nnz, GSL_SPMATRIX_CSR); /* triplet format */
    this->f = gsl_vector_alloc(this->matrix_length);        /* right hand side vector */
    this->u = gsl_vector_alloc(this->matrix_length);        /* solution vector */

    if(this->solver_method) // Deciding which solver element to use
        this->solver_function=&FDSolver::jacobi_element;
    else
        this->solver_function=&FDSolver::gauss_sidel_element;
};

FDSolver::~FDSolver(){ // Free memory for A, u, f
    gsl_vector_free(this->f);
    gsl_vector_free(this->u);
    gsl_spmatrix_free(this->A);
};

inline const double FDSolver::jacobi_element(const gsl_vector* u_prev, const int& j){
    return gsl_vector_get(u_prev, j);
};
inline const double FDSolver::gauss_sidel_element(const gsl_vector* u_prev, const int& j){
    return gsl_vector_get(this->u, j);
};

void FDSolver::system_solve(){//(int N_arg, gsl_spmatrix *M, gsl_vector *b, gsl_vector *x, bool jacOrGS){
    /* Some variables */
    this->construct_matrix(); // Construct A, f, u
    this->iterative_solve(); // Solve Au = f iteratively
}
void FDSolver::iterative_solve()
{
    std::cout << "Constructed problem, now solving..." << std::endl;
    // Iteration variables
    int row_idx{0}, i{0}, j{0}, k{0};
    double sigma{0}, residual{10000.00}, prev_residual{residual + 1}, Aii{0.0};
    gsl_vector *u_prev = gsl_vector_alloc(this->matrix_length); /* used only for iteration */
    /* initialiaze x at zero*/
    gsl_vector_set_zero(u_prev);
    /* iterate */
    const int nnz = this->A->nz;
    const int *row_idx_arr{this->A->p};
    const int *col_idx_arr{this->A->i};
    const double *data_A{this->A->data};
    printf("matrix is '%s' format.\n", gsl_spmatrix_type(this->A));
    while (this->tol < residual && k < max_iter)
    {
        for (row_idx = 0; row_idx < this->matrix_length; row_idx++)
        {
            Aii = gsl_spmatrix_get(this->A, row_idx, row_idx);
            sigma = -Aii * (this->*solver_function)(u_prev, row_idx);
            for (j = row_idx_arr[row_idx]; j < row_idx_arr[row_idx + 1]; j++)
                sigma += data_A[j] * (this->*solver_function)(u_prev, col_idx_arr[j]);
            gsl_vector_set(this->u, row_idx, (gsl_vector_get(this->f, row_idx) - sigma) / Aii);
        }
        // Update the residual and print it
        gsl_vector_sub(u_prev, this->u);
        residual = gsl_blas_dnrm2(u_prev);
        // std::cout << "residual: " << residual << "\n";

        // Make a copy of x into x_prev for the following iteration
        gsl_vector_memcpy(u_prev, this->u);

        // Only continue if the resiudal is decreasing
        if (prev_residual <= residual)
        {
            std::printf("Residual did not decrease for %d th iteration, was %f is now %f", k, prev_residual, residual);
            break;
        }

        // Update the previous residual
        prev_residual = residual;
        k++;
    }

    // Print u
    // for (i = 0; i < this->matrix_length; i++)
    // {
    //     std::cout << "u_" << i << "=" << gsl_vector_get(this->u, i) << "\n";
    // }
    gsl_vector_sub(f, this->u);
    std::cout<<"\nFinal L2 error : "<< gsl_blas_dnrm2(f)<<std::endl;
    printf("Done, in %d iterations, lets free up memory...\n", k);
    std::cout<<"Final residual is "<<residual<<std::endl;
    // Free space
    gsl_vector_free(u_prev);
};

// #endif // ITERATIVE_METHODS_H
