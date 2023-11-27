#include "fd_solver.h"

/*
** Default constructor
*/
FDSolver::FDSolver():FDSolver(0,0,true, 2, 0, 0, 1.0e-6, 1000000){
}; 

/*
** Parametrized constructor
*/
FDSolver::FDSolver(const size_t& num_nodes, const size_t& dim, const bool& solver_method, const size_t& order):FDSolver(num_nodes, dim, solver_method, order, 0, 0, 1.0e-6, 1000000){
};

FDSolver::FDSolver(const size_t& num_nodes, const size_t& dim, const bool& solver_method, const size_t& order, const bool& verify, const bool& debug, const double& tol, 
    const int& max_iter) : num_nodes(num_nodes),dim(dim), solver_method(solver_method), order(order), verify(verify), debug (debug), tol(tol), max_iter(max_iter), num_nodes_no_bndry(num_nodes-2){
        //Allocate memory for matrix and vector
    this->A = gsl_spmatrix_alloc(this->num_nodes_no_bndry ,this->num_nodes_no_bndry); /* triplet format */
    this->f = gsl_vector_alloc(this->num_nodes_no_bndry);        /* right hand side vector */
    this->u = gsl_vector_alloc(this->num_nodes_no_bndry);        /* solution vector */

    // testing!!
    std::cout << "Order of the method that we are using: " << order << std::endl;
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

const double FDSolver::jacobi_element(const gsl_vector* u_prev, const int& j){
    return gsl_vector_get(u_prev, j);
};
const double FDSolver::gauss_sidel_element(const gsl_vector* u_prev, const int& j){
    return gsl_vector_get(this->u, j);
};

// for all *_to_string() methods, might want to define the casted integers 
/*
    Convert the solver_method variable to a string to output L2 error to a file
*/
const std::string FDSolver::solver_method_to_string(){
    int solver_method_int = static_cast<int>(this->solver_method);

    if (this->solver_method_int){
        return "Jacobi";
    } else {
        return "Gauss-Seidel";
    }
}

const float FDSolver::output_L2_norm(){
    // f is the input, for the manufactured solution, input == solution
    gsl_vector_sub(this->f, this->u);

    // open output file
    ofstream outfile;
    outfile.open ("example.txt");
    outfile << "Writing this to a file.\n";
    outfile.close();
    return 0;

    // compute L2 norm of the difference and output to a file
    gsl_blas_dnrm2(this->f) << "../output/" + solver_method_to_string() + "_" + std::to_string(this->dim) + "D_" std::to_string(this->order) + "_order";
    
    // maybe print the L2 norm to the screen in debug mode

}

/*
    Convert finite difference method order to a string
*/

void FDSolver::system_solve(){//(int N_arg, gsl_spmatrix *M, gsl_vector *b, gsl_vector *x, bool jacOrGS){
    /* Some variables */
    this->construct_matrix(this->order); // Construct A, f, u
    // Iteration variables
    std::cout<<"Constructed problem, now solving..."<<std::endl;
    int i{0}, j{0}, k{0};
    double sigma{0.0}, residual{10000.00}, prev_residual{residual+1};
    gsl_vector *u_prev  = gsl_vector_alloc(this->num_nodes_no_bndry);       /* used only for iteration */
    /* initialiaze x at zero*/
    gsl_vector_set_zero(u_prev);
    /* iterate */
    while(this->tol < residual && k < max_iter){
        // check if we need to update this method
        for (i=0; i < this->num_nodes_no_bndry; i++){ // need to optimize this to nnz
            sigma = -gsl_spmatrix_get(this->A, i, i)*(this->*solver_function)(u_prev, i);
            for (j = 0; j < this->num_nodes_no_bndry; j++)
                sigma = sigma + gsl_spmatrix_get(this->A, i, j)*(this->*solver_function)(u_prev, j);
            gsl_vector_set(this->u, i, (gsl_vector_get(this->f, i) - sigma)/gsl_spmatrix_get(this->A, i, i));
        }
        // Update the residual and print it
        gsl_vector_sub(u_prev, this->u);
        residual = gsl_blas_dnrm2(u_prev);
        
        // only print residual in debug mode
        std::cout << "residual: " << residual << "\n";

        // Make a copy of x into x_prev for the following iteration
        gsl_vector_memcpy(u_prev, this->u);

        // Only continue if the resiudal is decreasing
        if (prev_residual <= residual){
            std::printf("Residual did not decrese for %d th iteration, was %f is now %f", k, prev_residual, residual);
            break;
        }

        // Update the previous residual
        prev_residual = residual;
        k++;
    }

    // Compute L2 Norm of the solution if we are in verification mode - perhaps turn this into a method
    if (this->verify){
        this->output_L2_norm();
    }

    // Print u and free up memory in debug mode only
    if (this->debug){
        for (i = 0; i < num_nodes_no_bndry; i++)
        {
            std::cout << "u_" << i << "=" << gsl_vector_get(this->u, i) << "\n";
        }
        std::cout<<"Done, lets free up memory..."<<std::endl;
    }
    // Free space
    gsl_vector_free(u_prev);
};


// #endif // ITERATIVE_METHODS_H