#include "oned_solver.h"

OneDSolver::OneDSolver():OneDSolver(0,true){
};

OneDSolver::OneDSolver(const size_t& num_nodes, const bool& solver_method):FDSolver::FDSolver(num_nodes,1,solver_method){
  std::cout<<"We are constructing oned solver with "<< num_nodes<<std::endl;
};

/*
** This function creates the matrices and vectors for the 1D mesh
*/
void OneDSolver::construct_matrix(){
    const double h = 1.0 / (this->num_nodes-1);                 /* grid spacing */
    // Print the grid spacing:
    std::cout << "grid spacing: " << h  << "\n";
    size_t i{0};
    double fi{0.0};
    /* construct the sparse matrix for the finite difference equation */

    /* construct first row */
    gsl_spmatrix_set(this->A, 0, 0, -2.0);
    gsl_spmatrix_set(this->A, 0, 1, 1.0);

    /* construct rows [1:n-2] */
    for (i = 1; i < this->num_nodes_no_bndry - 1; ++i)
    {
      gsl_spmatrix_set(A, i, i + 1, 1.0);
      gsl_spmatrix_set(A, i, i, -2.0);
      gsl_spmatrix_set(A, i, i - 1, 1.0);
    }

    /* construct last row */
    gsl_spmatrix_set(A, this->num_nodes_no_bndry - 1, this->num_nodes_no_bndry - 1, -2.0);
    gsl_spmatrix_set(A, this->num_nodes_no_bndry - 1, this->num_nodes_no_bndry - 2, 1.0);

    /* scale by h^2 */
    gsl_spmatrix_scale(A, -1.0 * scaling_constant / (h * h));

    /* construct right hand side vector */
    for (i = 0; i < this->num_nodes_no_bndry; ++i)
    {
      fi = sin(2.0 * M_PI * (i + 1) * h);
      gsl_vector_set(this->f, i, fi);
    }
    
    /* initialiaze u at zero*/
    gsl_vector_set_zero(this->u);
};
