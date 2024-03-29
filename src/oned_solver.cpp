#include "oned_solver.h"

OneDSolver::OneDSolver():OneDSolver(0,true){
};

// default to 1st dimension and 2nd order if nothing else is specified
OneDSolver::OneDSolver(const size_t& num_nodes, const bool& solver_method):FDSolver::FDSolver(num_nodes,1,solver_method,2, 3*num_nodes){
  std::cout<<"We are constructing oned solver with "<< num_nodes<<std::endl;
};
// including call without specification of verification and debug modes to allow for backward compatibility
OneDSolver::OneDSolver(const size_t& num_nodes, const int& dim, const bool& solver_method,const int& order):FDSolver::FDSolver(num_nodes,dim,solver_method,order, 0, 0,0,3*num_nodes){
  std::cout<<"We are constructing oned solver with "<< num_nodes<<std::endl;
};

OneDSolver::OneDSolver(const size_t& num_nodes, const int& dim, const bool& solver_method, const int& order, const bool& verify, const bool& debug, const bool& USE_PETSC):FDSolver::FDSolver(num_nodes,dim,solver_method,order, verify, debug, USE_PETSC, 3*num_nodes){
  std::cout<<"We are constructing oned solver entirely from libgrvy inputs now "<< num_nodes<<std::endl;
  /*
  if (order == 2){
    scaling_constant = 1.0 / (4 * M_PI * M_PI);
  } elif (order == 4){
    scaling_constant = 1.0 / (12 * M_PI * M_PI);
  } else{
      throw std::invalid_argument("order must be 2 or 4; 2nd or 4th order methods only!");
  }
  */
};

//OneDSolver constructor with all options passed from libgrvy
// ^ TODO

/*
** This function creates the matrices and vectors for the 1D mesh
*/
void OneDSolver::construct_matrix(){
    gsl_spmatrix* X;
    X = gsl_spmatrix_alloc_nzmax(this->num_nodes_no_bndry ,this->num_nodes_no_bndry, this->nnz, GSL_SPMATRIX_COO);
    const double h = 1.0 / (this->num_nodes-1);                 /* grid spacing */
    // Print the grid spacing:
    std::cout << "grid spacing: " << h  << "\n";
    size_t i{0};
    double fi{0.0};
    /* construct the sparse matrix for the finite difference equation */

    // change the method used based on the order of the system
    // might want to make this a function? not sure what is best
    if (this->order == 2){
      /* construct first row */
      gsl_spmatrix_set(X, 0, 0, -2.0);
      gsl_spmatrix_set(X, 0, 1, 1.0);

      /* construct rows [1:n-2] */
      for (i = 1; i < this->num_nodes_no_bndry - 1; ++i)
      {
        gsl_spmatrix_set(X, i, i + 1, 1.0);
        gsl_spmatrix_set(X, i, i, -2.0);
        gsl_spmatrix_set(X, i, i - 1, 1.0);
      }

      /* construct last row */
      gsl_spmatrix_set(X, this->num_nodes_no_bndry - 1, this->num_nodes_no_bndry - 1, -2.0);
      gsl_spmatrix_set(X, this->num_nodes_no_bndry - 1, this->num_nodes_no_bndry - 2, 1.0);

      /* scale by h^2 */
      gsl_spmatrix_scale(X, -1.0 * scaling_constant / (h * h));
    }
    else if (this->order == 4){
      /* construct first row */
      gsl_spmatrix_set(X, 0, 0, -24.0);
      gsl_spmatrix_set(X, 0, 1, 12.0);

      /* construct second row */
      gsl_spmatrix_set(X, 1, 0, 16.0);
      gsl_spmatrix_set(X, 1, 1, -30.0);
      gsl_spmatrix_set(X, 1, 2, 16.0);
      gsl_spmatrix_set(X, 1, 3, -1.0);

      /* construct rows [1:n-2] */
      for (i = 2; i < this->num_nodes_no_bndry - 2; ++i)
      {
          gsl_spmatrix_set(X, i, i + 2, -1.0);
          gsl_spmatrix_set(X, i, i + 1, 16.0);
          gsl_spmatrix_set(X, i, i, -30.0);
          gsl_spmatrix_set(X, i, i - 1, 16.0);
          gsl_spmatrix_set(X, i, i - 2, -1.0);
      }

      /* construct second to last row */
      gsl_spmatrix_set(X, this->num_nodes_no_bndry - 2, this->num_nodes_no_bndry-1, 16.0);
      gsl_spmatrix_set(X, this->num_nodes_no_bndry - 2, this->num_nodes_no_bndry-2, -30.0);
      gsl_spmatrix_set(X, this->num_nodes_no_bndry - 2, this->num_nodes_no_bndry-3, 16.0);
      gsl_spmatrix_set(X, this->num_nodes_no_bndry - 2, this->num_nodes_no_bndry-4, -1.0);

      /* construct last row */
      gsl_spmatrix_set(X, this->num_nodes_no_bndry - 1, this->num_nodes_no_bndry - 1, -24.0);
      gsl_spmatrix_set(X, this->num_nodes_no_bndry - 1, this->num_nodes_no_bndry - 2, 12.0);
      
      /* scale by h^2 */
      gsl_spmatrix_scale(X, -1.0 * scaling_constant / (12 * h * h));
    } 
    else {
      std::cout<< "order must be 2 or 4; 2nd or 4th order methods only!" << std::endl;
      exit(1);
    }


    /* construct right hand side vector */
    for (i = 0; i < this->num_nodes_no_bndry; ++i)
    {
      fi = sin(2.0 * M_PI * (i + 1) * h);
      gsl_vector_set(this->f, i, fi);
    }
    
    /* initialiaze u at zero*/
    gsl_vector_set_zero(this->u);
    const int status = gsl_spmatrix_csr(this->A, X);
    std::cout<<"Successfully copied X COO to A CSR"<<std::endl;
    gsl_spmatrix_free(X);
};
