#include "twod_solver.h"

TwoDSolver::TwoDSolver():TwoDSolver(0,true){
};

// default to 1st dimension and 2nd order if nothing else is specified
TwoDSolver::TwoDSolver(const size_t& num_nodes, const bool& solver_method):FDSolver::FDSolver(num_nodes,2,solver_method,2){
    std::cout<<"We are constructing oned solver with "<< num_nodes<<std::endl;
};

TwoDSolver::TwoDSolver(const size_t& num_nodes, const size_t& dim, const bool& solver_method, const size_t& order):FDSolver::FDSolver(num_nodes,dim,solver_method,order){
    std::cout<<"We are constructing oned solver entirely from libgrvy inputs now "<< num_nodes<<std::endl;
};

TwoDSolver::TwoDSolver(const size_t& num_nodes, const size_t& dim, const bool& solver_method, const size_t& order,const bool& verify, const bool& debug):FDSolver::FDSolver(num_nodes,dim,solver_method,order,verify,debug){
    std::cout<<"We are constructing oned solver entirely from libgrvy inputs now "<< num_nodes<<std::endl;
};

//TwoDSolver constructor with all options passed from libgrvy
// ^ TODO

/*
** This function creates the matrices and vectors for the 1D mesh
*/
void TwoDSolver::construct_matrix(){
    const double h = 1.0 / (this->num_nodes-1);                 /* grid spacing */
    // Print the grid spacing:
    std::cout << "grid spacing: " << h  << "\n";
    size_t i{0};
    double fi{0.0};
    /* construct the sparse matrix for the finite difference equation */
    double xi, yj, f_val;
    int centre{0}, left{0}, right{0}, up{0}, down{0};

    for (int i = 0; i < this->num_nodes_no_bndry; i++){
      for (int j = 0; j < this->num_nodes_no_bndry; j++){
        // Set the matrix
        centre = i*this->num_nodes_no_bndry+j;
        left = (i-1)*this->num_nodes_no_bndry+j;
        right = (i+1)*this->num_nodes_no_bndry+j;
        up = i*this->num_nodes_no_bndry+j+1;
        down = i*this->num_nodes_no_bndry+j-1;
        gsl_spmatrix_set(A, centre, centre, -4.0);
        if(i != 0)
            gsl_spmatrix_set(A, centre, left, 1.0);
        if(i != this->num_nodes_no_bndry-1)
            gsl_spmatrix_set(A, centre, right, 1.0);
        if(j != 0)
            gsl_spmatrix_set(A, centre, down, 1.0);
        if(j != this->num_nodes_no_bndry-1)
            gsl_spmatrix_set(A, centre, up, 1.0);
        // Set the RHS
        xi = (i + 1) * h;
        yj = (j + 1) * h;
        f_val = sin(2.0 * M_PI * xi)*sin(2.0 * M_PI * yj);
        gsl_vector_set(f, i*this->num_nodes_no_bndry+j, f_val);
      }
    }
    gsl_spmatrix_scale(A, (-1.0)*scaling_constant_2d/(h*h));
    /* initialiaze u at zero*/
    gsl_vector_set_zero(this->u);
};
