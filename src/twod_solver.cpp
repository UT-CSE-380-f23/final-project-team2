#include "twod_solver.h"

TwoDSolver::TwoDSolver():TwoDSolver(0,true){
};

// default to 1st dimension and 2nd order if nothing else is specified
TwoDSolver::TwoDSolver(const size_t& num_nodes, const bool& solver_method):FDSolver::FDSolver(num_nodes,2,solver_method,2, 9*num_nodes*num_nodes){
    std::cout<<"We are constructing oned solver with "<< num_nodes<<std::endl;
};

TwoDSolver::TwoDSolver(const size_t& num_nodes, const int& dim, const bool& solver_method, const int& order):FDSolver::FDSolver(num_nodes,dim,solver_method,order, 9*num_nodes*num_nodes){
    std::cout<<"We are constructing oned solver entirely from libgrvy inputs now "<< num_nodes<<std::endl;
};

TwoDSolver::TwoDSolver(const size_t& num_nodes, const int& dim, const bool& solver_method, const int& order,const bool& verify, const bool& debug):FDSolver::FDSolver(num_nodes,dim,solver_method,order,verify,debug, 9*num_nodes*num_nodes){
    std::cout<<"We are constructing oned solver entirely from libgrvy inputs now "<< num_nodes<<std::endl;
};

//TwoDSolver constructor with all options passed from libgrvy
// ^ TODO

/*
** This function creates the matrices and vectors for the 1D mesh
*/
void TwoDSolver::construct_matrix(){
    gsl_spmatrix* X;
    X = gsl_spmatrix_alloc_nzmax(this->matrix_length ,this->matrix_length, this->nnz, GSL_SPMATRIX_COO);
    const double h = 1.0 / (this->num_nodes-1);                 /* grid spacing */
    // Print the grid spacing:
    std::cout << "grid spacing: " << h  << "\n";
    /* construct the sparse matrix for the finite difference equation */
    double xi, yj, f_val;
    int centre{0}, left{0}, right{0}, up{0}, down{0}, up2{0}, down2{0}, left2{0}, right2{0};
    if (this->order == 2){
        for (int i = 0; i < this->num_nodes_no_bndry; i++){
        for (int j = 0; j < this->num_nodes_no_bndry; j++){
            // Set the matrix
            centre = i*this->num_nodes_no_bndry+j;
            left = (i-1)*this->num_nodes_no_bndry+j;
            right = (i+1)*this->num_nodes_no_bndry+j;
            up = i*this->num_nodes_no_bndry+j+1;
            down = i*this->num_nodes_no_bndry+j-1;
            gsl_spmatrix_set(X, centre, centre, -4.0);
            if(i != 0)
                gsl_spmatrix_set(X, centre, left, 1.0);
            if(i != this->num_nodes_no_bndry-1)
                gsl_spmatrix_set(X, centre, right, 1.0);
            if(j != 0)
                gsl_spmatrix_set(X, centre, down, 1.0);
            if(j != this->num_nodes_no_bndry-1)
                gsl_spmatrix_set(X, centre, up, 1.0);
            // Set the RHS
            xi = (i + 1) * h;
            yj = (j + 1) * h;
            f_val = sin(2.0 * M_PI * xi)*sin(2.0 * M_PI * yj);
            gsl_vector_set(f, i*this->num_nodes_no_bndry+j, f_val);
        }
        }
        gsl_spmatrix_scale(X, (-1.0)*scaling_constant_2d/(h*h));
    }
    else if (this->order == 4)
    {
        /*
        Fill the internal nodes
        */
        for (int i = 1; i < this->num_nodes_no_bndry-1; i++){
            for (int j = 1; j < this->num_nodes_no_bndry-1; j++){
            // Set the matrix
            centre = i*this->num_nodes_no_bndry+j;
            left = (i-1)*this->num_nodes_no_bndry+j;
            right = (i+1)*this->num_nodes_no_bndry+j;
            up = i*this->num_nodes_no_bndry+j+1;
            down = i*this->num_nodes_no_bndry+j-1;
            left2 = (i-2)*this->num_nodes_no_bndry+j;
            right2 = (i+2)*this->num_nodes_no_bndry+j;
            up2 = i*this->num_nodes_no_bndry+j+2;
            down2 = i*this->num_nodes_no_bndry+j-2;
            gsl_spmatrix_set(X, centre, centre, -60.0);
            gsl_spmatrix_set(X, centre, left, 16.0);
            gsl_spmatrix_set(X, centre, right, 16.0);

            gsl_spmatrix_set(X, centre, down, 16.0);
            gsl_spmatrix_set(X, centre, up, 16.0);

            // End points
            if(i != 1)  {gsl_spmatrix_set(X, centre, left2, -1.0);}
            if(i != this->num_nodes_no_bndry-2){gsl_spmatrix_set(X, centre, right2, -1.0);}

            if(j != 1)  {gsl_spmatrix_set(X, centre, down2, -1.0);}
            if(j != this->num_nodes_no_bndry-2){gsl_spmatrix_set(X, centre, up2, -1.0);}

            // Set the RHS
            xi = (i + 1) * h;
            yj = (j + 1) * h;
            f_val = sin(2.0 * M_PI * xi)*sin(2.0 * M_PI * yj);
            gsl_vector_set(f, centre, f_val);

            }
        }

        /*
        indices to be used
        */
        int i = 0;
        int j = 0;


        /*
        i = 0 case
        */

        i = 0;
        for (int j = 1; j < this->num_nodes_no_bndry-1; j++){
            centre = i*this->num_nodes_no_bndry+j;
            left = (i-1)*this->num_nodes_no_bndry+j;
            right = (i+1)*this->num_nodes_no_bndry+j;
            up = i*this->num_nodes_no_bndry+j+1;
            down = i*this->num_nodes_no_bndry+j-1;
            left2 = (i-2)*this->num_nodes_no_bndry+j;
            right2 = (i+2)*this->num_nodes_no_bndry+j;
            up2 = i*this->num_nodes_no_bndry+j+2;
            down2 = i*this->num_nodes_no_bndry+j-2;
            // Set the matrix
            gsl_spmatrix_set(X, centre, centre, -54.0);
            gsl_spmatrix_set(X, centre, right, 12.0);
            gsl_spmatrix_set(X, centre, down, 16.0);
            gsl_spmatrix_set(X, centre, up, 16.0);

            // End points
            // if(i != 1)  {gsl_spmatrix_set(A, centre, left2, -1.0);}
            // if(i != this->num_nodes_no_bndry-2){gsl_spmatrix_set(A, centre, right2, -1.0);}

            if(j != 1)  {gsl_spmatrix_set(X, centre, down2, -1.0);}
            if(j != this->num_nodes_no_bndry-2){gsl_spmatrix_set(X, centre, up2, -1.0);}

            // Set the RHS
            xi = (i + 1) * h;
            yj = (j + 1) * h;
            f_val = sin(2.0 * M_PI * xi)*sin(2.0 * M_PI * yj);
            gsl_vector_set(f, centre, f_val);

        }

        /*
        i = this->num_nodes_no_bndry-1 case
        */

        i = this->num_nodes_no_bndry-1;
        for (int j = 1; j < this->num_nodes_no_bndry-1; j++){
            centre = i*this->num_nodes_no_bndry+j;
            left = (i-1)*this->num_nodes_no_bndry+j;
            right = (i+1)*this->num_nodes_no_bndry+j;
            up = i*this->num_nodes_no_bndry+j+1;
            down = i*this->num_nodes_no_bndry+j-1;
            left2 = (i-2)*this->num_nodes_no_bndry+j;
            right2 = (i+2)*this->num_nodes_no_bndry+j;
            up2 = i*this->num_nodes_no_bndry+j+2;
            down2 = i*this->num_nodes_no_bndry+j-2;
            // Set the matrix
            gsl_spmatrix_set(X, centre, centre, -54.0);

            gsl_spmatrix_set(X, centre, left, 12.0);
            // gsl_spmatrix_set(A, centre, right, 12.0);

            gsl_spmatrix_set(X, centre, down, 16.0);
            gsl_spmatrix_set(X, centre, up, 16.0);

            // End points
            // if(i != 1)  {gsl_spmatrix_set(A, centre, left2, -1.0);}
            // if(i != this->num_nodes_no_bndry-2){gsl_spmatrix_set(A, centre, right2, -1.0);}

            if(j != 1)  {gsl_spmatrix_set(X, centre, down2, -1.0);}
            if(j != this->num_nodes_no_bndry-2){gsl_spmatrix_set(X, centre, up2, -1.0);}

            // Set the RHS
            xi = (i + 1) * h;
            yj = (j + 1) * h;
            f_val = sin(2.0 * M_PI * xi)*sin(2.0 * M_PI * yj);
            gsl_vector_set(f, centre, f_val);
        }

        /*
        j = 0 case
        */
        j = 0;
        for (int i = 1; i < this->num_nodes_no_bndry-1; i++){
            centre = i*this->num_nodes_no_bndry+j;
            left = (i-1)*this->num_nodes_no_bndry+j;
            right = (i+1)*this->num_nodes_no_bndry+j;
            up = i*this->num_nodes_no_bndry+j+1;
            down = i*this->num_nodes_no_bndry+j-1;
            left2 = (i-2)*this->num_nodes_no_bndry+j;
            right2 = (i+2)*this->num_nodes_no_bndry+j;
            up2 = i*this->num_nodes_no_bndry+j+2;
            down2 = i*this->num_nodes_no_bndry+j-2;

            // Set the matrix
            gsl_spmatrix_set(X, centre, centre, -54.0);

            gsl_spmatrix_set(X, centre, left, 16.0);
            gsl_spmatrix_set(X, centre, right, 16.0);

            //gsl_spmatrix_set(A, centre, down, 12.0);
            gsl_spmatrix_set(X, centre, up, 12.0);

            // End points
            if(i != 1)  {gsl_spmatrix_set(X, centre, left2, -1.0);}
            if(i != this->num_nodes_no_bndry-2){gsl_spmatrix_set(X, centre, right2, -1.0);}

            //if(j != 1)  {gsl_spmatrix_set(A, centre, down2, -1.0);}
            //if(j != this->num_nodes_no_bndry-2){gsl_spmatrix_set(A, centre, up2, -1.0);}

            // Set the RHS
            xi = (i + 1) * h;
            yj = (j + 1) * h;
            f_val = sin(2.0 * M_PI * xi)*sin(2.0 * M_PI * yj);
            gsl_vector_set(f, centre, f_val);

        }

        /*
        j = this->num_nodes_no_bndry-1 case
        */
        j = this->num_nodes_no_bndry-1;
        for (int i = 1; i < this->num_nodes_no_bndry-1; i++){
            centre = i*this->num_nodes_no_bndry+j;
            left = (i-1)*this->num_nodes_no_bndry+j;
            right = (i+1)*this->num_nodes_no_bndry+j;
            up = i*this->num_nodes_no_bndry+j+1;
            down = i*this->num_nodes_no_bndry+j-1;
            left2 = (i-2)*this->num_nodes_no_bndry+j;
            right2 = (i+2)*this->num_nodes_no_bndry+j;
            up2 = i*this->num_nodes_no_bndry+j+2;
            down2 = i*this->num_nodes_no_bndry+j-2;
            // Set the matrix
            gsl_spmatrix_set(X, centre, centre, -54.0);

            gsl_spmatrix_set(X, centre, left, 16.0);
            gsl_spmatrix_set(X, centre, right, 16.0);

            gsl_spmatrix_set(X, centre, down, 12.0);
            // gsl_spmatrix_set(A, centre, up, 12.0);

            // End points
            if(i != 1)  {gsl_spmatrix_set(X, centre, left2, -1.0);}
            if(i != this->num_nodes_no_bndry-2){gsl_spmatrix_set(X, centre, right2, -1.0);}

            //if(j != 1)  {gsl_spmatrix_set(A, centre, down2, -1.0);}
            //if(j != this->num_nodes_no_bndry-2){gsl_spmatrix_set(A, centre, up2, -1.0);}

            // Set the RHS
            xi = (i + 1) * h;
            yj = (j + 1) * h;
            f_val = sin(2.0 * M_PI * xi)*sin(2.0 * M_PI * yj);
            gsl_vector_set(f, centre, f_val);

        }

        /*
        Handle the four edge cases #1
        */
        i = 0;
        j = 0;
        centre = i*this->num_nodes_no_bndry+j;
        left = (i-1)*this->num_nodes_no_bndry+j;
        right = (i+1)*this->num_nodes_no_bndry+j;
        up = i*this->num_nodes_no_bndry+j+1;
        down = i*this->num_nodes_no_bndry+j-1;
        left2 = (i-2)*this->num_nodes_no_bndry+j;
        right2 = (i+2)*this->num_nodes_no_bndry+j;
        up2 = i*this->num_nodes_no_bndry+j+2;
        down2 = i*this->num_nodes_no_bndry+j-2;

        gsl_spmatrix_set(X, centre, centre, -48.0);

        //gsl_spmatrix_set(A, centre, left, 12.0);
        gsl_spmatrix_set(X, centre, right, 12.0);

        //gsl_spmatrix_set(A, centre, down, 12.0);
        gsl_spmatrix_set(X, centre, up, 12.0);

        // Set the RHS
        xi = (i + 1) * h;
        yj = (j + 1) * h;
        f_val = sin(2.0 * M_PI * xi)*sin(2.0 * M_PI * yj);
        gsl_vector_set(f, centre, f_val);

        // /*
        // Handle the four edge cases #2
        // */
        i = 0;
        j = this->num_nodes_no_bndry-1;
        centre = i*this->num_nodes_no_bndry+j;
        left = (i-1)*this->num_nodes_no_bndry+j;
        right = (i+1)*this->num_nodes_no_bndry+j;
        up = i*this->num_nodes_no_bndry+j+1;
        down = i*this->num_nodes_no_bndry+j-1;
        left2 = (i-2)*this->num_nodes_no_bndry+j;
        right2 = (i+2)*this->num_nodes_no_bndry+j;
        up2 = i*this->num_nodes_no_bndry+j+2;
        down2 = i*this->num_nodes_no_bndry+j-2;

        gsl_spmatrix_set(X, centre, centre, -48.0);

        //gsl_spmatrix_set(X, centre, left, 12.0);
        gsl_spmatrix_set(X, centre, right, 12.0);

        gsl_spmatrix_set(X, centre, down, 12.0);
        //gsl_spmatrix_set(X, centre, up, 12.0);

        // Set the RHS
        xi = (i + 1) * h;
        yj = (j + 1) * h;
        f_val = sin(2.0 * M_PI * xi)*sin(2.0 * M_PI * yj);
        gsl_vector_set(f, centre, f_val);

        /*
        Handle the four edge cases #3
        */
        i = this->num_nodes_no_bndry-1;
        j = 0;

        centre = i*this->num_nodes_no_bndry+j;
        left = (i-1)*this->num_nodes_no_bndry+j;
        right = (i+1)*this->num_nodes_no_bndry+j;
        up = i*this->num_nodes_no_bndry+j+1;
        down = i*this->num_nodes_no_bndry+j-1;
        left2 = (i-2)*this->num_nodes_no_bndry+j;
        right2 = (i+2)*this->num_nodes_no_bndry+j;
        up2 = i*this->num_nodes_no_bndry+j+2;
        down2 = i*this->num_nodes_no_bndry+j-2;

        gsl_spmatrix_set(X, centre, centre, -48.0);

        gsl_spmatrix_set(X, centre, left, 12.0);
        //gsl_spmatrix_set(X, centre, right, 12.0);

        //gsl_spmatrix_set(X, centre, down, 12.0);
        gsl_spmatrix_set(X, centre, up, 12.0);

        // Set the RHS
        xi = (i + 1) * h;
        yj = (j + 1) * h;
        f_val = sin(2.0 * M_PI * xi)*sin(2.0 * M_PI * yj);
        gsl_vector_set(f, centre, f_val);

        // /*
        // Handle the four edge cases #4
        // */
        i = this->num_nodes_no_bndry-1;
        j = this->num_nodes_no_bndry-1;
        centre = i*this->num_nodes_no_bndry+j;
        left = (i-1)*this->num_nodes_no_bndry+j;
        right = (i+1)*this->num_nodes_no_bndry+j;
        up = i*this->num_nodes_no_bndry+j+1;
        down = i*this->num_nodes_no_bndry+j-1;
        left2 = (i-2)*this->num_nodes_no_bndry+j;
        right2 = (i+2)*this->num_nodes_no_bndry+j;
        up2 = i*this->num_nodes_no_bndry+j+2;
        down2 = i*this->num_nodes_no_bndry+j-2;
        gsl_spmatrix_set(X, centre, centre, -48.0);

        gsl_spmatrix_set(X, centre, left, 12.0);
        //gsl_spmatrix_set(X, centre, right, 12.0);

        gsl_spmatrix_set(X, centre, down, 12.0);
        //gsl_spmatrix_set(X, centre, up, 12.0);

        // Set the RHS
        xi = (i + 1) * h;
        yj = (j + 1) * h;
        f_val = sin(2.0 * M_PI * xi)*sin(2.0 * M_PI * yj);
        gsl_vector_set(f, centre, f_val);

        // Scale the matrix
        gsl_spmatrix_scale(X, (-1.0)*scaling_constant_2d/(12*h*h)); /* Remember to comment/uncoment this */
    }
    else {
      std::cout<< "order must be 2 or 4; 2nd or 4th order methods only!" << std::endl;
      exit(1);
    }
    /* initialiaze u at zero*/
    gsl_vector_set_zero(this->u);
    const int status = gsl_spmatrix_csr(this->A, X);
    if(status)
        std::cout<<"Successfully copied X COO to A CSR"<<std::endl;
    gsl_spmatrix_free(X);
};
