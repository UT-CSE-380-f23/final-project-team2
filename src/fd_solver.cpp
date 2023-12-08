#include "fd_solver.h"

/*
** Default constructor
*/
FDSolver::FDSolver():FDSolver(0,0,true, 2, 0, 0, 0, 0, 1.0e-12, 1000000){
}; 

/*
** Parametrized constructor
*/
// this constructor might be deprecated, think it will return wrong result.
FDSolver::FDSolver(const size_t& num_nodes, const int& dim, const bool& solver_method, const int& order, const size_t& nnz):FDSolver(num_nodes, dim, solver_method, order, 0, 0, 0, nnz){
};

FDSolver::FDSolver(const size_t& num_nodes, const int& dim, const bool& solver_method, const int& order, const bool& verify, const bool& debug, const bool& USE_PETSC, const size_t& nnz):FDSolver(num_nodes, dim, solver_method, order, verify, debug, USE_PETSC, nnz, 1.0e-12, 1000000){
};

FDSolver::FDSolver(const size_t& num_nodes, const int& dim, const bool& solver_method, const int& order, const bool& verify, const bool& debug, const bool& USE_PETSC, const size_t& nnz, const double& tol, 
    const size_t& max_iter) : num_nodes(num_nodes),dim(dim), solver_method(solver_method), order(order), verify(verify), debug(debug), petsc_enabled(USE_PETSC), tol(tol), max_iter(max_iter), num_nodes_no_bndry(num_nodes-2), matrix_length(std::pow(this->num_nodes_no_bndry, this->dim)), nnz(nnz){
        //Allocate memory for matrix and vector
    // testing!!
    std::cout << "Order of the method that we are using: " << order << std::endl;
    this->A = gsl_spmatrix_alloc_nzmax(this->matrix_length ,this->matrix_length, this->nnz, GSL_SPMATRIX_CSR); /* triplet format */
    this->f = gsl_vector_alloc(this->matrix_length);        /* right hand side vector */
    this->f_temp = gsl_vector_alloc(this->matrix_length);   /* temporary vector to compute L2 norm */
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

/*
    Single function that outputs the solution to an hdf5 file
*/
void FDSolver::save_hdf5_data(const char* outfile){

    int dataset_size = int(pow(double(this->num_nodes),double(this->dim)));
    double h = 1.0 / (this->num_nodes-1); // this is not defined in this class...
    NodeData nodes[dataset_size];
    // might want to change nodes to be nodes[][] for the 2D case...

    // Initialize the data structure for the 1D case
    if (this->dim == 1){
        // interior nodes
        for (int i = 1; i < this->num_nodes - 1; ++i) {
            nodes[i].x = i*h;               
            nodes[i].numerical_temp = gsl_vector_get(this->u,i-1);
            nodes[i].analytical_temp = gsl_vector_get(this->f,i-1);
        }
        //boundary nodes
        nodes[0].x = 0.0;
        nodes[this->num_nodes - 1].x = 1.0;
        nodes[0].numerical_temp = 0.0;
        nodes[this->num_nodes - 1].numerical_temp = 0.0;
        nodes[0].analytical_temp = 0.0;
        nodes[this->num_nodes - 1].analytical_temp = 0.0;
        
    } else if(this->dim == 2){
        // interior/non-zero nodes
        for (int i = 1; i < this->num_nodes - 1; i++){
            for (int j=1; j < this->num_nodes - 1; j++){
                // nodal coordaintes
                nodes[i*(this->num_nodes) + j].x = i*h;
                nodes[i*(this->num_nodes) + j].y = j*h;
                // numerical and analytical temp
                nodes[i*(this->num_nodes) + j].numerical_temp = gsl_vector_get(this->u,(i-1)*(this->num_nodes_no_bndry) + j - 1);
                nodes[i*(this->num_nodes) + j].analytical_temp = gsl_vector_get(this->f,(i-1)*(this->num_nodes_no_bndry) + j - 1);
            }
        }
        // exterior nodes, where BC is applied
        for (int j = 0; j < this->num_nodes; j++){
            //data[0][j] = 0; // x = 0

            /*
                Nomenclature for "interior" boundaries

                u_mat = [
                    0, 0, 0, 0, 0
                    0, x, x, x, 0
                    0, x, x, x, 0
                    0, x, x, x, 0
                    0, 0, 0, 0, 0
                ]
                ->
                u_vec = [
                    0, 0, 0, 0, 0, 0, x, x, x, 0, 0, x, x, x, 0
                    *, *, *, *, *, !, ......., !, !, ........
                ]
                * : "exterior" boundary
                ! : "interior" boundary
            */

            // set values on the exterior boundaries (x = 0 and x = 1) to 0
            nodes[j].numerical_temp = 0.0;
            nodes[dataset_size - 1 - j].numerical_temp = 0.0;
            nodes[j].analytical_temp = 0.0;
            nodes[dataset_size - 1 - j].analytical_temp = 0.0;
            nodes[j].x = 0.0;
            nodes[dataset_size - 1 - j].x = 0.0;
            nodes[j].y = 0.0;
            nodes[dataset_size - 1 - j].y = 0.0;

            // set values on the "interior" boundaries (y = 0, y = 1 I believe) to 0.0
            nodes[j*this->num_nodes].numerical_temp = 0.0;
            nodes[(j+1)*this->num_nodes - 1].numerical_temp = 0.0;
            nodes[j*this->num_nodes].analytical_temp = 0.0;
            nodes[(j+1)*this->num_nodes - 1].analytical_temp = 0.0;
            nodes[j*this->num_nodes].x = 0.0;
            nodes[(j+1)*this->num_nodes - 1].x = 0.0;
            nodes[j*this->num_nodes].y = 0.0;
            nodes[(j+1)*this->num_nodes - 1].y = 0.0;

            // this will set the values on the corners to 0.0 twice

        }
    } else {

        std::cout<< "Can only solve a system in 1D or 2D!" << std::endl;
        exit(1);
    }

    // Initialize HDF5
    hid_t file_id = H5Fcreate(outfile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // Create a dataspace for the dataset
    hsize_t dims[1] = {dataset_size};
    hid_t dataspace_id = H5Screate_simple(1, dims, NULL);

    // Create a compound datatype for the node data
    hid_t datatype_id = H5Tcreate(H5T_COMPOUND, sizeof(NodeData));
    H5Tinsert(datatype_id, "x", HOFFSET(NodeData, x), H5T_NATIVE_DOUBLE);
    if (this->dim == 2){
        H5Tinsert(datatype_id, "y", HOFFSET(NodeData, y), H5T_NATIVE_DOUBLE);
    }
    H5Tinsert(datatype_id, "numerical_temperature", HOFFSET(NodeData, numerical_temp), H5T_NATIVE_DOUBLE);
    H5Tinsert(datatype_id, "analytical_temperature", HOFFSET(NodeData, analytical_temp), H5T_NATIVE_DOUBLE);

    // Create the dataset
    hid_t dataset_id = H5Dcreate2(file_id, "/data", datatype_id, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // Write the data to the dataset
    H5Dwrite(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, nodes);

    // Close resources
    H5Dclose(dataset_id);
    H5Tclose(datatype_id);
    H5Sclose(dataspace_id);
    H5Fclose(file_id);
}

// for all *_to_string() methods, might want to define the casted integers 
/*
    Convert the solver_method variable to a string to output L2 error to a file
*/
const std::string FDSolver::solver_method_to_string(){
    int solver_method_int = static_cast<int>(this->solver_method);

    if (solver_method_int){
        return "Jacobi";
    } else {
        return "Gauss-Seidel";
    }
}

void FDSolver::output_L2_norm(){

    // set f_temp equal to f so we can compute the L2 norm without changing f!
    gsl_vector_memcpy(this->f_temp, this->f);
    // subtract computed solution from analytical solution
    gsl_vector_sub(this->f_temp, this->u);

    // open output file
    std::ofstream outfile;
    outfile.open("../output/" + solver_method_to_string() + "_" + std::to_string(this->dim) + "D_" + std::to_string(this->order) + "_order_N_" + std::to_string(this->num_nodes));
    // compute L2 norm of the difference and output to outfile
    double l2norm{gsl_blas_dnrm2(this->f_temp)};
    std::cout<<"L2 norm of the error : "<<l2norm<<std::endl;
    outfile<<l2norm;
    outfile.close();
    std::cout << "Done outputing to file" << std::endl;
}

/*
    Convert finite difference method order to a string
*/

void FDSolver::system_solve(const char* outfile){//(int N_arg, gsl_spmatrix *M, gsl_vector *b, gsl_vector *x, bool jacOrGS){
    /* Some variables */
    this->construct_matrix(); // Construct A, f, u

    // solve the system with or without petsc
    if(this->petsc_enabled){
        this->petsc_solver();
    }else{
        this->iterative_solve(); 
    }
    // Solve Au = f iteratively
    // should probably define the function to use similarly to how we picked the iterative solver to use....check that out.
    /*
    if (this->dim == 1){
        this->save_hdf5_1d_data(outfile);
    } else if(this->dim == 2){
        this->save_hdf5_2d_data(outfile);
    } else{
        std::cout<< "Can only solve 1D or 2D system!" << std::endl;
        exit(1);
    }*/
    this->save_hdf5_data(outfile);
    //this->save_solution(outfile); // save the solution to an HDF5 file.
    
}

PetscErrorCode FDSolver::petsc_solver(){
      
  /*
    Some necessary variables
  */
  PetscErrorCode ierr;  // PETSC error
  Mat P;                // PETSC matrix
  Vec x,b;              // PETSC vectors, x is the solution, b the RHS
  PetscInt n, m;        // sizes -- used for both matrix and vector
  KSP ksp;

  const int nnz = this->A->nz;
  const int *row_idx_arr{this->A->p};
  const int *col_idx_arr{this->A->i};
  const double *data_A{this->A->data};

  /* Only call this *once* in a program */
  PetscInitialize(0,NULL,NULL,NULL);

  /*
    Build the PETS matrix
  */

  // Set it up
  m = this->matrix_length; n = this->matrix_length; // size of the matrix
  ierr = MatCreate(PETSC_COMM_WORLD, &P);CHKERRQ(ierr);
  ierr = MatSetType(P, MATAIJ);CHKERRQ(ierr);
  ierr = MatSetSizes(P, PETSC_DECIDE, PETSC_DECIDE, m, n);CHKERRQ(ierr);
  ierr = MatSetUp(P);CHKERRQ(ierr);

  const int size1 = this->A->size1;
  for (int i = 0; i < size1; i++){
    for (int j = row_idx_arr[i]; j < (row_idx_arr[i+1]); j++){
      ierr = MatSetValue(P, i, col_idx_arr[j], data_A[j], INSERT_VALUES);CHKERRQ(ierr);
    }
  }

  /*
  for (int i = 0; i < A_row->size1; i++){
    for (int j = A_row->p[i]; j < (A_row->p[i+1]); j++){
      ierr = MatSetValue(P, i, A_row->i[j], A_row->data[j], INSERT_VALUES);CHKERRQ(ierr);
    }
  }
  */

  // Populate it
  //PetscReal apq;
  PetscInt p,q;
  
  //for(int s = 0; s < (this->A)->nz; s++){
  //  p = (this->A)->i[s];        // Column
  //  q = (this->A)->p[s];        // Row
  //  apq = (this->A)->data[s];   // Value
  //  ierr = MatSetValue(P, p, q, apq, INSERT_VALUES);CHKERRQ(ierr); 
  //}

  // Assemble the matrix
  ierr = MatAssemblyBegin(P, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(P, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  if (this->debug){
    std::cout << "The right matrix (triplet representation): \n";
    ierr = MatView(P,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }


  /*
    Build the PETS vector
  */

  // Set up
  n = this->matrix_length;
  ierr = VecCreateSeq(PETSC_COMM_SELF,n,&b);CHKERRQ(ierr);
  ierr = VecSet(b,0.0);CHKERRQ(ierr);

  // Populate the RHS
  PetscReal bp;
  for (int s = 0; s < n; s++){
    p = s;
    bp = gsl_vector_get(this->f, s); 
    ierr = VecSetValue(b, p, bp, INSERT_VALUES);CHKERRQ(ierr);
  }

  // Assembly
  ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(b);CHKERRQ(ierr);

  // View the RHS (optional)
  if (this->debug){
    std::cout << "The right hand side: \n";
    ierr = VecView(b, 0);CHKERRQ(ierr);
  }

  /*
  
    Now we solve the system using GMRES
  
  */
  PetscReal rtol = 1e-10; // default is 1e-5
  PetscReal atol=1e-50; // default
  PetscReal dtol=1e5;   // default
  PetscReal maxits=1e4; // default
 
  /* The solution vector */
  ierr = VecCreateSeq(PETSC_COMM_SELF,n,&x);CHKERRQ(ierr);
  ierr = VecSet(x,0.0);CHKERRQ(ierr);
  /* solve linear system */
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,P,P);CHKERRQ(ierr);
  ierr = KSPSetType(ksp, KSPGMRES);
  // set the tolerance
  ierr = KSPSetTolerances(ksp, rtol, atol, dtol, maxits); CHKERRQ(ierr);
  ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);

  /* some post-processing */
  KSPConvergedReason reason;
  PetscInt nits;
  PetscReal rnorm;
  KSPGetConvergedReason(ksp, &reason);
  KSPGetIterationNumber(ksp, &nits);
  KSPGetResidualNorm(ksp, &rnorm);

  if(reason > 0){
    std::cout << "Converged after " << nits << " iterations, with a residual norm of " << rnorm << " \n";
  }else{
    std::cout << "Did not converge after " << nits << " iterations \n";
    std::cout << "residual = " << rnorm << "\n";
  };

  /* compute the L2 error */
  PetscScalar a = -1.0;
  ierr = VecAXPY(b, a, x);CHKERRQ(ierr);
  PetscReal   l2_error = 0.0;
  ierr = VecNorm(b, NORM_2,&l2_error);CHKERRQ(ierr);

  std::cout << "L_2 error: " << l2_error << "\n";

    /*

    Populate the GSL solution

    */

    PetscInt ix[1] = {0};
    PetscScalar y[1] = {0.0};

    for (int s = 0; s < n; s++){
        ix[0] = s;
        ierr = VecGetValues(x, 1, ix, y);CHKERRQ(ierr);
        gsl_vector_set(this->u, s, y[0]);
    }

  /* 
    Clean up 
  */
  
  // KSP
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  // Matrix
  ierr = MatDestroy(&P);CHKERRQ(ierr);
  // RHS
  ierr = VecDestroy(&b);CHKERRQ(ierr);
  // Solution vector
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  /* Only call this *once* in a program -- finalize petsc */
  ierr = PetscFinalize();

  return ierr;
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

    // Compute L2 Norm of the solution if we are in verification mode - perhaps turn this into a method
    if (this->verify){
        std::cout << "saving L2 norm" << std::endl;
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
