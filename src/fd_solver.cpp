#include "fd_solver.h"

/*
** Default constructor
*/
FDSolver::FDSolver():FDSolver(0,0,true, 2, 0, 0, 0, 1.0e-12, 1000000){
}; 

/*
** Parametrized constructor
*/
FDSolver::FDSolver(const size_t& num_nodes, const int& dim, const bool& solver_method, const int& order, const size_t& nnz):FDSolver(num_nodes, dim, solver_method, order, nnz, 1.0e-12, 1000000){
};

FDSolver::FDSolver(const size_t& num_nodes, const int& dim, const bool& solver_method, const int& order, const bool& verify, const bool& debug, const size_t& nnz):FDSolver(num_nodes, dim, solver_method, order, verify, debug, nnz, 1.0e-12, 1000000){
};

FDSolver::FDSolver(const size_t& num_nodes, const int& dim, const bool& solver_method, const int& order, const bool& verify, const bool& debug, const size_t& nnz, const double& tol, 
    const size_t& max_iter) : num_nodes(num_nodes),dim(dim), solver_method(solver_method), order(order), verify(verify), debug(debug), tol(tol), max_iter(max_iter), num_nodes_no_bndry(num_nodes-2), matrix_length(std::pow(this->num_nodes_no_bndry, this->dim)), nnz(nnz){
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

// properly initializes the objects for the 
void FDSolver::save_hdf5_1d_data(const char* outfile){
    //#define DATASETNAME "heatEqnOutput"
    #define DATASETNAME "heatEqnOutput"
    int RANK=this->dim;
    
    //#define RANK        this->dim
    hid_t   file, dataset;       /* file and dataset handles */
    hid_t   datatype, dataspace; /* handles */
    herr_t  status;
    hsize_t dimsf[1];
    double     data[this->num_nodes]; /* data to write */

    dimsf[0]  = this->num_nodes;

    // interior/non-zero nodes
    for (int j = 1; j < this->num_nodes - 1; j++){
        data[j] = gsl_vector_get(this->u,j-1);
    }

    // boundary conditions for now (ew!)
    data[0] = 0;
    data[this->num_nodes - 1] = 0;

    /*
     * Data  and output buffer initialization.
     */
    /*
    for (j = 0; j < this->grvy_parser.N; j++)
        for (i = 0; i < NY; i++)
            data[j][i] = i + j;
    */
    //data[0][0] += 1e-6;
    /*
     * 0+eps 1 2 3 4 5
     * 1 2 3 4 5 6
     * 2 3 4 5 6 7
     * 3 4 5 6 7 8
     * 4 5 6 7 8 9
     */
    file = H5Fcreate(outfile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /*
     * Describe the size of the array and create the data space for fixed
     * size dataset.
     */

    dataspace = H5Screate_simple(RANK, dimsf, NULL);
    //dataspace = H5Screate_scalar(RANK, dimsf, NULL);

    /*
     * Define datatype for the data in the file.
     * We will store little endian INT numbers.
     */
    datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
    status   = H5Tset_order(datatype, H5T_ORDER_LE);

    /*
     * Create a new dataset within the file using defined dataspace and
     * datatype and default dataset creation properties.
     */
    dataset = H5Dcreate2(file, DATASETNAME, datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /*
     * Write the data to the dataset using default transfer properties.
     */
    //std::cout << "checking value of data " << data[0] << std::endl;
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    /*
     * Close/release resources.
     */
    H5Sclose(dataspace);
    H5Tclose(datatype);
    H5Dclose(dataset);
    H5Fclose(file);

}

void FDSolver::save_hdf5_2d_data(const char* outfile){
    // variables for HDF5
    #define DATASETNAME "heatEqnOutput"
    #define RANK        this->dim
    hid_t   file, dataset;       /* file and dataset handles */
    hid_t   datatype, dataspace; /* handles */
    herr_t  status;
    hsize_t dimsf[2];
    double     data[this->num_nodes][this->num_nodes]; /* data to write */
    
    dimsf[0]  = this->num_nodes;
    dimsf[1]  = this->num_nodes;

    // interior/non-zero nodes
    for (int i = 1; i < this->num_nodes - 1; i++){
        for (int j=1; j < this->num_nodes - 1; j++){
            data[i][j] = gsl_vector_get(this->u,(i-1)*(this->num_nodes_no_bndry) + j - 1);
        }
    }
    // exterior nodes, where BC is applied
    for (int j = 0; j < this->num_nodes; j++){
        data[0][j] = 0; // x = 0
        data[this->num_nodes - 1][j] = 0; // x = 1 (x \in [0,1] for us)

        data[j][0] = 0; // y = 0
        data[j][this->num_nodes - 1] = 0; // y = 1

    }

    // generate the data we are going to output to the file
    //double data = (this->*hdf5_function)();
    /*
     * Create a new file using H5F_ACC_TRUNC access,
     * default file creation properties, and default file
     * access properties.
     */
    //file = H5Fcreate(H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    file = H5Fcreate(outfile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /*
     * Describe the size of the array and create the data space for fixed
     * size dataset.
     */

    dataspace = H5Screate_simple(RANK, dimsf, NULL);
    //dataspace = H5Screate_scalar(RANK, dimsf, NULL);

    /*
     * Define datatype for the data in the file.
     * We will store little endian INT numbers.
     */
    datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
    status   = H5Tset_order(datatype, H5T_ORDER_LE);

    /*
     * Create a new dataset within the file using defined dataspace and
     * datatype and default dataset creation properties.
     */
    dataset = H5Dcreate2(file, DATASETNAME, datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /*
     * Write the data to the dataset using default transfer properties.
     */
    //std::cout << "checking value of data " << data[0] << std::endl;
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    /*
     * Close/release resources.
     */
    H5Sclose(dataspace);
    H5Tclose(datatype);
    H5Dclose(dataset);
    H5Fclose(file);

    // this definition might have to change
    //double     data[this->grvy_parser.N]; /* data to write */
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

    // f is the input, for the manufactured solution, input == solution...compute error
    gsl_vector_sub(this->f, this->u);

    // open output file
    std::ofstream outfile;
    outfile.open("../output/" + solver_method_to_string() + "_" + std::to_string(this->dim) + "D_" + std::to_string(this->order) + "_order_N_" + std::to_string(this->num_nodes));
    // compute L2 norm of the difference and output to outfile
    double l2norm{gsl_blas_dnrm2(this->f)};
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
    this->iterative_solve(); // Solve Au = f iteratively
    // should probably define the function to use similarly to how we picked the iterative solver to use....check that out.
    if (this->dim == 1){
        this->save_hdf5_1d_data(outfile);
    } else if(this->dim == 2){
        this->save_hdf5_2d_data(outfile);
    } else{
        std::cout<< "Can only solve 1D or 2D system!" << std::endl;
        exit(1);
    }
    //this->save_solution(outfile); // save the solution to an HDF5 file.
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
