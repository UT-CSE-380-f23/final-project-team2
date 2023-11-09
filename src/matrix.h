#ifndef MATRIX_H
#define MATRIX_H

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

Performs the Jacobi iteration

*/

// void jacobi(int N_arg, gsl_spmatrix *M, gsl_vector *b, gsl_vector *x){
//     std::cout << "hello from the beggining of jacobi" << "\n";

//     /* Some variables */
//     const size_t N = N_arg;                                     /* number of grid points */
//     const size_t n = N-2;                                       /* subtract 2 to exclude boundaries */

//     gsl_spmatrix *D_inv     = gsl_spmatrix_alloc(n ,n);         /* triplet format */
//     gsl_spmatrix *H         = gsl_spmatrix_alloc(n ,n);         /* triplet format */

//     gsl_spmatrix *EF;       /* compressed format */
//     gsl_spmatrix *D_inv_c;  /* compressed format */
//     gsl_spmatrix *H_c;      /* compressed format */

//     gsl_vector *r       = gsl_vector_alloc(n);                    /* right hand side vector */
//     gsl_vector *x_prev  = gsl_vector_alloc(n);                    /* used for iteration */

//     size_t i;

//     /*  Update rule:
//             -> x_k+1 = H x_k + r
//         With:
//             -> H = D_inv * EF
//             -> r = D_inv * b
//         And:
//             -> EF = -(M - diag(M))
//     */

//     /* Create D_inv by iterating over the diagonal of M */
//     /* 
//         Update M by deleting its diagonal:
//         Note that this process will destroy the original matrix M
//     */

//     for (i = 0; i < n; ++i)
//     {
//       gsl_spmatrix_set(D_inv, i, i, 1.0/(gsl_spmatrix_get(M, i, i)));
//       gsl_spmatrix_set(M, i, i, 0);
//     }

//     /* Print these matrices */
//     for (int j = 0; j < n; ++j){
//         for (i = 0; i < n; ++i){
//             // std::cout << i << "," << j << " "; 
//             std::cout << gsl_spmatrix_get(M, i, j) << " ";
//         }
//         std::cout << "\n";
//     }

//     /* Print these matrices */
//     for (int j = 0; j < n; ++j){
//         for (i = 0; i < n; ++i){
//             // std::cout << i << "," << j << " "; 
//             std::cout << gsl_spmatrix_get(D_inv, i, j) << " ";
//         }
//         std::cout << "\n";
//     }

//     /* scale M by -1 */
//     gsl_spmatrix_scale(M, -1.0);

//     /* convert to compressed column format */
//     EF  = gsl_spmatrix_ccs(M);
//     D_inv_c = gsl_spmatrix_ccs(D_inv);
//     H_c = gsl_spmatrix_ccs(H);

//     /* Construct H and r using the sparse BLAS Support of GSL */
//     gsl_spblas_dgemm(1 , D_inv_c, EF, H_c);
//     gsl_spblas_dgemv(CblasNoTrans, 1, D_inv_c, b, 0.0, r);

//     /* Print this matrix */
//     for (int j = 0; j < n; ++j){
//         for (i = 0; i < n; ++i){
//             // std::cout << i << "," << j << " "; 
//             std::cout << gsl_spmatrix_get(H_c, i, j) << " ";
//         }
//         std::cout << "\n";
//     }
//     std::cout << "\n";

//     /* Print the RHS */
//     for (i = 0; i < n; ++i){
//         // std::cout << i << "," << j << " "; 
//         std::cout << gsl_vector_get(r, i) << " ";
//     }
//     std::cout << "\n";


//     /* Apply the Jacobi iteration */
//     const double tol    = 1.0e-6;     /* solution relative tolerance */
//     const size_t max_iter = 10000; /* maximum iterations */
//     int iter = 0;
//     double residual     = 100.0;

//     /* initialiaze x */
//     gsl_vector_set_zero(x);
//     /* initial guess x_prev = 0 */
//     gsl_vector_set_zero(x_prev);
    
//     /* iterate */
//     while (tol < residual && ++iter < max_iter){
//         // Perform the matrix vector multiplication (H x_k)
//         gsl_spblas_dgemv(CblasNoTrans, 1.0, H_c, x_prev, 0.0, x);
//         // Add r vector (H x_k + r)
//         gsl_vector_add(x, r);

//         // // Print x
//         // for (i = 0; i < n; i++)
//         // {
//         // std::cout << "x_" << i << "=" << gsl_vector_get (x_prev, i) << "\n";
//         // }

//         // Swap them to continue the iteration
//         gsl_vector_swap(x, x_prev);

//         // Update the residual and print it
//         gsl_vector_sub(x, x_prev);
//         residual = gsl_blas_dnrm2(x);
//         std::cout << "residual: " << residual << "\n";

//     }

//     if (tol >= residual){
//         std::cout << "The iteration coverged to the following solution: " << "\n";
//         // Print x
//         for (i = 0; i < n; i++)
//         {
//         std::cout << "x_" << i << "=" << gsl_vector_get (x_prev, i) << "\n";
//         }
//     }

//     /* Free variables */
//     gsl_spmatrix_free(D_inv);
//     gsl_spmatrix_free(H);
    
//     gsl_spmatrix_free(EF);
//     gsl_spmatrix_free(D_inv_c);
//     gsl_spmatrix_free(H_c);

//     gsl_vector_free(r);
//     gsl_vector_free(x_prev);

//     std::cout << "hello from the beggining of jacobi" << "\n";
// }

/*

One dimension - second order. 

*/

void oneDimSecOrd(int N_arg){
    std::cout << "hello from the beggining of oneDimSecOrd" << "\n";

    double k       = 1.0 / (4 * M_PI * M_PI);     /* constant in front of matrix */
    const size_t N = N_arg;                       /* number of grid points */
    const size_t n = N-2;                         /* subtract 2 to exclude boundaries */
    const double h = 1.0 / (N-1);                 /* grid spacing */
    
    // Print the grid spacing:
    std::cout << "grid spacing: " << h  << "\n";

    gsl_spmatrix *A = gsl_spmatrix_alloc(n ,n); /* triplet format */
    gsl_spmatrix *C;                            /* compressed format */
    gsl_vector *f = gsl_vector_alloc(n);        /* right hand side vector */
    gsl_vector *u = gsl_vector_alloc(n);        /* solution vector */
    size_t i;

    /* construct the sparse matrix for the finite difference equation */

    /* construct first row */
    gsl_spmatrix_set(A, 0, 0, -2.0);
    gsl_spmatrix_set(A, 0, 1, 1.0);

    /* construct rows [1:n-2] */
    for (i = 1; i < n - 1; ++i)
    {
      gsl_spmatrix_set(A, i, i + 1, 1.0);
      gsl_spmatrix_set(A, i, i, -2.0);
      gsl_spmatrix_set(A, i, i - 1, 1.0);
    }

    /* construct last row */
    gsl_spmatrix_set(A, n - 1, n - 1, -2.0);
    gsl_spmatrix_set(A, n - 1, n - 2, 1.0);

    /* scale by h^2 */
    gsl_spmatrix_scale(A, -1.0 * k / (h * h));

    /* construct right hand side vector */
    for (i = 0; i < n; ++i)
    {
      double xi = (i + 1) * h;
      double fi = sin(2.0 * M_PI * xi);
      gsl_vector_set(f, i, fi);
    }

    /* convert to compressed column format (May not be needed))*/ 
    C = gsl_spmatrix_ccs(A); 

    // Solve the system
    jacobi(N_arg, A, f, u);

    /* Free the space */
    gsl_spmatrix_free(A);
    gsl_spmatrix_free(C);
    gsl_vector_free(f);
    gsl_vector_free(u);

    std::cout << "hello from the end of oneDimSecOrd" << "\n";
}

#endif // MATRIX_H