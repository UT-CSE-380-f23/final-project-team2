#ifndef ITERATIVE_METHODS_H
#define ITERATIVE_METHODS_H

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

int jacobi(int size_of_system, gsl_spmatrix *M, gsl_vector *b, gsl_vector *x){
    const size_t n = size_of_system;                            /* subtract 2 to exclude boundaries */

    gsl_spmatrix *D_inv     = gsl_spmatrix_alloc(n ,n);         /* triplet format */
    gsl_spmatrix *H         = gsl_spmatrix_alloc(n ,n);         /* triplet format */

    gsl_spmatrix *EF;       /* compressed format */
    gsl_spmatrix *D_inv_c;  /* compressed format */
    gsl_spmatrix *H_c;      /* compressed format */

    gsl_vector *r       = gsl_vector_alloc(n);                    /* right hand side vector */
    gsl_vector *x_prev  = gsl_vector_alloc(n);                    /* used for iteration */

    size_t i;

    /*  Update rule:
            -> x_k+1 = H x_k + r
        With:
            -> H = D_inv * EF
            -> r = D_inv * b
        And:
            -> EF = -(M - diag(M))
    */

    /* Create D_inv by iterating over the diagonal of M */
    /* 
        Update M by deleting its diagonal:
        Note that this process will destroy the original matrix M
    */

    for (i = 0; i < n; ++i)
    {
      gsl_spmatrix_set(D_inv, i, i, 1.0/(gsl_spmatrix_get(M, i, i)));
      gsl_spmatrix_set(M, i, i, 0);
    }

    /* Print these matrices */
    /*
    for (int j = 0; j < n; ++j){
        for (i = 0; i < n; ++i){
            // std::cout << i << "," << j << " "; 
            std::cout << gsl_spmatrix_get(M, i, j) << " ";
        }
        std::cout << "\n";
    }
    */

    /* Print these matrices */
    /*
    for (int j = 0; j < n; ++j){
        for (i = 0; i < n; ++i){
            // std::cout << i << "," << j << " "; 
            std::cout << gsl_spmatrix_get(D_inv, i, j) << " ";
        }
        std::cout << "\n";
    }
    */

    /* scale M by -1 */
    gsl_spmatrix_scale(M, -1.0);

    /* convert to compressed column format */
    EF  = gsl_spmatrix_ccs(M);
    D_inv_c = gsl_spmatrix_ccs(D_inv);
    H_c = gsl_spmatrix_ccs(H);

    /* Construct H and r using the sparse BLAS Support of GSL */
    gsl_spblas_dgemm(1 , D_inv_c, EF, H_c);
    gsl_spblas_dgemv(CblasNoTrans, 1, D_inv_c, b, 0.0, r);

    /* Print this matrix */

    /*
    for (int j = 0; j < n; ++j){
        for (i = 0; i < n; ++i){
            // std::cout << i << "," << j << " "; 
            std::cout << gsl_spmatrix_get(H_c, i, j) << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
    */

    
    /* Print the RHS */
    /*
    for (i = 0; i < n; ++i){
        // std::cout << i << "," << j << " "; 
        std::cout << gsl_vector_get(r, i) << " ";
    }
    std::cout << "\n";
    */

    /* Apply the Jacobi iteration */
    const double tol    = 1.0e-6;       /* solution relative tolerance */
    const size_t max_iter = 100000;        /* maximum iterations */
    int iter = 0;
    double residual             = 100.0;
    double prevResidual         = residual + 1.0;

    /* initialiaze x */
    gsl_vector_set_zero(x);
    /* initial guess x_prev = 0 */
    gsl_vector_set_zero(x_prev);
    
    /* iterate */
    while (tol < residual && ++iter < max_iter ){
        // Perform the matrix vector multiplication (H x_k)
        gsl_spblas_dgemv(CblasNoTrans, 1.0, H_c, x_prev, 0.0, x);
        // Add r vector (H x_k + r)
        gsl_vector_add(x, r);

        // // Print x
        // for (i = 0; i < n; i++)
        // {
        // std::cout << "x_" << i << "=" << gsl_vector_get (x_prev, i) << "\n";
        // }

        // Swap them to continue the iteration
        gsl_vector_swap(x, x_prev);

        // Update the residual and print it
        gsl_vector_sub(x, x_prev);
        residual = gsl_blas_dnrm2(x);
        // std::cout << "residual: " << residual << "\n";

        // Only continue if the resiudal is decreasing
        if (prevResidual < residual){
            break;
        }

        // Update the previous residual
        prevResidual = residual;

    }

    std::cout << "The finalized with the following residual: " << residual << "\n";
    // Print x
    // for (i = 0; i < n; i++){ std::cout << "x_" << i << "=" << gsl_vector_get (x_prev, i) << "\n";  }

    /*
        Now we compute the error of the solution. 
    */

    // Compute the error vector
    gsl_vector_sub(b, x_prev);

    // Print error vector
    // for (i = 0; i < n; i++){ std::cout << "error_" << i << "=" << gsl_vector_get (b, i) << "\n";}
    
    // Print the l2 error
    std::cout << "\n \n" << "L2 error: " << gsl_blas_dnrm2(b) << "\n \n";

    /* Free variables */
    gsl_spmatrix_free(D_inv);
    gsl_spmatrix_free(H);
    
    gsl_spmatrix_free(EF);
    gsl_spmatrix_free(D_inv_c);
    gsl_spmatrix_free(H_c);

    gsl_vector_free(r);
    gsl_vector_free(x_prev);

    return 0;
}


/*

Jacobi iteration

*/

/*

Performs the following pseudo-code:


k = 0
while convergence not reached do
    for i := 1 step until n do
        σ = 0
        for j := 1 step until n do
            if j ≠ i then
                σ = σ + aij xj(k)
            end
        end
        xi(k+1) = (bi − σ) / aii
    end
    increment k
end

*/

int jacobi_element_based(int size_of_system, gsl_spmatrix *M, gsl_vector *b, gsl_vector *x){
    const size_t n = size_of_system;                                       /* The size of the system */

    /* Allocate memory for x_prev*/
    gsl_vector *x_prev  = gsl_vector_alloc(n);       /* used for iteration */


    /* Variables used for the Jacobi iteration */
    const double tol    = 1.0e-6;                   /* solution relative tolerance */
    const int max_iter = 10000;                     /* maximum iterations */
    int k = 0;
    double sigma = 0;
    double residual             = 100.0;
    double prevResidual         = residual + 1.0;

    /* initialiaze x */
    gsl_vector_set_zero(x);
    /* initial guess x_prev = 0 */
    gsl_vector_set_zero(x_prev);

    int i = 0;
    int j = 0;

    /* iterate */
    while (tol < residual && ++k < max_iter){

        /* Perform Jacobi*/
        for (i = 0; i < n; i++){
            sigma = 0;
            
            for (j = 0; j < n; j++){
                if(i != j){
                    sigma = sigma + gsl_spmatrix_get(M, i, j) * gsl_vector_get(x_prev, j);
                }
            }
            
            gsl_vector_set(x, i, (gsl_vector_get(b, i) - sigma)/gsl_spmatrix_get(M, i, i));
        }

        // Swap them to continue the iteration
        gsl_vector_swap(x, x_prev);

        // Update the residual and print it
        gsl_vector_sub(x, x_prev);
        residual = gsl_blas_dnrm2(x);
        std::cout << "residual: " << residual << "\n";

        // Only continue if the resiudal is decreasing
        if (prevResidual < residual){
            break;
        }

        // Update the previous residual
        prevResidual = residual;

    }

    /* Print x */
    for (i = 0; i < n; i++)
    {
        std::cout << "x_" << i << "=" << gsl_vector_get (x_prev, i) << "\n";
    }

    /* Free memory */

    return 0;
}


/*

Can perform both Jacobi and Gauss–Seidel iteration

*/

/*

For Gauss-Seidel, it performs the following pseudo-code:


algorithm Gauss–Seidel method is
    inputs: A, b
    output: φ

    Choose an initial guess φ to the solution
    repeat until convergence
        for i from 1 until n do
            σ ← 0
            for j from 1 until n do
                if j ≠ i then
                    σ ← σ + aijφj
                end if
            end (j-loop)
            φi ← (bi − σ) / aii
        end (i-loop)
        check if convergence is reached
    end (repeat)

*/

int iteration_element_based(int size_of_system, gsl_spmatrix *M, gsl_vector *b, gsl_vector *x, bool jacOrGS){
    std::cout << " Starting iteration solver: " << "\n";

    /* Size of the system */
    const size_t n = size_of_system;                                      

    /* Allocate memory for x_prev*/
    gsl_vector *x_prev  = gsl_vector_alloc(n);       /* used for iteration */


    /* Variables used for the Gauss-Seidel iteration */
    const double tol    = 1.0e-6;                   /* solution relative tolerance */
    const int max_iter = 1000000;                     /* maximum iterations */
    int k = 0;
    double sigma = 0;
    double residual             = 100.0;
    double prevResidual         = residual + 1.0;

    /* initialiaze x at zero*/
    gsl_vector_set_zero(x);

    // Set the initial value
    for (int i = 0; i < n; i++)
    {
        gsl_vector_set (x, i, 0.5);
    }

    /* initialiaze x_prev at zero*/
    gsl_vector_set_zero(x_prev);

    // Iteration variables
    int i = 0;
    int j = 0;

    /* iterate */
    while (tol < residual && ++k < max_iter){

        if (jacOrGS == true){
            /* Perform Jacobi */
            for (i = 0; i < n; i++){
                sigma = 0;
            
                for (j = 0; j < n; j++){
                    if(i != j){
                        sigma = sigma + gsl_spmatrix_get(M, i, j) * gsl_vector_get(x_prev, j);
                    }
                }
            
            gsl_vector_set(x, i, (gsl_vector_get(b, i) - sigma)/gsl_spmatrix_get(M, i, i));
            }
        }else {
            /* Perform Gauss-Seidel */
            for (i = 0; i < n; i++){
                sigma = 0;
            
                for (j = 0; j < n; j++){
                    if(i != j){
                        sigma = sigma + gsl_spmatrix_get(M, i, j) * gsl_vector_get(x, j);
                    }
                }
                gsl_vector_set(x, i, (gsl_vector_get(b, i) - sigma)/gsl_spmatrix_get(M, i, i));
            }
        }

        // Update the residual and print it
        gsl_vector_sub(x_prev, x);
        residual = gsl_blas_dnrm2(x_prev);
        std::cout << "residual: " << residual << "\n";

        // Make a copy of x into x_prev for the following iteration
        gsl_vector_memcpy(x_prev, x);

        // Only continue if the resiudal is decreasing
        if (prevResidual < residual){
            break;
        }

        // Update the previous residual
        prevResidual = residual;

    }

    // Print x
    for (i = 0; i < n; i++)
    {
        std::cout << "x_" << i << "=" << gsl_vector_get (x, i) << "  RHS_" << i << "=" << gsl_vector_get (b, i) << "\n";
    }

    // Compute the error vector
    gsl_vector_sub(b, x);

    std::cout << "\n";

    // Print error vector
    for (i = 0; i < n; i++)
    {
        std::cout << "error_" << i << "=" << gsl_vector_get (b, i) << "\n";
    }

    
    // Print the l2 error
    std::cout << "\n \n" << "L2 error: " << gsl_blas_dnrm2(b) << "\n \n";

    // Free space
    gsl_vector_free(x_prev);
    return 0;
}


#endif // ITERATIVE_METHODS_H