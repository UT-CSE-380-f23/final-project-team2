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