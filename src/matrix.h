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

    // Solve the system
    iteration_element_based(N_arg, A, f, u, true);

    /* Free the space */
    gsl_spmatrix_free(A);
    gsl_vector_free(f);
    gsl_vector_free(u);

    std::cout << "hello from the end of oneDimSecOrd" << "\n";
}

/*

One dimension - fourth order. 

*/

void oneDimFouOrd(int N_arg){
    std::cout << "hello from the beggining of oneDimFouOrd" << "\n";

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
    gsl_spmatrix_set(A, 0, 0, -24.0);
    gsl_spmatrix_set(A, 0, 1, 12.0);

    /* construct second row */
    gsl_spmatrix_set(A, 1, 0, 16.0);
    gsl_spmatrix_set(A, 1, 1, -30.0);
    gsl_spmatrix_set(A, 1, 2, 16.0);
    gsl_spmatrix_set(A, 1, 3, -1.0);

    /* construct rows [1:n-2] */
    for (i = 2; i < n - 2; ++i)
    {
        gsl_spmatrix_set(A, i, i + 2, -1.0);
        gsl_spmatrix_set(A, i, i + 1, 16.0);
        gsl_spmatrix_set(A, i, i, -30.0);
        gsl_spmatrix_set(A, i, i - 1, 16.0);
        gsl_spmatrix_set(A, i, i - 2, -1.0);
    }

    /* construct second to last row */
    gsl_spmatrix_set(A, n-2, n-1, 16.0);
    gsl_spmatrix_set(A, n-2, n-2, -30.0);
    gsl_spmatrix_set(A, n-2, n-3, 16.0);
    gsl_spmatrix_set(A, n-2, n-4, -1.0);

    /* construct last row */
    gsl_spmatrix_set(A, n - 1, n - 1, -24.0);
    gsl_spmatrix_set(A, n - 1, n - 2, 12.0);

    /* print the matrix */
    for (int j = 0; j < n; ++j){
        for (i = 0; i < n; ++i){
            // std::cout << i << "," << j << " "; 
            std::cout << gsl_spmatrix_get(A, j, i) << " ";
        }
        std::cout << "\n";
    }

    /* scale by h^2 */
    gsl_spmatrix_scale(A, -1.0 * k / (12 * h * h));

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
    iteration_element_based(N_arg, A, f, u, false);

    /* Free the space */
    gsl_spmatrix_free(A);
    gsl_spmatrix_free(C);
    gsl_vector_free(f);
    gsl_vector_free(u);

    std::cout << "hello from the end of oneDimFouOrd" << "\n";
}

/*

Two dimension - second order. 

*/

int TwoDimSecOrd(int N_arg){
    std::cout << "hello from the beggining of TwoDimSecOrd" << "\n";

    double constant       = 1.0 / (8 * M_PI * M_PI);     /* constant in front of matrix */
    const size_t N = N_arg;                       /* number of grid points */
    const size_t k = N-2;                         /* subtract 2 to exclude boundaries */
    const int    n = k*k;                         /* size of the sytem/matrix */
    const double h = 1.0 / (N-1);                 /* grid spacing */
    
    // Print the grid spacing:
    std::cout << "grid spacing: " << h  << "\n";

    gsl_spmatrix *A = gsl_spmatrix_alloc(n ,n); /* triplet format */
    gsl_vector *f = gsl_vector_alloc(n);        /* right hand side vector */
    gsl_vector *u = gsl_vector_alloc(n);        /* solution vector */

    /* construct the sparse matrix for the finite difference equation */
    double xi, yj, f_val;

    for (int i = 0; i < k; i++){
    for (int j = 0; j < k; j++){
      // Set the matrix
      gsl_spmatrix_set(A, i*k+j, i*k+j, -4.0);

      if(i != 0)  {gsl_spmatrix_set(A, i*k+j, (i-1)*k+j, 1.0);}
      if(i != k-1){gsl_spmatrix_set(A, i*k+j, (i+1)*k+j, 1.0);}

      if(j != 0)  {gsl_spmatrix_set(A, i*k+j, i*k+(j-1), 1.0);}
      if(j != k-1){gsl_spmatrix_set(A, i*k+j, i*k+(j+1), 1.0);}

      // Set the RHS
      xi = (i + 1) * h;
      yj = (j + 1) * h;
      f_val = sin(2.0 * M_PI * xi)*sin(2.0 * M_PI * yj);
      gsl_vector_set(f, i*k+j, f_val);

    }
    }

    // Scale RHS
    gsl_vector_scale(f, (-1.0)*h*h/constant);

    /* Print the matrix */
    double val = 0;

    for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      val = gsl_spmatrix_get(A, i, j);

      if (val > 0){
        std::cout << "+" << val << " ";
      }else if(val == 0){
        std::cout << " " << val << " ";
      }else{
        std::cout << val << " ";
      }
    }
      std::cout << "\n";
    }

    return 0;

    // // Solve the system
    // iteration_element_based(N_arg, A, f, u, true);

    /* Free the space */
    gsl_spmatrix_free(A);
    gsl_vector_free(f);
    gsl_vector_free(u);

    std::cout << "hello from the end of TwoDimSecOrd" << "\n";
}



#endif // MATRIX_H