#ifndef INTEGRATION_H
#define INTEGRATION_H

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

/*

Given a function and extra parameters passed as reference, it computes the integral of such function.

*/

// Computes the integral of the given function
double integral(double a, double b, double (*pf)(double, void *), void * p_parameters, int n){
    // Set up the integration
    gsl_integration_glfixed_table * w = gsl_integration_glfixed_table_alloc(n);
    double result;

    // Create the function and parameters
    gsl_function F;
    F.function = pf;
    F.params = p_parameters;

    // Perform the integration
    result = gsl_integration_glfixed(&F, a, b, w);

    // frees the memory assocated with the workspace
    gsl_integration_glfixed_table_free(w);

    // return the result
    return result;
}

/*

Given two functions passed as reference, it returns the product of these two functions
evaluated at point x.

*/ 

// Struct that contains the two input functions for the inner product, along with their parameters
struct TwoFunctions {
    // The two functions to whose inner product will be computer
    double (*f_1)(double, void *);
    double (*f_2)(double, void *);

    // The two set of parameters that the two functions may use in their computations
    void * params_1;
    void * params_2;
};

// Creates the function to be integrated
double product(double x, void * pointer_args){
    // Cast the void* back to the struct type
    TwoFunctions* structPtr = static_cast<TwoFunctions*>(pointer_args);

    // Unpack the two functions
    double (*f)(double, void *) = structPtr->f_1;
    double (*g)(double, void *) = structPtr->f_2;

    // Unpack the two parameters
    void * p1 = structPtr->params_1;
    void * p2 = structPtr->params_2;

    // Return the product of the two functions
    return f(x, p1)*g(x, p2);
}

// Computes the inner product of two functions passed as arguments
double inner_product (double a, double b, double (*f)(double, void *), double (*g)(double, void *), void * p_1, void * p_2, int n){ 
    
    // Create an instance of MyStruct, and fill it
    TwoFunctions functions;

    // The two functions
    functions.f_1 = f;
    functions.f_2 = g;

    // The two parameters that are passed to the functions
    functions.params_1 = p_1;
    functions.params_2 = p_2;

    // Compute the integral
    return integral(a, b, product, &functions, n);

}

/*

Given a functions passed as reference, it returns the L^2 norm of this function

*/ 

// Computes the inner product of two functions passed as arguments
double norm_l2 (double a, double b, double (*f)(double, void *), void * params, int n){ 
    return sqrt(inner_product(a, b, f, f, params, params, n));
}


#endif // INTEGRATION_H
