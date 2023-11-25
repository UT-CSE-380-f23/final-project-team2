// 11/6/23
// Kenneth Meyer, Rodrigo Gonzalez, Jenil Shah

#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <iostream>
#include "oned_solver.h"
//#include <boost/program_options.hpp>
//#include "parsing.h"

//#include "iterative_methods.h"
//#include "matrix.h"
#include "parsing.h"

// variables
//int n_samples;
//int DIM, ORDER, N, VERIFY;
//SOLVER, MODE;

// used to parse command-line arguments
//namespace po = boost::program_options;

/* input options

    --dim : 1, 2 (default 1)
        Choice between 1D or 2D problems
    --order : 2, 4 (defaul)
        Choice between 2nd or 4th-order finite differences
    --N : int
        Mesh sizing options
    --verify : bool ???
        Option to run in a verification mode
    --solver : jacobi, gauss
        Choice between Jacobi or Gauss-Seidel solvers
    --mode : standard, debug
        Option to control standard output mode (e.g. standard vs debug)
*/

//int main(int n_args,char *argv[]){
int main(int argc, char *argv[]){

    GrvyParser grvy_parser(argc,argv);

    // can add some tests/checks using this output; need to toggle verification/debug mode vs. normal mode.

    // do this in debug mode
    std::cout << "Checking variables from main " << std::endl;
    std::cout << "verify    = " << grvy_parser.verify << std::endl;
    std::cout << "mode      = " << grvy_parser.mode << std::endl;
    std::cout << "N         = " << grvy_parser.N << std::endl;
    return 0;

    //size_t num_nodes = atoi(argv[1]);
    //std::cout<<"We are in main with nodes "<<num_nodes<<std::endl;
    
    // initialize a 1-D solver with the libgrvy parser object
    //OneDSolver sl{num_nodes, false};
    //sl.system_solve(); // Example without libgrvy
}