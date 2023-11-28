// 11/6/23
// Kenneth Meyer, Rodrigo Gonzalez, Jenil Shah

#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <iostream>
#include "oned_solver.h"
#include "twod_solver.h"
//#include <boost/program_options.hpp>
// #include "parsing.h"

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

int main(int argc, char *argv[]){

    // initialize libgrvy parser and parse input file
    GrvyParser grvy_parser(argc,argv);

    // solve 1D or 2D system
    if (grvy_parser.DIM == 1){
        std::cout << "solving a system in 1D!! " << std::endl;
        OneDSolver sl{grvy_parser.N, grvy_parser.DIM, grvy_parser.solver, grvy_parser.ORDER, grvy_parser.verify, grvy_parser.mode};
        sl.system_solve();
    } else if (grvy_parser.DIM == 2){
        std::cout << "solving a system in 2D!! " << std::endl;
        TwoDSolver sl{grvy_parser.N, grvy_parser.DIM, grvy_parser.solver, grvy_parser.ORDER, grvy_parser.verify, grvy_parser.mode};
        sl.system_solve();
    }
    else {
        std::cout << "You don't have a " << grvy_parser.DIM << " solver implemented yet!" << std::endl;
    }
    return 0;
}
