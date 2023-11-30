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
    GrvyParser grvy_parser(argv[1]); // To take input file as an argument
    grvy_parser.parse(); // Parsing
    // can add some tests/checks using this output; need to toggle verification/debug mode vs. normal mode.

    // do this in debug mode
    if(grvy_parser.mode==true)
    {
        std::cout << "Checking variables from main " << std::endl;
        std::cout << "verify    = " << grvy_parser.verify << std::endl;
        std::cout << "mode      = " << grvy_parser.mode << std::endl;
        std::cout << "N         = " << grvy_parser.N << std::endl;
        // Adding verbosity for debug mode
    }
    
    // size_t num_nodes = atoi(argv[1]);
    //std::cout<<"We are in main with nodes "<<num_nodes<<std::endl;
    
    // initialize a 1-D solver with the libgrvy parser object
    // OneDSolver sl{num_nodes, false};
    // sl.system_solve(); // Example without libgrvy

    // ideally, we would want to pass the grvy_parser object to this, but I'm going to avoid that for now.
    // why?: because it would just make the code too complex/require too many changes at this point.
    if(grvy_parser.N < 4)
        std::cout << "num_nodes too low!"<<std::endl;
    else
    {
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
            std::cout << "You don't have a " << grvy_parser.DIM << "D solver implemented yet!" << std::endl;
        }
    }
    return 0;
}
