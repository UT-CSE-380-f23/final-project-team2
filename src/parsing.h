/*

    Header file for parsing command line input and parsing logfiles.

*/
#ifndef parsing
#define parsing

#include<string>
#include<cstdio>
#include<cstdlib>
#include<iostream>
#include<grvy.h>
//#include<hdf5.h>


// example from here: https://www.boost.org/doc/libs/1_63_0/doc/html/program_options/tutorial.html#idp523371328
// is used to help code the function below:
/**
 * @brief Parses command line arguments 
 * 
 * @param ac input from main program
 * @param av input from main program
 * @return int 
 */

class GrvyParser {
    public:
        //std::string solver;
        int N;
        int DIM, ORDER;
        bool verify, mode, solver, USE_PETSC;
        const char* input_file;

        // the constructor!
        GrvyParser();
        GrvyParser(const char* input_file);
        void parse();
};
/*
int deprecated_parsing(int ac, char** av){

    
        Parses command line arguments, which are as follows:
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
    

    // declare all supported options - some of these should be made optional
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "show help message")
        ("dim", po::value<int>(), "Set problem dimension, 1 or 2")
        ("order", po::value<int>(), "Finite Difference order, 2 or 4")
        ("N", po::value<int>(), "Set mesh size")
        ("verify", po::value<int>(), "Set verification mode, 0 or 1")
        ("solver", po::value<std::string>(), "Chose solver, 'jacobi' or 'gauss-seidel'")
        ("mode", po::value<std::string>(), "Chose mode, 'standard' or 'debug'")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);

    // return help message and do not run code if --help is given
    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 1;
    }

    // set default options for some variables if none are chosen
    if (vm.count("dim")) {
        std::cout << "Solving problem in " << vm["dim"].as<int>() << "D" << "\n";
    } else {
        std::cout << "Must specify problem dimension!!";
    }


    return 0;
}
*/
#endif