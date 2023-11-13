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

using namespace GRVY;
using namespace std;

/**
 * @brief Parses runtime options from input.txt file
 * 
 * @param argc input from main()
 * @param argv input from main()
 * @return int 
 */
class GrvyParser {
    public:
        std::string solver;
        int N, DIM, ORDER;
        int verify, mode;
        GrvyParser(int argc, char **argv){
    
            GRVY_Input_Class iparse; // Input parsing object

            // Initialize/read the file 
            if(! iparse.Open("../heat-input.txt"))
                exit(1);  

            // Read in option variables
            if( iparse.Read_Var("verify",&verify,0) )
                printf("--> %-11s = %i\n","verify",verify);
        
            if( iparse.Read_Var("mode",&mode,0) )
                printf("--> %-11s = %i\n","mode",mode);
        
            // Read in domain variables
            if( iparse.Read_Var("domain/N",&N,3) )
                printf("--> %-11s = %i\n","N",N);
        
            if( iparse.Read_Var("domain/DIM",&DIM,1) )
                printf("--> %-11s = %i\n","DIM",DIM);

            // Read in solver variables
            if( iparse.Read_Var("solver/solver",&solver,"gauss-seidel") )
                cout << "--> solver      = " << solver << endl;
            if( iparse.Read_Var("solver/ORDER",&ORDER,2) )
                printf("--> %-11s = %i\n","ORDER",ORDER);

            // Close file and exit
            iparse.Close();

        }
};


// example from here: https://www.boost.org/doc/libs/1_63_0/doc/html/program_options/tutorial.html#idp523371328
// is used to help code the function below:
/**
 * @brief Parses command line arguments 
 * 
 * @param ac input from main program
 * @param av input from main program
 * @return int 
 */
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