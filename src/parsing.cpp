#include "parsing.h"

// includes classes for libgrvy parsing and parsing with hdf5 files


using namespace GRVY;
//using namespace std;

// GrvyParser constructor (only method of the class)
GrvyParser::GrvyParser(int argc, char **argv){

    GRVY_Input_Class iparse; // Input parsing object

    // likely want to check for debug/verify mode when printing things...
    // also should compare speed using printf vs. std::cout

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
    if( iparse.Read_Var("solver/solver",&solver,0) )
        std::cout << "--> solver      = " << solver << std::endl;
    if( iparse.Read_Var("solver/ORDER",&ORDER,2) )
        printf("--> %-11s = %i\n","ORDER",ORDER);

    // Close file and exit
    iparse.Close();
};