#ifndef output
#define output


#include<string>
#include<cstdio>
#include<cstdlib>
#include<iostream>
#include<hdf5.h>
#include<parsing.h>

/**
 * @brief Parses command line arguments 
 * 
 * @param ac input from main program
 * @param av input from main program
 * @return int 
 */
class Output {
    public:
        const char* outfile;

        // the constructor!
        Output(GrvyParser& grvy_parser);
        Output(GrvyParser& grvy_parser, const char* outfile);
        void dump();

    private:
        GrvyParser &grvy_parser;
};

#endif