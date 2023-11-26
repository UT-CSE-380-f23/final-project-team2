# Building the Code

Below are some notes when building the code. Some of these will be irrelevant to the usage of the final version of the code, but are useful to understand if custom versions of GSL, libgrvy, petsc, and other libraries are used, and if these libraries are set up on a personal computer.

## CLANG vs. g++
MAC and Linux users have different C++ compilers installed by default. MAC users have CLANG installed, while Linux g++. There are a few minor differences between the two in how they build code:

1. The order of the libraries included in compilation/linking MATTER in g++ - but not in CLANG.
    - ex: ```g++ ... -lgslcblas -lgsl -lm``` will throw an error in g++, but not in CLANG. GSL needs to be loaded/linked before the GSL-specific CBLAS libraries can be referenced. `g++ ... -lgsl -lgslcblas -lm` will work in g++. 
2. Where the libraries are included in the compilation/linking step matters
    - it seems like `-lgslcblas -lgsl -lm` needs to be included at the **end** of both the compilation and linking steps in order for the libraries to be recognized when using g++, this does not appear to be an issue in CLANG.