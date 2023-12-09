[![CI](https://github.com/UT-CSE-380-f23/final-project-team2/actions/workflows/ci_testing.yaml/badge.svg)](https://github.com/UT-CSE-380-f23/final-project-team2/actions/workflows/ci_testing.yaml)

# final-project-team2
final-project-team2 created by GitHub Classroom
Authors:  Rodrigo Gonzalez, Kenneth L  Meyer, Jenil Shah

## Installation Notes
todo: check if the code should be able to be built and ran on *any* system, or just lonestar6. I.e. does the `Makefile` just need to work properly for any user on ls6? Or for users beyond ls6 as well?

## Running the code on docker
1. Clone the repo
2. Run the following command :
`docker run -v <cloned_folder>:/home/test jenilut/cse380_team2 csh <cloned_folder>/Docker/docker_make.sh`

### Installing libgrvy
Documentation of a local copy of [libgrvy](https://github.com/hpcsi/grvy) is shown below. The link to the version installed can be found [here](https://github.com/hpcsi/grvy/releases/tag/0.38.0)

```shell
wget https://github.com/hpcsi/grvy/releases/download/0.38.0/grvy-0.38.0.tar.gz
cd grvy-0.38.0/
module load boost # libgrvy is built on boost/it is a required module
./configure --prefix=$WORK/soft/bin/grvy # choose your installation location, might need to create it yourself? not 100% sure.
make
make check
make install
```


### Running a libgrvy application:
The following commands were used to successfully run an input example provided in libgrvy's documentation on Lonestar6:
```shell
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$WORK/soft/bin/grvy/lib
g++ -I $WORK/soft/bin/grvy/include main.cpp -L$WORK/soft/bin/grvy/lib -lgrvy -o test_parsing
./test_parsing
```

where `/usr/bin` will be replaced with wherever grvy was installed. On ls6, environment variables for the `include` and `lib` paths might already exist.

### GSL Usage on TACC

To load GSL on TACC and check if it loaded, type

```shell
module load gsl
module list
```

To properly link and compile code that uses gsl, the location of the GSL distribution, libraries, include files, and tools are required. They are `TACC_GSL_DIR`, `TACC_GSL_LIB`, `TACC_GSL_INC` and `TACC_GSL_BIN`, respectively. Including these paths in your Makefile would look something like `$(TACC_GSL_LIB)`

### Documentation
We will have a model document, but we will also use docstrings to properly comment and document our code. If time permits, we will use doxygen/sphinx in this effort.

To install a C++ docstring template in VSCode, type `ctrl+p` and run/paste the following line: `ext install cschlosser.doxdocgen`. See [here](https://marketplace.visualstudio.com/items?itemName=cschlosser.doxdocgen) for more information on how this works.
