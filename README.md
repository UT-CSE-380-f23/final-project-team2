# final-project-team2
final-project-team2 created by GitHub Classroom
Authors:  Rodrigo Gonzalez, Kenneth L  Meyer, Jenil Shah

## Installation Notes
todo: check if the code should be able to be built and ran on *any* system, or just lonestar6. I.e. does the `Makefile` just need to work properly for any user on ls6? Or for users beyond ls6 as well?

### Running a libgrvy application:
The following commands were used to successfully run an input example provided in libgrvy's documentation:
```shell
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/bin/grvy/lib
g++ -I /usr/bin/grvy/include cpp-input.cpp -L/usr/bin/grvy/lib -lgrvy -o cpp_run_grvy_test
```

where `/usr/bin` will be replaced with wherever grvy was installed. On ls6, environment variables for the `include` and `lib` paths might already exist.

### Documentation
We will have a model document, but we will also use docstrings to properly comment and document our code. If time permits, we will use doxygen/sphinx in this effort.

To install a C++ docstring template in VSCode, type `ctrl+p` and run/paste the following line: `ext install cschlosser.doxdocgen`. See [here](https://marketplace.visualstudio.com/items?itemName=cschlosser.doxdocgen) for more information on how this works.
