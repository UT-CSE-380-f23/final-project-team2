#include "output.h"

// should just make this a public function; need to change how this is all set up!!
// ^ this is becuase we cannot pass the abstract class FDSolver() to a potential hdf5 class


// default constructor if no output file is specified
Output::Output(GrvyParser& grvy_parser, FDSolver& fdsolver) : Output::Output(grvy_parser, FDSolver& fdsolver, "../output/heat-eqn-output-defaultfile.h5"){};

Output::Output(GrvyParser& grvy_parser, FDSolver& fdsolver, const char* outfile) : grvy_parser(grvy_parser), outfile(outfile), fdsolver(fdsolver){};
//#define H5FILE_NAME "SDS.h5"
#define DATASETNAME "heatEqnOutput"
//#define NX          5 /* dataset dimensions */
//#define NY          6
#define RANK        this->grvy_parser.DIM
// use the 
void Output::dump()
{


    hid_t   file, dataset;       /* file and dataset handles */
    hid_t   datatype, dataspace; /* handles */
    hsize_t dimsf[this->grvy_parser.DIM];            /* dataset dimensions */
    herr_t  status;
    // this definition might have to change
    double     data[this->grvy_parser.N]; /* data to write */
    //double     data[NX][NY]; /* data to write */
    int     i, j;

    std::cout << "testing hdf5 output grvy parser with " << this->grvy_parser.N << " nodes" << std::endl;
    std::cout << "testing hdf5 output grvy parser in " << this->grvy_parser.DIM << " dimension" << std::endl;

    // not sure if this is the best way to do this; data will be changing size
    if (this->grvy_parser.DIM == 1){
        //double     data[this->grvy_parser.N];
        dimsf[0]  = this->grvy_parser.N;
        // write our data - hardcoding boundary conditions for now (ew!)
        data[0] = 0;
        for (j = 1; j < this->grvy_parser.N - 1; j++){
            std::cout << j << std::endl;
            // we'll need a function that generates the full description of u; this is just the solution
            // on interior grid points
            data[j] = gsl_vector_get(this->fdsolver.u,j);
        }
        data[this->grvy_parser.N - 1] = 0;

        // check if the size of our data matches, throw an error otherwise

    } else if (this->grvy_parser.DIM == 2){
        //double     data[this->grvy_parser.N][this->grvy_parser.N];
        dimsf[0]  = this->grvy_parser.N;
        dimsf[1]  = this->grvy_parser.N;
    } else{
        //double data;
        std::cout << "Only implemented hdf5 file generator in 1D and 2D, invalid dimension" << std::endl;
        exit(1);
    }
    /*
     * Data  and output buffer initialization.
     */
    /*
    for (j = 0; j < this->grvy_parser.N; j++)
        for (i = 0; i < NY; i++)
            data[j][i] = i + j;
    */
    //data[0][0] += 1e-6;
    /*
     * 0+eps 1 2 3 4 5
     * 1 2 3 4 5 6
     * 2 3 4 5 6 7
     * 3 4 5 6 7 8
     * 4 5 6 7 8 9
     */

    /*
     * Create a new file using H5F_ACC_TRUNC access,
     * default file creation properties, and default file
     * access properties.
     */
    //file = H5Fcreate(H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    file = H5Fcreate(outfile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /*
     * Describe the size of the array and create the data space for fixed
     * size dataset.
     */

    dataspace = H5Screate_simple(RANK, dimsf, NULL);
    //dataspace = H5Screate_scalar(RANK, dimsf, NULL);

    /*
     * Define datatype for the data in the file.
     * We will store little endian INT numbers.
     */
    datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
    status   = H5Tset_order(datatype, H5T_ORDER_LE);

    /*
     * Create a new dataset within the file using defined dataspace and
     * datatype and default dataset creation properties.
     */
    dataset = H5Dcreate2(file, DATASETNAME, datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /*
     * Write the data to the dataset using default transfer properties.
     */
    std::cout << "checking value of data " << data[0] << std::endl;
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    /*
     * Close/release resources.
     */
    H5Sclose(dataspace);
    H5Tclose(datatype);
    H5Dclose(dataset);
    H5Fclose(file);

}