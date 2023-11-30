#include "output.h"


// default constructor if no output file is specified
Output::Output(GrvyParser& grvy_parser) : Output::Output(grvy_parser, "../output/heat-eqn-output-defaultfile.h5"){};

Output::Output(GrvyParser& grvy_parser, const char* outfile) : grvy_parser(grvy_parser), outfile(outfile){};
//#define H5FILE_NAME "SDS.h5"
#define DATASETNAME "IntArray"
#define NX          5 /* dataset dimensions */
#define NY          6
#define RANK        2
// use the 
void Output::dump()
{


    hid_t   file, dataset;       /* file and dataset handles */
    hid_t   datatype, dataspace; /* handles */
    hsize_t dimsf[2];            /* dataset dimensions */
    herr_t  status;
    double     data[NX][NY]; /* data to write */
    int     i, j;

    std::cout << "testing hdf5 output grvy parser with " << this->grvy_parser.N << " nodes" << std::endl;

    /*
     * Data  and output buffer initialization.
     */
    for (j = 0; j < NX; j++)
        for (i = 0; i < NY; i++)
            data[j][i] = i + j;

    data[0][0] += 1e-6;
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
    dimsf[0]  = NX;
    dimsf[1]  = NY;
    dataspace = H5Screate_simple(RANK, dimsf, NULL);

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
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    /*
     * Close/release resources.
     */
    H5Sclose(dataspace);
    H5Tclose(datatype);
    H5Dclose(dataset);
    H5Fclose(file);

}