/* c++ file for mutual functions and variables */

#include <string.h>
#include "cpp-wiring.h"



// set default values
long grid_size[3] = { -1, -1, -1 };
long nentries = -1;
long sheet_size = -1;
long row_size = -1;
long infinity = -1;

// create new default variables (will be overwritten)
std::unordered_map<long, char> segment = std::unordered_map<long, char>();
std::unordered_set<long> synapses = std::unordered_set<long>();




///////////////////////////////////////
//// POINT CLOUD UTILITY FUNCTIONS ////
///////////////////////////////////////

/* conventient I/O function */
void CppPopulatePointCloud(const char *prefix, const char *dataset, long label) {
    // read in the point cloud for this label
    char filename[4096];
    sprintf(filename, "%s/%s/%06ld.pts", dataset, prefix, label);

    FILE *fp = fopen(filename, "rb");
    if (!fp) { fprintf(stderr, "Failed to read %s.\n", filename); exit(-1); }

    long input_grid_size[3];
    long npoints;
    if (fread(&(input_grid_size[OR_Z]), sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", filename); exit(-1); }
    if (fread(&(input_grid_size[OR_Y]), sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", filename); exit(-1); }
    if (fread(&(input_grid_size[OR_X]), sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", filename); exit(-1); }
    if (fread(&npoints, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", filename); exit(-1); }

    // add padding around each segment (only way that populate offsets works!!)
    grid_size[OR_Z] = input_grid_size[OR_Z] + 2;
    grid_size[OR_Y] = input_grid_size[OR_Y] + 2;
    grid_size[OR_X] = input_grid_size[OR_X] + 2;
    
    // set global indexing parameters (do here since for loop calls IndicesToIndex)
    nentries = grid_size[OR_Z] * grid_size[OR_Y] * grid_size[OR_X];
    sheet_size = grid_size[OR_Y] * grid_size[OR_X];
    row_size = grid_size[OR_X];
    infinity = grid_size[OR_Z] * grid_size[OR_Z] + grid_size[OR_Y] * grid_size[OR_Y] + grid_size[OR_X] * grid_size[OR_X];

    for (long ip = 0; ip < npoints; ++ip) {
        long voxel_index;
        if (fread(&voxel_index, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", filename); exit(-1); }

        long iz = voxel_index / (input_grid_size[OR_Y] * input_grid_size[OR_X]);
        long iy = (voxel_index - iz * (input_grid_size[OR_Y] * input_grid_size[OR_X])) / input_grid_size[OR_X];
        long ix = voxel_index % input_grid_size[OR_X];

        //  pad the location by one
        iz += 1; iy += 1; ix += 1;

        // find the new voxel index
        long iv = IndicesToIndex(ix, iy, iz);

        if (!strcmp(dataset, "segmentations")) {
            segment[iv] = 1;
        }
        else if (!strcmp(dataset, "skeletons")) {
            segment[iv] = 1;
        }
        else if (!strcmp(dataset, "synapses")) {
            segment[iv] = 3;
            synapses.insert(iv);
        }
        else if (!strcmp(dataset, "somae")) {
            segment[iv] = 4;
        }
        else if (!strcmp(dataset, "volumetric_somae/surfaces")) {
            segment[iv] = 4;
        }
        else { fprintf(stderr, "Unrecognized point cloud: %s.\n", dataset); exit(-1); }
    }

    // close file
    fclose(fp);
}
