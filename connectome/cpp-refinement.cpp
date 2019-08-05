/* c++ file for running skeleton refinement algorithm */

#include "cpp-MinBinaryHeap.h"
#include "cpp-wiring.h"



struct DijkstraData {
    long iv;
    DijkstraData *prev;
    double distance;
    bool visited;
};



void CppSkeletonRefinement(const char *prefix, long label, double resolution[3])
{
    // start timing statistics
    clock_t start_time = clock();

    // clear the global variables
    segment = std::unordered_map<long, char>();
    synapses = std::unordered_set<long>();
    std::unordered_map<long, long> dijkstra_map = std::unordered_map<long, long>();

    // populate the point clouds with segment voxels and anchor points
    CppPopulatePointCloud(prefix, "skeletons", label);
    CppPopulatePointCloud(prefix, "synapses", label);
    CppPopulatePointCloud(prefix, "volumetric_somae/surfaces", label);
    
    // get the number of elements in the skeleton
    long nelements = segment.size();
    
    DijkstraData *voxel_data = new DijkstraData[nelements];
    if (!voxel_data) exit(-1);
    
    // initialize the priority queue
    DijkstraData tmp;
    MinBinaryHeap<DijkstraData *> voxel_heap(&tmp, (&tmp.distance), segment.size());

    // initialize all data
    long index = 0;
    for (std::unordered_map<long, char>::iterator it = segment.begin(); it != segment.end(); ++it, ++index) {
        voxel_data[index].iv = it->first;
        voxel_data[index].prev = NULL;
        voxel_data[index].distance = infinity;
        voxel_data[index].visited = false;
        dijkstra_map[it->first] = index;

        // this is the soma
        if (it->second == 4) {
            // insert the source into the heap
            voxel_data[index].distance = 0.0;
            voxel_data[index].visited = true;
            voxel_heap.Insert(index, &(voxel_data[index]));
        }
    }

    
    // visit all vertices
    long voxel_index;
    while (!voxel_heap.IsEmpty()) {
        DijkstraData *current = voxel_heap.DeleteMin();
        voxel_index = current->iv;

        // visit all 26 neighbors of this index
        long ix, iy, iz;
        IndexToIndices(voxel_index, ix, iy, iz);

        for (long iw = iz - 1; iw <= iz + 1; ++iw) {
            for (long iv = iy - 1; iv <= iy + 1; ++iv) {
                for (long iu = ix - 1; iu <= ix + 1; ++iu) {
                    // get the linear index for this voxel
                    long neighbor_index = IndicesToIndex(iu, iv, iw);

                    // skip if background
                    if (!segment[neighbor_index]) continue;

                    // get the corresponding neighbor data
                    long dijkstra_index = dijkstra_map[neighbor_index];
                    DijkstraData *neighbor_data = &(voxel_data[dijkstra_index]);

                    // find the distance between these voxels
                    long deltaz = resolution[OR_Z] * (iw - iz);
                    long deltay = resolution[OR_Y] * (iv - iy);
                    long deltax = resolution[OR_X] * (iu - ix);

                    // get the distance between (ix, iy, iz) and (iu, iv, iw)
                    double distance = sqrt(deltax * deltax + deltay * deltay + deltaz * deltaz);

                    // get the distance to get to this voxel through the current voxel (requires a penalty for visiting this voxel)
                    double distance_through_current = current->distance + distance;
                    double distance_without_current = neighbor_data->distance;

                    if (!neighbor_data->visited) {
                        neighbor_data->prev = current;
                        neighbor_data->distance = distance_through_current;
                        neighbor_data->visited = true;
                        voxel_heap.Insert(dijkstra_index, neighbor_data);
                    }
                    else if (distance_through_current < distance_without_current) {
                        neighbor_data->prev = current;
                        neighbor_data->distance = distance_through_current;
                        voxel_heap.DecreaseKey(dijkstra_index, neighbor_data);
                    }
                }
            }
        }
    }
  
    std::unordered_set<long> wiring_diagram = std::unordered_set<long>();

    // go through all of the synapses and add all of the skeleton points to the source
    for (std::unordered_set<long>::iterator it = synapses.begin(); it != synapses.end(); ++it) {
        // get the voxel and corresponding entry in the dijkstra data frame
        long voxel_index = *it;
        long dijkstra_index = dijkstra_map[voxel_index];

        DijkstraData *data = &(voxel_data[dijkstra_index]);

        while (data != NULL) {
            // add to the list of skeleton points
            long iv = data->iv;

            // convert to unpadded coordinates
            long ix, iy, iz;
            IndexToIndices(iv, ix, iy, iz);
            // unpad x, y, and z
            ix -= 1; iy -= 1; iz -= 1;
            // reconvert to linear coordinates
            iv = iz * (grid_size[OR_Y] - 2) * (grid_size[OR_X] - 2) + iy * (grid_size[OR_X] - 2) + ix;

            wiring_diagram.insert(iv);

            data = data->prev;
        }
    }
    
    char wiring_filename[4096];
    sprintf(wiring_filename, "connectomes/%s/%06ld.pts", prefix, label);
    char distance_filename[4096];
    sprintf(distance_filename, "distances/%s/%06ld.pts", prefix, label);

    FILE *wfp = fopen(wiring_filename, "wb"); 
    if (!wfp) { fprintf(stderr, "Failed to write to %s.\n", wiring_filename); exit(-1); }
    FILE *dfp = fopen(distance_filename, "wb");
    if (!dfp) { fprintf(stderr, "Failed to write to %s.\n", distance_filename); exit(-1); }

    // remove padding for file write
    grid_size[OR_Z] -= 2;
    grid_size[OR_Y] -= 2;
    grid_size[OR_X] -= 2;

    long nskeleton_points = wiring_diagram.size();    
    if (fwrite(&(grid_size[OR_Z]), sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s.\n", wiring_filename); exit(-1); }
    if (fwrite(&(grid_size[OR_Y]), sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s.\n", wiring_filename); exit(-1); }
    if (fwrite(&(grid_size[OR_X]), sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s.\n", wiring_filename); exit(-1); }
    if (fwrite(&nskeleton_points, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s.\n", wiring_filename); exit(-1); }

    if (fwrite(&(grid_size[OR_Z]), sizeof(long), 1, dfp) != 1) { fprintf(stderr, "Failed to write to %s.\n", distance_filename); exit(-1); }
    if (fwrite(&(grid_size[OR_Y]), sizeof(long), 1, dfp) != 1) { fprintf(stderr, "Failed to write to %s.\n", distance_filename); exit(-1); }
    if (fwrite(&(grid_size[OR_X]), sizeof(long), 1, dfp) != 1) { fprintf(stderr, "Failed to write to %s.\n", distance_filename); exit(-1); }
    if (fwrite(&nskeleton_points, sizeof(long), 1, dfp) != 1) { fprintf(stderr, "Failed to write to %s.\n", distance_filename); exit(-1); }


    for (std::unordered_set<long>::iterator it = wiring_diagram.begin(); it != wiring_diagram.end(); ++it) {
        long voxel_index = *it;
        if (fwrite(&voxel_index, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s.\n", wiring_filename); exit(-1); }

        // get the upsampled location
        long ix, iy, iz;
        iz = voxel_index / (grid_size[OR_X] * grid_size[OR_Y]);
        iy = (voxel_index - iz * grid_size[OR_X] * grid_size[OR_Y]) / grid_size[OR_X];
        ix = voxel_index % grid_size[OR_X];

        // need to repad everything
        iz += 1;
        iy += 1;
        ix += 1;

        long padded_index = iz * (grid_size[OR_X] + 2) * (grid_size[OR_Y] + 2) + iy * (grid_size[OR_X] + 2) + ix;

        // get the corresponding neighbor data
        long dijkstra_index = dijkstra_map[padded_index];
        DijkstraData *dijkstra_data = &(voxel_data[dijkstra_index]);
        double distance = dijkstra_data->distance;

        if (fwrite(&voxel_index, sizeof(long), 1, dfp) != 1) { fprintf(stderr, "Failed to write to %s.\n", wiring_filename); exit(-1); }
        if (fwrite(&distance, sizeof(double), 1, dfp) != 1) { fprintf(stderr, "Failed to write to %s.\n", wiring_filename); exit(-1); }
    }
    
    fclose(wfp);
    fclose(dfp);

    delete[] voxel_data;

    double total_time = (double) (clock() - start_time) / CLOCKS_PER_SEC;

    char time_filename[4096];
    sprintf(time_filename, "running_times/refinement/%s-%06ld.time", prefix, label);

    FILE *tfp = fopen(time_filename, "wb");
    if (!tfp) { fprintf(stderr, "Failed to write to %s.\n", time_filename); exit(-1); }

    // write the number of points and the total time to file
    if (fwrite(&nelements, sizeof(long), 1, tfp) != 1) { fprintf(stderr, "Failed to write to %s.\n", time_filename); exit(-1); }
    if (fwrite(&total_time, sizeof(double), 1, tfp) != 1) { fprintf(stderr, "Failed to write to %s.\n", time_filename); exit(-1); }

    // close file
    fclose(tfp);
}