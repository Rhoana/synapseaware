/* c++ file for the teaser skeletonization strategy */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unordered_set>
#include "cpp-MinBinaryHeap.h"
#include "cpp-teaser.h"




static unsigned char *segmentation = NULL;
static unsigned char *skeleton = NULL;
static float *DBF = NULL;
static float *penalties = NULL;
static float *PDRF = NULL;
static unsigned char *inside = NULL;



static long inside_voxels = 0;
static float scale = 1.3;
static long buffer = 10;



// set default values
static long grid_size[3] = { -1, -1, -1 };
static long nentries = -1;
static long sheet_size = -1;
static long row_size = -1;
static long infinity = -1;

// mask variables for bitwise operations

static long n26_offsets[26];

static void PopulateOffsets(void)
{
    n26_offsets[0] = -1 * grid_size[OR_Y] * grid_size[OR_X] - grid_size[OR_X] - 1;
    n26_offsets[1] = -1 * grid_size[OR_Y] * grid_size[OR_X] - grid_size[OR_X];
    n26_offsets[2] = -1 * grid_size[OR_Y] * grid_size[OR_X] - grid_size[OR_X] + 1;
    n26_offsets[3] = -1 * grid_size[OR_Y] * grid_size[OR_X] - 1;
    n26_offsets[4] = -1 * grid_size[OR_Y] * grid_size[OR_X];
    n26_offsets[5] = -1 * grid_size[OR_Y] * grid_size[OR_X] + 1;
    n26_offsets[6] = -1 * grid_size[OR_Y] * grid_size[OR_X] + grid_size[OR_X] - 1;
    n26_offsets[7] = -1 * grid_size[OR_Y] * grid_size[OR_X] + grid_size[OR_X];
    n26_offsets[8] = -1 * grid_size[OR_Y] * grid_size[OR_X] + grid_size[OR_X] + 1;

    n26_offsets[9] = -1 * grid_size[OR_X] - 1;
    n26_offsets[10] = -1 * grid_size[OR_X];
    n26_offsets[11] = -1 * grid_size[OR_X] + 1;
    n26_offsets[12] = -1;
    n26_offsets[13] = +1;
    n26_offsets[14] = grid_size[OR_X] - 1;
    n26_offsets[15] = grid_size[OR_X];
    n26_offsets[16] = grid_size[OR_X] + 1;

    n26_offsets[17] = grid_size[OR_Y] * grid_size[OR_X] - grid_size[OR_X] - 1;
    n26_offsets[18] = grid_size[OR_Y] * grid_size[OR_X] - grid_size[OR_X];
    n26_offsets[19] = grid_size[OR_Y] * grid_size[OR_X] - grid_size[OR_X] + 1;
    n26_offsets[20] = grid_size[OR_Y] * grid_size[OR_X] - 1;
    n26_offsets[21] = grid_size[OR_Y] * grid_size[OR_X];
    n26_offsets[22] = grid_size[OR_Y] * grid_size[OR_X] + 1;
    n26_offsets[23] = grid_size[OR_Y] * grid_size[OR_X] + grid_size[OR_X] - 1;
    n26_offsets[24] = grid_size[OR_Y] * grid_size[OR_X] + grid_size[OR_X];
    n26_offsets[25] = grid_size[OR_Y] * grid_size[OR_X] + grid_size[OR_X] + 1;
}



//////////////////////////////////////
//// COORDINATE UTILITY FUNCTIONS ////
//////////////////////////////////////

static void IndexToIndices(long iv, long &ix, long &iy, long &iz)
{
    iz = iv / sheet_size;
    iy = (iv - iz * sheet_size) / row_size;
    ix = iv % row_size;
}



static long IndicesToIndex(long ix, long iy, long iz)
{
    return iz * sheet_size + iy * row_size + ix;
}



static void ComputeDistanceFromBoundaryField(void)
{
    // start timing statistics
    clock_t start_time = clock();

    // allocate memory for bounday map and distance transform
    long *b = new long[nentries];
    for (long iz = 0; iz < grid_size[OR_Z]; ++iz) {
        for (long iy = 0; iy < grid_size[OR_Y]; ++iy) {
            for (long ix = 0; ix < grid_size[OR_X]; ++ix) {
                if (!segmentation[IndicesToIndex(ix, iy, iz)]) {
                    b[IndicesToIndex(ix, iy, iz)] = 0;
                    continue;
                }

                // inside voxels that are on the boundary have value 1 (based on TEASER paper figure 2)
                if ((ix == 0 or iy == 0 or iz == 0 or (ix == grid_size[OR_X] - 1) or (iy == grid_size[OR_Y] - 1) or (iz == grid_size[OR_Z] - 1)) ||
                    (ix > 0 and !segmentation[IndicesToIndex(ix - 1, iy, iz)]) ||
                    (iy > 0 and !segmentation[IndicesToIndex(ix, iy - 1, iz)]) ||
                    (iz > 0 and !segmentation[IndicesToIndex(ix, iy, iz - 1)]) ||
                    (ix < grid_size[OR_X] - 1 and !segmentation[IndicesToIndex(ix + 1, iy, iz)]) ||
                    (iy < grid_size[OR_Y] - 1 and !segmentation[IndicesToIndex(ix, iy + 1, iz)]) ||
                    (iz < grid_size[OR_Z] - 1 and !segmentation[IndicesToIndex(ix, iy, iz + 1)])) {
                    b[IndicesToIndex(ix, iy, iz)] = 1;
                }
                else {
                    b[IndicesToIndex(ix, iy, iz)] = infinity;
                }
            }
        }
    }

    // go along the z dimenion first for every (x, y) coordinate
    for (long ix = 0; ix < grid_size[OR_X]; ++ix) {
        for (long iy = 0; iy < grid_size[OR_Y]; ++iy) {

            long k = 0;
            long *v = new long[grid_size[OR_Z] + 1];
            float *z = new float[grid_size[OR_Z] + 1];

            v[0] = 0;
            z[0] = -1 * infinity;
            z[1] = infinity;

            for (long q = 1; q < grid_size[OR_Z]; ++q) {
                // label for jump statement
                zlabel:
                float s = ((b[IndicesToIndex(ix, iy, q)] + q * q) - (b[IndicesToIndex(ix, iy, v[k])] + v[k] * v[k])) / (float)(2 * q - 2 * v[k]);
                
                if (s <= z[k]) {
                    --k;
                    goto zlabel;
                }
                else {
                    ++k;
                    v[k] = q;
                    z[k] = s;
                    z[k + 1] = infinity;
                }
            }

            k = 0;
            for (long q = 0; q < grid_size[OR_Z]; ++q) {
                while (z[k + 1] < q)
                    ++k;

                DBF[IndicesToIndex(ix, iy, q)] = (q - v[k]) * (q - v[k]) + b[IndicesToIndex(ix, iy, v[k])];
            }

            // free memory 
            delete[] v;
            delete[] z;
        }
    }

    // update the boundary values with this distance
    for (long iz = 0; iz < grid_size[OR_Z]; ++iz) {
        for (long iy = 0; iy < grid_size[OR_Y]; ++iy) {
            for (long ix = 0; ix < grid_size[OR_X]; ++ix) {
                b[IndicesToIndex(ix, iy, iz)] = DBF[IndicesToIndex(ix, iy, iz)];
            }
        }
    }

    // go along the y dimension second for every (z, x) coordinate
    for (long iz = 0; iz < grid_size[OR_Z]; ++iz) {
        for (long ix = 0; ix < grid_size[OR_X]; ++ix) {

            long k = 0;
            long *v = new long[grid_size[OR_Y] + 1];
            float *z = new float[grid_size[OR_Y] + 1];

            v[0] = 0;
            z[0] = -1 * infinity;
            z[1] = infinity;

            for (long q = 1; q < grid_size[OR_Y]; ++q) {
                // label for jump statement
                ylabel:
                float s = ((b[IndicesToIndex(ix, q, iz)] + q * q) - (b[IndicesToIndex(ix, v[k], iz)] +  v[k] * v[k])) / (float)(2 * q - 2 * v[k]);
                
                if (s <= z[k]) {
                    --k;
                    goto ylabel;
                }
                else {
                    ++k; 
                    v[k] = q;
                    z[k] = s;
                    z[k + 1] = infinity;
                }
            }

            k = 0;
            for (long q = 0; q < grid_size[OR_Y]; ++q) {
                while (z[k + 1] < q)
                    ++k;
            
                DBF[IndicesToIndex(ix, q, iz)] = (q - v[k]) * (q - v[k]) + b[IndicesToIndex(ix, v[k], iz)];
            }

            // free memory
            delete[] v;
            delete[] z;
        }
    }

    // update the boundary values with this distance
    for (long iz = 0; iz < grid_size[OR_Z]; ++iz) {
        for (long iy = 0; iy < grid_size[OR_Y]; ++iy) {
            for (long ix = 0; ix < grid_size[OR_X]; ++ix) {
                b[IndicesToIndex(ix, iy, iz)] = DBF[IndicesToIndex(ix, iy, iz)];
            }
        }
    }


    // go along the x dimension last for every (y, z) coordinate
    for (long iy = 0; iy < grid_size[OR_Y]; ++iy) {
        for (long iz = 0; iz < grid_size[OR_Z]; ++iz) {

            long k = 0;
            long *v = new long[grid_size[OR_X] + 1];
            float *z = new float[grid_size[OR_X] + 1];

            v[0] = 0;
            z[0] = -1 * infinity;
            z[1] = infinity;

            for (long q = 1; q < grid_size[OR_X]; ++q) {
                // label for jump statement
                xlabel:
                float s = ((b[IndicesToIndex(q, iy, iz)] + q * q) - (b[IndicesToIndex(v[k], iy, iz)] + v[k] * v[k])) / (float)(2 * q - 2 * v[k]);

                if (s <= z[k]) {
                    --k;
                    goto xlabel;
                }
                else {
                    ++k;
                    v[k] = q;
                    z[k] = s;
                    z[k + 1] = infinity;
                }
            }

            k = 0;
            for (long q = 0;  q < grid_size[OR_X]; ++q) {
                while (z[k + 1] < q)
                    ++k;

                DBF[IndicesToIndex(q, iy, iz)] = (q - v[k]) * (q - v[k]) + b[IndicesToIndex(v[k], iy, iz)];
            }

            // free memory
            delete[] v;
            delete[] z;
        }
    }

    for (long iv = 0; iv < nentries; ++iv) {
        DBF[iv] = sqrt(DBF[iv]);
    }

    double total_time = (double) (clock() - start_time) / CLOCKS_PER_SEC;
    printf("Computed distance from boundary field in %0.2f seconds.\n", total_time);

    // free memory
    delete[] b;
}



struct DijkstraData {
    long iv;
    DijkstraData *prev;
    float voxel_penalty;
    float distance;
    bool visited;
};



long ComputeDistanceFromVoxelField(long source_index)
{
    DijkstraData *voxel_data = new DijkstraData[nentries];
    if (!voxel_data) exit(-1);

    // initialize all data
    for (int iv = 0; iv < nentries; ++iv) {
        voxel_data[iv].iv = iv;
        voxel_data[iv].prev = NULL;
        voxel_data[iv].voxel_penalty = penalties[iv];
        voxel_data[iv].distance = infinity;
        voxel_data[iv].visited = false;
    }

    // initialize the priority queue
    DijkstraData tmp;
    MinBinaryHeap<DijkstraData *> voxel_heap(&(tmp), (&tmp.distance), nentries);

    // insert the source into the heap
    voxel_data[source_index].distance = 0.0;
    voxel_data[source_index].visited = true;
    voxel_heap.Insert(source_index, &(voxel_data[source_index]));

    // visit all vertices
    long voxel_index = 0;
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
                    if (!segmentation[neighbor_index]) continue;
                    
                    // get the corresponding neighbor data
                    DijkstraData *neighbor_data = &(voxel_data[neighbor_index]);

                    // find the distance between these voxels
                    long deltaz = (iw - iz);
                    long deltay = (iv - iy);
                    long deltax = (iu - ix);

                    // get the distance between (ix, iy, iz) and (iu, iv, iw)
                    float distance = sqrt(deltax * deltax + deltay * deltay + deltaz * deltaz);

                    // get the distance to get to this voxel through the current voxel (requires a penalty for visiting this voxel)
                    float distance_through_current = current->distance + distance + neighbor_data->voxel_penalty;
                    float distance_without_current = neighbor_data->distance;

                    if (!neighbor_data->visited) {
                        neighbor_data->prev = current;
                        neighbor_data->distance = distance_through_current;
                        neighbor_data->visited = true;
                        voxel_heap.Insert(neighbor_index, neighbor_data);
                    }
                    else if (distance_through_current < distance_without_current) {
                        neighbor_data->prev = current;
                        neighbor_data->distance = distance_through_current;
                        voxel_heap.DecreaseKey(neighbor_index, neighbor_data);
                    }
                }
            }
        }
    }

    // first call to this function needs to return the root and does not compute the skeleton
    if (!PDRF) {
        // free memory
        delete[] voxel_data;

        // return the farthest voxel (to get the root voxel)
        return voxel_index;
    }

    // save the PDRF (only called when given root voxel)
    for (long iv = 0; iv < nentries; ++iv) {
        if (!segmentation[iv]) continue;
        PDRF[iv] = voxel_data[iv].distance;
    }

    // continue until there are no more inside voxels
    while (inside_voxels) {
        printf("  Inside voxels remaining: %ld\n", inside_voxels);
        float farthest_pdrf = -1;
        long starting_voxel = -1;

        // find the farthest PDRF that is still inside
        for (long iv = 0; iv < nentries; ++iv) {
            if (!inside[iv]) continue;
            if (PDRF[iv] > farthest_pdrf) {
                farthest_pdrf = PDRF[iv];
                starting_voxel = iv;
            }
        }

        for (long iv = 0; iv < nentries; ++iv) {
            if (!inside[iv]) continue;
            long ix, iy, iz;
            IndexToIndices(iv, ix, iy, iz);

            // get the skeleton path from this location to the root
            DijkstraData *current = &(voxel_data[starting_voxel]);

            while (!skeleton[current->iv]) {
                long ii, ij, ik;
                IndexToIndices(current->iv, ii, ij, ik);
                // what is the distance between this skeleton location and the inside location
                float deltax = (ii - ix);
                float deltay = (ij - iy);
                float deltaz = (ik - iz);

                float distance = sqrt(deltax * deltax + deltay * deltay + deltaz * deltaz);

                if (distance < scale * DBF[current->iv] + buffer) {
                    inside[iv] = 0;
                    inside_voxels--;
                    break;
                }

                // update skeleton pointer
                current = current->prev;
            }
        }

        DijkstraData *current = &(voxel_data[starting_voxel]);
        while (!skeleton[current->iv]) {
            skeleton[current->iv] = 1;
            current = current->prev;
        }
    }

    // free memory
    delete[] voxel_data;

    return -1;
}



void ComputePenalties(void)
{
    // get the maximum distance from the boundary
    float M = 0;
    for (long iv = 0; iv < nentries; ++iv) {
        if (DBF[iv] > M) M = DBF[iv]; 
    }

    // choose 5000 so that 3000 length voxel paths have correct floating point precision
    const float pdrf_scale = 5000;
    for (long iv = 0; iv < nentries; ++iv) {
        penalties[iv] = pdrf_scale * pow(1 - DBF[iv] / M, 16);
    }
}



static bool IsEndpoint(long iv)
{
    long ix, iy, iz;
    IndexToIndices(iv, ix, iy, iz);

    short nnneighbors = 0;
    for (long iw = iz - 1; iw <= iz + 1; ++iw) {
        for (long iv = iy - 1; iv <= iy + 1; ++iv) {
            for (long iu = ix - 1; iu <= ix + 1; ++iu) {
                long linear_index = IndicesToIndex(iu, iv, iw);
                if (segmentation[linear_index]) nnneighbors++;
            }
        }
    }

    // return if there is one neighbor (other than iv) that is 1
    if (nnneighbors <= 2) return true;
    else return false;
}



static int maximum(long a, long b, long c) {
    int max = (a < b) ? b : a;
    return ((max < c) ? c : max);
}



void CppTeaserSkeletonization(const char *prefix, long label, float resolution[3])
{
    // start timing statistics
    clock_t start_time = clock();

    // get the downsample factor
    float downsample_factor[3] = { 100 / resolution[OR_Z], 100 / resolution[OR_Y], 100 / resolution[OR_X] };

    // read in the point cloud for this label
    char filename[4096];
    sprintf(filename, "segmentations/%s/%06ld.pts", prefix, label);

    FILE *fp = fopen(filename, "rb");
    if (!fp) { fprintf(stderr, "Failed to read %s.\n", filename); exit(-1); }

    long input_grid_size[3];
    long npoints;
    if (fread(&(input_grid_size[OR_Z]), sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", filename); exit(-1); }
    if (fread(&(input_grid_size[OR_Y]), sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", filename); exit(-1); }
    if (fread(&(input_grid_size[OR_X]), sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", filename); exit(-1); }
    if (fread(&npoints, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", filename); exit(-1); }

    grid_size[OR_Z] = input_grid_size[OR_Z] / downsample_factor[OR_Z];
    grid_size[OR_Y] = input_grid_size[OR_Y] / downsample_factor[OR_Y];
    grid_size[OR_X] = input_grid_size[OR_X] / downsample_factor[OR_X];

    // pad the grid size by 2
    grid_size[OR_Z] += 2;
    grid_size[OR_Y] += 2;
    grid_size[OR_X] += 2;

    // set global indexing parameters
    nentries = grid_size[OR_Z] * grid_size[OR_Y] * grid_size[OR_X];
    sheet_size = grid_size[OR_Y] * grid_size[OR_X];
    row_size = grid_size[OR_X];
    infinity = grid_size[OR_Z] * grid_size[OR_Z] + grid_size[OR_Y] * grid_size[OR_Y] + grid_size[OR_X] * grid_size[OR_X];

    PopulateOffsets();

    // allocate memory for global variables
    segmentation = new unsigned char[nentries];
    skeleton = new unsigned char[nentries];
    penalties = new float[nentries];
    inside = new unsigned char[nentries];
    DBF = new float[nentries];
    for (long iv = 0; iv < nentries; ++iv) {
        segmentation[iv] = 0;
        skeleton[iv] = 0;
        penalties[iv] = 0;
        inside[iv] = 0;
        DBF[iv] = 0;
    }

    inside_voxels = 0;
    for (long ip = 0; ip < npoints; ++ip) {
        long voxel_index;
        if (fread(&voxel_index, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", filename); exit(-1); }

        long iz = voxel_index / (input_grid_size[OR_Y] * input_grid_size[OR_X]);
        long iy = (voxel_index - iz * (input_grid_size[OR_Y] * input_grid_size[OR_X])) / input_grid_size[OR_X];
        long ix = voxel_index % input_grid_size[OR_X];

        // get the location in the new grid
        iz = (int) (iz / downsample_factor[OR_Z]);
        iy = (int) (iy / downsample_factor[OR_Y]);
        ix = (int) (ix / downsample_factor[OR_X]);

        // pad the location by one
        iz += 1; iy += 1; ix += 1;

        long iv = IndicesToIndex(ix, iy, iz);
        // skip if already marked
        if (segmentation[iv]) continue;

        segmentation[iv] = 1;
        inside[iv] = 1;
        inside_voxels++;
    }
    printf("Input Size for label %ld: %ld\n", label, inside_voxels);
    fclose(fp);

    ComputeDistanceFromBoundaryField();

    // set any voxel as the source
    long source_voxel = -1;
    for (long iv = 0; iv < nentries; ++iv)
        if (inside[iv]) { source_voxel = iv; break; }

    // find a root voxel which is guaranteed to be at an extrema point
    long root_voxel = ComputeDistanceFromVoxelField(source_voxel);
    skeleton[root_voxel] = 1;
    inside[root_voxel] = 0;
    inside_voxels--;

    ComputePenalties();
    PDRF = new float[nentries];
    ComputeDistanceFromVoxelField(root_voxel);

    long num = 0;
    for (long iv = 0; iv < nentries; ++iv) {
        if (skeleton[iv]) num++;
    }

    // delete segmentation and point to skeleton to find which indices are endpoints
    delete[] segmentation;
    segmentation = skeleton;

    std::unordered_set<long> upsampled_skeleton = std::unordered_set<long>();
    for (long is = 0; is < nentries; ++is) {
        if (!skeleton[is]) continue;
        bool endpoint = false;
        if (IsEndpoint(is)) endpoint = true;
        else endpoint = false;

        // get the upsampled value
        long ix, iy, iz;
        IndexToIndices(is, ix, iy, iz);

        // unpad the coordinates
        --ix; --iy; --iz;

        // upsample to full resolution
        ix = (int)(ix * downsample_factor[OR_X]);
        iy = (int)(iy * downsample_factor[OR_Y]);
        iz = (int)(iz * downsample_factor[OR_Z]);

        long upsampled_index = iz * input_grid_size[OR_Y] * input_grid_size[OR_X] + iy * input_grid_size[OR_X] + ix;
        if (endpoint) upsampled_skeleton.insert(-1 * upsampled_index);
        else upsampled_skeleton.insert(upsampled_index);

        // look at all 26 neighbors of iv
        for (long ip = 0; ip < 26; ++ip) {
            long neighbor_index = is + n26_offsets[ip];
            if (!skeleton[neighbor_index]) continue;

            long iu, iv, iw;
            IndexToIndices(neighbor_index, iu, iv, iw);

            // unpad the coordinates
            --iu; --iv; --iw;

            // upsample to full resolution
            iu = (int)(iu * downsample_factor[OR_X]);
            iv = (int)(iv * downsample_factor[OR_Y]);
            iw = (int)(iw * downsample_factor[OR_Z]);

            while ((ix != iu) or (iy != iv) or (iz != iw)) {
                long xdiff = abs(iu - ix);
                long ydiff = abs(iv - iy);
                long zdiff = abs(iw - iz);
                long max_diff = maximum(xdiff, ydiff, zdiff);

                if (xdiff == max_diff) {
                    if (iu > ix) iu--;
                    else iu++;
                }
                if (ydiff == max_diff) {
                    if (iv > iy) iv--;
                    else iv++;
                }
                if (zdiff == max_diff) {
                    if (iw > iz) iw--;
                    else iw++;
                }

                long upsampled_path_index = iw * input_grid_size[OR_Y] * input_grid_size[OR_X] + iv * input_grid_size[OR_X] + iu;
                upsampled_skeleton.insert(upsampled_path_index);
            }
        }
    }

    char output_filename[4096];
    sprintf(output_filename, "baselines/teasers/%s/%06d.pts", prefix, label);

    FILE *wfp = fopen(output_filename, "wb");
    if (!wfp) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }

    long nskeleton_points = upsampled_skeleton.size();
    if (fwrite(&(input_grid_size[OR_Z]), sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
    if (fwrite(&(input_grid_size[OR_Y]), sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
    if (fwrite(&(input_grid_size[OR_X]), sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
    if (fwrite(&nskeleton_points, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
    
    printf("Skeleton size: %ld\n", nskeleton_points);
    for (std::unordered_set<long>::iterator it = upsampled_skeleton.begin(); it != upsampled_skeleton.end(); ++it) {
        long index = *it;
        if (fwrite(&index, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
    }

    // close file
    fclose(wfp);

    // free memory
    delete[] skeleton;
    delete[] penalties;
    delete[] PDRF;
    delete[] inside;
    delete[] DBF;

    // reset global variables
    segmentation = NULL;
    skeleton = NULL;
    penalties = NULL;
    PDRF = NULL;
    inside = NULL;
    DBF = NULL;

    double total_time = (double) (clock() - start_time) / CLOCKS_PER_SEC;

    char time_filename[4096];
    sprintf(time_filename, "running_times/teasers/%s-%06ld.time", prefix, label);

    FILE *tfp = fopen(time_filename, "wb");
    if (!tfp) { fprintf(stderr, "Failed to write to %s.\n", time_filename); exit(-1); }

    // write the number of points and the total time to file
    if (fwrite(&npoints, sizeof(long), 1, tfp) != 1) { fprintf(stderr, "Failed to write to %s.\n", time_filename); exit(-1); }
    if (fwrite(&total_time, sizeof(double), 1, tfp) != 1) { fprintf(stderr, "Failed to write to %s.\n", time_filename); exit(-1); }

    // close file
    fclose(tfp);
}
