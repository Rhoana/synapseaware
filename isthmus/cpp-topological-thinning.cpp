/* c++ file to generate skeletons with thinning methods */

#include <stdio.h>
#include <stdlib.h>
#include <unordered_map>
#include <unordered_set>
#include "cpp-topological-thinning.h"



// constant variables

static const int lookup_table_size = 1 << 23;
static const int NTHINNING_DIRECTIONS = 6;
static const int UP = 0;
static const int DOWN = 1;
static const int NORTH = 2;
static const int SOUTH = 3;
static const int EAST = 4;
static const int WEST = 5;



// lookup tables

static unsigned char *lut_simple;
static unsigned char *lut_isthmus;

// set default values
static long grid_size[3] = { -1, -1, -1 };
static long nentries = -1;
static long sheet_size = -1;
static long row_size = -1;
static long infinity = -1;
static std::unordered_map<long, char> segment;


// mask variables for bitwise operations

static long long_mask[26];
static unsigned char char_mask[8];
static long n26_offsets[26];
static long n6_offsets[6];


static void set_long_mask(void)
{
    long_mask[ 0] = 0x00000001;
    long_mask[ 1] = 0x00000002;
    long_mask[ 2] = 0x00000004;
    long_mask[ 3] = 0x00000008;
    long_mask[ 4] = 0x00000010;
    long_mask[ 5] = 0x00000020;
    long_mask[ 6] = 0x00000040;
    long_mask[ 7] = 0x00000080;
    long_mask[ 8] = 0x00000100;
    long_mask[ 9] = 0x00000200;
    long_mask[10] = 0x00000400;
    long_mask[11] = 0x00000800;
    long_mask[12] = 0x00001000;
    long_mask[13] = 0x00002000;
    long_mask[14] = 0x00004000;
    long_mask[15] = 0x00008000;
    long_mask[16] = 0x00010000;
    long_mask[17] = 0x00020000;
    long_mask[18] = 0x00040000;
    long_mask[19] = 0x00080000;
    long_mask[20] = 0x00100000;
    long_mask[21] = 0x00200000;
    long_mask[22] = 0x00400000;
    long_mask[23] = 0x00800000;
    long_mask[24] = 0x01000000;
    long_mask[25] = 0x02000000;
}

static void set_char_mask(void)
{
    char_mask[0] = 0x01;
    char_mask[1] = 0x02;
    char_mask[2] = 0x04;
    char_mask[3] = 0x08;
    char_mask[4] = 0x10;
    char_mask[5] = 0x20;
    char_mask[6] = 0x40;
    char_mask[7] = 0x80;
}

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

    // use this order to go UP, DOWN, NORTH, SOUTH, EAST, WEST
    // DO NOT CHANGE THIS ORDERING
    n6_offsets[0] = -1 * grid_size[OR_X];
    n6_offsets[1] = grid_size[OR_X];
    n6_offsets[2] = -1 * grid_size[OR_Y] * grid_size[OR_X];
    n6_offsets[3] = grid_size[OR_Y] * grid_size[OR_X];
    n6_offsets[4] = +1;
    n6_offsets[5] = -1;
}


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




// very simple double linked list data structure

typedef struct {
    long iv, ix, iy, iz;
    void *next;
    void *prev;
} ListElement;

typedef struct {
    void *first;
    void *last;
} List;

typedef struct {
    long iv, ix, iy, iz;
} Voxel;

typedef struct {
    Voxel v;
    ListElement *ptr;
    void *next;
} Cell;

typedef struct {
    Cell *head;
    Cell *tail;
    int length;
} PointList;

typedef struct {
    ListElement *first;
    ListElement *last;
} DoubleList;

List surface_voxels;



static void NewSurfaceVoxel(long iv, long ix, long iy, long iz)
{
    ListElement *LE = new ListElement();
    LE->iv = iv;
    LE->ix = ix;
    LE->iy = iy;
    LE->iz = iz;

    LE->next = NULL;
    LE->prev = surface_voxels.last;

    if (surface_voxels.last != NULL) ((ListElement *) surface_voxels.last)->next = LE;
    surface_voxels.last = LE;
    if (surface_voxels.first == NULL) surface_voxels.first = LE;
}



static void RemoveSurfaceVoxel(ListElement *LE) 
{
    ListElement *LE2;
    if (surface_voxels.first == LE) surface_voxels.first = LE->next;
    if (surface_voxels.last == LE) surface_voxels.last = LE->prev;

    if (LE->next != NULL) {
        LE2 = (ListElement *)(LE->next);
        LE2->prev = LE->prev;
    }
    if (LE->prev != NULL) {
        LE2 = (ListElement *)(LE->prev);
        LE2->next = LE->next;
    }
    delete LE;
}



static void CreatePointList(PointList *s) 
{
    s->head = NULL;
    s->tail = NULL;
    s->length = 0;
}



static void AddToList(PointList *s, Voxel e, ListElement *ptr) 
{
    Cell *newcell = new Cell();
    newcell->v = e;
    newcell->ptr = ptr;
    newcell->next = NULL;

    if (s->head == NULL) {
        s->head = newcell;
        s->tail = newcell;
        s->length = 1;
    }
    else {
        s->tail->next = newcell;
        s->tail = newcell;
        s->length++;
    }
}



static Voxel GetFromList(PointList *s, ListElement **ptr)
{
    Voxel V;
    Cell *tmp;
    V.iv = -1;
    V.ix = -1;
    V.iy = -1;
    V.iz = -1;
    (*ptr) = NULL;
    if (s->length == 0) return V;
    else {
        V = s->head->v;
        (*ptr) = s->head->ptr;
        tmp = (Cell *) s->head->next;
        delete s->head;
        s->head = tmp;
        s->length--;
        if (s->length == 0) {
            s->head = NULL;
            s->tail = NULL;
        }
        return V;
    }
}



static void DestroyPointList(PointList *s) {
    ListElement *ptr;
    while (s->length) GetFromList(s, &ptr);
}



static void InitializeLookupTables(const char *lookup_table_directory)
{
    char lut_filename[4096];
    FILE *lut_file;

    // read the simple lookup table
    sprintf(lut_filename, "%s/lut_simple.dat", lookup_table_directory);
    lut_simple = new unsigned char[lookup_table_size];
    lut_file = fopen(lut_filename, "rb");
    if (!lut_file) {
        fprintf(stderr, "Failed to read %s\n", lut_filename);
        exit(-1);
    }
    if (fread(lut_simple, 1, lookup_table_size, lut_file) != lookup_table_size) {
        fprintf(stderr, "Failed to read %s\n", lut_filename);
        exit(-1);
    }
    fclose(lut_file);

    // read the isthmus lookup table
    sprintf(lut_filename, "%s/lut_isthmus.dat", lookup_table_directory);
    lut_isthmus = new unsigned char[lookup_table_size];
    lut_file = fopen(lut_filename, "rb");
    if (!lut_file) {
        fprintf(stderr, "Failed to read %s\n", lut_filename);
        exit(-1);
    }
    if (fread(lut_isthmus, 1, lookup_table_size, lut_file) != lookup_table_size) {
        fprintf(stderr, "Failed to read %s\n", lut_filename);
        exit(-1);
    }
    fclose(lut_file);

    // set the mask variables
    set_char_mask();
    set_long_mask();
}



static void CollectSurfaceVoxels(void)
{
    // go through all voxels and check their six neighbors
    for (std::unordered_map<long, char>::iterator it = segment.begin(); it != segment.end(); ++it) {
        // all of these elements are either 1 or 3 and in the segment
        long index = it->first;

        long ix, iy, iz;
        IndexToIndices(index, ix, iy, iz);
        // check the 6 neighbors
        for (long iv = 0; iv < NTHINNING_DIRECTIONS; ++iv) {
            long neighbor_index = index + n6_offsets[iv];

            long ii, ij, ik;
            IndexToIndices(neighbor_index, ii, ij, ik);
            
            // skip the fake boundary elements
            if ((ii == 0) or (ii == grid_size[OR_X] - 1)) continue;
            if ((ij == 0) or (ij == grid_size[OR_Y] - 1)) continue;
            if ((ik == 0) or (ik == grid_size[OR_Z] - 1)) continue;
            
            if (segment.find(neighbor_index) == segment.end()) {
                // this location is a boundary so create a surface voxel and break
                it->second = 2;

                NewSurfaceVoxel(index, ix, iy, iz);

                // note this location as surface in erosions list
                break;
            }
        }
    }
}



static unsigned int Collect26Neighbors(long ix, long iy, long iz)
{
    unsigned int neighbors = 0;
    long index = IndicesToIndex(ix, iy, iz);
    
    // some of these lookups will create a new entry but the region is 
    // shrinking so memory overhead is minimal
    for (long iv = 0; iv < 26; ++iv) {
        if (segment[index + n26_offsets[iv]]) neighbors |= long_mask[iv];
    }

    return neighbors;
}



static bool Simple26_6(unsigned int neighbors)
{
    return lut_simple[(neighbors >> 3)] & char_mask[neighbors % 8];
}



static bool Isthmus(unsigned int neighbors)
{
    return lut_isthmus[(neighbors >> 3)] & char_mask[neighbors % 8];
}



static void DetectSimpleBorderPoints(PointList *deletable_points, int direction)
{
    ListElement *LE = (ListElement *)surface_voxels.first;
    while (LE != NULL) {
        long iv = LE->iv;
        long ix = LE->ix;
        long iy = LE->iy;
        long iz = LE->iz;

        // not an isthmus
        if (segment[iv] == 2) {
            long value = 0;
            switch (direction) {
            case UP: {
                value = segment[IndicesToIndex(ix, iy - 1, iz)];
                break;
            }
            case DOWN: {
                value = segment[IndicesToIndex(ix, iy + 1, iz)];
                break;
            }
            case NORTH: {
                value = segment[IndicesToIndex(ix, iy, iz - 1)];
                break;
            }
            case SOUTH: {
                value = segment[IndicesToIndex(ix, iy, iz + 1)];
                break;
            }
            case EAST: {
                value = segment[IndicesToIndex(ix + 1, iy, iz)];
                break;
            }
            case WEST: {
                value = segment[IndicesToIndex(ix - 1, iy, iz)];
                break;
            }
            }

            // see if the required point belongs to a different segment
            if (!value) {
                unsigned int neighbors = Collect26Neighbors(ix, iy, iz);

                // deletable point
                if (Simple26_6(neighbors)) {
                    Voxel voxel;
                    voxel.iv = iv;
                    voxel.ix = ix;
                    voxel.iy = iy;
                    voxel.iz = iz;
                    AddToList(deletable_points, voxel, LE);
                }
                else {
                    if (Isthmus(neighbors)) {
                        segment[iv] = 3;
                    }
                }
            }
        }
        LE = (ListElement *) LE->next;
    }
}



static long ThinningIterationStep(void)
{
    long changed = 0;

    // iterate through every direction
    for (int direction = 0; direction < NTHINNING_DIRECTIONS; ++direction) {
        PointList deletable_points;
        ListElement *ptr;

        CreatePointList(&deletable_points);
        DetectSimpleBorderPoints(&deletable_points, direction);

        while (deletable_points.length) {
            Voxel voxel = GetFromList(&deletable_points, &ptr);

            long iv = voxel.iv;
            long ix = voxel.ix;
            long iy = voxel.iy;
            long iz = voxel.iz;

            unsigned int neighbors = Collect26Neighbors(ix, iy, iz);
            if (Simple26_6(neighbors)) {
                // delete the simple point
                segment[iv] = 0;

                // add the new surface voxels
                if (segment[IndicesToIndex(ix - 1, iy, iz)] == 1) {
                    NewSurfaceVoxel(IndicesToIndex(ix - 1, iy, iz), ix - 1, iy, iz);
                    segment[IndicesToIndex(ix - 1, iy, iz)] = 2;
                }
                if (segment[IndicesToIndex(ix + 1, iy, iz)] == 1) {
                    NewSurfaceVoxel(IndicesToIndex(ix + 1, iy, iz), ix + 1, iy, iz);
                    segment[IndicesToIndex(ix + 1, iy, iz)] = 2;
                }
                if (segment[IndicesToIndex(ix, iy - 1, iz)] == 1) {
                    NewSurfaceVoxel(IndicesToIndex(ix, iy - 1, iz), ix, iy - 1, iz);
                    segment[IndicesToIndex(ix, iy - 1, iz)] = 2;
                }
                if (segment[IndicesToIndex(ix, iy + 1, iz)] == 1) {
                    NewSurfaceVoxel(IndicesToIndex(ix, iy + 1, iz), ix, iy + 1, iz);
                    segment[IndicesToIndex(ix, iy + 1, iz)] = 2;
                }
                if (segment[IndicesToIndex(ix, iy, iz - 1)] == 1) {
                    NewSurfaceVoxel(IndicesToIndex(ix, iy, iz - 1), ix, iy, iz - 1);
                    segment[IndicesToIndex(ix, iy, iz - 1)] = 2;
                }
                if (segment[IndicesToIndex(ix, iy, iz + 1)] == 1) {
                    NewSurfaceVoxel(IndicesToIndex(ix, iy, iz + 1), ix, iy, iz + 1);
                    segment[IndicesToIndex(ix, iy, iz + 1)] = 2;
                }

                // remove this from the surface voxels
                RemoveSurfaceVoxel(ptr);
                changed += 1;
            }
        }
        DestroyPointList(&deletable_points);
    }


    // return the number of changes
    return changed;
}



static void SequentialThinning(void)
{
    // create a vector of surface voxels
    CollectSurfaceVoxels();
    int iteration = 0;
    long changed = 0;
    do {
        changed = ThinningIterationStep();
        iteration++;
        printf("  Iteration %d deleted %ld points\n", iteration, changed);
    } while (changed);
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
                if (segment[linear_index]) nnneighbors++;
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

void CppTopologicalThinning(const char *prefix, const char *lookup_table_directory, long label, float resolution[3])
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

    // create (and clear) the global variables
    segment = std::unordered_map<long, char>(10000000);

    // initialize all of the lookup tables
    InitializeLookupTables(lookup_table_directory);

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
        segment[iv] = 1;
    }
    fclose(fp);

    // get the number of points
    long initial_points = segment.size();
    printf("Label %ld initial points: %ld\n", label, initial_points);

    // call the sequential thinning algorithm
    SequentialThinning();

    // count the number of remaining points
    long num = 0;
    ListElement *LE = (ListElement *) surface_voxels.first;
    std::unordered_set<long> skeleton = std::unordered_set<long>();
    while (LE != NULL) {
        num++;
        skeleton.insert(LE->iv);
        LE = (ListElement *)LE->next;
        
    }

    std::unordered_set<long> upsampled_skeleton = std::unordered_set<long>();
    for  (std::unordered_set<long>::iterator it = skeleton.begin(); it != skeleton.end(); ++it) {
        long is = *it;

        bool endpoint = false;
        if (IsEndpoint((is))) endpoint = true;

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
            
            if (skeleton.find(neighbor_index) == skeleton.end()) continue;

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
    sprintf(output_filename, "baselines/topological-thinnings/%s/%06ld.pts", prefix, label);

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

    delete[] lut_simple;
    delete[] lut_isthmus;

    double total_time = (double) (clock() - start_time) / CLOCKS_PER_SEC;

    char time_filename[4096];
    sprintf(time_filename, "running_times/topological-thinnings/%s-%06ld.time", prefix, label);

    FILE *tfp = fopen(time_filename, "wb");
    if (!tfp) { fprintf(stderr, "Failed to write to %s.\n", time_filename); exit(-1); }

    // write the number of points and the total time to file
    if (fwrite(&initial_points, sizeof(long), 1, tfp) != 1) { fprintf(stderr, "Failed to write to %s.\n", time_filename); exit(-1); }
    if (fwrite(&total_time, sizeof(double), 1, tfp) != 1) { fprintf(stderr, "Failed to write to %s.\n", time_filename); exit(-1); }

    // close file
    fclose(tfp);
}
