
## Synapse-Aware Skeleton Generation

The skeleton benchmark and code for the 2019 MICCAIR paper *Synapse-Aware Skeleton Generation for Neural Circuits* [1]. For more information: https://www.rhoana.org/synapseaware.

### Installation

````
git clone https://github.com/Rhoana/synapseaware.git .
cd synapseaware
conda create -n synapseaware_env python=3.7
conda install --file requirements.txt
cd connectome
python setup.py build_ext --inplace
cd ../isthmus 
python setup.py build_ext --inplace
cd ../teaser
python setup.py build_ext --inplace
````

Add the parent directory to this repository to your PYTHONPATH variable. 

### Meta Files

Each new dataset needs a meta file named meta/{PREFIX}.meta where {PREFIX} is a unique identifier for the dataset. All functions in this repository require as input this {PREFIX} identifier to find information such as resolution and grid size. An example meta file is provided in `examples/meta/Fib25.meta`. 

### Directory Structure

For each new {PREFIX}, you need to create the following directories:

````
baselines/teasers/{PREFIX}
baselines/topological-thinnings/{PREFIX}
connectomes/{PREFIX}
distances/{PREFIX}
segmentations/{PREFIX}
skeletons/{PREFIX}
somae/{PREFIX}
surfaces/{PREFIX} (optional)
synapses/{PREFIX}
volumetric_somae/segmentations/{PREFIX} (optional)
volumetric_somae/surfaces/{PREFIX}
widths/{PREFIX}
````

If the segment does not contain the somae, a central point should be chosen as the soma and be copied to the `somae` and `volumetric_somae/surfaces` directories. The algorithm also assumes the following directories:

````
meta
running_times/refinement
running_times/skeletons
width-errors (optional)
````

The script ``examples/scripts/create-folders.py`` creates the intended directory structure. 

### Example Script

The script ``examples/scripts/connectome.py`` extracts the skeleton for an example from the FIB-25 dataset publicly available online [2]. A central point was chosen as the soma location since this particular segment does not contain the cell body.

### Input Format

All of the point clouds have the format:

````
unsigned long: Z Grid Size
unsigned long: Y Grid Size
unsigned long: X Grid Size
unsigned long: number of points 
unsigned long *: points in linear index
````

The file ``utilities/dataIO.py`` contains example input/output operations for this file format. Use the following methods to convert between Cartesian coordinates and the linear index:

````
# iv is the linear index
# yres is the Y Grid Size
# xres is the X Grid Size
def IndexToIndices(iv):
	iz = iv // (yres * xres)
	iy = (iv - iz * yres * xres) // xres
	ix = iv % xres
	return iz, iy, ix

def IndicesToIndex(ix, iy, iz):
	return iz * yres * xres + iy * xres + ix
````

### Citations

[1] Matejek, B., Wei, D., Wang, X., Zhao, J., Palágyi, K., and Pfister, H., Synapse-Aware Skeleton Generation for Neural Circuits. In *International Conference on Medical Image Computing and Computer-Assisted Intervention, pp. YYY-YYY. Springer, Cham, 2019. 

[2] Takemura, S.Y., Xu, C.S., Lu, Z., Rivlin, P.K., Parag, T., Olbris, D.J., Plaza, S., Zhao, T., Katz, W.T., Umayam, L. and Weaver, C., 2015. Synaptic circuits and their variations within different columns in the visual system of Drosophila. _Proceedings of the National Academy of Sciences_, _112_(44), pp.13711-13716.
