# GeomeTRe 

GeomeTRe calculates geometrical properties of tandem repeat proteins. 
It requires a protein structure and the start and end positions 
of each repeat unit (with optional insertion positions) as input. 
If insertion positions are provided, those segments are excluded to improve accuracy. 

For most known STRPs, repeat unit and insertion coordinates are available 
from the manually curated [RepeatsDB](https://repeatsdb.org/) database.

It takes about 2 minutes to calculate 100 structures with 10 threads and a SSD disk.

## Algorithm
The algorithm computes the three Tait-Bryan angles - yaw, pitch, and roll - 
by simulating an airplane traversing the protein from its N-terminus to C-terminus. 
In this analogy, the airplane points to the centroid of the next repeat unit,
and the angles correspond to the maneuvers required to move from one unit to the next. 

The algorithm also determines handedness, defined by the roll direction of movement 
(clockwise/right-handed or anticlockwise/left-handed), and the sign of the pitch 
(positive for upward, negative for downward movement).


## Installation
The software can be installed with `pip` or you can just clone it and use it. 

In order to enable Pymol visualization (optional), we recommend to create a new environment and install the Pymol 
bundle first and after all the other dependencies. 

Installation with pip
```bash
pip install git+https://github.com/BioComputingUP/GeomeTRe.git
```

Installation with Conda and Pymol
```bash
# Create and activate a new conda environment
conda create -n geometre
conda activate geometre

# Install Pymol
conda install -c conda-forge -c schrodinger pymol-bundle

# Clone the GeomeTRe repository
git clone https://github.com/BioComputingUP/GeomeTRe.git

# Set path to the module in your environment
export PYTHONPATH="${PYTHONPATH}:/home/user/Desktop/GeomeTRe/src/"
```

## Dependencies
The following dependencies are required to run:
- Python 3.9 or higher
- Packages (installed automatically via `pip install`):
  - numpy
  - pandas
  - scikit-learn
  - scikit-image
  - biopython
  - tmtools
  
**PyMOL**: PyMOL must be installed via `conda` before running the package if you intend to enable visualization.
`conda install -c conda-forge pymol-open-source`

## Usage
GeomeTRe can be used in single mode to process a single structure, in batch mode to process an
entire dataset and it can be executed for just rendering its results in Pymol. 

Single structure execution 
```bash
# Single mode without pip. Same as above but invoking main.py directly
python3 __main__.py single 2xqh.pdb A result.csv 161_175,176_189,190_203,204_217,218_233,234_249,250_263,264_276,305_326,327_350,373_392,393_416 -ins_def 351_372

# Single mode with pip installation - pdb id, chain, output file, units, insertions (optional)
geometre single 2xqh.pdb A result.csv 161_175,176_189,190_203,204_217,218_233,234_249,250_263,264_276,305_326,327_350,373_392,393_416 -ins_def 351_372
```

Visualize output in Pymol 
```bash
# Visualize without pip
python3 __main__.py draw 2xqh.pdb result.npy

# Visualize with pip installation
geometre pymol 2xqh.pdb result.npy
```

Batch execution 
```bash
# Download structures. It extract PDB IDs from the first column of the TSV file and download them in the pdb_dir folder.
python3 __main__.py batch data/input_batch_short.tsv data/result_batch.csv -pdb_dir data/pdbs -threads 4

# If you don't provide the -pdb_dir argument the program don't download structures based on PDB id, 
# but it expects structures are available in the path provided in the first column of the TSV file (see format section below)
python3 __main__.py batch data/input_batch_short.tsv data/result_batch.csv -threads 4

# With pip installation
geometre batch data/input_batch_short.tsv data/result_batch.csv -threads 4
```


### Library

GeomeTRe can be used as module directly in a Python script:

```python
from geometre.process import compute

df = compute(filepath=input_file, 
			 chain=chain, 
			 units_ids='161_175,176_189,190_203,204_217', 
			 o_path=out_file, 
			 ins_ids='351_372')
```


## Formats 

### Output single mode
CSV table with the computed parameters (.csv)
 	- pdb_id: the PDB id of the molecule
 	- chain: chain of PDB structure
 	- unit_start: start position of repeat unit
 	- unit_end: end position of repeat unit
 	- curvature: the curvature, computed as the angle between the vectors connecting the rotation center to two consecutive units
	- twist: the twist, computed as the component of the rotation that aligns the two units w.r.t. the twist axis
	- twist_hand: computed as the handedness of the rotation, w.r.t. the twist axis
	- pitch: the pitch, computed the same way as the twist, but orthogonalizing w.r.t the pitch axis
	- pitch_hand: computed as the handedness of the rotation, w.r.t. the pitch axis.
	- tm-score: the tm-score of the structural alignment
	- yaw: the residual yaw rotation: in a perfect structure, this is 0, as we already compensate for the yaw when we align the axes of the two units to the standard reference axes. A high yaw 
	  means a bad performance on the algorithm for that unit pair.
 	- Additionally, the last 2 rows are showing mean and standard deviations of each parameter. The first column is all zeros, since the rows refer to the unit and the unit before it.

Pymol parameters for drawing (.npy)

### Pymol drawing
PyMOL drawing contains the following axis:
	- In red, twist axis of each repeat unit(RU) which is always parallel to the longest dimension of the protein
	- In green, pitch axis of each RU
	- In blue, yaw(curvature) axis of each RU
The example of PyMOL drawing in png format with explanation text is below

![Example of PyMOL drawing](/example_2xqh.png)

### Output batch mode

- pdb id: the PDB id of the molecule
- chain: chain of PDB structure
- curvature_mean: mean of curvature
- curvature_std: standard deviation of curvature
- twist_mean: mean of twist
- twist_std: standard deviation of twist
- twist_hand_mean: mean of twist handedness
- twist_hand_std: standard deviation of twist handedness
- pitch_mean: mean of pitch
- pitch_std: standard deviation of pitch
- pitch_hand_mean: mean of pitch handedness
- pitch_hand_std: standard deviation of pitch handedness
- tm-score_mean: mean of tmtool score
- tm-score_std: standard deviation of tmtool score
- yaw_mean: mean of yaw
- yaw_std: standard deviation of yaw


### data/input_batch_short.tsv
tsv file containing the proteins in the RepeatsDB database (version 4) for batch mode calculations.

