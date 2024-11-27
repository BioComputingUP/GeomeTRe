### TR_geometry.py
The base script, containing the algorithm


### TR_geometry_documentation.md
# TR_geometry.py documentation


## Command line syntax
python TR_geometry.py input_path chain units_def -o output_path -ins insertions --draw --batch
- Input_path: path to .pdb file to analyze
- Chain: chain to analyze
- Units_def: insertion codes delimiting the units, written as s1-e1,s2-e2,...
- -ins insertions: insertion codes delimiting insertions, written as the units
- -o output_path: optional argument. If used, the output of the program will be saved as a csv file within the specified directory. Any directory in the path must already exist.
- --draw: if used, the output will include a PyMOL drawing of the protein structure
- --batch: used by run_TR_geometry.sh (print statistics on the parameters as a tsv line)

Using the -h option displays this info.

EX: python3 TR_geometry.py Output/4cj9.pdb A 30_60,61_93,94_126,127_159,160_192,193_225,226_258,259_291,292_324,325_357,358_390,391_423,424_456,457_489,490_522,523_555,556_588,589_621,622_654,655_687,688_720,721_753 -o Outputs --draw
## Required Dependencies
This package requires the following dependencies to run:
- Python 3.6 or higher
- Required Python packages (installed automatically via `pip install`):
  - numpy
  - pandas
  - scikit-learn
  - scikit-image
  - biopython
  - tmtools
- **PyMOL**: PyMOL must be installed via `conda` before running the package if you intend to use the `--draw` option.
# Install these using conda
conda install -c conda-forge pymol-open-source

## What it does

### Rotation

The program first computes the rotation between units in the following way:
First we take sliding windows of 6 units, projecting them on a plane with PCA, and finding the rotation center as the center of the widest circular crown containining them all. We then reverse the PCA transformation to get the center of rotation in the 3D space.
Then for each pair of units, we find a center of rotation for that pair as the center of rotation found by a window including that pair that corresponded to the largest crown.

### Twist, pitch and handedness

We proceed selecting a sliding pair of consecutive units

Then for each unit of the pair, the program computes two reference axes (twist and pitch axis).

The pitch axis is the vector connecting the geometric center to the center of rotation relative to that pair.

The twist axis is the vector connecting the two barycenters, orthogonalized w.r.t. the two different pitch axes.

Then for each pair of units, we rotate them to bring their axes into correspondence to the standard axes (twist to (1,0,0), pitch to (0,1,0))

We then use TM-align to find the best rotation that overlaps the first unit on the second one.

Then we decompose the rotation of the alignment w.r.t. the reference axes of the first unit to find the pitch, twist and handedness, using Euler angles. 

## Output

There are two outputs: a table with the computed parameters, and a pymol drawing of the geometry of the molecule.
In the table, the first row is all zeros, since the rows refer to the unit and the unit before it.
There are also rows with mean and standard deviation of the parameters (ignoring the first row.)
The table has the following columns:
- pdb_id: the pdb id of the molecule, plus the region
- Chain: the chain
- Unit start: the start of the unit
- Unit end: the end of the unit
- Curvature: the curvature, computed as the angle between the vectors connecting the rotation center to two consecutive units
- Twist: the twist, computed as the component of the rotation that aligns the two units w.r.t. the twist axis
- Twisthandedness: computed as the handedness of the rotation, w.r.t. the twist axis.
- Pitch: the pitch, computed the same way as the twist, but orthogonalizing w.r.t the pitch axis.
- Pitchhandedness: computed as the handedness of the rotation, w.r.t. the pitch axis.
- TM-score: the tm-score of the structural alignment
- Yaw: the residual yaw rotation: in a perfect structure, this is 0, as we already compensate for the yaw when we align the axes of the two units to the standard reference axes. A high yaw means a bad performance on the algorithm for that unit pair.

The pymol drawing contains the following objects:
- In yellow, the line connecting the geometric centers of each unit to each other, and the lines connecting the rotation centers to each geometric center
- In green, the twist axis of each unit
- In blue, the pitch axis of each unit. Units that are at the edge of the 6 units window have two pitch axes (one for each rotation center)
- In red, orange and white, the components found by PCA.

### run_TR_geometry.sh
A bash script that uses GNU parallel to run GeomeTRe on a batch of proteins. Command line usage is as follows:
bash run_GeomeTRe.sh input_file.tsv njobs output_file.tsv.

Where :
- input_file.tsv is a tab separated file with 3 columns: the pdb id+chain (ex:ae4gA), the residue numbers delimiting the units (as in GeomeTRe_documentation.md, and the residues delimiting the insertions (NA if no insertions)
- njobs is the number of jobs to run in parallel
- output_file.tsv is the path to the output file, which will be in a tsv format, and contain the following columns: pdb id, chain, start of the repeated region, end of the repeated region, and mean,std for all the parameters: curvature, twist (value and sign), pitch (value and sign), TM-align similarity, residual yaw.
  
The .cif files corresponding to the pdb ids will be downloaded to a temporary directory, that is deleted at the end

### repeatsdb.tsv
tsv file containining the proteins in the RepeatsDB database (version 4), in the format needed for run_GeomeTRe.sh
