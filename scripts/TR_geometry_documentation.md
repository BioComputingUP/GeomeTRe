# TR_geometry.py documentation


## Command line syntax
python TR_geometry.py input_path chain units_def -o output_path -ins insertions--draw
- Input_path: path to .pdb file to analyze
- Chain: chain to analyze
- Units_def: insertion codes delimiting the units, written as s1-e1,s2-e2,...
- -ins insertions: insertion codes delimiting insertions, written as the units
- -o output_path: optional argument. If used, the output of the program will be saved as a csv file within the specified directory. Any directory in the path must already exist.
- --draw: if used, the output will include a PyMOL drawing of the protein structure

Using the -h option displays this info.

EX: python TR_geometry.py Data\4cj9.pdb A 30-60,61-93,94-126,127-159,160-192,193-225, 226-258,259-291,292-324,325-357,358-390,391-423,424-456,457-489,490-522,523-555,
556-588,589-621,622-654,655-687,688-720,721-753 -o Outputs --draw
## Required packages

- Argparse and pathlib for command line processing 
- Pymol for the final visualization of the structure
- Numpy and Pandas
- Biopython for the processing of the .cif file
- Scipy for its Rotation class
- sklearn for the Principal Component Analysis
- skimage for the circle fitting function

## What it does

### Rotation

The program first computes the rotation between units in the following way:
First we take sliding windows of 6 units, projecting them on a plane with PCA, and fitting a circle on the resulting points. We then reverse the PCA transformation to get the center of rotation in the 3D space.
Then for each pair of units, we find a center of rotation for that pair as the center of rotation that included that pair and gave the best fit in terms of RMSD.

### Twist, pitch and handedness

We proceed selecting a sliding pair of consecutive units

Then for each unit of the pair, the program computes two reference axes (twist and pitch axis).

The pitch axis is the vector connecting the geometric center to the center of rotation relative to that pair.

The twist axis is the vector connecting the two barycenters, orthogonalized w.r.t. the two different pitch axes.



Then for each pair of units, we rotate the reference axes of the second one to overlap those of the first one.
We then apply this rotation to the second unit, to bring it within the reference frame of the first one.

We then represent each unit as its components found by PCA.
Since sometimes the components may change in order between units, we match components by similarity, then make sure the verses are equal.
We then find the rotation that aligns the units as the rotation that aligns the matched components of their PCAs.

Then we decompose the rotation of the alignment w.r.t. the reference axes of the first unit to find the pitch, twist and handedness.
We do this by taking a vector, applying the rotation to it, orthogonalizing these two vectors w.r.t. the reference axis and looking at the resulting angle.

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
- Handedness: computed as the handedness of the rotation, w.r.t. the twist axis.
- Pitch: the pitch, computed the same way as the twist, but orthogonalizing w.r.t the pitch axis.

The pymol drawing contains the following objects:
- In yellow, the line connecting the geometric centers of each unit to each other, and the lines connecting the rotation centers to each geometric center
- In green, the twist axis of each unit
- In blue, the pitch axis of each unit. Units that are at the edge of the 6 units window have two pitch axes (one for each rotation center)
- In red, orange and white, the components found by PCA.
