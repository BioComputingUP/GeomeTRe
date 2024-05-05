- TR_geometry.py is the base script, containing the algorithm
- run_TR_geometry.sh is a bash script that uses GNU parallel to run TR_geometry on a batch of proteins. Command line usage is as follows:
bash run_TR_geometry.sh input_file.tsv njobs output_file.tsv
Where :
- input_file.tsv is a tab separated file with 3 columns: the pdb id+chain (ex:ae4gA), the residue numbers delimiting the units (as in TR_geometry_documentation.md, and the residues delimiting the insertions (NA if no insertions)
- njobs is the number of jobs to run in parallel
- output_file.tsv is the path to the output file, which will be in a tsv format, and contain the following columns: pdb id, chain, start of the repeated region, end of the repeated region, and mean,std for all the parameters: curvature, twist (value and sign), pitch (value and sign), TM-align similarity, residual yaw.
  
The .cif files corresponding to the pdb ids will be downloaded to a temporary directory, that is deleted at the end
