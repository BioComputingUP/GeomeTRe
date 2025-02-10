from geometre import GeomeTRe

# Initialize GeomeTRe with drawing disabled for batch mode
geometre = GeomeTRe(draw_enabled=False)

# Run batch processing
output_path = geometre.compute(
    "batch",
    tsv_path="/home/zarifa/Desktop/papers/Geometry/GeomeTRe/input_data/test/test.tsv",
    num_threads=4,
    file_format="pdb",
    pdb_dir="/home/zarifa/Desktop/papers/Geometry/GeomeTRe/input_data/test/pdb_dir",
    output_path="/home/zarifa/Downloads/cirkiiin.csv"
)

print(f"Batch processing completed! Results saved at: {output_path}")
