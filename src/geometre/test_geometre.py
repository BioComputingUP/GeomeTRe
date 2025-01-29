from geometre import GeomeTRe

# Initialize GeomeTRe with PyMOL drawing enabled
geometre = GeomeTRe(draw_enabled=True)

# Process a single PDB file
output_path = geometre.compute(
    "single",
    filepath="/home/zarifa/Downloads/2xqh.pdb",
    chain="A",
    units_ids="161_175,176_189,190_203,204_217,218_233,234_249,250_263,264_276,305_326,327_350,373_392,393_416",
    ins_ids="351_372"
)

# Display results from the output CSV file
geometre.display_results(output_path)
