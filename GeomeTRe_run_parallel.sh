#!/bin/bash
export files_dir=$(tr -dc A-Za-z0-9 </dev/urandom | head -c 6; echo)
mkdir "$files_dir"
do_GeomeTRe () {
    pdb="${1::-1}"
    chain="${1: -1}"
    units="$2"
    ins="$3"
    curl -s -o "$files_dir"/"$pdb".cif https://files.rcsb.org/download/"$pdb".cif 
    if [ "$ins" = 'NA' ]; then
        python3 -W ignore ./GeomeTRe.py "$files_dir"/"$pdb".cif "$chain" "$units" --batch
    else
        python3 -W ignore ./GeomeTRe.py "$files_dir"/"$pdb".cif "$chain" "$units" -ins "$ins" --batch
    fi
     rm -f "$files_dir"/"$pdb".cif  # Delete the CIF file immediately after processing
}

export -f do_GeomeTRe
echo -e "pdb\tchain\tstart\tend\tcurv_mean\tcurv_std\ttwist_mean\ttwist_std\t \
twist_sign_mean\ttwist_sign_std\tpitch_mean\t \
pitch_std\tpitch_sign_mean\tpitch_sign_std\t \
tmscores_mean\ttmscores_std\tyaw_mean\tyaw_std" > "$3"
parallel -a "$1" -j "$2"  --bar --colsep '\t' \
do_GeomeTRe {} >> "$3"
rm -rf "$files_dir"
