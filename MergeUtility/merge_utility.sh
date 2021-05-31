#!/bin/bash
# Defaults:
input_dir="./"
ext='001.fastq.gz'

# Usage:
echo " Merge Illumina Fastqs based on Lanes
 USAGE:
 $(basename $0) -o OUTPUT_DIR [-i INPUT_DIR]
 Options:
   -i    Input directory (default: $input_dir)
   -o    Output directory (will be created)
   -e    Extension without first dot (default: $ext)
"
while getopts o:i:e: option; do
        case "${option}" in

        i) input_dir=${OPTARG} ;;
        o) output_dir=${OPTARG} ;;
        e) ext=${OPTARG} ;;
        ?) echo " Unnknown parameter $OPTARG" ;;
        esac
done
shift "$(($OPTIND - 1))"

if [ -z ${output_dir+x} ]; then
        echo " FATAL ERROR: Please specify output directory:  -o OUTPUT_DIR"
        exit 9
fi

if [ -d "${output_dir}" ]; then
        echo " FATAL ERROR: Directory '$output_dir' was found. Please specify a new name"
        exit 7
fi

if test "$BASH" = "" || "$BASH" -uc "a=();true \"\${a[@]}\"" 2>/dev/null; then
        set -euo pipefail
else
        set -eo pipefail
fi
shopt -s nullglob globstar
IFS=$'\n\t'
mkdir "$output_dir"
# Loop files in {input_dir} with extension {ext}
for sample_file in ${input_dir}*_*.${ext}; do
        sample_name=$(basename "$sample_file" | cut -f 1 -d "_")
        sample_index=$(basename "$sample_file" | cut -f 2 -d "_")
        sample_strand=$(basename "$sample_file" | cut -f 4 -d "_")

        echo " > Adding $sample_file to ${sample_name}_${sample_strand}.${ext}"
        cat $sample_file >>${output_dir}/${sample_name}_${sample_index}_L001_${sample_strand}_${ext}
done
cd $output_dir/
echo "Renaming ..."

for sample_file in ./*.${ext}; do
        sample_name=$(basename "$sample_file" | cut -f 1 -d "_")
        sample_index=$(basename "$sample_file" | cut -f 2 -d "_")
        sample_strand=$(basename "$sample_file" | cut -f 4 -d "_")
        find . -type f -name "${sample_name}_${sample_index}_L001_${sample_strand}_001.${ext}" \
                -exec mv {} ./${sample_name}_${sample_index}_${sample_strand}_${ext} \;
done
echo "Done!"
