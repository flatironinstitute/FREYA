
"cmwf_csv.sh" is the processing pipeline that takes as input HTS
*.fastq files as well as various reference files, and produces as
output various summary files and of *.vcf files. Please see block
comments with the script for more information.

This pipeline has several dependencies, which are, in effect,
described by the Dockerfile.

"VCF_mutation_picker.0.6.py" is invoked as the last step of the
pipeline. It produces a per-gene summary of results.


