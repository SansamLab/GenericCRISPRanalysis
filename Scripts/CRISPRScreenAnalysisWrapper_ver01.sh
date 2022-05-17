#!/bin/bash -l

# Chris Sansam and Kevin Boyd
# 2020-06-14; version 01
# Script for running MAGeCK pipeline

################################################################################
# function to provide usage help
################################################################################

usage()
{
  cat <<EOF


  Usage: $0 
            [ -o Object_store_paths ]
            [ -s SRA_IDs ]
            [ -b sample_labels ]
            [ -q fastq filenames ] 
            [ -f fastq_directory ] 
            [ -l sgRNA_library_file ]
            [ -d design_matrix_file ]
            [ -n negative_control_sgRNA_file ]

  note:  When multiple arguments are provided for a single option, they must
         all be enclosed in a set of quotes. For example:
         -s "SRR4880321 SRR4880322"


EOF
  exit 2
}

################################################################################
# print all arguments to log file
################################################################################

cat <<EOF
$0
Arguments:
$@


EOF

################################################################################
# load necessary modules
################################################################################

ml ncbi_sra/2.8.2-1 o3-utils mageck fastqc multiqc aspera

################################################################################
# set arguments to variables
################################################################################

unset SRA_IDs sgRNA_library_file design_matrix_file sample_labels \
negative_control_sgRNA_file Object_store_paths fastq_directory fastq_filenames \

while getopts 'o:s:b:q:f:l:d:n:' OPTION; do
  case "$OPTION" in
    o) Object_store_paths="$OPTARG" ;;
    s) SRA_IDs="$OPTARG" ;;
    b) sample_labels="$OPTARG" ;;
    q) fastq_filenames="$OPTARG" ;;
    f) fastq_directory="$OPTARG" ;;
    l) sgRNA_library_file="$OPTARG" ;;
    d) design_matrix_file="$OPTARG" ;;
    n) negative_control_sgRNA_file="$OPTARG" ;;
    \?) usage ;;
  esac
done
shift "$(($OPTIND -1))"

################################################################################
# if it is specified but doesn't exist, make the fastq directory
################################################################################

if [ -n "$fastq_directory" ]; then
  mkdir -p $fastq_directory
else
  echo "The fastq_directory was not specified, so fastq files were stored in scratch."
fi

################################################################################
# if fastq directory is not specified; set it to a scratch directory
################################################################################

if [ -z "$fastq_directory" ]; then
  fastq_directory=$(mktemp -d)
  cleanup() {
    err=$?
    rm -rf "${fastq_directory}"
    trap '' EXIT INT TERM
    exit $err
  }
  sig_cleanup() {
    err=$?
    trap '' EXIT
    cleanup $err
  }
  trap cleanup EXIT
  trap sig_cleanup INT QUIT TERM
fi

################################################################################
# transfer and unzip any "-o"-specified object store fastq files to the 
# fastq_directory
################################################################################

if [ -n "$Object_store_paths" ]; then
  for path in $Object_store_paths; do
      bsenme=$(basename -s ".gz" $path)
      bsenme="$fastq_directory""/""$bsenme"        
      o3-do.rb $path "gunzip -c {1} > $bsenme"
  done
else
  echo "The -o is missing, so fastq files were not transferred from object store."
fi

################################################################################
# download any "-s"-specified SRA fastq files to the fastq_directory
################################################################################

if [ -n "$SRA_IDs" ]; then
  for SRA_ID in $SRA_IDs; do
    prefetch "$SRA_ID"
    fastq-dump --outdir "$fastq_directory" "$SRA_ID"
  done
else
  echo "The -s is missing, so fastq files were not transferred from SRA."
fi

rm /s/sansam-lab/*.sra

################################################################################
# run fastqc and multiqc on all fastq files
################################################################################

if [[ -n "$SRA_IDs"  || -n "$Object_store_paths" ]]; then
  mkdir -p fastqc_results
  mkdir -p fastqc_results/multiqc_results
  echo "Running fastqc."
  fastqc -o fastqc_results ${fastq_directory}/*.fastq
  echo "Running multiqc."
  multiqc -o fastqc_results/multiqc_results fastqc_results/.
else
  echo "No new fastq files were transferred, so fastqc was not run."
fi

################################################################################
# make count table from all fastq files with MAGeCK count
################################################################################

if [ -n "$sgRNA_library_file" ]; then
  # prepend all filenames with scratch directory
    allFastQFiles=$(printf "${fastq_directory}/%s\n" $fastq_filenames)
    allFastQFiles=$(echo $allFastQFiles|sed "s|,|,$fastq_directory/|g")
  echo "Running mageck count."
  mageck count \
  -l $sgRNA_library_file \
  --fastq $allFastQFiles \
  --sample-label $sample_labels
else
  echo "The -l argument is missing, so MAGeCK count was not run."
fi

################################################################################
# compare samples with MAGeCK mle
################################################################################

if [ -n "$design_matrix_file" ]; then
  #day0Label=$(echo $sample_labels | sed 's/,.*$//')
  echo "Running mageck mle."
  mageck mle \
  -k sample1.count.txt \
  -d $design_matrix_file \
  --threads 12 \
  --control-sgrna $negative_control_sgRNA_file
else
  echo "The -d argument is missing, so MAGeCK mle was not run."
fi


