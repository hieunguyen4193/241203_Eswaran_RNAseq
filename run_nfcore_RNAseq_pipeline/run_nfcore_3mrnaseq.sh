#!/bin/bash

nfcore_version=3.17.0

# Check if the --test option is provided
if [ "$#" -gt 0 ] && [ "$1" = "--test" ]; then
  ################## test script #########################
  nextflow run nf-core/rnaseq -r $nfcore_version \
  -profile test,docker --outdir results
else
  ################## GPM samplesheet #####################
  # gpm samplesheet-rnaseq -st 'forward' -sn true -si 2 samplesheet.csv /data/fastq/241203_NB501289_0892_AHJVVTBGXV/241203_Eswaran_KJM_3mRNAseq/

  ################## Run nfcore pipeline #################
  nextflow run nf-core/rnaseq -r ${nfcore_version} -profile docker -c nextflow.config -resume \
  --input samplesheet.csv --outdir results \
  --genome GENCODE_GRCm39_v35_ERCC --skip_biotype_qc --gencode --featurecounts_group_type gene_type \
  --extra_salmon_quant_args="--noLengthCorrection" \
  --extra_star_align_args="--alignIntronMax 1000000 --alignIntronMin 20 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --outFilterMismatchNmax 999 --outFilterMultimapNmax 20 --outFilterType BySJout --outFilterMismatchNoverLmax 0.1 --clip3pAdapterSeq AAAAAAAA" \
  --with_umi --umitools_extract_method="regex" --umitools_bc_pattern="^(?P<umi_1>.{6})(?P<discard_1>.{4}).*"
fi

###### For rerun the pipeline #################################
# use: -resume

###### Options for genome: ####################################
# GENCODE_GRCh38_v46, GENCODE_GRCh38_v46_ERCC, 
# GENCODE_GRCm39_v35, GENCODE_GRCm39_v35_ERCC
# --gencode --featurecounts_group_type gene_type

###### QuantSeq 3’mRNA-Seq V2 Library Prep Kit with UDI ########
# --extra_star_align_args="--alignIntronMax 1000000 \
# --alignIntronMin 20 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 \
# --outFilterMismatchNmax 999 --outFilterMultimapNmax 20 \
# --outFilterType BySJout --outFilterMismatchNoverLmax 0.1 \
# --clip3pAdapterSeq AAAAAAAA" \
# --with_umi --umitools_extract_method="regex" \
# --umitools_bc_pattern="^(?P<umi_1>.{6})(?P<discard_1>.{4}).*"

###### For any 3'mRNAseq #######################################
# --extra_salmon_quant_args="--noLengthCorrection"

###### For runs with ERCC spike-ins ############################
# --genome GENCODE_GRCh38_v46_ERCC --skip_biotype_qc

