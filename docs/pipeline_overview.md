# Pipeline overview

This pipeline is designed for RNA-seq based multi-omics integration.

## Step 1: FASTQ download

The pipeline downloads paired-end FASTQ files from SRA using the provided SRR accession ID.

## Step 2: Read trimming

Raw reads are cleaned using `fastp`.

## Step 3: Alignment

Cleaned reads are aligned to the reference genome using `HISAT2`.

## Step 4: Read counting

Gene-level read counts are generated using `featureCounts` and a GTF annotation file.

## Step 5: Functional annotation

Protein sequences are annotated using `eggNOG-mapper`.

## Step 6: Multi-omics integration

The Python integration script combines expression, gene annotation, KEGG pathway information, and optional metabolic model-derived information into a final CSV file.

## Final output

```text
final_complete_multiomics.csv
```
