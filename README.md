# RNA-seq Multi-omics Integration Pipeline using Nextflow

This repository contains a Nextflow-based RNA-seq analysis and multi-omics integration pipeline.

The pipeline processes RNA-seq data, performs functional annotation using eggNOG-mapper, and integrates gene expression, annotation, KEGG pathway information, and optional metabolic model-derived information into a final multi-omics CSV file.

## Workflow overview

The pipeline includes the following steps:

1. Download paired-end FASTQ files from SRA
2. Trim and quality-filter reads using `fastp`
3. Align reads to a reference genome using `HISAT2`
4. Generate gene-level counts using `featureCounts`
5. Annotate protein sequences using `eggNOG-mapper`
6. Integrate expression, annotation, KEGG, and optional metabolic model information using Python

## Repository structure

```text
rna-seq-multiomics-nextflow/
│
├── README.md
├── main.nf
├── nextflow.config
├── requirements.txt
├── .gitignore
├── LICENSE
│
├── bin/
│   └── integrate.py
│
├── config/
│   └── params.config
│
├── assets/
│   └── workflow_diagram.png
│
├── data/
│   └── README.md
│
├── reference/
│   └── README.md
│
├── results/
│   └── README.md
│
├── docs/
│   ├── pipeline_overview.md
│   └── tools_and_dependencies.md
│
└── example/
    ├── sample_params.config
    └── example_run_command.txt
```

## Requirements

The following tools are required:

- Nextflow
- Java
- SRA Toolkit
- fastp
- HISAT2
- SAMtools
- Subread / featureCounts
- eggNOG-mapper
- Python 3
- pandas

## Install Python dependencies

```bash
pip install -r requirements.txt
```

## Conda installation example

```bash
conda create -n rnaseq_nextflow -c bioconda -c conda-forge \
    nextflow sra-tools fastp hisat2 samtools subread eggnog-mapper pandas

conda activate rnaseq_nextflow
```

## Input files

| Input | Description |
|---|---|
| SRR ID | SRA accession ID for RNA-seq data |
| `annotation.gtf` | Genome annotation file |
| HISAT2 index | Pre-built HISAT2 reference index |
| `protein.faa` | Protein FASTA file for eggNOG annotation |
| `important_genes.csv` | Expression/count table used for integration |
| Metabolic model CSV | Optional metabolic model table |

## Configuration

Edit `nextflow.config` before running the pipeline:

```nextflow
params {
    srr = 'SRR4101836'

    gtf     = 'reference/annotation.gtf'
    index   = 'reference/hisat2_index/genome'
    protein = 'reference/protein.faa'
    counts  = 'data/important_genes.csv'

    eggnog_db = '/path/to/eggnog_db'
    metabolic_model = ''
    locus_prefix = 'DJ41_'
    outdir = 'results'
}
```

If you are not using metabolic model integration, keep:

```nextflow
metabolic_model = ''
```

## Running the pipeline

```bash
nextflow run main.nf
```

Run with custom parameters:

```bash
nextflow run main.nf \
    --srr SRR4101836 \
    --gtf reference/annotation.gtf \
    --index reference/hisat2_index/genome \
    --protein reference/protein.faa \
    --counts data/important_genes.csv \
    --eggnog_db /path/to/eggnog_db
```

Run with optional metabolic model integration:

```bash
nextflow run main.nf \
    --metabolic_model data/msystems.00157-18-sd001_converted.csv \
    --locus_prefix DJ41_
```

## Output files

The main outputs are written inside the `results/` directory.

| Folder | Output |
|---|---|
| `results/fastq` | Raw FASTQ files downloaded from SRA |
| `results/fastp` | Trimmed FASTQ files |
| `results/bam` | Aligned BAM file |
| `results/counts` | Gene-level count table |
| `results/eggnog` | eggNOG annotation output |
| `results/integration` | Final integrated multi-omics CSV |

Final output:

```text
results/integration/final_complete_multiomics.csv
```

## Final integrated output columns

Depending on input availability, the final CSV may contain:

- RNA locus tag
- Metabolic locus tag
- Expression value
- Gene symbol
- Functional description
- COG annotation
- GO annotation
- KEGG ID
- KEGG pathway ID
- Metabolite name
- Metabolite pathway ID
- Metabolite pathway name

## Important GitHub note

Do not upload large files such as:

- FASTQ files
- BAM/SAM files
- SRA files
- HISAT2 indexes
- eggNOG databases
- Large reference genomes
- Full result folders

Those files are ignored using `.gitignore`.

## Author

Mohamed Siddiq
