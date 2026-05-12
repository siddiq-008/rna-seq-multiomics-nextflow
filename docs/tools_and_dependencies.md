# Tools and dependencies

## Required tools

| Tool | Purpose |
|---|---|
| Nextflow | Workflow management |
| Java | Required by Nextflow |
| SRA Toolkit | Download FASTQ files from SRA |
| fastp | Read trimming and quality filtering |
| HISAT2 | RNA-seq read alignment |
| SAMtools | BAM/SAM processing |
| featureCounts | Gene-level read counting |
| eggNOG-mapper | Functional annotation |
| Python 3 | Integration script |
| pandas | Data processing in Python |

## Conda installation example

```bash
conda create -n rnaseq_nextflow -c bioconda -c conda-forge \
    nextflow sra-tools fastp hisat2 samtools subread eggnog-mapper pandas
```

Activate the environment:

```bash
conda activate rnaseq_nextflow
```

## Python-only dependency

```bash
pip install -r requirements.txt
```
