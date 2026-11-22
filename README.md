# 16S amplicon workflow in Nextflow (QIIME 2 + Docker)

Nextflow port of the first stages of a 16S rRNA amplicon (V4) OTU workflow.
Every process runs inside a pinned QIIME 2 Docker image.

## What it does

1. **IMPORT**: builds a manifest from the paired FASTQ files and imports them into QIIME 2.
2. **MERGE**: merges read pairs into contigs (vsearch) and discards contigs containing ambiguous bases (`--maxns`, default 0).
3. **DEREP**: per-sample dereplication, keeping abundance and sample identity.
4. **CHIMERA**: UCHIME de novo chimera detection (vsearch), then filters the feature table and sequences.

Final output: a chimera-free feature table and its sequences.

## Scope

Only these four stages are ported to Nextflow. The remaining stages of the full workflow
(closed-reference OTU clustering against SILVA, taxonomy filtering, tree, taxonomic profile
and alpha rarefaction) were run as a separate bash workflow and are not part of this repository.

## Requirements

| Software | Version used |
|---|---|
| Nextflow | 26.10.7 (build 15252) |
| Docker Engine | 29.10.0 |
| Java | 17 or later (Nextflow requirement) |
| QIIME 2 image | `quay.io/qiime2/qiime2:2026.10` (pulled automatically on first run; several GB) |

Tested on Ubuntu 24.04 under WSL2.

## Data

The example data is the mothur MiSeq SOP dataset (mouse gut samples and a mock community,
paired-end 2x250 V4 reads). Download it from the MiSeq SOP page: https://mothur.org/wiki/miseq_sop/

Extract it and make the FASTQ files available in a `data/` folder next to `main.nf`, for example:

```
ln -s /path/to/MiSeq_SOP data
```

The pipeline expects files named `<sample>_<...>_R1_001.fastq` and `<sample>_<...>_R2_001.fastq`.
The FASTQ data is not included in this repository.

## Run

```
nextflow run main.nf
```

Re-run only what changed with `-resume`:

```
nextflow run main.nf -resume
```

Parameters:

| Parameter | Default | Meaning |
|---|---|---|
| `--reads` | `data/*_R{1,2}_001.fastq` | Glob for the paired FASTQ files |
| `--maxns` | `0` | Maximum ambiguous bases allowed per merged contig |

The output folder is set to `results/` in `nextflow.config` (override with `-output-dir`).

## Outputs

| File | Contents |
|---|---|
| `results/import/demux.qza` | Imported paired-end reads |
| `results/merge/merged.qza` | Merged contigs |
| `results/derep/derep_table.qza`, `derep_seqs.qza` | Dereplicated feature table and sequences |
| `results/chimera/nochim_table.qza`, `nochim_seqs.qza` | Chimera-free feature table and sequences |

## Reproducibility check

On the MiSeq SOP example data, the chimera-free feature table from this pipeline was compared
with the table produced by the equivalent step-by-step bash workflow. The exported tables
were identical.

## Notes and limitations

- `bin/make_manifest.sh` assumes Illumina-style file names and uses the text before the first
  underscore as the sample ID. Other naming schemes need a different manifest step.
- No primer trimming is included, because the MiSeq SOP reads contain no primers. Data that
  still contains primers needs a trimming step added before merging.
- The mock community sample (`Mock`) is processed together with the other samples.
- After downloading the repository, make the helper script executable with `chmod +x bin/make_manifest.sh` (uploading files through the GitHub web page does not keep the executable permission).
