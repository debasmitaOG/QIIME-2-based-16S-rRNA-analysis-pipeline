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
| QIIME 2 image | `quay.io/qiime2/qiime2:2026.10` |

Tested on Ubuntu 24.04 under WSL2.

## Data

The example data is the mothur MiSeq SOP dataset (mouse gut samples and a mock community,
paired-end 2x250 V4 reads). Download it from the MiSeq SOP page: https://mothur.org/wiki/miseq_sop/

The FASTQ data is not included in this repository. The pipeline expects files named
`<sample>_<...>_R1_001.fastq` and `<sample>_<...>_R2_001.fastq`.

## How to rerun

**1. Get the repository, the data and the image**

```
git clone <https://github.com/debasmitaOG/QIIME-2-based-16S-rRNA-analysis-pipeline.git>
cd <QIIME-2-based-16S-rRNA-analysis-pipeline>
chmod +x bin/make_manifest.sh
docker pull quay.io/qiime2/qiime2:2026.10
```

Extract the MiSeq SOP data and make the FASTQ files available in a `data/` folder next to `main.nf`,
for example:

```
ln -s /path/to/MiSeq_SOP data
```

**2. Run**

```
nextflow run main.nf
```

**3. Rerun**

Reuse finished work (only processes that changed run again):

```
nextflow run main.nf -resume
```

Start from scratch: from inside the project folder, delete the `work` and `results` folders and the
`.nextflow` cache, then run again.

### Parameters

| Parameter | Default | Meaning |
|---|---|---|
| `--reads` | `data/*_R{1,2}_001.fastq` | Glob for the paired FASTQ files |
| `--maxns` | `0` | Maximum ambiguous bases allowed per merged contig |

The output folder is set to `results/` in `nextflow.config` (override with `-output-dir`).

### Outputs

| File | Contents |
|---|---|
| `results/import/demux.qza` | Imported paired-end reads |
| `results/merge/merged.qza` | Merged contigs |
| `results/derep/derep_table.qza`, `derep_seqs.qza` | Dereplicated feature table and sequences |
| `results/chimera/nochim_table.qza`, `nochim_seqs.qza` | Chimera-free feature table and sequences |

## Findings

All numbers are from the MiSeq SOP example data: 20 samples (19 mouse samples and 1 mock community).

**Nextflow vs bash.** The chimera-free feature table from this Nextflow workflow was exported and
compared with the table from the equivalent step-by-step bash workflow. The two exported tables were
identical.

**Full workflow.** The numbers below come from the bash workflow (stages after chimera removal), not from
this repository.

| Point in the workflow | OTUs | Reads | Samples |
|---|---|---|---|
| After closed-reference OTU clustering (97%) | 393 | 112,383 | 20 |
| After removing chloroplast, mitochondria, archaea and eukaryotes | 387 | 112,209 | 20 |
| After removing the mock community sample (final table) | 359 | 107,928 | 19 |

- The lineage filter removed 6 OTUs and 174 reads.
- Reads per sample in the final table: minimum 2,263, median 4,664, maximum 15,096.
- Rarefaction depth was set to 2,260.
- The phylogenetic tree contains 359 sequences, equal to the number of OTUs in the final table.
- Reference: SILVA 138.2 SSURef NR99, trimmed to alignment positions 13862-23445 (V4), re-clustered to 97% with a
  consensus taxonomy.

Interpretation of the rarefaction curves and the taxonomic barplot: to be added.

## QC decisions and why

| Decision | Reason |
|---|---|
| Discard merged contigs containing any ambiguous base (`--maxns 0`) | Every sequence carried forward has fully called bases |
| Dereplicate within each sample | Collapses identical sequences while keeping their counts, so later steps work on unique sequences |
| UCHIME de novo chimera removal | Chimeras are artificial sequences formed during PCR; UCHIME de novo flags them using the abundance of the dereplicated sequences |
| Closed-reference OTU clustering at 97% against SILVA 138.2; reads without a match are discarded | Each OTU is a SILVA reference sequence with known taxonomy |
| Reference trimmed to alignment positions 13862-23445 | These are the V4 positions given in a QIIME 2 forum tutorial, so the reference covers the same region as the reads |
| Remove chloroplast, mitochondria, archaeal and eukaryotic sequences | They are not part of the bacterial community being profiled |
| Remove the mock community sample before the taxonomic profile and alpha diversity analysis | It is a control of known composition, not a biological sample |
| Rarefaction depth 2,260 | Just below the smallest sample (2,263 reads), so all 19 samples are kept |

**Not done:**
- Primer trimming. The mothur MiSeq SOP does not trim primers from the reads; its primer-region step is applied to the SILVA reference alignment.
- Quality-score filtering, homopolymer screening and length screening.
- Denoising.
- A mock-community error rate.

## What I would not conclude from this data

This workflow groups reads into OTUs at 97% identity against a reference. It does not denoise reads or
infer amplicon sequence variants (ASVs). Because of that:

- **No ASV-level or strain-level conclusions.** Sequences that differ by a few bases fall into the same OTU, so
  any conclusion that depends on exact sequence variants cannot be drawn from this data.
- **No direct comparison with ASV-based results.** OTU counts and observed-OTU values are not the same quantity
  as ASV counts.
- **No conclusions about novel or unclassified organisms, or about total diversity.** Closed-reference
  clustering discards reads that have no match to the reference, so they are excluded by design.
- **No claims about the accuracy of the pipeline.** The mock community was not used to measure an error rate.
- **No statistical conclusions.** Only alpha-diversity rarefaction curves were produced and no statistical tests were run,
  so any difference between groups (for example early and late samples) is descriptive. Beta diversity and
  differential abundance were not analysed.
- **No generalisation beyond this dataset.** It is a small example dataset (19 samples) used to build and
  demonstrate the workflow.

## Notes

- `bin/make_manifest.sh` assumes Illumina-style file names and uses the text before the first underscore as
  the sample ID. Other naming schemes need a different manifest step.
- The mock community sample (`Mock`) is processed together with the other samples through the Nextflow stages.
- Uploading files through the GitHub web page does not keep the executable permission, which is why the
  rerun steps include `chmod +x bin/make_manifest.sh`.

## Data and tools

- Data and example workflow: mothur MiSeq SOP, https://mothur.org/wiki/miseq_sop/
- Tools: Nextflow, Docker, QIIME 2 (2026.10 release), VSEARCH, SILVA 138.2
