#!/usr/bin/env bash
# =============================================================================
# 16S rRNA amplicon pipeline (final): QIIME 2 + VSEARCH, closed-reference OTUs (SILVA)
#
# STATUS: run end to end on the mothur MiSeq SOP test data with QIIME 2 2026.10
#         (conda env: rachis-qiime2-2026.10). Per-sample read counts (step 9)
#         were obtained separately with biom summarize-table, outside this script.
#
# Changes from the earlier v2 script (config values only, no logic changes):
#   THREADS=12                                  (was 8; set to this machine's nproc)
#   SILVA_REF_SEQS="refs/silva-97-seqs.qza"     (was a placeholder name)
#   SILVA_REF_TAX="refs/silva-97-tax.qza"       (was a placeholder name)
#   SILVA_REF_ALN="refs/silva-aln-V4.qza"       (was a placeholder name)
#   MAX_DEPTH=2260                              (was 2400)
#
# The three SILVA reference files are made BEFORE this script runs; the commands
# are in the runbook file (16s_pipeline_runbook.sh), section "Reference preparation".
#
# Usage:  bash 16s_otu_pipeline_final.sh all
#         bash 16s_otu_pipeline_final.sh prep import merge        (run chosen steps)
#
# Step numbers match the workflow table:
#   0 inputs | 1 import | 2-3 merge + ambiguous bases | 4 dereplicate | 5 chimeras
#   6 closed-ref OTUs | 7 taxonomy + filtering | 8 alignment + tree
#   9 table summary (biom summarize-table, run separately) | 10 taxonomic profile
#   11 alpha rarefaction | 12 export
#
# Not included: primer trimming (the MiSeq SOP reads contain no primers).
# =============================================================================
set -euo pipefail

# ============================= CONFIG ========================================
# Edit this block only. Nothing below should be hardcoded.
DATA_DIR="$HOME/16s_project/data/MiSeq_SOP"   # folder with the FASTQs and mothur design files
OUT="$HOME/16s_project/results"
THREADS=12

MANIFEST="$OUT/00_inputs/manifest.tsv"        # created by step 0
METADATA="$OUT/00_inputs/metadata.tsv"        # created by step 0

# SILVA reference files: all three must come from the SAME release (matching IDs)
# Paths are relative to the folder the script is run from (~/16s_project)
SILVA_REF_SEQS="refs/silva-97-seqs.qza"       # FeatureData[Sequence]        (97% V4 reference set)
SILVA_REF_TAX="refs/silva-97-tax.qza"         # FeatureData[Taxonomy]        (taxonomy for those IDs)
SILVA_REF_ALN="refs/silva-aln-V4.qza"         # FeatureData[AlignedSequence] (V4 reference alignment)

OTU_ID=0.97                                   # closed-reference identity threshold
STRAND="plus"                                 # use "both" if some reads may be reverse-oriented

# Lineages to remove. These labels were checked against the SILVA 138.2 taxonomy used here:
# chloroplast = o__Chloroplast, mitochondria = f__Mitochondria
EXCLUDE_TAXA="mitochondria,chloroplast,d__Archaea,d__Eukaryota"

# Optional: drop samples (e.g. the mock community) before diversity analysis.
# Leave empty to keep everything.
EXCLUDE_SAMPLES_WHERE="[time]='Mock'"

MAX_DEPTH=2260                                # just under the smallest sample (2263 reads)

mkdir -p "$OUT" "$OUT/00_inputs" "$OUT/export"

# ============================= 0. PREPARE INPUTS =============================
# Builds the manifest and the metadata file from the mothur design files.
# (For your own data, write manifest.tsv and metadata.tsv by hand instead.)
step_prep() {
  # gzip copies of the FASTQs so QIIME 2 import is safe; originals untouched
  mkdir -p "$OUT/00_inputs/fastq_gz"
  for f in "$DATA_DIR"/*.fastq; do
    gzip -c "$f" > "$OUT/00_inputs/fastq_gz/$(basename "$f").gz"
  done

  # Manifest: sample-id, forward path, reverse path (absolute paths required)
  awk -v d="$OUT/00_inputs/fastq_gz" 'BEGIN{OFS="\t";
        print "sample-id","forward-absolute-filepath","reverse-absolute-filepath"}
      {print $1, d"/"$2".gz", d"/"$3".gz"}' \
    "$DATA_DIR/stability.files" > "$MANIFEST"

  # Metadata: first column must be named sample-id
  # Joins mouse.time.design (group, time) with mouse.dpw.metadata (group, dpw)
  {
    printf "sample-id\ttime\tdpw\n"
    awk 'NR==FNR{if(FNR>1) dpw[$1]=$2; next}
         FNR>1{print $1"\t"$2"\t"dpw[$1]}' \
      "$DATA_DIR/mouse.dpw.metadata" "$DATA_DIR/mouse.time.design"
    printf "Mock\tMock\t\n"      # the mock community has no design entry
  } > "$METADATA"

  echo "Manifest rows: $(($(wc -l < "$MANIFEST") - 1)) | Metadata rows: $(($(wc -l < "$METADATA") - 1))"
}

# ============================= 1. IMPORT =====================================
step_import() {
  qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-format PairedEndFastqManifestPhred33V2 \
    --input-path "$MANIFEST" \
    --output-path "$OUT/01_demux.qza"
}

# ============================= 2-3. MERGE + AMBIGUOUS BASES ==================
# vsearch reverse-complements R2 internally when merging.
# --p-maxns 0 discards any contig containing an ambiguous base (step 3).
step_merge() {
  qiime vsearch merge-pairs \
    --i-demultiplexed-seqs "$OUT/01_demux.qza" \
    --p-maxns 0 \
    --p-threads "$THREADS" \
    --o-merged-sequences "$OUT/02_merged.qza" \
    --o-unmerged-sequences "$OUT/02_unmerged.qza"
}

# ============================= 4. DEREPLICATE ================================
# Per-sample dereplication; the output table keeps abundance and sample identity.
step_derep() {
  qiime vsearch dereplicate-sequences \
    --i-sequences "$OUT/02_merged.qza" \
    --o-dereplicated-table "$OUT/04_table.qza" \
    --o-dereplicated-sequences "$OUT/04_repseqs.qza"
}

# ============================= 5. CHIMERA REMOVAL ============================
step_chimera() {
  qiime vsearch uchime-denovo \
    --i-table "$OUT/04_table.qza" \
    --i-sequences "$OUT/04_repseqs.qza" \
    --output-dir "$OUT/05_chimera"

  qiime feature-table filter-features \
    --i-table "$OUT/04_table.qza" \
    --m-metadata-file "$OUT/05_chimera/nonchimeras.qza" \
    --o-filtered-table "$OUT/05_table_nochim.qza"

  qiime feature-table filter-seqs \
    --i-data "$OUT/04_repseqs.qza" \
    --m-metadata-file "$OUT/05_chimera/nonchimeras.qza" \
    --o-filtered-data "$OUT/05_repseqs_nochim.qza"
}

# ============================= 6. CLOSED-REFERENCE OTUs ======================
# Reads are matched to SILVA at 97%. Each OTU IS the SILVA reference sequence it
# matched, so OTU IDs = SILVA IDs. Reads with no match are discarded.
step_closedref() {
  qiime vsearch cluster-features-closed-reference \
    --i-table "$OUT/05_table_nochim.qza" \
    --i-sequences "$OUT/05_repseqs_nochim.qza" \
    --i-reference-sequences "$SILVA_REF_SEQS" \
    --p-perc-identity "$OTU_ID" \
    --p-strand "$STRAND" \
    --p-threads "$THREADS" \
    --o-clustered-table "$OUT/06_otu_table.qza" \
    --o-clustered-sequences "$OUT/06_otu_seqs.qza" \
    --o-unmatched-sequences "$OUT/06_unmatched.qza"
}

# ============================= 7. TAXONOMY + FILTERING =======================
# Because OTUs are SILVA sequences, taxonomy comes straight from the SILVA
# taxonomy file (no separate classifier needed). Then unwanted lineages are removed.
step_taxfilter() {
  qiime taxa filter-table \
    --i-table "$OUT/06_otu_table.qza" \
    --i-taxonomy "$SILVA_REF_TAX" \
    --p-exclude "$EXCLUDE_TAXA" \
    --o-filtered-table "$OUT/07_table_lineage.qza"

  # Optional: remove selected samples (e.g. the mock community)
  if [[ -n "$EXCLUDE_SAMPLES_WHERE" ]]; then
    qiime feature-table filter-samples \
      --i-table "$OUT/07_table_lineage.qza" \
      --m-metadata-file "$METADATA" \
      --p-where "$EXCLUDE_SAMPLES_WHERE" \
      --p-exclude-ids \
      --o-filtered-table "$OUT/07_table_final.qza"
  else
    cp "$OUT/07_table_lineage.qza" "$OUT/07_table_final.qza"
  fi

  # Keep only the OTU sequences that are still in the table
  qiime feature-table filter-seqs \
    --i-data "$OUT/06_otu_seqs.qza" \
    --i-table "$OUT/07_table_final.qza" \
    --o-filtered-data "$OUT/07_otu_seqs_final.qza"
}

# ============================= 8. ALIGNMENT + TREE ===========================
# OTU sequences are already in the SILVA reference alignment, so their aligned
# versions are pulled out of it (this worked: alignment had the same number of
# sequences as the final table has OTUs).
step_tree() {
  qiime feature-table filter-seqs \
    --i-data "$SILVA_REF_ALN" \
    --i-table "$OUT/07_table_final.qza" \
    --o-filtered-data "$OUT/08_aligned.qza"

  qiime alignment mask \
    --i-alignment "$OUT/08_aligned.qza" \
    --o-masked-alignment "$OUT/08_masked.qza"

  qiime phylogeny fasttree \
    --i-alignment "$OUT/08_masked.qza" \
    --p-n-threads "$THREADS" \
    --o-tree "$OUT/08_tree_unrooted.qza"

  qiime phylogeny midpoint-root \
    --i-tree "$OUT/08_tree_unrooted.qza" \
    --o-rooted-tree "$OUT/08_tree_rooted.qza"
}

# ============================= 10. TAXONOMIC PROFILE =========================
step_profile() {
  qiime taxa barplot \
    --i-table "$OUT/07_table_final.qza" \
    --i-taxonomy "$SILVA_REF_TAX" \
    --m-metadata-file "$METADATA" \
    --o-visualization "$OUT/10_taxa_barplot.qzv"
}

# ============================= 11. ALPHA RAREFACTION =========================
# OTU-based metrics only (they do not use the tree).
step_alpha() {
  qiime diversity alpha-rarefaction \
    --i-table "$OUT/07_table_final.qza" \
    --m-metadata-file "$METADATA" \
    --p-max-depth "$MAX_DEPTH" \
    --p-metrics observed_features \
    --p-metrics shannon \
    --o-visualization "$OUT/11_alpha_rarefaction.qzv"
}

# ============================= 12. EXPORT ====================================
# Plain files for plotting and imaging software.
# The barplot export contains per-taxonomic-level CSV tables.
step_export() {
  qiime tools export --input-path "$OUT/11_alpha_rarefaction.qzv" --output-path "$OUT/export/alpha_rarefaction"
  qiime tools export --input-path "$OUT/10_taxa_barplot.qzv"      --output-path "$OUT/export/taxa_barplot"
  qiime tools export --input-path "$OUT/07_table_final.qza"       --output-path "$OUT/export/otu_table"
  qiime tools export --input-path "$OUT/08_tree_rooted.qza"       --output-path "$OUT/export/tree"

  biom convert \
    -i "$OUT/export/otu_table/feature-table.biom" \
    -o "$OUT/export/otu_table/otu_table.tsv" \
    --to-tsv
}

# ============================= RUNNER ========================================
ALL_STEPS=(prep import merge derep chimera closedref taxfilter tree profile alpha export)

if [[ $# -eq 0 ]]; then
  echo "Usage: $0 all | <step> [<step> ...]"
  echo "Steps: ${ALL_STEPS[*]}"
  exit 1
fi

steps=("$@")
[[ "$1" == "all" ]] && steps=("${ALL_STEPS[@]}")

for s in "${steps[@]}"; do
  echo "=== Running step: $s ($(date '+%F %T')) ==="
  "step_$s"
done
