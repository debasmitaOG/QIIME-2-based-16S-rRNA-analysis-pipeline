params.reads = 'data/*_R{1,2}_001.fastq'
params.maxns = 0

process IMPORT {
    input:
    path reads

    output:
    path 'demux.qza'

    script:
    """
    make_manifest.sh manifest.tsv
    qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-format PairedEndFastqManifestPhred33V2 \
        --input-path manifest.tsv \
        --output-path demux.qza
    """
}

process MERGE {
    cpus 4

    input:
    path demux

    output:
    path 'merged.qza'

    script:
    """
    qiime vsearch merge-pairs \
        --i-demultiplexed-seqs ${demux} \
        --p-maxns ${params.maxns} \
        --p-threads ${task.cpus} \
        --o-merged-sequences merged.qza \
        --o-unmerged-sequences unmerged.qza
    """
}

process DEREP {
    input:
    path merged

    output:
    path 'derep_table.qza', emit: table
    path 'derep_seqs.qza', emit: seqs

    script:
    """
    qiime vsearch dereplicate-sequences \
        --i-sequences ${merged} \
        --o-dereplicated-table derep_table.qza \
        --o-dereplicated-sequences derep_seqs.qza
    """
}

process CHIMERA {
    input:
    path table
    path seqs

    output:
    path 'nochim_table.qza', emit: table
    path 'nochim_seqs.qza', emit: seqs

    script:
    """
    qiime vsearch uchime-denovo \
        --i-table ${table} \
        --i-sequences ${seqs} \
        --output-dir chimera

    qiime feature-table filter-features \
        --i-table ${table} \
        --m-metadata-file chimera/nonchimeras.qza \
        --o-filtered-table nochim_table.qza

    qiime feature-table filter-seqs \
        --i-data ${seqs} \
        --m-metadata-file chimera/nonchimeras.qza \
        --o-filtered-data nochim_seqs.qza
    """
}

workflow {
    main:
    ch_reads = channel.fromPath(params.reads, checkIfExists: true).collect()
    ch_demux = IMPORT(ch_reads)
    ch_merged = MERGE(ch_demux)
    DEREP(ch_merged)
    CHIMERA(DEREP.out.table, DEREP.out.seqs)

    publish:
    demux = ch_demux
    merged = ch_merged
    derep_table = DEREP.out.table
    derep_seqs = DEREP.out.seqs
    nochim_table = CHIMERA.out.table
    nochim_seqs = CHIMERA.out.seqs
}

output {
    demux {
        path 'import'
    }
    merged {
        path 'merge'
    }
    derep_table {
        path 'derep'
    }
    derep_seqs {
        path 'derep'
    }
    nochim_table {
        path 'chimera'
    }
    nochim_seqs {
        path 'chimera'
    }
}
