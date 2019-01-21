#!/usr/bin/env nextflow

params.db = "${baseDir}/db/hug"
params.in = "${baseDir}/data/many"
params.x = 'fasta'


genomes = Channel
    .fromPath("${params.in}/*.${params.x}")
    .map { file -> [file.baseName, file] }

hmmlist = Channel
    .fromPath("${params.db}/*.hmm")
    .collect()


process predict {
    tag "${genome_id}"
    publishDir "output/1/${genome_id}"

    input:
    set val(genome_id), file(genome_seq) from genomes

    output:
    set genome_id, "${genome_id}.genes.faa" into genes

    """
    prodigal -a ${genome_id}.genes.faa -i ${genome_seq}
    """
}

process search {
    tag "${genome_id}, ${hmm.baseName}"
    publishDir "output/1/${genome_id}"

    input:
    set val(genome_id), file(gene_seqs) from genes
    each hmm from hmmlist

    output:
    set val("${hmm.baseName}"), genome_id, file(gene_seqs), "${genome_id}.${hmm.baseName}.tbl" into marker_genes

    """
    hmmsearch --tblout ${genome_id}.${hmm.baseName}.tbl ${hmm} ${gene_seqs} > ${genome_id}.${hmm.baseName}.out
    """
}

process extract {
    tag "${genome_id}, ${marker_id}"
    publishDir "output/1/${genome_id}"

    input:
    set val(marker_id), val(genome_id), file(gene_seqs), file(tblout) from marker_genes

    output:
    set marker_id, "${genome_id}.${marker_id}.faa" into marker_genes_2

    """
    grep -v "^#" ${tblout} | awk 'FNR == 1 {print \$1}' > gene_id
    if [ -s gene_id ]
    then
      seqkit grep -f gene_id ${gene_seqs} | seqkit replace -p '.+' -r ${genome_id} > ${genome_id}.${marker_id}.faa
    else
      printf ">${genome_id}\nX\n" > ${genome_id}.${marker_id}.faa
    fi
    """
}

marker_genes_2.groupTuple(sort: true).set{ gene_blocks }

process align {
    tag "${marker_id}"
    publishDir "output/2"

    input:
    set val(marker_id), file(marker_seqs) from gene_blocks

    output:
    file "${marker_id}.aln" into raw_alignment

    """
    cat ${marker_seqs} | muscle > ${marker_id}.aln
    """
}

process concat {
    publishDir "output/3"

    input:
    file alignment_files from raw_alignment.toSortedList()

    output:
    file "concat.aln" into concat_alignment

    """
    seqkit concat -w 0 -t protein ${alignment_files} | seqkit replace -s -p X -r - > concat.aln
    """
}

process trim {
    publishDir "output/3"

    input:
    file alignment from concat_alignment

    output:
    file "trimmed.aln" into trimmed_alignment

    """
    trimal -in ${alignment} -gt 0.5 -keepseqs > trimmed.aln
    """
}

process build {
    publishDir "output/3"

    input:
    file alignment from trimmed_alignment

    output:
    file "tree.nwk" into tree

    """
    fasttree ${alignment} > tree.nwk
    """
}
