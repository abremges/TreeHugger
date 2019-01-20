#!/usr/bin/env nextflow

params.db = "${baseDir}/db/hug"
params.in = "${baseDir}/data/one"
params.x = 'fasta'


genomes = Channel
    .fromPath("${params.in}/*.${params.x}")
    .map { file -> [file.baseName, file] }

hmmlist = Channel
    .fromPath("${params.db}/*.hmm")
    .collect()


process prodigal {
    tag "${id}"
    publishDir "output/prodigal/${id}"

    input:
    set val(id), file(genome) from genomes

    output:
    set id, "${id}.genes.faa" into genes

    """
    prodigal -a ${id}.genes.faa -i ${genome}
    """
}


process hmmsearch {
    tag "${id}, ${hmm.baseName}"
    publishDir "output/hmmsearch/${id}"

    input:
    set val(id), file(genes) from genes
    each hmm from hmmlist

    output:
    set val("${hmm.baseName}"), id, file(genes), "${id}.${hmm.baseName}.gid" into marker_genes

    """
    hmmsearch --tblout ${id}.${hmm.baseName}.tbl ${hmm} ${genes} > ${id}.${hmm.baseName}.out
    grep -v "^#" ${id}.${hmm.baseName}.tbl | awk 'FNR == 1 {print \$1}' > ${id}.${hmm.baseName}.gid
    """
}

process extract {
    tag "${genome_id}, ${marker_id}"
    publishDir "output/extract/${genome_id}"

    input:
    set val(marker_id), val(genome_id), file(genes), file(gene_id) from marker_genes

    output:
    set marker_id, genome_id, "${genome_id}.${marker_id}.faa" into marker_genes_2

    """
    if [ -s ${gene_id} ]
    then
      seqkit seq -i ${genes} | seqkit grep -p \$(< ${gene_id}) | seqkit replace -p \$(< ${gene_id}) -r ${genome_id} > ${genome_id}.${marker_id}.faa
    else
      printf ">${genome_id}\nX\n" > ${genome_id}.${marker_id}.faa
    fi
    """
}

marker_genes_2.groupTuple(sort: true).set{ gene_blocks }

process muscle {
    tag "${marker_id}"
    publishDir "output/muscle"

    input:
    set val(marker_id), val(genome_id), file(gene_sequence) from gene_blocks

    output:
    set marker_id, genome_id, "${marker_id}.muscle.fasta" into raw_alignment

    """
    cat ${gene_sequence} | muscle > ${marker_id}.muscle.fasta
    """
}

process trimal {
    tag "${marker_id}"
    publishDir "output/trimal"

    input:
    set val(marker_id), val(genome_id), file(aligned_sequence) from raw_alignment

    output:
    file "${marker_id}.trimal.fasta" into trimmed_alignment

    """
    trimal -in ${aligned_sequence} -gt 0.5 -keepseqs > ${marker_id}.trimal.fasta
    """
}

trimmed_alignment.toSortedList().view()

// process concat {
//     tag "..."
//     publishDir "output/concat"
//
//     input:
//     set file(alignments) from trimmed_alignment.collect()
//
//     output:
//     file "alignment.fasta" into concat_alignment
//
//     """
//     for summary in ${alignments}
//     do
//         echo seqkit concat \${summary} >> alignment.fasta
//     done
//
//     """
// }
