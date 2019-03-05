#!/usr/bin/env nextflow
treehugger_version = 'v0.3'
println """
▄▄▄█████▓ ██▀███  ▓█████ ▓█████  ██░ ██  █    ██   ▄████   ▄████ ▓█████  ██▀███
▓  ██▒ ▓▒▓██ ▒ ██▒▓█   ▀ ▓█   ▀ ▓██░ ██▒ ██  ▓██▒ ██▒ ▀█▒ ██▒ ▀█▒▓█   ▀ ▓██ ▒ ██▒
▒ ▓██░ ▒░▓██ ░▄█ ▒▒███   ▒███   ▒██▀▀██░▓██  ▒██░▒██░▄▄▄░▒██░▄▄▄░▒███   ▓██ ░▄█ ▒
░ ▓██▓ ░ ▒██▀▀█▄  ▒▓█  ▄ ▒▓█  ▄ ░▓█ ░██ ▓▓█  ░██░░▓█  ██▓░▓█  ██▓▒▓█  ▄ ▒██▀▀█▄
  ▒██▒ ░ ░██▓ ▒██▒░▒████▒░▒████▒░▓█▒░██▓▒▒█████▓ ░▒▓███▀▒░▒▓███▀▒░▒████▒░██▓ ▒██▒
  ▒ ░░   ░ ▒▓ ░▒▓░░░ ▒░ ░░░ ▒░ ░ ▒ ░░▒░▒░▒▓▒ ▒ ▒  ░▒   ▒  ░▒   ▒ ░░ ▒░ ░░ ▒▓ ░▒▓░
    ░      ░▒ ░ ▒░ ░ ░  ░ ░ ░  ░ ▒ ░▒░ ░░░▒░ ░ ░   ░   ░   ░   ░  ░ ░  ░  ░▒ ░ ▒░
  ░        ░░   ░    ░      ░    ░  ░░ ░ ░░░ ░ ░ ░ ░   ░ ░ ░   ░    ░     ░░   ░
            ░        ░  ░   ░  ░ ░  ░  ░   ░           ░       ░    ░  ░   ░
                                                ${treehugger_version?:''}
"""

// TODO Expose command-line options for each tool in TreeHugger

// TODO Provide flags for commonly used marker gene sets
// amphora2_arc amphora2_bac checkm gtdb_arc122 gtdb_bac120 hug phylosift speci ubcg
params.db = "${baseDir}/db/hug"
params.in = "${baseDir}/data/one"
params.x = 'fasta'


genomes = Channel
    .fromPath("${params.in}/*.${params.x}")
    .map { file -> [file.baseName, file] }

datadir = Channel
    .fromPath("${params.db}", type: 'dir')
    .first()

hmmlist = Channel
    .fromPath("${params.db}/*.hmm")
    .collect()

// TODO Optionally skip gene prediction --from-genes
process predict {
    tag "${genome_id}"

    input:
    set val(genome_id), file(genome_seq) from genomes

    output:
    set genome_id, "${genome_id}.genes.faa" into genes

    """
    prodigal -m -a ${genome_id}.genes.faa -i ${genome_seq}
    """ // -p meta ?
}

process search {
    tag "${genome_id}, ${hmm.baseName}"
    publishDir "output/1/${genome_id}"

    input:
    set val(genome_id), file(gene_seqs) from genes
    each hmm from hmmlist
    file datadir // Docker workaround

    output:
    set val("${hmm.baseName}"), "${genome_id}.${hmm.baseName}.faa" into marker_genes

    """
    hmmsearch -E 1e-10 --tblout ${genome_id}.${hmm.baseName}.tbl ${hmm} ${gene_seqs}
    grep -v "^#" ${genome_id}.${hmm.baseName}.tbl | awk 'FNR == 1 {print \$1}' > ${genome_id}.${hmm.baseName}.gid
    if [ -s ${genome_id}.${hmm.baseName}.gid ]
    then
      seqkit grep -f ${genome_id}.${hmm.baseName}.gid ${gene_seqs} | seqkit replace -p '.+' -r "${genome_id} ${hmm.baseName}" > ${genome_id}.${hmm.baseName}.faa
    else
      printf ">${genome_id} ${hmm.baseName}\nX\n" > ${genome_id}.${hmm.baseName}.faa
    fi
    """
}

process align {
    tag "${marker_id}"
    publishDir "output/2"

    input:
    set val(marker_id), file(marker_seqs) from marker_genes.groupTuple()

    output:
    file "${marker_id}.aln" into alignment

    """
    cat ${marker_seqs} | muscle > ${marker_id}.aln
    """ // -maxiters 8 ?
}

process build {
    publishDir "output/3"

    input:
    file alignment_files from alignment.toSortedList( { a, b -> a.baseName <=> b.baseName } )

    output:
    file "concat.aln"
    file "trimmed.aln"
    file "tree.nwk"

    """
    seqkit concat -w 0 -t protein ${alignment_files} | seqkit replace -s -p X -r - > concat.aln
    trimal -automated1 -in concat.aln > trimmed.aln
    fasttree -lg -gamma trimmed.aln > tree.nwk
    """ // trimal -keepseqs ?
}
