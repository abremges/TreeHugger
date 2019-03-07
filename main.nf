#!/usr/bin/env nextflow
treehugger_version = 'v0.3.1'
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

// Marker gene set
params.custom = ''
if (params.custom) {
    params.marker = 'custom'
    db = "${params.custom}"
} else {
    params.marker = 'hug'
    db = "${baseDir}/db/hug"
    if (params.marker =~ /^a/) {
        if (params.marker =~ /arc/)
            db = "${baseDir}/db/amphora2_arc"
        else
            db = "${baseDir}/db/amphora2_bac"
    } else if (params.marker =~ /^p/)
        db = "${baseDir}/db/phylosift"
    else if (params.marker =~ /^s/)
        db = "${baseDir}/db/speci"
    else if (params.marker =~ /^c/)
        db = "${baseDir}/db/checkm"
    else if (params.marker =~ /^u/)
        db = "${baseDir}/db/ubcg"
    else if (params.marker =~ /^g/) {
        if (params.marker =~ /arc/)
            db = "${baseDir}/db/gtdb_arc122"
        else
            db = "${baseDir}/db/gtdb_bac120"
    }
}
datadir = Channel
    .fromPath("${db}", type: 'dir')
    .first()
hmmlist = Channel
    .fromPath("${db}/*.hmm")
    .collect()
println("[config] Marker gene set:\n${params.marker} - ${db}")

// Input folder (and file extension)
params.in = 'input'
params.x = 'fasta'
(genomes, input) = Channel
    .fromPath("${params.in}/*.${params.x}")
    .map { file -> [file.baseName, file] }
    .into(2)
println("[config] Input genomes: ")
input.println { it[0] + " - " + it[1] }

// TODO Optionally skip gene prediction --from-genes
params.prodigal = '-m -p meta'
process predict {
    tag "${genome_id}"

    input:
    set val(genome_id), file(genome_seq) from genomes

    output:
    set genome_id, "${genome_id}.genes.faa" into genes

    """
    prodigal ${params.prodigal} -a ${genome_id}.genes.faa -i ${genome_seq}
    """
}

params.hmmsearch = '-E 1e-10'
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
    hmmsearch ${params.hmmsearch} --tblout ${genome_id}.${hmm.baseName}.tbl ${hmm} ${gene_seqs}
    grep -v "^#" ${genome_id}.${hmm.baseName}.tbl | awk 'FNR == 1 {print \$1}' > ${genome_id}.${hmm.baseName}.gid
    if [ -s ${genome_id}.${hmm.baseName}.gid ]
    then
      seqkit grep -f ${genome_id}.${hmm.baseName}.gid ${gene_seqs} | seqkit replace -p '.+' -r "${genome_id} ${hmm.baseName}" > ${genome_id}.${hmm.baseName}.faa
    else
      printf ">${genome_id} ${hmm.baseName}\nX\n" > ${genome_id}.${hmm.baseName}.faa
    fi
    """
}

params.muscle = '-maxiters 16' // -maxiters 8
process align {
    tag "${marker_id}"
    publishDir "output/2"

    input:
    set val(marker_id), file(marker_seqs) from marker_genes.groupTuple()

    output:
    file "${marker_id}.aln" into alignment

    """
    cat ${marker_seqs} | muscle ${params.muscle} > ${marker_id}.aln
    """
}

params.trimal = '-automated1' // -keepseqs
params.fasttree = '-lg -gamma'
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
    trimal ${params.trimal} -in concat.aln > trimmed.aln
    fasttree ${params.fasttree} trimmed.aln > tree.nwk
    """
}
