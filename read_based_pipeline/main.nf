#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
===============================================================================
Soil Metagenome Analysis Pipeline
===============================================================================
*/

workflow {

    // Gather input fasta files
    read_ch_input = Channel.fromFilePairs(params.input_files_input)
    read_ch_treatment = Channel.fromFilePairs(params.input_files_treatment)

    genome_ch = Channel.fromPath(params.input_genomes)
    .map { file -> tuple(file.baseName, file) }

    plasmid_ch = Channel.fromPath(params.input_plasmids)
    .map { file -> tuple(file.baseName, file) }

    // Combined fasta file input channel
    read_ch = read_ch_input.concat(read_ch_treatment)

    // Define Gene databases 
    def CARD_tup = tuple('CARD', file(params.CARD_DB))
    def Res_tup = tuple('ResFinderFG', file(params.RES_DB))
    def MGE_tup = tuple('mobileOG', file(params.MGE_DB))

    // Combined ARG fasta file
    Create_ARG_DB(CARD_tup, Res_tup)

    // Channel with all databases
    db_ch = arg_db_ch.mix(Channel.of(MGE_tup, Res_tup, MGE_tup))

    // Screen for genes in db
    Diamond_index(db_ch)
    Diamond_mapping_genes(read_ch.combine(Diamond_index.out.db))
    Diamond_mapping_plasmids(plasmid_ch.combine(Diamond_index.out.db))

    // Screen for point mutations
    Mumame(read_ch, params.mumame_db)   

    //cluster databases
    Cluster_db(db_ch)

    // Screen for genes in GI
    input_channel = Channel.fromPath(params.island_viewer_data)
        .splitCsv(sep: '\t', header: true)
        .map { row -> row.accession }
        .unique()

    download_genomes(input_channel)
    extract_genomic_islands(download_genomes.out.genomes, params.island_viewer_data)
    screen_GI(extract_genomic_islands.out.genomic_islands, Diamond_index.out.db)
}

process Create_ARG_DB {
    conda 'cd-hit=4.8.1'
    publishDir "${params.directory_out}/CD_HIT/", mode:"copy"
    cpus 1

    input:
    tuple val(type), path(CARD)
    tuple val(type), path(RF)

    output:
    tuple val("ARG_db"), path("ARG_db_clstr_100.faa"), emit:arg_db
    path "ARG_db_clstr_100.clstr", emit:cluster
   
    script:
    """
    cat ${CARD} ${RF} > comb_ARG.aa
    
    # Cluster and  create aa file
    cd-hit -i comb_ARG.aa -o ARG_db_clstr_100 -n 5 -c 1  -d 200 -T ${task.cpus} 
    mv ARG_db_clstr_100  ARG_db_clstr_100.faa

    # replace spaces in ResFinderFG header names
    sed -i 's/ /_/g' ARG_db_clstr_100.faa
    """
}

process Cluster_db {
    conda 'cd-hit=4.8.1'
    publishDir "${params.directory_out}/CD_HIT/", mode:"copy"
    cpus 1

    input:
    tuple val(type), path(db)

    output:
    tuple val("type"), path("${type}_95*"), emit:clusters
   
    script:
    """
    cd-hit -i ${db} -o ${type}_95  -n 5 -c 0.95  -d 200 -T ${task.cpus} 
    mv ${type}_95  ${type}_95.faa
    """
}

process Diamond_index {
    cpus 1
    conda "diamond=2.1.9"

    input:
    tuple val(type), path(db)

    output:
    tuple val(type), path("${type}.dmnd"), emit: db

    shell:
    """
    # Create index
    diamond makedb --in ${db} -d ${type}
    """
}

process Diamond_mapping_genes {
    cpus 10
    conda "diamond=2.1.9"
    publishDir "${params.directory_out}/Diamond/ARG/", mode: 'copy'

    input:
    tuple val(sample_id), path(reads), val(type), path(db)

    output:
    path "${sample_id}_${type}.tsv" , emit: diamond_out

    shell:
    """
    diamond blastx --db ${db} --query ${reads[0]} --out ${sample_id}_${type}.tsv --threads ${task.cpus} --id 80 --query-cover 80 --outfmt 6 -k0
    """
}

process Mumame {
    conda 'bengtssonpalme::mumame' 
    publishDir "${params.directory_out}/Mumame/", mode:"copy"
    cpus 20

    input:
    tuple val(sample_id), path(reads) 
    val db

    output:
    path "${sample_id}.per-mutation.txt", emit: mutation_txt
    path "${sample_id}.table.txt", emit: tablels
    path "${sample_id}*", emit: other

    script:
    """
    mumame -d ${db} -o ${sample_id} *.fq.gz
    """
}

process Diamond_mapping_plasmids {
    cpus 10
    conda "diamond=2.1.9 prodigal=2.6.3"
    publishDir "${params.directory_out}/Diamond/ARG_plasmids/", mode: 'copy'

    input:
    tuple val(sample_id), path(plasmid), val(type), path(db)

    output:
    path "${sample_id}_${type}.tsv" , emit: diamond_out

    shell:
    """
    gunzip  -c ${plasmid} > file.fna
    prodigal -i file.fna -a protein.translations.faa -p single
    diamond blastp --db ${db} --query protein.translations.faa --out ${sample_id}_${type}.tsv --threads ${task.cpus} --id 90 --subject-cover 80 --query-cover 80  --outfmt 6 -k0
    rm file.fna
    """
}

process download_genomes {
    cpus 1
    conda "bioconda::entrez-direct=22.4 gzip"
    errorStrategy 'ignore'
    
    input:
    val(accession)

    output:
    tuple val(accession), path("${accession}.fasta.gz"), emit: genomes

    script:
    """
     #!/bin/bash
    set -e
    
    max_attempts=3
    attempt=1
    
    download_and_check() {
        efetch -db nucleotide -id ${accession} -format fasta > ${accession}.fasta
        
        if grep -q '^>' ${accession}.fasta && grep -q '^[ACGTNacgtn]' ${accession}.fasta; then
            gzip ${accession}.fasta
            return 0
        else
            rm ${accession}.fasta
            return 1
        fi
    }
    
    until download_and_check || [ \$attempt -eq \$max_attempts ]
    do
        attempt=\$((attempt + 1))
    done
    
    if [ \$attempt -eq \$max_attempts ] && [ ! -f "${accession}.fasta.gz" ]; then
        exit 1
    fi
    """
}

process extract_genomic_islands {
    cpus 1
    conda "python=3.9 biopython=1.81 pandas=2.2.0 "
    publishDir "${params.directory_out}/Diamond/GI_genomes/"
    
    input:
    tuple val(accession), path(genome_file)
    path input_file

    output:
    tuple val(accession), path("${accession}_islands.fasta.gz"), emit: genomic_islands

    script:
    """
    #!/usr/bin/env python3
    import gzip
    from Bio import SeqIO
    import pandas as pd

    # Read the input TSV file (assuming it's not gzipped)
    df = pd.read_csv("${input_file}", sep='\t')

    # Filter the dataframe for the current accession
    islands = df[df['accession'] == "${accession}"]

    # Open the input gzipped genome file and output gzipped fasta file
    with gzip.open("${genome_file}", "rt") as input_handle, gzip.open("${accession}_islands.fasta.gz", "wt") as output_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            for idx, row in islands.iterrows():
                start = int(row['start'])
                end = int(row['end'])
                island = record.seq[start:end+1]
                SeqIO.write(SeqIO.SeqRecord(island, id=f"${accession}_island_{idx}", description=f"start:{start} end:{end}"), output_handle, "fasta")
    """
}

process screen_GI {
    cpus 1
    conda "diamond=2.1.9 prodigal=2.6.3"
    publishDir "${params.directory_out}/Diamond/ARG_GI/", mode: 'copy'
    
    input:
    tuple val(accession), path(island_file)
    tuple val(type), path(db)

    output:
    path "${accession}_GI.tsv", emit: diamond_out

    script:
    """
    gunzip  -c  ${island_file} > islands.fasta
    prodigal -i islands.fasta -a protein.translations.faa -p meta
    diamond blastp --db ${db} --query protein.translations.faa --out ${accession}_GI.tsv --threads ${task.cpus} --id 90 --subject-cover 80 --query-cover 80  --outfmt 6 -k0
    rm islands.fasta
    """
}