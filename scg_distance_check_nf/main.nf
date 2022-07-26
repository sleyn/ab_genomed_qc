nextflow.enable.dsl = 2

/*
 * Define deault parameters:
 * 1. CSV file defining single copy marker genes from CheckM. Columns:
 *    - Gene: Gene name
 *    - Function: Gene function
 *    - PGFam: Protein Global Family from PATRIC
 *    - ATCC17978_sequence: Protein sequence from the reference ATCC17978 genome
 */

params.gene_table = "sc_genes.csv"

/*
 * Processes:
 */

/*
 * STEP 1: Make temporary file with family ID
 */

 process EMIT_FAMILY {
    memory '2 GB'
    container 'semenleyn/patric_cli_1.039_ubuntu_20.04'

    input:
        tuple val(gene), val(pgfam), val(atcc_sequence)
    
    output:
        path "${pgfam}_temp.family", emit: family_id_tbl_ch

    script:
    """
    p3-echo -t family ${pgfam} > ${pgfam}_temp.family
    """
 }

 /*
  * STEP 2: Download a table with each protein family member:
  *     - genome_id
  *     - genome_name
  *     - patric_id
  *     - plfam_id
  *     - pgfam_id
  *     - product,gene
  */

process DOWNLOAD_FAMILY_MEMBER_TABLE {
    memory '2 GB'
    container 'semenleyn/patric_cli_1.039_ubuntu_20.04'
    containerOptions '-v "$(pwd):/temp"'
    publishDir "${launchDir}/output/family_tables/", mode: 'copy'

    input:
        tuple val(gene), val(pgfam), val(atcc_sequence)
        path family_table_in
    
    output:
        path "${pgfam}.tsv", emit: family_table_ch
    
    script:
    """
    p3-get-family-features \
        --input /temp/${family_table_in} \
        --ftype=global \
		--equal "genome_name,Acinetobacter baumannii" \
        --attr genome_id,genome_name,patric_id,plfam_id,pgfam_id,product,gene,aa_sequence > ${pgfam}.tsv
    """
}

/*
 * STEP 3: Convert table of proteins to the FASTA format.
 */

process TBL_TO_FASTA {
    container 'semenleyn/patric2fasta:latest'
    containerOptions '-v "$(pwd):/temp"'
    publishDir "${launchDir}/output/family_fasta/", mode: 'copy'

    input:
        tuple val(gene), val(pgfam), val(atcc_sequence)
        path family_table

    output:
        path "${pgfam}.fasta", emit: family_fasta_ch

    script:
    """
    patric2fasta.py -t /temp/${family_table} -o ${pgfam}.fasta --ref_sequence ${atcc_sequence}
    """
}

workflow {
    sc_gene_ch = Channel
                    .fromPath("${params.gene_table}")
                    .splitCsv(header: true)
                    .map { row -> tuple(row.Gene, row.PGFam, row.ATCC17978_sequence)}

    
    EMIT_FAMILY(sc_gene_ch)
    DOWNLOAD_FAMILY_MEMBER_TABLE(
        sc_gene_ch,
        EMIT_FAMILY.out.family_id_tbl_ch
        )
    TBL_TO_FASTA(
        sc_gene_ch, 
        DOWNLOAD_FAMILY_MEMBER_TABLE.out.family_table_ch
        )
//    ALIGN_FASTA()
//    GENERATE_DISTANCE()
}