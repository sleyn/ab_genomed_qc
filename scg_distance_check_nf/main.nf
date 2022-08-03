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
        tuple val(pgfam), path("${pgfam}_temp.family"), val(atcc_sequence), emit: family_id_tbl_ch

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
        tuple val(pgfam), path(family_table_in), val(atcc_sequence)
    
    output:
        path "${pgfam}.tsv"
        tuple val(pgfam), path("${pgfam}.tsv"), val(atcc_sequence), emit: family_table_ch
    
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
    container 'semenleyn/ab-gen-qual-py-env:latest'
    containerOptions '-v "$(pwd):/temp" -w "/temp"'
    publishDir "${launchDir}/output/family_fasta/", mode: 'copy'

    input:
        tuple val(pgfam), path(family_table), val(atcc_sequence)

    output:
        path "${pgfam}.fasta"
        tuple val(pgfam), path("${pgfam}.fasta"), emit: family_fasta_ch

    script:
    """
    patric2fasta.py -t /temp/${family_table} -o ${pgfam}.fasta --ref_sequence ${atcc_sequence}
    """
}

/*
 * STEP 4. Align fasta sequences
 */

process ALIGN_FASTA {
    memory '4 GB'
    container 'semenleyn/mafft:latest'
    containerOptions '-v "$(pwd):/temp" -w "/temp"'
    publishDir "${launchDir}/output/family_aln/", mode: 'copy'

    input:
        tuple val(pgfam), path(family_fasta)

    output:
        path "${pgfam}.aln"
        tuple val(pgfam), path("${pgfam}.aln"), emit: family_aln_ch

    script:
    """
    mafft --thread -1 /temp/${family_fasta} > ${pgfam}.aln
    """
}

/*
 * STEP 5. Calculate distance matrix
 */

process CALCULATE_DISTANCE {
    memory '8 GB'
    cpus 1
    container 'semenleyn/ab-gen-qual-r-env:latest'
    containerOptions '-v "$(pwd):/temp" -w "/temp"'
    publishDir "${launchDir}/output/family_distance/", mode: 'copy'

    input:
        tuple val(pgfam), path(family_aln)

    output:
        path "${pgfam}.dist" , emit: family_dist

    script:
    """
    calculate_distances.R ${family_aln} ${pgfam}.dist
    """
}

/*
 * STEP 6. Combine distances for each PGFam, convert them to P-values assuming normal 
 * distribution of distances and filter by average p-value.
 */

process COMBINE_DISTANCES {
    container 'semenleyn/ab-gen-qual-r-env:latest'
    containerOptions '-v "$(pwd):/temp" -w "/temp"'
    publishDir "${launchDir}/output/family_combined/", mode: 'copy'

    input:
        path "*"

    output:
        path "distance_combined.tsv"
        path "distance_pvalue_combined.tsv"
        path "dist_mean_vs_sd.pdf"
        path "dist_pvalue.pdf"

    script:
    """
    compare_distances.R . .
    """
}

workflow {
    sc_gene_ch = Channel
                    .fromPath("${params.gene_table}")
                    .splitCsv(header: true)
                    .map { row -> tuple(row.Gene, row.PGFam, row.ATCC17978_sequence)}

    EMIT_FAMILY(sc_gene_ch)
    DOWNLOAD_FAMILY_MEMBER_TABLE(EMIT_FAMILY.out.family_id_tbl_ch)
    TBL_TO_FASTA(DOWNLOAD_FAMILY_MEMBER_TABLE.out.family_table_ch)
    ALIGN_FASTA(TBL_TO_FASTA.out.family_fasta_ch)
    CALCULATE_DISTANCE(ALIGN_FASTA.out.family_aln_ch)
    COMBINE_DISTANCES(CALCULATE_DISTANCE.out.family_dist.collect())
}

