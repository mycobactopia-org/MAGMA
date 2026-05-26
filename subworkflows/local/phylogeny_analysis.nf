include { GATK_SELECT_VARIANTS as GATK_SELECT_VARIANTS_PHYLOGENY } from '../../modules/local/gatk/select_variants'
include { GATK_VARIANTS_TO_TABLE                                  } from '../../modules/local/gatk/variants_to_table'
include { SNPSITES                                                } from '../../modules/local/snpsites/snpsites'
include { SNPDISTS                                                } from '../../modules/local/snpdists/snpdists'
include { IQTREE                                                  } from '../../modules/local/iqtree/iqtree'


workflow PHYLOGENY_ANALYSIS {

    take:
    prefix_ch      // val: prefix string used in output filenames (e.g. 'IncComplex')
    arg_files_ch   // channel: paths to exclusion intervals / resource files
    vcf_ch         // channel: [ meta, tbi, vcf ]

    main:

    // Build the -XL (exclude-intervals) argument string from interval/list files
    args_ch = arg_files_ch
        .filter { it instanceof java.nio.file.Path && (it.getExtension() == "gz" || it.getExtension() == "list") }
        .reduce('') { a, b -> "$a -XL ${b.getName()}" }

    resources_files_ch = arg_files_ch
        .filter { it instanceof java.nio.file.Path && (it.getExtension() == "gz" || it.getExtension() == "list") }
        .collect()
        .ifEmpty([])

    resources_file_indexes_ch = arg_files_ch
        .filter { it instanceof java.nio.file.Path && it.getExtension() == "tbi" }
        .collect()
        .ifEmpty([])

    GATK_SELECT_VARIANTS_PHYLOGENY(
        'SNP',
        prefix_ch,
        vcf_ch,
        args_ch,
        resources_files_ch,
        resources_file_indexes_ch,
        params.ref_fasta,
        [params.ref_fasta_fai, params.ref_fasta_dict]
    )

    GATK_VARIANTS_TO_TABLE(prefix_ch, GATK_SELECT_VARIANTS_PHYLOGENY.out.variantsVcfTuple)

    SNPSITES(prefix_ch, GATK_VARIANTS_TO_TABLE.out)

    SNPDISTS(prefix_ch, SNPSITES.out)

    IQTREE(prefix_ch, SNPSITES.out)

    emit:
    snpsites_tree_tuple = SNPSITES.out.join(IQTREE.out.tree_tuple) // [ meta, fasta, treefile ]
    snp_dists_ch        = SNPDISTS.out.snp_dists_file
}
