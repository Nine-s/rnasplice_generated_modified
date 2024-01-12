nextflow.enable.dsl = 2


include { FASTQC  } from '../../../ninon/description_prototype/rnasplice_modules/fastqc.nf'
include { TRIMGALORE  } from '../../../ninon/description_prototype/rnasplice_modules/trimgalore.nf'
include { SALMON_GENOMEGENERATE  } from '../../../ninon/description_prototype/rnasplice_modules/salmon_genome_generate.nf'
include { SALMON_QUANT  } from '../../../ninon/description_prototype/rnasplice_modules/salmon.nf'
include { HISAT2_INDEX  } from '/path/Hisat2_index'
include { HISAT2_ALIGN  } from '/path/Hisat2'
include { SAMTOOLS  } from '../../../ninon/description_prototype/rnasplice_modules/samtools.nf'
include { CUSTOM_GETCHROMSIZES  } from '../../../ninon/description_prototype/rnasplice_modules/getchromsizes.nf'
include { BEDTOOLS_GENOMECOV  } from '../../../ninon/description_prototype/rnasplice_modules/bedtoolsgenomecov.nf'
include { BEDCLIP as BEDCLIP_FORWARD; BEDCLIP as BEDCLIP_REVERSE } from '../../../ninon/description_prototype/rnasplice_modules/bedclip.nf'
include { BEDGRAPHTOBIGWIG as BEDGRAPH_TO_BIGWIG_FORWARD; BEDGRAPHTOBIGWIG as BEDGRAPH_TO_BIGWIG_REVERSE } from '../../../ninon/description_prototype/rnasplice_modules/bedgraphtobigwig.nf'
include { DEXSEQ_ANNOTATION  } from '../../../ninon/description_prototype/rnasplice_modules/dexseq_annotation.nf'
include { DEXSEQ_COUNT  } from '../../../ninon/description_prototype/rnasplice_modules/dexseq_count.nf'
include { MERGE_RESULTS_DEXSEQ  } from '../../../ninon/description_prototype/rnasplice_modules/merge_results_dexseq.nf'
include { DEXSEQ_EXON  } from '../../../ninon/description_prototype/rnasplice_modules/dexseq_exon.nf'
include { GFFREAD_TX2GENE  } from '../../../ninon/description_prototype/rnasplice_modules/gffread_tx2gene.nf'
include { MERGE_RESULTS_SALMON  } from '../../../ninon/description_prototype/rnasplice_modules/merge_results.nf'
include { TXIMPORT  } from '../../../ninon/description_prototype/rnasplice_modules/tximport.nf'
include { DRIMSEQ_FILTER  } from '../../../ninon/description_prototype/rnasplice_modules/drimseq_filter.nf'
include { DEXSEQ_DTU  } from '../../../ninon/description_prototype/rnasplice_modules/dexseq_dtu.nf'
include { MULTIQC  } from '../../../ninon/description_prototype/rnasplice_modules/multiqc.nf'
include { FASTQSPLIT as FASTQSPLIT_HISAT2 } from '/path/FASTQSPLIT.nf'
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_HISAT2 } from '/path/SAMTOOLS.nf'

workflow{
        read_pairs_ch = Channel
            .fromPath( params.csv_input )
            .splitCsv(header: true, sep: ',')
            .map {row -> tuple(row.sample, [row.path_r1, row.path_r2], row.condition)}
            .view()
        
SALMON_GENOMEGENERATE(params.genome, params.transcripts_fasta)
GFFREAD_TX2GENE(params.annotation_gtf)
FASTQC(read_pairs_ch)
DEXSEQ_ANNOTATION(params.annotation_gtf)
TRIMGALORE(read_pairs_ch)
SALMON_QUANT(TRIMGALORE.out.preprocessed_reads, SALMON_GENOMEGENERATE.out.index)
HISAT2_INDEX(params.genome, params.annotation_gtf)
CUSTOM_GETCHROMSIZES(params.genome)
MERGE_RESULTS_SALMON(SALMON_QUANT.out.transcripts.collect())
FASTQSPLIT_HISAT2(TRIMGALORE.out.preprocessed_reads)
TXIMPORT(MERGE_RESULTS_SALMON.out.gathered_bam, GFFREAD_TX2GENE.out.tx2gene)
DRIMSEQ_FILTER(TXIMPORT.out.txi_dtu, TXIMPORT.out.tximport_tx2gene, params.csv_input, params.min_samps_gene_expr, params.min_samps_feature_expr, params.min_samps_feature_prop, params.min_feature_expr, params.min_feature_prop, params.min_gene_expr)
HISAT2_ALIGN(FASTQSPLIT_HISAT2.out.split_reads, HISAT2_INDEX.out.index, params.annotation_gtf)
MULTIQC(SALMON_QUANT.out.json_info.collect(), TRIMGALORE.out.log.collect(), HISAT2_ALIGN.out.log_final.collect(), FASTQC.out.zip.collect())
SAMTOOLS_MERGE_HISAT2(HISAT2_ALIGN.out.sam.collect())
SAMTOOLS(SAMTOOLS_MERGE_HISAT2.out.merged)
DEXSEQ_COUNT(SAMTOOLS.out.bam, DEXSEQ_ANNOTATION.out.gff, params.alignment_quality)
DEXSEQ_DTU(DRIMSEQ_FILTER.out.drimseq_samples_tsv, DRIMSEQ_FILTER.out.drimseq_counts_tsv, params.csv_contrastsheet, params.n_dexseq_plot)
MERGE_RESULTS_DEXSEQ(DEXSEQ_COUNT.out.dexseq_clean_txt.collect())
BEDTOOLS_GENOMECOV(SAMTOOLS.out.bam)
BEDCLIP_FORWARD(BEDTOOLS_GENOMECOV.out.bedgraph_forward, CUSTOM_GETCHROMSIZES.out.sizes)
BEDGRAPH_TO_BIGWIG_FORWARD(BEDCLIP_FORWARD.out.bedgraph, CUSTOM_GETCHROMSIZES.out.sizes)
BEDCLIP_REVERSE(BEDTOOLS_GENOMECOV.out.bedgraph_reverse, CUSTOM_GETCHROMSIZES.out.sizes)
DEXSEQ_EXON(MERGE_RESULTS_DEXSEQ.out.clean_counts, DEXSEQ_ANNOTATION.out.gff, params.csv_input, params.csv_contrastsheet, params.n_dexseq_plot)
BEDGRAPH_TO_BIGWIG_REVERSE(BEDCLIP_REVERSE.out.bedgraph, CUSTOM_GETCHROMSIZES.out.sizes)

}