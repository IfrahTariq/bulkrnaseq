version 1.0

import "../tasks/fastqc.wdl" as fastqc_task
import "../tasks/fastp.wdl" as fastp_task
import "../tasks/star.wdl" as star_task
import "../tasks/rnaseqc2.wdl" as rnaseqc2_task
import "../tasks/rsem.wdl" as rsem_task


workflow RNA_preprocessing_pipeline {

  input {

    File fastq1
    File fastq2
    String sample_id

    #annotation GTF
    File genes_gtf="gs://gtex-resources/GENCODE/gencode.v34.GRCh38.ERCC.genes.collapsed_only.gtf"

    #star_task index
    File star_index="gs://gtex-resources/STAR_genomes/STARv275a_genome_GRCh38_noALT_noHLA_noDecoy_ERCC_v34_oh100.tar.gz"

    #rsem index
    File rsem_reference="gs://gtex-resources/RSEM_references/rsem_reference_GRCh38_gencode34_ercc.tar.gz"

  }

  call fastqc_task.fastqc as raw_fastqc1 {
    input:
      seqFile=fastq1
  }

  call fastqc_task.fastqc as raw_fastqc2 {
    input:
      seqFile=fastq2
  }

  
  call star_task.star_workflow as star {
    input:
      prefix=sample_id,
      fastq1=fastq1,
      fastq2=fastq2,
      star_index=star_index
  }

  call rnaseqc2_task.rnaseqc2 as rnaseqc2 {
    input:
      bam_file=star.bam_file,
      genes_gtf=genes_gtf,
      sample_id=sample_id
  }

  call rsem_task.rsem as rsem {
    input:
      transcriptome_bam=star.transcriptome_bam,
      prefix=sample_id,
      rsem_reference=rsem_reference
  }


  output {
    #fastqc raw data
    File raw_htmlReport1 = raw_fastqc1.htmlReport
    File raw_reportZip1 = raw_fastqc1.reportZip
    File raw_htmlReport2 = raw_fastqc2.htmlReport
    File raw_reportZip2 = raw_fastqc2.reportZip

    #star
    File bam_file=star.bam_file
    File bam_index=star.bam_index
    File transcriptome_bam=star.transcriptome_bam
    File chimeric_junctions=star.chimeric_junctions
    File chimeric_bam_file=star.chimeric_bam_file
    File read_counts=star.read_counts
    File junctions=star.junctions
    File junctions_pass1=star.junctions_pass1
    Array[File] logs=star.logs

    #rnaseqc2
    File gene_tpm=rnaseqc2.gene_tpm
    File gene_counts=rnaseqc2.gene_counts
    File exon_counts=rnaseqc2.exon_counts
    File metrics=rnaseqc2.metrics
    File insertsize_distr=rnaseqc2.insertsize_distr
    File gc_content=rnaseqc2.gc_content

    #rsem
    File genes=rsem.genes
    File isoforms=rsem.isoforms
  }

}

