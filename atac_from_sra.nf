#!/usr/bin/env nextflow
// Copyright (C) 2020 Kade Pettie

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// Nextflow pipeline to map ATAC-seq reads from SRA fastqs.

out = file(params.outdir)
if( !out.exists() ) out.mkdir()

import org.jsoup.Jsoup;
import org.jsoup.helper.Validate;
import org.jsoup.nodes.Document;
import org.jsoup.nodes.Element;
import org.jsoup.select.Elements;

///// SCRAPE NCBI FOR SRA SAMPLE METADATA /////

if (params.fastq_local_glob) {

  Channel.fromPath( params.fastq_local_glob )
    .map{ itg = ( it.getName() =~ /([a-z]{3})([1-2]{1})_(R[1-2]{1})\.fastq\.gz/ )[0]
          return [ '_', itg[1].toUpperCase(), "rep${itg[2]}", it ] }
    .groupTuple(by: [0,1,2])
    .into{ SRA; COUNT_RAW }

} else {

  sampkey = []

  Document doc = Jsoup.connect( params.urlfull ).get();
  Elements links = doc.select("a[href]");
  for (Element link : links) {
      if ( link.attr("href") =~ /sra\/SRX/ ) {
          Document doc2 = Jsoup.connect( link.attr("abs:href") ).get()
          Elements links2 = doc2.select("a[href]")
          for (Element link2 : links2) {
              if ( link2.text() =~ /SRR/ ) {
                  ax = link2.text()
              }
          }
          (full, pop, rep) = ( link.text() =~ /.+([A-Z]{3}).+(\d)$/ )[0]
          sampkey.add( [ax, pop, "rep$rep"] )
      }
  }

  Channel.fromList(sampkey)
      .combine( Channel.fromSRA( params.sra_id, apiKey: params.sra_api_key ), by: 0 )
      .into{ SRA; COUNT_RAW }

}

///////// FASTQS FROM SRA //////////

process cutadapt {
  label 'cutadapt'
  publishDir "${params.outdir}/fastq/"

  when:
  params.mode =~ /(all)/

  input:
  tuple id, pop, rep, file(fq) from SRA

  output:
  tuple pop, rep, path("*.fastq1.cutadapt.gz"), path("*.fastq2.cutadapt.gz") into BAM_INIT, COUNT_CUTADAPT

  script:
  """
  cutadapt \
  --cores ${params.mapcores} \
  -e 0.20 \
  -a CTGTCTCTTATACACATCT \
  -A CTGTCTCTTATACACATCT \
  -m 5 \
  -o ${pop}_${rep}.fastq1.cutadapt.gz \
  -p ${pop}_${rep}.fastq2.cutadapt.gz \
  ${fq[0]} \
  ${fq[1]}
  """

}

process map_initial {
  label 'map'
  publishDir "${params.outdir}/bams/initial/"

  when:
  params.mode =~ /(all)/

  input:
  tuple pop, rep, path(fq1), path(fq2) from BAM_INIT

  output:
  tuple pop, rep, path(outbam) into RMDUP, COUNT_INITIAL

  script:
  outunsort = "${pop}_${rep}.initial.uncoordsort.bam"
  outbam = "${pop}_${rep}.initial.bam"
  """
  bowtie2 \
  -p ${params.mapcores} \
  ${params.bt2_opts} \
  -x ${params.bt2_idx_prefix} \
  -1 $fq1 \
  -2 $fq2 \
  | samtools view -b - \
  > $outunsort; \
  samtools sort \
  --threads ${params.mapcores} \
  -o $outbam \
  $outunsort
  """

}

process rmdup {
  label 'rmdup'
  publishDir "${params.outdir}/bams/initial/"

  when:
  params.mode =~ /(all)/

  input:
  tuple pop, rep, path(inbam) from RMDUP

  output:
  tuple pop, rep, path(outbam) into SNPS, COUNT_RMDUP
  path(outlog)

  script:
  outbam = "${pop}_${rep}.initial.rmdup.bam"
  outlog = "${outbam}.log"
  """
  picard MarkDuplicates \
  ASSUME_SORT_ORDER=coordinate \
  SORTING_COLLECTION_SIZE_RATIO=.05 \
  MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
  MAX_RECORDS_IN_RAM=2500000 \
  OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 \
  REMOVE_DUPLICATES=true \
  DUPLICATE_SCORING_STRATEGY=RANDOM \
  INPUT=$inbam \
  OUTPUT=$outbam \
  METRICS_FILE=$outlog
  """

}

Channel.fromPath( params.sample_info )
  .combine( Channel.fromList( params.vcf_pops ) )
  .set{ LD_IDS }

process get_vcf_inds {
  label 'r'
  publishDir "${params.outdir}/inds/"

  when:
  params.mode =~ /(all)/

  input:
  tuple path(si), pops from LD_IDS

  output:
  tuple pops, path(name) into VCF_INDS

  script:
  name = "${pops}.individuals.txt"
  """
  Rscript ${params.bin}/get_1kg_pop_individuals.R \
  --sample_info $si \
  --pops $pops \
  --name $name
  """

}

Channel.fromPath( params.vcf_dir + '/' + params.vcf_stem )
    .combine(VCF_INDS)
    .combine(params.maf)
    .combine(params.mac)
    .set{ VCFS }

process subset_vcf {
  label 'snps'
  publishDir "${params.outdir}/vcfs/$pops/maf${maf}mac${mac}/"

  when:
  params.mode =~ /(all)/

  input:
  tuple file(vcf), pops, path(inds), maf, mac from VCFS

  output:
  tuple pops, maf, mac, path(outname) into VCFS_SUB

  script:
  (vcfname, chrom) = ( vcf.getName() =~ /ALL\.(chr[1-9,X,Y]{1,2})\b.+/ )[0]
  outname = "${chrom}.vcf.gz"
  if (maf) {
    allele_filt_string = "--maf $maf"
  } else {
    allele_filt_string = "--mac $mac"
  }
  """
  vcftools \
  --gzvcf $vcf \
  --keep $inds \
  $allele_filt_string \
  --recode \
  --stdout \
  | bgzip -c \
  > $outname
  """

}

process get_snps {
  label 'snps'
  publishDir "${params.outdir}/snps/maf${maf}mac${mac}/"

  when:
  params.mode =~ /(all)/

  input:
  tuples pops, maf, mac, path(vcf) from VCFS_SUB

  output:
  path(snpout) into SNPIDS // ${chr}.snps.txt.gz

  script:
  chrom = vcf.simpleName
  snpout = "${chrom}.snps.txt.gz"
  // for now throw out copy number variants other insertion types (preceded by '<')
  // TODO: check escape characters working properly
  """
  pigz -dc $vcf \
  | grep -v "^#" \
  | awk '{{printf ("%s\t%s\t%s\n", \$2, \$4, \$5)}}' \
  | grep -v -e \< \
  | pigz \
  > $snpout
  """

}

process find_intersecting_snps {
  label 'intersecting'
  publishDir "${params.outdir}/hornet"

  when:
  params.mode =~ /(all)/

  input:
  tuple pop, rep, path(inbam), path(snps) from SNPS.combine( SNPIDS.collect() )

  output:
  tuple pop, rep, path("*.fq1.gz"), path("*.fq2.gz") into REMAP
  tuple pop, rep, path("*.to.remap.bam") into FILT_REMAP, COUNT_REMAP
  tuple pop, rep, path("*.keep.bam") into KEEP_MERGE, COUNT_KEEP

  script:
  """
  mkdir snps; \
  mv *.snps.txt.gz snps; \
  python ${params.hornet_code}/find_intersecting_snps.py \
  -p $inbam \
  snps
  """

}

process remap {
  label 'map'
  publishDir "${params.outdir}/hornet"

  when:
  params.mode =~ /(all)/

  input:
  tuple pop, rep, path(fq1), path(fq2) from REMAP

  output:
  tuple pop, rep, path(outbam) into REMAPPED

  script:
  outunsort = "${pop}_${rep}.remap.unnamesort.bam"
  outbam = "${pop}_${rep}.remap.bam"
  """
  bowtie2 \
  -p ${params.mapcores} \
  ${params.bt2_opts} \
  -x ${params.bt2_idx_prefix} \
  -1 $fq1 \
  -2 $fq2 \
  | samtools view -b - \
  > $outunsort; \
  samtools sort \
  -n \
  --threads ${params.mapcores} \
  -o $outbam \
  $outunsort
  """

}

process filter_remapped {
  label 'filter'
  publishDir "${params.outdir}/hornet"

  when:
  params.mode =~ /(all)/

  input:
  tuple pop, rep, path(orig), path(remapped) from FILT_REMAP.combine( REMAPPED, by: [0,1] )

  output:
  tuple pop, rep, path(outbam) into KEPT_MERGE, COUNT_KEPT

  script:
  outbam = "${pop}_${rep}.kept.bam"
  """
  python ${params.hornet_code}/filter_remapped_reads.py \
  -p \
  $orig \
  $remapped \
  $outbam
  """

}

process hornet_merge_mapq {
  label 'samtools'
  publishDir "${params.outdir}/hornet", pattern: "*.unmapq.bam"
  publishDir "${params.outdir}/bams/hornet", pattern: "*.final.bam"

  when:
  params.mode =~ /(all)/

  input:
  tuple pop, rep, path(keep), path(kept) from KEEP_MERGE.combine( KEPT_MERGE, by: [0,1] )

  output:
  tuple pop, rep, path(outbam) into ANAL, COUNT_FINAL
  tuple pop, rep, path(unfiltbam) into COUNT_UNMAPQ

  script:
  unfiltbam = "${pop}_${rep}.unmapq.bam"
  outunsort = "${pop}_${rep}.uncoordsort.bam"
  outbam = "${pop}_${rep}.final.bam"
  """
  samtools merge \
  --threads ${params.samcores} \
  $unfiltbam \
  $keep \
  $kept; \
  samtools view \
  --threads ${params.samcores} \
  -bq ${params.mapq} \
  $unfiltbam \
  > $outunsort; \
  samtools sort \
  --threads ${params.samcores} \
  -o $outbam \
  $outunsort
  """

}

COUNT_RAW.map{ it -> ['raw', it[1], it[2], it[3] ] }
  .concat( COUNT_CUTADAPT.map{ it -> [ 'trimmed', it[0], it[1], [it[2], it[3]] ] } )
  .concat( COUNT_INITIAL.map{ it -> [ 'initial', it[0], it[1], it[2] ] } )
  .concat( COUNT_RMDUP.map{ it -> [ 'rmdup', it[0], it[1], it[2] ] } )
  .concat( COUNT_KEEP.map{ it -> [ 'no_var', it[0], it[1], it[2] ] } )
  .concat( COUNT_REMAP.map{ it -> [ 'var', it[0], it[1], it[2] ] } )
  .concat( COUNT_KEPT.map{ it -> [ 'kept_var', it[0], it[1], it[2] ] } )
  .concat( COUNT_UNMAPQ.map{ it -> [ 'unmapq', it[0], it[1], it[2] ] } )
  .concat( COUNT_FINAL.map{ it -> [ 'final', it[0], it[1], it[2] ] } )
  .set{ COUNT }

process count_reads {
  label 'samtools'
  publishDir "${params.outdir}/counts"

  cpus { rtype in ['raw', 'trimmed'] ? 1 : params.samcores }
  memory { rtype in ['raw', 'trimmed'] ? '4G' : '48G' }

  when:
  params.mode =~ /(all)/

  input:
  tuple rtype, pop, rep, path(reads) from COUNT

  output:
  path(outname) into CONCAT_COUNTS

  shell:
  outname = "${pop}_${rep}_${rtype}.PE_read_count.txt"
  if ( rtype in ['raw', 'trimmed'] ) {
    '''
    l1=$(zcat !{reads[0]} | wc -l); \
    l2=$(zcat !{reads[1]} | wc -l); \
    r1=$( echo $l1/4 | bc); \
    r2=$( echo $l2/4 | bc); \
    rt=$( echo $r1+$r2 | bc); \
    p=$( echo $rt/2 | bc); \
    echo !{pop} $'\t' !{rep} $'\t' !{rtype} $'\t' $p > !{outname}
    '''
  } else {
    '''
    r=$(samtools flagstat --threads !{params.samcores} !{reads} | head -5 | tail -1 | cut -f 1 -d ' '); \
    p=$( echo $r/2 | bc); \
    echo !{pop} $'\t' !{rep} $'\t' !{rtype} $'\t' $p > !{outname}
    '''
  }

}

process concat_counts {
  publishDir "${params.outdir}/counts"

  time = '10m'

  when:
  params.mode =~ /(all)/

  input:
  path(countlist) from CONCAT_COUNTS.collect()

  output:
  path(outname) into PLOT_COUNTS

  shell:
  counts = countlist.findAll{ it.toString().endsWith('.txt') }.sort()
  colnames = ['pop','rep','map_step','PE_reads']
  outname = "PE_read_counts.txt"
  '''
  echo !{colnames.join(" $'\\t' ")} > !{outname}; \
  cat !{counts.join(' ')} >> !{outname}
  '''

}
