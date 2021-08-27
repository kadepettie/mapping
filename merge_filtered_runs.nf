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

// Nextflow pipeline to merge mapped, filtered reads from separate sequencing
// runs of the same library and remove any duplicate reads between the two runs.

out = file(params.outdir)
if( !out.exists() ) out.mkdir()

Channel.fromPath(params.bam_runs)
  .map{ itg = ( it.toString() =~ /.*(run[1-9]).*\/(?!.*\/)(.*)\.final\.bam$/ )[0]
        return [ itg[1], itg[2], it ]
   }
  .set{ RENAME }

// rename to avoid filename conflicts in merge staging directory

process rename {
  stageInMode 'symlink'

  when:
  params.mode =~ /(all)/

  input:
  tuple run, samp, bam from RENAME

  output:
  tuple samp, outbam into TO_MERGE

  script:
  outbam = "${run}_${samp}.bam"
  """
  cp $bam $outbam
  """

}

TO_MERGE
  .groupTuple( by: [0] )
  .tap{ MERGE_CHECK }
  .set{ MERGE }

MERGE_CHECK.view()

process merge_runs {
  label 'samtools'
  publishDir "${params.outdir}/merged"
  stageInMode "symlink"

  when:
  params.mode =~ /(all)/

  input:
  tuple samp, path(bams) from MERGE

  output:
  tuple samp, path(outbam) into RMDUP, COUNT_MERGED

  script:
  bamstring = bams.findAll{ it.toString().endsWith('.bam') }.sort()
  outunsort = "${samp}.uncoordsort.bam"
  outbam = "${samp}.dups.bam"
  """
  samtools merge \
  --threads ${params.samcores} \
  $outunsort \
  ${bamstring.join(' ')}; \
  samtools sort \
  --threads ${params.samcores} \
  -o $outbam \
  $outunsort
  """

}

// rmdup then count

process rmdup {
  label 'rmdup'
  publishDir "${params.outdir}/rmdup/", pattern: "*.bam"
  publishDir "${params.outdir}/rmdup/logs/", pattern: "*.log"

  when:
  params.mode =~ /(all)/

  input:
  tuple samp, path(inbam) from RMDUP

  output:
  tuple samp, path(outbam) into ANAL, COUNT_RMDUP
  path(outlog)

  script:
  outbam = "${samp}.bam"
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

COUNT_MERGED.map{ it -> ['merged', it[0], it[1]] }
  .concat( COUNT_RMDUP.map{ it -> ['rmdup', it[0], it[1]] } )
  .set{ COUNT }

process count_reads {
  label 'samtools'
  publishDir "${params.outdir}/counts/separate"

  when:
  params.mode =~ /(all)/

  input:
  tuple rtype, samp, path(reads) from COUNT

  output:
  path(outname) into CONCAT_COUNTS

  shell:
  outname = "${samp}_${rtype}.PE_read_count.txt"
  '''
  r=$(samtools flagstat --threads !{params.samcores} !{reads} | head -5 | tail -1 | cut -f 1 -d ' '); \
  p=$( echo $r/2 | bc); \
  echo !{samp} $'\t' !{rtype} $'\t' $p > !{outname}
  '''

}

process concat_counts {
  publishDir "${params.outdir}/counts/"

  time = '10m'

  when:
  params.mode =~ /(all)/

  input:
  path(countlist) from CONCAT_COUNTS.collect()

  output:
  path(outname) into PLOT_COUNTS

  shell:
  counts = countlist.findAll{ it.toString().endsWith('.txt') }.sort()
  colnames = ['samp','map_step','PE_reads']
  outname = "PE_read_counts.txt"
  '''
  echo !{colnames.join(" \\$'\\t' ")} > !{outname}; \
  cat !{counts.join(' ')} >> !{outname}
  '''

}
