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

// Nextflow pipeline to convert bam file of PE reads back to separate fastqs.

out = file(params.outdir)
if( !out.exists() ) out.mkdir()

Channel.fromPath( params.bwt2_global_glob )
  .map{ it -> [ it.getName() - ~/(_hg19\.bwt2glob\.)(bam)?(unmap\.fastq)?(\.gz)?/, it ] }
  .groupTuple(by: [0])
  .set{ READS }

process bam_to_fastqs {
  tag "$fqout = $bam + $fq"
  storeDir "${params.outdir}/$sample/"

  input:
  tuple prefix, path(reads) from READS

  output:
  tuple prefix, path(fq1), path(fq2) into FASTQS

  shell:
  sample = prefix - ~/(_L00[0-9]{1})?(_R1|_R2)?(_001)?$/
  bamfq = reads.collect{ it.getName() }.sort()
  bam = bamfq[0]
  fq = bamfq[1] // must be uncompressed for now (default in hicpro)
  fqout = "${prefix}.fastq.gz"
  // -n don't append '/1' or '/2' to the end of read names
  // -i add Illumina Casava 1.8 format entry to header (eg 1:N:0:ATCACG)
  // ^doesn't work without additional option specifying the barcode
  '''
  bctag=$(head -1 !{fq} | awk -F $' ' '{ print $2 }');
  samtools fastq \
  -n \
  --threads !{task.cpus} \
  !{bam} \
  | paste - - - - \
  | awk -F $'\t' -v bctag="$bctag" '{ print $1" "bctag, $2, $3, $4 }' OFS='\t' \
  | cat - <(cat !{fq} | paste - - - - ) \
  | sort \
  -k 1,1 \
  -S !{task.memory} \
  -T !{params.tmpdir}/!{prefix} \
  --parallel=!{task.cpus} \
  | tr "\t" "\n" \
  | pigz -c \
  -p !{task.cpus} \
  > !{fqout}
  '''

}
