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

Channel.fromPath( params.bam_glob )
  .map{ it -> [ it.getName() - ~/(_001_hg19\.bwt2pairs\.bam)/, it ] }
  .set{ BAM }

process bam_to_fastqs {

  storeDir "${params.outdir}/$sample/"

  input:
  tuple prefix, path(bam) from BAM

  output:
  tuple prefix, path(fq1), path(fq2) into FASTQS

  script:
  sample = prefix - ~/(_L00[0-9]{1}$)/
  fq1 = "${prefix}_R1_001.fastq.gz"
  fq2 = "${prefix}_R2_001.fastq.gz"
  // -n don't append '/1' or '/2' to the end of read names
  // -i add Illumina Casava 1.8 format entry to header (eg 1:N:0:ATCACG)
  // ^doesn't work without additional option specifying the barcode
  """
  samtools fastq \
  -n \
  -1 $fq1 \
  -2 $fq2 \
  --threads ${task.cpus} \
  $bam
  """

}
