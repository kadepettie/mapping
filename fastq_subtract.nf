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

// Nextflow pipeline to subtract reads from one fastq file from another
// combined fastq file.

out = file(params.outdir)
if( !out.exists() ) out.mkdir()

Channel.fromPath( params.fastq_comb_glob )
  .map{ itg = ( it.getName() =~ /([a-z]{3})([1-2]{1})_(R[1-2]{1})\.fastq\.gz/ )[0]
        return [ "comb", itg[1].toUpperCase(), "rep${itg[2]}", itg[3], it ] }
  .set{ COMB }

Channel.fromPath( params.fastq_sub_glob )
  .map{ itg = ( it.getName() =~ /([a-z]{3})([1-2]{1})_(R[1-2]{1})\.fastq\.gz/ )[0]
        return [ "sub", itg[1].toUpperCase(), "rep${itg[2]}", itg[3], it ] }
  .set{ SUB }

process rename_comb {

  stageInMode "symlink"

  when:
  params.mode =~ /(all)/

  input:
  tuple t, pop, rep, rd, path(fq) from COMB.concat(SUB)

  output:
  tuple fname, path(fqout) into RENAMED

  script:
  fname = fq.getName()
  fqout = "${t}_${fname}"
  """
  mv $fq $fqout
  """

}

RENAMED
  .groupTuple( by: 0 )
  .set{ COMBSUB }

process subtract {

  publishDir "${params.outdir}/run2/"

  when:
  params.mode =~ /(all)/

  input:
  tuple fname, path(fq) from COMBSUB

  output:
  path(fname) into SUBBED

  script:
  """
  gunzip -c \
  ${fq[0]} \
  ${fq[1]} \
  | paste  - - - - \
  | sort \
  | uniq -u \
  | tr "\\t" "\\n" \
  > $fname
  """

}
