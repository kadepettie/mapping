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

params{
  outdir = '190721output'
  mode = "all" // all
  fastq_comb_glob = "/godot/shared_data/atacseq/fastqs/combined/*.fastq.gz"
  fastq_sub_glob = "/godot/shared_data/ashley_ATACseqData/fastqs/combined_runs/*.fastq.gz"
}

process{
  errorStrategy = 'finish'
  maxForks = 200
  stageInMode = 'rellink'
  cpus = 1
  memory = '4G'
  time  = '8h'
}

executor{
  name = 'slurm'
  submitRateLimit = '1 sec'
  queueSize = 500
}

notification{
  enabled = true
  to = 'kpettie@stanford.edu'
}

Channel.fromPath( params.fastq_comb_glob )
  .map{ itg = ( it.getName() =~ /([a-z]{3})([1-2]{1})_(R[1-2]{1})\.fastq\.gz/ )[0]
        return [ "comb", itg[1].toUpperCase(), "rep${itg[2]}", itg[3], it ] }
  .into{ COMB }

Channel.fromPath( params.fastq_sub_glob )
  .map{ itg = ( it.getName() =~ /([a-z]{3})([1-2]{1})_(R[1-2]{1})\.fastq\.gz/ )[0]
        return [ "sub", itg[1].toUpperCase(), "rep${itg[2]}", itg[3], it ] }
  .into{ SUB }

process rename_comb {

  stageInMode "symlink"

  when:
  params.mode =~ /(all)/

  input:
  tuple t, pop, rep, rd, path(fq) from RENAME.concat(SUB)

  output:
  tuple fname, path(fqout) into RENAMED

  script:
  fname = fq.getName
  fqout = "${t}_${fname}"
  """
  mv $fq $fqout
  """

}

RENAMED
  .groupTuple( by: 0 )
  .set{ COMBSUB }

process subtract {

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
