#!/bin/bash

if [[ $# == 5 ]]; then

   which trimmomatic
   if [[ $? != 0 ]]; then
      echo "Could not locate trimmomatic"
      echo "Please install by one of the following means"
      echo "'sudo apt-get install trimmomatic'"
      echo "'conda install trimmomatic -c bioconda'"
      echo "Add bioconda channel if it doesn't already exist: 'conda config --add chennels bioconda'"
      exit $?
   else
      id=$1
      dname=$2 # strip trailing forward slash
      leadx=$3
      trailx=$4
      t=$5
      mkdir -p ${dname}paired
      mkdir -p ${dname}unpaired
      pdr="${dname}paired/"
      udr="${dname}unpaired/"
      n="$((50/$t))"
      while read -r line; do
          trimmomatic PE -phred33 $line ILLUMINACLIP:$HOME/bioTools/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:0 TRAILING:38 SLIDINGWINDOW:4:15 MINLEN:36 -threads $t
      done < $id
      mv ${dname}*_fp.fastq.gz ${dname}*_rp.fastq.gz $pdr
      mv ${dname}*_fu.fastq.gz ${dname}*_ru.fastq.gz $udr
   fi
else
        echo """
	Usage: ./trimFastq.sh <idlist> <fpath> <leading> <trailing> <threads>

	 idlist: File containing SRA run accessions. (NB: Must be in same path as fastq files)
	  fpath: Path to fastq files. (NB: Must end with a forward '/' slash)
	threads: Optional [default = 1]
	leading: Number of bases ot crop from the 5' end
       trailing: Number of bases to crop from the 3' end

	e.g.: 1_S1_L001_R1_001.fastq.gz 1_S1_L001_R2_001.fastq.gz
	
    """
fi
