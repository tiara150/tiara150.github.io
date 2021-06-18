#!/bin/bash

#download transcriptoom
`curl -O $1`
#folder for data

 transcriptoom=$(echo $1 | sed 's/.*\///')
echo $transcriptoom

 index=$(echo $transcriptoom | sed 's/rna.fna.gz/index/g')

#build an index on the transcriptome
salmon index -t $transcriptoom -i $index

mkdir -p data
cd data

#first run accession
firstrac=$2
first=${firstrac:0:6}
lastfirst=${firstrac:9}

#last run accession
lastrac=$3	
lastlast=${lastrac:9}


for i in `seq $lastfirst $lastlast`;
do
  inb=$(printf %03d ${i})
  mkdir -p ${firstrac:0:9}${i};
  cd ${firstrac:0:9}${i};
  if [ "$4" -eq "1" ]
    then

      wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$first/$inb/${firstrac:0:9}${i}/${firstrac:0:9}${i}.fastq.gz
      fi
    if [ "$4" -eq "2" ]

    then
      wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$first/$inb/${firstrac:0:9}${i}/${firstrac:0:9}${i}_1.fastq.gz
      wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$first/$inb/${firstrac:0:9}${i}/${firstrac:0:9}${i}_2.fastq.gz
  fi
  cd ..;
done
cd ..

#Quantify our samples
for ((fn=${lastfirst};fn<=${lastlast};fn++));

do
samp=`basename ${firstrac:0:9}${fn}`
echo "Processing sample ${samp}"
  if [ "$4" -eq "1" ]
    then
    salmon quant -i $index -l A -r data/${samp}/${samp}.fastq.gz -o ./quants/${samp}
    fi
  if [ "$4" -eq "2" ] 
    then
    salmon quant -i $index -l A \
         -1 data/${samp}/${samp}_1.fastq.gz \
         -2 data/${samp}/${samp}_2.fastq.gz \
         -p 8 --validateMappings -o quants/${samp}
    fi
done
