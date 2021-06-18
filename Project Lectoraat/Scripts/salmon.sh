conda activate salmon
#download transcriptoom
`curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/339/765/GCF_003339765.1_Mmul_10/GCF_003339765.1_Mmul_10_rna.fna.gz`
#folder for data
mkdir data
cd data

#build an index on the transcriptome
salmon index -t athal.fa.gz -i athal_index

salmon index -t GCF_003339765.1_Mmul_10_rna.fna.gz -i GCF_003339765.1_Mmul_10_index

#Obtaining sequencing data
for i in `seq 16 31`; 
do 
  mkdir SRR120105${i}; 
  cd SRR120105${i}; 
  wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR120/0${i}/SRR120105${i}/SRR120105${i}.fastq.gz
  cd ..; 
done
cd ..

 
#Quantify our samples
for fn in data/SRR120105{16..31};
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
salmon quant -i GCF_003339765.1_Mmul_10_index -l A \
         -1 ${fn}/${samp}.fastq
         -p 8 --validateMappings -o quants/${samp}_quant
done 
