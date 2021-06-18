#https://combine-lab.github.io/salmon/getting_started/
#conda create --name salmon

#terminal opnieuw opstarten
#conda install -c bioconda salmon
#conda install tbb=2020.2
#conda install -c conda-forge -c bioconda salmon=1.5.0

#activate environment
conda activate salmon

#folder for tutorial
mkdir salmon_tutorial
cd salmon_tutorial
 
#download transcriptoom
curl https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/339/765/GCF_003339765.1_Mmul_10/GCF_003339765.1_Mmul_10_rna.fna.gz -o athal.fa.gz

#folder for data
mkdir data
cd data

#build an index on the transcriptome
salmon index -t athal.fa.gz -i athal_index

#obtain the raw data and  the corresponding read files 
for i in `seq 25 40`; 
do 
  mkdir DRR0161${i}; 
  cd DRR0161${i}; 
  wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR016/DRR0161${i}/DRR0161${i}_1.fastq.gz; 
  wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR016/DRR0161${i}/DRR0161${i}_2.fastq.gz; 
  cd ..; 
done
cd .. 

#Quantify our samples
for fn in data/DRR0161{25..40};
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
salmon quant -i athal_index -l A \
         -1 ${fn}/${samp}_1.fastq.gz \
         -2 ${fn}/${samp}_2.fastq.gz \
         -p 8 --validateMappings -o quants/${samp}_quant
done 
