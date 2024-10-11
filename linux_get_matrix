cd $1

mkdir qc
mkdir star

/home/bingxu/software/TrimGalore-0.6.10/trim_galore --paired -q 30 -j 18 -o $1/qc/ --gzip --stringency 4 --fastqc *_1.fq.gz *_2.fq.gz

STAR \
  --genomeDir /home/bingxu/RNAseq/arabidopsis_database_index \
  --runThreadN 20 \
  --readFilesCommand zcat \
  --runMode alignReads \
  --quantMode GeneCounts \
  --limitOutSJcollapsed 4000000 \
  --outSAMtype BAM SortedByCoordinate \
  --limitBAMsortRAM 16000000000 \
  --readFilesIn $1/qc/*_1.fq.gz $1/qc/*_2.fq.gz \
  --outFileNamePrefix $1/star/
featureCounts -a "/home/bingxu/RNAseq/Arabidopsis_data/genomic.gtf" -p -T 8 -o $1/matrix.txt $1/star/Aligned.sortedByCoord.out.bam

# bash ./this_script ./folder_your_raw_data
