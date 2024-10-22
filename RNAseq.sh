#!/bin/bash

### Quality filter (optional)
gunzip *.gz
files=$(ls *.fastq|grep -e ".fastq$")
for file in ${files[*]}
do
	echo $file "processing..."
	/master/xli/software/fastx_toolkit/bin/fastq_quality_filter -Q33 -q 20 -p 80 -i $file -o $file-HQ
	echo $file "done!"
done


### rename file, simplfy (optional)
mv MAL12-10a_S87_L008_R1_001.fastq-HQ MAL12-10a.R1

### only use pair-end reads (optional)
export PATH=$PATH:/master/xli/software/bin
files=(ls ./)
for file in ${files[*]}
do
	echo "processing $file..."
	echo "files: -f $file.R1 $file.R2"
	common-find.pl -f $file.R1 -b $file.R2
	echo "$file done!"
done

### mapping & reads count
export PATH=$PATH:/master/xli/software/STAR-2.5.3a/bin/Linux_x86_64/:/master/xli/software/samtools-1.3/:/master/xli/software/HTSeq-0.7.2/scripts
#Build a genome index
STAR  --runMode genomeGenerate --runThreadN 2 --genomeDir /path/to/genomeDir --genomeFastaFiles /path/to/genome/fasta

for file in ${files[*]}
do
	echo "STAR processing $file..."
	STAR --runThreadN 12 --genomeDir /master/xli/Index/Pfal_STAR_index/ --readFilesIn $file.R1-filter $file.R2-filter --outFileNamePrefix $file --outReadsUnmapped Fasta >1.log 2>&1
	echo "STAR $file done!"
	echo "samtools $file..."
	samtools view -bS ${file}Aligned.out.sam > $file.bam
	samtools sort -n $file.bam -o $file.sorted
	echo "samtools $file done!"
	echo "htseq-count $file..."
	htseq-count -f bam -r name -s no --mode=intersection-nonempty -i Parent $file.sorted /master/xli/Index/Pfal32/PlasmoDB-32_Pfalciparum3D7.gff > $file.count
	echo "htseq-count $file done!"
done

### merge results
paste MAL12-1a.count MAL12-1b.count MAL12-2a.count MAL12-2b.count MAL12-3a.count MAL12-3b.count MAL12-4a.count MAL12-4b.count MAL12-5a.count MAL12-5b.count MAL12-6a.count MAL12-6b.count MAL12-7a.count MAL12-7b.count MAL12-8a.count MAL12-8b.count MAL12-9a.count MAL12-9b.count MAL12-10a.count MAL12-10b.count MAL39-1a.count MAL39-1b.count MAL39-2a.count MAL39-2b.count MAL39-3a.count MAL39-3b.count MAL39-4a.count MAL39-4b.count MAL39-5a.count MAL39-5b.count MAL39-6a.count MAL39-6b.count MAL39-7a.count MAL39-7b.count MAL39-8a.count MAL39-8b.count MAL39-9a.count MAL39-9b.count MAL39-10a.count MAL39-10b.count MAL47-1a.count MAL47-1b.count MAL47-2a.count MAL47-2b.count MAL47-3a.count MAL47-3b.count MAL47-4a.count MAL47-4b.count MAL47-5a.count MAL47-5b.count MAL47-6a.count MAL47-6b.count MAL47-7a.count MAL47-7b.count MAL47-8a.count MAL47-8b.count MAL47-9a.count MAL47-9b.count MAL47-10a.count MAL47-10b.count MAL54-1a.count MAL54-1b.count MAL54-2a.count MAL54-2b.count MAL54-3a.count MAL54-3b.count MAL54-4a.count MAL54-4b.count MAL54-5a.count MAL54-5b.count MAL54-6a.count MAL54-6b.count MAL54-7a.count MAL54-7b.count MAL54-8a.count MAL54-8b.count MAL54-9a.count MAL54-9b.count MAL54-10a.count MAL54-10b.count | cut -f 1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100,102,104,106,108,110,112,114,116,118,120,122,124,126,128,130,132,134,136,138,140,142,144,146,148,150,152,154,156,158,160 >MA2

### RPKM caculation
library(edgeR)
library(limma)
x<- read.delim("C:/Users/xli/Desktop/count/MA2_gene_length.txt", header=TRUE, row.names=1,stringsAsFactors=FALSE)
targets<-data.frame(Lane=c(1:80),Treatment=c(rep("emb",80)),Lable=c(paste("emb",1:80,sep="")))
y<-DGEList(counts=x[,1:80],group=targets$Treatment,genes=data.frame(Length=x[,81]))
colnames(y)<-targets$Label
y$samples$libsize<-colSums(y$counts)
y<-calcNormFactors(y)
write.csv(rpkm(y,normalized.lib.sizes=T,gene.length=x[,81]),"MA2_RPKM.csv")

