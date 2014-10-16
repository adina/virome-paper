#example of generating bowtie scripts
for x in *R1*gz; do echo "bowtie2 --sensitive -x ../VLP-final-contigs.renamed.fa -1 ${x%*_R*}_R1_*gz -2 ${x%*_R*}_R2*gz -S ${x%*_R*}.VLP-VLP.sam 1> ${x%*_R*}.VLP-VLP.stdout 2> ${x%*_R*}.VLP-VLP.errout"; done > bowtie1.sh
