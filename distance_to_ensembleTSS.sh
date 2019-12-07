awk '{FS=OFS="\t"}{if($6=="+")print $1,$2,$2;else print $1,$3,$3}' /project/ibilab/library/annotation/hg19/hg19_ensembl.bed12 | sort -k1,1 -k2,2n | uniq  > hg19_TSS.bed
bedtools closest -a multiMTL.bed -b hg19_TSS.bed -d -t first | cut -f7 | sed "s|-1|1000000|g" > t1
R --no-save < log/distance_to_ensembleTSS.R
rm t1 hg19_TSS.bed





