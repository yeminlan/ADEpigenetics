## get multiMTL table1 with H3K27ac/H3K9ac/H3K122ac
for i in `ls ../maskedMTL/MTL.* | grep -v "gain\|loss\|Y\|O\|A"`; do ln -s $i .;done
cat MTL.H3K27ac.bed MTL.H3K9ac.bed MTL.H3K122ac.bed | sort -k1,1 -k2,2n | bedtools merge > multiMTL.bed
R --no-save < log/peak_size.R
rm MTL.*.bed

## get presence/absence in each
for i in `ls ../maskedMTL/MTL.* | grep "Y.bed\|O.bed\|A.bed"`; do ln -s $i .;done
for i in `ls MTL.*.bed|sed "s|.bed||g"`; do bedtools intersect -a multiMTL.bed -b $i.bed -c | cut -f4 > In.$i;done
for i in `ls In.*`; do echo "$i" | cat - $i > $i.tmp; done
paste In.MTL.CTCF.*.tmp In.MTL.Rad21.*.tmp In.MTL.H3K27ac.*.tmp In.MTL.H3K9ac.*.tmp In.MTL.H3K122ac.*.tmp In.MTL.H3K4me1.*.tmp | sed "s|In.MTL.|In.|g" > InMTL.txt
rm MTL.*.bed In.MTL.*

## parse AUC from individual-mark-MTL AUC.txt files
# H3K27ac
echo -e "#chr\tstart\tend" |cat - ../run_H3K27ac/peak_calling.v2/MTL.bed > t1
paste t1 ../run_H3K27ac/peak_calling.v2/AUC.txt > t2
bedtools intersect -a multiMTL.bed -b t2 -wa -wb > t3
R --no-save < log/summarize_AUC.R
mv AUC.txt AUC.H3K27ac.txt
rm t1 t2 t3
# H3K9ac
echo -e "#chr\tstart\tend" |cat - ../run_H3K9ac/peak_calling.v2/MTL.bed > t1
paste t1 ../run_H3K9ac/peak_calling.v2/AUC.txt > t2
bedtools intersect -a multiMTL.bed -b t2 -wa -wb > t3
R --no-save < log/summarize_AUC.R
mv AUC.txt AUC.H3K9ac.txt
rm t1 t2 t3
# H3K122ac
echo -e "#chr\tstart\tend" |cat - ../run_H3K122ac/peak_calling.v2/MTL.bed > t1
paste t1 ../run_H3K122ac/peak_calling.v2/AUC.txt > t2
bedtools intersect -a multiMTL.bed -b t2 -wa -wb > t3
R --no-save < log/summarize_AUC.R
mv AUC.txt AUC.H3K122ac.txt
rm t1 t2 t3
#
paste AUC.H3K27ac.txt AUC.H3K9ac.txt AUC.H3K122ac.txt > AUC.txt
rm AUC.H3K27ac.txt AUC.H3K9ac.txt AUC.H3K122ac.txt

## get MTL.result.csv
R --no-save < log/get_csv.R

## add result of pairwise comparison
R --no-save < log/pairwise_comparison.R

## add nearest gene
annotatePeaks.pl multiMTL.bed hg19 > multiMTL.anno.txt 
R --no-save < log/add_nearest_gene.R

## add rna result
ln -s ../brain.rna/brain.rna.csv .
R --no-save < log/add_rna.R

## has CTCF/Rad21 within 5kb
bedtools closest -a multiMTL.bed -b ../run_CTCF/peak_calling.v2/MTL.bed -t "first" -d | awk '{print($7<=5000)?1:0}' > MTL.hasCTCF_5kb.bed
bedtools closest -a multiMTL.bed -b ../run_Rad21/peak_calling.v2/MTL.bed -t "first" -d | awk '{print($7<=5000)?1:0}' > MTL.hasRad21_5kb.bed

## collapse MTL table so that each transcript appear once (record only the nearest peak to its TSS)
R --no-save < log/MTL.result.collapsed.R

## count multiMTLs for venn diagram
R --no-save < log/get_venn.R

## overlap with PrefrontalCortex K9ac peaks
bedtools intersect -a multiMTL.bed -b ../run_ROSMAP_PrefrontalCortex_H3K9ac/differential/Again.bed -c | awk '{print ($4>0)?"gain":"-"}' > t1
bedtools intersect -a multiMTL.bed -b ../run_ROSMAP_PrefrontalCortex_H3K9ac/differential/Aloss.bed -c | awk '{print ($4>0)?"loss":"-"}' > t2
paste -d- t1 t2 | sed "s|---|nonsig|g" | sed "s|--||g" > PrefrontalCortex.sig.txt
rm t1 t2

## overlap with PrefrontalCortex K9ac Tau-burden peaks
bedtools intersect -a multiMTL.bed -b ../run_ROSMAP_PrefrontalCortex_H3K9ac/ROSMAP/peaks_correlated_with_tau.pos.bed -c | awk '{print ($4>0)?"gain":"-"}' > t1
bedtools intersect -a multiMTL.bed -b ../run_ROSMAP_PrefrontalCortex_H3K9ac/ROSMAP/peaks_correlated_with_tau.neg.bed -c | awk '{print ($4>0)?"loss":"-"}' > t2
paste -d- t1 t2 | sed "s|---|nonsig|g" | sed "s|--||g" > PrefrontalCortex_TauBurden.sig.txt
rm t1 t2

