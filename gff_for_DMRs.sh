#!/bin/bash
#convert a GFF file into an annotation file for DMR calling
#example /home/etienne/scripts/allDMRs/gff_for_DMRs.sh -i /media/etienne/save/lab/jbrowse/genomes/gala_v4/myGalaAnnotation/Gala4_annotated_EB_chr0.gff -o Gala4_DMR_gene_regions.gff -f gene -d 2000
while [[ "$#" -gt 0 ]]; do case $1 in
	-i|--input) input="$2"; shift;;
	-o|--output) output="$2"; shift;;
	-f|--feature) feature="$2"; shift;;
	-d|--distance) distance="$2"; shift;;
	 *) echo "Unknown parameter passed: $1"; exit 1;;
esac; shift; done

awk -F"\t" -v feature=$feature '$3==feature && $7=="+"' $input > $output"_plus"
awk -F"\t" -v feature=$feature '$3==feature && $7=="-"' $input > $output"_minus"

awk -F"\t" -v s=$distance '{$5=$4;$4-=s;$3="5p"}1' OFS='\t' $output"_plus" > $output"_plus_promoter"
awk -F"\t" -v s=$distance '{$4=$5;$5+=s;$3="3p"}1' OFS='\t' $output"_plus" > $output"_plus_terminator"
awk -F"\t" -v s=$distance '{$4=$5;$5+=s;$3="5p"}1' OFS='\t' $output"_minus" > $output"_minus_promoter"
awk -F"\t" -v s=$distance '{$5=$4;$4-=s;$3="3p"}1' OFS='\t' $output"_minus" > $output"_minus_terminator"
awk -F"\t" -v s=$distance '{$3="body"}1' OFS='\t' $output"_plus" > $output"_plus_body"
awk -F"\t" -v s=$distance '{$3="body"}1' OFS='\t' $output"_minus" > $output"_minus_body"
cat $output"_plus_body" $output"_minus_body" $output"_plus_promoter" $output"_plus_terminator" $output"_minus_promoter" $output"_minus_terminator" > $output"_cat"

awk '{split($9,a,"Name="); print a[2]}' $output"_cat" | awk '{split($1,b,";"); print b[1]}' > $output"_names.txt"

paste $output"_names.txt" $output"_cat" > $output"_cat_names"

awk '{print $2,$2"_"$1"_"$4,$5,$6}' OFS='\t' $output"_cat_names" > $output

rm $output"_plus"
rm $output"_minus"
rm $output"_minus_promoter"
rm $output"_minus_terminator"
rm $output"_plus_promoter"
rm $output"_plus_terminator"
rm $output"_plus_body"
rm $output"_minus_body"
rm $output"_cat"
rm $output"_names.txt"
rm $output"_cat_names"

