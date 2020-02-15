#!/bin/bash
#requires mmv
#example: /home/etienne/scripts/allDMRs/DMR_caller_parallel.sh -g /media/etienne/save/lab/jbrowse/genomes/gala_v4/myGalaAnnotation/Gala4_EB.fasta -c /media/etienne/save/lab/methylomes/gala/WGBS/2.conf -o /media/etienne/now/parallel/sample2_gff -p 4 -f /media/etienne/save/lab/methylomes/gala/WGBS/Gala4_DMR_gene_regions.gff

#parse command line arguments
gff_file="none"
while [[ "$#" -gt 0 ]]; do case $1 in
	-g|--genome) ref="$2"; shift;;
	-c|--conf) conf_tvg="$2"; shift;;
	-o|--output)  output_directory="$2"; shift;;
	-p|--processors) processors="$2"; shift;;
	-f|--gff) gff_file="$2"; shift;;
	 *) echo "Unknown parameter passed: $1"; exit 1;;
esac; shift; done


START=$(date +%s.%N)
echo "preparing pipeline..."
script_location="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

if [ $gff_file == "none" ]; then
	pipeline=$script_location"/dmr_pipeline.py"
else 
	pipeline=$script_location"/dmr_pipeline_gff.py"
fi

echo "pipeline is here: "$pipeline

mkdir -p $output_directory
mkdir -p $output_directory/genome
mkdir -p $output_directory/controls
mkdir -p $output_directory/samples
mkdir -p $output_directory/dmrs
echo "splitting fasta file..."
if [ -f $output_directory/genome/chr1.fasta ]
then
	echo "genome already processed"
else
	seqkit split --quiet -f -i $ref -O $output_directory/genome -o ""
    mmv "$output_directory/genome/*.id_*.fasta" "$output_directory/genome/#2.fasta" #rename files to nice chr nomenclature
fi
#parse configuration
var1=$(sed '1q;d' $conf_tvg)
var2=$(sed '2q;d' $conf_tvg)

controls_tmp=$(echo $var1 | cut -f2 -d=)
samples_tmp=$(echo $var2 | cut -f2 -d=)

IFS=',' # comma is set as delimiter
read -ra controls <<< "$controls_tmp" # str is read into an array as tokens separated by IFS
if [ "$(ls -A $output_directory"/controls/")" ]; then
    echo "controls methylation file splitting already done"
	else
		for i in "${controls[@]}"; do # access each element of array
    		data=${i##*/}
    		outpath=$output_directory"/controls/"$data"_"
    		echo "working on" $data
    		awk -v var="$outpath" 'FNR >1 {if (last != $1) close(last); print >> var$1; last = $1}' $i
		done
fi

read -ra samples <<< "$samples_tmp" # str is read into an array as tokens separated by IFS
if [ "$(ls -A $output_directory"/samples/")" ]; then
    echo "samples methylation file splitting already done"
	else
		for i in "${samples[@]}"; do # access each element of array
    		data=${i##*/}
    		outpath=$output_directory"/samples/"$data"_"
    		echo "working on" $data
    		awk -v var="$outpath" 'FNR >1 {if (last != $1) close(last); print >> var$1; last = $1}' $i
		done
fi

#call the DMRs with Nicola's script chromosome by chromosome, first prepare scripts to run with parallel

if [ -f $output_directory/dmrs/commands.txt ]
	then 
		rm $output_directory/dmrs/commands.txt
fi

for chromosome in $output_directory/genome/*.fasta
do
	#((j=j%N)); ((j++==0)) && wait
	info=$(basename $chromosome .fasta)
	#echo "processing chromosome" $info
	#echo "output:" $output_directory
	mkdir -p $output_directory/dmrs/$info
	mkdir -p $output_directory/dmrs/$info/data
	mkdir -p $output_directory/dmrs/$info/results
	controlsMeth=($output_directory/controls/*_$info)
	samplesMeth=($output_directory/samples/*_$info)
	varC=$(echo "${controlsMeth[@]}")
	varS=$(echo "${samplesMeth[@]}")
	echo "echantillon1="$varC > $output_directory/dmrs/$info/$info.conf
	echo "echantillon2="$varS >> $output_directory/dmrs/$info/$info.conf
	sed -i 's/[\t ]/,/g' $output_directory/dmrs/$info/$info.conf
	if [ -f $output_directory/dmrs/$info/results/dmrs_not_merged.sign.gff3 ]
	then 
		echo "dmr calling already done for" $info
	else
		echo "python3 $pipeline $chromosome 200 50 200 $output_directory/dmrs/$info/$info.conf $output_directory/dmrs/$info $info $gff_file" >> $output_directory/dmrs/commands.txt
	fi
done

echo "calculating DMRs, this will take a while..."

parallel -j $processors :::: $output_directory/dmrs/commands.txt

#check if all dmrs have been called
incomplete=0
for chromosome in $output_directory/genome/*.fasta
do
	info=$(basename $chromosome .fasta)
	if [ ! -f $output_directory/dmrs/$info/results/dmrs_not_merged.sign.gff3 ]; then
    	echo "DMR calling failed for" $info "try re-running this script with less processors (use -p 1 option).."
    	incomplete=1
	fi
done

	if [ $incomplete == 0 ]
	then
		echo "all dmrs called, generating gff3..."
		conf=$(basename $conf_tvg .conf)
		find $output_directory/dmrs -name "dmrs_not_merged.sign.gff3" |xargs cat > $output_directory/$conf.tab #fetch all results into one file
		#filter DMRs and get nice formatting
		Rscript $script_location/filterDMR.R $output_directory/$conf.tab $output_directory
		echo "cleaning up..."
		rm -r $output_directory/controls
		rm -r $output_directory/samples
		rm -r $output_directory/genome
		#rm -r $output_directory/dmrs
	fi

END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)
echo "Analysis time:" $DIFF " seconds"
#read -p "Press [Enter] key continue..."