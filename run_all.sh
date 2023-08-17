# my version

num_threads=10
topdir="/public/groups/sanfordlab/people/juphilip/results/primates/andrew_pipeline/remap_01_2020/"

## make species outdirs

human_out="$topdir""/human_out/"
chimpanzee_out="$topdir""/chimpanzee_out/"
orangutan_out="$topdir""/orangutan_out/"


mkdir -p $human_out
mkdir -p $chimpanzee_out
mkdir -p $orangutan_out



## paths to genome fasta, fai files

genome_base="/public/groups/sanfordlab/people/juphilip/datasets/references/"

human_genome_fasta="$genome_base""/human/hg38/hg38.fa"
human_genome_fai="$genome_base""/human/hg38/hg38.fa.fai"

chimpanzee_genome_fasta="$genome_base""/chimp/panTro6/panTro6.fa"
chimpanzee_genome_fai="$genome_base""/chimp/panTro6/panTro6.fa.fai"

orangutan_genome_fasta="$genome_base""/orang/ponAbe3/ponAbe3.fa"
orangutan_genome_fai="$genome_base""/orang/ponAbe3/ponAbe3.fa.fai"


## fastq prefixes

human_fastq_dir=/public/groups/sanfordlab/people/juphilip/datasets/primate_fastq/human/
chimpanzee_fastq_dir=/public/groups/sanfordlab/people/juphilip/datasets/primate_fastq/chimp/
orangutan_fastq_dir=/public/groups/sanfordlab/people/juphilip/datasets/primate_fastq/orang/

human_fastq_prefix=$(ls /public/groups/sanfordlab/people/juphilip/datasets/primate_fastq/human/ | awk -F"_[12].fastq" '{ print $1}' | sort | uniq)
chimpanzee_fastq_prefix=$(ls /public/groups/sanfordlab/people/juphilip/datasets/primate_fastq/chimp/ | awk -F"_[12].fastq" '{ print $1}' | sort | uniq)
orangutan_fastq_prefix=$(ls /public/groups/sanfordlab/people/juphilip/datasets/primate_fastq/orang/ | awk -F"_[12].fastq" '{ print $1}' | sort | uniq)


## paths to bam file directories
human_bam="$topdir""/human_bam/"
chimpanzee_bam="$topdir""/chimp_bam/"
orangutan_bam="$topdir""/orang_bam/"



## make sample_info.tsv

#cat <( echo -e "sample_name\tfraction" ) <( paste <(ls $human_bam/*.bam | awk -F\/ '{ print $NF}' | awk -F"_" '{ print $1"_"$2}') <(echo -e "2\n2\n3\n3\n4\n4\n5\n5\n6\n6\n7\n7\n8\n8\n1\n1\n1")) > $human_out/sample_info.tsv &

#cat <( echo -e "sample_name\tfraction" ) <( paste <(ls $chimpanzee_bam/*.bam | awk -F\/ '{ print $NF}' | awk -F"_" '{ print $1"_"$2"_"$3}') <(echo -e "1\n1\n2\n2\n3\n3\n4\n4\n5\n5\n6\n6")) > $chimpanzee_out/sample_info.tsv &

#cat <( echo -e "sample_name\tfraction" ) <( paste <(ls $orangutan_bam/*.bam | awk -F\/ '{ print $NF}' | awk -F"_" '{ print $1"_"$2}') <(echo -e "1\n2\n3\n4\n5\n6\n1\n2\n3\n4\n5\n6")) > $orangutan_out/sample_info.tsv &



## paths to gencode v27 and derived CAT GTFs
annotation_base="$genome_base"

awk '$3=="gene" {OFS="\t"; print $14, $12}' "$annotation_base""/human/hg38/gencode.v32.annotation.gtf" | sort | uniq | tr -d '\";' > "$annotation_base""/human/hg38/gencode.v32_gene_names_types.tsv"
sed -i '1i gene\ttype' "$annotation_base""/human/hg38/gencode.v32_gene_names_types.tsv"

human_full_gtf="$annotation_base""/human/hg38/hg38_cat_hg38_V32_Lv3_v2.gtf"
chimpanzee_full_gtf="$annotation_base""/chimp/panTro6/panTro6_cat_hg38_V32_Lv3_v2.gtf"
orangutan_full_gtf="$annotation_base""/orang/ponAbe3/ponAbe3_cat_hg38_V32_Lv3_v2.gtf"


## ??
human_tx_gtf="$annotation_base""/gencode.v32.basic.annotation.gtf"
chimpanzee_tx_gtf="$annotation_base""/chimp/panTro6/panTro6_cat_gencode.v32.basic.gtf"
orangutan_tx_gtf="$annotation_base""/orang/ponAbe3/ponAbe3_cat_gencode.v32.basic.gtf"




## rename transcripts in GTFs

renamed_gtf_dir="$annotation_base""/renamed_transcripts/"

#mkdir -p "$renamed_gtf_dir"

rename_tx_script="/public/home/anjowall/repos/multi_species_as_events/rename_cat_transcripts.py"

#/public/home/anjowall/anaconda2/bin/python2.7 "$rename_tx_script" \
#	--reference_gtf "$human_tx_gtf" \
#	--cat_gtf "$chimpanzee_full_gtf" \
#	--outdir "$renamed_gtf_dir" \
#	--prefix Chimp_gencode.v32.basic

#/public/home/anjowall/anaconda2/bin/python2.7 "$rename_tx_script" \
#	--reference_gtf "$human_tx_gtf" \
#	--cat_gtf "$orangutan_full_gtf" \
#	--outdir "$renamed_gtf_dir" \
#	--prefix Orangutan_gencode.v32.basic



chimpanzee_tx_gtf="$renamed_gtf_dir""Chimp_gencode.v32.basic_renamed_to_ref.gtf"
orangutan_tx_gtf="$renamed_gtf_dir""Orangutan_gencode.v32.basic_renamed_to_ref.gtf"


## common tx ids

#join <(join <(awk '$3=="transcript" { print $14}' $chimpanzee_tx_gtf | sort | uniq ) <(awk '$3=="transcript" { print $14}' $orangutan_tx_gtf | sort | uniq )) <(awk '$3=="transcript" { print $12}' $human_tx_gtf | sort | uniq ) | sort | uniq | tr -d '\";' > $renamed_gtf_dir/common_tx_ids.txt



human_common_tx_gtf=$renamed_gtf_dir/Human_gencode.v32.basic_common_tx.annotation.gtf
chimpanzee_common_tx_gtf=$renamed_gtf_dir/Chimp_gencode.v32.basic_common_tx.annotation.gtf
orangutan_common_tx_gtf=$renamed_gtf_dir/Orangutan_gencode.v32.basic_common_tx.annotation.gtf

#fgrep -f $renamed_gtf_dir/common_tx_ids.txt $human_full_gtf | grep -v "PAR_Y" > $human_common_tx_gtf
#fgrep -f $renamed_gtf_dir/common_tx_ids.txt $chimpanzee_full_gtf > $chimpanzee_common_tx_gtf
#fgrep -f $renamed_gtf_dir/common_tx_ids.txt $orangutan_full_gtf > $orangutan_common_tx_gtf


## run stringtie

#stringtie_fun () {

#	stringtie "$1" -G "$2" -o "$3"/$(echo "$1" | awk -F\/ '{ print $NF}').gtf

#}

human_stringtie_out="$human_out""/stringtie/"
chimpanzee_stringtie_out="$chimpanzee_out""/stringtie/"
orangutan_stringtie_out="$orangutan_out""/stringtie/"

#unset i

#(
#for bam in $human_bam/*.bam
#do
#((i=i%num_threads)); ((i++==0)) && wait
#stringtie_fun "$bam" "$human_full_gtf" "$human_stringtie_out" &
#done
#wait
#)


#stringtie --merge -G $human_full_gtf -o $human_stringtie_out/human_merged.gtf  $(ls $human_stringtie_out/*.gtf | grep -v "merged") 

#unset i

#(
#for bam in $chimpanzee_bam/*.bam
#do
#((i=i%num_threads)); ((i++==0)) && wait
#stringtie_fun "$bam" "$chimpanzee_full_gtf" "$chimpanzee_stringtie_out" &
#done
#wait
#)

#stringtie --merge -G $chimpanzee_full_gtf -o $chimpanzee_stringtie_out/chimpanzee_merged.gtf  $(ls $chimpanzee_stringtie_out/*.gtf | grep -v "merged") 

#unset i

#(
#for bam in $orangutan_bam/*.bam
#do
#((i=i%num_threads)); ((i++==0)) && wait
#stringtie_fun "$bam" "$orangutan_full_gtf" "$orangutan_stringtie_out" &
#done
#wait
#)

#stringtie --merge -G $orangutan_full_gtf -o $orangutan_stringtie_out/orangutan_merged.gtf  $(ls $orangutan_stringtie_out/*.gtf | grep -v "merged") 



## run gffread to get fasta files for transcriptomes

human_tx_merged=$human_stringtie_out/human_merged.gtf
chimpanzee_tx_merged=$chimpanzee_stringtie_out/chimpanzee_merged.gtf
orangutan_tx_merged=$orangutan_stringtie_out/orangutan_merged.gtf


human_tx_merged_fasta=$human_stringtie_out/human_merged.fa
chimpanzee_tx_merged_fasta=$chimpanzee_stringtie_out/chimpanzee_merged.fa
orangutan_tx_merged_fasta=$orangutan_stringtie_out/orangutan_merged.fa


#gffread $human_tx_merged -g $human_genome_fasta -w $human_tx_merged_fasta &
#gffread $chimpanzee_tx_merged -g $chimpanzee_genome_fasta -w $chimpanzee_tx_merged_fasta &
#gffread $orangutan_tx_merged -g $orangutan_genome_fasta -w $orangutan_tx_merged_fasta &



## run kallisto
kallisto="/public/groups/sanfordlab/people/anjowall/bin/kallisto_linux-v0.45.0/kallisto"

human_kallisto_index=$human_stringtie_out/human_merged.idx
chimpanzee_kallisto_index=$chimpanzee_stringtie_out/chimpanzee_merged.idx
orangutan_kallisto_index=$orangutan_stringtie_out/orangutan_merged.idx


human_kallisto_out=$human_out/kallisto/
chimpanzee_kallisto_out=$chimpanzee_out/kallisto/
orangutan_kallisto_out=$orangutan_out/kallisto/

mkdir -p $human_kallisto_out
mkdir -p $chimpanzee_kallisto_out
mkdir -p $orangutan_kallisto_out

#$kallisto index -i $human_kallisto_index $human_tx_merged_fasta & 
#$kallisto index -i $chimpanzee_kallisto_index $chimpanzee_tx_merged_fasta &
#$kallisto index -i $orangutan_kallisto_index $orangutan_tx_merged_fasta &


#for prefix in $human_fastq_prefix
#do
#$kallisto quant -i $human_kallisto_index -o $human_kallisto_out/$prefix/ -t $num_threads -b 50 --bias --rf-stranded "$human_fastq_dir""$prefix""_1.fastq.gz" "$human_fastq_dir""$prefix""_2.fastq.gz"
#done

for prefix in $chimpanzee_fastq_prefix
do
$kallisto quant -i $chimpanzee_kallisto_index -o $chimpanzee_kallisto_out/$prefix/ -t $num_threads -b 50 --bias --rf-stranded "$chimpanzee_fastq_dir""$prefix""_1.fastq.gz" "$chimpanzee_fastq_dir""$prefix""_2.fastq.gz"
done

for prefix in $orangutan_fastq_prefix
do
$kallisto quant -i $orangutan_kallisto_index -o $orangutan_kallisto_out/$prefix/ -t $num_threads -b 50 --bias --rf-stranded "$orangutan_fastq_dir""$prefix""_1.fastq.gz" "$orangutan_fastq_dir""$prefix""_2.fastq.gz"
done





## run kallisto with common transcripts only

### common tx fasta

human_common_tx_fasta=$renamed_gtf_dir/Human_gencode.v27.basic_common_tx.annotation.fa
chimpanzee_common_tx_fasta=$renamed_gtf_dir/Chimp_gencode.v27.basic_common_tx.annotation.fa
orangutan_common_tx_fasta=$renamed_gtf_dir/Orangutan_gencode.v27.basic_common_tx.annotation.fa

gffread $human_common_tx_gtf -g $human_genome_fasta -w $human_common_tx_fasta &
gffread $chimpanzee_common_tx_gtf -g $chimpanzee_genome_fasta -w $chimpanzee_common_tx_fasta &
gffread $orangutan_common_tx_gtf -g $orangutan_genome_fasta -w $orangutan_common_tx_fasta &


human_common_tx_kallisto_index=$renamed_gtf_dir/Human_gencode.v27.basic_common_tx.annotation.idx
chimpanzee_common_tx_kallisto_index=$renamed_gtf_dir/Chimp_gencode.v27.basic_common_tx.annotation.idx
orangutan_common_tx_kallisto_index=$renamed_gtf_dir/Orangutan_gencode.v27.basic_common_tx.annotation.idx

human_common_tx_kallisto_out=$human_out/kallisto_common_tx/
chimpanzee_common_tx_kallisto_out=$chimpanzee_out/kallisto_common_tx/
orangutan_common_tx_kallisto_out=$orangutan_out/kallisto_common_tx/

mkdir -p $human_common_tx_kallisto_out
mkdir -p $chimpanzee_common_tx_kallisto_out
mkdir -p $orangutan_common_tx_kallisto_out

$kallisto index -i $human_common_tx_kallisto_index $human_common_tx_fasta & 
$kallisto index -i $chimpanzee_common_tx_kallisto_index $chimpanzee_common_tx_fasta &
$kallisto index -i $orangutan_common_tx_kallisto_index $orangutan_common_tx_fasta &


for prefix in $human_fastq_prefix
do
mkdir -p $human_common_tx_kallisto_out/$prefix/
$kallisto quant -i $human_common_tx_kallisto_index -o $human_common_tx_kallisto_out/$prefix/ -t $num_threads -b 50 --bias --rf-stranded "$human_fastq_dir""$prefix""_1.fastq.gz" "$human_fastq_dir""$prefix""_2.fastq.gz"
done

for prefix in $chimpanzee_fastq_prefix
do
mkdir -p $chimpanzee_common_tx_kallisto_out/$prefix/
$kallisto quant -i $chimpanzee_common_tx_kallisto_index -o $chimpanzee_common_tx_kallisto_out/$prefix/ -t $num_threads -b 50 --bias --rf-stranded "$chimpanzee_fastq_dir""$prefix""_1.fastq.gz" "$chimpanzee_fastq_dir""$prefix""_2.fastq.gz"
done

for prefix in $orangutan_fastq_prefix
do
mkdir -p $orangutan_common_tx_kallisto_out/$prefix/
$kallisto quant -i $orangutan_common_tx_kallisto_index -o $orangutan_common_tx_kallisto_out/$(echo $prefix | awk -F"_" '{ print $1"_"$2}')/ -t $num_threads -b 50 --bias --rf-stranded "$orangutan_fastq_dir""$prefix""_1.fastq.gz" "$orangutan_fastq_dir""$prefix""_2.fastq.gz"
done




## run junctioncounts

human_junctioncounts_out=$human_out/junctioncounts/
chimpanzee_junctioncounts_out=$chimpanzee_out/junctioncounts/
orangutan_junctioncounts_out=$orangutan_out/junctioncounts/

infer_pairwise_events=/public/home/anjowall/repos/junctionCounts/junctionCounts/infer_pairwise_events.py
junctioncounts=/public/home/anjowall/repos/junctionCounts/junctionCounts/junctionCounts.py

python2.7 $infer_pairwise_events \
	--transcript_gtf $human_tx_merged \
	--outdir $human_junctioncounts_out &

python2.7 $infer_pairwise_events \
	--transcript_gtf $chimpanzee_tx_merged \
	--outdir $chimpanzee_junctioncounts_out &

python2.7 $infer_pairwise_events \
	--transcript_gtf $orangutan_tx_merged \
	--outdir $orangutan_junctioncounts_out &


gffread $human_junctioncounts_out/splice_lib_events.gtf -g $human_genome_fasta -w $human_junctioncounts_out/splice_lib_events.fa &
gffread $chimpanzee_junctioncounts_out/splice_lib_events.gtf -g $chimpanzee_genome_fasta -w $chimpanzee_junctioncounts_out/splice_lib_events.fa &
gffread $orangutan_junctioncounts_out/splice_lib_events.gtf -g $orangutan_genome_fasta -w $orangutan_junctioncounts_out/splice_lib_events.fa &


junctioncounts_fun () {
python2.7 $junctioncounts \
--event_gtf $1 \
--bam $2 \
--forward_read R2 \
--outdir $3 \
--sample_name $4 \
--event_ioe $5 \
--calc_gene_frac 
}



unset i

(
for bam in $human_bam/*.bam
do
((i=i%num_threads)); ((i++==0)) && wait
junctioncounts_fun $human_junctioncounts_out/splice_lib_events.gtf $bam $human_junctioncounts_out $(echo $bam | awk -F\/ '{ print $NF}' | awk -F"_" '{ print $1"_"$2}') $human_junctioncounts_out/splice_lib_events.ioe &
done
wait
)

unset i

(
for bam in $chimpanzee_bam/*.bam
do
((i=i%num_threads)); ((i++==0)) && wait
junctioncounts_fun $chimpanzee_junctioncounts_out/splice_lib_events.gtf $bam $chimpanzee_junctioncounts_out $(echo $bam | awk -F\/ '{ print $NF}' | awk -F"_" '{ print $1"_"$2"_"$3}') $chimpanzee_junctioncounts_out/splice_lib_events.ioe &
done
wait
)

unset i

(
for bam in $orangutan_bam/*.bam
do
((i=i%num_threads)); ((i++==0)) && wait
junctioncounts_fun $orangutan_junctioncounts_out/splice_lib_events.gtf $bam $orangutan_junctioncounts_out $(echo $bam | awk -F\/ '{ print $NF}' | awk -F"_" '{ print $1"_"$2}') $orangutan_junctioncounts_out/splice_lib_events.ioe &
done
wait
)




### concatenate junctioncounts output

for i in $human_junctioncounts_out/*count_psi_outfile.tsv
do
awk 'NR>1' "$i" >> $human_junctioncounts_out/temp
done

cat <( ls $human_junctioncounts_out/*count_psi_outfile.tsv | head -n 1 | xargs head -n 1 ) $human_junctioncounts_out/temp > $human_junctioncounts_out/human_all_junctioncounts.tsv

rm $human_junctioncounts_out/temp
rm $human_junctioncounts_out/*count_psi_outfile.tsv


for i in $chimpanzee_junctioncounts_out/*count_psi_outfile.tsv
do
awk 'NR>1' "$i" >> $chimpanzee_junctioncounts_out/temp
done

cat <( ls $chimpanzee_junctioncounts_out/*count_psi_outfile.tsv | head -n 1 | xargs head -n 1 ) $chimpanzee_junctioncounts_out/temp > $chimpanzee_junctioncounts_out/chimpanzee_all_junctioncounts.tsv

rm $chimpanzee_junctioncounts_out/temp
rm $chimpanzee_junctioncounts_out/*count_psi_outfile.tsv


for i in $orangutan_junctioncounts_out/*count_psi_outfile.tsv
do
awk 'NR>1' "$i" >> $orangutan_junctioncounts_out/temp
done

cat <( ls $orangutan_junctioncounts_out/*count_psi_outfile.tsv | head -n 1 | xargs head -n 1 ) $orangutan_junctioncounts_out/temp > $orangutan_junctioncounts_out/orangutan_all_junctioncounts.tsv

rm $orangutan_junctioncounts_out/temp
rm $orangutan_junctioncounts_out/*count_psi_outfile.tsv



## summarize by species (run calc_dpsi.R)

/public/home/anjowall/anaconda2/bin/Rscript /public/home/anjowall/repos/junctionCounts/scripts/calc_dpsi.R \
--all_junctioncounts $human_junctioncounts_out/human_all_junctioncounts.tsv \
--sample_info $human_out/sample_info.tsv \
--outdir $human_junctioncounts_out \
--comparisons fraction \
--levels fraction:1-2-3-4-5 \
--num_threads $num_threads


/public/home/anjowall/anaconda2/bin/Rscript /public/home/anjowall/repos/junctionCounts/scripts/calc_dpsi.R \
--all_junctioncounts $chimpanzee_junctioncounts_out/chimpanzee_all_junctioncounts.tsv \
--sample_info $chimpanzee_out/sample_info.tsv \
--outdir $chimpanzee_junctioncounts_out \
--comparisons fraction \
--levels fraction:1-2-3-4-5 \
--num_threads $num_threads


/public/home/anjowall/anaconda2/bin/Rscript /public/home/anjowall/repos/junctionCounts/scripts/calc_dpsi.R \
--all_junctioncounts $orangutan_junctioncounts_out/orangutan_all_junctioncounts.tsv \
--sample_info $orangutan_out/sample_info.tsv \
--outdir $orangutan_junctioncounts_out \
--comparisons fraction \
--levels fraction:1-2-3-4-5 \
--num_threads $num_threads


## get genes for events

python2.7 /public/home/anjowall/repos/multi_species_as_events/assign_events_to_genes.py \
--event_gtf $human_junctioncounts_out/splice_lib_events.gtf \
--gene_gtf $human_tx_gtf \
--transcript_gtf $human_tx_merged \
--prefix human \
--outdir $human_junctioncounts_out &


python2.7 /public/home/anjowall/repos/multi_species_as_events/assign_events_to_genes.py \
--event_gtf $chimpanzee_junctioncounts_out/splice_lib_events.gtf \
--gene_gtf $chimpanzee_tx_gtf \
--transcript_gtf $chimpanzee_tx_merged \
--prefix chimpanzee \
--outdir $chimpanzee_junctioncounts_out &


python2.7 /public/home/anjowall/repos/multi_species_as_events/assign_events_to_genes.py \
--event_gtf $orangutan_junctioncounts_out/splice_lib_events.gtf \
--gene_gtf $orangutan_tx_gtf \
--transcript_gtf $orangutan_tx_merged \
--prefix orangutan \
--outdir $orangutan_junctioncounts_out &


## run cds insertion

#### get CCDS tx ids from gencode



awk '$3=="transcript"' $annotation_base/gencode.v32.annotation.gtf | \
grep "CCDS" | awk '{ print $12}' | tr -d '\";' | sort | uniq > \
$annotation_base/ccds_tx_ids.txt

fgrep -f $annotation_base/ccds_tx_ids.txt $annotation_base/gencode.v32.annotation.gtf > $annotation_base/gencode.v32.ccds.gtf
fgrep -f $annotation_base/ccds_tx_ids.txt $annotation_base/chimp/panTro6/panTro6_cat_hg38_V32_Lv3_v2.gtf > $annotation_base/Chimp.ccds.gtf
fgrep -f $annotation_base/ccds_tx_ids.txt $annotation_base/orang/ponAbe3/ponAbe3_cat_hg38_V32_Lv3_v2.gtf > $annotation_base/Orangutan.ccds.gtf

### cannot set CCDS since CAT GTFs do not retain tags

python2.7 /public/home/anjowall/repos/cds_insertion/cds_insertion/cds_insertion.py \
	--transcript_gtf $human_tx_merged \
	--transcript_fasta $human_tx_merged_fasta \
	--annotation_gtf $annotation_base/gencode.v32.ccds.gtf \
	--bigGenePred_as_path /public/groups/sanfordlab/people/anjowall/bin/bigGenePred.as \
	--gtfToGenePred_path /public/groups/sanfordlab/people/anjowall/bin/gtfToGenePred \
	--genePredToBigGenePred_path /public/groups/sanfordlab/people/anjowall/bin/genePredToBigGenePred \
	--bedToBigBed_path /public/groups/sanfordlab/people/anjowall/bin/bedToBigBed \
	--bedtools_path /public/home/anjowall/bedtools2/bin/bedtools \
	--chrNameLength_path /public/groups/sanfordlab/people/juphilip/datasets/references/human/hg39/star_genome/chrNameLength.txt \
	--make_bigBed \
	--outdir $human_stringtie_out > $human_stringtie_out/cds_ins_stdout.txt 2> $human_stringtie_out/cds_ins_stderr.txt &


python2.7 /public/home/anjowall/repos/cds_insertion/cds_insertion/cds_insertion.py \
	--transcript_gtf $chimpanzee_tx_merged \
	--transcript_fasta $chimpanzee_tx_merged_fasta \
	--annotation_gtf $annotation_base/Chimp.ccds.gtf \
	--bigGenePred_as_path /public/groups/sanfordlab/people/anjowall/bin/bigGenePred.as \
	--gtfToGenePred_path /public/groups/sanfordlab/people/anjowall/bin/gtfToGenePred \
	--genePredToBigGenePred_path /public/groups/sanfordlab/people/anjowall/bin/genePredToBigGenePred \
	--bedToBigBed_path /public/groups/sanfordlab/people/anjowall/bin/bedToBigBed \
	--bedtools_path /public/home/anjowall/bedtools2/bin/bedtools \
	--chrNameLength_path /public/groups/sanfordlab/people/juphilip/datasets/references/chimp/panTro6/star_genome/chrNameLength.txt \
	--make_bigBed \
	--outdir $chimpanzee_stringtie_out > $chimpanzee_stringtie_out/cds_ins_stdout.txt 2> $chimpanzee_stringtie_out/cds_ins_stderr.txt &

python2.7 /public/home/anjowall/repos/cds_insertion/cds_insertion/cds_insertion.py \
	--transcript_gtf $orangutan_tx_merged \
	--transcript_fasta $orangutan_tx_merged_fasta \
	--annotation_gtf $annotation_base/Orangutan.ccds.gtf \
	--bigGenePred_as_path /public/groups/sanfordlab/people/anjowall/bin/bigGenePred.as \
	--gtfToGenePred_path /public/groups/sanfordlab/people/anjowall/bin/gtfToGenePred \
	--genePredToBigGenePred_path /public/groups/sanfordlab/people/anjowall/bin/genePredToBigGenePred \
	--bedToBigBed_path /public/groups/sanfordlab/people/anjowall/bin/bedToBigBed \
	--bedtools_path /public/home/anjowall/bedtools2/bin/bedtools \
	--chrNameLength_path /public/groups/sanfordlab/people/juphilip/datasets/references/orang/ponAbe3/star_genome_ponAbe2/chrNameLength.txt \
	--make_bigBed \
	--outdir $orangutan_stringtie_out > $orangutan_stringtie_out/cds_ins_stdout.txt 2> $orangutan_stringtie_out/cds_ins_stderr.txt &

## run find switch events


python2.7 ~/repos/cds_insertion/cds_insertion/find_switch_events.py \
	--outdir $human_junctioncounts_out \
	--ioe_file $human_junctioncounts_out/splice_lib_events.ioe \
	--event_gtf $human_junctioncounts_out/splice_lib_events.gtf \
	--transcript_dict_pkl $human_stringtie_out/cds_insertion_transcript_dict.pkl \
	 > $human_junctioncounts_out/find_switch_stdout.txt 2> $human_junctioncounts_out/find_switch_stderr.txt &

python2.7 ~/repos/cds_insertion/cds_insertion/find_switch_events.py \
	--outdir $chimpanzee_junctioncounts_out \
	--ioe_file $chimpanzee_junctioncounts_out/splice_lib_events.ioe \
	--event_gtf $chimpanzee_junctioncounts_out/splice_lib_events.gtf \
	--transcript_dict_pkl $chimpanzee_stringtie_out/cds_insertion_transcript_dict.pkl \
	 > $chimpanzee_junctioncounts_out/find_switch_stdout.txt 2> $chimpanzee_junctioncounts_out/find_switch_stderr.txt &

python2.7 ~/repos/cds_insertion/cds_insertion/find_switch_events.py \
	--outdir $orangutan_junctioncounts_out \
	--ioe_file $orangutan_junctioncounts_out/splice_lib_events.ioe \
	--event_gtf $orangutan_junctioncounts_out/splice_lib_events.gtf \
	--transcript_dict_pkl $orangutan_stringtie_out/cds_insertion_transcript_dict.pkl \
	 > $orangutan_junctioncounts_out/find_switch_stdout.txt 2> $orangutan_junctioncounts_out/find_switch_stderr.txt &


## summarize by species

## assess event support

## cross species annotation


paftools_path=/public/home/anjowall/minimap2/misc/paftools_mod.js

human_chimpanzee_paf=$annotation_base/human_chimpanzee_genome_alignment/human_to_chimpanzee_ref.paf
human_orangutan_paf=$annotation_base/human_orangutan_genome_alignment/human_to_orangutan_ref.paf


chimpanzee_human_paf=$annotation_base/chimpanzee_human_genome_alignment/chimpanzee_to_human_ref.paf
chimpanzee_orangutan_paf=$annotation_base/chimpanzee_orangutan_genome_alignment/chimpanzee_to_orangutan_ref.paf


orangutan_human_paf=$annotation_base/orangutan_human_genome_alignment/orangutan_to_human_ref.paf
orangutan_chimpanzee_paf=$annotation_base/orangutan_chimpanzee_genome_alignment/orangutan_to_chimpanzee_ref.paf


multi_species_events_standalone="/public/home/anjowall/repos/multi_species_as_events/multi_species_events_standalone.py"



python2.7 $multi_species_events_standalone \
--event_file_prefixes $human_junctioncounts_out/splice_lib_events,$chimpanzee_junctioncounts_out/splice_lib_events \
--species_list human,chimpanzee \
--indir / \
--genome_fasta_paths $human_genome_fasta,$chimpanzee_genome_fasta \
--paf_paths $human_chimpanzee_paf,$chimpanzee_human_paf \
--paftools_lift_outdir $topdir/human_chimpanzee_cross_events/ \
--paftools_path $paftools_path \
--add_gene_symbols \
--gene_gtfs $human_full_gtf,$chimpanzee_full_gtf \
--transcript_gtfs $human_tx_merged,$chimpanzee_tx_merged &



python2.7 $multi_species_events_standalone \
--event_file_prefixes $human_junctioncounts_out/splice_lib_events,$orangutan_junctioncounts_out/splice_lib_events \
--species_list human,orangutan \
--indir / \
--genome_fasta_paths $human_genome_fasta,$orangutan_genome_fasta \
--paf_paths $human_orangutan_paf,$orangutan_human_paf \
--paftools_lift_outdir $topdir/human_orangutan_cross_events/ \
--paftools_path $paftools_path \
--add_gene_symbols \
--gene_gtfs $human_full_gtf,$orangutan_full_gtf \
--transcript_gtfs $human_tx_merged,$orangutan_tx_merged &



python2.7 $multi_species_events_standalone \
--event_file_prefixes $chimpanzee_junctioncounts_out/splice_lib_events,$orangutan_junctioncounts_out/splice_lib_events \
--species_list chimpanzee,orangutan \
--indir / \
--genome_fasta_paths $chimpanzee_genome_fasta,$orangutan_genome_fasta \
--paf_paths $chimpanzee_orangutan_paf,$orangutan_chimpanzee_paf \
--paftools_lift_outdir $topdir/chimpanzee_orangutan_cross_events/ \
--paftools_path $paftools_path \
--add_gene_symbols \
--gene_gtfs $chimpanzee_full_gtf,$orangutan_full_gtf \
--transcript_gtfs $chimpanzee_tx_merged,$orangutan_tx_merged &


python2.7 $multi_species_events_standalone \
--event_file_prefixes $human_junctioncounts_out/splice_lib_events,$chimpanzee_junctioncounts_out/splice_lib_events,$orangutan_junctioncounts_out/splice_lib_events \
--species_list human,chimpanzee,orangutan \
--indir / \
--genome_fasta_paths $human_genome_fasta,$chimpanzee_genome_fasta,$orangutan_genome_fasta \
--paf_paths $human_chimpanzee_paf,$human_orangutan_paf,$chimpanzee_human_paf,$chimpanzee_orangutan_paf,$orangutan_human_paf,$orangutan_chimpanzee_paf \
--paftools_lift_outdir $topdir/all_species_cross_events/ \
--paftools_path $paftools_path \
--add_gene_symbols \
--gene_gtfs $human_full_gtf,$chimpanzee_full_gtf,$orangutan_full_gtf \
--transcript_gtfs $human_tx_merged,$chimpanzee_tx_merged,$orangutan_tx_merged


python2.7 $multi_species_events_standalone \
--event_file_prefixes $chimpanzee_junctioncounts_out/splice_lib_events,$orangutan_junctioncounts_out/splice_lib_events,$macaque_junctioncounts_out/splice_lib_events \
--species_list chimpanzee,orangutan,macaque \
--indir / \
--genome_fasta_paths $chimpanzee_genome_fasta,$orangutan_genome_fasta,$macaque_genome_fasta \
--paf_paths $chimpanzee_orangutan_paf,$chimpanzee_macaque_paf,$orangutan_chimpanzee_paf,$orangutan_macaque_paf,$macaque_chimpanzee_paf,$macaque_orangutan_paf \
--paftools_lift_outdir $topdir/chimpanzee_orangutan_macaque_cross_events/ \
--paftools_path $paftools_path \
--add_gene_symbols \
--gene_gtfs $chimpanzee_full_gtf,$orangutan_full_gtf,$macaque_full_gtf \
--transcript_gtfs $chimpanzee_tx_merged,$orangutan_tx_merged,$macaque_tx_merged

python2.7 $multi_species_events_standalone \
--event_file_prefixes $human_junctioncounts_out/splice_lib_events,$chimpanzee_junctioncounts_out/splice_lib_events,$orangutan_junctioncounts_out/splice_lib_events \
--species_list human,chimpanzee,orangutan \
--indir / \
--genome_fasta_paths $human_genome_fasta,$chimpanzee_genome_fasta,$orangutan_genome_fasta \
--paf_paths $human_chimpanzee_paf,$human_orangutan_paf,$chimpanzee_human_paf,$chimpanzee_orangutan_paf,$orangutan_human_paf,$orangutan_chimpanzee_paf \
--paftools_lift_outdir $topdir/human_chimpanzee_orangutan_cross_events/ \
--paftools_path $paftools_path \
--add_gene_symbols \
--gene_gtfs $human_full_gtf,$chimpanzee_full_gtf,$orangutan_full_gtf \
--transcript_gtfs $human_tx_merged,$chimpanzee_tx_merged,$orangutan_tx_merged




### run junctionCounts with new IDs from all species run

all_species_event_annot_dir=$topdir/all_species_cross_events/
all_species_jc_out=$all_species_event_annot_dir/junctioncounts/

mkdir -p $all_species_jc_out

#### redefine jc fun to eliminate --calc_gene_frac, which would required IOE files

junctioncounts_fun () {
python2.7 $junctioncounts \
--event_gtf $1 \
--bam $2 \
--forward_read R2 \
--outdir $3 \
--sample_name $4
}

## generate event GTFs restricted to events annotated in all species

human_restricted_events=$all_species_event_annot_dir/human_splice_lib_events_restricted.gtf
chimpanzee_restricted_events=$all_species_event_annot_dir/chimpanzee_splice_lib_events_restricted.gtf
orangutan_restricted_events=$all_species_event_annot_dir/orangutan_splice_lib_events_restricted.gtf

grep "human_chimpanzee_orangutan" $all_species_event_annot_dir/human_splice_lib_events.gtf \
> $human_restricted_events
grep "human_chimpanzee_orangutan" $all_species_event_annot_dir/chimpanzee_splice_lib_events.gtf \
> $chimpanzee_restricted_events
grep "human_chimpanzee_orangutan" $all_species_event_annot_dir/orangutan_splice_lib_events.gtf \
> $orangutan_restricted_events


unset i

(
for bam in $human_bam/*.bam
do
((i=i%num_threads)); ((i++==0)) && wait
junctioncounts_fun $human_restricted_events $bam $all_species_jc_out $(echo $bam | awk -F\/ '{ print $NF}' | awk -F"_" '{ print $1"_"$2}') &
done
wait
)

unset i

(
for bam in $chimpanzee_bam/*.bam
do
((i=i%num_threads)); ((i++==0)) && wait
junctioncounts_fun <(grep "human_chimpanzee_orangutan" $all_species_event_annot_dir/chimpanzee_splice_lib_events.gtf) $bam $all_species_jc_out $(echo $bam | awk -F\/ '{ print $NF}' | awk -F"_" '{ print $1"_"$2"_"$3}') &
done
wait
)

unset i

(
for bam in $orangutan_bam/*.bam
do
((i=i%num_threads)); ((i++==0)) && wait
junctioncounts_fun <(grep "human_chimpanzee_orangutan" $all_species_event_annot_dir/orangutan_splice_lib_events.gtf) $bam $all_species_jc_out $(echo $bam | awk -F\/ '{ print $NF}' | awk -F"_" '{ print $1"_"$2}') &
done
wait
)



### concatenate junctioncounts output

for i in $all_species_jc_out/*count_psi_outfile.tsv
do
awk 'NR>1' "$i" >> $all_species_jc_out/temp
done

cat <( ls $all_species_jc_out/*count_psi_outfile.tsv | head -n 1 | xargs head -n 1 ) $all_species_jc_out/temp > $all_species_jc_out/all_species_junctioncounts.tsv

rm $all_species_jc_out/temp
rm $all_species_jc_out/*count_psi_outfile.tsv



### create sample_info for all samples with species

cat <(echo -e "sample_name\tfraction\tspecies") \
<(awk 'NR>1 && $2 <= 6 {OFS="\t"; print $1,$2,"human"}' $human_out/sample_info.tsv) \
<(awk 'NR>1 {OFS="\t"; print $1,$2,"chimpanzee"}' $chimpanzee_out/sample_info.tsv) \
<(awk 'NR>1 {OFS="\t"; print $1,$2,"orangutan"}' $orangutan_out/sample_info.tsv) > \
$all_species_jc_out/all_species_sample_info.tsv


## summarize by species (run calc_dpsi.R)

cat <( awk '{ print $10}' human_splice_lib_events_restricted.gtf | sort | uniq ) <(awk '{ print $10}' chimpanzee_splice_lib_events_restricted.gtf | sort | uniq ) <(awk '{ print $10}' orangutan_splice_lib_events_restricted.gtf | sort | uniq ) | sort | uniq -c | tr -d '\";' | awk '$1 == 4 {print $2}' > common_event_ids.txt


fgrep -f ../common_event_ids.txt all_species_junctioncounts.tsv > temp

cat <(head -n 1 all_species_junctioncounts.tsv ) temp > all_species_junctioncounts_restricted.tsv
rm temp

/public/home/anjowall/anaconda2/bin/Rscript /public/home/anjowall/repos/junctionCounts/scripts/calc_dpsi.R \
--all_junctioncounts $all_species_jc_out/all_species_junctioncounts_restricted.tsv \
--sample_info $all_species_jc_out/all_species_sample_info.tsv \
--outdir $all_species_jc_out \
--comparisons species:fraction,fraction:species \
--levels species:human-chimpanzee-orangutan,fraction:1-2-3-4-5 \
--num_threads $num_threads \
--calc_distance \
--calc_ddpsi


python2.7 /public/home/anjowall/repos/splice_lib_utils/get_competing_splice_sites.py \
--outdir $all_species_jc_out \
--species human,chimpanzee,orangutan \
--event_gtfs $human_restricted_events,$chimpanzee_restricted_events,$orangutan_restricted_events \
--genome_fastas $human_genome_fasta,$chimpanzee_genome_fasta,$orangutan_genome_fasta \
--region_bed \
--exonic_length 50 \
--intronic_length 250 \
--bedtools_path /public/home/anjowall/bedtools2/bin/bedtools

bwtool="/public/home/anjowall/bin/bwtool"

$bwtool extract bed $all_species_jc_out/human_splice_site_regions.bed \
/public/groups/sanfordlab/people/anjowall/genomes/hg38/hg38.phastCons100way.bw \
$all_species_jc_out/human_splice_site_regions_phastcons.tsv -tabs

$bwtool extract bed human_splice_site_regions.bed \
/public/groups/sanfordlab/people/anjowall/genomes/hg38/hg38.phyloP100way.bw \
$all_species_jc_out/human_splice_site_regions_phylop.tsv -tabs




python2.7 /public/home/anjowall/repos/multi_species_as_events/assign_events_to_genes.py \
--event_gtf $all_species_event_annot_dir/human_splice_lib_events.gtf \
--gene_gtf $human_full_gtf \
--transcript_gtf $human_tx_merged \
--prefix human \
--outdir $all_species_event_annot_dir &


