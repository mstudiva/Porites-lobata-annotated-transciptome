# Porites lobata Transcriptome Annotation, version December 18, 2024
# Created by Misha Matz (matz@utexas.edu), modified by Michael Studivan (studivanms@gmail.com) for use on FAU's HPC (KoKo)


#------------------------------
# BEFORE STARTING, replace, in this whole file:
#	- studivanms@gmail.com by your actual email;
#	- mstudiva with your KoKo user name.

# The idea is to copy the chunks separated by empty lines below and paste them into your cluster
# terminal window consecutively.

# The lines beginning with hash marks (#) are explanations and additional instructions -
# please make sure to read them before copy-pasting.


#------------------------------
# setup

# To install Bioperl as a conda environment
conda create -y -n bioperl perl-bioperl

# getting scripts
cd ~/bin
git clone https://github.com/z0on/annotatingTranscriptomes.git
mv annotatingTranscriptomes/* .
rm -rf annotatingTranscriptomes
rm launcher_creator.py

git clone https://github.com/z0on/emapper_to_GOMWU_KOGMWU.git
mv emapper_to_GOMWU_KOGMWU/* .
rm -rf emapper_to_GOMWU_KOGMWU

git clone https://github.com/mstudiva/Porites-lobata-annotated-transcriptome.git
mv Porites-lobata-annotated-transcriptome/* .
rm -rf Porites-lobata-annotated-transcriptome

# creating backup directory
mkdir backup

# creating annotation directory
cd
mkdir annotate
cd annotate


#------------------------------
# getting transcriptomes

# Noel (2023) https://doi.org/10.1186/s13059-023-02960-7
# cds from genome at https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_942486035.1/
# protein translations at https://www.genoscope.cns.fr/corals/genomes.html

# Renaming gene identifiers for ease
sed -i 's/lcl|CALNXK/Plobata/' Plobata.fasta
sed -i 's/PLOB_/Plobata/' Plobata.fasta
sed -i 's/Plob_/Plobata/' Plobata_pro.fasta

# transcriptome statistics
conda activate bioperl
echo "seq_stats.pl Plobata.fasta > seqstats_Plobata.txt" > seq_stats
launcher_creator.py -j seq_stats -n seq_stats -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch seq_stats.slurm

Plobata.fasta
-------------------------
42872 sequences.
1384 average length.
63084 maximum length.
105 minimum length.
N50 = 1953
59.3 Mb altogether (59326467 bp).
0 ambiguous Mb. (0 bp, 0%)
0 Mb of Ns. (0 bp, 0%)
-------------------------


#------------------------------
# Extracting contig and isogroup IDs into a lookup table

grep ">" Plobata.fasta | awk '{
    header = substr($1, 2);                  # Remove the leading ">" from the first field
    match($0, /locus_tag=([^]]+)/, gene);    # Extract the locus_tag (gene name)
    print header "\t" gene[1];               # Print the full header and locus_tag
}' > Plobata_seq2iso.tab


#------------------------------
# GO annotation
# updated based on Misha Matz's new GO and KOG annotation steps on GitHub: https://github.com/z0on/emapper_to_GOMWU_KOGMWU

# selecting the longest contig per isogroup
fasta2SBH.pl Plobata_pro.fasta >Plobata_pro_out.fasta

# scp *_pro_out.fasta to laptop, submit it to
http://eggnog-mapper.embl.de
cd /path/to/local/directory
scp mstudiva@koko-login.hpc.fau.edu:~/path/to/HPC/directory/\*_pro_out.fasta .

# copy link to job ID status and output file, paste it below instead of current link:
# http://eggnog-mapper.embl.de/job_status?jobname=MM_u7kolx5l

# once it is done, download results to HPC:
wget http://eggnog-mapper.embl.de/MM_u7kolx5l/out.emapper.annotations

# GO:
awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$10 }' out.emapper.annotations | grep GO | perl -pe 's/,/;/g' >Plobata_iso2go.tab
# gene names:
awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$8 }' out.emapper.annotations | grep -Ev "\tNA" >Plobata_iso2geneName.tab


#------------------------------
# KOG annotation
# updated based on Misha Matz's new GO and KOG annotation steps on github: https://github.com/z0on/emapper_to_GOMWU_KOGMWU

cp ~/bin/kog_classes.txt .

#  KOG classes (single-letter):
awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$7 }' out.emapper.annotations | grep -Ev "\tNA" >Plobata_iso2kogClass1.tab
# converting single-letter KOG classes to text understood by KOGMWU package (must have kog_classes.txt file in the same dir):
awk 'BEGIN {FS=OFS="\t"} NR==FNR {a[$1] = $2;next} {print $1,a[$2]}' kog_classes.txt Plobata_iso2kogClass1.tab > Plobata_iso2kogClass.tab


#------------------------------
# KEGG annotations:

# first rename fasta headers as gene ID (isogroup) rather than contig ID
awk 'BEGIN {
    while (getline < "Plobata_seq2iso.tab") {
        map[$1] = $2
    }
}
/^>/ {
    gene_id = substr($0, 2, index($0, " ") - 2)
    if (gene_id in map) {
        sub(gene_id, map[gene_id])
    }
}
{ print }' Plobata.fasta > Plobata_iso.fasta

# selecting the longest contig per isogroup
fasta2SBH.pl Plobata_iso.fasta >Plobata_iso_out.fasta

# Sanity check: How many unique isogroups do we have?
grep ">" Plobata.fasta | sort | uniq | wc -l            # 42872
grep ">" Plobata_pro.fasta | sort | uniq | wc -l        # 42872
grep ">" Plobata_pro_out.fasta | sort | uniq | wc -l    # 42871
grep ">" Plobata_iso.fasta | sort | uniq | wc -l        # 42872
grep ">" Plobata_iso_out.fasta | sort | uniq | wc -l    # 42871
# The two _out files should roughly match

# scp *_iso_out.fasta to your laptop
cd /path/to/local/directory
scp mstudiva@koko-login.hpc.fau.edu:~/path/to/HPC/directory/\*_iso_out.fasta .
# use web browser to submit _iso.fasta file to KEGG's KAAS server (http://www.genome.jp/kegg/kaas/)
# select SBH method, upload nucleotide query
https://www.genome.jp/kaas-bin/kaas_main?mode=user&id=1734545226&key=zsuB4sJ4

# Once it is done, download to HPC - it is named query.ko by default
wget https://www.genome.jp/tools/kaas/files/dl/1734545226/query.ko

# selecting only the lines with non-missing annotation:
cat query.ko | awk '{if ($2!="") print }' > Plobata_iso2kegg.tab

# the KEGG mapping result can be explored for completeness of transcriptome in terms of genes found,
# use 'html' output link from KAAS result page, see how many proteins you have for conserved complexes and pathways, such as ribosome, spliceosome, proteasome etc


#------------------------------
# file transfer

# copy all files to local machine
cd /path/to/local/directory
scp mstudiva@koko-login.fau.edu:~/path/to/HPC/directory/\* .
