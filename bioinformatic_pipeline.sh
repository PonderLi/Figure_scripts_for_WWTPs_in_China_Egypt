#! /bin/bash
#### qc
KNEADDATA_DB_HUMAN_GENOME=/share/soft/miniconda3/envs/kneaddata/db/
kneaddata -i $fq1 -i $fq2 --output-prefix $id \
      -o ./ -v -p 16 --remove-intermediate-output \
      --trimmomatic /share/soft/miniconda3/envs/kneaddata/bin/../share/trimmomatic-0.39-2/ \
      --trimmomatic-options 'SLIDINGWINDOW:4:20 MINLEN:50' \
      --bowtie2-options '--very-sensitive --dovetail' \
      -db $KNEADDATA_DB_HUMAN_GENOME \
      --fastqc /share/soft/miniconda3/envs/kneaddata/bin/ \
      --run-fastqc-start --run-fastqc-end --reorder 
      
#### assembly
megahit -1 $fq1 \
        -2 $fq2 \
        --min-contig-len 1000\
        -m 30000000000 \
        -t 16  \
        --presets meta-sensitive  \
        --out-dir $name
        
#### reads analysis for ARGs with ARGs-OAP
source /share/soft/miniconda3/bin/activate /share/soft/miniconda3/envs/argoaps
args_oap stage_one -i data  -o stage_1 -f fastq.gz -t 16  && \
args_oap stage_two -i stage_1   -o stage_2 -t 16

#### orf predict 
$seqkit replace -p .+ -r "${name}_{nr}" $name/final.contigs.fa -o $name/$name.contig.fa
$pprodigal -i $name/$name.contig.fa -g 11 -a $name/${name}.gene.faa \
    -d $name/${name}.gene.fna -o $name/${name}.gene.gff -f gff -p meta -T 16 -C 30000
    
#### non-redudant gene catolog
seqkit='/share/soft/miniconda3/envs/denovo/bin/seqkit'
mmseqs easy-linclust all_gene.faa all_uniq tmp\
       -e 0.001 --min-seq-id 0.95 -c 0.80 \
       --threads 16 

sed -i 's/ .*//' all_uniq_rep_seq.fasta
$seqkit seq -m 33 all_uniq_rep_seq.fasta -o all_uniq_rep_seq_100.fasta
$seqkit seq -n all_uniq_rep_seq_100.fasta |$seqkit grep -f - all_gene.fna > all_uniq_gene_100.fna 

sed -i 's/ .*//' all_uniq_gene_100.fna
sed -i 's/\*//' all_uniq_rep_seq_100.fasta
mv all_uniq_rep_seq_100.fasta all_uniq_gene_100.faa

#### identify ARGs
diamond blastp \
        --db /home/database/sarg/Full_database/full_sarg.dmnd \
        --query tes.faa \
        -o sarg_result2/sarg_result.tsv --more-sensitive \
        --evalue 1e-5 --id 80 --query-cover 70 \
        --max-target-seqs 1 --header -f 6 --threads 32
rgi main --input_sequence  tes.faa\
    --output_file card_result/card.txt \
    --input_type protein \
    --clean -a DIAMOND -n 32
deeparg predict --model LS \
    -i tes.faa \
    -o deeparg_result/deeparg.txt \
    -d /home/HuAnyi/database/deeparg/ \
    --type prot --min-prob 0.8 \
    --arg-alignment-identity 50 \
    --arg-alignment-evalue 1e-10 \
    --arg-num-alignments-per-entry 1000

#### coverm
coverm contig   -r ./index/all_uniq_gene_100.fna \
                -1 data/*paired_1.fastq.gz \
                -2 data/*paired_2.fastq.gz\
                --output-file abundance_gene.csv\
                -p bwa-mem \
                --min-read-percent-identity 95 \
                --min-read-aligned-percent 75  \
                --methods count mean trimmed_mean rpkm tpm\
                -t 16 --no-zeros \
                --exclude-supplementary --min-covered-fraction 10 
                
#### identify pathogenic contigs
clarkpath='/home/contig_pathogen/'
CLARK -n 16 -k 31 -m 1 \
        -T $clarkpath/targets_addresses.txt \
        -D $clarkpath/DBD/ \
        -O all_contigs_rep_seq.fasta \
        -R ./results_k31_m1


#### identify plasmid contig
genomad end-to-end --threads 32 --cleanup \
        --splits 0 -s 8  \
        arg_contig.fa genomad_7_output \
        /home/HuAnyi/database/genomad_2/genomad_db 
python /home/database/PLASme/PLASMe/PLASMe.py\
       ../arg_all_contig.fa result \
      -d /home/database/PLASme/PLASMe/DB \
      -c 0.9 -i 0.9 -p 0.5 -t 32
DeepMicroClass predict -i ../genomad/all_contig_uniq_rep_seq.fasta -o DMC_new
#### identify viral contigs
genomad end-to-end --threads 32 --cleanup \
        --splits 0 -s 8  \
        arg_contig.fa genomad_7_output \
        /home/HuAnyi/database/genomad_2/genomad_db 
phabox2 --task phamer --dbdir ~/database/phabox/phabox_db_v2 \
        --outpth phamer --contigs ../genomad/all_contig_uniq_rep_seq.fasta \
        --len 990 --threads 32
virsorter run -i ../genomad/all_contig_uniq_rep_seq.fasta \ 
              --include-groups dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae \
              -j 32 --min-score 0.5
#### kairos
nextflow ~/database/kairos/kairos-main/kairos-dd.nf \
                 --max_cpus 8 --max_overlap 0.5 \
                 --input_contigs $1 --taxa_df /home/kairos/taxadf_all.tsv \
                 --outdir $2 \
                 --target_database ~/kairos/all_arg.dmnd \
                 --MGE_database ~/database/MobileOG/mobileOG-db_beatrix-1.6.All.faa.dmnd \
                 --num_chunks 8


