DATI 

mongo "mongodb://172.21.2.89:27017/disprot8"

la più aggiornata è entries_2020_12_c

oppure sull'interfaccia di Disprot


>> TRASFORMARE IL JSON IN FASTA

cat disprot_2020_12.json | jq -cr '.data[] | [.acc, .sequence] | @csv' | tr -d '"' | awk -F"," '{printf(">%s\n%s\n",$1,$2)}'

>> SPLITTARE IL FASTA IN SINGOLI FASTA

fasta file is in data/disprot
awk -F "|" '/^>/ {close(F); ID=$1; gsub("^>", "", ID); F=ID".fasta"} {print >> F}' yourfile.fa
rm .fasta # gets generated in the process, don't know what it is but I don't need it

>> PER L'ALLINEAMENTO GLOBALE
>> CREARE LA MATRICE TRIANGOLARE SENZA DIAGONALE
python pairwise.py # to create the pairwise couples


>> I have changed folders, YOU NEED TO CHANGE PATHS BEFORE EXECUTING THE BASH SCRIPTS
# dividi global_triang_matrix in liste da 5k coppie e crea lista di liste
mkdir -p lists sge-out results
./split-list.sh data/alignments_global_needle/global_triangular_matrix.csv
#rsync -avhz --delete lists/ marbev@echidna:IDPanalysis/lists/
# rsync -avhz  list_of_lists.dat marbev@echidna:IDPanalysis/

# lancia l'array job su sge locale
rm -rf sge-out/submit-needle-on-lists.sh* results/*
rm -f needle-summary.dat
qsub -t 1-$(cat list_of_lists.dat|wc -l) -o sge-out/ -e sge-out/ -cwd ./submit-needle-on-lists.sh list_of_lists.dat
# ... internamente chiamerà:
# ./needle-on-list.sh lists/list_*.dat

# per vedere come cambia la coda in tempo reale
 watch 'qstat -f -q super'

# colleziona i risultati
cat results/list_*.out

#>>>>>> LANCIARE EMBOSS NEEDLE SULLE COPPIE DI SEQUENZE

# WATER
rm -rf sge-out/submit-water-on-lists.sh* results/*
rm -f water-summary.dat
qsub -t 1-$(cat list_of_lists.dat|wc -l) -o sge-out/ -e sge-out/ -cwd ./submit-water-on-lists.sh list_of_lists.dat

# SPLITTARE IL FILE DELLE REGIONI DISORDINATE IN FILE MULTIPLI

awk '/^>/{if(file){close(file)};split($1,a,"|");split($2,b,"=");file=a[2] "_" b[2] ".fasta"} {print > file}'  inputFile.txt



>> NEEDLE SEQUENCES ANALYSIS

run needle
enter the results folder and copy results with scp marbev@echidna:IDPanalysis/results_needle/list* .
concat all results with cat * > all-global-needle.txt


>> PIPELINE DISORDER ANALYSIS AFTER NEEDLE UNIPROTS
# filter >= 30% identity
1. filter_identity(30) filter.py 
# generate clusters
2. cluster_file() cluster.py # starts to be useful only when the filtering is large, otherwise data are all connected in a big cluster
# keep uniprots with a structural state region associated
3. filter_structural() filter.py
# filter uniprots from clusters
4. filter_cluster() cluster.py
# rearrange data into dataframes
5. create_tables() rearrange.py
# calculate union overlap and visualize a graph
7. visualize_overlap_identity() visualize.py
# visualize histograms
8. visualize_overlap_identity_hist() visualize.py
#filter to visualize a subset of alignments
>> study of cases with 80-100% identity
8. visualize_alignments() visualize.py
