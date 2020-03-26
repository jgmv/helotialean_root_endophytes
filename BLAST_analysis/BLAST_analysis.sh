 #!/usr/bin/env bash


### prepare environment
# set variables
EMAIL="your@email" # tell NCBI who you are
BLAST_DB="fungalITS" # BLAST Database for comparisons
SRA_DB="path/to/SRA" # path to the folder where the SRA databases are stored

# create folders
mkdir -p output


### download necessary scripts
# custom python scripts
if test -f ./scripts/python_seq_analysis/subset_fasta_seqs.sh
then
  echo "python_seq_analysis.sh found"
else
  git clone https://github.com/jgmv/python_seq_analysis.git temp
  mv temp/* scripts
  rm -rf temp
  for i in scripts/*.sh
  do
    . $i
  done
fi

# custom scripts to interact with NCBI data
if test -f ./scripts/ncbi_data_analysis/fetch_gb.py
then
  echo "ncbi_data_analysis found"
else
  git clone https://github.com/jgmv/ncbi_data_analysis.git temp
  mv temp/* scripts
  rm -rf temp
  for i in scripts/*.py
  do
    chmod +x $i
  done
fi

# system-wide access to scripts
export PATH="$PATH:scripts"


### select isolates to include in the analysis
declare -a arr=("P1101" "P1162" "P1176" "P1280" "P1281" "P1313" "P1323" "P1331" "P1344" "P1348" "P1381" "P1396" "P1442" "P1443" "P1450" "P1451" "P1518" "P1572" "P1615" "P1686" "P1749" "P1751" "P1772" "P1823" "P1850" "P1854" "P1866" "P1893" "P1900" "P1909" "P1924" "P1935" "P1940" "P1944" "P1950" "P1959" "P1963" "P1973" "P1974" "P1978" "P1984" "P1992" "P1993" "P2010" "P2092" "P2118" "P2166" "P2298" "P2323" "P2414" "P2437" "P2440" "P2450" "P2477" "P2482" "P2485" "P2500" "P2509" "P2510" "P2511" "P2514" "P2515" "P2535" "P2641" "P2646" "P2739" "P2766" "P2775" "P2785" "P2792" "P2793" "P2794" "P2810" "P2816" "P2828" "P2832" "P2851" "P2867" "P2899" "P2925" "P2932" "P2962" "P2966" "P2973" "P2995" "P6033" "P6045" "P6077" "P6081" "P6089" "P6101" "P6112" "P6587")

for i in "${arr[@]}"
do
  echo $i >> strain_selection.txt
done

subset_fasta_seqs.py data/ITS.fasta strain_selection.txt -o \
  data/strain_selection.fasta
rm strain_selection.txt


### Run standalone BLAST against ITS sequences in NCBI GenBank
blastn -query data/strain_selection.fasta -db $BLAST_DB \
  -out output/blast_result.txt \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sgi"
sed 's/^[^|]*|\([^|]*\)|.*/\1/g' output/blast_result.txt > output/gi_list.txt
fetch_gb.py output/gi_list.txt -o output/gb_records.gb -email $EMAIL
get_metadata_from_gb.py output/gb_records.gb -o output/gb_records.csv


### Run standalone BLAST against selected SRA objects
ls $SRA_DB/*.sra > db.txt
cut -d'/' -f6 db.txt | cut -d'.' -f1 > output/databases.txt
rm db.txt

echo -n > output/sra_blast_output.csv
IFS=$'\n'
set -f
for i in $(cat < output/databases.txt)
do
  echo "BLAST against "$i
  blastn_vdb -db  $i -query data/strain_selection.fasta \
    -out $i_blast_result.csv \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
  cat $i_blast_result.csv >> output/sra_blast_output.csv
  rm $i_blast_result.csv
done

# keep only matches with pident > 98
cat output/sra_blast_output.csv | awk '$3 > 98' > \
  output/sra_blast_output_filtered.csv


### Process output data
Rscript --no-save scripts/process_BLAST_results.R


### end
