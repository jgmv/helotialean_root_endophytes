# get scripts from GitHub
mkdir -p scripts
if test -f ./scripts/process_fungal_ITS-LSU/fungal_ITS-LSU_process.sh
then
  echo "fungal_ITS-LSU_process.sh found"
else
  git clone https://github.com/jgmv/process_fungal_ITS-LSU.git \
    scripts/process_fungal_ITS-LSU
fi

if test -f ./scripts/python_seq_analysis/subset_fasta_seqs.sh
then
  echo "python_seq_analysis.sh found"
else
  git clone https://github.com/jgmv/python_seq_analysis.git \
    scripts/python_seq_analysis
  for i in scripts/python_seq_analysis/*.py
  do
    chmod +x $i
  done
fi

if test -f ./scripts/ncbi_data_analysis/fetch_gb.py
then
  echo "ncbi_data_analysis found"
else
  git clone https://github.com/jgmv/ncbi_data_analysis.git \
    scripts/ncbi_data_analysis
  for i in scripts/ncbi_data_analysis/*.py
  do
    chmod +x $i
  done
fi


export PATH="$PATH:scripts/process_fungal_ITS-LSU:scripts/python_seq_analysis:scripts/ncbi_data_analysis"

# identify and perform phylogenetic analysis (if not done yet)
if test -f ./multilocus_phylogeny/RAxML_bipartitions.ITS_LSU
then
  echo "analysis done"
else
  bash fungal_ITS-LSU_process.sh ITS.fasta LSU.fasta
fi

### code above this can be omitted

# manually select representative strains
declare -a arr=("P1162" "P1176" "P1280" "P1331" "P1323" "P1348" "P1381" \
 "P1396" "P1615" "P1854" "P1963" "P1984" "P2010" "P2414" "P2437" "P2739" \
 "P2794" "P2816" "P6045" "P6089" "P6112" "P6587")

for i in "${arr[@]}"
do
  echo $i >> strain_selection.txt
  grep $i otu_clustering/otu_list.csv
done

subset_fasta_seqs.py ITS.fasta strain_selection.txt -o strain_selection.fasta
rm strain_selection.txt

# standalone BLAST of ITS sequences against NCBI
blastn -query strain_selection.fasta -db fungalITS -out blast_result.txt \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sgi"
sed 's/^[^|]*|\([^|]*\)|.*/\1/g' blast_result.txt > gi_list.txt
fetch_gb.py gi_list.txt -o gb_records.gb -email "MaciaVicente@em.uni-frankfurt.de"
get_metadata_from_gb.py gb_records.gb -o gb_records.csv

# standalone BLAST of ITS sequences against selected SRA objects
ls /home/jose/databases/sra/*.sra > db.txt
cut -d'/' -f6 db.txt | cut -d'.' -f1 > databases.txt
rm db.txt

echo -n > sra_blast_output.csv
IFS=$'\n'
set -f
for i in $(cat < databases.txt)
do
  echo "BLAST against "$i
  blastn_vdb -db  $i -query strain_selection.fasta -out $i_blast_result.csv \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
  cat $i_blast_result.csv >> sra_blast_output.csv
  rm $i_blast_result.csv
done

# keep only matches with pident > 98
cat sra_blast_output.csv | awk '$3 > 98' > sra_blast_output_filtered.csv







