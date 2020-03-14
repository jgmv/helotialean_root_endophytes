### download necessary scripts
# custom bash scripts
if test -f ./scripts/bash_seq_analysis/otuList.sh
then
  echo "./scripts/bash_seq_analysis found"
else
  git clone https://github.com/jgmv/bash_seq_analysis.git scripts/bash_seq_analysis
fi

for i in $(ls scripts/bash_seq_analysis/*.sh)
do
  . $i
done


### group ITS sequences into OTUs
# 97% similarity
blastclust -i its.fasta -o otus_97.csv -p F -b T -S 97
otuList otus_97.csv otu_list_97.csv

# 98% similarity
blastclust -i its.fasta -o otus_98.csv -p F -b T -S 98
otuList otus_98.csv otu_list_98.csv

# 99% similarity
blastclust -i its.fasta -o otus_99.csv -p F -b T -S 99
otuList otus_99.csv otu_list_99.csv

### end
