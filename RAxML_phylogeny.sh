### RAxML tree
raxmlHPC-SSE3 -s Alignments/Cadophora_ILRT_modified_2020.fasta -n ITS_LSU_RPB2_TEF -f a -m GTRGAMMA -q Alignments/partitions.txt -x 12345 -p 12345 -# 1000
mv RAxML_* RAxML_tree/

