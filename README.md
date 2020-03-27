# helotialean_root_endophytes
Data and code for the data analyses used in [Maciá-Vicente, Piepenbring, Koukol (2020) Brassicaceous roots as an unexpected diversity hot-spot of helotialean endophytes. Preprint in Research Square, 10.21203/rs.3.rs-17858/v1](https://doi.org/10.21203/rs.3.rs-17858/v1).

## Contents
### BLAST_analysis
Code to compare the ITS sequence from the focal isolates of the study with records in NCBI GenBank and with selected Sequence Read Archive (SRA) objects. The script generates files `blast_nt.csv` and `blast_sra.csv` that are required for the data analyses detailed below.

In order to run this script, it is necessary to have a local BLAST database and the selected SRA objects for comparison, which can be obtained using the code available at [https://github.com/jgmv/local_BLAST_db](https://github.com/jgmv/local_BLAST_db).

This analysis can take several days to run given the number of comparisons, and the intermediate files generated have sizes >2~Gb. Therefore, the relevant files for downstream analyses are provided in [data_analysis](https://github.com/jgmv/helotialean_root_endophytes/data_analysis)

### OTU_clustering
Data and code to clasify to ITS sequences of the focal isolates into Operational Taxonomic Units (OTUs) at 97, 98, and 99 % sequence similarity thresholds. This results in files `otu_list_97.csv`, `otu_list_98.csv`, and `otu_list_99.csv` that are required during the data analyses to generate Fig.~S1 of the manuscript.
These output files are provided in [data_analysis](https://github.com/jgmv/helotialean_root_endophytes/data_analysis)

### data_analysis
Data and code for the statistical analyses in R. The following R packages are required: `ape`,`ggplot2`,`ggthemes`,`Hmisc`,`maps`,`PCAmixdata`, `rnaturalearth`, `sp`, and `vegan`.

## Additional files
* The alignments and phylogenetic trees generated in this study are available at TreeBASE under accession [S25942](https://www.treebase.org/treebase-web/search/study/summary.html?id=25942).
* The code for creating local BLAST files with the NCBI nt database and selected SRA objects is provided at [https://github.com/jgmv/local_BLAST_db](https://github.com/jgmv/local_BLAST_db)
