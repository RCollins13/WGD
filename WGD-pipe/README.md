
This repository includes scripts to apply the 6-F correction to bincov matrix.  To apply the correction, a 4 column tab-deliminated file including the following information is should be prepared for each sample:
```
1. name of chromosome 
2. full path to raw binCov file
3. desired full path to output corrected binCov file
4. full path to 6F metadata matrix for that chromosome
```

File should be named as `6F_metadata_ref/RD_multiCorrection.{sample}.txt`, here's an example:
```
1	bincov_per_chr/CMC-HBCC-DNA-ACC-5669.1.bed.gz	bincov_corrected/CMC-HBCC-DNA-ACC-5669.1.bed.gz	6F_metadata_matrix/hg37/6F_correction_matrix.1.100bp_matrix.bed.gz
2	bincov_per_chr/CMC-HBCC-DNA-ACC-5669.2.bed.gz	bincov_corrected/CMC-HBCC-DNA-ACC-5669.2.bed.gz	6F_metadata_matrix/hg37/6F_correction_matrix.2.100bp_matrix.bed.gz
3	bincov_per_chr/CMC-HBCC-DNA-ACC-5669.3.bed.gz	bincov_corrected/CMC-HBCC-DNA-ACC-5669.3.bed.gz	6F_metadata_matrix/hg37/6F_correction_matrix.3.100bp_matrix.bed.gz
4	bincov_per_chr/CMC-HBCC-DNA-ACC-5669.4.bed.gz	bincov_corrected/CMC-HBCC-DNA-ACC-5669.4.bed.gz	6F_metadata_matrix/hg37/6F_correction_matrix.4.100bp_matrix.bed.gz
5	bincov_per_chr/CMC-HBCC-DNA-ACC-5669.5.bed.gz	bincov_corrected/CMC-HBCC-DNA-ACC-5669.5.bed.gz	6F_metadata_matrix/hg37/6F_correction_matrix.5.100bp_matrix.bed.gz
6	bincov_per_chr/CMC-HBCC-DNA-ACC-5669.6.bed.gz	bincov_corrected/CMC-HBCC-DNA-ACC-5669.6.bed.gz	6F_metadata_matrix/hg37/6F_correction_matrix.6.100bp_matrix.bed.gz
7	bincov_per_chr/CMC-HBCC-DNA-ACC-5669.7.bed.gz	bincov_corrected/CMC-HBCC-DNA-ACC-5669.7.bed.gz	6F_metadata_matrix/hg37/6F_correction_matrix.7.100bp_matrix.bed.gz
8	bincov_per_chr/CMC-HBCC-DNA-ACC-5669.8.bed.gz	bincov_corrected/CMC-HBCC-DNA-ACC-5669.8.bed.gz	6F_metadata_matrix/hg37/6F_correction_matrix.8.100bp_matrix.bed.gz
9	bincov_per_chr/CMC-HBCC-DNA-ACC-5669.9.bed.gz	bincov_corrected/CMC-HBCC-DNA-ACC-5669.9.bed.gz	6F_metadata_matrix/hg37/6F_correction_matrix.9.100bp_matrix.bed.gz
10	bincov_per_chr/CMC-HBCC-DNA-ACC-5669.10.bed.gz	bincov_corrected/CMC-HBCC-DNA-ACC-5669.10.bed.gz	6F_metadata_matrix/hg37/6F_correction_matrix.10.100bp_matrix.bed.gz
11	bincov_per_chr/CMC-HBCC-DNA-ACC-5669.11.bed.gz	bincov_corrected/CMC-HBCC-DNA-ACC-5669.11.bed.gz	6F_metadata_matrix/hg37/6F_correction_matrix.11.100bp_matrix.bed.gz
12	bincov_per_chr/CMC-HBCC-DNA-ACC-5669.12.bed.gz	bincov_corrected/CMC-HBCC-DNA-ACC-5669.12.bed.gz	6F_metadata_matrix/hg37/6F_correction_matrix.12.100bp_matrix.bed.gz
13	bincov_per_chr/CMC-HBCC-DNA-ACC-5669.13.bed.gz	bincov_corrected/CMC-HBCC-DNA-ACC-5669.13.bed.gz	6F_metadata_matrix/hg37/6F_correction_matrix.13.100bp_matrix.bed.gz
14	bincov_per_chr/CMC-HBCC-DNA-ACC-5669.14.bed.gz	bincov_corrected/CMC-HBCC-DNA-ACC-5669.14.bed.gz	6F_metadata_matrix/hg37/6F_correction_matrix.14.100bp_matrix.bed.gz
15	bincov_per_chr/CMC-HBCC-DNA-ACC-5669.15.bed.gz	bincov_corrected/CMC-HBCC-DNA-ACC-5669.15.bed.gz	6F_metadata_matrix/hg37/6F_correction_matrix.15.100bp_matrix.bed.gz
16	bincov_per_chr/CMC-HBCC-DNA-ACC-5669.16.bed.gz	bincov_corrected/CMC-HBCC-DNA-ACC-5669.16.bed.gz	6F_metadata_matrix/hg37/6F_correction_matrix.16.100bp_matrix.bed.gz
17	bincov_per_chr/CMC-HBCC-DNA-ACC-5669.17.bed.gz	bincov_corrected/CMC-HBCC-DNA-ACC-5669.17.bed.gz	6F_metadata_matrix/hg37/6F_correction_matrix.17.100bp_matrix.bed.gz
18	bincov_per_chr/CMC-HBCC-DNA-ACC-5669.18.bed.gz	bincov_corrected/CMC-HBCC-DNA-ACC-5669.18.bed.gz	6F_metadata_matrix/hg37/6F_correction_matrix.18.100bp_matrix.bed.gz
19	bincov_per_chr/CMC-HBCC-DNA-ACC-5669.19.bed.gz	bincov_corrected/CMC-HBCC-DNA-ACC-5669.19.bed.gz	6F_metadata_matrix/hg37/6F_correction_matrix.19.100bp_matrix.bed.gz
20	bincov_per_chr/CMC-HBCC-DNA-ACC-5669.20.bed.gz	bincov_corrected/CMC-HBCC-DNA-ACC-5669.20.bed.gz	6F_metadata_matrix/hg37/6F_correction_matrix.20.100bp_matrix.bed.gz
21	bincov_per_chr/CMC-HBCC-DNA-ACC-5669.21.bed.gz	bincov_corrected/CMC-HBCC-DNA-ACC-5669.21.bed.gz	6F_metadata_matrix/hg37/6F_correction_matrix.21.100bp_matrix.bed.gz
22	bincov_per_chr/CMC-HBCC-DNA-ACC-5669.22.bed.gz	bincov_corrected/CMC-HBCC-DNA-ACC-5669.22.bed.gz	6F_metadata_matrix/hg37/6F_correction_matrix.22.100bp_matrix.bed.gz
X	bincov_per_chr/CMC-HBCC-DNA-ACC-5669.X.bed.gz	bincov_corrected/CMC-HBCC-DNA-ACC-5669.X.bed.gz	6F_metadata_matrix/hg37/6F_correction_matrix.X.100bp_matrix.bed.gz
Y	bincov_per_chr/CMC-HBCC-DNA-ACC-5669.Y.bed.gz	bincov_corrected/CMC-HBCC-DNA-ACC-5669.Y.bed.gz	6F_metadata_matrix/hg37/6F_correction_matrix.Y.100bp_matrix.bed.gz
```


To run the correction through snakemake:
```
snakemake	
```

