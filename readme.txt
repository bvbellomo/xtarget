X-Target is a Python package is a simple, extensible Python module identify miRNA target sites and predict there effectiveness. X-Target was specifically designed to enable easy testing of potential new features.


Local Directories:

X-Target is configured to use a working directory of 'c:\bvb\school\' however this is easily changed.  X-Target requires 'fasta' and 'cds' directories for intermediate data inside the working directory.



Data Sources:

X-Target uses transfection data from:
McGeary, S. E., Lin, K. S., Shi, C. Y., Pham, T. M., Bisaria, N., Kelley, G. M., & Bartel, D. P. (2019). ‘The biochemical basis of microRNA target-ing efficacy’, Science, 366/6472: eaav1741. American Association for the Advancement of Science. DOI: 10.1126/science.aav1741

GSE140217_HeLa_transfection_logtpm.txt should be downloaded from NCBI Geo (https://www.ncbi.nlm.nih.gov/geo)


X-Target uses MANE v1.4 data from:
Morales, J., Pujar, S., Loveland, J. E., Astashyn, A., Bennett, R., Berry, A., Cox, E., et al. (2022). ‘A joint NCBI and EMBL-EBI transcript set for clinical genomics and research’, Nature, 604/7905: 310–5. DOI: 10.1038/s41586-022-04558-8

MANE.GRCh38.v1.4.summary.txt should be downloaded from NCBI (https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.4/)


X-Target uses Ensembl Biomart mapping data.

Create an unfiltered query that eliminates duplicates for 'Gene stable ID', 'Gene stable ID version', 'Transcript stable ID', 'Transcript stable ID version', 'HGNC ID', 'RefSeq mRNA ID' 'RefSeq mRNA predicted ID' and save it as biomart_refseq_map.txt.  This will be a large file.


X-Target requires a handed edited mir.csv created by X-Target's authors and included with the source code.  Sequences were taken from https://www.mirbase.org.



Python modules:

fasta_download.py

Downloads the latest fasta formatted data for MANE sequences from NCBI.  Requires MANE.GRCh38.v1.4.summary.txt


cds_get.py

Annotates 5'UTR, coding and 3'UTR based on data from NCBI.  Requires MANE.GRCh38.v1.4.summary.txt


bioutils.py

Module of biology related utility functions.


scan.py

Scans messenger RNAs for microRNA target sites.


regression.py

Regression analysis of scanned target sites.