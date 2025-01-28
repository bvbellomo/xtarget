import numpy as np
import os
import requests

path = 'C:\\bvb\\school\\'
mane_filename = 'MANE.GRCh38.v1.4.summary.txt'
efetch_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=fasta&id='

with open(path + mane_filename, 'r') as mane_file:
    mane_data = np.genfromtxt(mane_file, skip_header=1, dtype=str, delimiter='\t', usecols=(5))

for id_version in mane_data:
    fasta_filename = id_version + '.fasta'
    fasta_path = path + 'fasta\\' + fasta_filename
    if not os.path.exists(fasta_path):
        url = efetch_url + id_version
        response = requests.get(url)
        with open(fasta_path, 'w') as fasta_file:
            fasta_file.write(response.text)


print('done')
