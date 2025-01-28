import os
import requests
from pathlib import Path


def read_cds(id):
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=fasta_cds_na&id=' + id
    response = requests.get(url)
    lines = response.text.splitlines()
    for line in lines:
        if line.startswith('>'):
            data = line.split('|')[1]
            id = data[0:data.find('_cds')]
            loc = '[location='
            index = data.find(loc) + len(loc)
            sloc = data[index:data.find(']', index)]
            ary = sloc.split('..')
            if ary[0] == '<1':
                print('partial open reading frame for ' + id)
            elif 'join' in ary[0]:
                print('multiple open reading frames for ' + id)
            else:
                return (id, int(ary[0]), int(ary[1]))


path = 'C:\\bvb\\school\\'
cds_path = path + 'cds\\'
fasta_path = path + 'fasta\\'
files = os.listdir(cds_path)
cds = {}
max_filename = 0
for f in files:
    max_filename = max(max_filename, int(Path(f).stem))
    with open(cds_path + f) as fi:
        line = fi.readline()
        while line:
            ary = line.split(',')
            cds[ary[0]] = (ary[1], ary[2])
            line = fi.readline()

with open(cds_path + str(max_filename + 1) + '.csv', 'w') as fo:
    files = os.listdir(fasta_path)
    for f in files:
        id = Path(f).stem
        if id.startswith('NM_032790'):
            pass
        if id not in cds:
            cd = read_cds(id)
            if cd is not None:
                fo.write(cd[0] + ',' + str(cd[1]) + ',' + str(cd[2]) + '\n')
                fo.flush()

print('done')
