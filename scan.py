import csv
import os
import re
from collections import OrderedDict
from datetime import datetime

import numpy as np

import bioutils


# constants
PATH = 'C:\\bvb\\school\\'
MCGEARY_LIN_FILENAME = 'GSE140217_HeLa_transfection_logtpm.txt'
BIOMART_REFSEQ_MAP_FILENAME = 'biomart_refseq_map.txt'
MANE_FILENAME = 'MANE.GRCh38.v1.3.summary.txt'
MIR_FILENAME = 'mir.csv'
REGIONS = 3
FEATURES_PER_REGION = 4

# module-level variables
ref_to_mane_map = {}
mane_ref_ids = set()
mir_strands = {}
cds_data = {}


def get_fasta_filename(accession):
    dir = os.path.join(PATH, 'fasta')
    files = [file for file in os.listdir(dir) if file.startswith(accession) and file.endswith('.fasta')]
    valid_filenames = []
    for filename in files:
        ary = filename.split('.')
        if ary[0] == accession:
            valid_filenames.append((filename, int(ary[1])))
    sorted_filenames = sorted(valid_filenames, key=lambda x: 0 - x[1])
    if len(sorted_filenames) > 0:
        first_filename = sorted_filenames[0]
        return os.path.join(dir, first_filename[0]), first_filename[1]
    return None, None


def scan_mrna(scan_writer, accession: str, mirnas: list) -> bool:
    # Always use the latest MANE transcript
    if accession in ref_to_mane_map:
        mapped_id = ref_to_mane_map[accession]
    else:
        mapped_id = accession
    if mapped_id not in mane_ref_ids:
        print(mapped_id + '\tcould not be mapped to MANE')
        return False

    fasta_path, version = get_fasta_filename(mapped_id)
    if fasta_path is None:
        print(mapped_id + '\tno fasta file')
        return False

    id_version = mapped_id + '.' + str(version)

    if not os.path.exists(fasta_path):
        print(mapped_id + '\tNo fasta file for ' + id_version)
        return False
    _, seq = bioutils.readFasta(fasta_path)

    if id_version not in cds_data:
        print(mapped_id + '\tno cds')
        return False
    cds = cds_data[id_version]

    coding_start = cds[0] - 1
    coding_end = cds[1] - 1

    for mirna in mirnas:
        mir_data = mir_strands[mirna]
        if mir_data is None:
            print('No mir data for ' + mirna)
            return False
        # coding coordinates are 1 based, we want to pass 0 based coordinates
        scan_mrna_mirna(scan_writer, [accession, id_version, mirna, '5p'], seq, mir_data[1], len(mir_data[0]), coding_start, coding_end)
        scan_mrna_mirna(scan_writer, [accession, id_version, mirna, '3p'], seq, mir_data[3], len(mir_data[2]), coding_start, coding_end)
    return True


def get_seed(mirna):
    seed = list(mirna[0:len(mirna)-1])
    for i, c in enumerate(seed):
        seed[i] = bioutils.dna_complement(c)
    seed.reverse()
    return np.array(seed)


def is_canonical(mirna: list, mrna: list, start: int, pos: int):
    # 1 = 6mer = position 2–7 match seed match
    # 2 = a1 = position 2–7 match + A at position 1
    # 3 = m8 = position 2–8 match
    # 4 = full = 2–8 with an A opposite position
    mrna_len = pos - start
    for i in range(mrna_len - 7, mrna_len - 1):
        if mrna[start + i] != mirna[i]:
            return 0

    a = int(mrna[pos - 1] == 'A')
    n1 = int(mrna[pos - 8] == mirna[-8])
    site_classification = 1 + a + n1 * 2
    # g13 – g16 supplemental
    # supplemental = 0
    # for i in range(13, 17):
    #     if mrna[-i] == mirna[-i]:
    #         supplemental += 1
    # if supplemental >= 2:
    #     site_classification += 1
    return site_classification


def get_col(region: int, canonical: int) -> int:
    match region:
        case 1:
            i = 0
        case 3:
            i = FEATURES_PER_REGION
        case 5:
            i = FEATURES_PER_REGION * 2
    i += canonical - 1
    return i


def scan_mrna_mirna(scan_writer, label: list, mrna: list, seed: list, mir_len: int, coding_start, coding_end):
    data = [0 for i in range(REGIONS * FEATURES_PER_REGION)]

    pos = mir_len
    while pos < len(mrna):
        start = pos - mir_len + 1
        canonical = is_canonical(seed, mrna, start, pos)

        if canonical >= 1:
            if pos < coding_start:
                region = 1  # 5'UTR
            elif start > coding_end:
                region = 5  # 3'UTR
            else:
                region = 3  # Coding
            data[get_col(region, canonical)] += 1
        pos += 1
    scan_writer.writerow(label + data)


def main():
    start_time = datetime.now()

    # read mir strand data
    with open(PATH + MIR_FILENAME, encoding='utf-8') as file_in:
        reader = csv.reader(file_in, delimiter=',')
        next(reader, None)  # skip the headers
        for row in reader:
            mir_strands[row[1]] = (row[2], get_seed(row[2]), row[3], get_seed(row[3]))

    # read cds data
    cds_path = PATH + 'cds\\'
    files = os.listdir(cds_path)
    for file_in in files:
        with open(cds_path + file_in, encoding='utf-8') as file_in:
            line = file_in.readline()
            while line:
                ary = line.split(',')
                cds_data[ary[0]] = (int(ary[1]), int(ary[2]))
                line = file_in.readline()

    # read MANE data
    transcript_to_ref = {}
    hgnc_to_ref = {}
    with open(PATH + MANE_FILENAME, 'r', encoding='utf-8') as file_in:
        reader = csv.reader(file_in, delimiter='\t')
        next(reader, None)  # skip the headers
        for row in reader:
            hgnc = row[2]
            ref_seq_id = row[5].split('.')[0]
            transcript_stable_id = row[7].split('.')[0]
            transcript_to_ref[transcript_stable_id] = ref_seq_id
            hgnc_to_ref[hgnc] = ref_seq_id
            mane_ref_ids.add(ref_seq_id)

    # read biomart_refseq_map
    mane_ref_id = ''
    with open(PATH + BIOMART_REFSEQ_MAP_FILENAME, 'r', encoding='utf-8') as file_in:
        reader = csv.reader(file_in, delimiter='\t')
        next(reader, None)  # skip the headers
        for row in reader:
            transcript_stable_id = row[2]
            hgnc = row[4]
            ref_seq_id = row[5]
            ref_seq_pred_id = row[6]
            if hgnc in hgnc_to_ref:
                mane_ref_id = hgnc_to_ref[hgnc]
            elif transcript_stable_id in transcript_to_ref:
                mane_ref_id = transcript_to_ref[transcript_stable_id]
            else:
                mane_ref_id = None

            if mane_ref_id is not None:
                if ref_seq_id != '' and ref_seq_id != mane_ref_id:
                    ref_to_mane_map[ref_seq_id] = mane_ref_id
                if ref_seq_pred_id != '' and ref_seq_pred_id != mane_ref_id:
                    ref_to_mane_map[ref_seq_pred_id] = mane_ref_id

    with open(PATH + MCGEARY_LIN_FILENAME, 'r', encoding='utf-8') as file_in:
        header = file_in.readline()
        #   list mirna's from header removing duplicates but preserving order
        mirnas = list(OrderedDict.fromkeys(re.sub('_(rep|batch)\\d', '', header).split()))
        line = file_in.readline()
        row_count = 0
        success_count = 0
        with open(PATH + 'scan.csv', 'w', newline='', encoding='utf-8') as scan_file_out:
            scan_writer = csv.writer(scan_file_out)
            while line:
                accession = line.split()[0].split('.')[0]
                if scan_mrna(scan_writer, accession, mirnas):
                    success_count += 1
                row_count += 1
                line = file_in.readline()
    print('done: ' + str(success_count) + ' of ' + str(row_count) + ' succeeded in ' + str(datetime.now() - start_time))


main()
