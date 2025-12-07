# import cProfile
import csv
import os
import re
import numpy as np
from numba import njit
from collections import OrderedDict
from datetime import datetime

import bioutils


# constants
PATH = 'C:\\bvb\\school\\'
MCGEARY_LIN_FILENAME = 'GSE140217_HeLa_transfection_logtpm.txt'
BIOMART_REFSEQ_MAP_FILENAME = 'biomart_refseq_map.txt'
MANE_FILENAME = 'MANE.GRCh38.v1.3.summary.txt'
MIR_FILENAME = 'mir.csv'
THRESHOLD = 50

# 1 = noncanonical NW match
# 2 = 6mer = position 2–7 match seed match
# 3 = a1 = position 2–7 match + A at position 1
# 4 = m8 = position 2–8 match
# 5 = full = 2–8 with an A opposite position
FEATURES: dict[int, str] = {1: 'noncanonical', 2: '6mer', 3: 'a1', 4: 'm8', 5: 'full'}
FEATURE_COUNT: int = len(FEATURES)
FEATURE_TO_ID = {name: id for id, name in FEATURES.items()}

REGION_5UTR: int = 1
REGION_CODING: int = 3
REGION_3UTR: int = 5
REGION_COUNT: int = 3
PRINT_ROW_COUNT = 10

# module-level variables
ref_to_mane_map = {}
mane_ref_ids = set()
mir_strands: dict[str, tuple[str, str, str, str]] = {}
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
    fasta_result = bioutils.readFasta(fasta_path)
    if not fasta_result:
        print(mapped_id + '\tBad fasta file for ' + id_version)
        return False

    _, seq = fasta_result

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
        label = [accession, id_version, mirna, '5p']
        row = scan_mrna_mirna(seq, mir_data[1], coding_start, coding_end)
        scan_writer.writerow(label + row)
        label = [accession, id_version, mirna, '3p']
        row = scan_mrna_mirna(seq, mir_data[3], coding_start, coding_end)
        scan_writer.writerow(label + row)
    return True


def get_seed(mirna) -> str:
    seed = ''
    for _, c in enumerate(mirna):
        seed += bioutils.dna_complement(c)
    return seed[::-1]


@njit
def is_6mer(mirna: str, mrna: str, start: int, pos: int) -> bool:
    mrna_len = pos - start
    for i in range(mrna_len - 7, mrna_len - 1):
        if mrna[start + i] != mirna[i]:
            return False
    return True


@njit
def classify_site(seed: str, mrna: str, pos: int, matrix):
    start = pos - len(seed) + 1
    if is_6mer(seed, mrna, start, pos):
        a = int(mrna[pos - 1] == 'A')
        n1 = int(mrna[pos - 8] == seed[-8])
        site_classification = 2 + a + n1 * 2
        return site_classification
    else:
        score = align(seed, mrna[pos - len(seed):pos], matrix)
        if score >= THRESHOLD:
            return 1
        return 0


@njit
def get_col(region_id: int, feature_id: int) -> int:
    match region_id:
        case 1:
            i = 0
        case 3:
            i = FEATURE_COUNT
        case 5:
            i = FEATURE_COUNT * 2
        case _:
            raise Exception('Invalid region')
    i += feature_id - 1
    return i


@njit
def align(mirna: str, mrna: str, matrix) -> int:
    score_match_cg = 4
    score_match_at = 3
    score_mismatch = -2
    score_gap = -2

    if len(mirna) != len(mrna):
        raise Exception("unsupported")

    for i, c1 in enumerate(mrna):
        for j, c2 in enumerate(mirna):
            diag = matrix[i][j]
            if c1 == c2:
                if c1 == "C" or c1 == "G":
                    diag += score_match_cg
                else:
                    diag += score_match_at
                matrix[i + 1][j + 1] = diag
            else:
                diag += score_mismatch
                left = matrix[i + 1][j] + score_gap
                up = matrix[i][j + 1] + score_gap
                matrix[i + 1][j + 1] = max(diag, up, left)
    
    return matrix[len(mirna)][len(mrna)]


def scan_mrna_mirna(mrna: str, seed: str, coding_start, coding_end) -> list:
    data = [0 for i in range(REGION_COUNT * FEATURE_COUNT)]
    matrix = np.zeros((len(mrna) + 1, len(mrna) + 1), dtype=np.int8)

    pos = len(seed)
    while pos < len(mrna):
        feature_id = classify_site(seed, mrna, pos, matrix)

        if feature_id >= 1:
            if pos < coding_start:
                region = REGION_5UTR
            elif pos - len(seed) + 1 > coding_end:
                region = REGION_3UTR
            else:
                region = REGION_CODING
            data[get_col(region, feature_id)] += 1
        pos += 1
    return data


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
            header_left = ['accession', 'id_version', 'mirna', 'strand']
            header = header_left.copy()
            header.extend([''] * REGION_COUNT * FEATURE_COUNT)
            REGIONS: dict[int, str] = {REGION_5UTR: "5'UTR", REGION_CODING: "Coding", REGION_3UTR: "3'UTR"}
            for region_id, region_name in REGIONS.items():
                for feature_id, feature_name in FEATURES.items():
                    header[4 + get_col(region_id, feature_id)] = region_name + '_' + feature_name
            scan_writer.writerow(header)
            while line:
                accession = line.split()[0].split('.')[0]
                if scan_mrna(scan_writer, accession, mirnas):
                    success_count += 1
                row_count += 1
                line = file_in.readline()
                if row_count % PRINT_ROW_COUNT == 0:
                    print(str(row_count) + ' rows done')
    print('done: ' + str(success_count) + ' of ' + str(row_count) + ' succeeded in ' + str(datetime.now() - start_time))


# profiler = cProfile.Profile()
# profiler.enable()
main()
# profiler.disable()
# profiler.dump_stats('c:\\bvb\\school\\scan.stats')
