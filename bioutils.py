import os.path


def dna_complement(nucleotide):
    match nucleotide:
        case 'A':
            return 'T'
        case 'T' | 'U':
            return 'A'
        case 'C':
            return 'G'
        case 'G':
            return 'C'
    raise Exception('unknown nucleotide ' + nucleotide)


def rna_complement(nucleotide):
    match nucleotide:
        case 'A':
            return 'U'
        case 'T' | 'U':
            return 'A'
        case 'C':
            return 'G'
        case 'G':
            return 'C'
    raise Exception('unknown nucleotide ' + nucleotide)


def readFasta(path: str) -> tuple[str, str] | None:
    if not os.path.exists(path):
        return None

    with open(path) as f:
        line = f.readline()
        transcript_label_version = line.split()[0][1::]
        seq = ''
        while line:
            if not line.startswith('>'):
                seq = seq + line
            line = f.readline()
    return transcript_label_version, ''.join(seq.splitlines())
