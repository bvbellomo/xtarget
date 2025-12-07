# import cProfile
from dataclasses import dataclass
import math
import typing
from typing import Self

import numpy as np
import numpy.typing as npt


# constants
PATH = 'C:\\bvb\\school\\'
MCGEARY_LIN_FILENAME = 'GSE140217_HeLa_transfection_logtpm.txt'
SCAN_FILENAME = 'scan.csv'
TPSITE_FILENAME = 'tpsites.csv'
MCGEARY_BATCH_COUNT = 5
MCGEARY_REPLICATION_COUNT = 2
MAX_WEIGHT_ITERATIONS = 10000
HALT_GRADIENT = 0.001  # 0.001 works well
BASELINE_MAX_ITERATIONS = 100
BASELINE_CONVERGENCE = 0.1  # 0.01
LEARNING_RATE = 0.01  # 0.001 recommended by Adam paper;0.013 seems good, 0.014 oscillates. 0.01 works well
STRAND_COEF_ITERATIONS = 10
STRAND_CONVERGENCE = 0.0001
MAX_MIR_LEN = 23

FEATURES: dict[int, str] = {1: 'noncan', 2: '6mer', 3: 'a1', 4: 'm8', 5: 'full'}
REGIONS: dict[int, str] = {1: "5'UTR", 2: "Coding", 3: "3'UTR"}
SCAN_COL_COUNT = len(FEATURES) * len(REGIONS)

EXCLUDE_MRNA = {'NM_001037293', 'NM_001136562'}  # Discrepencies between what MANE data considers the same gene and McGeary considered the same gene make these difficult to work with

MAX_BASELINE_SCALE = 1.02


# classes and functions
class RegressionException(Exception):
    """Raise for general exceptions in regression code"""


@dataclass
class BatchMaxValueItem:
    v: np.float64
    values: list[int]

    @classmethod
    def empty(cls) -> Self:
        return cls(np.float64(0), [])


def sigmoid(x: typing.Any, clip_min=-500, clip_max=500) -> typing.Any:
    return 1 / (1 + np.exp(-np.clip(x, clip_min, clip_max)))


def squared_sign(x):
    return x * x if x > 0 else -x * x


def sqrt_sign(x):
    return math.sqrt(x) if x > 0 else 0 - math.sqrt(0 - x)


class BaselineData:
    def __init__(self):
        self._baseline_by_id = {}
        self._max_initial_values_by_id = {}
        self._baseline_values = []
        self.baseline_by_row: npt.NDArray[np.float64] | None = None

    def load(self, accession: str, batch_values: list[BatchMaxValueItem]) -> None:
        """_summary_

        Args:
            id (str): _description_
            batch_values (typing.List[tuple[float, typing.List[int]]]): _description_
        """
        new_baseline_array = []
        max_mrna_value = 0
        for _, v in enumerate(batch_values):
            val_index = len(self._baseline_values)
            self._baseline_values.append(v.v)
            max_mrna_value = max(max_mrna_value, v.v)
            new_baseline_array.append((val_index, v.values))
        self._baseline_by_id[accession] = new_baseline_array
        self._max_initial_values_by_id[accession] = max_mrna_value

    def append_row_no(self, accession: str, batch: int, row: int):
        self._baseline_by_id[accession][batch][1].append(row)

    def update_row_data(self):
        max_row_index = 0
        rows_to_delete = {}
        for k, v in self._baseline_by_id.items():
            for bi, bv in enumerate(v):
                if bv[1]:  # some rows may not have loaded successfully
                    max_row_index = max(max_row_index, *bv[1])
                elif k in rows_to_delete:
                    rows_to_delete[k].append(bi)
                else:
                    rows_to_delete[k] = [bi]
        for dk, dv in rows_to_delete.items():
            b = self._baseline_by_id[dk]
            if len(dv) == len(b):
                del self._baseline_by_id[dk]
            else:
                raise RegressionException("Removing some, but not all, rows for an mRNA can be implemented, but should not be needed.")
        self.baseline_by_row = np.zeros(max_row_index + 1, np.float64)
        for k, v in self._baseline_by_id.items():
            for bi, bv in enumerate(v):
                for _, rv in enumerate(bv[1]):
                    self.baseline_by_row[rv] = self._baseline_values[bv[0]]

    def adjust_baseline(self, X, y, weights) -> float:
        for row in self._baseline_by_id:
            if row in EXCLUDE_MRNA:
                print('error excluding mRNA')

        max_change = 0
        repression_all = np.dot(X, weights)
        # for mrna_val in self._baseline_by_id.values():
        for _, (mrna_key, mrna_val) in enumerate(self._baseline_by_id.items()):
            for batch_index, batch_val in enumerate(mrna_val):
                new_baseline_accumulation = 0
                new_baseline_count = len(batch_val[1])
                for i in batch_val[1]:
                    repression = repression_all[i]
                    max_baseline_value = self._max_initial_values_by_id[mrna_key] * MAX_BASELINE_SCALE
                    if repression < 1:
                        a = y[i] / (1 - repression)
                        if a < max_baseline_value:
                            new_baseline_accumulation += a
                        else:
                            new_baseline_accumulation += max_baseline_value
                    else:
                        new_baseline_accumulation += max_baseline_value

                new_baseline_value = new_baseline_accumulation / new_baseline_count
                old = self._baseline_values[batch_val[0]]
                max_change = max(abs(old - new_baseline_value), max_change)
                self._baseline_values[batch_val[0]] = new_baseline_value
                if self.baseline_by_row is not None:
                    self.baseline_by_row[batch_val[1]] = new_baseline_value
                else:
                    raise Exception('bad baseline_by_row')
        return max_change

    def set_ids_and_values(self, baseline_by_id, baseline_values) -> None:
        self._baseline_by_id = baseline_by_id
        self._baseline_values = baseline_values

    def copy(self) -> Self:
        new_data = type(self)()
        if self.baseline_by_row:
            new_data.baseline_by_row = np.copy(self.baseline_by_row)
        else:
            raise Exception('bad baseline_by_row')
        new_data.set_ids_and_values(self._baseline_by_id, self._baseline_values.copy())
        return new_data

    def print_exclusion_errors(self) -> None:
        for row in self._baseline_by_id:
            if row in EXCLUDE_MRNA:
                print('error excluding mRNA')


def get_col(canonical: int, bound_state: bool) -> int:
    i = canonical - 1
    if bound_state:
        i += len(FEATURES)
    return i


def r_squared(y, yhat, col_count: int) -> tuple[float, float, int]:
    rss = sum((y-yhat)**2)
    tss = sum((y-np.mean(y))**2)
    result = 1 - (float(rss))/tss
    samples = len(y)
    adjusted_r_squared = 1 - (1 - result) * (samples - 1) / (samples - col_count - 1)
    return adjusted_r_squared, result, rss


def check_duplicates(scan, exclude):
    map_to_mane = {}
    reverse_map = {}
    mrna = set()
    for row in scan:
        id_version = str(row[1])
        accession = id_version.split('.')[0]
        original = str(row[0])
        if original not in exclude:
            map_to_mane[original] = id_version
            if accession in reverse_map and reverse_map[accession] != original:
                print('duplicate: both ' + reverse_map[accession] + ' and ' + original + ' map to ' + accession)
            else:
                reverse_map[accession] = original
                mrna.add(accession)


def load_data() -> tuple[BaselineData,
                         typing.List[tuple[str, str, int, int]],    # labels
                         npt.NDArray[np.float64],                   # X
                         npt.NDArray[np.float64],                   # X_strand
                         typing.Dict[str, tuple[float, float]],     # strand_coefficients
                         npt.NDArray[np.float64]]:                  # y

    with open(PATH + SCAN_FILENAME, 'r', encoding='utf-8') as file_in1:
        scan = np.genfromtxt(file_in1, skip_header=0, dtype=None, delimiter=',', encoding='utf-8')

    check_duplicates(scan, EXCLUDE_MRNA)

    with open(PATH + MCGEARY_LIN_FILENAME, 'r', encoding='utf-8') as file_in2:
        cols = file_in2.readline().split()
        mcgeary_lin: np.ndarray = np.genfromtxt(file_in2, skip_header=0, dtype=None, delimiter='\t', encoding='utf-8')

    # initialize strand coefficients
    strand_coefficients = {}
    for col in cols:
        mir = col.split('_')[0]
        if mir not in strand_coefficients:
            strand_coefficients[mir] = (0.5, 0.5)

    cols.insert(0, 'id')

    # initialize mrna expression for baseline and after transfection
    transfected = {}
    baseline_data = BaselineData()
    for row in mcgeary_lin:
        accession = str(row[0]).split('.')[0]
        if accession not in EXCLUDE_MRNA:
            max_v_batch = [BatchMaxValueItem.empty() for _ in range(MCGEARY_BATCH_COUNT)]
            for i in range(1, len(cols)):
                v: np.float64 = 2**row[i]  # transform values from log(2) to linear
                col_ary = cols[i].split('_')
                batch = int(col_ary[2].replace('batch', ''))
                max_v_batch[batch - 1].v = max(max_v_batch[batch - 1].v, v)
                transfected[accession + '_' + cols[i]] = v
            baseline_data.load(accession, max_v_batch)

    baseline_data.print_exclusion_errors()

    X_strand = []
    labels = []
    y = []
    x_strand_row = []
    row_number = 0
    for row in scan:
        original = str(row[0])
        if original not in EXCLUDE_MRNA:
            mir = str(row[2])
            strand = str(row[3])
            if original == 'NM_000031' and mir == 'mir143' and strand == '5p':
                pass
            # This assumes strands in the scan are always in order from 5p to 3p
            if strand == '5p':
                x_strand_row = []
                for i in range(SCAN_COL_COUNT):
                    val = int(row[4 + i])
                    x_strand_row.append(val)
            elif strand == '3p':
                for i in range(SCAN_COL_COUNT):
                    val = int(row[4 + i])
                    x_strand_row.append(val)
                for bi in range(1, MCGEARY_BATCH_COUNT + 1):
                    for ri in (range(1, MCGEARY_REPLICATION_COUNT + 1)):
                        key = original + '_' + mir + '_rep' + str(ri) + '_batch' + str(bi)
                        if key in transfected:
                            labels.append((original, mir, ri, bi))
                            y.append(transfected[key])
                            X_strand.append(x_strand_row)
                            baseline_data.append_row_no(original, bi-1, row_number)
                            row_number += 1
    baseline_data.update_row_data()
    X = np.zeros((len(y), SCAN_COL_COUNT))
    X_strand = np.asarray(X_strand)
    y = np.asarray(y)

    update_x_data(labels, X, X_strand, strand_coefficients)

    return baseline_data, labels, X, X_strand, strand_coefficients, y


def update_x_data(labels, X, X_strand, strand_coefficients):
    for li, lv in enumerate(labels):
        s = strand_coefficients[lv[1]]
        X[li] = s[0] * X_strand[li, :SCAN_COL_COUNT] + s[1] * X_strand[li, SCAN_COL_COUNT:]


def adjust_weights_gd(baseline_by_row: npt.NDArray[np.float64], weights, X, y):
    # debug_change_history = []
    # debug_weight_history = []
    oscillation = 0
    last_dw = []
    halt = HALT_GRADIENT * LEARNING_RATE
    t = 0
    for t in range(1, MAX_WEIGHT_ITERATIONS + 1):
        rep = np.clip(np.dot(X, weights), 0, 1)
        pred = np.multiply(baseline_by_row, 1 - rep)
        # mask = np.where(rep > 1, 0, 1)
        # delta = np.multiply(2 * (pred - y), mask)
        delta = 2 * (pred - y)
        dw = 0 - np.dot(X.T, delta) / len(y)
        old = weights
        new_weights = weights - dw * LEARNING_RATE

        idx = np.flatnonzero(new_weights < 0)
        if len(idx) > 0:
            X2 = np.delete(X, idx, axis=1)
            weights2 = np.delete(weights, idx)
            weights2 = adjust_weights_gd(baseline_by_row, weights2, X2, y)
            return replace_zero_vals(idx, weights2)
        weights = np.clip(new_weights, 0, 1)
        change = old - weights

        if np.linalg.norm(change) <= halt:
            return weights
        if t > 1:
            if np.dot(last_dw, dw) < 0:
                oscillation += 1
        last_dw = dw

    print('Warning, MAX_WEIGHT_ITERATIONS reached with Gradient Descent, t=' + str(t) + ', oscillation=' + str(oscillation))
    return weights


def remove_zero_vals(X):
    # remove any column that is all zeros
    idx = np.argwhere(np.all(X[..., :] == 0, axis=0))
    if len(idx) > 0:
        return idx, np.delete(X, idx, axis=1)
    else:
        return idx, X


def replace_zero_vals(idx, w):
    if len(idx) == 0:
        return w

    w2 = np.zeros(len(idx) + len(w), np.float64)
    j = 0
    for i in range(len(w2)):
        if i not in idx:
            w2[i] = w[j]
            j += 1
    return w2


def adjust_weights_analytic(baseline_by_row: npt.NDArray[np.float64], X, y):
    x2 = (X * baseline_by_row[:, None])
    xt = x2.T
    xtx = np.matmul(xt, x2)
    det = np.linalg.det(xtx)
    if det != 0:
        i = np.linalg.inv(xtx)
    else:
        i = np.linalg.pinv(xtx)
    y2 = baseline_by_row - y
    xty = np.matmul(xt, y2)
    w = np.matmul(i, xty)
    return np.clip(w, 0, 1)


def adjust_weights(baseline_by_row: npt.NDArray[np.float64], X, y):
    weights = adjust_weights_analytic(baseline_by_row, X, y)
    weights = adjust_weights_gd(baseline_by_row, weights, X, y)
    return weights


def adjust_baselines_and_weights(baseline_data: BaselineData, X, y):
    if baseline_data.baseline_by_row is None:
        raise Exception('bad baseline_data.baseline_by_row')
    weights = adjust_weights(baseline_data.baseline_by_row, X, y)
    for i in range(BASELINE_MAX_ITERATIONS):
        change = baseline_data.adjust_baseline(X, y, weights)
        if change < BASELINE_CONVERGENCE:
            return weights
        weights = adjust_weights(baseline_data.baseline_by_row, X, y)
    print('Warning, BASELINE_MAX_ITERATIONS reached')
    return weights


def adjust_strands_baselines_and_weights(baseline_data: BaselineData, X, X_strand, strand_coefficients, y, labels):
    weights = None
    for _ in range(STRAND_COEF_ITERATIONS):
        change = adjust_strand_coef(baseline_data.baseline_by_row, labels, X, X_strand, strand_coefficients, y)
        # weights and baselines are much faster and more important than strand coefficients
        # train them last so they have a good fit to our final strand coefficients
        weights = adjust_baselines_and_weights(baseline_data, X, y)
        if change < STRAND_CONVERGENCE:
            return weights
    print('Warning, STRAND_COEF_ITERATIONS reached')
    return weights


def adjusted_strand_ratio(fp: np.float64, tp: np.float64) -> tuple[np.float64, np.float64]:
    if tp > 0 and fp > 0:
        if fp > tp:
            r = tp / fp
            return (2 - r) / 2, r / 2
        else:
            r = fp / tp
            return r / 2, (2 - r) / 2
    elif tp <= 0 and fp <= 0:
        return np.float64(0), np.float64(0)
    elif fp > 0:
        return np.float64(1), np.float64(0)
    elif tp > 0:
        return np.float64(0), np.float64(1)
    raise Exception('bad strand ratio')


def adjust_strand_coef(baseline_by_row, labels, X, X_strand, strand_coefficients, y) -> float:
    max_change = 0
    for mir_name in strand_coefficients.keys():  # for each mirna
        # build a strand-seperated X feature array filtered for 1 miRNA only
        X2 = []
        y2 = []
        baseline_by_row_by_strand = []
        for label_index, label_val in enumerate(labels):  # for each row
            if label_val[1] == mir_name:
                X2.append(X_strand[label_index])
                y2.append(y[label_index])
                baseline_by_row_by_strand.append(baseline_by_row[label_index])
        X2 = np.asarray(X2)
        y2 = np.asarray(y2)
        baseline_by_row_by_strand = np.asarray(baseline_by_row_by_strand)

        # calculate strand coefficients
        weights2 = adjust_weights(baseline_by_row_by_strand, X2, y2)
        tm_power = np.float64(0)
        fm_power = np.float64(0)
        for row_i in range(0, len(y2)):
            b = baseline_by_row_by_strand[row_i]
            if b > 0:
                for col_i in range(0, SCAN_COL_COUNT):
                    fm_power += weights2[col_i] * X2[row_i, col_i] * b
                for col_i in range(SCAN_COL_COUNT, SCAN_COL_COUNT * 2):
                    tm_power += weights2[col_i] * X2[row_i, col_i] * b
        fm_power, tm_power = adjusted_strand_ratio(fm_power, tm_power)
        max_change = max(max_change, abs(fm_power - strand_coefficients[mir_name][0]))  # fm_power + tm_power will always be 1, it doesn't matter which we use
        strand_coefficients[mir_name] = (fm_power, tm_power)

    update_x_data(labels, X, X_strand, strand_coefficients)
    return max_change


def print_strand_coef(strand_coefficients):
    print('strand coefficients')
    for k, v in strand_coefficients.items():
        print(f'{k}\t5p = {v[0]:.2f}\t3p = {v[1]:.2f}')


def predict(baseline_data: BaselineData, X, weights):
    if baseline_data.baseline_by_row is None:
        raise Exception('bad baseline')
    return np.multiply(baseline_data.baseline_by_row, 1 - np.clip(np.dot(X, weights), 0, 1))

# def print_weights2(weights):
#     weight_labels = [''] * (len(REGIONS) * len(FEATURES))
#     for region_id, region_name in REGIONS.items():
#         for feature_id, feature_name in FEATURES.items():
#             weight_labels[region_id + feature_id] = region_name + '_' + feature_name

#     weight_labels = ["5'UTR 6mer", "5'UTR a1", "5'UTR m8", "5'UTR full",
#                      "Coding 6mer", "Coding a1", "Coding m8", "Coding full"]
#     if len(weights) == 12:
#         weight_labels += ["3'UTR 6mer", "3'UTR a1", "3'UTR m8", "3'UTR full"]
#     elif len(weights) == 16:
#         weight_labels += ["3'UTR 6mer N", "3'UTR a1 N", "3'UTR m8 N", "3'UTR full N",
#                           "3'UTR 6mer Y", "3'UTR a1 Y", "3'UTR m8 Y", "3'UTR full Y"]
#     else:
#         raise ValueError('could not label unknown number of weights')

#     for i, w in enumerate(weights):
#         print(f'{weight_labels[i]}\t{w:4f}')


def format_feature(val, precision) -> str:
    return f'\t{val:.{precision}f}' if val else '\t0'


def format_feature_region_matrix(data, precision) -> str:
    region_labels = ["5'UTR", "Coding"]
    if len(data) == SCAN_COL_COUNT:
        region_labels += ["3'UTR"]
    elif len(data) == len(FEATURES) * (len(REGIONS) + 1):
        region_labels += ["3'UTR N", "3'UTR Y"]
    else:
        return '\t'.join(f'{w:4f}' for w in data)

#    class_labels = ["6mer", "a1", "m8", "full"]
    # class_labels = ["6mer", "6mer+s", "a1", "a1+s", "m8", "m8+s", "full", "full+s"]

    s = '\t' + '\t'.join(FEATURES.values()) + '\n'
    for i, l in enumerate(region_labels):
        s += l

        for j in range(len(FEATURES)):
            s += format_feature(data[i * len(FEATURES) + j], precision)
        s += '\n'
    return s


def print_status(text: str, y, y_pred, col_count):
    ar2, r2, rss = r_squared(y, y_pred, col_count)
    print(f'{text}\t1-adjusted r2 = {1 - ar2:.5f}\t1-r squared = {1- r2:.5f}\trss = {rss:.0f}')


def main():
    baseline_data, labels, X, X_strand, strand_coefficients, y = load_data()
    print(format_feature_region_matrix(np.sum(X, 0), 0))

    weights = np.tile(.01, SCAN_COL_COUNT)
    weights = adjust_baselines_and_weights(baseline_data, X, y)
    y_pred = predict(baseline_data, X, weights)
    print_status('After training initial weights and baselines: ', y, y_pred, SCAN_COL_COUNT)
    print(format_feature_region_matrix(weights, 4))

    weights = adjust_strands_baselines_and_weights(baseline_data, X, X_strand, strand_coefficients, y, labels)
    print(format_feature_region_matrix(np.sum(X, 0), 0))
    y_pred = predict(baseline_data, X, weights)
    print_status('After training strand coefficients: ', y, y_pred, SCAN_COL_COUNT)
    print()
    print_strand_coef(strand_coefficients)
    print()
    print(format_feature_region_matrix(weights, 4))
    print('done')


# profiler = cProfile.Profile()
# profiler.enable()
main()
# profiler.disable()
# profiler.dump_stats(PATH + 'regression.stats')
