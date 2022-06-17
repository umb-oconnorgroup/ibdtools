#! /usr/bin/env python3

import numpy as np
from subprocess import run


def read_ibdtools_mat(ibdtools_mat_fn):
    """read ibd matrix and turn it to a numpy array"""

    # decompress the matrix file
    flat_mat_fn = "tmp_ibdtools_mat_flat"
    cmd = f"gzip -dc {ibdtools_mat_fn} > {flat_mat_fn}"
    run(cmd, shell=True)

    # uint8: is_hap,                         offset 0
    # uint32: d                              offset 1
    # uint64: subpop_array_sz_bytes (xx)     offset 5
    # xx bytes: subpop array                 offset 13
    # uint64: ibd array sz in bytes (yy)     offset 13 + xx
    # yy bytes: array                        offset 21 + xx

    # dimension
    is_hap = np.fromfile(flat_mat_fn, dtype=np.uint8, offset=0, count=1)[0]
    d = np.fromfile(flat_mat_fn, dtype=np.uint32, offset=1, count=1)[0]
    xx = np.fromfile(flat_mat_fn, dtype=np.uint64, offset=5, count=1)[0]
    yy = np.fromfile(flat_mat_fn, dtype=np.uint64, offset=np.uint64(13 + xx), count=1)[
        0
    ]

    arr_sz = yy // 2  # each array item is uint16
    subpop_sz = np.uint64(xx // 4)  # each supop id is uint32

    print(f"- is_hap      : {is_hap}")
    print(f"- arr_sz      : {arr_sz}")
    print(f"- d           : {d}")
    print(f"- #subpop ids : {subpop_sz}")

    # load subpop id
    if xx == 0:
        subpop_ids = np.array([])
    else:
        subpop_ids = np.fromfile(
            flat_mat_fn, dtype=np.uint32, offset=13, count=subpop_sz
        )

    # load the whole arr
    cnt = d * (d - 1) // 2
    arr = np.fromfile(
        flat_mat_fn, dtype=np.uint16, offset=np.uint64(21 + xx), count=cnt
    )

    # array to lower triangular matrix
    M = np.zeros(shape=(d, d), dtype=np.float32)
    tril_idx = np.tril_indices(d, k=-1)

    # In ibdtools, IBD total is coded as 16 bit unsigned integer after
    # multiplied by 10
    M[tril_idx] = arr / 10

    # Make symmetricall matrix
    M += M.T

    return is_hap, subpop_ids, M


res = read_ibdtools_mat("./gw.mat.mat")
print(res, res[2].sum())
