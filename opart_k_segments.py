import itertools
import numpy as np

def get_mean(sequence, chpnt):
    mean = np.zeros(len(sequence))
    chpnt = np.append(chpnt, len(sequence)-1)
    chpnt = np.append(-1, chpnt)
    chpnt = chpnt + 1
    for i in range(len(chpnt)-1):
        mean[chpnt[i]:chpnt[i+1]] = np.mean(sequence[chpnt[i]:chpnt[i+1]])
    return mean


def square_loss(seq1, seq2):
    return np.mean(np.square(seq1 - seq2))


def opart(sequence, k):
    if k > 1:
        combinations = list(itertools.combinations(np.arange(len(sequence)-1), k-1))
    else:
        return np.array([], dtype=np.int64)
    best_square_loss = np.inf
    for chpnt in combinations:
        mean = get_mean(sequence, chpnt)
        if square_loss(sequence, mean) < best_square_loss:
            best_square_loss = square_loss(sequence, mean)
            best_chpnt = chpnt
    return np.array(best_chpnt, dtype=np.int64)