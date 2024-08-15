import numpy as np

# Get cumulative sum vectors
def get_cumsum(sequence):
    y = np.cumsum(sequence)
    z = np.cumsum(np.square(sequence))
    return np.append([0], y), np.append([0], z)


# function to create loss value from 'start' to 'end' given cumulative sum vector y (data) and z (square)
def L(start, end, y, z):
    _y = y[end+1] - y[start]
    _z = z[end+1] - z[start]
    return _z - np.square(_y)/(end-start+1)


# function to get the list of changepoint from vector tau_star
def trace_back(tau_star):
    tau = tau_star[-1]
    chpnt = np.array([len(tau_star)], dtype=int)
    while tau > 0:
        chpnt = np.append(tau, chpnt)
        tau = tau_star[tau-1]
    return np.append(0, chpnt)


# solve opart
def opart(sequence, lda):
    y, z = get_cumsum(sequence)             # cumsum vector
    sequence_length = len(sequence)-1       # length of sequence

    # Set up
    C = np.zeros(sequence_length + 1)
    C[0] = -lda

    # Get tau_star
    tau_star = np.zeros(sequence_length+1, dtype=int)
    for t in range(1, sequence_length+1):
        V = C[:t] + lda + L(1 + np.arange(t), t, y, z)  # calculate set V
        last_chpnt = np.argmin(V)                       # get optimal tau from set V
        C[t] = V[last_chpnt]                            # update C_i
        tau_star[t] = last_chpnt                        # update tau_star

    set_of_chpnt = trace_back(tau_star[1:])             # get set of changepoints
    return set_of_chpnt[1:-1]