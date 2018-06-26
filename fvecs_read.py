# Read a set of vectors stored in the fvec format (int + n * float)
# The function returns a set of output vector (one vector per column)
#
# Syntax:
#   v = fvecs_read (filename)     -> read all vectors
#   v = fvecs_read (filename, n)  -> read n vectors
#   v = fvecs_read (filename, [a b]) -> read the vectors from a to b (indices starts from 1)

import numpy as np

def fvecs_read(filename):
    fv = np.fromfile(filename, dtype=np.float32)
    if fv.size == 0:
        return np.zeros((0,0))
    dim = fv.view(np.int32)[0]
    assert dim > 0
    fv = fv.reshape(-1, 1 + dim)
    if not all(fv.view(np.int32)[:, 0] == dim):
        raise IOError("Non-uniform vector sizes in " + filename)
    fv = fv[:, 1:]

    return fv


print (fvecs_read("/Users/Stella/Documents/CSFYP/GQR/gqr/data/audio/audio_base.fvecs"))
