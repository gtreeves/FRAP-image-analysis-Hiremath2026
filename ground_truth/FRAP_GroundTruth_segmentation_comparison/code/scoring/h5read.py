"""Read a variable from a MATLAB v7.3 (HDF5) .mat file, as saved by the MATLAB part of this package.
Arrays come back in MATLAB's index order (rows, columns, frames)."""
import h5py, numpy as np

def read(fn, name):
    with h5py.File(fn, 'r') as f:
        a = np.array(f[name][()])
    return a.T   # HDF5 stores MATLAB arrays column-major; transpose back to MATLAB order
