"""Reading MATLAB result and mask files."""
import numpy as np, scipy.io as sio


def load_soln(path):
    """Variable 'soln' of a *_results.mat file (MATLAB v7) as an object with attribute access."""
    return sio.loadmat(path, squeeze_me=True, struct_as_record=False)['soln']


def read_h5(path, name):
    """One variable of a MATLAB v7.3 (HDF5) file, e.g. the pipeline mask files, in MATLAB orientation."""
    import h5py
    with h5py.File(path, 'r') as f:
        return np.array(f[name]).T
