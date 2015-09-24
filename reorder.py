import pyximport
pyximport.install()

import utils

import numpy as np
import array

def sortDatasetByMass(spectra):
    """
    Orders all data points by m/z for fast image generation.
    Does so in a memory-efficient fashion.

    Parameters
    ----------
    spectra: iterator over (index, mzs, intensities) tuples

    Returns
    -------
    ret : (ndarray, ndarray, ndarray)
        Arrays containing ordered m/z values, intensities and 
        spectrum indices.
    """

    mzs_list = None
    count_list = None
    indices = []
    reps = []
    mzs_dtype = intensities_dtype = None

    def typecode(numpy_array):
        dtype = numpy_array.dtype
        if dtype == np.float64:
            return 'd'
        elif dtype == np.float32:
            return 'f'
        else:
            raise "unsupported data type"

    for index, mzs, intensities in spectra:
        indices.append(index)
        reps.append(len(mzs))
        if not mzs_list:
            mzs_dtype = mzs.dtype
            intensities_dtype = intensities.dtype
            mzs_list = array.array(typecode(mzs))
            count_list = array.array(typecode(intensities))
        mzs_list.extend(mzs)
        count_list.extend(intensities)

    # convert array.arrays into numpy arrays (without copying)
    mzs_list = np.frombuffer(mzs_list, mzs_dtype)
    count_list = np.frombuffer(count_list, intensities_dtype)

    # np.int64 branch will probably never be tested...
    order_dtype = np.int64 if len(mzs_list) >= 2**31 else np.int32
    order = np.arange(len(mzs_list), dtype=order_dtype)

    # sort mzs_list and set the permutation accordingly
    utils.keysort(mzs_list, order)

    # rearrange intensities and indices as well
    indices = np.array(indices, np.int32)
    reps = np.array(reps, np.uint64)
    utils.permute(count_list, order)
    utils.permuterepeat(reps, indices, order)
    idx_list = order

    return indices, mzs_list, count_list, idx_list
