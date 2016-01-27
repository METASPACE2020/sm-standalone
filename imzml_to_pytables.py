from pyimzml.ImzMLParser import *
import tables as t
import numpy as np
from progressbar import ProgressBar
import pyximport
pyximport.install()
from utils import keysort

import heapq
import os
import gc

class Peak(t.IsDescription):
    mz = t.Float64Col()
    intensity = t.Float32Col()
    x = t.Int16Col()
    y = t.Int16Col()

def convert(imzml_filename, output_filename, complevel=5):
    filters = t.Filters(complib='blosc', complevel=complevel)

    imzml = ImzMLParser(imzml_filename)
        
    chunk_size = int(10e6)
    convert.mzs_tmp = np.zeros(chunk_size, dtype=np.float64)
    convert.ints_tmp = np.zeros(chunk_size, dtype=np.float32)
    convert.xs_tmp = np.zeros(chunk_size, dtype=np.uint16)
    convert.ys_tmp = np.zeros(chunk_size, dtype=np.uint16)

    convert.n_rows = 0
    convert.part_number = 0

    convert.tmp_filenames = []

    convert.order_tmp = np.zeros(chunk_size, dtype=np.uint32)

    def write(peak, order):
        for j in xrange(len(order)):
            peak['mz'] = convert.mzs_tmp[order[j]]
            peak['intensity'] = convert.ints_tmp[order[j]]
            peak['x'] = convert.xs_tmp[order[j]]
            peak['y'] = convert.ys_tmp[order[j]]
            peak.append()

    def fill_arange(array, n):
        for j in xrange(n):
            array[j] = j

    def append_peaks(mzs):
        convert.tmp_filenames.append(output_filename + ".part" + str(convert.part_number))
        h5file = t.open_file(convert.tmp_filenames[-1], mode = "w")
        table = h5file.create_table("/", 'peaks', Peak, "Peaks",
                                    filters=t.Filters(complib='blosc', complevel=1),
                                    expectedrows=chunk_size)
        n = len(mzs)
        fill_arange(convert.order_tmp, n)
        keysort(mzs, convert.order_tmp[:n])
        order = convert.order_tmp[:n]

        write(table.row, order)

        table.flush()
        convert.part_number += 1
        convert.n_rows += len(mzs)
        h5file.close()

    print "Writing sorted temporary files..."

    convert.k = 0

    def add_spectrum(spectrum, coords):
        k = convert.k
        for j in xrange(len(spectrum[0])):
            mz, intensity = spectrum[0][j], spectrum[1][j]
            convert.mzs_tmp[k] = mz
            convert.ints_tmp[k] = intensity
            convert.xs_tmp[k] = coords[0]
            convert.ys_tmp[k] = coords[1]
            convert.k += 1
            if convert.k == chunk_size:
                append_peaks(convert.mzs_tmp)
                convert.k = 0

    def run(imzml):
        pbar = ProgressBar(maxval=len(imzml.coordinates)).start()
        i = 0
        for coords in imzml.coordinates:
            #spectrum = map(np.array, imzml.getspectrum(i))
            spectrum = imzml.getspectrum(i)
            add_spectrum(spectrum, coords)
            i += 1
            if i % 100 == 0:
                pbar.update(i)
                gc.collect()
        append_peaks(convert.mzs_tmp[:convert.k])
        pbar.finish()
                    
    run(imzml)

    files = []
    records = []
    merged = t.open_file(output_filename, mode='w')
    table = merged.create_table("/", 'peaks', Peak, "Peaks", filters=filters, expectedrows=convert.n_rows)
    peak = table.row

    nx = imzml.imzmldict['max count of pixels x']
    ny = imzml.imzmldict['max count of pixels y']
    merged.create_array("/", 'dimensions', np.array([nx, ny]), "Scan dimensions")

    for fn in convert.tmp_filenames:
        files.append(t.open_file(fn))
        iterator = files[-1].root.peaks.iterrows()
        records.append(((x['mz'], x['intensity'], x['x'], x['y']) for x in iterator))

    print "Merging temporary files"
    pbar = ProgressBar(maxval=convert.n_rows).start()
    k = 0
    for mz, intensity, x, y in heapq.merge(*records):
        peak['mz'] = mz
        peak['intensity'] = intensity
        peak['x'] = x
        peak['y'] = y
        peak.append()
        k += 1
        if k % 8192 == 0:
            pbar.update(k)
            gc.collect()
    pbar.finish()

    table.flush()
    table.cols.mz.create_index(filters=filters)
    merged.close()

    for f in files:
        f.close()
    for fn in convert.tmp_filenames:
        os.unlink(fn)

def read_mz_image(imh5_filename, mz, ppm, hotspot_removal=True):
    h5file = t.open_file(imh5_filename, mode = "r")

    peaks = h5file.root.peaks
    lmz = mz * (1 - ppm * 1e-6)
    rmz = mz * (1 + ppm * 1e-6)

    img = np.zeros((250, 280))

    # uses pyTables indexing for fast retrieval
    for x in peaks.where("(mz > lmz) & (mz < rmz)"):
        img[x['x'] - 1, x['y'] - 1] += x['intensity'] 
    h5file.close()
    
    if hotspot_removal:
        perc = np.percentile(img, 99)
        img[img > perc] = perc
    return img

if __name__ == '__main__':
    import sys
    convert(sys.argv[1], sys.argv[2])
