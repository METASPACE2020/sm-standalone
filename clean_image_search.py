import matplotlib
matplotlib.use('Agg')
import h5py
import numpy as np
import matplotlib.pyplot as plt

from pyMS.pyisocalc import pyisocalc
from pyIMS.inMemoryIMS import inMemoryIMS
from pyIMS.image_measures.isotope_image_correlation import *
from pyimzml.ImzMLParser import ImzMLParser

def read_mean_spectrum(raw_h5):
    for k in raw_h5['Regions'].keys():
        try:
            mzs = raw_h5['Regions'][k]['SamplePositions/SamplePositions'][:]
            total_intensities = raw_h5['Regions'][k]['Intensities'][:]
            break
        except:
            pass
    return mzs, total_intensities

def read_random_spectrum(raw_h5):
    spots = raw_h5['Spots'].keys()
    spot = np.random.choice(spots, 1)[0]
    mzs = raw_h5['Spots/' + spot + '/InitialMeasurement/SamplePositions/SamplePositions'][:]
    intensities = raw_h5['Spots/' + spot + '/InitialMeasurement/Intensities'][:]
    return mzs, intensities

def generate_summary_spectrum_lowmem(ppm, imzml):
    mz_min, mz_max = imzml.get_mz_range()

    mz_min = round(mz_min) - 1.0 
    mz_max = round(mz_max) + 1.0 
    mz_axis = mz_min * (1.0+ppm*1e-6) ** np.arange(np.log(mz_max / mz_min) / np.log(1.0+ppm*1e-6) + 1)

    mean_spectrum = np.zeros(mz_axis.shape)

    for i, _ in enumerate(imzml.coordinates):
        mzs, intensities = imzml.getspectrum(i)
        bins = np.floor((np.log(mzs) - np.log(mz_min)) / np.log(1.0+ppm*1e-6)).astype(int)
        mean_spectrum[bins] += intensities
    return mz_axis, mean_spectrum

def generate_summary_spectrum3(mz_axis, imzml):
    mean_spectrum = np.zeros(mz_axis.shape[0] + 1)
    frequencies = np.zeros(mz_axis.shape[0] + 1, dtype=np.int32)

    for i, _ in enumerate(imzml.coordinates):
        mzs, intensities = imzml.getspectrum(i)
        bins = mz_axis.searchsorted(mzs)
        mean_spectrum[bins] += intensities
        frequencies[bins] += 1
    mean_spectrum[-2] += mean_spectrum[-1]
    return mz_axis, mean_spectrum[:-1], frequencies[:-1]

#### computing resolution estimates ####

def fwhm_interval(peak_pos, mzs, intensities):
    l = peak_pos
    r = peak_pos + 1
    while l >= 0 and intensities[l] > 0.5 * intensities[peak_pos]: l -= 1
    while r + 1 < len(intensities) and intensities[r] > 0.5 * intensities[peak_pos]: r += 1
    return l, r
    
def resolution_at_peak(peak_pos, mzs, intensities):
    l, r = fwhm_interval(peak_pos, mzs, intensities)
    l_mz, r_mz = (mzs[l] + mzs[l+1]) / 2, (mzs[r] + mzs[r-1]) / 2
    if intensities[l+1] != intensities[l]:
        l_ratio = (0.5 * intensities[peak_pos] - intensities[l]) / (intensities[l+1] - intensities[l])
        l_mz = mzs[l] * l_ratio + mzs[l+0] * (1.0 - l_ratio)
    if intensities[r-1] != intensities[r]:
        r_ratio = (0.5 * intensities[peak_pos] - intensities[r]) / (intensities[r-1] - intensities[r])
        r_mz = mzs[r] * r_ratio + mzs[r-1] * (1.0 - r_ratio)
    mz_diff = r_mz - l_mz
    return np.average(mzs[peak_pos-2:peak_pos+3], weights=intensities[peak_pos-2:peak_pos+3]) / mz_diff

from pyMS.centroid_detection import gradient
from sklearn.linear_model import RANSACRegressor

def resolution_estimate(raw_h5, n_spectra=25):
    slopes = []
    intercepts = []
    for i in range(n_spectra):
        mzs, intensities = read_random_spectrum(raw_h5)
        peak_positions = np.array(gradient(mzs, intensities)[-1])
        intensities_at_peaks = intensities[peak_positions]
        high_intensity_threshold = np.percentile(intensities_at_peaks, 40)
        peak_positions = peak_positions[intensities[peak_positions] > high_intensity_threshold]
        resolutions = []
        for i, peak_pos in enumerate(peak_positions):
            resolutions.append(resolution_at_peak(peak_pos, mzs, intensities))
        resolutions = np.array(resolutions)
        mzs = mzs[peak_positions]
        mzs = mzs[resolutions > 0]
        resolutions = resolutions[resolutions > 0]
        ransac = RANSACRegressor()
        ransac.fit(np.log(mzs).reshape((-1,1)), np.log(resolutions).reshape((-1,1)))
        slope = ransac.estimator_.coef_[0][0]
        intercept = ransac.estimator_.intercept_[0]
        slopes.append(slope)
        intercepts.append(intercept)
    slope = np.median(slopes)
    intercept = np.median(intercepts)
    return lambda mz: np.exp(intercept + slope * np.log(mz)) 

def isotope_pattern_score(emp_intensities, theoretical_intensities):
    ints = theoretical_intensities[:len(emp_intensities)]
    emp_intensities = emp_intensities/np.linalg.norm(emp_intensities)
    ints = ints / np.linalg.norm(ints)
    return 1 - np.mean(abs(emp_intensities - ints))

def get_peak_indices(mz_axis, isotope_mzs):
    return mz_axis.searchsorted(isotope_mzs)

class Peak(object):
    def __init__(self, l, r, top, mz_axis, intensities):
        self.leftmost_bin = l
        self.rightmost_bin = r
        self.peak_bin = top
        self.top_intensity = intensities[top]
        self.top_mass = mz_axis[top]
        self.mz_interval = [mz_axis[l], mz_axis[r + 1]]
        self.total_intensity = intensities[l : r + 1].sum()

def get_peaks(mz_axis, total_intensities, indices):
    prev_top_intensity = None
    peaks = []
    for k, index in enumerate(indices):
        max_delta = (k + 2) / 2
        size = total_intensities.shape[0]
        if index < 1 or index > size - 1: break
        go_to_left = total_intensities[index - 1] > total_intensities[index]
        go_to_right = total_intensities[index + 1] > total_intensities[index]
        if go_to_left and go_to_right:
            break
        step = -1 if go_to_left else 1
        top = index
        delta = 0
        while 0 <= top + step < size and total_intensities[top] < total_intensities[top + step]:
            top += step
            delta += 1
        if delta > max_delta:
            break
        r = top
        while r + 1 < size and total_intensities[r + 1] <= total_intensities[r] and total_intensities[r + 1] > 0:
            r += 1
        l = top
        while l > 0 and total_intensities[l - 1] <= total_intensities[l] and total_intensities[l - 1] > 0:
            l -= 1
        if total_intensities[top] == 0.0:
            break
        u, v = total_intensities[l] / total_intensities[top], total_intensities[r] / total_intensities[top]
        if max(u, v) > 0.2:
            break

        #if prev_top_intensity and total_intensities[top] > prev_top_intensity:
        #    break
        prev_top_intensity = total_intensities[top]
        peaks.append(Peak(l, r, top, mz_axis, total_intensities))
    return peaks

#### main function of this notebook ####

from collections import defaultdict

def find_clean_molecules(mzs, total_intensities, patterns, min_peaks=3,
                         min_iso_corr=0.95, min_intensity_share=0.99):
    counts = defaultdict(int)

    all_masses = np.concatenate([patterns[key][0] for key in patterns])
    all_indices = mzs.searchsorted(all_masses)

    result = []
    offset = 0
    for key in patterns:
        f, a = key
        masses, intensities = patterns[key]
        k = 0
        indices = all_indices[offset:offset + len(masses)]
        offset += len(masses)
        peaks = get_peaks(mzs, total_intensities, indices)
        k = len(peaks)
        
        if k >= min_peaks:
            emp_intensities = np.array([peak.top_intensity for peak in peaks])
            emp_intensities /= np.sum(emp_intensities)
            ips = isotope_pattern_score(emp_intensities, intensities)
            # the desired properties of isotope peaks are:
            # a) good agreement with theoretical pattern
            # b) catching significant part of the total intensity
            if ips >= min_iso_corr and intensities[:k].sum() / intensities.sum() > min_intensity_share:
                counts[k] += 1
                theor_ints = intensities[:k] / intensities.sum()
                if (emp_intensities / theor_ints).max() > 1.2: # suggests that the corresponding peak is not clean
                    continue

                result.append((f, a))
    return result

def generate_patterns(formulas_fn, resolution_func, mz_range):
    mz_min, mz_max = mz_range
    patterns = {}
    adducts = ['H', 'K', 'Na']
    formulae = [s.strip() for s in open(formulas_fn).readlines()]
    for f in formulae:
        for a in adducts:
            sf = pyisocalc.SumFormulaParser.parse_string(f + a)
            raw_pattern = pyisocalc.isodist(sf, cutoff=1e-4, charge=1)
            mz = raw_pattern.get_spectrum()[0][0]
            if mz < mz_min or mz > mz_max:
                continue
            fwhm = mz / resolution_func(mz)
            mzs, intensities = pyisocalc.apply_gaussian(raw_pattern, fwhm, exact=False).get_spectrum(source="centroids")
            mzs = np.array(mzs)
            intensities = np.array(intensities)
            order = np.argsort(intensities)[::-1]
            patterns[(f, a)] = (mzs[order], intensities[order])
    return patterns

class SpectralMatch(object):
    def __init__(self, f, a, patterns, mzs, intensities, resolution_func):
        self.formula = f
        self.adduct = a
        self.theor_mzs, self.theor_ints = patterns[(f, a)]
        indices = get_peak_indices(mzs, self.theor_mzs)
        self.peaks = get_peaks(mzs, intensities, indices)
        self.emp_intensities = np.array([p.top_intensity for p in self.peaks])
        self.ips = isotope_pattern_score(self.emp_intensities, self.theor_ints)
        self._resolution = resolution_func
        
    @property
    def peak_count(self):
        return len(self.peaks)
    
    @property
    def resolution(self):
        return self._resolution(self.theor_mzs[0])
    
    @property
    def isotope_pattern_score(self):
        return self.ips
    
    def __str__(self):
        return " + ".join((self.formula, self.adduct))

def get_images(matches, imzml, n_peak_images, mz_axis, n_bins=9):
    intervals = []
    for m in matches:
        assert len(m.theor_mzs) >= n_peak_images
        k = 0
        for p in m.peaks:
            intervals.append(p.mz_interval)
            k += 1
            if k == n_peak_images: # enough!
                break
        while k < n_peak_images:
            mz = m.theor_mzs[k]
            central_bin = mz_axis.searchsorted(mz)
            intervals.append([mz_axis[central_bin - n_bins/2], mz_axis[central_bin + n_bins/2]])
            k += 1
    #intervals = [p.mz_interval for m in matches for p in m.peaks]
    lower, upper = map(np.array, zip(*intervals))
    
    nr = imzml.imzmldict["max count of pixels x"]
    nc = imzml.imzmldict["max count of pixels y"]

    images = np.zeros((len(intervals), nr * nc))
    for i, (row, col, _) in enumerate(imzml.coordinates):
        mzs, intensities = map(np.array, imzml.getspectrum(i))
        cumul_ints = np.concatenate(([0.0], np.cumsum(intensities)))
        k = (col - 1) * nr + (row - 1)
        lidx = mzs.searchsorted(lower, 'l')
        ridx = mzs.searchsorted(upper, 'r')
        images[:, k] = cumul_ints[ridx] - cumul_ints[lidx]
    return images, nc, nr
 

# main program

import logging
reload(logging)
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] %(message)s', datefmt='%H:%M:%S')

import sys
h5 = h5py.File(sys.argv[1])
imzml = ImzMLParser(sys.argv[2])
formulas_fn = sys.argv[3]
pdf_filename = sys.argv[4]

n_spectra = 25
logging.info("estimating resolution from %d random raw spectra..." % n_spectra)
resolution_func = resolution_estimate(h5)
logging.info("resolution is %d @ 200" % round(resolution_func(200)))

mz_range = imzml.get_mz_range()
logging.info("m/z range: %f .. %f" % mz_range)

logging.info("generating isotope patterns...")
patterns = generate_patterns(formulas_fn, resolution_func, mz_range)

logging.info("computing mean spectrum from centroided data...")
mzs, _ = read_mean_spectrum(h5)
mzs, intensities, frequencies = generate_summary_spectrum3(mzs, imzml)

logging.info("finding candidate molecules")
matches = []
molecules = find_clean_molecules(mzs, intensities, patterns, min_peaks=3, min_intensity_share=0.97,\
                                 min_iso_corr=0.9)

for f, a in sorted(molecules, key=lambda key:patterns[key][0][0]):
    matches.append(SpectralMatch(f, a, patterns, mzs, intensities, resolution_func))
    
logging.info("extracting molecular images from imzml...")
n_bins = 15
n = 5
matches = [m for m in matches if len(m.theor_mzs) >= n]
images, nrow, ncol = get_images(matches, imzml, n, mzs, n_bins)

from matplotlib import gridspec
from colourmaps import viridis_colormap

def set_axis_color(ax, c):
    ax.tick_params(color=c, labelcolor=c)
    for spine in ax.spines.values():
        spine.set_edgecolor(c)
        
_, mean_intensities = read_mean_spectrum(h5)

viridis = viridis_colormap()
offset = 0

from matplotlib.backends.backend_pdf import PdfPages

logging.info("saving results to pdf...")

plt.ioff()

with PdfPages(pdf_filename) as pdf:
    for match in sorted(matches, key = lambda m: m.theor_mzs[0]):
        k = match.peak_count
        if k < 3:
            offset += n
            continue

        l = min(k, len(match.theor_ints))
        try:
            ims = isotope_image_correlation(images[offset : offset + l, :], weights=match.theor_ints[1:l])
        except TypeError:
            ims = -1 # TODO

        #if ims < 0.7:
        #    offset += n
        #    continue

        #fig = plt.figure(figsize=(3.3 * n, 13))
        fig = plt.figure(figsize=(8.27, 11.69), dpi=100) # A4 format
        for k in range(n):
            gs = gridspec.GridSpec(6, n, height_ratios=[3, 5, 2.5, 2.5, 2.5, 1.5])

            img = images[offset + k,:].reshape((nrow, ncol))
            perc = np.percentile(img, 99)
            img[img > perc] = perc
            plt.subplot(gs[n + k])
            plt.imshow(img, cmap=viridis)
            plt.title("%.4f" % (round(match.theor_mzs[k], 4)))
            plt.axis('off', frameon=False)

            if k < match.peak_count:
                l, r = match.peaks[k].leftmost_bin, match.peaks[k].rightmost_bin
                ax_handler = lambda ax: None
            else:
                central_bin = mzs.searchsorted(match.theor_mzs[k])
                l = central_bin - n_bins/2
                r = central_bin + n_bins/2 - 1
                ax_handler = lambda ax: set_axis_color(ax, 'red')

            # plot mean spectrum
            ax = plt.subplot(gs[2 * n + k])
            plt.plot(mzs[l:r+1], mean_intensities[l:r+1], '.-')
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            ax.get_xaxis().set_visible(False)
            ax.axvline(match.theor_mzs[k], color='green')
            ax_handler(ax)

            # plot spectrum around the peaks
            ax = plt.subplot(gs[3 * n + k])
            plt.plot(mzs[l:r+1], intensities[l:r+1], '.-')
            ax.get_xaxis().set_visible(False)
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            ax.axvline(match.theor_mzs[k], color='green')
            ax_handler(ax)

            # plot pixel counts
            ax = plt.subplot(gs[4 * n + k])
            plt.plot(mzs[l:r+1], frequencies[l:r+1], '.-')
            ax.get_xaxis().set_visible(False)
            ax.axvline(match.theor_mzs[k], color='green')
            ax_handler(ax)

        offset += n
        #plt.show()
        
        emp_ints = match.emp_intensities / match.emp_intensities.sum()
        k = match.peak_count
        txt = '''
        {0} + {1}
        First {2} peaks are clean
        m/z's: {3}
        intensities (theoretical): {4}
        intensities (empirical): {5}
        isotope pattern score: {6}
        image correlation score (only for clean peaks): {7}
        clean peak intensity share: {8}
        '''.format(match.formula, match.adduct,
                   match.peak_count,
                   match.theor_mzs[:k].round(4),
                   match.theor_ints[:k].round(1),
                   (emp_ints * 100 / emp_ints[0]).round(1),
                   match.isotope_pattern_score,
                   ims,
                   match.theor_ints[:k].sum() / match.theor_ints.sum())
        
        fig.text(.1, .85, txt)
        fig.tight_layout()

        pdf.savefig(fig)
        plt.close(fig)

logging.info("done!")
