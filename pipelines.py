from pyIMS.image_measures.isotope_pattern_match import isotope_pattern_match
from pyIMS.image_measures.isotope_image_correlation import isotope_image_correlation
from pyIMS.image_measures.level_sets_measure import measure_of_chaos
from pyIMS.hdf5.inMemoryIMS_hdf5 import inMemoryIMS_hdf5
from pyIMS.ion_datacube import ion_datacube

from pyMS.pyisocalc import pyisocalc

from itertools import product
import os
import sys
import cPickle
import numpy as np
import matplotlib.image
import logging
import h5py

import scipy.signal as signal
from matplotlib.colors import Normalize

logging.basicConfig(level=logging.INFO, format='[%(asctime)s] %(message)s', datefmt='%H:%M:%S')

class Pipeline(object):
    def __init__(self, config):
        self.config = config
        self.q = config['image_generation']['q'] 
        self.ppm = config['image_generation']['ppm'] #parts per million -  a measure of how accurate the mass spectrometer is
        self.nlevels = config['image_generation']['nlevels']  # parameter for measure of chaos
        self.data_file = config['file_inputs']['data_file']

        self.measure_value_score = {}
        self.iso_correlation_score = {}
        self.iso_ratio_score = {}

        self.chunk_size = 200

        from colourmaps import viridis_colormap
        self.cmap = viridis_colormap()

        self.measure_tol = config['results_thresholds']['measure_tol']
        self.iso_corr_tol = config['results_thresholds']['iso_corr_tol']
        self.iso_ratio_tol = config['results_thresholds']['iso_ratio_tol']

    def run(self):
        logging.info("==== computing/loading isotope patterns")
        self.load_queries()
        logging.info("==== loading data")
        self.load_data()
        logging.info("==== computing scores for all formulae")
        self.compute_scores()
        logging.info("==== outputting results")
        self.print_results()

    def make2DImage(self, img):
        # Regions with no data are denoted as -1.
        # This allows to apply measure of chaos to the image
        # and quickly filter the pixels out when producing PNGs
        result = np.zeros(self.nrows * self.ncols) - 1
        result[self.pixel_indices] = img
        return result.reshape((self.nrows, self.ncols))

    def save_image(self, img, filename_out):
        mask = img >= 0
        values = img[mask]
        norm = Normalize(vmin=np.min(values), vmax=np.max(values))
        colorized_img = np.zeros((self.nrows, self.ncols, 4))
        colorized_img[mask] = self.cmap(norm(values))
        # set alpha channel to 0 for pixels with no data
        colorized_img[img < 0, -1] = 0
        matplotlib.image.imsave(filename_out, colorized_img)

    def print_images(self, imgs, sum_formula, adduct):
        total_img = np.zeros((self.nrows, self.ncols))

        img_output_dir = self.output_directory()
        for i, mz in enumerate(self.mz_list[sum_formula][adduct][0]):
            filename_out = "{img_output_dir}/{sum_formula}_{adduct}_{mz}.png".format(**locals())
            with open(filename_out,'w') as f_out:
                img = self.make2DImage(imgs[i])
                total_img += img
                self.save_image(img, filename_out)
        
        filename_out = "{img_output_dir}/_{sum_formula}_{adduct}.png".format(**locals())
        self.save_image(total_img, filename_out)

    # template method
    def compute_scores(self):
        ### Runs the main pipeline
        # Get sum formula and predicted m/z peaks for molecules in database
        # Parse dataset
        raise NotImplementedError

    def load_data(self):
        raise NotImplementedError

    def algorithm_name(self):
        raise NotImplementedError

    # creates the output directory if it doesn't exist
    def output_directory(self):
        output_dir = self.config['file_inputs']['results_folder']
        if os.path.isdir(output_dir) == False:
            os.mkdir(output_dir)
        return output_dir

    # while all_images can have whatever shape, first_image is used for chaos measure
    def process_query(self, sum_formula, adduct, all_images, first_image, intensities):
        self.score_chaos(sum_formula, adduct, first_image)
        self.score_corr(sum_formula, adduct, all_images, intensities[1:])
        self.score_ratio(sum_formula, adduct, all_images, intensities)

    def hot_spot_removal(self, xics):
        for xic in xics:
            xic_q = np.percentile(xic, self.q)
            xic[xic > xic_q] = xic_q
        return xics

    def score_chaos(self, sum_formula, adduct, img):
        if not sum_formula in self.measure_value_score:
            self.measure_value_score[sum_formula] = {}
        result = 1 - measure_of_chaos(img, self.nlevels, interp=False)[0]
        if result == 1:
            result = 0
        self.measure_value_score[sum_formula][adduct] = result

    def score_corr(self, sum_formula, adduct, imgs, weights):
        if not sum_formula in self.iso_correlation_score:
            self.iso_correlation_score[sum_formula] = {}
        result = 1.0        # return 1 if there's a single peak
        if len(weights) > 1:
            result = isotope_image_correlation(imgs, weights=weights)
        self.iso_correlation_score[sum_formula][adduct] = result

    def score_ratio(self, sum_formula, adduct, imgs, intensities):
        if not sum_formula in self.iso_ratio_score:
            self.iso_ratio_score[sum_formula] = {}
        self.iso_ratio_score[sum_formula][adduct] = isotope_pattern_match(imgs, intensities)

    # config.file_inputs.database_file must contain one formula per line
    def load_queries(self):
        def calculate_isotope_patterns(sum_formulae,adducts='',isocalc_sig=0.01,isocalc_resolution = 200000.,isocalc_do_centroid = True, charge=1):
            ### Generate a mz list of peak centroids for each sum formula with the given adduct
            # todo - parse sum formula and adduct properly so that subtractions (losses) can be utilised (this code already exists somewhere)
            mz_list={}
            for n, sum_formula in enumerate(sum_formulae):   
                isotope_ms = pyisocalc.isodist(sum_formula + adduct,
                                               plot=False,
                                               sigma=isocalc_sig,
                                               charges=charge,
                                               resolution=isocalc_resolution,
                                               do_centroid=isocalc_do_centroid)
                if not sum_formula in mz_list:
                     mz_list[sum_formula] = {}
                mz_list[sum_formula][adduct] = isotope_ms.get_spectrum(source='centroids')
            return mz_list

        # Extract variables from config dict
        config = self.config
        db_filename = config['file_inputs']['database_file']
        db_dump_folder = config['file_inputs']['database_load_folder']  
        isocalc_sig = config['isotope_generation']['isocalc_sig']  
        isocalc_resolution = config['isotope_generation']['isocalc_resolution']  
        if len(config['isotope_generation']['charge']) > 1:
            print 'Warning: only first charge state currently accepted'
        charge = int('{}{}'.format(config['isotope_generation']['charge'][0]['polarity'], config['isotope_generation']['charge'][0]['n_charges'])) #currently only supports first charge!!
        self.adducts=[a['adduct'] for a in config['isotope_generation']['adducts']]
      
        # Read in molecules
        self.sum_formulae = [l.strip() for l in open(db_filename).readlines()]
        # Check if already generated and load if possible, otherwise calculate fresh   
        db_name =  os.path.splitext(os.path.basename(db_filename))[0] 
        mz_list={}
        nmz = 0
        nf = 0
        for adduct in self.adducts:
            load_file = '{}/{}_{}_{}_{}.dbasedump'.format(db_dump_folder,db_name,adduct,isocalc_sig,isocalc_resolution)
            if os.path.isfile(load_file):
                logging.info("loading cached isotope patterns for adduct %s" % adduct)
                mz_list_tmp = cPickle.load(open(load_file,'r'))
            else:
                logging.info("calculating isotope patterns for adduct %s" % adduct)
                mz_list_tmp = calculate_isotope_patterns(self.sum_formulae,adducts=(adduct,),isocalc_sig=isocalc_sig,isocalc_resolution=isocalc_resolution,charge=charge)
                if db_dump_folder != "":
                    if os.path.isdir(db_dump_folder)==False:
                        os.mkdir(db_dump_folder)
                    with open(load_file, 'w') as f:
                        cPickle.dump(mz_list_tmp, f)
            # add patterns to total list
            for sum_formula in mz_list_tmp:
                if not sum_formula in mz_list:
                    mz_list[sum_formula] = {}
                mzs, ints = mz_list_tmp[sum_formula][adduct]
                order = ints.argsort()[::-1]
                mz_list[sum_formula][adduct] = (mzs[order], ints[order])
                nmz += len(ints)
                nf += 1 
        self.mz_list = mz_list
        logging.info("all isotope patterns generated and loaded (#formula+adduct pairs: {nf}, #peaks: {nmz})".format(**locals()))

    def passes_filters(self, sum_formula, adduct):
        def ok(dictionary, threshold):
            return dictionary[sum_formula][adduct] > threshold

        return ok(self.measure_value_score, self.measure_tol)\
           and ok(self.iso_correlation_score, self.iso_corr_tol)\
           and ok(self.iso_ratio_score, self.iso_ratio_tol)

    def print_results(self):
        filename_in = self.config['file_inputs']['data_file']
        output_dir = self.config['file_inputs']['results_folder']
        # Save the processing results
        if os.path.isdir(output_dir)==False:
            os.mkdir(output_dir)
        filename_out = '{}{}{}_full_results_{}.txt'.format(output_dir,os.sep,os.path.splitext(os.path.basename(filename_in))[0], self.algorithm_name())
        with open(filename_out,'w') as f_out:
            f_out.write('sf,adduct,mz,moc,spec,spat,pass\n'.format())
            for sum_formula, adduct in product(self.sum_formulae, self.adducts):
                moc_pass = self.passes_filters(sum_formula, adduct)
                f_out.write('{},{},{},{},{},{},{}\n'.format(
                        sum_formula,
                        adduct,
                        self.mz_list[sum_formula][adduct][0][0],
                        self.measure_value_score[sum_formula][adduct],
                        self.iso_correlation_score[sum_formula][adduct],
                        self.iso_ratio_score[sum_formula][adduct],
                        moc_pass)) 

        filename_out = '{}{}{}_pass_results_{}.txt'.format(output_dir,os.sep,os.path.splitext(os.path.basename(filename_in))[0], self.algorithm_name())
        with open(filename_out,'w') as f_out:
            f_out.write('sf,adduct,mz,moc,spec,spat\n'.format())
            for sum_formula, adduct in product(self.sum_formulae, self.adducts):
                if self.passes_filters(sum_formula, adduct):
                    f_out.write('{},{},{},{},{},{}\n'.format(
                        sum_formula, adduct,
                        self.mz_list[sum_formula][adduct][0][0],
                        self.measure_value_score[sum_formula][adduct],
                        self.iso_correlation_score[sum_formula][adduct],
                        self.iso_ratio_score[sum_formula][adduct]))

    def report_scoring_progress(self, n_processed):
        #logging.info(str(round(float(n_processed) * 100.0 / len(self.sum_formulae), 2)) + "% sum formulae processed")
        logging.info(str(n_processed) + "/" + str(len(self.sum_formulae)) + " sum formulae processed")

    def _calculate_dimensions(self):
        cube = ion_datacube()
        cube.add_coords(self.coords)
        self.nrows = cube.nRows
        self.ncols = cube.nColumns
        self.pixel_indices = cube.pixel_indices

class ReferencePipeline(Pipeline):
    def __init__(self, config):
        super(ReferencePipeline, self).__init__(config)

    def algorithm_name(self):
        return "reference"

    def load_data(self):
        self.IMS_dataset = inMemoryIMS_hdf5(self.data_file)
        self.coords = self.IMS_dataset.coords
        self._calculate_dimensions()

    def compute_scores(self):
        for i, sum_formula in enumerate(self.sum_formulae):
            if i % self.chunk_size == 0 and i > 0:
                self.report_scoring_progress(i)
            for adduct in self.adducts:
                ion_datacube = self.IMS_dataset.get_ion_image(self.mz_list[sum_formula][adduct][0], self.ppm) #for each spectrum, sum the intensity of all peaks within tol of mz_list
                ion_datacube.xic = self.hot_spot_removal(ion_datacube.xic)

                img = ion_datacube.xic_to_image(0)
                intensities = self.mz_list[sum_formula][adduct][1]

                self.process_query(sum_formula, adduct, ion_datacube.xic, img, intensities)

                if self.passes_filters(sum_formula, adduct):
                    assert np.allclose(ion_datacube.xic_to_image(0), self.make2DImage(ion_datacube.xic[0]))
                    self.print_images(ion_datacube.xic, sum_formula, adduct)

        self.report_scoring_progress(len(self.sum_formulae))

class inMemoryIMS_low_mem(inMemoryIMS_hdf5):
    def __init__(self, filename):
        self.load_file(filename)

    def load_file(self, filename, min_mz=0, max_mz=np.inf, min_int=0, index_range=[]):
        if filename.endswith('.hdf5'):
            self.hdf = h5py.File(filename, 'r')
            keys = index_range or map(int, self.hdf['/spectral_data'].keys())
        else:
            self.imzml = ImzMLParser.ImzMLParser(filename)
            keys = index_range or range(len(self.imzml.coordinates))

        self.coords = np.zeros((len(keys), 3))

        def spectra_iter_hdf5(keys):
            for i in keys:
                tmp_str = "/spectral_data/" + str(i)
                mzs = self.hdf[tmp_str + '/centroid_mzs/'][()]
                counts = self.hdf[tmp_str + '/centroid_intensities/'][()]
                self.coords[i, :] = self.hdf[tmp_str + '/coordinates']
                valid = np.where((mzs > min_mz) & (mzs < max_mz) & (counts > min_int))
                counts = counts[valid]
                mzs = mzs[valid]
                yield (i, mzs, counts)

        def spectra_iter_imzml(keys):
            for i in keys:
                coords = self.imzml.coordinates[i]
                mzs, counts = map(np.array, self.imzml.getspectrum(i))
                if len(coords) == 2:
                    coords = (coords[0], coords[1], 0)
                self.coords[i, :] = coords
                yield (i, mzs, counts)

        sp_iter = spectra_iter_hdf5 if filename.endswith('.hdf5') else spectra_iter_imzml

        import reorder
        data = reorder.sortDatasetByMass(sp_iter(keys))
        self.index_list, self.mz_list, self.count_list, self.idx_list = data
        self.max_index = max(self.index_list)

        cube = ion_datacube()
        cube.add_coords(self.coords)
        self.cube_pixel_indices = cube.pixel_indices
        self.cube_n_row, self.cube_n_col = cube.nRows, cube.nColumns

        self.mz_min = self.mz_list[0]
        self.mz_max = self.mz_list[-1]
        self.histogram_mz_axis = {}

        # split binary searches into two stages for better locality
        self.window_size = 1024
        self.mz_sublist = self.mz_list[::self.window_size].copy()

class LowMemReferencePipeline(ReferencePipeline):
    def __init__(self, config):
        super(LowMemReferencePipeline, self).__init__(config)

    def algorithm_name(self):
        return "reference_low_mem"

    def load_data(self):
        self.IMS_dataset = inMemoryIMS_low_mem(self.data_file)
        self.coords = self.IMS_dataset.coords
        self._calculate_dimensions()

#####################################################################################################################

# a few helper functions used by the NewPipeline

from pyMS import centroid_detection
def prepare(mzs, ints, centroids=True):
    if centroids == True:
        mzs_list, intensity_list = mzs, ints
    else:    
        ints=signal.savgol_filter(ints, 5, 2)
        mzs_list, intensity_list, indices_list = \
            centroid_detection.gradient(np.asarray(mzs), np.asarray(ints), max_output=-1, weighted_bins=3)
    mzs_list = np.asarray(mzs_list).astype(np.float64)
    intensity_list = np.asarray(intensity_list).astype(np.float32)
    intensity_list[intensity_list < 0] = 0
    return mzs_list, intensity_list

from collections import namedtuple
Spectrum = namedtuple('Spectrum', ['index', 'mzs', 'cumsum_int', 'coords'])

class Spectrum(object):
    def __init__(self, i, mzs, intensities, coords):
        self.index = int(i)
        self.mzs = np.asarray(mzs)
        self.cumsum_int = np.cumsum(np.concatenate(([0], intensities)))
        self.coords = np.asarray(coords)

from pyimzml import ImzMLParser
def readImzML(filename, centroids=True):
    f_in = ImzMLParser.ImzMLParser(filename)       
    for i, coords in enumerate(f_in.coordinates):
        mzs, ints = prepare(*f_in.getspectrum(i), centroids=centroids)
        if len(coords) == 2:
            coords = (coords[0], coords[1], 0)
        yield Spectrum(i, mzs, ints, map(lambda x: x-1, coords))

def readHDF5(filename, centroids=True):
    hdf = h5py.File(filename, 'r')
    for i in hdf['/spectral_data'].keys():
        tmp_str = "/spectral_data/" + i
        mzs = hdf[tmp_str + '/centroid_mzs/']
        ints = hdf[tmp_str + '/centroid_intensities/']
        coords = hdf[tmp_str + '/coordinates/']
        mzs, ints = prepare(mzs, ints, centroids=centroids)
        yield Spectrum(i, mzs, ints, map(lambda x: x-1, coords))

class NewPipeline(Pipeline):
    def __init__(self, config):
        super(NewPipeline, self).__init__(config)

    def algorithm_name(self):
        return "new"

    def load_data(self):
        if self.data_file.endswith(".imzML"):
            spectra = readImzML(self.data_file)
        elif self.data_file.endswith(".hdf5"):
            spectra = readHDF5(self.data_file)
        else:
            raise "the input format is unsupported"

        self.spectra = list(spectra)

        self.coords = np.zeros((len(self.spectra), 3))
        for i, sp in enumerate(self.spectra):
            self.coords[i, :] = sp.coords
        self._calculate_dimensions()

    def compute_scores(self):
        chunk_size = self.chunk_size
        n_chunks = len(self.sum_formulae) / chunk_size + 1
        for offset in xrange(0, len(self.sum_formulae), chunk_size):
            formulae = self.sum_formulae[offset : offset+chunk_size]
            r = self.get_ion_images(self.spectra, formulae)

            for xic, (sum_formula, adduct) in zip(r, product(formulae, self.adducts)):
                imgs = self.hot_spot_removal(xic)
                img = self.make2DImage(imgs[0])
                intensities = self.mz_list[sum_formula][adduct][1]
                self.process_query(sum_formula, adduct, imgs, img, intensities)

                if self.passes_filters(sum_formula, adduct):
                    self.print_images(imgs, sum_formula, adduct)
            del r
            import gc
            gc.collect()
            self.report_scoring_progress(offset + len(formulae))

    def process_spectra_multiple_queries(self, mol_mz_intervals, spectra):
        from numba import njit
        @njit
        def numba_multiple_queries(lower, upper, lperm, uperm, mzs, cumsum_int,
                                   result, pixel, tmp1, tmp2):
            m = len(mzs)
            n = len(lower)

            i = j = 0
            while i < n:
                x = lower[i]
                while j < m and x > mzs[j]:
                    j += 1
                tmp1[lperm[i]] = j
                i += 1

            while i > 0:
                i -= 1
                x = upper[i]
                while j > 0 and mzs[j - 1] > x:
                    j -= 1
                tmp2[uperm[i]] = j
                    
            i = 0
            while i < n:
                r = cumsum_int[tmp2[i]] - cumsum_int[tmp1[i]]
                result[i, pixel] += r
                i += 1

        lower, upper, lperm, uperm = mol_mz_intervals
        n = len(lower)
        
        result = np.zeros((n, len(spectra)))

        query_ids = np.zeros(n, dtype=np.int)
        intensities = np.zeros(n)
        tmp1 = np.zeros(n, dtype=np.int)
        tmp2 = np.zeros(n, dtype=np.int)

        for i, sp in enumerate(spectra):
            numba_multiple_queries(
                    lower, upper, lperm, uperm, sp.mzs, sp.cumsum_int,
                    result, i, tmp1, tmp2)

        return result

    def get_ion_images(self, spectra, formulae):
        peaks = [self.mz_list[f][a][0] for f in formulae for a in self.adducts]
        query_lens = np.array(map(len, peaks))
        mzs = np.array([s for _q in peaks for s in _q])
        tols = mzs * self.ppm / 1e6
        lower = mzs - tols
        upper = mzs + tols
        lower_order = lower.argsort()
        lower_sorted = lower[lower_order]
        upper_order = upper.argsort()
        upper_sorted = upper[upper_order]
        bounds = (lower_sorted, upper_sorted, lower_order, upper_order)
        qres = self.process_spectra_multiple_queries(bounds, spectra)
        qres = np.split(qres, np.cumsum(query_lens)[:-1])
        return qres

if __name__ == '__main__':
    import json
    import sys
    config = json.loads(open(sys.argv[1]).read())
    if config['method'] == 'reference':
        pipeline = ReferencePipeline(config)
    elif config['method'] == 'reference_low_mem':
        pipeline = LowMemReferencePipeline(config)
    elif config['method'] == 'new':
        pipeline = NewPipeline(config)
    else:
        print "method not recognized"
        sys.exit(1)
    pipeline.run()

