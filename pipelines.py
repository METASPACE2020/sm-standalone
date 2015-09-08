import sys

from pyIMS.image_measures.isotope_pattern_match import isotope_pattern_match
from pyIMS.image_measures.isotope_image_correlation import isotope_image_correlation
from pyIMS.image_measures.level_sets_measure import measure_of_chaos
from pyMS.pyisocalc import pyisocalc
from itertools import product
import os
import sys
import cPickle
import numpy as np

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

        self.measure_tol = config['results_thresholds']['measure_tol']
        self.iso_corr_tol = config['results_thresholds']['iso_corr_tol']
        self.iso_ratio_tol = config['results_thresholds']['iso_ratio_tol']

    def run(self):
        print "==== computing/loading isotope patterns"
        self.load_queries()
        print "==== computing scores for all formulae"
        self.compute_scores()
        print "==== outputting results"
        self.print_results()

    # template method
    def compute_scores(self):
        ### Runs the main pipeline
        # Get sum formula and predicted m/z peaks for molecules in database
        # Parse dataset
        raise NotImplementedError

    def algorithm_name(self):
        raise NotImplementedError

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
        for adduct in self.adducts:
            load_file = '{}/{}_{}_{}_{}.dbasedump'.format(db_dump_folder,db_name,adduct,isocalc_sig,isocalc_resolution)
            if os.path.isfile(load_file):
                print "loading cached isotope patterns for adduct", adduct
                mz_list_tmp = cPickle.load(open(load_file,'r'))
            else:
                print "calculating isotope patterns for adduct", adduct
                mz_list_tmp = calculate_isotope_patterns(sum_formulae,adducts=(adduct,),isocalc_sig=isocalc_sig,isocalc_resolution=isocalc_resolution,charge=charge)
                if db_dump_folder != "":
                    cPickle.dump(mz_list_tmp,open(load_file,'w'))
            # add patterns to total list
            for sum_formula in mz_list_tmp:
                if not sum_formula in mz_list:
                    mz_list[sum_formula] = {}
                mzs, ints = mz_list_tmp[sum_formula][adduct]
                order = ints.argsort()[::-1]
                mz_list[sum_formula][adduct] = (mzs[order], ints[order])
        self.mz_list = mz_list
        print 'all isotope patterns generated and loaded'

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


class ReferencePipeline(Pipeline):
    def __init__(self, config):
        super(ReferencePipeline, self).__init__(config)

    def algorithm_name(self):
        return "reference"

    def compute_scores(self):
        from pyIMS.hdf5.inMemoryIMS_hdf5 import inMemoryIMS_hdf5
        IMS_dataset=inMemoryIMS_hdf5(self.data_file)
            
        for i, sum_formula in enumerate(self.sum_formulae):
            if i % 1000 == 0 and i > 0:
                print str(round(float(i) * 100.0 / len(self.sum_formulae), 2))+"%"
            for adduct in self.adducts:
                # 1. Generate ion images
                ion_datacube = IMS_dataset.get_ion_image(self.mz_list[sum_formula][adduct][0], self.ppm) #for each spectrum, sum the intensity of all peaks within tol of mz_list
                ion_datacube.xic = self.hot_spot_removal(ion_datacube.xic)

                # 2. Spatial Chaos 
                img = ion_datacube.xic_to_image(0)
                intensities = self.mz_list[sum_formula][adduct][1]

                self.score_chaos(sum_formula, adduct, img)
                # only compare pixels with values in the monoisotopic (otherwise high correlation for large empty areas)
                notnull_monoiso = ion_datacube.xic[0] > 0 
                #for xic in ion_datacube.xic:
                #    xic = xic[notnull_monoiso]

                # 3. Score correlation with monoiso
                self.score_corr(sum_formula, adduct, ion_datacube.xic, intensities[1:])
                # 4. Score isotope ratio
                self.score_ratio(sum_formula, adduct, ion_datacube.xic, intensities)

class NewPipeline(Pipeline):
    def __init__(self, config):
        super(NewPipeline, self).__init__(config)
        self.nrows = config['image_generation']['rows']
        self.ncols = config['image_generation']['columns']

    def algorithm_name(self):
        return "new"

    def compute_scores(self):
        def txt_to_spectrum(s):
            arr = s.strip().split("|")
            intensities = np.fromstring("0 " + arr[2], sep=' ')
            return (int(arr[0]), np.fromstring(arr[1], sep=' '), np.cumsum(intensities))
        print "reading spectra"
        spectra_str = open(self.data_file).readlines()
        spectra = map(txt_to_spectrum, spectra_str)
        
        chunk_size = 500
        n_chunks = len(self.sum_formulae) / chunk_size + 1

        for offset in xrange(0, len(self.sum_formulae), chunk_size):
            print("processing chunk #%d/%d" % (offset / chunk_size + 1, n_chunks))
            print "\tgetting ion images"
            formulae = self.sum_formulae[offset:offset + chunk_size]
            r = self.get_ion_images(spectra, formulae)
            print "\tcomputing scores"

            for xic, (sum_formula, adduct) in zip(r, product(formulae, self.adducts)):
                imgs = self.hot_spot_removal(xic)
                img = imgs[0].reshape((self.nrows, self.ncols)).copy()
                intensities = self.mz_list[sum_formula][adduct][1]
                self.process_query(sum_formula, adduct, imgs, img, intensities)

            del r
            import gc
            gc.collect()

    def process_spectra_multiple_queries(self, mol_mz_intervals, spectra):
        from numba import njit
        @njit
        def numba_multiple_queries(lower, upper, lperm, uperm, mzs, cumsum_int,
                                   result, pixel, tmp1, tmp2):
            ''' Equivalent to numpy_multiple_queries when lower and upper are sorted '''
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
        
        result = np.zeros((n, self.nrows * self.ncols))

        query_ids = np.zeros(n, dtype=np.int)
        intensities = np.zeros(n)
        tmp1 = np.zeros(n, dtype=np.int)
        tmp2 = np.zeros(n, dtype=np.int)

        for sp in spectra:
            pixel, mzs, cumsum_int = sp
            numba_multiple_queries(
                    lower, upper, lperm, uperm, mzs, cumsum_int,
                    result, pixel, tmp1, tmp2)

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
    elif config['method'] == 'new':
        pipeline = NewPipeline(config)
    else:
        print "method not recognized"
        sys.exit(1)
    pipeline.run()

