# for PNG generation
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import re

from colourmaps import viridis_colormap

from pipelines import inMemoryIMS_low_mem

from collections import OrderedDict

# webserver
import bottle

# isotope pattern generation
from pyMS.pyisocalc import pyisocalc

from pyIMS.ion_datacube import ion_datacube

import cffi
ffi = cffi.FFI()
ffi.cdef(open("/home/lomereiter/github/ims-cpp/cffi/ims.h").read())
imzb = ffi.dlopen("/home/lomereiter/github/ims-cpp/build/libimzb_cffi.so")

class ImzbReader(object):
    def __init__(self, filename):
        self._reader = ffi.gc(imzb.imzb_reader_new(filename),
                              imzb.imzb_reader_free)

    def height(self):
        return imzb.imzb_reader_height(self._reader)

    def width(self):
        return imzb.imzb_reader_width(self._reader)

    def get_mz_image(self, mz, ppm):
        data = np.zeros(self.height() * self.width(), dtype=np.float32)
        imzb.imzb_reader_image(self._reader, ffi.cast("double", mz), ffi.cast("double", ppm),
                               ffi.from_buffer(data))
        return data.reshape((self.height(), self.width()))

    def get_datacube(self, mzs, ppm):
        cube = ion_datacube()
        cube.xic = []
        cube.nRows = self.height()
        cube.nColumns = self.width()
        cube.pixel_indices = None

        for mz in mzs:
            img = self.get_mz_image(mz, ppm)
            if cube.pixel_indices is None:
                cube.pixel_indices = np.where(img.ravel() >= 0)[0]
            img = img.ravel()[cube.pixel_indices]
            img[img < 0] = 0.0
            cube.xic.append(img)
        return cube

class ImageWebserver(bottle.Bottle):
    def __init__(self, *args, **kwargs):
        super(ImageWebserver, self).__init__(*args, **kwargs)

    def run(self, filenames, **kwargs):
        self.load_data(filenames)
        print "running webserver..."
        super(ImageWebserver, self).run(**kwargs)

    def load_data(self, filenames):
        def prettify_fn(fn):
            import os
            return os.path.splitext(os.path.basename(fn))[0]

        if len(filenames) == 0:
            print "usage: python simple_webserver.py <file.imzML>"
            print "       python simple_webserver.py <file.hdf5>"
            print "       python simple_webserver.py <file1.imzb> [<file2.imzb> ...]"
            sys.exit(0)
        if len(filenames) > 1 and not all(fn.endswith(".imzb") for fn in filenames):
            print "multiple-file mode is supported only for .imzb files"
            sys.exit(2)
        if len(filenames) == 1 and not filenames[0].endswith(".imzb"):
            filename = filenames[0]
            if filename.endswith(".imzML") or filename.endswith(".hdf5"):
                print "loading data..."
                self.data = inMemoryIMS_low_mem(filename)
                self.in_memory = True
                self.paths = { prettify_fn(filename) : filename }
            else:
                print "unsupported format"
                sys.exit(3)
        else:
            self.paths = OrderedDict()
            for fn in filenames:
                if os.path.exists(fn):
                    self.paths[prettify_fn(fn)] = ImzbReader(fn)
                else:
                    print "WARNING: file " + fn + " doesn't exist, skipping"
            self.in_memory = False

    def get_ion_image(self, dataset, mz, tol):
        if self.in_memory is True:
            return self.data.get_ion_image(np.array([mz]), tol).xic_to_image(0)
        else:
            return self.paths[dataset].get_mz_image(mz, tol)

    def get_datacube(self, dataset, mzs, tol):
        if self.in_memory is True:
            return self.data.get_ion_image(np.array(mzs), tol)
        else:
            return self.paths[dataset].get_datacube(mzs, tol)

app = ImageWebserver()

@app.route('/', method='GET')
def show_form():
    return bottle.template('show_images', hs_removal=True, selected_dataset=app.paths.iterkeys().next(),
                           isotope_patterns={}, formula="", selected_adduct='H', pretty_formula="", tol=5,
                           resolution=100000, datasets=app.paths.keys())

import io
import os
import numpy as np
from matplotlib.colors import Normalize

cmap = viridis_colormap()

@app.route("/show_image/<dataset>/<mz>/<tol>")
def generate_image(dataset, mz, tol):
    mz, tol = float(mz), float(tol)
    img = app.get_ion_image(dataset, mz, tol)
    if img.shape[0] > img.shape[1]:
        img = img.T
    buf = io.BytesIO()
    mask = img >= 0
    if bottle.request.query.remove_hotspots:
        pc = np.percentile(img[mask], 99)
        img[img > pc] = pc
    values = img[mask]
    norm = Normalize(vmin=np.min(values), vmax=np.max(values))
    colorized_img = np.zeros((img.shape[0], img.shape[1], 4))
    colorized_img[mask] = cmap(norm(values))
    # set alpha channel to 0 for pixels with no data
    colorized_img[img < 0, -1] = 0
    plt.imsave(buf, colorized_img, format='png')
    bottle.response.content_type = 'image/png'
    buf.seek(0, os.SEEK_END)
    bottle.response.content_length = buf.tell()
    buf.seek(0)
    return buf

@app.route("/correlation_plot/<dataset>/<formula>/<adduct>/<mzs>/<intensities>/<tol>")
def generate_correlation_plot(dataset, formula, adduct, mzs, intensities, tol):
    mzs = np.array(map(float, mzs.split(",")))
    intensities = np.array(map(float, intensities.split(",")))
    order = intensities.argsort()[::-1]
    mzs = mzs[order]
    intensities = intensities[order]
    tol = float(tol)

    datacube = app.get_datacube(dataset, np.array(mzs), tol)
    images = datacube.xic

    buf = io.BytesIO()
    transform = np.sqrt
    base_intensities = images[0]

    plt.figure(figsize=(16, 8))
    ax1 = plt.subplot(1, 2, 1)
    plt.title("per-pixel isotope pattern agreement (higher is better)")
    n = min(len(datacube.xic), len(intensities))
    full_images = np.array([transform(datacube.xic_to_image(i)) for i in xrange(n)])
    full_images /= np.linalg.norm(full_images, ord=2, axis=0)
    normalized_ints = transform(intensities[:n])
    normalized_ints /= np.linalg.norm(normalized_ints)
    #correlations = np.einsum("ijk,i", full_images, normalized_ints)
    #plt.imshow(correlations, vmin=0, vmax=1)
    deviations = 1 - np.amax(np.abs(np.transpose(full_images, (1, 2, 0)) - normalized_ints), axis=2)
    if deviations.shape[0] > deviations.shape[1]:
        deviations = deviations.T
    plt.imshow(deviations, vmin=0, vmax=1, cmap="gnuplot", interpolation='none')
    plt.axis('off')

    # http://stackoverflow.com/questions/26034777/matplotlib-same-height-for-colorbar-as-for-plot
    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes("right", size="5%", pad="3%")

    cbar = plt.colorbar(cax = cax1)

    markersize = min(20, (10000.0 / (1 + np.sum(images[1] > 0))) ** 0.5)

    plt.subplot(1, 2, 2)

    plt.xlabel("sqrt( principal peak intensities )")
    plt.ylabel("sqrt( other peak intensities )")
    plt.title(formula + " + " + adduct + " (m/z={:.4f})".format(mzs[0]) +\
              "\n(lines are based on the predicted isotope pattern)")

    colors = ['blue', 'red', 'green', 'purple', 'black']

    for i in xrange(1, min(5, len(images))):
        ratio = intensities[i] / intensities[0]
        observed = images[i]

        mask = base_intensities > 0
        label = "m/z={0:.4f} {1:.1%}".format(mzs[i], intensities[i] / 100.0)
        plt.plot(transform(base_intensities[mask]), transform(observed[mask]), '.', markersize=markersize,
                 color = colors[i-1], label=label)

        xs = transform(base_intensities[mask])
        ys = transform(base_intensities[mask] * ratio)
        order = xs.argsort()
        plt.plot(xs[order], ys[order], color=colors[i-1], linewidth=0.5)
    lgnd = plt.legend(loc='upper left', numpoints=10)
    # http://stackoverflow.com/questions/24706125/setting-a-fixed-size-for-points-in-legend
    for handle in lgnd.legendHandles:
        handle._legmarker.set_markersize(6)
    plt.tight_layout(w_pad=5.0)
    plt.savefig(buf)
    plt.close()

    bottle.response.content_type = 'image/png'
    buf.seek(0, os.SEEK_END)
    bottle.response.content_length = buf.tell()
    buf.seek(0)
    return buf

from pyIMS.image_measures.isotope_pattern_match import isotope_pattern_match
from pyIMS.image_measures.isotope_image_correlation import isotope_image_correlation
from pyIMS.image_measures.level_sets_measure import measure_of_chaos

@app.route("/show")
def show_images_get():
    dataset = bottle.request.params.get('dataset', app.paths.iterkeys().next())
    formula = bottle.request.params.get('formula', '')
    tolerance = float(bottle.request.params.get('tolerance', 5.0))
    resolution = float(bottle.request.params.get('resolution', 1e5))
    selected_adduct = bottle.request.params.get('adduct', 'H')
    hs_removal = bottle.request.GET.get('hs_removal', False)
    if hs_removal == 'on':
        hs_removal = True
    pts = float(bottle.request.params.get('pts', 10))
    cutoff = float(bottle.request.params.get('pyisocalc_cutoff', 1e-3))

    adducts = ['H', 'K', 'Na']
    isotope_patterns = {}
    for adduct in adducts:
        sf = pyisocalc.SumFormulaParser.parse_string(formula + adduct)
        raw_pattern = pyisocalc.isodist(sf, cutoff)
        fwhm = raw_pattern.get_spectrum()[0][0] / resolution
        pattern = pyisocalc.apply_gaussian(raw_pattern, fwhm, pts, exact=True)

        mzs, intensities = map(np.array, pattern.get_spectrum(source='centroids'))
        k = 5
        if len(mzs) > k:
            order = intensities.argsort()[::-1]
            mzs = mzs[order][:k]
            intensities = intensities[order][:k]
            order = mzs.argsort()
            mzs = mzs[order]
            intensities = intensities[order]

        datacube = app.get_datacube(dataset, mzs, tolerance)
        if hs_removal:
            for img in datacube.xic:
                if len(img) > 0:
                    pc = np.percentile(img, 99)
                    img[img > pc] = pc

        chaos = measure_of_chaos(datacube.xic_to_image(0), 30, interp=False)[0]

        iso_corr = isotope_pattern_match(datacube.xic, intensities)

        img_corr = 1.0 # return 1 if there's a single peak
        if len(intensities[1:]) > 1:
            img_corr = isotope_image_correlation(datacube.xic, weights=intensities[1:])

        stats = {'measure of chaos': chaos,
                 'image correlation score': img_corr,
                 'isotope pattern score': iso_corr}

        isotope_patterns[adduct] = (mzs, intensities, stats)
    return bottle.template('show_images', hs_removal=hs_removal,
                           isotope_patterns=isotope_patterns, formula=formula, selected_adduct=selected_adduct,
                           pretty_formula=re.sub(r"(\d+)", r"<sub>\1</sub>", formula),
                           resolution=resolution, tol=tolerance, datasets=app.paths.keys(),
                           selected_dataset=dataset)

import sys
app.run(sys.argv[1:], port=8080)
