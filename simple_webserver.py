# for PNG generation
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import re

from colourmaps import viridis_colormap

from pipelines import inMemoryIMS_low_mem

# webserver
import bottle

# isotope pattern generation
from pyMS.pyisocalc import pyisocalc

from imzml_to_pytables import read_mz_image, read_ion_datacube

class ImageWebserver(bottle.Bottle):
    def __init__(self, *args, **kwargs):
        super(ImageWebserver, self).__init__(*args, **kwargs)

    def run(self, filename, **kwargs):
        print "loading data..."
        self.load_data(filename)
        print "running webserver..."
        super(ImageWebserver, self).run(**kwargs)

    def load_data(self, filename):
        if filename.endswith(".imzML") or filename.endswith(".hdf5"):
            self.data = inMemoryIMS_low_mem(filename)
            self.in_memory = True
        elif filename.endswith(".imh5"):
            self.filename = filename
            self.in_memory = False

    def get_ion_image(self, mz, tol):
        if self.in_memory is True:
            return self.data.get_ion_image(np.array([mz]), tol).xic_to_image(0)
        else:
            return read_mz_image(self.filename, mz, tol, hotspot_removal=False)

    def get_datacube(self, mzs, tol):
        if self.in_memory is True:
            return self.data.get_ion_image(np.array(mzs), tol)
        else:
            return read_ion_datacube(self.filename, mzs, tol)

app = ImageWebserver()

@app.route('/', method='GET')
def show_form():
    return bottle.template('show_images', hs_removal=True,
                           isotope_patterns={}, formula="", selected_adduct='H', pretty_formula="", tol=5,
                           resolution=100000)

@app.route('/', method='POST')
def show_images():
    formula = bottle.request.forms.get('formula')
    tolerance = float(bottle.request.forms.get('tolerance'))
    hs_removal = bottle.request.forms.get('hs_removal')
    resolution = float(bottle.request.forms.get('resolution'))
    pts = int(bottle.request.forms.get('pts'))
    cutoff = float(bottle.request.forms.get('pyisocalc_cutoff'))
    adduct = bottle.request.forms.get('adduct')
    return _generate_page(formula, adduct, tolerance, hs_removal, resolution, pts, cutoff)

import io
import os
import numpy as np
from matplotlib.colors import Normalize

cmap = viridis_colormap()

@app.route("/show_image/<mz>/<tol>")
def generate_image(mz, tol):
    mz, tol = float(mz), float(tol)
    img = app.get_ion_image(mz, tol)
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

@app.route("/correlation_plot/<formula>/<adduct>/<mzs>/<intensities>/<tol>")
def generate_correlation_plot(formula, adduct, mzs, intensities, tol):
    mzs = np.array(map(float, mzs.split(",")))
    intensities = np.array(map(float, intensities.split(",")))
    order = intensities.argsort()[::-1]
    mzs = mzs[order]
    intensities = intensities[order]
    tol = float(tol)

    datacube = app.get_datacube(np.array(mzs), tol)
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
    plt.imshow(deviations, vmin=0, vmax=1, cmap="gnuplot")
    plt.axis('off')

    # http://stackoverflow.com/questions/26034777/matplotlib-same-height-for-colorbar-as-for-plot
    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes("right", size="5%", pad="3%")

    cbar = plt.colorbar(cax = cax1)

    markersize = min(20, (10000.0 / (1 + np.sum(images[1] > 0))) ** 0.5)

    plt.subplot(1, 2, 2)

    plt.xlabel("sqrt( principal peak intensities )")
    plt.ylabel("sqrt( other peak intensities )")
    plt.title(formula + " + " + adduct + " (m/z={:.4f})".format(mzs[0]))

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

@app.route("/show/<formula>")
def show_images_get(formula):
    tolerance = float(bottle.request.params.get('ppm', 5.0))
    resolution = float(bottle.request.params.get('resolution', 1e5))
    return _generate_page(formula, adduct, tolerance, True, resolution, 10, 0.001)

def _generate_page(formula, selected_adduct, tolerance, hs_removal, resolution, pts, cutoff):
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

        datacube = app.get_datacube(mzs, tolerance)
        if hs_removal:
            for img in datacube.xic:
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
                           resolution=resolution, tol=tolerance)

import sys
app.run(sys.argv[1], port=8080)
