# for PNG generation
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from pipelines import inMemoryIMS_low_mem

# webserver
import bottle

# isotope pattern generation
from pyMS.pyisocalc import pyisocalc

class ImageWebserver(bottle.Bottle):
    def __init__(self, *args, **kwargs):
        super(ImageWebserver, self).__init__(*args, **kwargs)

    def run(self, filename, **kwargs):
        print "loading data..."
        self.data = inMemoryIMS_low_mem(filename)
        print "running webserver..."
        super(ImageWebserver, self).run(**kwargs)

app = ImageWebserver()

@app.route('/')
def show_form():
    return '''
        <form action="/show_images" method="post">
            Formula: <input name="formula" type="text"><br/>
            Tolerance (ppm): <input name="tolerance" type="number" step="any" value="1.0"><br/>
            Pyisocalc cutoff: <input name="pyisocalc_cutoff" type="number" step="any" value="1e-5"><br/>
            <input value="Show images" type="submit"><br/>
            <input type="checkbox" checked="true" name="hs_removal">Remove hotspots
        </form>'''

import io
import os
import numpy as np

@app.route("/show_image/<mz>/<tol>")
def generate_image(mz, tol):
    mz, tol = float(mz), float(tol)
    img = app.data.get_ion_image(np.array([mz]), tol).xic_to_image(0)
    buf = io.BytesIO()
    if bottle.request.query.remove_hotspots:
        pc = np.percentile(img, 99)
        img[img > pc] = pc
    plt.imsave(buf, img, format='png')
    bottle.response.content_type = 'image/png'
    buf.seek(0, os.SEEK_END)
    bottle.response.content_length = buf.tell()
    buf.seek(0)
    return buf

from pyIMS.image_measures.isotope_pattern_match import isotope_pattern_match
from pyIMS.image_measures.isotope_image_correlation import isotope_image_correlation
from pyIMS.image_measures.level_sets_measure import measure_of_chaos

@app.route("/show_images", method="post")
def show_images():
    formula = bottle.request.forms.get('formula')
    tolerance = bottle.request.forms.get('tolerance')
    hs_removal = bottle.request.forms.get('hs_removal')
    cutoff = float(bottle.request.forms.get('pyisocalc_cutoff'))
    adducts = ['H', 'K', 'Na']
    isotope_patterns = {}
    for adduct in adducts:
        pattern = pyisocalc.isodist(formula + adduct, plot=False, sigma=0.01,
                                    charges=1, resolution=200000, cutoff=cutoff, do_centroid=True)

        mzs, intensities = pattern.get_spectrum(source='centroids')
        datacube = app.data.get_ion_image(np.array(mzs), float(tolerance))
        if hs_removal:
            for img in datacube.xic:
                pc = np.percentile(img, 99)
                img[img > pc] = pc

        chaos = 1 - measure_of_chaos(datacube.xic_to_image(0), 30, interp=False)[0]
        if chaos == 1: chaos = 0

        img_corr = 1.0 # return 1 if there's a single peak
        if len(intensities[1:]) > 1:
            img_corr = isotope_image_correlation(datacube.xic, weights=intensities[1:])

        iso_corr = isotope_pattern_match(datacube.xic, intensities)
        stats = {'measure of chaos': chaos,
                 'image correlation score': img_corr,
                 'isotope pattern score': iso_corr}
        
        isotope_patterns[adduct] = (mzs, intensities, stats)
    return bottle.template('show_images', hs_removal=hs_removal,
                           isotope_patterns=isotope_patterns, formula=formula, tol=tolerance)

import sys
app.run(sys.argv[1], port=8080)
