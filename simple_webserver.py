# for PNG generation
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from colourmaps import viridis_colormap

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
            Resolution: <input name="resolution" type="number" step="any" value="200000"><br/>
            Points per FWHM: <input name="pts" type="number" step="1" value="10"></br>
            <input value="Show images" type="submit"><br/>
            <input type="checkbox" checked="true" name="hs_removal">Remove hotspots
        </form>'''

import io
import os
import numpy as np
from matplotlib.colors import Normalize

cmap = viridis_colormap()

@app.route("/show_image/<mz>/<tol>")
def generate_image(mz, tol):
    mz, tol = float(mz), float(tol)
    img = app.data.get_ion_image(np.array([mz]), tol).xic_to_image(0)
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

from pyIMS.image_measures.isotope_pattern_match import isotope_pattern_match
from pyIMS.image_measures.isotope_image_correlation import isotope_image_correlation
from pyIMS.image_measures.level_sets_measure import measure_of_chaos

@app.route("/show_images", method="post")
def show_images():
    formula = bottle.request.forms.get('formula')
    tolerance = float(bottle.request.forms.get('tolerance'))
    hs_removal = bottle.request.forms.get('hs_removal')
    resolution = float(bottle.request.forms.get('resolution'))
    pts = int(bottle.request.forms.get('pts'))
    cutoff = float(bottle.request.forms.get('pyisocalc_cutoff'))
    adducts = ['H', 'K', 'Na']
    isotope_patterns = {}
    for adduct in adducts:
        sf = pyisocalc.SumFormulaParser.parse_string(formula + adduct)
        raw_pattern = pyisocalc.isodist(sf, cutoff)
        fwhm = raw_pattern.get_spectrum()[0][0] / resolution
        pattern = pyisocalc.apply_gaussian(raw_pattern, fwhm, pts, exact=True)

        mzs, intensities = pattern.get_spectrum(source='centroids')
        datacube = app.data.get_ion_image(np.array(mzs), tolerance)
        if hs_removal:
            for img in datacube.xic:
                pc = np.percentile(img, 99)
                img[img > pc] = pc

        chaos = 1 - measure_of_chaos(datacube.xic_to_image(0), 30, interp=False)[0]
        if chaos == 1: chaos = 0

        iso_corr = isotope_pattern_match(datacube.xic, intensities)

        img_corr = 1.0 # return 1 if there's a single peak
        if len(intensities[1:]) > 1:
            img_corr = isotope_image_correlation(datacube.xic, weights=intensities[1:])

        stats = {'measure of chaos': chaos,
                 'image correlation score': img_corr,
                 'isotope pattern score': iso_corr}
        
        isotope_patterns[adduct] = (mzs, intensities, stats)
    return bottle.template('show_images', hs_removal=hs_removal,
                           isotope_patterns=isotope_patterns, formula=formula, tol=tolerance)

import sys
app.run(sys.argv[1], port=8080)
