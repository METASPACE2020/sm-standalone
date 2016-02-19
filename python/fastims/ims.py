import cffi

from .utils import shared_lib, full_filename

ffi = cffi.FFI()
ffi.cdef(open(full_filename("ims.h")).read())
ims = ffi.dlopen(full_filename(shared_lib("ims_cffi")))

class IsotopePattern(object):
    def __init__(self, sum_formula, threshold=1e-4, fft_threshold=1e-8):
        self.ptr = ffi.gc(ims.isotope_pattern_new_from_sf(sum_formula, threshold, fft_threshold),
                          ims.isotope_pattern_free)

    def centroids(self, resolution, min_abundance=1e-4, points_per_fwhm=25):
        obj = object.__new__(IsotopePattern)
        obj.ptr = ffi.gc(ims.isotope_pattern_centroids(self.ptr, resolution,
                                                       min_abundance, points_per_fwhm),
                         ims.isotope_pattern_free)
        return obj

    def size(self):
        return ims.isotope_pattern_size(self.ptr)

    def masses(self):
        buf = ffi.new("double[]", self.size())
        ims.isotope_pattern_masses(self.ptr, buf)
        return list(buf)

    def abundances(self):
        buf = ffi.new("double[]", self.size())
        ims.isotope_pattern_abundances(self.ptr, buf)
        return list(buf)

    def addCharge(self, charge):
        ims.isotope_pattern_add_charge(self.ptr, charge)

    def trim(self, n_peaks):
        ims.isotope_pattern_trim(self.ptr, n_peaks)

    def envelope(self, resolution):
        def envelopeFunc(mz):
            return ims.isotope_pattern_envelope(self.ptr, resolution, mz)
        return envelopeFunc

if __name__ == '__main__':
    import sys
    pattern = IsotopePattern(sys.argv[1])
    pattern = pattern.centroids(int(sys.argv[2]))
    pattern.addCharge(1)
    pattern.trim(5)
    print(zip(pattern.masses(), pattern.abundances()))
