import cffi
import numpy as np

from .utils import shared_lib, full_filename

ffi = cffi.FFI()
ffi.cdef(open(full_filename("ims.h")).read())
ims = ffi.dlopen(full_filename(shared_lib("ims_cffi")))

class ImzbReader(object):
    def __init__(self, filename):
        self.ptr = ffi.gc(ims.imzb_reader_new(filename),
                          ims.imzb_reader_free)

    def width(self):
        return ims.imzb_reader_width(self.ptr)

    def height(self):
        return ims.imzb_reader_height(self.ptr)

    def image(self, mz, ppm):
        buf = np.zeros(self.width() * self.height(), dtype=np.float32)
        ptr = ffi.cast("float*", buf.__array_interface__['data'][0])
        ims.imzb_reader_image(self.ptr, mz, ppm, ptr)
        return buf.reshape((self.height(), -1))

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import sys
    imzb = ImzbReader(sys.argv[1])
    img = imzb.image(float(sys.argv[2]), 5)
    plt.imshow(img)
    plt.show()
