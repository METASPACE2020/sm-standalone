#include "cffi/common.hpp"
#include "imzb/reader.hpp"

extern "C" {
  DLL_PUBLIC void* imzb_reader_new(const char* filename) {
    return new imzb::ImzbReader(filename);
  }

  DLL_PUBLIC void imzb_reader_free(void* reader) {
    delete static_cast<imzb::ImzbReader*>(reader);
  }

  DLL_PUBLIC int imzb_reader_height(void* reader) {
    return static_cast<imzb::ImzbReader*>(reader)->height();
  }

  DLL_PUBLIC int imzb_reader_width(void* reader) {
    return static_cast<imzb::ImzbReader*>(reader)->width();
  }

  DLL_PUBLIC void imzb_reader_image(void* reader, double mz, double ppm, float* outbuf) {
    static_cast<imzb::ImzbReader*>(reader)->readImage(mz, ppm, outbuf);
  }
}
