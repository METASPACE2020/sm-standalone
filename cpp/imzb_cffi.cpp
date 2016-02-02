#include "imzb_reader.h"

extern "C" {
  void* imzb_reader_new(const char* filename) {
    return new imzb::ImzbReader(filename);
  }

  void imzb_reader_free(void* reader) {
    delete static_cast<imzb::ImzbReader*>(reader);
  }

  int imzb_reader_height(void* reader) {
    return static_cast<imzb::ImzbReader*>(reader)->height();
  }

  int imzb_reader_width(void* reader) {
    return static_cast<imzb::ImzbReader*>(reader)->width();
  }

  void imzb_reader_image(void* reader, double mz, double ppm, float* outbuf) {
    static_cast<imzb::ImzbReader*>(reader)->readImage(mz, ppm, outbuf);
  }
}
