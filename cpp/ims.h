#pragma once

#include "ims_image.h"
#include "ims_measures.h"

#include <vector>
#include <cstdint>

namespace ims {

struct Position {
  uint32_t x;
  uint32_t y;
  uint32_t z;
};

struct Spectrum {
  std::vector<double> mzs;
  std::vector<float> intensities;

  Position coords;
};

struct Peak {
  Position coords;
  double mz;
  float intensity;
};

}
