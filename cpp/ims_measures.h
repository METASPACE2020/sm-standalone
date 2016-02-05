#pragma once

#include "ims_image.h"
#include "isotope_pattern.h"

#include <vector>
#include <valarray>
#include <cassert>
#include <iostream>

namespace ims {

  double isotopeImageCorrelation(const ims::ImageF* images, size_t n,
                                 const ms::IsotopePattern& pattern);
  
  double isotopeImageCorrelation(const std::vector<ims::ImageF>& images,
                                 const ms::IsotopePattern& pattern);

  double isotopePatternMatch(const ims::ImageF* images, size_t n,
                             const ms::IsotopePattern& pattern);

  double isotopePatternMatch(const std::vector<ims::ImageF>& images,
                             const ms::IsotopePattern& pattern);
}
