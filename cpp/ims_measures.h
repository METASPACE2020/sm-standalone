#pragma once

#include "ims_image.h"
#include "isotope_pattern.h"

#include <vector>
#include <valarray>
#include <cassert>
#include <iostream>

namespace ims {

  double isotopeImageCorrelation(const ims::ImageF* images, size_t n,
                                 const ms::IsotopePattern& pattern)
  {
    assert(n <= pattern.size());
    if (n < 2)
      return 0.0;

    auto n_empty_pixels = images[0].countEmptyPixels();
    auto cov = [&](const ims::ImageF& i1, const ims::ImageF& i2) -> float {
      // skip division as it cancels out in corr. coef. calculation
      return (i1.intensities() * i2.intensities()).sum() - n_empty_pixels;
    };

    auto principle_peak_norm = std::sqrt(cov(images[0], images[0]));
    constexpr double eps = 1e-6;
    if (std::fabs(principle_peak_norm) < eps)
      return 0.0;

    std::valarray<float> correlations(n - 1);
    for (size_t i = 1; i < n; ++i) {
      assert(images[i].shape() == images[0].shape());
      auto peak_norm = std::sqrt(cov(images[i], images[i]));
      if (std::fabs(peak_norm) < eps)
        correlations[i - 1] = 0.0;
      else
        correlations[i - 1] = cov(images[0], images[i]) / principle_peak_norm / peak_norm;
    }

    double weighted_sum = 0.0, total_weight = 0.0;
    for (size_t i = 0; i < correlations.size(); i++) {
      weighted_sum += correlations[i] * pattern.abundances[i + 1];
      total_weight += pattern.abundances[i + 1];
    }
    return weighted_sum / total_weight;
  }

  double isotopeImageCorrelation(const std::vector<ims::ImageF>& images,
                                 const ms::IsotopePattern& pattern)
  {
    return isotopeImageCorrelation(&images[0], images.size(), pattern);
  }

  double isotopePatternMatch(const ims::ImageF* images, size_t n,
                             const ms::IsotopePattern& pattern) {
    assert(n <= pattern.size());
    std::valarray<float> total_intensities(n), expected(n);
    auto n_empty_pixels = images[0].countEmptyPixels();
    double norm_squared1 = 0.0, norm_squared2 = 0.0;
    for (size_t i = 0; i < n; ++i) {
      total_intensities[i] = images[i].intensities().sum() + n_empty_pixels;
      norm_squared1 += std::pow(total_intensities[i], 2);
      expected[i] = pattern.abundances[i];
      norm_squared2 += std::pow(expected[i], 2);
    }
    if (std::fabs(norm_squared1) < 1e-6)
      return 0.0;
    total_intensities /= std::sqrt(norm_squared1);
    expected /= std::sqrt(norm_squared2);
    double score = 1.0;
    for (size_t i = 0; i < n; i++)
      score -= std::abs(expected[i] - total_intensities[i]) / n;
    return score;
  }
}
