#include "ims/image_measures.hpp"

template <typename T>
double measure_of_chaos(T* image, int width, int height, int n_levels) {
  ims::ImageF img(height, width);
  auto p = img.rawPtr();
  for (size_t i = 0, n = height * width; i < n; i++)
    p[i] = image[i];
  return measureOfChaos(img, n_levels);
}

extern "C" {
  double measure_of_chaos_f(float* image, int width, int height, int n_levels) {
    return measure_of_chaos<float>(image, width, height, n_levels);
  }

  double measure_of_chaos_d(double* image, int width, int height, int n_levels) {
    return measure_of_chaos<double>(image, width, height, n_levels);
  }
}
