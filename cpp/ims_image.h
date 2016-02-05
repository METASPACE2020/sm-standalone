#pragma once

#include <valarray>
#include <algorithm>

namespace ims {

  inline size_t pixelIndex(size_t x, size_t y, size_t width) {
    return x * width + y;
  }

  template <typename T>
    class Image {
      std::valarray<T> intensities_;
      size_t height_, width_;
      public:
      Image(size_t height, size_t width) :
        intensities_(height * width), height_(height), width_(width)
      {
      }

      size_t height() const { return height_; }
      size_t width() const { return width_; }

      const T& intensity(size_t x, size_t y) const {
        return intensities_[pixelIndex(x, y, width_)];
      }

      T& intensity(size_t x, size_t y) {
        return intensities_[pixelIndex(x, y, width_)];
      }

      T* rawPtr() { return &intensities_[0]; }

      size_t countEmptyPixels() const {
        size_t n = 0;
        for (auto& x: intensities_) if (x < 0) ++n;
        return n;
      }

      const std::valarray<T>& intensities() const { return intensities_; }

      std::pair<size_t, size_t> shape() const {
        return std::make_pair(height(), width());
      }
    };

  typedef Image<float> ImageF;
  typedef Image<double> ImageD;
}
