#pragma once

#include "ims.h"
#include "blosc.h"

#include <string>
#include <fstream>
#include <vector>
#include <cassert>

namespace imzb {

class ImzbWriter {
public:
  static const uint32_t BLOCK_SIZE = 4096;

private:
  std::ofstream out_, index_;

  std::vector<ims::Peak> in_buf_;
  std::vector<char> out_buf_;
  uint32_t filled_;
  uint64_t bytes_written_;

  void dump() {
    int n = blosc_compress_ctx(5, 1, sizeof(ims::Peak),
                               filled_ * sizeof(ims::Peak), &in_buf_[0],
                               &out_buf_[0], out_buf_.size(),
                               "blosclz", 0, 1);
    assert(n > 0);
    filled_ = 0;
    out_.write(&out_buf_[0], n);
    binary_write(index_, in_buf_[0].mz);
    binary_write(index_, bytes_written_);
    bytes_written_ += n;
  }

  template <typename T>
  void binary_write(std::ofstream& stream, T& value) {
    stream.write(reinterpret_cast<char*>(&value), sizeof(value));
  }

public:
  ImzbWriter(const std::string& filename) :
    out_(filename), index_(filename + ".idx"),
    in_buf_(BLOCK_SIZE),
    out_buf_(BLOCK_SIZE * sizeof(ims::Peak) + BLOSC_MAX_OVERHEAD * 2),
    filled_(0), bytes_written_(0)
  {
  }

  void writePeak(const ims::Peak& peak) {
    if (filled_ == in_buf_.size())
      dump();
    in_buf_[filled_++] = peak;
  }

  void close() {
    dump();
    out_.close();
    index_.close();
  }

  ~ImzbWriter() {
    if (out_.is_open())
      close();
  }
};

}
