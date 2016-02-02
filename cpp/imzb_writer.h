#pragma once

#include "ims.h"
#include "imzb_fileutils.h"
#include "imzb_fileheader.h"
#include "imzb_index.h"
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
  std::ofstream out_;
  imzb::Index index_;

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
    index_.mzs.push_back(in_buf_[0].mz);
    index_.offsets.push_back(bytes_written_);
    bytes_written_ += n;
  }

  std::string filename_;

public:
  ImzbWriter(const std::string& filename) :
    out_(filename),
    in_buf_(BLOCK_SIZE),
    out_buf_(BLOCK_SIZE * sizeof(ims::Peak) + BLOSC_MAX_OVERHEAD * 2),
    filled_(0), bytes_written_(0), filename_(filename)
  {
  }

  void setMask(const imzb::Mask& mask) {
    index_.header.mask = mask;
  }

  void writePeak(const ims::Peak& peak) {
    if (filled_ == in_buf_.size())
      dump();
    in_buf_[filled_++] = peak;
  }

  void close() {
    dump();
    out_.close();

    std::ofstream out_idx(filename_ + ".idx");
    index_.write(out_idx);
    out_idx.close();
  }

  ~ImzbWriter() {
    if (out_.is_open())
      close();
  }
};

}
