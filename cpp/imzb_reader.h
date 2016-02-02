#pragma once

#include "ims.h"
#include "blosc.h"
#include "imzb_index.h"
#include "imzb_writer.h"

#include <string>
#include <fstream>
#include <vector>
#include <cstdint>
#include <cassert>
#include <memory>
#include <algorithm>

namespace imzb {

class ImzbReader {
  std::string fn_;
  std::ifstream in_;

  IndexPtr index_;
  size_t block_idx_;

  std::vector<char> buffer_;
  std::vector<ims::Peak> peaks_;
  size_t n_peaks_;
  size_t pos_;

  size_t decompressBlock(size_t block_idx,
                         std::ifstream& in,
                         std::vector<char>& inbuf,
                         std::vector<ims::Peak>& outbuf) const
  {
    uint64_t start = index_->offsets[block_idx];
    uint64_t end = index_->offsets[block_idx + 1];
    inbuf.resize(end - start);
    in.seekg(start);
    in.read(&inbuf[0], end - start);
    return blosc_decompress_ctx(&inbuf[0], &outbuf[0],
        outbuf.size() * sizeof(ims::Peak), 1) / sizeof(ims::Peak);
  }

  bool readNextBlock() {
    ++block_idx_;
    if (block_idx_ == index_->mzs.size()) {
      n_peaks_ = 0;
      return false;
    }

    n_peaks_ = decompressBlock(block_idx_, in_, buffer_, peaks_);
    pos_ = 0;
    return true;
  }

  bool empty_;

public:
  ImzbReader(const std::string& filename) :
    fn_(filename),
    in_(filename), block_idx_(0), peaks_(imzb::ImzbWriter::BLOCK_SIZE),
    n_peaks_(0), pos_(0), empty_(false)
  {
    std::ifstream in_idx(filename + ".idx");

    index_ = std::make_shared<Index>();
    index_->read(in_idx);

    in_.seekg(0, in_.end);
    index_->offsets.push_back(in_.tellg());
    in_.seekg(0, in_.beg);
  }

  bool readNext(ims::Peak& peak) {
    if (empty_)
      return false;

    if (pos_ >= n_peaks_) {
      if (!readNextBlock()) {
        empty_ = true;
        return false;
      }
    }

    peak = peaks_[pos_];
    ++pos_;
    return true;
  }

  void reset() {
    in_.seekg(0, in_.beg);
    n_peaks_ = pos_ = block_idx_ = 0;
    empty_ = false;
  }

  std::vector<ims::Peak> slice(double min_mz, double max_mz) const {
    assert(min_mz < max_mz);
    std::vector<char> inbuf;
    std::vector<ims::Peak> result, outbuf{imzb::ImzbWriter::BLOCK_SIZE};
    size_t start_block = index_->startBlock(min_mz);
    size_t end_block = index_->endBlock(max_mz);
    std::ifstream in{fn_};

    for (size_t i = start_block; i < end_block; ++i) {
      size_t n = decompressBlock(i, in, inbuf, outbuf);
      auto beg = outbuf.cbegin(), end = outbuf.cbegin() + n;
      if (outbuf.front().mz < min_mz)
        beg = std::lower_bound(beg, end, min_mz,
            [](const ims::Peak& p, double mz) { return p.mz < mz; });
      if (outbuf.back().mz > max_mz)
        end = std::upper_bound(beg, end, max_mz,
            [](double mz, const ims::Peak& p) { return mz < p.mz; });
      result.insert(result.end(), beg, end);
    }
    return result;
  }

  uint32_t height() const { return index_->header.mask.height; }
  uint32_t width() const { return index_->header.mask.width; }
};

}
