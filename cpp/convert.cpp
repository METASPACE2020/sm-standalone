#include "imzml_reader.h"
#include "imzb_reader.h"
#include "imzb_writer.h"

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <memory>
#include <queue>
#include <cstdio>

class Sorter {
  std::string fn_;
  std::vector<ims::Peak> buffer_;
  std::vector<std::string> tmp_filenames_;
  size_t filled_;
  bool closed_;

  void dump() {
    std::sort(buffer_.begin(), buffer_.begin() + filled_,
        [](const ims::Peak& a, const ims::Peak& b) { return a.mz < b.mz; }
    );
    std::stringstream tmp_fn;
    tmp_fn << fn_ << "." << tmp_filenames_.size();
    tmp_filenames_.push_back(tmp_fn.str());
    imzb::ImzbWriter writer(tmp_filenames_.back());
    for (size_t i = 0; i < filled_; ++i)
      writer.writePeak(buffer_[i]);
    writer.close();

    filled_ = 0;
  }

  struct PeakAndFile {
    ims::Peak peak;
    size_t file_index;

    bool operator<(const PeakAndFile& other) const {
      return peak.mz > other.peak.mz;
    }
  };

  void merge() {
    std::priority_queue<PeakAndFile> queue;

    std::vector<std::shared_ptr<imzb::ImzbReader>> readers;
    imzb::ImzbWriter writer(fn_);
    writer.setMask(mask_);

    ims::Peak peak;
    for (const auto& tmp_fn: tmp_filenames_) {
      readers.push_back(std::make_shared<imzb::ImzbReader>(tmp_fn));
      if (readers.back()->readNext(peak))
        queue.push(PeakAndFile{peak, readers.size() - 1});
    }

    while (!queue.empty()) {
      auto item = queue.top();
      writer.writePeak(item.peak);
      queue.pop();

      if (readers[item.file_index]->readNext(peak))
        queue.push(PeakAndFile{peak, item.file_index});
    }

    writer.close();

    for (const auto& fn: tmp_filenames_) {
      std::remove(fn.c_str());
      auto idx_fn = fn + ".idx";
      std::remove(idx_fn.c_str());
    }
  }

  imzb::Mask mask_;
public:
  Sorter(const std::string& filename, size_t buffer_size=10000000) :
    fn_(filename),buffer_(buffer_size), filled_(0), closed_(false)
  {
  }

  void setMask(const imzb::Mask& mask) {
    mask_ = mask;
  }

  void addPeak(const ims::Peak& peak) {
    if (filled_ == buffer_.size())
      dump();
    buffer_[filled_++] = peak;
  }

  void close() {
    dump();
    merge();
    closed_ = true;
  }

  ~Sorter() {
    if (!closed_)
      close();
  }
};

int main(int argc, char** argv) {
  if (argc == 1) {
    std::cout << "usage: ./convert <file.imzML> <out.imzb>";
    return 0;
  }

  imzml::ImzmlReader imzml(argv[1]);
  Sorter sorter(argv[2]);
  imzb::Mask mask{imzml.height(), imzml.width()};
  
  std::vector<std::string> tmp_filenames_;
  ims::Spectrum sp;
  while (imzml.readNextSpectrum(sp)) {
    for (size_t i = 0; i < sp.mzs.size(); ++i) {
      // skip invalid (NaN) and zero peaks
      if (sp.mzs[i] > 0 && sp.intensities[i] > 0)
        sorter.addPeak(ims::Peak{sp.coords, sp.mzs[i], sp.intensities[i]});
    }
    mask.set(sp.coords.x, sp.coords.y);
  }
  sorter.setMask(mask);
  sorter.close();
}
