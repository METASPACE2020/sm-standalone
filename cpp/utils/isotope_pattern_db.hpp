#pragma once

#include "ms/isocalc.hpp"
#include "msgpack.hpp"

#include <map>
#include <string>
#include <vector>
#include <mutex>
#include <fstream>
#include <cassert>
#include <iostream>

namespace utils {

class InstrumentProfile {
public:
  virtual double resolutionAt(double mz) const = 0;
  virtual ~InstrumentProfile() {}
};

class OrbitrapProfile : public InstrumentProfile {
  double resolution_200_;
public:
  OrbitrapProfile(double resolution_at_200) :
    resolution_200_(resolution_at_200)
  {}

  double resolutionAt(double mz) const {
    return resolution_200_ * (200.0 / mz);
  }

  OrbitrapProfile() {
  }
};

class IsotopePatternDB {
  typedef std::pair<std::string, std::string> SFAdductPair;
  std::vector<SFAdductPair> pairs_;

  // use a simple structure that doesn't require special msgpack methods
  std::map<SFAdductPair, std::map<std::string, std::vector<double>>> patterns_;

  static std::vector<SFAdductPair> makePairs(
    const std::vector<std::string>& sum_formulas,
    const std::vector<std::string>& adducts)
  {
    std::vector<SFAdductPair> pairs;
    for (auto& f: sum_formulas)
      for (auto& a: adducts)
        pairs.push_back(std::make_pair(f, a));
    return pairs;
  }

public:
  IsotopePatternDB(const std::vector<std::string>& sum_formulas,
                   const std::vector<std::string>& adducts)
    : IsotopePatternDB(makePairs(sum_formulas, adducts))
  {
  }

  IsotopePatternDB(const std::vector<std::pair<std::string, std::string>>& sf_adduct_pairs) :
    pairs_(sf_adduct_pairs)
  {
  }

  IsotopePatternDB(const std::string& dump_filename) {
    load(dump_filename);
  }

  void computeIsotopePatterns(const utils::InstrumentProfile& instrument,
                              size_t max_peaks=5)
  {
    std::mutex map_mutex;

#pragma omp parallel for
    for (size_t i = 0; i < pairs_.size(); i++) {
      const auto& f = pairs_[i].first;
      auto adduct = pairs_[i].second;
      int charge = adduct[0] == '-' ? -1 : +1;
      if (adduct[0] != '-' && adduct[0] != '+')
        adduct = "+" + adduct;
      auto full_formula = f + adduct;
      auto res = instrument.resolutionAt(ms::monoisotopicMass(full_formula));
      auto pattern = ms::computeIsotopePattern(full_formula)\
                         .centroids(res).charged(charge).trimmed(max_peaks);
      auto key = std::make_pair(f, adduct);
      std::lock_guard<std::mutex> lock(map_mutex);
      patterns_[key]["mzs"] = pattern.masses;
      patterns_[key]["abundances"] = pattern.abundances;
    }
  }

  void save(const std::string& output_filename) {
    // TODO: write instrument profile parameters?
    msgpack::sbuffer sbuf;
    msgpack::pack(sbuf, patterns_);
    std::ofstream db(output_filename, std::ios::binary);
    db.write(sbuf.data(), sbuf.size());
  }

  void load(const std::string& input_filename) {
    // load precomputed isotope patterns from msgpack file
    std::ifstream in(input_filename, std::ios::binary | std::ios::ate);
    size_t bufsize = in.tellg();
    msgpack::sbuffer sbuf(bufsize);
    in.seekg(0, std::ios::beg);
    in.read(sbuf.data(), bufsize);
    in.close();

    msgpack::unpacked unpacked;
    msgpack::unpack(&unpacked, sbuf.data(), bufsize);

    msgpack::object obj = unpacked.get();
    obj.convert(patterns_);
    for (auto& item: patterns_) pairs_.push_back(item.first);
  }

  ms::IsotopePattern operator()(const std::string& formula, const std::string& adduct) const {
    auto& entry = patterns_.at(std::make_pair(formula, adduct));
    return ms::IsotopePattern{entry.at("mzs"), entry.at("abundances")};
  }

  size_t size() const {
    assert(patterns_.size() == pairs_.size());
    return patterns_.size();
  }

  const std::vector<std::pair<std::string, std::string>>& keys() const {
    return pairs_;
  }
};

}
