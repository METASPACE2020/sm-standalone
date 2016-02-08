#pragma once

#include "ms/isocalc.hpp"
#include "msgpack.hpp"

#include <map>
#include <string>
#include <vector>
#include <mutex>
#include <fstream>

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

  void computeIsotopePatterns(const utils::InstrumentProfile& instrument,
                              const std::string& output_filename,
                              size_t max_peaks=5)
  {
    // use a simple structure that doesn't require special msgpack methods
    std::map<SFAdductPair, std::map<std::string, std::vector<double>>> patterns;

    std::mutex map_mutex;

#pragma omp parallel for shared(patterns)
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
      patterns[key]["mzs"] = pattern.masses;
      patterns[key]["abundances"] = pattern.abundances;
    }

    // TODO: write instrument profile parameters?
    msgpack::sbuffer sbuf;
    msgpack::pack(sbuf, patterns);
    std::ofstream db(output_filename, std::ios::binary);
    db.write(sbuf.data(), sbuf.size());
  }
};

}
