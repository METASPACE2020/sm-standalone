#include "ms/isocalc.hpp"
#include "imzb/reader.hpp"

#include "utils/isotope_pattern_db.hpp"

#include "boost/program_options.hpp"
#include "boost/algorithm/string/split.hpp"

#include <iostream>
#include <iomanip>
#include <vector>
#include <cassert>
#include <string>

struct Metrics {
  double img_corr;
  double iso_corr;
  double moc;
};

int main(int argc, char** argv) {
  std::string db_filename, imzb_filename;
  std::string output_filename;
  double ppm;
  bool remove_hotspots;

  namespace po = boost::program_options;

  po::options_description desc("Options");
  desc.add_options()
    ("ppm", po::value<double>(&ppm)->default_value(3.0),
     "m/z-window half-width in ppm")
    ("out", po::value<std::string>(&output_filename)->default_value("/dev/stdout"))
    ("remove-hotspots", po::value<bool>(&remove_hotspots)->default_value(true),
     "apply hotspot removal prior to computing metrics")
    ("help", "Print help");

  po::options_description hidden;
  hidden.add_options()
    ("db_filename", po::value<std::string>(&db_filename))
    ("input_file", po::value<std::string>(&imzb_filename));

  po::positional_options_description p;
  p.add("db_filename", 1);
  p.add("input_file", 1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv)
      .options(po::options_description().add(desc).add(hidden)).positional(p).run(), vm);
  po::notify(vm);

  if (vm.count("help") || imzb_filename.empty() || db_filename.empty()) {
    std::cout << "Usage: detect [options] <isotope_patterns.db> <input.imzb>\n";
    std::cout << desc << std::endl;
    return 0;
  }

  utils::IsotopePatternDB db(db_filename);
  imzb::ImzbReader reader(imzb_filename);
  std::ofstream out(output_filename);

  auto& keys = db.keys();
  std::vector<Metrics> metrics(keys.size());

  thread_local std::vector<ims::ImageF> images;
  thread_local std::vector<float> hotspot_removal_buf;

#pragma omp parallel for
  for (size_t i = 0; i < keys.size(); i++) {
    std::string f = keys[i].first;
    std::string adduct = keys[i].second;

    auto p = db(f, adduct);
    if (p.size() == 0)
      continue;

    if (images.size() < p.size()) {
      for (size_t i = images.size(); i < p.size(); i++)
        images.emplace_back(reader.height(), reader.width());
    }

    if (remove_hotspots && hotspot_removal_buf.empty())
      hotspot_removal_buf.resize(reader.width() * reader.height());

    for (size_t j = 0; j < p.size(); ++j) {
      reader.readImage(p.masses[j], ppm, images[j].rawPtr());
      if (remove_hotspots)
        images[j].removeHotspots(99.0, &hotspot_removal_buf[0]);
    }

    auto img_corr = ims::isotopeImageCorrelation(&images[0], p.size(), p);
    auto iso_corr = ims::isotopePatternMatch(&images[0], p.size(), p);
    auto moc = ims::measureOfChaos(images[0], 30);
    metrics[i] = Metrics{img_corr, iso_corr, moc};
  }

  for (size_t i = 0; i < keys.size(); ++i) {
    std::string f = keys[i].first;
    std::string adduct = keys[i].second;
    out << f << "\t" << adduct << "\t" <<
      metrics[i].img_corr << "\t" << metrics[i].iso_corr << "\t" << metrics[i].moc << "\n";
  }

  return 0;
}
