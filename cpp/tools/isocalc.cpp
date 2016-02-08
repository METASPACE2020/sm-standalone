#include "utils/isotope_pattern_db.hpp"
#include "boost/program_options.hpp"
#include "boost/algorithm/string/split.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>

int main(int argc, char** argv) {

  double resolution;
  unsigned max_peaks;

  namespace po = boost::program_options;

  std::string input_file, output_file;

  std::string adducts_str;

  po::options_description desc("Options");
  desc.add_options()
    ("resolution", po::value<double>(&resolution)->default_value(140000.0),
     "Resolving power at m/z=200")
    ("adducts", po::value<std::string>(&adducts_str)->default_value("+H,+K,+Na"),
     "Comma-separated list of adducts")
    ("max-peaks", po::value<unsigned>(&max_peaks)->default_value(5),
     "Maximum number of peaks to store")
    ("help", "Print help");

  po::options_description hidden;
  hidden.add_options()
    ("input-file", po::value<std::string>(&input_file),
     "List of sum formulas, one per line")
    ("output-file", po::value<std::string>(&output_file),
     "Output file");

  po::positional_options_description p;
  p.add("input-file", 1);
  p.add("output-file", 1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv)
      .options(po::options_description().add(desc).add(hidden)).positional(p).run(), vm);
  po::notify(vm);

  if (vm.count("help") || input_file.empty() || output_file.empty()) {
    std::cout << "Usage: isocalc [options] <input.txt> <output.db>\n";
    std::cout << "       where input contains one sum formula per line.\n";
    std::cout << desc << std::endl;
    return 0;
  }

  std::vector<std::string> target_adducts;
  boost::split(target_adducts, adducts_str, [](char c) { return c == ','; });
  std::sort(target_adducts.begin(), target_adducts.end());
  auto it = std::unique(target_adducts.begin(), target_adducts.end());
  target_adducts.erase(it, target_adducts.end());

  std::cout << "Resolution @ m/z=200: " << resolution << std::endl;
  std::cout << "Generating isotope pattern database for adducts ";
  for (auto& a: target_adducts) {
    std::cout << a;
    if (a != target_adducts.back()) std::cout << ", ";
  }
  std::cout << "..." << std::endl;

  std::ifstream in(input_file);

  std::vector<std::string> sum_formulas;

  std::string f;
  while (!in.eof()) {
    std::getline(in, f);
    if (ms::monoisotopicMass(f) > 100 && ms::monoisotopicMass(f) < 2000) {
      sum_formulas.push_back(f);
    }
  }

  utils::IsotopePatternDB db{sum_formulas, target_adducts};
  utils::OrbitrapProfile orbitrap{resolution};
  db.computeIsotopePatterns(orbitrap, output_file, max_peaks);

  std::cout << "Isotope patterns have been saved to " << output_file << std::endl;

  return 0;
}
