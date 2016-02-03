#pragma once

#include <map>
#include <vector>
#include <string>

namespace ms {
  struct IsotopePattern {
    std::vector<double> masses;
    std::vector<float> abundances;

    size_t size() const { return masses.size(); }
  };

  struct Element;
  extern const std::map<std::string, ms::Element> periodic_table;

  struct Element {
    std::string abbr;
    unsigned number;
    ms::IsotopePattern isotope_pattern;

    static bool isKnown(const std::string& abbr) {
      return periodic_table.find(abbr) != periodic_table.end();
    }

    static Element getByName(const std::string& abbr) {
      return periodic_table.find(abbr)->second;
    }
  };
}
