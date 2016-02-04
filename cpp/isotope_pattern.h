#pragma once

#include <vector>
#include <cstdlib>
#include <initializer_list>

namespace ms {
  struct IsotopePattern {
    std::vector<double> masses;
    std::vector<double> abundances;

    IsotopePattern() {}
    IsotopePattern(double mass) : masses{mass}, abundances{1.0} {}
    IsotopePattern(std::initializer_list<double> masses,
                   std::initializer_list<double> abundances) : masses{masses}, abundances{abundances}
    {}

    bool isUnit() const { return masses.size() == 1 && masses[0] == 0.0; }

    ms::IsotopePattern multiply(const IsotopePattern& other, double threshold = 0.0) const;

    // sorts by decreasing intensity and scales max intensity to 1
    void normalize();

    size_t size() const { return masses.size(); }
  };
}
