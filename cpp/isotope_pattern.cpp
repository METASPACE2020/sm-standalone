#include "isotope_pattern.h"
#include "periodic_table.h"

#include <vector>
#include <algorithm>

namespace ms {

void IsotopePattern::addCharge(int charge) {
  for (auto& m: masses)
    m -= charge * ms::electronMass;
}

IsotopePattern IsotopePattern::multiply(const IsotopePattern& other, double threshold) const {
  if (this->isUnit()) return other;
  if (other.isUnit()) return *this;

  const auto& p1 = *this, p2 = other;
  size_t n1 = p1.size(), n2 = p2.size();
  IsotopePattern result;
  for (size_t i = 0; i < n1; i++)
    for (size_t j = 0; j < n2; j++) {
      auto abundance = p1.abundances[i] * p2.abundances[j];
      if (abundance > threshold) {
        result.masses.push_back(p1.masses[i] + p2.masses[j]);
        result.abundances.push_back(abundance);
      }
    }
  result.normalize();
  return result;
}

void IsotopePattern::normalize() {
  for (size_t i = 0; i < size(); i++)
    for (size_t j = i + 1; j < size(); j++)
      if (abundances[i] < abundances[j]) {
        std::swap(abundances[i], abundances[j]);
        std::swap(masses[i], masses[j]);
      }
  auto top = abundances[0];
  for (auto& item: abundances)
    item /= top;
}

}
