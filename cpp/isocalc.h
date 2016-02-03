#pragma once

#include "periodic_table.h"

#include <iostream>

#include <string>
#include <sstream>
#include <map>
#include <vector>
#include <cstdint>

namespace ms {
  typedef std::map<std::string, uint16_t> ElementCounter;
}

namespace sf_parser {

  typedef ms::ElementCounter ElementCounter;

  class ParseError : public std::exception {
    std::string s_;
    size_t pos_;
  public:
    ParseError(const std::string& str, size_t pos) : s_(str), pos_(pos) {}
    virtual const char* what() const noexcept {
      std::ostringstream ss;
      ss << s_ << " at position " << pos_;
      return ss.str().c_str();
    }
  };

  class NegativeTotalError : public std::exception {
    std::string element_;
    int total_;
  public:
    NegativeTotalError(const std::string& element, int total) : element_(element), total_(total) {}
    virtual const char* what() const noexcept {
      std::ostringstream ss;
      ss << "total number of " << element_ << " elements (" << total_ << ") is less than zero";
      return ss.str().c_str();
    }
  };

  class SumFormulaParser {
    std::string s;
    size_t n;
    ElementCounter counter_;
  public:
    SumFormulaParser(const std::string& input) : s(input), n(0) {
      parseSumFormula(counter_);
    }

    const ElementCounter& elementCounts() const {
      return counter_;
    }

  private:
    bool eof() const { return n >= s.length(); }

    void checkAvailable() const {
      if (eof()) throw new ParseError("unexpected end of input", n);
    }

    char nextChar() { checkAvailable(); return s[n++]; }

    char peek() const { checkAvailable(); return s[n]; }

    uint16_t parseOptionalNumber() {
      if (eof() || !std::isdigit(peek()))
        return 1;
      auto pos = n;
      uint16_t result = 0;
      while (!eof() && std::isdigit(peek())) {
        if (result >= 6553)
          throw new ParseError("the number is too big", pos);
        result = result * 10 + (nextChar() - '0');
      }
      return result;
    }

    void parseElement(ElementCounter& counter) {
      std::string element;
      element.push_back(nextChar());
      auto pos = n;
      if (!std::isupper(element.back())) throw new ParseError("expected an element", n);
      while (!eof() && std::islower(peek()))
        element.push_back(nextChar());
      uint16_t num = parseOptionalNumber();
      if (!ms::Element::isKnown(element))
        throw new ParseError("unknown element " + element, pos);
      counter[element] += num;
    }

    void parseSimpleFragment(ElementCounter& counter) {
      while (!eof() && std::isupper(peek())) {
        parseElement(counter);
      }
    }

    void parseFragment(ElementCounter& counter) {
      while (!eof()) {
        if (peek() == '(') {
          nextChar();
          ElementCounter tmp;
          parseFragment(tmp);
          if (nextChar() != ')')
            throw new ParseError("expected closing parenthesis", n - 1);
          auto repeats = parseOptionalNumber();
          for (auto& item: tmp)
            counter[item.first] += item.second * repeats;
        } else if (std::isupper(peek())) {
          parseSimpleFragment(counter);
        } else {
          break;
        }
      }
    }

    void parseMolecularComplex(ElementCounter& counter) {
      ElementCounter tmp;
      auto repeats = parseOptionalNumber();
      while (!eof()) {
        if (peek() == '.' || peek() == '-' || peek() == '+')
          break;
        parseFragment(tmp);
      }
      for (auto& item: tmp)
        counter[item.first] += repeats * item.second;
    }

    void parseSumFormula(ElementCounter& counter) {
      parseMolecularComplex(counter);

      while (!eof()) {
        if (peek() == '.') {
          nextChar();
          parseMolecularComplex(counter);
        } else {
          break;
        }
      }

      if (!eof()) {
        ElementCounter adduct;
        char sign = nextChar();
        int mult;
        if (sign == '-') mult = -1;
        else if (sign == '+') mult = 1;
        else throw new ParseError("expected +/-", n-1);
        parseMolecularComplex(adduct);
        for (auto& item: adduct) {
          auto total = int(counter[item.first]) + mult * int(item.second);
          if (total < 0)
            throw new NegativeTotalError(item.first, total);
          counter[item.first] = total;
        }
      }
    }
  };

  ElementCounter parseSumFormula(const std::string& formula) {
    sf_parser::SumFormulaParser parser{formula};
    return parser.elementCounts();
  }
}
