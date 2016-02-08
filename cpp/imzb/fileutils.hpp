#pragma once

#include <fstream>

template <typename T>
void binary_read(std::ifstream& stream, T& value) {
  stream.read(reinterpret_cast<char*>(&value), sizeof(value));
}

template <typename T>
void binary_write(std::ofstream& stream, T& value) {
  stream.write(reinterpret_cast<char*>(&value), sizeof(value));
}
