#pragma once

#include "ims.h"
#include "blosc.h"
#include "imzb_writer.h"

#include <string>
#include <fstream>
#include <vector>
#include <cstdint>

namespace imzb {
	
class ImzbReader {
	std::ifstream in_;

	std::vector<double> mzs_;
	std::vector<uint64_t> offsets_;
	size_t block_idx_;

	std::vector<char> buffer_;
	std::vector<ims::Peak> peaks_;
	size_t n_peaks_;
	size_t pos_;

	template <typename T>
	void binary_read(std::ifstream& stream, T& value) {
		stream.read(reinterpret_cast<char*>(&value), sizeof(value));
	}

	bool readNextBlock() {
		++block_idx_;
		if (block_idx_ == mzs_.size()) {
			n_peaks_ = 0;
			return false;
		}
		uint64_t start = offsets_[block_idx_];
		uint64_t end = offsets_[block_idx_ + 1];
		buffer_.resize(end - start);
		in_.seekg(start);
		in_.read(&buffer_[0], end - start);
		n_peaks_ = blosc_decompress_ctx(&buffer_[0], &peaks_[0],
				peaks_.size() * sizeof(ims::Peak), 1) / sizeof(ims::Peak);
		pos_ = 0;
		return true;
	}

	bool empty_;

public:
	ImzbReader(const std::string& filename) :
		in_(filename), block_idx_(0), peaks_(imzb::ImzbWriter::BLOCK_SIZE),
		n_peaks_(0), pos_(0), empty_(false)
	{
		std::ifstream in_idx(filename + ".idx");
		double mz;
		uint64_t offset;
		while (!in_idx.eof()) {
			binary_read(in_idx, mz);
			binary_read(in_idx, offset);
			mzs_.push_back(mz);
			offsets_.push_back(offset);
		}

		in_.seekg(0, in_.end);
		offsets_.push_back(in_.tellg());
		in_.seekg(0, in_.beg);
	}

	bool readNext(ims::Peak& peak) {
		if (empty_)
			return false;

		if (pos_ >= n_peaks_) {
			if (!readNextBlock()) {
				empty_ = true;
				return false;
			}
		}

		peak = peaks_[pos_];
		++pos_;
		return true;
	}
};

}
