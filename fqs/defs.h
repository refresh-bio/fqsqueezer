#pragma once
// *******************************************************************************************
// This file is a part of FQSqueezer software distributed under GNU GPL 3 licence.
// The homepage of the FQSqueezer project is http://sun.aei.polsl.pl/REFRESH/fqsqueezer
//
// Author: Sebastian Deorowicz
// Version: 1.1
// Date   : 2020-06-16
// *******************************************************************************************

#include <array>
#include <string>
using namespace std;

typedef array<uint32_t, 4> stats_t;
typedef array<uint64_t, 4> stats64_t;

//#define ENABLE_STATS

#define KMER_MARGIN_2_SYMBOLS

const uint32_t BIN_PREFIX_LEN = 4;
const uint32_t NO_BINS = (1u << (2 * BIN_PREFIX_LEN));
const uint32_t BUF_SIZE = 128 << 10;
const uint32_t SIV_FIELD_SIZE = 2;
const uint32_t SMER_COUNTER_BITS = 12;
const uint32_t BMER_COUNTER_BITS = 6;
const uint32_t IMER_COUNTER_BITS = 6;
const uint32_t HMER_COUNTER_BITS = 6;
const uint32_t SMER_MAX_COUNTER_VALUE = (1u << SMER_COUNTER_BITS) - 1u;
const uint32_t BMER_MAX_COUNTER_VALUE = (1u << BMER_COUNTER_BITS) - 1u;
const uint32_t PMER_MAX_COUNTER_VALUE = (1u << SIV_FIELD_SIZE) - 1u;

const size_t NO_STREAMS     = 4;
const size_t STREAM_META    = 0;
const size_t STREAM_ID      = 1;
const size_t STREAM_DNA     = 2;
const size_t STREAM_QUALITY = 3;

const std::string FQS_VERSION = std::string("FQSqueezer (k-mer-based compression of sequencing data) 1.1");

const unsigned char dna_alphabet[] = "ACGTN";

typedef uint64_t context_t;

enum class counts_level_t {none = 0, pmer = 1, smer = 2, bmer = 3, mixed = 4, bmer_unc = 5};

template<typename T> T NormalizeValue(T val, T min_val, T max_val)
{
	if (val < min_val)
		return min_val;
	if (val > max_val)
		return max_val;
	return val;
}

struct read_desc_t
{
	uint8_t *id;
	uint8_t *dna;
	uint8_t *plus;
	uint8_t *quality;

	read_desc_t() : id(nullptr), dna(nullptr), plus(nullptr), quality(nullptr)
	{};

	read_desc_t(uint8_t* _id, uint8_t* _dna, uint8_t *_plus, uint8_t* _quality) : id(_id), dna(_dna), plus(_plus), quality(_quality)
	{};

	uint32_t id_len() {
		return (uint32_t) (dna - id);
	}

	// Length of the DNA symbols
	uint32_t read_len() {
		return (uint32_t)(plus - dna - 1);
	}

	// Size of all fields of the read
	uint32_t read_size() {
		return (uint32_t)(quality - id + read_len() + 1);
	}
};

// EOF
