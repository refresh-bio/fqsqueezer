#pragma once
// *******************************************************************************************
// This file is a part of FQSqueezer software distributed under GNU GPL 3 licence.
// The homepage of the MSAC project is http://sun.aei.polsl.pl/REFRESH/fqsqueezer
//
// Author: Sebastian Deorowicz
// Version: 1.0
// Date   : 2019-02-22
// *******************************************************************************************

#include "defs.h"
#include <mutex>
#include <nmmintrin.h>
#include <string>
#include <random>
#include <vector>
#include <algorithm>

using namespace std;

#include <chrono>
#include <cstdint>

// TimeHelpers from https://codereview.stackexchange.com/questions/132852/easy-to-use-c-class-for-asking-current-time-stamp-in-milli-micro-and-nanose/132863
namespace TimeHelpers
{
	namespace
	{
		inline std::chrono::high_resolution_clock Clock()
		{
			static std::chrono::high_resolution_clock clock;
			return clock;
		}

		template<typename Out, typename In>
		Out TimeAs()
		{
			return std::chrono::duration_cast<In>
				(Clock().now().time_since_epoch()).count();
		}
	}

	template<typename Out>
	inline Out TimeFromEpochInMilliSeconds()
	{
		return TimeAs<Out, std::chrono::milliseconds>();
	}

	template<typename Out>
	inline Out TimeFromEpochInMicroSeconds()
	{
		return TimeAs<Out, std::chrono::microseconds>();
	}

	template<typename Out>
	inline Out TimeFromEpochInNanoSeconds()
	{
		return TimeAs<Out, std::chrono::nanoseconds>();
	}
}

template<typename T>
constexpr T dna_code(const T x) 
{
	switch (x)
	{
	case (T) 'A':	return (T) 0; break;
	case (T) 'C':	return (T) 1; break;
	case (T) 'G':	return (T) 2; break;
	case (T) 'T':	return (T) 3; break;
	}

	return (T)4;
}


template<typename T>
constexpr T reverse_complement(const T x)
{
	switch (x)
	{
	case dna_code('A'):	return (T)dna_code('T'); break;
	case dna_code('C'):	return (T)dna_code('G'); break;
	case dna_code('G'):	return (T)dna_code('C'); break;
	case dna_code('T'):	return (T)dna_code('A'); break;
	}

	return (T) 4;
}

template<typename T>
constexpr T complement(const T x)
{
	switch (x)
	{
	case dna_code('A'):	return (T)dna_code('T'); break;
	case dna_code('C'):	return (T)dna_code('G'); break;
	case dna_code('G'):	return (T)dna_code('C'); break;
	case dna_code('T'):	return (T)dna_code('A'); break;
	}

	return (T)4;
}

constexpr uint8_t reverse_complement_alhpa(const uint8_t x)
{
	switch (x)
	{
	case 'A': return 'T'; break;
	case 'C': return 'G'; break;
	case 'G': return 'C'; break;
	case 'T': return 'A'; break;
	}

	return 'N';
}

template<typename T>
constexpr uint32_t pop_count(T x)
{
	uint32_t r = 0;

	for (; x; ++r)
		x &= x - 1;

	return r;
}


// ************************************************************************************
class CProgress
{
	mutex mtx;
	int show_comment;

	int64_t max_value;
	int64_t current_value;
	string comment;

	string prev_text;

	void show_progress(void);

public:
	CProgress();
	~CProgress();

	void Init(int64_t _max_value, bool _show_comment);
	void SetComment(string _comment);

	void Step(int64_t increment);
};

constexpr uint64_t ilog2(uint64_t x)
{
	uint64_t r = 0;

	for (; x; ++r)
		x >>= 1;

	return r;
}

constexpr uint64_t no_bytes(uint64_t x)
{
	uint64_t r = 1;

	x >>= 8;

	for (; x; ++r)
		x >>= 8;

	return r;
}

uint64_t popcnt(uint64_t x);
uint64_t notz(uint64_t x);

string trim(string s);

void sort_copy_stats(stats_t &dest, const stats_t &src);

constexpr uint64_t modulo_divisor(uint64_t x, int mod)
{
	switch (mod)
	{
	case 0: return 0;
	case 1: return 0;
	case 2: return x % 2;
	case 3: return x % 3;
	case 4: return x % 4;
	case 5: return x % 5;
	case 6: return x % 6;
	case 7: return x % 7;
	case 8: return x % 8;
	case 9: return x % 9;
	case 10: return x % 10;
	case 11: return x % 11;
	case 12: return x % 12;
	case 13: return x % 13;
	case 14: return x % 14;
	case 15: return x % 15;
	case 16: return x % 16;
	case 17: return x % 17;
	case 18: return x % 18;
	case 19: return x % 19;
	case 20: return x % 20;
	case 21: return x % 21;
	case 22: return x % 22;
	case 23: return x % 23;
	case 24: return x % 24;
	case 25: return x % 25;
	case 26: return x % 26;
	case 27: return x % 27;
	case 28: return x % 28;
	case 29: return x % 29;
	case 30: return x % 30;
	case 31: return x % 31;
	case 32: return x % 32;
	case 33: return x % 33;
	case 34: return x % 34;
	case 35: return x % 35;
	case 36: return x % 36;
	case 37: return x % 37;
	case 38: return x % 38;
	case 39: return x % 39;
	case 40: return x % 40;
	case 41: return x % 41;
	case 42: return x % 42;
	case 43: return x % 43;
	case 44: return x % 44;
	case 45: return x % 45;
	case 46: return x % 46;
	case 47: return x % 47;
	case 48: return x % 48;
	case 49: return x % 49;
	case 50: return x % 50;
	case 51: return x % 51;
	case 52: return x % 52;
	case 53: return x % 53;
	case 54: return x % 54;
	case 55: return x % 55;
	case 56: return x % 56;
	case 57: return x % 57;
	case 58: return x % 58;
	case 59: return x % 59;
	case 60: return x % 60;
	case 61: return x % 61;
	case 62: return x % 62;
	case 63: return x % 63;
	case 64: return x % 64;
	default: return x % mod;
	}
}

class CCounterIncrementer {
	mt19937 mt;
	uint32_t thr;
	uint32_t mult;
	uint32_t max_val;

	vector<uint32_t> v_mapping;

	uint32_t decode_real_value(uint32_t value)
	{
		if (value <= thr)
			return value;

		return (v_mapping[(uint64_t)value] + v_mapping[(uint64_t) value+1]) / 2;
	}

	uint32_t encode_real_value(uint32_t real_value)
	{
		if (real_value <= thr)
			return real_value;

		auto end_dist = min(max_val, real_value) + 1;
		auto p = upper_bound(v_mapping.begin() + thr, v_mapping.begin() + end_dist, real_value) - 1;

		uint32_t pos = (uint32_t) (p - v_mapping.begin());

		if (pos >= max_val)
			return max_val;

		uint32_t rest = real_value - *p;

		if (mt() % (v_mapping[(uint64_t)pos+1] - v_mapping[(uint64_t)pos]) < rest)
			++pos;

		return pos;
	}

public:
	CCounterIncrementer() : thr(numeric_limits<unsigned int>::max()), mult(1), max_val(1) {};
		
	void Reset(uint32_t _thr, uint32_t _mult, uint32_t _max_val)
	{
		mt.seed(5481);
		max_val = _max_val;
		thr = min(_thr, max_val);
		mult = _mult;

		v_mapping.resize((uint64_t)max_val + 2);

		for (uint32_t i = 0; i <= thr; ++i)
			v_mapping[i] = i;

		uint32_t inc = mult;
		for (uint32_t i = thr + 1; i <= max_val; ++i, inc += mult)
			v_mapping[i] = v_mapping[i - 1] + inc;

		v_mapping[(uint64_t)max_val + 1] = v_mapping[(uint64_t)max_val];		// guard
	};

	uint32_t Increment(uint32_t counter)
	{
		if (counter <= thr)
			return counter + 1;
		else
		{
			if (mt() % (mult * (counter - thr)) == 0)
				return counter + 1;
			else
				return counter;
		}
	};

	uint32_t Increment(uint32_t counter, uint32_t increment)
	{
		uint32_t real_counter = decode_real_value(counter);
		uint32_t real_inc = decode_real_value(increment);

		return encode_real_value(real_counter + real_inc);
	};
};

// EOF
