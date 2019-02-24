#pragma once
// *******************************************************************************************
// This file is a part of FQSqueezer software distributed under GNU GPL 3 licence.
// The homepage of the MSAC project is http://sun.aei.polsl.pl/REFRESH/fqsqueezer
//
// Author: Sebastian Deorowicz
// Version: 1.0
// Date   : 2019-02-22
// *******************************************************************************************

#include <algorithm>
#include "xmmintrin.h"
#include "utils.h"

using namespace std;

template<unsigned FIELD_SIZE> class TSmallIntVector {
	uint32_t key_size;
	uint64_t *data;
	uint64_t data_size;
	const uint64_t divider = 64ull / FIELD_SIZE;
	const uint64_t mask = (1ull << FIELD_SIZE) - 1ull;
	const uint32_t idx_shift = (uint32_t) (ilog2(64 / FIELD_SIZE) - 1u);
	const uint32_t idx_mask = 64u / FIELD_SIZE - 1u;
	atomic<uint64_t> no_updates;
	atomic<uint64_t> no_filled;

public:
	TSmallIntVector(const uint32_t _key_size)
	{
		key_size = _key_size;

		data_size = (1ull << key_size) / divider;

		data = new uint64_t[data_size];
		fill_n(data, data_size, 0ull);

		no_updates = 0;
		no_filled = 0;
	}

	~TSmallIntVector()
	{
		delete[] data;
	}

	void prefetch(const uint64_t idx) const
	{
		uint64_t word = idx >> idx_shift;
		_mm_prefetch((const char *)(data+word), _MM_HINT_T0);
	}

	bool increment(const uint64_t idx)
	{
		uint64_t word = idx >> idx_shift;
		uint64_t reminder = idx & idx_mask;
		uint64_t cur_mask = mask << (FIELD_SIZE * reminder);

		if ((data[word] & cur_mask) == cur_mask)
			return false;		// max field value

		bool new_item = (data[word] & cur_mask) == 0;

		data[word] += 1ull << (FIELD_SIZE * reminder);

		return new_item;
	}

	uint64_t test(const uint64_t idx) const
	{
		uint64_t word = idx >> idx_shift;
		uint64_t reminder = idx & idx_mask;

		if (word >= data_size)
		{
			cout << "test ERROR:   data_size: " + to_string(data_size) + "   word: " + to_string(word) + "\n";
			fflush(stdout);
		}

		return (data[word] >> (FIELD_SIZE * reminder)) & mask;
	}

	void counts(const uint64_t idx, stats_t &counts) const
	{
		uint64_t word = idx >> idx_shift;
		uint64_t reminder = (idx & idx_mask) & ~3ull;
		
		uint32_t rem_FS = (uint32_t) (FIELD_SIZE * reminder);
		counts[0] = (uint32_t) ((data[word] >> rem_FS) & mask);
		rem_FS += FIELD_SIZE;
		counts[1] = (uint32_t) ((data[word] >> rem_FS) & mask);
		rem_FS += FIELD_SIZE;
		counts[2] = (uint32_t) ((data[word] >> rem_FS) & mask);
		rem_FS += FIELD_SIZE;
		counts[3] = (uint32_t) ((data[word] >> rem_FS) & mask);
	}

	void counts_update(const uint64_t idx, stats_t &counts) const
	{
		uint64_t word = idx >> idx_shift;
		uint64_t reminder = (idx & idx_mask) & ~3ull;

		uint32_t rem_FS = (uint32_t) (FIELD_SIZE * reminder);
		counts[0] += (uint32_t) ((data[word] >> rem_FS) & mask);
		rem_FS += FIELD_SIZE;
		counts[1] += (uint32_t)((data[word] >> rem_FS) & mask);
		rem_FS += FIELD_SIZE;
		counts[2] += (uint32_t)((data[word] >> rem_FS) & mask);
		rem_FS += FIELD_SIZE;
		counts[3] += (uint32_t)((data[word] >> rem_FS) & mask);
	}

	uint64_t test_shorter(const uint64_t idx, const uint64_t size) const
	{
		uint64_t shift = key_size - size;
		uint64_t start = idx << shift;
		uint64_t end = (idx + 1) << shift;

		if (end - start >= 128ull / FIELD_SIZE)
			return sum_shorter_fast(start, end);
		else
			return sum_shorter(start, end);
	}

	uint64_t sum_shorter(const uint64_t start, const uint64_t end) const
	{
		uint64_t r = 0;

		for (uint64_t x = start; x < end; ++x)
			r += test(x);

		return r;
	}

	uint64_t sum_shorter_fast(const uint64_t start, const uint64_t end) const
	{
		uint64_t r = 0;

		uint64_t word_start = start >> idx_shift;
		uint64_t reminder = start & idx_mask;

		if (reminder)
		{
			for (uint64_t d = data[word_start] >> (FIELD_SIZE * reminder); d; d >>= FIELD_SIZE)
				r += d & mask;
			++word_start;
		}

		uint64_t word_end = end >> idx_shift;
		reminder = end & idx_mask;

		if (reminder)
		{
			uint64_t m = (~0ull) << (FIELD_SIZE * reminder);
			for (uint64_t d = data[word_end] & ~m; d; d >>= FIELD_SIZE)
				r += d & mask;

			--word_end;
		}

		for (uint64_t word = word_start; word < word_end; ++word)
			for (uint64_t d = data[word]; d; d >>= FIELD_SIZE)
				r += d & mask;

		return r;
	}

	uint64_t conditionally_clear(const uint64_t idx, const uint64_t thr)
	{
		uint64_t word = idx >> idx_shift;
		uint64_t reminder = idx & idx_mask;
		uint64_t cur_mask = mask << (FIELD_SIZE * reminder);
		uint64_t val = (data[word] & cur_mask);

		if (val == 0 || val > thr)
			return 0;		

		data[word] -= val << (FIELD_SIZE * reminder);

		return val;
	}

	// Clear every element smaller or equal to thr
	void prune(const uint64_t thr = 1)
	{
		uint64_t no_items = 1ull << key_size;
		uint64_t del_items = 0;
		uint64_t del_counts = 0;

		for (uint64_t i = 0; i < no_items; ++i)
		{
			uint64_t x = conditionally_clear(i, thr);
			if (x)
			{
				++del_items;
				del_counts += x;
			}
		}

		no_filled -= del_items;
		no_updates -= del_counts;
	}

	double avg_filling_factor() const
	{
		if (no_filled)
			return (double)(no_updates) / (no_filled);

		return 0.0;
	}

	void update_no_filled(uint64_t _no_filled) 
	{
		no_filled += _no_filled;
	}

	void update_no_updates(uint64_t _no_updates)
	{
		no_updates += _no_updates;
	}

	uint64_t get_no_pmers() const
	{
		return no_filled;
	}

	uint64_t get_no_updates() const
	{
		return no_updates;
	}
};

// EOF
