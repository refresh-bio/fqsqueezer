// *******************************************************************************************
// This file is a part of FQSqueezer software distributed under GNU GPL 3 licence.
// The homepage of the FQSqueezer project is http://sun.aei.polsl.pl/REFRESH/fqsqueezer
//
// Author: Sebastian Deorowicz
// Version: 1.1
// Date   : 2020-06-16
// *******************************************************************************************

#include "ht_kmer.h"

using namespace std;

// ************************************************************************************
// CHT_pair_kmers
// ************************************************************************************
CHT_pair_kmers::CHT_pair_kmers(uint32_t _kmer_size, uint64_t _no_parts)
{
	no_parts = _no_parts;
	kmer_size = _kmer_size;

	ht_desc = new ht_desc_t[no_parts];

	for (size_t i = 0; i < no_parts; ++i)
	{
		ht_desc[i].size = HT_INIT_SIZE;
		ht_desc[i].filled = 0;
		ht_desc[i].max_filled = (uint64_t)(HT_INIT_SIZE * max_load_must);
		ht_desc[i].ht = new item_t[ht_desc[i].size + 1];
	}

	counter_shift = 2 * kmer_size;
	value_mask = (1ull << counter_shift) - 1ull;
	max_counter = (~0ull) >> counter_shift;

	for (auto i = 0ull; i < no_parts; ++i)
		_clear(i);
}

// ************************************************************************************
CHT_pair_kmers::~CHT_pair_kmers()
{	
	for(auto i = 0ull; i < no_parts; ++i)
		delete[] ht_desc[i].ht;

	delete[] ht_desc;
}

// ************************************************************************************
void CHT_pair_kmers::restruct(uint64_t id_part)
{
	auto old_ht = ht_desc[id_part].ht;
	auto old_ht_size = ht_desc[id_part].size;

	ht_desc[id_part].size *= 2;
	ht_desc[id_part].max_filled = (uint64_t)(ht_desc[id_part].size * max_load_must);
	ht_desc[id_part].filled = 0;

	ht_desc[id_part].ht = new item_t[ht_desc[id_part].size + 1];
	_clear(id_part);

	const uint64_t step_pe = 32;

	for (uint64_t i = 0; i < old_ht_size; ++i)
	{
		if (i + step_pe < old_ht_size)
			if (old_ht[i + step_pe] != EMPTY)
				prefetch(id_part, old_ht[i + step_pe].first);

		if (old_ht[i] != EMPTY)
			insert_raw(old_ht[i].first, old_ht[i].second);
	}

	delete[] old_ht;
}

// ************************************************************************************
void CHT_pair_kmers::_clear(uint64_t id_part)
{
	fill_n(ht_desc[id_part].ht, ht_desc[id_part].size + 1, EMPTY);

	ht_desc[id_part].filled = 0;
}

// ************************************************************************************
void CHT_pair_kmers::insert_raw(uint64_t key, uint64_t value_counter)
{
	uint64_t hash_key = hash(key);
	uint64_t id_part = get_part_id_from_hash(hash_key);
	uint64_t pos;
	   
	auto loc_ht_mask = ht_desc[id_part].size - 1ull;
	auto loc_ht = ht_desc[id_part].ht;

	for (pos = hash_key & loc_ht_mask; loc_ht[pos] != EMPTY; ++pos)
		;
	if(pos > loc_ht_mask)
		for (pos = 0; loc_ht[pos] != EMPTY; ++pos)
			;

	loc_ht[pos] = make_pair(key, value_counter);
	++ht_desc[id_part].filled;
}

// ************************************************************************************
void CHT_pair_kmers::prefetch(uint64_t part_id, uint64_t key)
{
	if (key == value_mask)
		return;

	uint64_t hash_key = hash(key);

#ifdef _WIN32
	_mm_prefetch((const char*)(ht_desc[part_id].ht + hash_key), _MM_HINT_T0);
#else
	__builtin_prefetch(ht_desc[part_id].ht + hash_key);
#endif
}

// ************************************************************************************
void CHT_pair_kmers::insert(uint64_t key, uint64_t value, uint64_t count)
{
	if (key == value_mask || value == value_mask)
		return;

	uint64_t hash_key = hash(key);
	uint64_t id_part = get_part_id_from_hash(hash_key);

	if (ht_desc[id_part].filled == ht_desc[id_part].max_filled)
		restruct(id_part);

	uint64_t pos;

	auto loc_ht = ht_desc[id_part].ht;
	auto loc_ht_mask = ht_desc[id_part].size - 1ull;

	for (pos = hash_key & loc_ht_mask; ; ++pos)
	{
		if (loc_ht[pos] == EMPTY)
		{
			if(pos > loc_ht_mask)
				break;

			if (count > max_counter)
				count = max_counter;

			loc_ht[pos].first = key;
			loc_ht[pos].second = value + (count << counter_shift);
			++ht_desc[id_part].filled;
			break;
		}

		if (loc_ht[pos].first == key && (loc_ht[pos].second & value_mask) == value)
		{
			uint64_t cur_counter = loc_ht[pos].second >> counter_shift;
			if (cur_counter + count < max_counter)
				loc_ht[pos].second += count << counter_shift;
			else
				loc_ht[pos].second += (max_counter - cur_counter) << counter_shift;
			break;
		}
	}

	if(pos > loc_ht_mask)
		for (pos = 0; ; ++pos)
		{
			if (loc_ht[pos] == EMPTY)
			{
				if (count > max_counter)
					count = max_counter;

				loc_ht[pos].first = key;
				loc_ht[pos].second = value + (count << counter_shift);
				++ht_desc[id_part].filled;
				break;
			}

			if (loc_ht[pos].first == key && (loc_ht[pos].second & value_mask) == value)
			{
				uint64_t cur_counter = loc_ht[pos].second >> counter_shift;
				if (cur_counter + count < max_counter)
					loc_ht[pos].second += count << counter_shift;
				else
					loc_ht[pos].second += (max_counter - cur_counter) << counter_shift;
				break;
			}
		}
}

// ************************************************************************************
uint64_t CHT_pair_kmers::count(uint64_t key, uint64_t value)
{
	auto hash_key = hash(key);
	uint64_t id_part = get_part_id_from_hash(hash_key);

	auto loc_ht = ht_desc[id_part].ht;
	auto loc_ht_mask = ht_desc[id_part].size - 1ull;
	uint64_t pos;

	for (pos = hash_key & loc_ht_mask; loc_ht[pos] != EMPTY; ++pos)
		if (loc_ht[pos].first == key && (loc_ht[pos].second & value_mask) == value)
			return loc_ht[pos].second >> counter_shift;

	if(pos > loc_ht_mask)
		for (pos = 0; loc_ht[pos] != EMPTY; ++pos)
			if (loc_ht[pos].first == key && (loc_ht[pos].second & value_mask) == value)
				return loc_ht[pos].second >> counter_shift;

	return 0;
}

// ************************************************************************************
void CHT_pair_kmers::find(uint64_t key, vector<uint64_t> &v_value_counter)
{
	auto hash_key = hash(key);
	uint64_t id_part = get_part_id_from_hash(hash_key);

	auto loc_ht = ht_desc[id_part].ht;
	auto loc_ht_mask = ht_desc[id_part].size - 1ull;
	uint64_t pos;

	for (pos = hash_key & loc_ht_mask; loc_ht[pos] != EMPTY; ++pos)
		if (loc_ht[pos].first == key)
			v_value_counter.emplace_back(loc_ht[pos].second);

	if (pos > loc_ht_mask)
		for (pos = 0; loc_ht[pos] != EMPTY; ++pos)
			if (loc_ht[pos].first == key)
				v_value_counter.emplace_back(loc_ht[pos].second);
}

// EOF
