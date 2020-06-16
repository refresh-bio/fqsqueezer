#pragma once
// *******************************************************************************************
// This file is a part of FQSqueezer software distributed under GNU GPL 3 licence.
// The homepage of the FQSqueezer project is http://sun.aei.polsl.pl/REFRESH/fqsqueezer
//
// Author: Sebastian Deorowicz
// Version: 1.1
// Date   : 2020-06-16
// *******************************************************************************************

#include <algorithm>
#include <atomic>
#include <tuple>
#include <iostream>
#include <map>
#include <vector>

#include "defs.h"
#include "kmer.h"
#include "xmmintrin.h"
#include "utils.h"
#include <cassert>

using namespace std;

// ************************************************************************************
// CHT_kmer
// ************************************************************************************
template<typename T> class CHT_kmer
{
	typedef T item_t;

	const item_t EMPTY = 0;
	const uint32_t HT_INIT_SIZE = 16ull;

	uint64_t kmer_len;
	uint64_t counter_size;
	uint64_t prefix_len;
	uint64_t tab_size;
	double max_load_factor;

	item_t counter_mask;
	item_t counter_max;

	uint64_t total_memory;
	atomic<uint64_t> no_kmers;
	atomic<uint64_t> no_inserts;

	typedef struct {
		item_t *ht;
		uint32_t size;
		uint32_t filled;
		uint32_t max_filled;
	} ht_desc_t;
	
	ht_desc_t *ht_desc;

	// ************************************************************************************
	void init_ht()
	{
		ht_desc = new ht_desc_t[tab_size];

		total_memory += sizeof(ht_desc_t) * tab_size;

		for (uint64_t i = 0; i < tab_size; ++i)
		{
			ht_desc[i].size = HT_INIT_SIZE;
			ht_desc[i].filled = 0;
			ht_desc[i].max_filled = (uint32_t)(ht_desc[i].size * max_load_factor);

			ht_desc[i].ht = new item_t[ht_desc[i].size + 1ull];
			total_memory += sizeof(item_t) * ht_desc[i].size;

			fill_n(ht_desc[i].ht, ht_desc[i].size + 1ull, EMPTY);
		}
	}

	// ************************************************************************************
	void free_ht()
	{
		for (uint64_t i = 0; i < tab_size; ++i)
			delete[] ht_desc[i].ht;

		delete[] ht_desc;
	}

	// ************************************************************************************
	void restruct(const uint32_t idx)
	{
		item_t *old_ht = ht_desc[idx].ht;
		uint64_t old_size = ht_desc[idx].size;

		ht_desc[idx].size *= 2;
		ht_desc[idx].filled = 0;
		ht_desc[idx].max_filled = (uint32_t)(ht_desc[idx].size * max_load_factor);

		ht_desc[idx].ht = new item_t[ht_desc[idx].size + 1ull];
		fill_n(ht_desc[idx].ht, ht_desc[idx].size + 1, EMPTY);

		for (uint64_t i = 0; i < old_size; ++i)
			if (old_ht[i] != EMPTY)
			{
				uint64_t kmer = join_kmer(idx, old_ht[i] >> (8 * sizeof(item_t) - 2 * (kmer_len - prefix_len)));

				--no_kmers;
				insert(kmer, nullptr, get_count(old_ht[i]));
			}

		delete[] old_ht;

		total_memory += sizeof(item_t) * old_size;
	}

	// ************************************************************************************
	inline uint64_t hash(uint64_t reminder, const uint64_t idx) const
	{
#ifdef KMER_MARGIN_2_SYMBOLS
		reminder >>= 8;
#else
		reminder >>= 4;
#endif

		reminder ^= reminder >> 33;
		reminder *= 0xff51afd7ed558ccdL;
		reminder ^= reminder >> 33;
		reminder *= 0xc4ceb9fe1a85ec53L;
		reminder ^= reminder >> 33;

		return reminder & (ht_desc[idx].size - 1);
	}

	// ************************************************************************************
	inline uint32_t get_idx(const uint64_t kmer) const
	{
#ifdef KMER_MARGIN_2_SYMBOLS
		return (uint32_t)(prefix_len ? ((kmer << 4) >> (64 - 2 * prefix_len)) : 0);
#else
		return (uint32_t) (prefix_len ? ((kmer << 2) >> (64 - 2 * prefix_len)) : 0);
#endif
	}

	// ************************************************************************************
	inline pair<uint32_t, item_t> split_kmer(const uint64_t kmer) const
	{
#ifdef KMER_MARGIN_2_SYMBOLS
		uint64_t p1 = kmer >> 60;
		uint64_t p2 = prefix_len ? ((kmer << 4) >> (64 - 2 * prefix_len)) : 0;
		uint64_t p3 = (kmer << (2 * (prefix_len + 2))) >> (64 - 2 * (kmer_len - prefix_len - 2) - 4);
#else
		uint64_t p1 = kmer >> 62;
		uint64_t p2 = prefix_len ? ((kmer << 2) >> (64 - 2 * prefix_len)) : 0;
		uint64_t p3 = (kmer << (2 * (prefix_len + 1))) >> (64 - 2 * (kmer_len - prefix_len - 1) - 2);
#endif

		return make_pair((uint32_t) p2, (item_t)(p3 + p1));
	}

	// ************************************************************************************
	inline uint64_t join_kmer(const uint32_t idx, const item_t reminder) const
	{
#ifdef KMER_MARGIN_2_SYMBOLS
		uint64_t p1 = reminder & 0xf;
		uint64_t p2 = idx;
		uint64_t p3 = reminder >> 4;

		return (p1 << 60) + (p2 << (64 - 2 * (prefix_len + 2))) + (p3 << (64 - 2 * kmer_len));
#else
		uint64_t p1 = reminder & 3;
		uint64_t p2 = idx;
		uint64_t p3 = reminder >> 2;

		return (p1 << 62) + (p2 << (64 - 2 * (prefix_len + 1))) + (p3 << (64 - 2 * kmer_len));
#endif
	}

	// ************************************************************************************
	bool non_empty(const stats_t &counts) const
	{
		return counts[0] || counts[1] || counts[2] || counts[3];
	}

	// ************************************************************************************
	inline uint32_t get_count(const item_t item) const
	{
		return (uint32_t)(item & counter_mask);
	}

	// ************************************************************************************
	bool find_full(const CKmer &kmer, stats_t &counts) const
	{
		counts[0] = counts[1] = counts[2] = counts[3] = 0;

		update_counts_full(kmer, counts);

		return non_empty(counts);
	}

	// ************************************************************************************
	void update_counts_full(const CKmer &kmer, stats_t &counts) const
	{
		_update_counts_full(kmer.data_normalized(), counts, kmer.is_normalized_dir());
	}

	// ************************************************************************************
	void _update_counts_full(const uint64_t kmer_data, stats_t &counts, bool is_normalized_dir) const
	{
		uint32_t idx;
		item_t reminder;
		item_t reminder_left_aligned;

#ifdef KMER_MARGIN_2_SYMBOLS
		uint32_t last_letter_shift = (uint32_t)(8 * sizeof(item_t) - 2 * (kmer_len - prefix_len) + (is_normalized_dir ? 4 : 2));
#else
		uint32_t last_letter_shift = (uint32_t)(8 * sizeof(item_t) - 2 * (kmer_len - prefix_len) + (is_normalized_dir ? 2 : 0));
#endif

		item_t last_letter_mask = ((item_t)3u) << last_letter_shift;
		item_t match_mask = ~(counter_mask | last_letter_mask);

		tie(idx, reminder) = split_kmer(kmer_data);
		reminder_left_aligned = reminder << (8 * sizeof(item_t) - 2 * (kmer_len - prefix_len));
		item_t reminder_left_aligned_masked = reminder_left_aligned & match_mask;

		uint64_t pos = hash(reminder, idx);
		auto cur_ht = ht_desc[idx].ht;

		if (is_normalized_dir)
		{
			for (; cur_ht[pos] != EMPTY; ++pos)
				if ((cur_ht[pos] & match_mask) == reminder_left_aligned_masked)
				{
					uint32_t sym = (cur_ht[pos] >> last_letter_shift) & 3u;
					counts[sym] += get_count(cur_ht[pos]);
				}

			if(pos == ht_desc[idx].size)
				for (pos = 0; cur_ht[pos] != EMPTY; ++pos)
					if ((cur_ht[pos] & match_mask) == reminder_left_aligned_masked)
					{
						uint32_t sym = (cur_ht[pos] >> last_letter_shift) & 3u;
						counts[sym] += get_count(cur_ht[pos]);
					}
		}
		else
		{
			for (; cur_ht[pos] != EMPTY; ++pos)
				if ((cur_ht[pos] & match_mask) == reminder_left_aligned_masked)
				{
					uint32_t sym = (cur_ht[pos] >> last_letter_shift) & 3u;
					sym = reverse_complement(sym);
					counts[sym] += get_count(cur_ht[pos]);
				}

			if (pos == ht_desc[idx].size)
				for (pos = 0; cur_ht[pos] != EMPTY; ++pos)
					if ((cur_ht[pos] & match_mask) == reminder_left_aligned_masked)
					{
						uint32_t sym = (cur_ht[pos] >> last_letter_shift) & 3u;
						sym = reverse_complement(sym);
						counts[sym] += get_count(cur_ht[pos]);
					}
		}
	}

	// ************************************************************************************
	bool find_partial(const CKmer &kmer, stats_t &counts, CCounterIncrementer &cinc) const
	{
		counts[0] = counts[1] = counts[2] = counts[3] = 0;

		stats_t loc_counts;

		CKmer kmer_tmp = kmer;
		uint32_t no_missing_letters = kmer.get_max_size() - kmer.get_cur_size();
		uint32_t no_trials = 1;

		for (uint32_t i = 0; i < no_missing_letters; ++i)
		{
			kmer_tmp.insert_front(0);
			no_trials *= 4;
		}

		const uint32_t pf_step_idx = 8;
		const uint32_t pf_step_data = 4;

		pair<uint64_t, bool> *v_kmers_to_check = (pair<uint64_t, bool> *) alloca(((uint64_t) no_trials + pf_step_idx) * sizeof(pair<uint64_t, bool>));

		v_kmers_to_check[0] = make_pair(kmer_tmp.data_normalized(), kmer_tmp.is_normalized_dir());
		for(uint32_t i = 0; i < pf_step_idx; ++i)
			v_kmers_to_check[no_trials + i] = v_kmers_to_check[0];

		for (uint32_t i = 1; i < no_trials; ++i)
		{
			uint32_t j = 0;
			while (true)
			{
				auto c = kmer_tmp.get_symbol(j);
				if (c < 3)
				{
					kmer_tmp.replace(c + 1, j);
					break;
				}
				++j;
			}

			for (uint32_t k = 0; k < j; ++k)
				kmer_tmp.replace(0, k);

			v_kmers_to_check[i].first = kmer_tmp.data_normalized();
			v_kmers_to_check[i].second = kmer_tmp.is_normalized_dir();
		}

		for (uint32_t i = 0; i < no_trials; ++i)
		{
			prefetch_idx(v_kmers_to_check[i + pf_step_idx].first);
			prefetch_data(v_kmers_to_check[i + pf_step_data].first);

			loc_counts[0] = loc_counts[1] = loc_counts[2] = loc_counts[3] = 0;

			_update_counts_full(v_kmers_to_check[i].first, loc_counts, v_kmers_to_check[i].second);

			for (int i = 0; i < 4; ++i)
				if (loc_counts[i])
					counts[i] = cinc.Increment(counts[i], loc_counts[i]);
		}

		return non_empty(counts);
	}

	// ************************************************************************************
	pair<uint32_t, uint64_t> _find(const uint64_t kmer, const bool insert_mode) 
	{
		uint32_t idx;
		item_t reminder;
		item_t reminder_left_aligned;

		tie(idx, reminder) = split_kmer(kmer);
		reminder_left_aligned = reminder << (8 * sizeof(item_t) - 2 * (kmer_len - prefix_len));

		if (insert_mode && ht_desc[idx].filled == ht_desc[idx].max_filled)
			restruct(idx);

		uint64_t pos = hash(reminder, idx);
		auto cur_ht = ht_desc[idx].ht;

		for (; cur_ht[pos] != EMPTY; ++pos)
			if ((cur_ht[pos] & ~counter_mask) == reminder_left_aligned)
				break;
		if(pos == ht_desc[idx].size)
			for (pos = 0; cur_ht[pos] != EMPTY; ++pos)
				if ((cur_ht[pos] & ~counter_mask) == reminder_left_aligned)
					break;

		if (insert_mode && cur_ht[pos] == EMPTY)
		{
			cur_ht[pos] = reminder_left_aligned;
			++ht_desc[idx].filled;

			++no_kmers;
		}

		return make_pair(idx, pos);
	}

public:
	// ************************************************************************************
	CHT_kmer(uint64_t _kmer_len, uint64_t _counter_size, double _max_load_factor)
	{
		kmer_len = _kmer_len;
		counter_size = _counter_size;
		max_load_factor = _max_load_factor;

		counter_mask = (((item_t)1u) << counter_size) - (item_t)1;
		counter_max = counter_mask;

		// Adjust prefix_len for best performance
		uint64_t max_suffix_len = (8 * sizeof(item_t) - counter_size) / 2;

		if (sizeof(item_t) == 4)
		{
			if (kmer_len < max_suffix_len + 7)
				prefix_len = 7;
			else
				prefix_len = kmer_len - max_suffix_len;
		}
		else if (sizeof(item_t) == 8)
			prefix_len = 0;
		else
			prefix_len = 0;

		tab_size = 1ull << (prefix_len * 2);
		no_kmers = 0;
		no_inserts = 0;

		total_memory = 0;

		ht_desc = nullptr;

		init_ht();
	}

	// ************************************************************************************
	~CHT_kmer()
	{
		map<uint64_t, uint64_t> hist_ht_kmer;

	/*	for (uint64_t i = 0; i < tab_size; ++i)
			hist_ht_kmer[ht_desc[i].size] += 1;*/

		free_ht();
	}

	// ************************************************************************************
	void clear(const uint32_t idx)
	{
		ht_desc[idx].filled = 0;
		fill_n(ht_desc[idx].ht, ht_desc[idx].size + 1, EMPTY);
	}

	// ************************************************************************************
	void insert(const uint64_t kmer, CCounterIncrementer *cinc, const uint32_t add = 0)
	{
		uint32_t idx;
		uint64_t pos;

		tie(idx, pos) = _find(kmer, true);

		auto &cur_ht = ht_desc[idx].ht;

		if(add)		// restruct
			cur_ht[pos] += add;			
		else
		{
			auto cnt = get_count(cur_ht[pos]);

			if(cnt < counter_max && cinc->Increment(cnt) != cnt)
				++cur_ht[pos];
		}
	}

	// ************************************************************************************
	int count(const uint64_t kmer) 
	{
		uint32_t idx;
		uint64_t pos;

		tie(idx, pos) = _find(kmer, false);

		auto &cur_ht = ht_desc[idx].ht;

		if (cur_ht[pos] != EMPTY)
			return get_count(cur_ht[pos]);

		return 0;
	}

	// ************************************************************************************
	void prefetch_idx(const uint64_t kmer) const
	{
		_mm_prefetch((const char *) (ht_desc + get_idx(kmer)), _MM_HINT_T0);
	}

	// ************************************************************************************
	void prefetch_idx_kernel(const uint64_t kernel) const
	{
#ifdef KMER_MARGIN_2_SYMBOLS
		_mm_prefetch((const char *)(ht_desc + (kernel >> (2 * kmer_len - 8 - 2 * prefix_len))), _MM_HINT_T0);
#else
		_mm_prefetch((const char *)(ht_desc + (kernel >> (2 * kmer_len - 4 - 2 * prefix_len))), _MM_HINT_T0);
#endif
	}
	
	// ************************************************************************************
	void prefetch_data(const uint64_t kmer) const
	{
		uint32_t idx;
		item_t reminder;

		tie(idx, reminder) = split_kmer(kmer);

		uint64_t pos = hash(reminder, idx);

		_mm_prefetch((const char *) (ht_desc[idx].ht + pos), _MM_HINT_T0);
	}

	// ************************************************************************************
	void prefetch_data_kernel(const uint64_t kernel) const
	{
#ifdef KMER_MARGIN_2_SYMBOLS
		uint64_t kmer = kernel << (64 - (2 * kmer_len - 8) - 4);
#else
		uint64_t kmer = kernel << (64 - (2 * kmer_len - 4) - 2);
#endif
		uint32_t idx;
		item_t reminder;
		
		tie(idx, reminder) = split_kmer(kmer);

		uint64_t pos = hash(reminder, idx);

		_mm_prefetch((const char *)(ht_desc[idx].ht + pos), _MM_HINT_T0);
	}

	// ************************************************************************************
	bool find(const CKmer &kmer, stats_t &counts, CCounterIncrementer &cinc) const
	{
		if (kmer.is_full())
			return find_full(kmer, counts);
		else
			return find_partial(kmer, counts, cinc);
	}

	// ************************************************************************************
	int count(const CKmer &kmer) 
	{
		if (kmer.is_full())
			return count(kmer.data_normalized());
		else
		{
			assert(0);
			return 0;
		}
	}

	// ************************************************************************************
	uint64_t get_total_memory() const
	{
		return total_memory;
	}

	// ************************************************************************************
	uint64_t get_no_kmers() const
	{
		return no_kmers;
	}

	// ************************************************************************************
	uint64_t get_no_inserts() const
	{
		return no_inserts;
	}

	// ************************************************************************************
	uint32_t get_ht_prefix_len() const
	{
		return (uint32_t) prefix_len;
	}

	// ************************************************************************************
	double get_kmer_coverage() const
	{
		return 1.0;
		return (double)no_inserts / no_kmers;
	}
};

// ************************************************************************************
//
// ************************************************************************************
class CHT_pair_kmers {
	typedef pair<uint64_t, uint64_t> item_t;

	const double max_load_must = 0.6;
	const double max_load_need = 0.55;
	const double max_load_maybe = 0.5;

	typedef struct {
		item_t *ht;
		uint32_t size;
		uint32_t filled;
		uint32_t max_filled;
	} ht_desc_t;

	ht_desc_t *ht_desc;

	uint32_t kmer_size;
	uint64_t value_mask;
	uint32_t counter_shift;
	uint64_t max_counter;
	uint64_t no_parts;

	const item_t EMPTY{ 0ull, 0ull };
	const uint64_t HT_INIT_SIZE = 128ull << 10;

	void restruct(uint64_t id_part);
	void _clear(uint64_t id_part);
	uint64_t hash(uint64_t key) const
	{
		key ^= key >> 33;
		key *= 0xff51afd7ed558ccdL;
		key ^= key >> 33;
		key *= 0xc4ceb9fe1a85ec53L;
		key ^= key >> 33;

		return key;
	}

	void insert_raw(uint64_t key, uint64_t value_counter);

	uint64_t get_part_id_from_hash(const uint64_t hash_val) const
	{
		return modulo_divisor(hash_val >> 48, (int) no_parts);
	}

public:
	CHT_pair_kmers(uint32_t _kmer_size, uint64_t _no_parts);
	~CHT_pair_kmers();

	uint64_t get_part_id(const uint64_t key) const
	{
		return get_part_id_from_hash(hash(key));
	}

	void prefetch(uint64_t part_id, uint64_t key);
	void insert(uint64_t key, uint64_t value, uint64_t count);
	uint64_t count(uint64_t key, uint64_t value);
	void find(uint64_t key, vector<uint64_t> &v_value_counter);

	uint64_t get_value_mask() const
	{
		return value_mask;
	}

	uint32_t get_counter_shift() const
	{
		return counter_shift;
	}

	uint64_t get_no_pairs() const
	{
		uint64_t r = 0;

		for (uint64_t i = 0; i < no_parts; ++i)
			r += ht_desc[i].filled;

		return r;
	};

	bool need_restruct_soon() const
	{
		double r = 0.0;

		for (uint64_t i = 0; i < no_parts; ++i)
		{
			double ff = (double)ht_desc[i].filled / ht_desc[i].size;
			if (ff > r)
				r = ff;
		}

		return r > max_load_need;
	}

	void RestructIfNecessary(uint64_t id_part) 
	{
		if ((double) ht_desc[id_part].filled / ht_desc[id_part].size > max_load_maybe)
			restruct(id_part);
	}

	void clear()
	{
		for (uint32_t i = 0; i < no_parts; ++i)
			_clear(i);
	}
};

// EOF
