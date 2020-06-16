#pragma once
// *******************************************************************************************
// This file is a part of FQSqueezer software distributed under GNU GPL 3 licence.
// The homepage of the FQSqueezer project is http://sun.aei.polsl.pl/REFRESH/fqsqueezer
//
// Author: Sebastian Deorowicz
// Version: 1.1
// Date   : 2020-06-16
// *******************************************************************************************

#include "defs.h"
#include "utils.h"
#include "params.h"

#include <array>

typedef array<context_t, 15> context_levels_t;

class CCodeContext 
{
	friend class CDNACompressor;

	CParams *params;

public:
	const uint32_t no_ctx_levels = 7u;

	const uint64_t ctx_r_sym_mask = 0b11111111ull;

	// Definitions of parts of context for codes
	const uint32_t size_pos = 14;
	const uint32_t size_counts_level = 3;
	const uint32_t size_counts = 7;
	const uint32_t size_rsym = 4;
	const uint32_t size_let_max = 3;
	const uint32_t size_cor_zone = 3;

	const context_t mask_pos = (1ull << size_pos) - 1ull;
	const context_t mask_counts_level = (1ull << size_counts_level) - 1ull;
	const context_t mask_counts = (1ull << size_counts) - 1ull;
	const context_t mask_rsym = (1ull << size_rsym) - 1ull;
	const context_t mask_let_max = (1ull << size_let_max) - 1ull;
	const context_t mask_cor_zone = (1ull << size_cor_zone) - 1ull;

	const uint32_t shift_pos = 0;
	const uint32_t shift_counts_level = shift_pos + size_pos;
	const uint32_t shift_counts[4] = {
		shift_counts_level + size_counts_level,
		shift_counts_level + size_counts_level + size_counts,
		shift_counts_level + size_counts_level + 2 * size_counts,
		shift_counts_level + size_counts_level + 3 * size_counts };
	const uint32_t shift_rsym = shift_counts[3] + size_counts;
	const uint32_t shift_let_max = shift_rsym + size_rsym;
	const uint32_t shift_cor_zone = shift_let_max + size_let_max;

	const context_t enable_pos = mask_pos << shift_pos;
	const context_t enable_counts_level = mask_counts_level << shift_counts_level;
	const context_t enable_counts_all =
		(mask_counts << shift_counts[0]) +
		(mask_counts << shift_counts[1]) +
		(mask_counts << shift_counts[2]) +
		(mask_counts << shift_counts[3]);
	const context_t enable_counts[4] = {
		(mask_counts << shift_counts[0]),
		(mask_counts << shift_counts[1]),
		(mask_counts << shift_counts[2]),
		(mask_counts << shift_counts[3]) };
	const context_t enable_rsym = mask_rsym << shift_rsym;
	const context_t enable_let_max = mask_let_max << shift_let_max;
	const context_t enable_cor_zone = mask_cor_zone << shift_cor_zone;

	const uint32_t eor_size = 5;			// number of trailing symbols in read for which we store separate stats

private:
	int no_of_zeros(const stats_t &counts) const;
	int no_of_values(const stats_t &counts, uint32_t value) const;
	int no_of_more_than_ones(const stats_t &counts) const;
	int max_element(const stats_t &counts) const;

	uint64_t convert_count(uint64_t count, counts_level_t counts_level, uint32_t cnt_lev);
	uint64_t convert_lev_1_count(uint64_t count, uint32_t cnt_lev);
	uint64_t convert_lev_2_4_count(uint64_t count, uint32_t cnt_lev);
	uint64_t convert_lev_3_count(uint64_t count, uint32_t cnt_lev);

	uint32_t let_max_element(const stats_t &counts, const stats64_t &s_letters) const;

	context_t transform_r_sym(context_t r_sym) const;

public:
	CCodeContext(CParams *_params) : params(_params) {};

	void determine_ctx_codes(context_levels_t &a_code_ctx, const stats_t &counts, const stats64_t &s_letters,
		uint32_t pos, counts_level_t counts_level, uint32_t corrected_zone, context_t ctx_r_sym, uint32_t read_len);
	std::string decode_ctx(context_t ctx, int &lev);
};

class CLettersContext
{
	friend class CDNACompressor;

	CParams *params;

	context_t decode_letter(uint64_t letters, uint32_t idx) const;

public:
	const uint32_t no_ctx_levels = 9u;

	// Definitions of parts of context for codes
	const uint32_t size_pos = 14;
	const uint32_t size_letter = 3;

	const context_t mask_pos = (1ull << size_pos) - 1ull;
	const context_t mask_letter = (1ull << size_letter) - 1ull;
	
	const uint32_t shift_pos = 0;
	const uint32_t shift_letters[10] = {
		size_pos,
		size_pos + 1 * size_letter,
		size_pos + 2 * size_letter,
		size_pos + 3 * size_letter,
		size_pos + 4 * size_letter,
		size_pos + 5 * size_letter,
		size_pos + 6 * size_letter,
		size_pos + 7 * size_letter,
		size_pos + 8 * size_letter,
		size_pos + 9 * size_letter,
	};

	const context_t enable_pos = mask_pos << shift_pos;
	const context_t enable_letters[10] = {
		mask_letter << (size_pos),
		mask_letter << (size_pos + 1 * size_letter),
		mask_letter << (size_pos + 2 * size_letter),
		mask_letter << (size_pos + 3 * size_letter),
		mask_letter << (size_pos + 4 * size_letter),
		mask_letter << (size_pos + 5 * size_letter),
		mask_letter << (size_pos + 6 * size_letter),
		mask_letter << (size_pos + 7 * size_letter),
		mask_letter << (size_pos + 8 * size_letter),
		mask_letter << (size_pos + 9 * size_letter),
	};

	const uint32_t eor_size = 5;			// number of trailing symbols in read for which we store separate stats

private:
public:
	CLettersContext(CParams *_params) : params(_params) {};

	void determine_ctx_letters(context_levels_t &a_letter_ctx, uint32_t pos, uint64_t letters, uint32_t read_len);

	std::string decode_ctx(context_t ctx);
};

// EOF