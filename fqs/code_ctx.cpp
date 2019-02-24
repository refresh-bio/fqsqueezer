// *******************************************************************************************
// This file is a part of FQSqueezer software distributed under GNU GPL 3 licence.
// The homepage of the MSAC project is http://sun.aei.polsl.pl/REFRESH/fqsqueezer
//
// Author: Sebastian Deorowicz
// Version: 1.0
// Date   : 2019-02-22
// *******************************************************************************************

#include "code_ctx.h"

#include <iostream>

//*****************************************************************************************************
uint64_t CCodeContext::convert_count(uint64_t count, counts_level_t counts_level, uint32_t cnt_lev)
{
	if (counts_level == counts_level_t::pmer)
		return convert_lev_1_count(count, cnt_lev);
	else if (counts_level == counts_level_t::bmer)
		return convert_lev_3_count(count, cnt_lev);
	else
		return convert_lev_2_4_count(count, cnt_lev);
}

//*****************************************************************************************************
uint64_t CCodeContext::convert_lev_1_count(uint64_t count, uint32_t cnt_lev)
{
	uint64_t flag = (uint64_t) cnt_lev << 5;

/*	if (cnt_lev == 0)
	{
		if (count < 3)	return flag + count;
		if (count < 5)	return flag + 3;
		return flag + 4;
	}*/

	if (cnt_lev == 0)
	{
		if (count < 5)	return flag + count;
		if (count < 8)	return flag + 5;
		if (count < 16)	return flag + 6;
		if (count < 32)	return flag + 7;
		if (count < 64)	return flag + 8;
		return flag + 9;
	}

/*	if (cnt_lev == 2 || cnt_lev == 1)
	{
		if (count < 8)	return flag + count;
		if (count < 16)	return flag + 8;
		if (count < 24)	return flag + 9;
		if (count < 32)	return flag + 10;
		if (count < 48)	return flag + 11;
		if (count < 64)	return flag + 12;
		if (count < 96)	return flag + 13;
		return flag + 14;
	}*/

	if (cnt_lev == 3 || cnt_lev == 1 || cnt_lev == 2)
	{
		if (count < 8)	return flag + count;
		if (count < 16)	return flag + 8;
		if (count < 24)	return flag + 9;
		if (count < 32)	return flag + 10;
		if (count < 40)	return flag + 11;
		if (count < 48)	return flag + 12;
		if (count < 56)	return flag + 13;
		if (count < 64)	return flag + 14;
		if (count < 80)	return flag + 15;
		if (count < 96)	return flag + 16;
		if (count < 112)	return flag + 17;
		if (count < 128)	return flag + 18;
		if (count < 144)	return flag + 19;
		if (count < 160)	return flag + 20;
		if (count < 176)	return flag + 21;
		if (count < 192)	return flag + 22;
		if (count < 224)	return flag + 23;
		if (count < 288)	return flag + 24;
		if (count < 384)	return flag + 25;
		if (count < 512)	return flag + 26;
		if (count < 1024)	return flag + 27;
		if (count < 2048)	return flag + 28;
		return flag + 29;
	}

	return 0;
}

//*****************************************************************************************************
uint64_t CCodeContext::convert_lev_2_4_count(uint64_t count, uint32_t cnt_lev)
{
	uint64_t flag = (uint64_t) cnt_lev << 5;

	if (cnt_lev == 0)
	{
		if (count < 3)	return flag + count;
		if (count < 5)	return flag + 3;
		return flag + 4;
	}

	if (cnt_lev == 1)
	{
		if (count < 5)	return flag + count;
		if (count < 8)	return flag + 5;
		if (count < 13)	return flag + 6;
		if (count < 20)	return flag + 7;
		if (count < 30)	return flag + 8;
		return flag + 9;
	}

	if (cnt_lev == 2)
	{
		if (count < 10)	return flag + count;
		if (count < 15)	return flag + 10;
		if (count < 20)	return flag + 11;
		if (count < 30)	return flag + 12;
		if (count < 50)	return flag + 13;
		return flag + 14;
	}

/*	if (cnt_lev == 3)
	{
		if (count < 10)	return flag + count;
		if (count < 16)	return flag + 10;
		if (count < 24)	return flag + 11;
		if (count < 32)	return flag + 12;
		if (count < 48)	return flag + 13;
		if (count < 64)	return flag + 14;
		if (count < 128)	return flag + 15;
		if (count < 256)	return flag + 16;
		if (count < 1024)	return flag + 17;
		if (count < 4096)	return flag + 18;
		if (count < 4096 * 2)	return flag + 19;
		if (count < 4096 * 4)	return flag + 20;
		if (count < 4096 * 8)	return flag + 21;
		if (count < 4096 * 8 + 64)	return flag + 22;
		if (count < 4096 * 8 + 128)	return flag + 23;
		if (count < 4096 * 8 + 256)	return flag + 24;
		if (count < 4096 * 8 + 512)	return flag + 25;
		if (count < 4096 * 8 + 1024)	return flag + 26;
		if (count < 4096 * 8 + 2048)	return flag + 27;
		if (count < 4096 * 8 + 4096)	return flag + 28;
		return flag + 29;
	}*/

	if (cnt_lev == 3)
	{
		if (count < 10)	return flag + count;
		if (count < 16)	return flag + 10;
		if (count < 24)	return flag + 11;
		if (count < 32)	return flag + 12;
		if (count < 48)	return flag + 13;
		if (count < 64)	return flag + 14;
		if (count < 128)	return flag + 15;
		if (count < 256)	return flag + 16;
		if (count < 512)	return flag + 17;
		if (count < 1024)	return flag + 18;
		if (count < 2048)	return flag + 19;
		if (count < 2048 + 32)	return flag + 20;
		if (count < 2048 + 64)	return flag + 21;
		if (count < 2048 + 128)	return flag + 22;
		if (count < 2048 + 192)	return flag + 23;
		if (count < 2048 + 256)	return flag + 24;
		if (count < 2048 + 384)	return flag + 25;
		if (count < 2048 + 512)	return flag + 26;
		if (count < 2048 + 768)	return flag + 27;
		if (count < 2048 + 1024)	return flag + 28;
		return flag + 29;
	}

	return 0;
}

//*****************************************************************************************************
uint64_t CCodeContext::convert_lev_3_count(uint64_t count, uint32_t cnt_lev)
{
	uint64_t flag = (uint64_t) cnt_lev << 5;

	if (cnt_lev == 0)
	{
		if (count < 3)	return flag + count;
		if (count < 5)	return flag + 3;
		return flag + 4;
	}

	if (cnt_lev == 1)
	{
		if (count < 5)	return flag + count;
		if (count < 8)	return flag + 5;
		if (count < 13)	return flag + 6;
		if (count < 20)	return flag + 7;
		if (count < 30)	return flag + 8;
		return flag + 9;
	}

	if (cnt_lev == 2)
	{
		if (count < 10)	return flag + count;
		if (count < 13)	return flag + 10;
		if (count < 20)	return flag + 11;
		if (count < 30)	return flag + 12;
		if (count < 50)	return flag + 13;
		return flag + 14;
	}

	if (cnt_lev == 3)
	{
/*		if (count < 16)	return flag + count;
		if (count < 18)	return flag + 16;
		if (count < 20)	return flag + 17;
		if (count < 24)	return flag + 18;
		if (count < 28)	return flag + 19;
		if (count < 32)	return flag + 21;
		if (count < 38)	return flag + 22;
		if (count < 44)	return flag + 23;
		if (count < 50)	return flag + 24;
		if (count < 56)	return flag + 25;
		if (count < 64)	return flag + 26;
		return flag + 27;*/
		if (count < 15)	return flag + count;
		if (count < 18)	return flag + 15;
		if (count < 20)	return flag + 16;
		if (count < 25)	return flag + 17;
		if (count < 30)	return flag + 18;
		if (count < 40)	return flag + 19;
		if (count < 50)	return flag + 20;
		if (count < 60)	return flag + 21;
		if (count < 64)	return flag + 22;
		if (count < 68)	return flag + 23;
		if (count < 72)	return flag + 24;
		if (count < 76)	return flag + 25;
		if (count < 80)	return flag + 26;
		if (count < 84)	return flag + 27;
		if (count < 88)	return flag + 28;
		return flag + 29;
	}

	return 0;
}

//*****************************************************************************************************
// Prepare s-contexts:
//  0 - nothing
//  1 + counts_level
//      counts[0] and counts[1] at 10 levels
//      counts[2] and counts[3] at 5 levels
//	    round(pos, 8)
//  2 + cor_zone
//      r_sym
//  3 + counts[0] and counts[1] at 15 levels
//  4 + counts[0] and counts[1] at 20 levels
//      counts[2] and counts[3] at 10 levels
//  5 + let_max
//  6 + round(pos, 4)

//*****************************************************************************************************
void CCodeContext::determine_ctx_codes(context_levels_t &a_code_ctx, const stats_t &counts, const stats64_t &s_letters,
	uint32_t pos, counts_level_t counts_level, uint32_t corrected_zone, context_t ctx_r_sym, uint32_t read_len)
{
	context_t mask = ~(context_t)0ull;
	context_t ctx = 0ull;

	uint32_t pos_limit[] = { 0, params->pmer_len, params->smer_len, params->bmer_len, params->bmer_len, params->bmer_len };


	stats_t counts_srt;
	sort_copy_stats(counts_srt, counts);

	// level 0
	a_code_ctx[0] = ctx | mask;

	// level 1
	mask ^= enable_counts_level | enable_counts_all | enable_pos;
	ctx += ((context_t)counts_level) << shift_counts_level;

	for (int i = 0; i < 2; ++i)
		ctx += ((context_t)convert_count(counts_srt[i], counts_level, 1)) << shift_counts[i];
	for (int i = 2; i < 4; ++i)
		ctx += ((context_t)convert_count(counts_srt[i], counts_level, 0)) << shift_counts[i];

	if (pos < pos_limit[(uint32_t) counts_level])
		ctx += ((context_t)pos) << shift_pos;
	else if (pos + eor_size >= read_len)
		ctx += ((context_t)(mask_pos - (read_len - pos))) << shift_pos;
	else
		ctx += ((context_t)pos_limit[(uint32_t) counts_level] + pos / 16) << shift_pos;

	a_code_ctx[1] = ctx | mask;

	// level 2
	mask ^= enable_cor_zone | enable_rsym;
	ctx += ((context_t)corrected_zone) << shift_cor_zone;
	ctx += ((context_t)transform_r_sym(ctx_r_sym)) << shift_rsym;
	a_code_ctx[2] = ctx | mask;

	// level 3
	ctx &= ~(enable_counts[0] | enable_counts[1]);
	ctx += ((context_t)convert_count(counts_srt[0], counts_level, 2)) << shift_counts[0];
	ctx += ((context_t)convert_count(counts_srt[1], counts_level, 2)) << shift_counts[1];
	a_code_ctx[3] = ctx | mask;

	// level 4
	ctx &= ~(enable_counts[0] | enable_counts[1] | enable_counts[2] | enable_counts[3]);
	ctx += ((context_t)convert_count(counts_srt[0], counts_level, 3)) << shift_counts[0];
	ctx += ((context_t)convert_count(counts_srt[1], counts_level, 3)) << shift_counts[1];
	ctx += ((context_t)convert_count(counts_srt[2], counts_level, 1)) << shift_counts[2];
	ctx += ((context_t)convert_count(counts_srt[3], counts_level, 1)) << shift_counts[3];
	a_code_ctx[4] = ctx | mask;

	// level 5
	mask ^= enable_let_max;
	ctx += ((context_t)let_max_element(counts, s_letters)) << shift_let_max;
	a_code_ctx[5] = ctx | mask;

	// level 6
	ctx &= ~enable_pos;
	if (pos < pos_limit[(uint32_t) counts_level])
		ctx += ((context_t)pos + (1u << 13)) << shift_pos;
	else if (pos + eor_size >= read_len)
		ctx += ((context_t)(mask_pos - (read_len - pos))) << shift_pos;
	else
		ctx += ((context_t)pos_limit[(uint32_t) counts_level] + pos / 8 + (1u << 13)) << shift_pos;
	a_code_ctx[6] = ctx | mask;
}

//*****************************************************************************************************
uint32_t CCodeContext::let_max_element(const stats_t &counts, const stats64_t &s_letters) const
{
	uint32_t r = 0;

	for (int i = 1; i < 4; ++i)
		if (counts[i] > counts[r])
			r = i;
		else if (counts[i] == counts[r] && s_letters[i] > s_letters[r])
			r = i;

	return r;
}

//*****************************************************************************************************
int CCodeContext::max_element(const stats_t &counts) const
{
	uint32_t r = 0;

	for (int i = 1; i < 4; ++i)
		if (counts[i] > counts[r])
			r = i;

	return r;
}

//*****************************************************************************************************
int CCodeContext::no_of_zeros(const stats_t &counts) const
{
	return (counts[0] == 0) + (counts[1] == 0) + (counts[2] == 0) + (counts[3] == 0);
}

//*****************************************************************************************************
int CCodeContext::no_of_more_than_ones(const stats_t &counts) const
{
	return (counts[0] > 1) + (counts[1] > 1) + (counts[2] > 1) + (counts[3] > 1);
}

//*****************************************************************************************************
int CCodeContext::no_of_values(const stats_t &counts, uint32_t value) const
{
	return (counts[0] == value) + (counts[1] == value) + (counts[2] == value) + (counts[3] == value);
}

//*****************************************************************************************************
context_t CCodeContext::transform_r_sym(context_t r_sym) const
{
	return popcnt(r_sym);

	uint32_t r = popcnt(r_sym);

	uint32_t x = (~r_sym) & ((~r_sym) - 1u);

	x &= 0xff;

	r = r + 8 - popcnt(x ^ r_sym);

	return min(r, 14u);
}

//*****************************************************************************************************
string CCodeContext::decode_ctx(context_t ctx, int &lev)
{
	string s = "[";
	context_t t;

	t = (ctx >> shift_counts_level) & mask_counts_level;
	lev = (int) t;

	if (t != mask_counts_level)
		s += "lev:" + to_string(t);
	else
		s += "lev:---";

	t = (ctx >> shift_pos) & mask_pos;
	if (t != mask_pos)
	{
		s += "  pos:";
		if (t > (1u << 13))
			s += "1:";
		else
			s += "0:";

		s += to_string(t & 0x1fff);
	}
	else
		s += "  pos:---";

	t = (ctx >> shift_cor_zone) & mask_cor_zone;
	if (t != mask_cor_zone)
		s += "  cz:" + to_string(t);
	else
		s += "  cz:---";
	
	t = (ctx >> shift_rsym) & mask_rsym;
	if (t != mask_rsym)
		s += "  r_sym:" + to_string(t);
	else
		s += "  r_sym:---";

	t = (ctx >> shift_let_max) & mask_let_max;
	if (t != mask_let_max)
		s += "  let_mc:" + to_string(t);
	else
		s += "  let_mc:---";

	t = (ctx >> shift_counts[0]) & mask_counts;
	if (t != mask_counts)
		s += string("  C0:") + ("0123"[t >> 5]) + ":" + to_string(t & 31);
	else
		s += "  C0:---";

	t = (ctx >> shift_counts[1]) & mask_counts;
	if (t != mask_counts)
		s += string("  C1:") + ("0123"[t >> 5]) + ":" + to_string(t & 31);
	else
		s += "  C1:---";

	t = (ctx >> shift_counts[2]) & mask_counts;
	if (t != mask_counts)
		s += string("  C2:") + ("0123"[t >> 5]) + ":" + to_string(t & 31);
	else
		s += "  C2:---";

	t = (ctx >> shift_counts[3]) & mask_counts;
	if (t != mask_counts)
		s += string("  C3:") + ("0123"[t >> 5]) + ":" + to_string(t & 31);
	else
		s += "  C3:---";

	s += "]";

	return s;
}


//*****************************************************************************************************
//
//*****************************************************************************************************
void CLettersContext::determine_ctx_letters(context_levels_t &a_letter_ctx, uint32_t pos, uint64_t letters, uint32_t read_len)
{
	context_t mask = ~(context_t)0u;
	context_t ctx = 0u;

	mask ^= enable_pos;
	
	if (pos < params->pmer_len)
		ctx += ((context_t)pos) << shift_pos;
	else if (pos + eor_size > read_len)
		ctx += ((context_t)(mask_pos - (read_len - pos))) << shift_pos;

	for (uint32_t i = 0; i < 10; ++i)
	{
		ctx += decode_letter(letters, i) << shift_letters[i];
		mask ^= enable_letters[i];

		a_letter_ctx[i] = ctx | mask;
	}
}

//*****************************************************************************************************
context_t CLettersContext::decode_letter(uint64_t letters, uint32_t idx) const
{
	return (letters >> (4 * idx)) & 0x7ull;
}

//*****************************************************************************************************
string CLettersContext::decode_ctx(context_t ctx)
{
	string s = "[";
	context_t t;

	t = (ctx >> shift_pos) & mask_pos;

	if (t != mask_pos)
		s += "pos:" + to_string(t);
	else
		s += "pos:---";

	s += "  letters:";

	for (int32_t i = 9; i >= 0; --i)
	{
		t = (ctx >> shift_letters[i]) & mask_letter;
		s.push_back("ACGTN**-"[t]);
	}

	s += "]";

	return s;
}

// EOF