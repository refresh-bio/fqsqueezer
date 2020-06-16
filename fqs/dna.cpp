// *******************************************************************************************
// This file is a part of FQSqueezer software distributed under GNU GPL 3 licence.
// The homepage of the FQSqueezer project is http://sun.aei.polsl.pl/REFRESH/fqsqueezer
//
// Author: Sebastian Deorowicz
// Version: 1.1
// Date   : 2020-06-16
// *******************************************************************************************

#include "dna.h"
#include "utils.h"

//*****************************************************************************************************
CDNACompressor::CDNACompressor() : CBasicCompressor()
{
	no_threads = 1;

	for (int i = 0; i < 256; ++i)
		dna_convert[i] = 4;
	dna_convert['A'] = dna_code('A');
	dna_convert['C'] = dna_code('C');
	dna_convert['G'] = dna_code('G');
	dna_convert['T'] = dna_code('T');

	tpl_ctx_rc_flags = nullptr;
	tpl_ctx_rc_letters = nullptr;
	tpl_ctx_rc_codes = nullptr;

	tpl_ctx_rc_prefix_sorted_flags = nullptr;	
	tpl_ctx_rc_prefix_sorted_nibbles = nullptr;
	tpl_ctx_rc_prefix_sorted_no_bytes = nullptr;
	tpl_ctx_rc_prefix_sorted_bytes = nullptr;
	tpl_ctx_rc_prefix_Ns = nullptr;

	ctx_rc_pe_minimizer_id = nullptr;
	tpl_ctx_rc_minimizer_pos = nullptr;
	
	siv_prune_thr = 7.0;

	ht_bmer_local = nullptr;
	ht_smer_local = nullptr;
	ht_pe_mers_local = nullptr;
	pmer_mod_shift = 0;

	fill_n(a_stat_minim_id.begin(), a_stat_minim_id.size(), 0u);

	ctx_flags = 0ull;

	code_ctx = new CCodeContext(&params);
	letters_ctx = new CLettersContext(&params);
}

//*****************************************************************************************************
CDNACompressor::~CDNACompressor()
{
	if (thread_id == 0 && params.verbosity == 2)
		show_contexts();

	if (tpl_ctx_rc_flags)
		delete tpl_ctx_rc_flags;
	if (tpl_ctx_rc_letters)
		delete tpl_ctx_rc_letters;
	if (tpl_ctx_rc_codes)
		delete tpl_ctx_rc_codes;

	if (tpl_ctx_rc_prefix_sorted_flags)
		delete tpl_ctx_rc_prefix_sorted_flags;
	if (tpl_ctx_rc_prefix_sorted_no_bytes)
		delete tpl_ctx_rc_prefix_sorted_no_bytes;
	if (tpl_ctx_rc_prefix_sorted_nibbles)
		delete tpl_ctx_rc_prefix_sorted_nibbles;
	if (tpl_ctx_rc_prefix_sorted_bytes)
		delete tpl_ctx_rc_prefix_sorted_bytes;
	if (tpl_ctx_rc_prefix_Ns)
		delete tpl_ctx_rc_prefix_Ns;

	if (ctx_rc_pe_minimizer_id)
		delete ctx_rc_pe_minimizer_id;
	if (tpl_ctx_rc_minimizer_pos)
		delete tpl_ctx_rc_minimizer_pos;

	if (ht_bmer_local)
		delete ht_bmer_local;
	if (ht_smer_local)
		delete ht_smer_local;

	if (ht_pe_mers_local)
		delete ht_pe_mers_local;

	delete code_ctx;
	delete letters_ctx;
}

//*****************************************************************************************************
bool CDNACompressor::Init()
{
	if (ht_bmer_local)
		delete ht_bmer_local;
	ht_bmer_local = new CHT_kmer<uint64_t>(params.bmer_len, BMER_COUNTER_BITS, params.ht_max_filling_factor);

	if (ht_smer_local)
		delete ht_smer_local;
	ht_smer_local = new CHT_kmer<uint64_t>(params.smer_len, SMER_COUNTER_BITS, params.ht_max_filling_factor);

	if (ht_pe_mers_local)
		delete ht_pe_mers_local;
	ht_pe_mers_local = new CHT_pair_kmers(params.bmer_len, 1);

	if (tpl_ctx_rc_flags)
		delete tpl_ctx_rc_flags;
	tpl_ctx_rc_flags = new CRangeCoderModelFixedSize<CVectorIOStream, 2>(rc, 12, 1 << 12, nullptr, 4, params.work_mode == work_mode_t::compress);

	if (tpl_ctx_rc_letters)
		delete tpl_ctx_rc_letters;
	int init_letters_stats[] = { 10, 10, 10, 10, 1 };
	tpl_ctx_rc_letters = new CRangeCoderModelFixedSize<CVectorIOStream, 5>(rc, 15, 1 << 15, init_letters_stats, 4, params.work_mode == work_mode_t::compress);

	if (tpl_ctx_rc_codes)
		delete tpl_ctx_rc_codes;
	int init_codes_stats[] = { 20, 6, 3, 2, 1 };
	tpl_ctx_rc_codes = new CRangeCoderModelFixedSize<CVectorIOStream, 5>(rc, 15, 1 << 15, init_codes_stats, 4, params.work_mode == work_mode_t::compress);

	max_prefix_sorted_flag_value = (1u << SIV_FIELD_SIZE);
	if (tpl_ctx_rc_prefix_sorted_flags)
		delete tpl_ctx_rc_prefix_sorted_flags;
	tpl_ctx_rc_prefix_sorted_flags = new CRangeCoderModel<CVectorIOStream>(rc, 1 + max_prefix_sorted_flag_value, 12, 1 << 12, nullptr, 4, params.work_mode == work_mode_t::compress);

	max_prefix_sorted_no_bytes_value = (2 * params.pmer_len + 7) / 8;
	if (tpl_ctx_rc_prefix_sorted_no_bytes)
		delete tpl_ctx_rc_prefix_sorted_no_bytes;
	tpl_ctx_rc_prefix_sorted_no_bytes = new CRangeCoderModel<CVectorIOStream>(rc, max_prefix_sorted_no_bytes_value, 12, 1 << 12, nullptr, 4, params.work_mode == work_mode_t::compress);

	if (tpl_ctx_rc_prefix_sorted_nibbles)
		delete tpl_ctx_rc_prefix_sorted_nibbles;
	tpl_ctx_rc_prefix_sorted_nibbles = new CRangeCoderModelFixedSize<CVectorIOStream, 16>(rc, 15, 1 << 15, nullptr, 4, params.work_mode == work_mode_t::compress);

	if (tpl_ctx_rc_prefix_sorted_bytes)
		delete tpl_ctx_rc_prefix_sorted_bytes;
	tpl_ctx_rc_prefix_sorted_bytes = new CRangeCoderModelFixedSize<CVectorIOStream, 256>(rc, 15, 1 << 15, nullptr, 4, params.work_mode == work_mode_t::compress);

	if (tpl_ctx_rc_prefix_Ns)
		delete tpl_ctx_rc_prefix_Ns;
	tpl_ctx_rc_prefix_Ns = new CRangeCoderModelFixedSize<CVectorIOStream, 2>(rc, 12, 1 << 12, nullptr, 4, params.work_mode == work_mode_t::compress);

	ctx_rc_pe_minimizer_id = new CRangeCoderModelFixedSize<CVectorIOStream, 16>(rc, 15, 1 << 15, nullptr, 4, params.work_mode == work_mode_t::compress);
	tpl_ctx_rc_minimizer_pos = new CRangeCoderModel<CVectorIOStream>(rc, 256, 15, 1 << 15, nullptr, 4, params.work_mode == work_mode_t::compress);

	reset_context(ctx_flags);
	reset_context(ctx_letters);

	reset_context(ctx_prefix_sorted_flags);

	pmer_can.Reset(params.pmer_len, kmer_mode_t::canonical);
	smer_can.Reset(params.smer_len, kmer_mode_t::canonical);
	bmer_can.Reset(params.bmer_len, kmer_mode_t::canonical);

	pmer_can_unc.Reset(params.pmer_len, kmer_mode_t::canonical);
	smer_can_unc.Reset(params.smer_len, kmer_mode_t::canonical);
	bmer_can_unc.Reset(params.bmer_len, kmer_mode_t::canonical);

	mt.seed(5481);
	cinc_b.Reset(7, 2, BMER_MAX_COUNTER_VALUE);
	cinc_s.Reset(SMER_MAX_COUNTER_VALUE / 2, 1, SMER_MAX_COUNTER_VALUE);
	cinc_lb.Reset(7, 2, BMER_MAX_COUNTER_VALUE);
	cinc_ls.Reset(SMER_MAX_COUNTER_VALUE / 2, 1, SMER_MAX_COUNTER_VALUE);

	pmer_can_prev = pmer_can;

	fill(s_letters.begin(), s_letters.end(), 0ull);

	n_rc_ctx = 0;

	return true;
}

//*****************************************************************************************************
uint8_t CDNACompressor::rank(const stats_t &counts, uint8_t sym) const
{
	if (sym == 4)
		return 4;

	uint8_t r = 0;

	for (int i = 0; i < 4; ++i)
		if (counts[i] != counts[sym])
			r += counts[sym] < counts[i];
		else if (s_letters[i] != s_letters[sym])
			r += s_letters[sym] < s_letters[i];
		else
			r += sym > i;

	return r;
}

//*****************************************************************************************************
// !!! Can be improved
uint8_t CDNACompressor::un_rank(stats_t &counts, uint8_t sym)
{
	if (sym == 4)
		return 4;

	for (int i = 0; i < 4; ++i)
		if (rank(counts, i) == sym)
			return i;

	return 4;
}

//*****************************************************************************************************
bool CDNACompressor::find_counts_p(const CKmer &kmer, stats_t &counts)
{
	if (!kmer.is_full())
	{
		CKmer pmer_tmp = kmer;

		for (int j = 0; j < 4; ++j)
		{
			pmer_tmp.replace_last(j);
			counts[j] = (uint32_t)siv_pmer->test_shorter(pmer_tmp.data_aligned_rc(), 2 * (uint64_t) pmer_tmp.get_cur_size());
		}
	}
	else
		siv_pmer->counts(kmer.data_aligned_dir(), counts);

	return non_empty(counts);
}

//*****************************************************************************************************
bool CDNACompressor::find_counts_rough_p(stats_t &counts)
{
	clear_counts(counts);

	CKmer pmer_tmp, pmer_pf;

	for (uint32_t i = 0; i < params.pmer_len - 1; ++i)
	{
		pmer_tmp.ResetFromCan(pmer_can, kmer_mode_t::direct);
		pmer_pf = pmer_tmp;

		for (int j = 0; j < 4; ++j)
		{
			if (i + 2 < params.pmer_len - 1)
			{
				pmer_pf.replace(j, i + 2);
				siv_pmer->prefetch(pmer_pf.data_aligned());
			}

			pmer_tmp.replace(j, i);
			siv_pmer->counts_update(pmer_tmp.data_aligned(), counts);
		}
	}

	return non_empty(counts);
}

//*****************************************************************************************************
bool CDNACompressor::find_counts_rough_s(stats_t &counts)
{
	clear_counts(counts);

	CKmer smer_tmp, smer_pf1, smer_pf2;
	stats_t loc_counts;

	for (uint32_t i = 0; i < params.smer_len - 1; ++i)
	{
		smer_tmp = smer_can;
		smer_pf1 = smer_tmp;
		smer_pf2 = smer_tmp;

		for (int j = 0; j < 4; ++j)
		{
			if (i + 2 < params.smer_len - 1)
			{
				smer_pf1.replace(j, i + 2);
				ht_smer->prefetch_idx(smer_pf1.kernel_canonical());
			}
			if (i + 1 < params.smer_len - 1)
			{
				smer_pf2.replace(j, i + 1);
				ht_smer->prefetch_data(smer_pf2.kernel_canonical());
			}

			smer_tmp.replace(j, i);
			
			if (ht_smer->find(smer_tmp, loc_counts, cinc_s))
				for (int i = 0; i < 4; ++i)
					counts[i] = cinc_s.Increment(counts[i], loc_counts[i]);			
		}
	}

	return non_empty(counts);
}

//*****************************************************************************************************
bool CDNACompressor::find_counts_rough_b(stats_t &counts)
{
	clear_counts(counts);

	CKmer bmer_tmp, bmer_pf1, bmer_pf2;
	stats_t loc_counts;

	for (uint32_t i = 0; i < params.bmer_len - 1; ++i)
	{
		bmer_tmp = bmer_can;
		bmer_pf1 = bmer_tmp;
		bmer_pf2 = bmer_tmp;

		for (int j = 0; j < 4; ++j)
		{
			if (i + 2 < params.bmer_len - 1)
			{
				bmer_pf1.replace(j, i + 2);
				ht_bmer->prefetch_idx(bmer_pf1.kernel_canonical());
			}
			if (i + 1 < params.bmer_len - 1)
			{
				bmer_pf2.replace(j, i + 1);
				ht_bmer->prefetch_data(bmer_pf2.kernel_canonical());
			}

			bmer_tmp.replace(j, i);

			if (ht_bmer->find(bmer_tmp, loc_counts, cinc_b))
				for (int i = 0; i < 4; ++i)
					counts[i] = cinc_b.Increment(counts[i], loc_counts[i]);
		}
	}

	return non_empty(counts);
}

//*****************************************************************************************************
bool CDNACompressor::repair_kmers_existing(uint8_t *p_cur, uint32_t pos, stats_t &counts, uint8_t sym, uint32_t prev_count)
{
	uint32_t max_sym = 0;
	for (int i = 1; i < 4; ++i)
		if (counts[i] > counts[max_sym])
			max_sym = i;
		else if(counts[i] == counts[max_sym] && s_letters[i] > s_letters[max_sym])
			max_sym = i;

	uint32_t rpl_thr[][2] = { { 0, 3 }};

	if (sym != 4)
	{
		if (max_sym == sym)
			return false;

		bool was_verified = false;
		for (uint64_t i = 0; i < sizeof(rpl_thr) / sizeof(rpl_thr[0]); ++i)
			if (counts[sym] == rpl_thr[i][0])
			{
				if (counts[max_sym] <= rpl_thr[i][1])
					return false;
				was_verified = true;
				break;
			}

		if (!was_verified)
			return false;
	}

	pmer_can.replace_last(max_sym);
	smer_can.replace_last(max_sym);
	bmer_can.replace_last(max_sym);

	cor_pos = pos;

	return true;
}

//*****************************************************************************************************
// Try correction of b-mer
bool CDNACompressor::repair_kmers_missing(uint8_t *p_cur, uint32_t pos, uint8_t sym, uint32_t existing_count)
{
	if (siv_pmer->avg_filling_factor() < 7.0)
		return false;

	int best_c = 4;
	int best_count = 0;
	int best_pos = 0;
	int max_dist = 6;
//	int max_dist = 4;

	uint32_t min_count = 2;
//	uint32_t min_count = 1;

	if (existing_count > 2)
		min_count = existing_count;

	for (int j = 1; j < max_dist; ++j)
	{
		CKmer bmer_tmp = bmer_can;
		CKmer bmer_tmp_p1 = bmer_can;
		CKmer bmer_tmp_p2 = bmer_can;

		bool is_pos_ok = false;

		for (uint32_t c = 0; c < 4; ++c)
		{
			if (j + 2 < max_dist)
			{
				if (bmer_can.get_symbol(bmer_can.get_cur_size() - 1 - (j + 2)) != c)
				{
					bmer_tmp_p1.replace(c, bmer_tmp_p1.get_cur_size() - 1 - (j + 2));
					ht_bmer->prefetch_idx(bmer_tmp_p1.data_normalized());
				}
			}

			if (j + 1 < max_dist)
			{
				if (bmer_can.get_symbol(bmer_can.get_cur_size() - 1 - (j + 1)) != c)
				{
					bmer_tmp_p2.replace(c, bmer_tmp_p2.get_cur_size() - 1 - (j + 1));
					ht_bmer->prefetch_data(bmer_tmp_p2.data_normalized());
				}
			}

			if (bmer_can.get_symbol(bmer_can.get_cur_size() - 1 - j) == c)
				continue;

			bmer_tmp.replace(c, bmer_tmp.get_cur_size() - 1 - j);
			int count = ht_bmer->count(bmer_tmp);

			if (count)
			{
				if (!is_pos_ok)
					is_pos_ok = true;
			}

			if (count >= best_count && count >= (int) min_count)
			{
				best_c = c;
				best_count = count;
				best_pos = j;
			}
		}
	}

	if (best_pos)
	{
		bmer_can.replace(best_c, bmer_can.get_cur_size() - 1 - best_pos);
		if(best_pos < (int) smer_can.get_cur_size())
			smer_can.replace(best_c, smer_can.get_cur_size() - 1 - best_pos);
		if (best_pos < (int) pmer_can.get_cur_size())
			pmer_can.replace(best_c, pmer_can.get_cur_size() - 1 - best_pos);

		cor_pos = max(cor_pos, pos - best_pos);

		return true;
	}

	return false;
}

//*****************************************************************************************************
counts_level_t CDNACompressor::find_counts(stats_t &counts)
{
	clear_counts(counts);

	uint32_t bmer_margin = params.bmer_len - params.smer_len - 1;
	uint32_t smer_margin = params.smer_len - params.pmer_len;

	smer_margin += 1;

	if (bmer_can.is_almost_full(bmer_margin))
	{
		if (find_counts_b(bmer_can, counts))
		{
			if (code_ctx->no_of_values(counts, BMER_MAX_COUNTER_VALUE) > 1)
			{
				stats_t counts2;
				find_counts_s(smer_can, counts2);

				for (int i = 0; i < 4; ++i)
					counts[i] += counts2[i];

				return counts_level_t::mixed;
			}

			return counts_level_t::bmer;
		}
		else
		{
			if (find_counts_b_local(bmer_can, counts))
				return counts_level_t::bmer;

			if (bmer_can != bmer_can_unc && find_counts_b(bmer_can_unc, counts))
				return counts_level_t::bmer_unc;
		}
	}

	if (smer_can.is_almost_full(smer_margin))
	{
		if (find_counts_s(smer_can, counts) || find_counts_s_local(smer_can, counts))
			return counts_level_t::smer;
	}
	else if (find_counts_p(pmer_can, counts))
		return counts_level_t::pmer;

	return counts_level_t::none;
}

//*****************************************************************************************************
// Direct encoding of a prefix
void CDNACompressor::compress_prefix_direct(const uint8_t *p, uint32_t size, uint32_t start_pos)
{
	uint8_t sym;
	reset_ctx_letters(ctx_letters);

	if(start_pos)
		for (uint32_t i = 0; i < start_pos; ++i)
		{
			sym = dna_convert[p[i]];
			update_ctx_letters(ctx_letters, sym);
		}

	for (uint32_t i = start_pos; i < params.prefix_len; ++i)
	{
		sym = dna_convert[p[i]];

		context_levels_t ctx_lev_codes;
		letters_ctx->determine_ctx_letters(ctx_lev_codes, i, ctx_letters, 0);

		context_t found_ctx;
		uint64_t v_counter;
		auto p_rc5 = find_rc_letters_context(m_ctx_rc_letters, ctx_lev_codes, tpl_ctx_rc_letters, found_ctx, v_counter);
		p_rc5->Encode(sym);

		update_ctx_letters(ctx_letters, sym);

		if (sym == 4)
		{
			sym = 0;	// N symbol
			cor_pos = i;
		}

		pmer_can.insert(sym);
		smer_can.insert(sym);
		bmer_can.insert(sym);

		pmer_can_unc.insert(sym);
		smer_can_unc.insert(sym);
		bmer_can_unc.insert(sym);
	}
}

//*****************************************************************************************************
void CDNACompressor::compress_prefix_sorted(const uint8_t *p, uint32_t size)
{
	uint8_t sym;
	reset_ctx_letters(ctx_letters);
	bool was_N = false;

	for (uint32_t i = 0; i < params.pmer_len; ++i)
	{
		sym = dna_convert[p[i]];
		update_ctx_letters(ctx_letters, sym);

		if (sym == 4)
		{
			sym = 3;	// N -> T
			was_N = true;
			N_run_len++;
		}
		else
			N_run_len = 0;

		pmer_can.insert(sym);
		smer_can.insert(sym);
		bmer_can.insert(sym);

		pmer_can_unc.insert(sym);
		smer_can_unc.insert(sym);
		bmer_can_unc.insert(sym);
	}

	uint64_t dif;
	uint64_t dif_no_bytes;
	uint64_t flag;

	// If prefix contains Ns we need to store their locations
	context_t ctx_Ns = 0;
	auto p_rc_Ns = find_rc_context(m_ctx_rc_prefix_Ns, ctx_Ns, tpl_ctx_rc_prefix_Ns);
	p_rc_Ns->Encode(was_N);

	ctx_prefix_sorted_flags = ((ctx_prefix_sorted_flags << 1) + (uint64_t) was_N) & 0xffff;

	if (pmer_can == pmer_can_prev)
		flag = max_prefix_sorted_flag_value;
	else
		flag = siv_pmer->test(pmer_can.data_aligned_dir());

	auto p_rc = find_rc_context(m_ctx_rc_prefix_sorted_flags, ctx_prefix_sorted_flags, tpl_ctx_rc_prefix_sorted_flags);
	p_rc->Encode((int) flag);
	ctx_prefix_sorted_flags = ((ctx_prefix_sorted_flags << 3) + flag) & 0xffff;

	if (flag < max_prefix_sorted_flag_value)
	{
		dif = 0;
		uint64_t max_i = pmer_can.data_aligned_dir();

		for (uint64_t i = pmer_can_prev.data_aligned_dir() + 1; i < max_i; ++i)
			if (siv_pmer->test(i) == flag)
				++dif;
		
		dif_no_bytes = no_bytes(dif);

		context_t ctx_prefix_sorted_no_bytes = ctx_prefix_sorted_flags;

		auto p_rc = find_rc_context(m_ctx_rc_prefix_sorted_no_bytes, ctx_prefix_sorted_no_bytes, tpl_ctx_rc_prefix_sorted_no_bytes);
		p_rc->Encode((int) dif_no_bytes - 1);

		if (dif_no_bytes == 1)
		{
			uint32_t hi_nibble = (uint32_t)(dif >> 4);
			uint32_t lo_nibble = (uint32_t) (dif & 0xf);

			context_t ctx = (1ull << 24) + flag;
			auto p_rc = find_rc_context(m_ctx_rc_prefix_sorted_nibbles, ctx, tpl_ctx_rc_prefix_sorted_nibbles);
			p_rc->Encode(hi_nibble);

			ctx = (2ull << 24) + flag * 256 + hi_nibble;
			p_rc = find_rc_context(m_ctx_rc_prefix_sorted_nibbles, ctx, tpl_ctx_rc_prefix_sorted_nibbles);
			p_rc->Encode(lo_nibble);
		}
		else
		{
			context_t ctx = flag * 65536 + dif_no_bytes * 256 + dif_no_bytes;
			auto p_rc = find_rc_context(m_ctx_rc_prefix_sorted_bytes, ctx, tpl_ctx_rc_prefix_sorted_bytes);
			uint64_t hi_byte = dif >> (dif_no_bytes * 8 - 8);
			p_rc->Encode((int) hi_byte);

			for (int i = 0; i < (int)dif_no_bytes - 1; ++i)
			{
				context_t ctx = (flag << 24) + (dif_no_bytes << 16) + (hi_byte << 8) + i;
				auto p_rc = find_rc_context(m_ctx_rc_prefix_sorted_bytes, ctx, tpl_ctx_rc_prefix_sorted_bytes);
				p_rc->Encode(dif & 0xff);
				dif >>= 8;
			}
		}
	}	

	// Encode flags to distinguish between Ts and Ns
	if (was_N)
	{
		for (uint32_t i = 0; i < params.pmer_len; ++i)
			if (p[i] == 'T' || p[i] == 'N')
			{
				auto p_rc_Ns = find_rc_context(m_ctx_rc_prefix_Ns, (uint64_t) i + 1, tpl_ctx_rc_prefix_Ns);
				p_rc_Ns->Encode(p[i] == 'N');
			}
	}
	
	pmer_can_prev = pmer_can;

	auto x = pmer_can.data_aligned_dir();
	(*my_pmers_to_add)[modulo_divisor(x >> pmer_mod_shift, no_threads)].push_back(x);
	x = pmer_can.data_aligned_rc();
	(*my_pmers_to_add)[modulo_divisor(x >> pmer_mod_shift, no_threads)].push_back(x);
}

//*****************************************************************************************************
void CDNACompressor::update_ctx_r_sym(uint64_t &ctx_r_sym, const int r_sym)
{
	ctx_r_sym <<= 1;
	if (r_sym == 0)
		ctx_r_sym += 1ull;

	ctx_r_sym &= code_ctx->ctx_r_sym_mask;
}

//*****************************************************************************************************
void CDNACompressor::compress_suffix(uint8_t *p, uint32_t size, bool original_order, uint32_t start_pos, bool reversed_pe)
{
	stats_t counts{0, 0, 0, 0};	
	context_t ctx_r_sym = 0;
	uint32_t prev_count = 0;
	int no_prev_letters = 0;

	for (uint32_t i = start_pos ? start_pos : original_order ? params.prefix_len : params.pmer_len; i < size; ++i)
	{
		uint8_t sym = dna_convert[p[i]];
		auto sym_to_kmers = (sym == 4) ? 0 : sym;

		// Insert empty symbol (technically 'A', but it does not matter here)
		pmer_can.insert_zero();
		smer_can.insert_zero();
		bmer_can.insert_zero();

		pmer_can_unc.insert_zero();
		smer_can_unc.insert_zero();
		bmer_can_unc.insert_zero();

		counts_level_t counts_level = find_counts(counts);

		if (counts_level == counts_level_t::bmer_unc)
		{
			bmer_can = bmer_can_unc;
			smer_can = smer_can_unc;
			pmer_can = pmer_can_unc;

			cor_pos = 0;
			counts_level = counts_level_t::bmer;
		}

		bool rough_counts = false;

		if (counts_level == counts_level_t::none)
		{
			if (bmer_can.is_full())
			{
				if (find_counts_rough_b(counts))
				{
					counts_level = counts_level_t::pmer;
					rough_counts = true;
				}
			}
			else if (smer_can.is_full())
			{
				if (find_counts_rough_s(counts))
				{
					counts_level = counts_level_t::pmer;
					rough_counts = true;
				}
			}
			else if(pmer_can.is_full() && counts_level == counts_level_t::none)
			{
				if (find_counts_rough_p(counts))
				{
					counts_level = counts_level_t::pmer;
					rough_counts = true;
				}
			}
		}
		
		if (counts_level != counts_level_t::none && N_run_len < 2)
		{
			int cor_dist = (counts_level == counts_level_t::pmer) ? params.pmer_len : (counts_level == counts_level_t::smer) ? params.smer_len : params.bmer_len;

			uint32_t cor_zone = ((int) i - (int) cor_pos < cor_dist) ? (1 + 2 * ((int) cor_dist - ((int) i - (int) cor_pos)) / (int)cor_dist) : 0;

			if (rough_counts)
				cor_zone = 3;

			context_levels_t ctx_lev_codes;
			if(!reversed_pe)
				code_ctx->determine_ctx_codes(ctx_lev_codes, counts, s_letters,
				i, counts_level, cor_zone, ctx_r_sym, size);
			else
				code_ctx->determine_ctx_codes(ctx_lev_codes, counts, s_letters,
					size - i - 1, counts_level, cor_zone, ctx_r_sym, ~0u);

			// Encode using counts
			context_t found_ctx;
			uint64_t v_counter;
			auto p_rc5 = find_rc_code_context(m_ctx_rc_codes, ctx_lev_codes, tpl_ctx_rc_codes, found_ctx, v_counter, &counts);

			uint8_t r_sym = rank(counts, sym);
			p_rc5->Encode(r_sym);

			if(bmer_can.is_almost_full(1))
			{
				CKmer bmer_tmp = bmer_can;
				bmer_tmp.replace_last(sym);
				bmer_tmp.insert_zero();
				ht_bmer->prefetch_data_kernel(bmer_tmp.kernel_canonical());
				bmer_tmp.insert_zero();
				ht_bmer->prefetch_idx_kernel(bmer_tmp.kernel_canonical());
			}

			no_prev_letters = 0;

			update_ctx_r_sym(ctx_r_sym, r_sym);
		}
		else
		{
			// Encode plain
			context_levels_t ctx_lev_codes;
			letters_ctx->determine_ctx_letters(ctx_lev_codes, i, ctx_letters, size);

			context_t found_ctx;
			uint64_t v_counter;
			auto p_rc5 = find_rc_letters_context(m_ctx_rc_letters, ctx_lev_codes, tpl_ctx_rc_letters, found_ctx, v_counter);
			p_rc5->Encode(sym);

			if (bmer_can.is_almost_full(1))
			{
				CKmer bmer_tmp = bmer_can;
				bmer_tmp.replace_last(sym);
				bmer_tmp.insert_zero();
				ht_bmer->prefetch_data_kernel(bmer_tmp.kernel_canonical());
				bmer_tmp.insert_zero();
				ht_bmer->prefetch_idx_kernel(bmer_tmp.kernel_canonical());
			}

			prev_count = 0;
			no_prev_letters++;

			update_ctx_r_sym(ctx_r_sym);
		}

		update_ctx_letters(ctx_letters, sym);

		if (sym == 4)
			++N_run_len;
		else
			N_run_len = 0;

		pmer_can.replace_last(sym_to_kmers);
		smer_can.replace_last(sym_to_kmers);
		bmer_can.replace_last(sym_to_kmers);
		
		pmer_can_unc.replace_last(sym_to_kmers);
		smer_can_unc.replace_last(sym_to_kmers);
		bmer_can_unc.replace_last(sym_to_kmers);

		if (sym < 4)
		{
			bool pmer_insert = true;

			if (bmer_can.is_full())
			{
				auto x = bmer_can.data_normalized();
				(*my_bmers_to_add)[modulo_divisor((x >> bsmer_mod_shift) & bsmer_mod_mask, no_threads)].push_back(x);
				ht_bmer_local->insert(x, &cinc_lb);

				if ((counts_level == counts_level_t::smer || counts_level == counts_level_t::bmer || counts_level == counts_level_t::mixed || counts_level == counts_level_t::bmer_unc) 
					&& counts[sym] >= (1u << SIV_FIELD_SIZE) - 1)
					pmer_insert = false;
			}

			if (smer_can.is_full())
			{
				auto x = smer_can.data_normalized();
				(*my_smers_to_add)[modulo_divisor((x >> bsmer_mod_shift) & bsmer_mod_mask, no_threads)].push_back(x);
				ht_smer_local->insert(x, &cinc_ls);
			}

			if (pmer_can.is_full() && i - cor_pos >= params.pmer_len-1)
			{
				if (pmer_insert)
				{
					auto x = pmer_can.data_aligned_dir();
					(*my_pmers_to_add)[modulo_divisor(x >> pmer_mod_shift, no_threads)].push_back(x);
					x = pmer_can.data_aligned_rc();
					(*my_pmers_to_add)[modulo_divisor(x >> pmer_mod_shift, no_threads)].push_back(x);
				}
				else
					no_pmer_hidden_updates += 2;
			}
		}

		if (bmer_can.is_full())
		{
			if (counts_level == counts_level_t::bmer || counts_level == counts_level_t::mixed || counts_level == counts_level_t::bmer_unc)
			{
				if (repair_kmers_existing(p + i, i, counts, sym, prev_count))
				{
					auto x = bmer_can.data_normalized();
					(*my_bmers_to_add)[modulo_divisor((x >> bsmer_mod_shift) & bsmer_mod_mask, no_threads)].push_back(x);
					ht_bmer_local->insert(x, &cinc_lb);
				}
			}
			else if ((counts_level == counts_level_t::none || counts_level == counts_level_t::pmer)
				&& bmer_can.is_full())
			{
				if (repair_kmers_missing(p + i, i, sym))
				{
					auto x = bmer_can.data_normalized();
					(*my_bmers_to_add)[modulo_divisor((x >> bsmer_mod_shift) & bsmer_mod_mask, no_threads)].push_back(x);
					ht_bmer_local->insert(x, &cinc_lb);
				}
			}
		}
	}
}

//*****************************************************************************************************
bool CDNACompressor::is_valid_minimizer(const uint64_t x) const
{
	uint64_t first_3_symbols = x >> (2 * params.bmer_len - 6);

	if (first_3_symbols == 0)		// AAA...
		return false;
	if (first_3_symbols == 1)		// AAC...
		return false;

	return true;
}

//*****************************************************************************************************
bool CDNACompressor::is_valid_maximizer(const uint64_t x) const
{
	uint64_t first_3_symbols = x >> (2 * params.bmer_len - 6);

	if (first_3_symbols == 0b111110)		// TTG...
		return false;
	if (first_3_symbols == 0b111111)		// TTT...
		return false;

	return true;
}

//*****************************************************************************************************
void CDNACompressor::merge_minim_results(vector<uint64_t> &v_cand, vector<uint64_t> &v_top_cand)
{
	if (v_cand.size() == 1)
	{
		swap(v_cand, v_top_cand);
		return;
	}

	uint64_t counter_shift = 2 * (uint64_t) params.bmer_len;
	uint64_t mask_value = (1ull << counter_shift) - 1ull;
	uint64_t max_counter = (~0ull) >> counter_shift;

	if (v_minim_cand.size() > 3 * (uint64_t) no_examined_pe_minim)
	{
		partial_sort(v_minim_cand.begin(), v_minim_cand.begin() + 3 * (uint64_t) no_examined_pe_minim, v_minim_cand.end(), [&](uint64_t x, uint64_t y) {
			uint64_t x_counter = x >> counter_shift;
			uint64_t y_counter = y >> counter_shift;

			if (x_counter != y_counter)
				return x_counter > y_counter;
			return (x & mask_value) < (y & mask_value);
		});
		v_minim_cand.resize(3 * (uint64_t)no_examined_pe_minim);
	}

	sort(v_minim_cand.begin(), v_minim_cand.end(), [mask_value](uint64_t x, uint64_t y) {
		return (x & mask_value) < (y & mask_value);
	});

	v_top_cand.clear();
	v_top_cand.reserve(v_cand.size());
	v_top_cand.push_back(v_cand.front());

	for (size_t i = 1; i < v_cand.size(); ++i)
	{
		if ((v_top_cand.back() & mask_value) != (v_cand[i] & mask_value))
			v_top_cand.push_back(v_cand[i]);
		else
		{
			uint64_t cnt_x = v_top_cand.back() >> counter_shift;
			uint64_t cnt_y = v_cand[i] >> counter_shift;
			if (cnt_x + cnt_y > max_counter)
				cnt_y = max_counter - cnt_x;
			v_top_cand.back() += cnt_y << counter_shift;
		}
	}

	if(v_top_cand.size() <= no_examined_pe_minim)
		sort(v_top_cand.begin(), v_top_cand.end(), [&](uint64_t x, uint64_t y) {
			uint64_t x_counter = x >> counter_shift;
			uint64_t y_counter = y >> counter_shift;

			if (x_counter != y_counter)
				return x_counter > y_counter;
			return (x & mask_value) < (y & mask_value);
		});
	else
		partial_sort(v_top_cand.begin(), v_top_cand.begin() + no_examined_pe_minim, v_top_cand.end(), [&](uint64_t x, uint64_t y) {
			uint64_t x_counter = x >> counter_shift;
			uint64_t y_counter = y >> counter_shift;

			if (x_counter != y_counter)
				return x_counter > y_counter;
			return (x & mask_value) < (y & mask_value);
		});
}

//*****************************************************************************************************
void CDNACompressor::generate_read_bmers(const uint8_t *p, uint32_t size)
{
	CKmer bmer;
	uint32_t i;
	uint8_t sym;

	v_read_bmers.clear();
	v_read_bmers.reserve(size);

	bmer.Reset(params.bmer_len, kmer_mode_t::direct);

	for (i = 0; i < size; ++i)
	{
		sym = dna_convert[p[i]];

		if (sym == 4)
			bmer.Reset();
		else
		{
			bmer.insert(sym);
			auto cand = bmer.data_aligned_dir();
			if (bmer.is_full() && is_valid_minimizer(cand))
				v_read_bmers.emplace_back(make_pair(cand, i - (params.bmer_len - 1)));
		}
	}
}

//*****************************************************************************************************
uint64_t CDNACompressor::find_minimizer(const uint8_t *p, uint32_t size) const
{
	CKmer bmer;
	uint64_t bmer_min_al = ht_pe_mers->get_value_mask();
	uint32_t i;
	uint8_t sym;

	bmer.Reset(params.bmer_len, kmer_mode_t::direct);

	for (i = 0; i < size; ++i)
	{
		sym = dna_convert[p[i]];

		if (sym == 4)
			bmer.Reset();
		else
		{
			bmer.insert(sym);
			auto cand = bmer.data_aligned_dir();
			if (bmer.is_full() && cand < bmer_min_al && is_valid_minimizer(cand))
				bmer_min_al = cand;
		}
	}

	return bmer_min_al;
}

//*****************************************************************************************************
uint64_t CDNACompressor::find_maximizer(const uint8_t *p, uint32_t size) const
{
	CKmer bmer;
	uint64_t bmer_max_al = 0ull;
	int32_t i;
	uint8_t sym;

	bmer.Reset(params.bmer_len, kmer_mode_t::direct);

	for (i = size-1; i >= 0; --i)
	{
		sym = dna_convert[p[i]];

		if (sym == 4)
			bmer.Reset();
		else
		{
			bmer.insert(sym);
			auto cand = bmer.data_aligned_dir();
			if (bmer.is_full() && cand > bmer_max_al && is_valid_maximizer(cand))
				bmer_max_al = cand;
		}
	}

	return bmer_max_al;
}

//*****************************************************************************************************
void CDNACompressor::append_pe_mers3(const uint8_t *p1, uint32_t size1, const uint8_t *p2, uint32_t size2)
{
	int minim_start_size = size1 - params.bmer_len + 1;
	int start_point_1 = minim_start_size / 3;
	int start_point_2 = 2 * minim_start_size / 3;

	auto minim_1_1st = find_minimizer(p1, start_point_1 + params.bmer_len - 1);
	auto minim_1_2nd = find_minimizer(p1 + start_point_1, start_point_2 - start_point_1 + params.bmer_len - 1);
	auto minim_1_3rd = find_minimizer(p1 + start_point_2, size1 - start_point_2);

	minim_start_size = size2 - params.bmer_len + 1;
	start_point_1 = minim_start_size / 3;
	start_point_2 = 2 * minim_start_size / 3;

	auto minim_2_1st = find_minimizer(p2, start_point_1 + params.bmer_len - 1);
	auto minim_2_2nd = find_minimizer(p2 + start_point_1, start_point_2 - start_point_1 + params.bmer_len - 1);
	auto minim_2_3rd = find_minimizer(p2 + start_point_2, size2 - start_point_2);

	auto minim_1_1st_dest = ht_pe_mers->get_part_id(minim_1_1st);
	auto minim_1_2nd_dest = ht_pe_mers->get_part_id(minim_1_2nd);
	auto minim_1_3rd_dest = ht_pe_mers->get_part_id(minim_1_3rd);
	auto minim_2_1st_dest = ht_pe_mers->get_part_id(minim_2_1st);
	auto minim_2_2nd_dest = ht_pe_mers->get_part_id(minim_2_2nd);
	auto minim_2_3rd_dest = ht_pe_mers->get_part_id(minim_2_3rd);

	int mid_point1 = (size1 + params.bmer_len) / 2;
	int mid_point2 = (size2 + params.bmer_len) / 2;

	auto maxim_1_second = find_maximizer(p1 + mid_point1 - params.bmer_len + 1, size1 - (mid_point1 - params.bmer_len + 1));
	auto maxim_2_second = find_minimizer(p2 + mid_point2 - params.bmer_len + 1, size2 - (mid_point2 - params.bmer_len + 1));
	maxim_1_second = (~maxim_1_second) & ht_pe_mers->get_value_mask();
	maxim_2_second = (~maxim_2_second) & ht_pe_mers->get_value_mask();

	(*my_pe_mers_to_add)[minim_1_1st_dest].emplace_back(make_tuple(minim_1_1st, minim_2_1st, 2));
//	(*my_pe_mers_to_add)[minim_1_1st_dest].emplace_back(make_tuple(minim_1_1st, minim_2_2nd, 3));
	(*my_pe_mers_to_add)[minim_1_1st_dest].emplace_back(make_tuple(minim_1_1st, minim_2_3rd, 4));
	(*my_pe_mers_to_add)[minim_1_1st_dest].emplace_back(make_tuple(minim_1_1st, maxim_1_second, 1));

	(*my_pe_mers_to_add)[minim_1_2nd_dest].emplace_back(make_tuple(minim_1_2nd, minim_2_1st, 3));
//	(*my_pe_mers_to_add)[minim_1_2nd_dest].emplace_back(make_tuple(minim_1_2nd, minim_2_2nd, 4));
	(*my_pe_mers_to_add)[minim_1_2nd_dest].emplace_back(make_tuple(minim_1_2nd, minim_2_3rd, 3));

	(*my_pe_mers_to_add)[minim_1_3rd_dest].emplace_back(make_tuple(minim_1_3rd, minim_2_1st, 4));
//	(*my_pe_mers_to_add)[minim_1_3rd_dest].emplace_back(make_tuple(minim_1_3rd, minim_2_2nd, 3));
	(*my_pe_mers_to_add)[minim_1_3rd_dest].emplace_back(make_tuple(minim_1_3rd, minim_2_3rd, 2));

	(*my_pe_mers_to_add)[minim_2_1st_dest].emplace_back(make_tuple(minim_2_1st, minim_1_1st, 2));
//	(*my_pe_mers_to_add)[minim_2_1st_dest].emplace_back(make_tuple(minim_2_1st, minim_1_2nd, 3));
	(*my_pe_mers_to_add)[minim_2_1st_dest].emplace_back(make_tuple(minim_2_1st, minim_1_3rd, 4));
	(*my_pe_mers_to_add)[minim_2_1st_dest].emplace_back(make_tuple(minim_2_1st, maxim_2_second, 1));

	(*my_pe_mers_to_add)[minim_2_2nd_dest].emplace_back(make_tuple(minim_2_2nd, minim_1_1st, 3));
//	(*my_pe_mers_to_add)[minim_2_2nd_dest].emplace_back(make_tuple(minim_2_2nd, minim_1_2nd, 4));
	(*my_pe_mers_to_add)[minim_2_2nd_dest].emplace_back(make_tuple(minim_2_2nd, minim_1_3rd, 4));

	(*my_pe_mers_to_add)[minim_2_3rd_dest].emplace_back(make_tuple(minim_2_3rd, minim_1_1st, 4));
//	(*my_pe_mers_to_add)[minim_2_3rd_dest].emplace_back(make_tuple(minim_2_3rd, minim_1_2nd, 3));
	(*my_pe_mers_to_add)[minim_2_3rd_dest].emplace_back(make_tuple(minim_2_3rd, minim_1_3rd, 2));

	ht_pe_mers_local->insert(minim_1_1st, minim_2_1st, 2);
	ht_pe_mers_local->insert(minim_1_1st, minim_2_3rd, 4);
	ht_pe_mers_local->insert(minim_1_1st, maxim_1_second, 1);

	ht_pe_mers_local->insert(minim_1_2nd, minim_2_1st, 3);
	ht_pe_mers_local->insert(minim_1_2nd, minim_2_3rd, 3);

	ht_pe_mers_local->insert(minim_1_3rd, minim_2_1st, 4);
	ht_pe_mers_local->insert(minim_1_3rd, minim_2_3rd, 2);

	ht_pe_mers_local->insert(minim_2_1st, minim_1_1st, 2);
	ht_pe_mers_local->insert(minim_2_1st, minim_1_3rd, 4);
	ht_pe_mers_local->insert(minim_2_1st, maxim_2_second, 1);

	ht_pe_mers_local->insert(minim_2_2nd, minim_1_1st, 3);
	ht_pe_mers_local->insert(minim_2_2nd, minim_1_3rd, 4);

	ht_pe_mers_local->insert(minim_2_3rd, minim_1_1st, 4);
	ht_pe_mers_local->insert(minim_2_3rd, minim_1_3rd, 2);
}

//*****************************************************************************************************
void CDNACompressor::decompress_suffix(uint8_t *p, uint32_t size, bool original_order, uint32_t start_pos, bool reversed_pe)
{
	stats_t counts{ 0,0,0,0 };
	uint8_t sym;
	context_t ctx_r_sym = 0;
	uint32_t prev_count = 0;
	int no_prev_letters = 0;

	for (uint32_t i = start_pos ? start_pos : original_order ? params.prefix_len : params.pmer_len; i < size; ++i)
	{
		// Insert empty symbol (technically 'A', but it does not matter here)
		pmer_can.insert_zero();
		smer_can.insert_zero();
		bmer_can.insert_zero();

		pmer_can_unc.insert_zero();
		smer_can_unc.insert_zero();
		bmer_can_unc.insert_zero();

		counts_level_t counts_level = find_counts(counts);

		if (counts_level == counts_level_t::bmer_unc)
		{
			bmer_can = bmer_can_unc;
			smer_can = smer_can_unc;
			pmer_can = pmer_can_unc;

			cor_pos = 0;
			counts_level = counts_level_t::bmer;
		}

		bool rough_counts = false;

		if (counts_level == counts_level_t::none)
		{
			if (bmer_can.is_full())
			{
				if (find_counts_rough_b(counts))
				{
					counts_level = counts_level_t::pmer;
					rough_counts = true;
				}
			}
			else if (smer_can.is_full())
			{
				if (find_counts_rough_s(counts))
				{
					counts_level = counts_level_t::pmer;
					rough_counts = true;
				}
			}
			else if (pmer_can.is_full())
			{
				if (find_counts_rough_p(counts))
				{
					counts_level = counts_level_t::pmer;
					rough_counts = true;
				}
			}
		}

		if (counts_level != counts_level_t::none && N_run_len < 2)
		{
			int cor_dist = (counts_level == counts_level_t::pmer) ? params.pmer_len : (counts_level == counts_level_t::smer) ? params.smer_len : params.bmer_len;

			uint32_t cor_zone = ((int)i - (int)cor_pos < cor_dist) ? (1 + 2 * ((int)cor_dist - ((int)i - (int)cor_pos)) / (int)cor_dist) : 0;

			if (rough_counts)
				cor_zone = 3;

			context_levels_t ctx_lev_codes;
			if (!reversed_pe)
				code_ctx->determine_ctx_codes(ctx_lev_codes, counts, s_letters,
					i, counts_level, cor_zone, ctx_r_sym, size);
			else
				code_ctx->determine_ctx_codes(ctx_lev_codes, counts, s_letters,
					size - i - 1, counts_level, cor_zone, ctx_r_sym, ~0u);

			// Encode using counts
			context_t found_ctx;
			uint64_t v_counter;
			auto p_rc5 = find_rc_code_context(m_ctx_rc_codes, ctx_lev_codes, tpl_ctx_rc_codes, found_ctx, v_counter, &counts);

			// Decode using counts
			uint8_t r_sym = p_rc5->Decode();
			sym = un_rank(counts, r_sym);

			if (bmer_can.is_almost_full(1))
			{
				CKmer bmer_tmp = bmer_can;
				bmer_tmp.replace_last(sym);
				bmer_tmp.insert_zero();
				ht_bmer->prefetch_data_kernel(bmer_tmp.kernel_canonical());
				bmer_tmp.insert_zero();
				ht_bmer->prefetch_idx_kernel(bmer_tmp.kernel_canonical());
			}

			no_prev_letters = 0;

			update_ctx_r_sym(ctx_r_sym, r_sym);
		}
		else
		{
			// Decode plain
			context_levels_t ctx_lev_codes;
			letters_ctx->determine_ctx_letters(ctx_lev_codes, i, ctx_letters, size);

			context_t found_ctx;
			uint64_t v_counter;
			auto p_rc5 = find_rc_letters_context(m_ctx_rc_letters, ctx_lev_codes, tpl_ctx_rc_letters, found_ctx, v_counter);
			sym = p_rc5->Decode();

			if (bmer_can.is_almost_full(1))
			{
				CKmer bmer_tmp = bmer_can;
				bmer_tmp.replace_last(sym);
				bmer_tmp.insert_zero();
				ht_bmer->prefetch_data_kernel(bmer_tmp.kernel_canonical());
				bmer_tmp.insert_zero();
				ht_bmer->prefetch_idx_kernel(bmer_tmp.kernel_canonical());
			}

			prev_count = 0;
			no_prev_letters++;

			update_ctx_r_sym(ctx_r_sym);
		}

		p[i] = dna_alphabet[sym];

		auto sym_to_kmers = (sym == 4) ? 0 : sym;

		update_ctx_letters(ctx_letters, sym);

		if (sym == 4)
			++N_run_len;
		else
			N_run_len = 0;

		pmer_can.replace_last(sym_to_kmers);
		smer_can.replace_last(sym_to_kmers);
		bmer_can.replace_last(sym_to_kmers);

		pmer_can_unc.replace_last(sym_to_kmers);
		smer_can_unc.replace_last(sym_to_kmers);
		bmer_can_unc.replace_last(sym_to_kmers);

		if (sym < 4)
		{
			bool pmer_insert = true;

			if (bmer_can.is_full())
			{
				auto x = bmer_can.data_normalized();
				(*my_bmers_to_add)[modulo_divisor((x >> bsmer_mod_shift) & bsmer_mod_mask, no_threads)].push_back(x);
				ht_bmer_local->insert(x, &cinc_lb);

				if ((counts_level == counts_level_t::smer || counts_level == counts_level_t::bmer || counts_level == counts_level_t::mixed || counts_level == counts_level_t::bmer_unc)
					&& counts[sym] >= (1u << SIV_FIELD_SIZE) - 1)
					pmer_insert = false;
			}

			if (smer_can.is_full())
			{
				auto x = smer_can.data_normalized();
				(*my_smers_to_add)[modulo_divisor((x >> bsmer_mod_shift) & bsmer_mod_mask, no_threads)].push_back(x);
				ht_smer_local->insert(x, &cinc_ls);
			}

			if (pmer_can.is_full() && i - cor_pos >= params.pmer_len - 1)
			{
				if (pmer_insert)
				{
					auto x = pmer_can.data_aligned_dir();
					(*my_pmers_to_add)[modulo_divisor(x >> pmer_mod_shift, no_threads)].push_back(x);
					x = pmer_can.data_aligned_rc();
					(*my_pmers_to_add)[modulo_divisor(x >> pmer_mod_shift, no_threads)].push_back(x);
				}
				else
					no_pmer_hidden_updates += 2;
			}
		}

		if (bmer_can.is_full())
		{
			if (counts_level == counts_level_t::bmer || counts_level == counts_level_t::mixed || counts_level == counts_level_t::bmer_unc)
			{
				if (repair_kmers_existing(p + i, i, counts, sym, prev_count))
				{
					auto x = bmer_can.data_normalized();
					(*my_bmers_to_add)[modulo_divisor((x >> bsmer_mod_shift) & bsmer_mod_mask, no_threads)].push_back(x);
					ht_bmer_local->insert(x, &cinc_lb);
				}
			}
			else if ((counts_level == counts_level_t::none || counts_level == counts_level_t::pmer)
				&& bmer_can.is_full())
			{
				if (repair_kmers_missing(p + i, i, sym))
				{
					auto x = bmer_can.data_normalized();
					(*my_bmers_to_add)[modulo_divisor((x >> bsmer_mod_shift) & bsmer_mod_mask, no_threads)].push_back(x);
					ht_bmer_local->insert(x, &cinc_lb);
				}
			}
		}
	}
}

//*****************************************************************************************************
void CDNACompressor::decompress_prefix_direct(uint8_t *p, uint32_t size)
{
	uint8_t sym;
	reset_ctx_letters(ctx_letters);

	for (uint32_t i = 0; i < params.prefix_len; ++i)
	{
		context_levels_t ctx_lev_codes;
		letters_ctx->determine_ctx_letters(ctx_lev_codes, i, ctx_letters, 0);

		context_t found_ctx;
		uint64_t v_counter;
		auto p_rc5 = find_rc_letters_context(m_ctx_rc_letters, ctx_lev_codes, tpl_ctx_rc_letters, found_ctx, v_counter);
		sym = p_rc5->Decode();

		p[i] = dna_alphabet[sym];

		update_ctx_letters(ctx_letters, sym);

		if (sym == 4)
		{
			sym = 0;
			cor_pos = i;
		}

		pmer_can.insert(sym);
		smer_can.insert(sym);
		bmer_can.insert(sym);

		pmer_can_unc.insert(sym);
		smer_can_unc.insert(sym);
		bmer_can_unc.insert(sym);
	}
}

//*****************************************************************************************************
void CDNACompressor::decompress_prefix_sorted(uint8_t *p, uint32_t size)
{
	pmer_can.Reset();
	smer_can.Reset();
	bmer_can.Reset();

	pmer_can_unc.Reset();
	smer_can_unc.Reset();
	bmer_can_unc.Reset();

	// We need to decode Ns in prefix (if any exist)
	context_t ctx_Ns = 0;
	auto p_rc_Ns = find_rc_context(m_ctx_rc_prefix_Ns, ctx_Ns, tpl_ctx_rc_prefix_Ns);
	bool was_N = p_rc_Ns->Decode();

	ctx_prefix_sorted_flags = ((ctx_prefix_sorted_flags << 1) + (uint64_t)was_N) & 0xffff;

	auto p_rc = find_rc_context(m_ctx_rc_prefix_sorted_flags, ctx_prefix_sorted_flags, tpl_ctx_rc_prefix_sorted_flags);
	uint64_t flag = p_rc->Decode();

	ctx_prefix_sorted_flags = ((ctx_prefix_sorted_flags << 3) + flag) & 0xffff;

	uint64_t dif;

	if (flag < max_prefix_sorted_flag_value)
	{
		context_t ctx_prefix_sorted_no_bytes = ctx_prefix_sorted_flags;
		auto p_rc = find_rc_context(m_ctx_rc_prefix_sorted_no_bytes, ctx_prefix_sorted_no_bytes, tpl_ctx_rc_prefix_sorted_no_bytes);
		int dif_no_bytes = p_rc->Decode() + 1;

		if (dif_no_bytes == 1)
		{
			context_t ctx = (1ull << 24) + flag;
			auto p_rc = find_rc_context(m_ctx_rc_prefix_sorted_nibbles, ctx, tpl_ctx_rc_prefix_sorted_nibbles);
			uint32_t hi_nibble = p_rc->Decode();

			ctx = (2ull << 24) + flag * 256 + hi_nibble;
			p_rc = find_rc_context(m_ctx_rc_prefix_sorted_nibbles, ctx, tpl_ctx_rc_prefix_sorted_nibbles);
			uint32_t lo_nibble = p_rc->Decode();

			dif = ((uint64_t) hi_nibble << 4) + (uint64_t) lo_nibble;
		}
		else
		{
			context_t ctx = (uint64_t) flag * 65536 + (uint64_t)dif_no_bytes * 256 + (uint64_t)dif_no_bytes;
			auto p_rc = find_rc_context(m_ctx_rc_prefix_sorted_bytes, ctx, tpl_ctx_rc_prefix_sorted_bytes);

			uint64_t hi_byte = (uint64_t) p_rc->Decode();
			dif = hi_byte << (dif_no_bytes * 8 - 8);

			for (int i = 0; i < (int)dif_no_bytes - 1; ++i)
			{
				context_t ctx = ((uint64_t) flag << 24) + ((uint64_t)dif_no_bytes << 16) + ((uint64_t)hi_byte << 8) + (uint64_t) i;
				auto p_rc = find_rc_context(m_ctx_rc_prefix_sorted_bytes, ctx, tpl_ctx_rc_prefix_sorted_bytes);
				uint64_t x = p_rc->Decode();
				dif += x << (i * 8);
			}
		}

		uint64_t k = pmer_can_prev.data_aligned_dir() + 1;

		for (uint64_t j = 0; j <= dif; ++k)
			if (siv_pmer->test(k) == flag)
				++j;
		--k;

		for (uint32_t i = 0; i < params.pmer_len; ++i)
		{
			pmer_can.insert_front(k & 3);
			k >>= 2;
		}
	}
	else
	{
		if(pmer_can_prev.is_full())
			pmer_can = pmer_can_prev;
		else
		{
			// It can happen than pmer_can_prev is empty and its binary representation is 0, which can be equal to some valid k-mer
			for (uint32_t i = 0; i < params.pmer_len; ++i)
				pmer_can.insert(0);
			pmer_can_prev = pmer_can;
		}
	}

	uint8_t sym;
	reset_ctx_letters(ctx_letters);

	for (uint32_t i = 0; i < params.pmer_len; ++i)
	{
		sym = (uint8_t) pmer_can.get_symbol(i);
		p[i] = dna_alphabet[sym];

		if (was_N)
		{
			if(sym == 3)
			{
				auto p_rc_Ns = find_rc_context(m_ctx_rc_prefix_Ns, (uint64_t) i + 1, tpl_ctx_rc_prefix_Ns);
				if (p_rc_Ns->Decode())
				{
					p[i] = 'N';
					sym = 4;
				}
			}
		}

		update_ctx_letters(ctx_letters, sym);

		if (sym == 4)
		{
			sym = 3;	// N -> T
			N_run_len++;
		}
		else
			N_run_len = 0;

		smer_can.insert(sym);
		bmer_can.insert(sym);

		pmer_can_unc.insert(sym);
		smer_can_unc.insert(sym);
		bmer_can_unc.insert(sym);
	}

	pmer_can_prev = pmer_can;

	auto x = pmer_can.data_aligned_dir();
	(*my_pmers_to_add)[modulo_divisor(x >> pmer_mod_shift, no_threads)].push_back(x);
	x = pmer_can.data_aligned_rc();
	(*my_pmers_to_add)[modulo_divisor(x >> pmer_mod_shift, no_threads)].push_back(x);
}

//*****************************************************************************************************
bool CDNACompressor::CompressDirect(uint8_t *p, uint32_t size, uint8_t *q, bool first_read_of_pair)
{
	reset_context(ctx_letters);

	read_cur.assign(p, p + size);

	if (params.duplicates_check && first_read_of_pair)
	{
		bool same_read = read_cur == read_prev;

		auto p_rc = find_rc_context(m_ctx_rc_flags, ctx_flags, tpl_ctx_rc_flags);
		p_rc->Encode(same_read);
		ctx_flags = ((ctx_flags << 1) + same_read) & 0xff;

		if (same_read)
			return true;
	}

	// Kodowanie reada
	pmer_can.Reset();
	smer_can.Reset();
	bmer_can.Reset();

	pmer_can_unc.Reset();
	smer_can_unc.Reset();
	bmer_can_unc.Reset();

	cor_pos = 0;
	N_run_len = 0;

	compress_prefix_direct(p, size, 0);
	compress_suffix(p, size, true, 0);

	if(first_read_of_pair)
		read_prev.swap(read_cur);

	update_s_letters(p, size);

	return true;
}

//*****************************************************************************************************
bool CDNACompressor::CompressDirectWithMinim(uint8_t *p, uint32_t size, uint32_t minim_pos)
{
	reset_context(ctx_letters);

	read_cur.assign(p, p + size);

	pmer_can.Reset();
	smer_can.Reset();
	bmer_can.Reset();

	pmer_can_unc.Reset();
	smer_can_unc.Reset();
	bmer_can_unc.Reset();

	cor_pos = 0;
	N_run_len = 0;

	reset_ctx_letters(ctx_letters);

	for (int i = minim_pos; i < (int) (minim_pos + params.bmer_len); ++i)
	{
		uint8_t c = dna_convert[p[i]];

		update_ctx_letters(ctx_letters, c);

		if (c == 4)
			c = 0;

		pmer_can.insert(c);
		bmer_can.insert(c);
		smer_can.insert(c);
		pmer_can_unc.insert(c);
		bmer_can_unc.insert(c);
		smer_can_unc.insert(c);
	}

	compress_suffix(p, size, true, params.bmer_len + minim_pos);
		
	uint8_t *rc_p = (uint8_t*) alloca((uint64_t)minim_pos + params.bmer_len);
	int rc_size = 0;

	for(int i = (int) minim_pos + (int) params.bmer_len - 1; i >= 0; --i)
		rc_p[rc_size++] = reverse_complement_alhpa(p[i]);

	pmer_can.Reset();
	smer_can.Reset();
	bmer_can.Reset();

	pmer_can_unc.Reset();
	smer_can_unc.Reset();
	bmer_can_unc.Reset();

	cor_pos = 0;
	N_run_len = 0;

	reset_ctx_letters(ctx_letters);

	for (int i = 0; i < (int) params.bmer_len; ++i)
	{
		uint8_t c = dna_convert[rc_p[i]];

		update_ctx_letters(ctx_letters, c);

		if (c == 4)
			c = 0;

		pmer_can.insert(c);
		bmer_can.insert(c);
		smer_can.insert(c);
		pmer_can_unc.insert(c);
		bmer_can_unc.insert(c);
		smer_can_unc.insert(c);
	}

	compress_suffix(rc_p, rc_size, true, params.bmer_len, true);
	
	update_s_letters(p, size);

	return true;
}

//*****************************************************************************************************
bool CDNACompressor::DecompressDirectWithMinim(uint8_t *p, uint32_t size, uint32_t minim_pos, uint64_t minim)
{
	int32_t i;

	reset_context(ctx_letters);

	pmer_can.Reset();
	smer_can.Reset();
	bmer_can.Reset();
	pmer_can_unc.Reset();
	smer_can_unc.Reset();
	bmer_can_unc.Reset();

	for (i = 0; i < (int) params.bmer_len; ++i)
	{
		uint64_t sym = (minim >> 62);

		update_ctx_letters(ctx_letters, sym);

		minim <<= 2;
		pmer_can.insert(sym);
		smer_can.insert(sym);
		bmer_can.insert(sym);
		pmer_can_unc.insert(sym);
		smer_can_unc.insert(sym);
		bmer_can_unc.insert(sym);
		p[minim_pos + i] = dna_alphabet[sym];
	}

	decompress_suffix(p, size, true, params.bmer_len + minim_pos);

	uint8_t *rc_p = (uint8_t*)alloca(minim_pos + params.bmer_len);
	int rc_size = 0;

	for (i = (int) minim_pos + (int) params.bmer_len - 1; i >= (int) minim_pos; --i)
		rc_p[rc_size++] = reverse_complement_alhpa(p[i]);

	pmer_can.Reset();
	smer_can.Reset();
	bmer_can.Reset();
	pmer_can_unc.Reset();
	smer_can_unc.Reset();
	bmer_can_unc.Reset();

	cor_pos = 0;
	N_run_len = 0;

	reset_ctx_letters(ctx_letters);

	for (i = 0; i < rc_size; ++i)
	{
		uint8_t c = dna_convert[rc_p[i]];

		update_ctx_letters(ctx_letters, c);

		if (c == 4)
			c = 0;

		pmer_can.insert(c);
		smer_can.insert(c);
		bmer_can.insert(c);
		pmer_can_unc.insert(c);
		smer_can_unc.insert(c);
		bmer_can_unc.insert(c);
	}

	decompress_suffix(rc_p, params.bmer_len + minim_pos, true, params.bmer_len, true);

	for (i = (int) minim_pos - 1; i >= 0; --i)
		p[i] = reverse_complement_alhpa(rc_p[rc_size++]);

	return true;
}

//*****************************************************************************************************
bool CDNACompressor::CompressSorted(uint8_t *p, uint32_t size, bool first_read_of_pair)
{
	reset_context(ctx_letters);

	read_cur.assign(p, p + size);

	if (params.duplicates_check && first_read_of_pair)
	{
		bool same_read = read_cur == read_prev;

		auto p_rc = find_rc_context(m_ctx_rc_flags, ctx_flags, tpl_ctx_rc_flags);
		p_rc->Encode(same_read);
		ctx_flags = ((ctx_flags << 1) + same_read) & 0xff;

		if (same_read)
			return true;
	}
	
	pmer_can.Reset();
	smer_can.Reset();
	bmer_can.Reset();

	pmer_can_unc.Reset();
	smer_can_unc.Reset();
	bmer_can_unc.Reset();

	cor_pos = 0;
	N_run_len = 0;

	compress_prefix_sorted(p, size);
	compress_suffix(p, size, false, 0);

	if(first_read_of_pair)
		read_prev.swap(read_cur);

	update_s_letters(p, size);

	return true;
}

//*****************************************************************************************************
bool CDNACompressor::find_minim_cand(const uint8_t *p, uint32_t size) 
{
	v_minim_cand.clear();

	int minim_start_size = size - params.bmer_len + 1;
	int start_point_1 = minim_start_size / 4;
	int start_point_2 = 2 * minim_start_size / 4;
	int start_point_3 = 3 * minim_start_size / 4;

	auto minim_1st = find_minimizer(p, start_point_1 + params.bmer_len - 1);
	auto minim_2nd = find_minimizer(p + start_point_1, start_point_2 - start_point_1 + params.bmer_len - 1);
	auto minim_3rd = find_minimizer(p + start_point_2, start_point_3 - start_point_2 + params.bmer_len - 1);
	auto minim_4th = find_minimizer(p + start_point_3, size - start_point_3);

	ht_pe_mers->find(minim_1st, v_minim_cand);
	ht_pe_mers->find(minim_2nd, v_minim_cand);
	ht_pe_mers->find(minim_3rd, v_minim_cand);
	ht_pe_mers->find(minim_4th, v_minim_cand);

	ht_pe_mers_local->find(minim_1st, v_minim_cand);
	ht_pe_mers_local->find(minim_2nd, v_minim_cand);
	ht_pe_mers_local->find(minim_3rd, v_minim_cand);
	ht_pe_mers_local->find(minim_4th, v_minim_cand);

	if (v_minim_cand.empty())
		return false;

	merge_minim_results(v_minim_cand, v_minim_top_cand);

	return true;
}

//*****************************************************************************************************
bool CDNACompressor::CompressPE(uint8_t *p1, uint32_t size1, uint8_t *p2, uint32_t size2)
{
	// 1st read of a pair
	if (params.dna_mode == dna_mode_t::pe_original)
		CompressDirect(p1, size1, nullptr, true);
	else
		CompressSorted(p1, size1, true);
	
	bool minim_found = find_minim_cand(p1, size1);

	uint32_t minim2_pos;
	int minim2_id = -1;
	
	if(!minim_found)
		a_stat_minim_id[20]++;
	else
	{
		uint64_t minim_mask = ht_pe_mers->get_value_mask();
		
		generate_read_bmers(p2, size2);

		for (size_t i = 0; i < v_minim_top_cand.size(); ++i)
		{
			auto m = v_minim_top_cand[i] & minim_mask;
			auto p = find_if(v_read_bmers.begin(), v_read_bmers.end(), [m](pair<uint64_t, uint32_t> x) {
				return x.first == m;
			});

			if (p != v_read_bmers.end())
			{
				minim2_id = (int) i;
				minim2_pos = p->second;
				break;
			}
		}

		if (minim2_id < 0)
			a_stat_minim_id[21]++;
		else if (minim2_id < 15)
			a_stat_minim_id[minim2_id]++;
		else
			a_stat_minim_id[15]++;

		if (minim2_id < 0 || minim2_id > 14)
			minim2_id = 15;
	}

	if (minim2_id < 0)
		CompressDirect(p2, size2, nullptr, false);
	else
	{
		ctx_rc_pe_minimizer_id->Encode(minim2_id);
		if (minim2_id == 15)
			CompressDirect(p2, size2, nullptr, false);
		else
		{
			// Encode minimizer pos
			if (minim2_pos < 254)
			{
				auto p_rc = find_rc_context(m_ctx_rc_minimizer_pos, minim2_id, tpl_ctx_rc_minimizer_pos);
				p_rc->Encode(minim2_pos);
			}
			else if(minim2_pos < 256 * 256)
			{
				auto p_rc = find_rc_context(m_ctx_rc_minimizer_pos, minim2_id, tpl_ctx_rc_minimizer_pos);
				p_rc->Encode(254);
				p_rc = find_rc_context(m_ctx_rc_minimizer_pos, (uint64_t) minim2_id + 0x100, tpl_ctx_rc_minimizer_pos);
				p_rc->Encode((int) (minim2_pos >> 8));
				p_rc = find_rc_context(m_ctx_rc_minimizer_pos, (uint64_t) minim2_id + 0x200, tpl_ctx_rc_minimizer_pos);
				p_rc->Encode((int) (minim2_pos & 0xff));
			}
			else
			{
				auto p_rc = find_rc_context(m_ctx_rc_minimizer_pos, (uint64_t) minim2_id, tpl_ctx_rc_minimizer_pos);
				p_rc->Encode(255);
				p_rc = find_rc_context(m_ctx_rc_minimizer_pos, (uint64_t)minim2_id + 0x300, tpl_ctx_rc_minimizer_pos);
				p_rc->Encode((int) (minim2_pos >> 16));
				p_rc = find_rc_context(m_ctx_rc_minimizer_pos, (uint64_t)minim2_id + 0x400, tpl_ctx_rc_minimizer_pos);
				p_rc->Encode((int) ((minim2_pos >> 8) & 0xff));
				p_rc = find_rc_context(m_ctx_rc_minimizer_pos, (uint64_t)minim2_id + 0x500, tpl_ctx_rc_minimizer_pos);
				p_rc->Encode((int) (minim2_pos & 0xff));
			}

			CompressDirectWithMinim(p2, size2, minim2_pos);
		}
	}

	append_pe_mers3(p1, size1, p2, size2);

	return true;
}

//*****************************************************************************************************
bool CDNACompressor::DecompressSE(uint8_t *p, uint32_t size)
{
	reset_context(ctx_letters);

	// Flag decoding
	if (params.duplicates_check)
	{
		auto p_rc = find_rc_context(m_ctx_rc_flags, ctx_flags, tpl_ctx_rc_flags);
		bool same_read = p_rc->Decode();
		ctx_flags = ((ctx_flags << 1) + same_read) & 0xff;

		if (same_read)
		{
			copy_n(read_prev.begin(), size, p);
			return true;
		}
	}

	pmer_can.Reset();
	smer_can.Reset();
	bmer_can.Reset();

	pmer_can_unc.Reset();
	smer_can_unc.Reset();
	bmer_can_unc.Reset();

	cor_pos = 0;
	N_run_len = 0;

	if (params.dna_mode == dna_mode_t::se_original)
	{
		decompress_prefix_direct(p, size);
		decompress_suffix(p, size, true, 0);
	}
	else
	{
		decompress_prefix_sorted(p, size);
		decompress_suffix(p, size, false, 0);
	}

	read_prev.assign(p, p + size);

	update_s_letters(p, size);

	return true;
}

//*****************************************************************************************************
bool CDNACompressor::DecompressPE(uint8_t *p1, uint32_t size1, uint8_t *p2, uint32_t size2)
{
	// 1st read of a pair
	reset_context(ctx_letters);
	bool same_read = false;

	// Flag decoding
	if (params.duplicates_check)
	{
		auto p_rc = find_rc_context(m_ctx_rc_flags, ctx_flags, tpl_ctx_rc_flags);
		same_read = p_rc->Decode();
		ctx_flags = ((ctx_flags << 1) + same_read) & 0xff;

		if (same_read)
			copy_n(read_prev.begin(), size1, p1);
	}

	if (!same_read)
	{
		pmer_can.Reset();
		smer_can.Reset();
		bmer_can.Reset();

		pmer_can_unc.Reset();
		smer_can_unc.Reset();
		bmer_can_unc.Reset();

		cor_pos = 0;
		N_run_len = 0;

		if (params.dna_mode == dna_mode_t::pe_original)
		{
			decompress_prefix_direct(p1, size1);
			decompress_suffix(p1, size1, true, 0);
		}
		else
		{
			decompress_prefix_sorted(p1, size1);
			decompress_suffix(p1, size1, false, 0);
		}

		read_prev.assign(p1, p1 + size1);

		update_s_letters(p1, size1);
	}

	// 2nd read of a pair
	reset_context(ctx_letters);

	bool minim_found = find_minim_cand(p1, size1);
	bool decompress_direct = false;
	uint64_t minim2_id;
	uint32_t minim2_pos;

	if (!minim_found)
	{
		decompress_direct = true;
	}
	else
	{
		minim2_id = ctx_rc_pe_minimizer_id->Decode();
		if (minim2_id == 15)
			decompress_direct = true;
		else
		{
			auto p_rc = find_rc_context(m_ctx_rc_minimizer_pos, minim2_id, tpl_ctx_rc_minimizer_pos);
			minim2_pos = p_rc->Decode();

			if(minim2_pos == 254)		// 2-byte value
			{
				minim2_pos = 0;
				p_rc = find_rc_context(m_ctx_rc_minimizer_pos, minim2_id + 0x100, tpl_ctx_rc_minimizer_pos);
				minim2_pos += p_rc->Decode() << 8;
				p_rc = find_rc_context(m_ctx_rc_minimizer_pos, minim2_id + 0x200, tpl_ctx_rc_minimizer_pos);
				minim2_pos += p_rc->Decode();
			}
			else if(minim2_pos == 255)	// 3-byte value
			{
				minim2_pos = 0;
				p_rc = find_rc_context(m_ctx_rc_minimizer_pos, minim2_id + 0x300, tpl_ctx_rc_minimizer_pos);
				minim2_pos += p_rc->Decode() << 16;
				p_rc = find_rc_context(m_ctx_rc_minimizer_pos, minim2_id + 0x400, tpl_ctx_rc_minimizer_pos);
				minim2_pos += p_rc->Decode() << 8;
				p_rc = find_rc_context(m_ctx_rc_minimizer_pos, minim2_id + 0x500, tpl_ctx_rc_minimizer_pos);
				minim2_pos += p_rc->Decode();
			}
		}
	}

	pmer_can.Reset();
	smer_can.Reset();
	bmer_can.Reset();

	pmer_can_unc.Reset();
	smer_can_unc.Reset();
	bmer_can_unc.Reset();

	cor_pos = 0;
	N_run_len = 0;

	if (decompress_direct)
	{
		decompress_prefix_direct(p2, size2);
		decompress_suffix(p2, size2, true, 0);
	}
	else
		DecompressDirectWithMinim(p2, size2, minim2_pos, v_minim_top_cand[minim2_id] << (64 - 2 * params.bmer_len));

	update_s_letters(p2, size2);

	append_pe_mers3(p1, size1, p2, size2);

	return true;
}

//*****************************************************************************************************
void CDNACompressor::update_s_letters(const uint8_t *p, uint32_t size)
{
	for (uint32_t i = 0; i < size; ++i)
	{
		if (dna_convert[p[i]] == 4)
			continue;

		s_letters[dna_convert[p[i]]]++;
		s_letters[reverse_complement(dna_convert[p[i]])]++;
	}
}

//*****************************************************************************************************
CDNACompressor::ctx_map_t::value_type CDNACompressor::find_rc_context(ctx_map_t &m_ctx_rc, context_t ctx, CRangeCoderModel<CVectorIOStream> *tpl, const stats_t *counts)
{
	auto p = m_ctx_rc.find(ctx);

	if (p == nullptr)
	{
		m_ctx_rc.insert(ctx, p = new CRangeCoderModel<CVectorIOStream>(*tpl));

		if (counts)
		{
			int init[5];
			partial_sort_copy(counts->begin(), counts->end(), init, init + 4, greater<uint32_t>());
			init[4] = 1;
			transform(init, init + 4, init, [](auto x) {return 4 * x + 1; });

			p->Init(init);

			++n_rc_ctx;
		}
	}

	return p;
}

//*****************************************************************************************************
CDNACompressor::ctx_map_2_t::value_type CDNACompressor::find_rc_context(ctx_map_2_t &m_ctx_rc, context_t ctx, CRangeCoderModelFixedSize<CVectorIOStream, 2> *tpl, const stats_t *counts)
{
	auto p = m_ctx_rc.find(ctx);

	if (p == nullptr)
		m_ctx_rc.insert(ctx, p = new CRangeCoderModelFixedSize<CVectorIOStream, 2>(*tpl));

	return p;
}

//*****************************************************************************************************
CDNACompressor::ctx_map_5_t::value_type CDNACompressor::find_rc_context(ctx_map_5_t &m_ctx_rc, context_t ctx, CRangeCoderModelFixedSize<CVectorIOStream, 5> *tpl, const stats_t *counts)
{
	auto p = m_ctx_rc.find(ctx);

	if (p == nullptr)
		m_ctx_rc.insert(ctx, p = new CRangeCoderModelFixedSize<CVectorIOStream, 5>(*tpl));

	return p;
}

//*****************************************************************************************************
CDNACompressor::ctx_map_5_t::value_type CDNACompressor::find_rc_code_context(ctx_map_5_t &m_ctx_rc, context_levels_t &ctx_lev, 
	CRangeCoderModelFixedSize<CVectorIOStream, 5> *tpl, 
	context_t &found_ctx,
	uint64_t &v_counter,
	const stats_t *counts)
{
	int i;
	ctx_map_5_t::value_type p, q;
	size_t *p_counter, *q_counter;

	p = nullptr;

	int start_lev = (int)(avg_code_ctx_level + 0.49);
	p = m_ctx_rc.find_ext(ctx_lev[start_lev], p_counter);
	if (p != nullptr && *p_counter < new_code_ctx_thr[start_lev])
	{
		*p_counter += 1;
		v_counter = *p_counter;
		update_avg(avg_code_ctx_level, avg_mult, start_lev);

		found_ctx = ctx_lev[start_lev];

		return p;
	}

	if (p == nullptr)
	{
		for (i = start_lev-1; i >= 0; --i)
		{
			p = m_ctx_rc.find_ext(ctx_lev[i], p_counter);
			if (p != nullptr)
			{
				found_ctx = ctx_lev[i];
				break;
			}
		}
	}
	else
	{
		for (i = start_lev + 1; i < (int) code_ctx->no_ctx_levels; ++i)
		{
			q = m_ctx_rc.find_ext(ctx_lev[i], q_counter);
			if (q == nullptr)
			{
				found_ctx = ctx_lev[(uint64_t)i-1];
				break;
			}
			if (*q_counter < new_code_ctx_thr[i])
			{
				update_avg(avg_code_ctx_level, avg_mult, i);
				*q_counter += 1;
				found_ctx = ctx_lev[i];
				v_counter = *q_counter;

				return q;
			}
			p = q;
			p_counter = q_counter;
		}
		--i;
	}

	if (p == nullptr)
	{
		m_ctx_rc.insert(ctx_lev[0], new CRangeCoderModelFixedSize<CVectorIOStream, 5>(*tpl));
		p = m_ctx_rc.find_ext(ctx_lev[0], p_counter);
		*p_counter += 1;
		v_counter = *p_counter;

		i = 0;
		found_ctx = ctx_lev[0];
	}

	if (*p_counter >= new_code_ctx_thr[i] && i + 1 < (int) code_ctx->no_ctx_levels)
	{
		m_ctx_rc.insert(ctx_lev[(uint64_t)i + 1], new CRangeCoderModelFixedSize<CVectorIOStream, 5>(*p));
		p = m_ctx_rc.find_ext(ctx_lev[(uint64_t)i + 1], p_counter);
		found_ctx = ctx_lev[(uint64_t)i+1];
		*p_counter += 1;
		v_counter = *p_counter;
	}
	else
	{
		*p_counter += 1;
		v_counter = *p_counter;
	}

	update_avg(avg_code_ctx_level, avg_mult, i);

	return p;
}

//*****************************************************************************************************
CDNACompressor::ctx_map_5_t::value_type CDNACompressor::find_rc_letters_context(ctx_map_5_t &m_ctx_rc, context_levels_t &ctx_lev, CRangeCoderModelFixedSize<CVectorIOStream, 5> *tpl,
	context_t &found_ctx, uint64_t &v_counter)
{
	int i;
	ctx_map_5_t::value_type p, q;
	size_t *p_counter, *q_counter;

	p = nullptr;

	int start_lev = (int)(avg_letters_ctx_level + 0.49);
	p = m_ctx_rc.find_ext(ctx_lev[start_lev], p_counter);
	if (p != nullptr && *p_counter < new_letters_ctx_thr[start_lev])
	{
		*p_counter += 1;
		v_counter = *p_counter;
		update_avg(avg_letters_ctx_level, avg_mult, start_lev);

		found_ctx = ctx_lev[start_lev];

		return p;
	}

	if (p == nullptr)
	{
		for (i = start_lev - 1; i >= 0; --i)
		{
			p = m_ctx_rc.find_ext(ctx_lev[i], p_counter);
			if (p != nullptr)
			{
				found_ctx = ctx_lev[i];
				break;
			}
		}
	}
	else
	{
		for (i = start_lev + 1; i < (int)letters_ctx->no_ctx_levels; ++i)
		{
			q = m_ctx_rc.find_ext(ctx_lev[i], q_counter);
			if (q == nullptr)
			{
				found_ctx = ctx_lev[(uint64_t)i - 1];
				break;
			}
			if (*q_counter < new_code_ctx_thr[i])
			{
				update_avg(avg_letters_ctx_level, avg_mult, i);
				*q_counter += 1;
				found_ctx = ctx_lev[i];
				v_counter = *q_counter;

				return q;
			}
			p = q;
			p_counter = q_counter;
		}
		--i;
	}

	if (p == nullptr)
	{
		m_ctx_rc.insert(ctx_lev[0], new CRangeCoderModelFixedSize<CVectorIOStream, 5>(*tpl));
		p = m_ctx_rc.find_ext(ctx_lev[0], p_counter);
		*p_counter += 1;
		v_counter = *p_counter;

		i = 0;
		found_ctx = ctx_lev[0];
	}

	if (*p_counter >= new_code_ctx_thr[i] && i + 1 < (int)letters_ctx->no_ctx_levels)
	{
		m_ctx_rc.insert(ctx_lev[(uint64_t)i + 1], new CRangeCoderModelFixedSize<CVectorIOStream, 5>(*p));
		p = m_ctx_rc.find_ext(ctx_lev[(uint64_t)i + 1], p_counter);
		found_ctx = ctx_lev[(uint64_t)i + 1];
		*p_counter += 1;
		v_counter = *p_counter;
	}
	else
	{
		*p_counter += 1;
		v_counter = *p_counter;
	}
	update_avg(avg_letters_ctx_level, avg_mult, i);

	return p;
}

//*****************************************************************************************************
CDNACompressor::ctx_map_16_t::value_type CDNACompressor::find_rc_context(ctx_map_16_t &m_ctx_rc, context_t ctx, CRangeCoderModelFixedSize<CVectorIOStream, 16> *tpl, const stats_t *counts)
{
	auto p = m_ctx_rc.find(ctx);

	if (p == nullptr)
		m_ctx_rc.insert(ctx, p = new CRangeCoderModelFixedSize<CVectorIOStream, 16>(*tpl));

	return p;
}

//*****************************************************************************************************
CDNACompressor::ctx_map_256_t::value_type CDNACompressor::find_rc_context(ctx_map_256_t &m_ctx_rc, context_t ctx, CRangeCoderModelFixedSize<CVectorIOStream, 256> *tpl, const stats_t *counts)
{
	auto p = m_ctx_rc.find(ctx);

	if (p == nullptr)
		m_ctx_rc.insert(ctx, p = new CRangeCoderModelFixedSize<CVectorIOStream, 256>(*tpl));

	return p;
}

//*****************************************************************************************************
void CDNACompressor::create_rc_context(ctx_map_t &m_ctx_rc, context_t new_ctx, ctx_map_t::value_type base_p_rc)
{
	auto p = m_ctx_rc.find(new_ctx);

	if (p == nullptr)
		m_ctx_rc.insert(new_ctx, new CRangeCoderModel<CVectorIOStream>(*base_p_rc));
}

//*****************************************************************************************************
void CDNACompressor::create_rc_context(ctx_map_2_t &m_ctx_rc, context_t new_ctx, ctx_map_2_t::value_type base_p_rc)
{
	auto p = m_ctx_rc.find(new_ctx);

	if (p == nullptr)
		m_ctx_rc.insert(new_ctx, new CRangeCoderModelFixedSize<CVectorIOStream, 2>(*base_p_rc));
}

//*****************************************************************************************************
void CDNACompressor::create_rc_context(ctx_map_5_t &m_ctx_rc, context_t new_ctx, ctx_map_5_t::value_type base_p_rc)
{
	auto p = m_ctx_rc.find(new_ctx);

	if (p == nullptr)
		m_ctx_rc.insert(new_ctx, new CRangeCoderModelFixedSize<CVectorIOStream, 5>(*base_p_rc));
}

//*****************************************************************************************************
void CDNACompressor::create_rc_context(ctx_map_16_t &m_ctx_rc, context_t new_ctx, ctx_map_16_t::value_type base_p_rc)
{
	auto p = m_ctx_rc.find(new_ctx);

	if (p == nullptr)
		m_ctx_rc.insert(new_ctx, new CRangeCoderModelFixedSize<CVectorIOStream, 16>(*base_p_rc));
}

//*****************************************************************************************************
void CDNACompressor::create_rc_context(ctx_map_256_t &m_ctx_rc, context_t new_ctx, ctx_map_256_t::value_type base_p_rc)
{
	auto p = m_ctx_rc.find(new_ctx);

	if (p == nullptr)
		m_ctx_rc.insert(new_ctx, new CRangeCoderModelFixedSize<CVectorIOStream, 256>(*base_p_rc));
}

//*****************************************************************************************************
void CDNACompressor::SetKmerDS(TSmallIntVector<SIV_FIELD_SIZE> *_siv_pmer, 
	CHT_kmer<uint32_t> *_ht_smer, CHT_kmer<uint32_t> *_ht_bmer, CHT_pair_kmers *_ht_pe_mers, 
	vector<vector<vector<uint64_t>>> *_pmers_to_add, vector<vector<vector<uint64_t>>> *_smers_to_add, 
	vector<vector<vector<uint64_t>>> *_bmers_to_add, vector<vector<vector<tuple<uint64_t, uint64_t, uint64_t>>>> *_pe_mers_to_add,
	uint64_t _thread_id)
{
	siv_pmer = _siv_pmer;
	ht_smer = _ht_smer;
	ht_bmer = _ht_bmer;
	ht_pe_mers = _ht_pe_mers;

	pmers_to_add = _pmers_to_add;
	smers_to_add = _smers_to_add;
	bmers_to_add = _bmers_to_add;
	pe_mers_to_add = _pe_mers_to_add;

	thread_id = _thread_id;

	my_pmers_to_add = &(*pmers_to_add)[thread_id];
	my_smers_to_add = &(*smers_to_add)[thread_id];
	my_bmers_to_add = &(*bmers_to_add)[thread_id];
	my_pe_mers_to_add = &(*pe_mers_to_add)[thread_id];

	no_threads = (uint32_t) pmers_to_add->size();

	pmer_mod_shift = 2 * params.pmer_len - 12;
	bsmer_mod_shift = 46;
	no_pmer_hidden_updates = 0;

#ifdef KMER_MARGIN_2_SYMBOLS
	bsmer_mod_mask = (1ull << (64 - bsmer_mod_shift - 4)) - 1ull;
#else
	bsmer_mod_mask = (1ull << (64 - bsmer_mod_shift - 2)) - 1ull;
#endif
}

//*****************************************************************************************************
void CDNACompressor::InsertKmersToHT()
{
	// p-mers
	uint64_t siv_no_filled = 0;
	uint64_t siv_no_updates = 0;

	uint64_t siv_prefetch_step = 32;

	for (uint64_t i = 0; i < no_threads; ++i)
	{
		auto &cur_vec = (*pmers_to_add)[i][thread_id];
		uint64_t v_size = cur_vec.size();

		for (uint64_t j = 0; j < v_size; ++j)
		{
			if (j + siv_prefetch_step < v_size)
				siv_pmer->prefetch(cur_vec[j + siv_prefetch_step]);

			siv_no_filled += (uint64_t)siv_pmer->increment(cur_vec[j]);
			++siv_no_updates;
		}
	}

	siv_pmer->update_no_filled(siv_no_filled);
	siv_pmer->update_no_updates(siv_no_updates + no_pmer_hidden_updates);
	no_pmer_hidden_updates = 0;

	// s-mers and b-mers
	uint64_t step_idx = 16 * 2;
	uint64_t step_data = 8 * 2;
	uint64_t step_pe = 16 * 2;

	for (auto kmer_data : { make_tuple(ht_smer, smers_to_add, &cinc_s), make_tuple(ht_bmer, bmers_to_add, &cinc_b) })
	{
		uint64_t my_inserts = 0;

		for (uint64_t i = 0; i < no_threads; ++i)
		{
			uint64_t max_j = (*get<1>(kmer_data))[i][thread_id].size();
			auto &vec = (*get<1>(kmer_data))[i][thread_id];
			
			for (uint64_t j = 0; j < max_j; ++j)
			{
				if (j + step_idx < max_j)
					get<0>(kmer_data)->prefetch_idx(vec[j + step_idx]);

				if (j + step_data < max_j)
					get<0>(kmer_data)->prefetch_data(vec[j + step_data]);

				get<0>(kmer_data)->insert(vec[j], get<2>(kmer_data));
				my_inserts++;
			}
		}
	}

	for (uint64_t i = 0; i < no_threads; ++i)
	{
		uint64_t no_pe_mers = (*pe_mers_to_add)[i][thread_id].size();
		auto &pe_mers_to_do = (*pe_mers_to_add)[i][thread_id];
		
		uint64_t prev_prefetch_key = 0ull;

		for (uint64_t j = 0; j < no_pe_mers; ++j)
		{
			if (j + step_pe < no_pe_mers)
			{
				uint64_t prefetch_key = get<0>(pe_mers_to_do[j + step_pe]);

				if (prefetch_key != prev_prefetch_key)
				{
					prev_prefetch_key = prefetch_key;
					ht_pe_mers->prefetch(thread_id, prefetch_key);
				}
			}

			auto &x = pe_mers_to_do[j];
			ht_pe_mers->insert(get<0>(x), get<1>(x), get<2>(x));
		}
	}
}

//*****************************************************************************************************
void CDNACompressor::ClearKmersToHT()
{
	for (int i = 0; i < (int) no_threads; ++i)
	{
		(*my_pmers_to_add)[i].clear();
		(*my_smers_to_add)[i].clear();
		(*my_bmers_to_add)[i].clear();
		(*my_pe_mers_to_add)[i].clear();
	}

	ht_bmer_local->clear(0);
	ht_smer_local->clear(0);
	ht_pe_mers_local->clear();
}

//*****************************************************************************************************
bool CDNACompressor::NeedRestructHT()
{
	return ht_pe_mers->need_restruct_soon();
}

//*****************************************************************************************************
void CDNACompressor::DoRestruct()
{
	ht_pe_mers->RestructIfNecessary(thread_id);
}

//*****************************************************************************************************
void CDNACompressor::show_contexts()
{
	vector<CContextHM<CRangeCoderModelFixedSize<CVectorIOStream, 5u>>::item_t> v_ctx;

	m_ctx_rc_codes.debug_list(v_ctx);

	FILE *ctx_log = fopen("ctx.log", "wb");
	fprintf(ctx_log, "Symbols - No contexts in map: %lld\n", (long long) v_ctx.size());
	vector<float> v_freq;
	size_t no_updates;

	double est_len = 0;
	double est_tot_len = 0;
	double est_lev_len[6] = { 0,0,0,0,0,0 };
	int lev;

	for (auto x : v_ctx)
	{
		x.rcm->GetSimpleModel()->GetLogStats(no_updates, v_freq);
		est_len = x.rcm->GetEstTotLen();
		est_tot_len += est_len;

		no_updates = x.counter;
		fprintf(ctx_log, "%16.3lf   %16llx  - %s  - %llu    0:%7.3f     1:%7.3f     2:%7.3f     3:%7.3f     4:%7.3f\n",
			est_len / 8,
			(unsigned long long) x.ctx, 
			code_ctx->decode_ctx(x.ctx, lev).c_str(),
			(long long unsigned int) no_updates, 100.0*v_freq[0], 100.0*v_freq[1], 100.0*v_freq[2], 100.0*v_freq[3], 100.0*v_freq[4]);

		est_lev_len[lev] += est_len;
	}

	fprintf(ctx_log, "Tot. est. len: %lld\n", (long long) est_tot_len);
	for(int i = 0; i < 6; ++i)
		fprintf(ctx_log, "Lev %d est. len: %lld\n", i, (long long) est_lev_len[i]);

	fclose(ctx_log);

	est_len = 0;
	est_tot_len = 0;

	ctx_log = fopen("ctx_let.log", "wb");
	m_ctx_rc_letters.debug_list(v_ctx);

	uint64_t tot_no_updates = 0;

	for (auto x : v_ctx)
	{
		x.rcm->GetSimpleModel()->GetLogStats(no_updates, v_freq);
		est_len = x.rcm->GetEstTotLen();
		est_tot_len += est_len;
		fprintf(ctx_log, "%16.3lf   %16llx  - %s   - %lld   0:%7.3f     1:%7.3f     2:%7.3f     3:%7.3f     4:%7.3f\n",
			(est_len / 8), (unsigned long long) x.ctx, letters_ctx->decode_ctx(x.ctx).c_str(), (long long) no_updates, 100.0*v_freq[0], 100.0*v_freq[1], 100.0*v_freq[2], 100.0*v_freq[3], 100.0*v_freq[4]);
		tot_no_updates += no_updates;
	}

	fprintf(ctx_log, "Total no. of updates: %lld\n", (long long) tot_no_updates);
	fprintf(ctx_log, "Total est_len : %lld\n", (long long) est_tot_len);

	fclose(ctx_log);
}

//*****************************************************************************************************
std::string CDNACompressor::decode_ctx_let(context_t ctx)
{
	string s = "[";

	s += "pos:" + to_string((ctx >> 50) & 0x3ffull) + "  ";
	for (int i = ctx_letters_len-1; i >= 0; --i)
		s.push_back("ACGTN              ."[(ctx >> (4 * i)) & 0xf]);

	s += "]";

	return s;
}

// EOF
