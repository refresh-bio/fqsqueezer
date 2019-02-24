#pragma once
// *******************************************************************************************
// This file is a part of FQSqueezer software distributed under GNU GPL 3 licence.
// The homepage of the MSAC project is http://sun.aei.polsl.pl/REFRESH/fqsqueezer
//
// Author: Sebastian Deorowicz
// Version: 1.0
// Date   : 2019-02-22
// *******************************************************************************************

#include <random>

#include "defs.h"
#include "compressor.h"
#include "params.h"
#include "io.h"
#include "rc.h"
#include "sub_rc.h"
#include "kmer.h"
#include "ht_kmer.h"
#include "bit_vec.h"
#include "code_ctx.h"

#include <array>
#include <map>
#include <unordered_map>

class CDNACompressor : public CBasicCompressor
{
	const uint32_t no_ctx_symbols = 1;
	typedef uint64_t context_t;

	const size_t new_code_ctx_thr[12] = { 1, 32, 64, 64, 128, 512, 1024, 32, 128, 256, 2048, 4 };					// a22
	double avg_code_ctx_level = 0.0;

	const size_t new_letters_ctx_thr[12] = { 1, 32, 64, 128, 256, 512, 2048, 4096, 8192, 16384, 16384, 4 };					// a22
	double avg_letters_ctx_level = 0.0;

	double avg_mult = 0.999;

	uint32_t no_threads;
	uint32_t pmer_mod_shift;
	uint32_t bsmer_mod_shift;
	uint32_t bsmer_mod_mask;
	uint64_t no_pmer_hidden_updates;
	double siv_prune_thr;

	ctx_map_2_t m_ctx_rc_flags;
	ctx_map_5_t m_ctx_rc_letters;
	ctx_map_5_t m_ctx_rc_codes;

	ctx_map_t m_ctx_rc_prefix_sorted_flags;
	ctx_map_t m_ctx_rc_prefix_sorted_no_bytes;
	ctx_map_16_t m_ctx_rc_prefix_sorted_nibbles;
	ctx_map_256_t m_ctx_rc_prefix_sorted_bytes;
	ctx_map_2_t m_ctx_rc_prefix_Ns;
	ctx_map_t m_ctx_rc_minimizer_pos;

	array<uint64_t, 23> a_stat_minim_id;

	CRangeCoderModelFixedSize<CVectorIOStream, 2> *tpl_ctx_rc_flags;
	CRangeCoderModelFixedSize<CVectorIOStream, 5> *tpl_ctx_rc_letters;
	CRangeCoderModelFixedSize<CVectorIOStream, 5> *tpl_ctx_rc_codes;

	CRangeCoderModel<CVectorIOStream> *tpl_ctx_rc_prefix_sorted_flags;
	CRangeCoderModel<CVectorIOStream> *tpl_ctx_rc_prefix_sorted_no_bytes;
	CRangeCoderModelFixedSize<CVectorIOStream, 16> *tpl_ctx_rc_prefix_sorted_nibbles;
	CRangeCoderModelFixedSize<CVectorIOStream, 256> *tpl_ctx_rc_prefix_sorted_bytes;
	CRangeCoderModelFixedSize<CVectorIOStream, 2> *tpl_ctx_rc_prefix_Ns;

	CRangeCoderModelFixedSize<CVectorIOStream, 16> *ctx_rc_pe_minimizer_id;
	CRangeCoderModel<CVectorIOStream> *tpl_ctx_rc_minimizer_pos;

	context_t ctx_flags;
	context_t ctx_letters;

	context_t ctx_prefix_sorted_flags;

	map<pair<uint64_t, uint64_t>, uint64_t> m_no_ctx_card;
	map<uint64_t, uint64_t> m_min_ctx_card;

	const uint32_t ctx_letters_len = 6;
	const context_t ctx_letters_mask = (1ull << (4*ctx_letters_len)) - 1ull;
	const uint32_t no_examined_pe_minim = 16;

	mt19937 mt;
	CCounterIncrementer cinc_b;
	CCounterIncrementer cinc_s;
	CCounterIncrementer cinc_lb;
	CCounterIncrementer cinc_ls;

	stats64_t s_letters;

	uint64_t n_rc_ctx;

	void reset_context(context_t &ctx) const {
		ctx = (context_t)0;
	}

	void update_avg(double &avg, double avg_mult, double cur_val) const {
		avg = avg_mult * avg + (1.0 - avg_mult) * cur_val;
	}

	void reset_ctx_letters(context_t &ctx) 	{
		ctx = ~(context_t)0;
	}

	void update_ctx_letters(context_t &ctx, context_t sym)	{
		ctx = (ctx << 4) + sym;
	}
	
	ctx_map_t::value_type find_rc_context(ctx_map_t &m_ctx_rc, context_t ctx, CRangeCoderModel<CVectorIOStream> *tpl, const stats_t *counts = nullptr);
	ctx_map_2_t::value_type find_rc_context(ctx_map_2_t &m_ctx_rc, context_t ctx, CRangeCoderModelFixedSize<CVectorIOStream, 2> *tpl, const stats_t *counts = nullptr);
	ctx_map_5_t::value_type find_rc_context(ctx_map_5_t &m_ctx_rc, context_t ctx, CRangeCoderModelFixedSize<CVectorIOStream, 5> *tpl, const stats_t *counts = nullptr);
	ctx_map_16_t::value_type find_rc_context(ctx_map_16_t &m_ctx_rc, context_t ctx, CRangeCoderModelFixedSize<CVectorIOStream, 16> *tpl, const stats_t *counts = nullptr);
	ctx_map_256_t::value_type find_rc_context(ctx_map_256_t &m_ctx_rc, context_t ctx, CRangeCoderModelFixedSize<CVectorIOStream, 256> *tpl, const stats_t *counts = nullptr);

	ctx_map_5_t::value_type find_rc_code_context(ctx_map_5_t &m_ctx_rc, context_levels_t &ctx_lev, CRangeCoderModelFixedSize<CVectorIOStream, 5> *tpl, 
		context_t &found_ctx, uint64_t &v_counter,
		const stats_t *counts = nullptr);

	ctx_map_5_t::value_type find_rc_letters_context(ctx_map_5_t &m_ctx_rc, context_levels_t &ctx_lev, CRangeCoderModelFixedSize<CVectorIOStream, 5> *tpl,
		context_t &found_ctx, uint64_t &v_counter);

	void create_rc_context(ctx_map_t &m_ctx_rc, context_t new_ctx, ctx_map_t::value_type base_p_rc);
	void create_rc_context(ctx_map_2_t &m_ctx_rc, context_t new_ctx, ctx_map_2_t::value_type base_p_rc);
	void create_rc_context(ctx_map_5_t &m_ctx_rc, context_t new_ctx, ctx_map_5_t::value_type base_p_rc);
	void create_rc_context(ctx_map_16_t &m_ctx_rc, context_t new_ctx, ctx_map_16_t::value_type base_p_rc);
	void create_rc_context(ctx_map_256_t &m_ctx_rc, context_t new_ctx, ctx_map_256_t::value_type base_p_rc);

	char dna_convert[256];

	vector<uint8_t> read_prev;
	vector<uint8_t> read_cur;
	vector<uint8_t> qual_cur;

	// Hash tables
	TSmallIntVector<SIV_FIELD_SIZE> *siv_pmer;
	CHT_kmer<uint32_t> *ht_smer;
	CHT_kmer<uint32_t> *ht_bmer;
	CHT_pair_kmers *ht_pe_mers;

	CHT_kmer<uint64_t> *ht_bmer_local;
	CHT_kmer<uint64_t> *ht_smer_local;
	CHT_pair_kmers *ht_pe_mers_local;

	CCodeContext *code_ctx;
	CLettersContext *letters_ctx;

	vector<vector<vector<uint64_t>>> *pmers_to_add;
	vector<vector<vector<uint64_t>>> *smers_to_add;
	vector<vector<vector<uint64_t>>> *bmers_to_add;
	vector<vector<vector<tuple<uint64_t, uint64_t, uint64_t>>>> *pe_mers_to_add;

	vector<vector<uint64_t>> *my_pmers_to_add;
	vector<vector<uint64_t>> *my_smers_to_add;
	vector<vector<uint64_t>> *my_bmers_to_add;
	vector<vector<tuple<uint64_t, uint64_t, uint64_t>>> *my_pe_mers_to_add;

	uint64_t thread_id;

	typedef struct {
		CKmer p_dir, p_rc;
		CKmer s_can;
		CKmer b_can;
	} kmer_pack_t;

	typedef struct {
		uint32_t counts_level;
		kmer_pack_t kmers;
		stats_t counts;
		uint32_t sug_fix;
	} history_counts_t;

	vector<kmer_pack_t> v_kmers, v_kmers_new;
	vector<history_counts_t> v_history_counts;

	CKmer pmer_can;
	CKmer smer_can;
	CKmer bmer_can;
	vector<uint64_t> v_minim_cand, v_minim_top_cand;
	vector<pair<uint64_t, uint32_t>> v_read_bmers;

	CKmer pmer_can_unc;
	CKmer smer_can_unc;
	CKmer bmer_can_unc;

	uint32_t cor_pos;
	uint32_t N_run_len;
	uint32_t max_prefix_sorted_flag_value;
	uint32_t max_prefix_sorted_no_bytes_value;

	CKmer pmer_can_prev;

	uint8_t rank(const stats_t &counts, uint8_t sym) const;
	uint8_t un_rank(stats_t &counts, uint8_t sym);

	bool empty(const stats_t &counts) const {
		return !counts[0] && !counts[1] && !counts[2] && !counts[3];
	}
	bool non_empty(const stats_t &counts) const {
		return counts[0] || counts[1] || counts[2] || counts[3];
	}
	bool deterministic(const stats_t &counts) const {
		return (!counts[0] + !counts[1] + !counts[2] + !counts[3]) == 3;
	}
	void update_ctx_r_sym(uint64_t &ctx_r_sym, const int r_sym = 4);
	void clear_counts(stats_t &counts) const {
		counts[0] = counts[1] = counts[2] = counts[3] = 0;
	}

	counts_level_t find_counts(stats_t &counts);
	bool find_counts_p(const CKmer &kmer, stats_t &counts);
	bool find_counts_s(const CKmer &kmer, stats_t &counts) {
		return ht_smer->find(kmer, counts, cinc_s);
	}

	bool find_counts_b(const CKmer &kmer, stats_t &counts) {
		return ht_bmer->find(kmer, counts, cinc_b);
	}

	bool find_counts_b_local(const CKmer &kmer, stats_t &counts) {
		return ht_bmer_local->find(kmer, counts, cinc_lb);
	}

	bool find_counts_s_local(const CKmer &kmer, stats_t &counts) {
		return ht_smer_local->find(kmer, counts, cinc_ls);
	}

	bool find_counts_rough_p(stats_t &counts);
	bool find_counts_rough_s(stats_t &counts);
	bool find_counts_rough_b(stats_t &counts);

	bool repair_kmers_existing(uint8_t *p_cur, uint32_t pos, stats_t &counts, uint8_t sym, uint32_t prev_count);
	bool repair_kmers_missing(uint8_t *p_cur, uint32_t pos, uint8_t sym, uint32_t existing_count = 0);

	void compress_prefix_direct(const uint8_t *p, uint32_t size, uint32_t start_pos);
	void compress_prefix_sorted(const uint8_t *p, uint32_t size);
	void compress_suffix(uint8_t *p, uint32_t size, bool original_order, uint32_t start_pos, bool reversed_pe = false);
	void append_pe_mers3(const uint8_t *p1, uint32_t size1, const uint8_t *p2, uint32_t size2);
	uint64_t find_minimizer(const uint8_t *p, uint32_t size) const;
	uint64_t find_maximizer(const uint8_t *p, uint32_t size) const;
	bool find_minim_cand(const uint8_t *p, uint32_t size);

	bool is_valid_minimizer(const uint64_t x) const;
	bool is_valid_maximizer(const uint64_t x) const;
	void merge_minim_results(vector<uint64_t> &v_cand, vector<uint64_t> &v_top_cand);
	void generate_read_bmers(const uint8_t *p, uint32_t size);
		
	void decompress_prefix_direct(uint8_t *p, uint32_t size);
	void decompress_prefix_sorted(uint8_t *p, uint32_t size);
	void decompress_suffix(uint8_t *p, uint32_t size, bool original_order, uint32_t start_pos, bool reversed_pe = false);

	void update_s_letters(const uint8_t *p, uint32_t size);

	void show_contexts();
	string decode_ctx_let(context_t ctx);

public:
	CDNACompressor();
	~CDNACompressor();

	bool Init();

	void SetKmerDS(TSmallIntVector<SIV_FIELD_SIZE> *_siv_pmer, 
		CHT_kmer<uint32_t> *_ht_smer, CHT_kmer<uint32_t> *_ht_bmer, CHT_pair_kmers *_ht_pe_mers, 
		vector<vector<vector<uint64_t>>> *_pmers_to_add, vector<vector<vector<uint64_t>>> *_smers_to_add, 
		vector<vector<vector<uint64_t>>> *_bmers_to_add, vector<vector<vector<tuple<uint64_t, uint64_t, uint64_t>>>> *_pe_mers_to_add,
		uint64_t _thread_id);

	bool CompressDirect(uint8_t *p, uint32_t size, uint8_t *q, bool first_read_of_pair);
	bool CompressDirectWithMinim(uint8_t *p, uint32_t size, uint32_t minim_pos);
	bool CompressSorted(uint8_t *p, uint32_t size, bool first_read_of_pair);
	bool CompressPE(uint8_t *p1, uint32_t size1, uint8_t *p2, uint32_t size2);
	bool DecompressDirectWithMinim(uint8_t *p, uint32_t size, uint32_t minim_pos, uint64_t minim);

	bool DecompressSE(uint8_t *p, uint32_t size);
	bool DecompressPE(uint8_t *p1, uint32_t size1, uint8_t *p2, uint32_t size2);

	void ResetReadPrev() { read_prev.clear(); };

	void InsertKmersToHT();
	void ClearKmersToHT();
	bool NeedRestructHT();
	void DoRestruct();
};

// EOF
