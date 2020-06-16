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
#include "compressor.h"
#include "params.h"
#include "io.h"
#include "rc.h"
#include "sub_rc.h"
#include "mtf.h"

#include <vector>
#include <array>

class CIdCompressor : public CBasicCompressor
{
	CMTF mtf;

	enum class token_type_t {literal = 0, numeric = 1};
	array<bool, 256> a_numeric, a_literal;

	typedef vector<tuple<token_type_t, uint8_t, uint32_t, uint32_t>> v_tokens_t;
	v_tokens_t v_tokens_prev, v_tokens_cur;
	vector<int64_t> v_deltas;
	uint8_t *id_prev;
	uint32_t size_prev;

	CRangeCoderModel<CVectorIOStream> *rc_id_mtf_flag;
	CRangeCoderModel<CVectorIOStream> *rc_id_mtf_code[7];
	CRangeCoderModel<CVectorIOStream> *rc_id_mtf_byte[4];

	// lossless
	context_t ctx_flags;
	context_t ctx_pe_flags;

	ctx_map_2_t m_ctx_rc_flags;
	ctx_map_2_t m_ctx_rc_pe_flags;
	ctx_map_256_t m_ctx_rc_numeric_size;
	ctx_map_4_t m_ctx_rc_numeric_small;
	ctx_map_t m_ctx_rc_numeric;
	ctx_map_t m_ctx_rc_literal;
	ctx_map_2_t m_ctx_rc_literal_same;
	ctx_map_2_t m_ctx_rc_literal_same_length;
	ctx_map_128_t m_ctx_rc_plain;

	CRangeCoderModelFixedSize<CVectorIOStream, 2> *tpl_ctx_rc_flags;
	CRangeCoderModelFixedSize<CVectorIOStream, 2> *tpl_ctx_rc_pe_flags;

	CRangeCoderModelFixedSize<CVectorIOStream, 256> *tpl_ctx_rc_numeric_size;
	CRangeCoderModel<CVectorIOStream> *tpl_ctx_rc_numeric;
	CRangeCoderModelFixedSize<CVectorIOStream, 4> *tpl_ctx_rc_numeric_small;
	CRangeCoderModel<CVectorIOStream> *tpl_ctx_rc_literal;
	CRangeCoderModelFixedSize<CVectorIOStream, 2> *tpl_ctx_rc_literal_same;
	CRangeCoderModelFixedSize<CVectorIOStream, 2> *tpl_ctx_rc_literal_same_length;
	CRangeCoderModelFixedSize<CVectorIOStream, 128> *tpl_ctx_rc_plain;

	ctx_map_2_t::value_type find_rc_context(ctx_map_2_t &m_ctx_rc, context_t ctx, CRangeCoderModelFixedSize<CVectorIOStream, 2> *tpl, const stats_t *counts = nullptr);
	ctx_map_4_t::value_type find_rc_context(ctx_map_4_t &m_ctx_rc, context_t ctx, CRangeCoderModelFixedSize<CVectorIOStream, 4> *tpl, const stats_t *counts = nullptr);
	ctx_map_128_t::value_type find_rc_context(ctx_map_128_t &m_ctx_rc, context_t ctx, CRangeCoderModelFixedSize<CVectorIOStream, 128> *tpl, const stats_t *counts = nullptr);
	ctx_map_256_t::value_type find_rc_context(ctx_map_256_t &m_ctx_rc, context_t ctx, CRangeCoderModelFixedSize<CVectorIOStream, 256> *tpl, const stats_t *counts = nullptr);
	ctx_map_t::value_type find_rc_context(ctx_map_t &m_ctx_rc, context_t ctx, CRangeCoderModel<CVectorIOStream> *tpl, const stats_t *counts = nullptr);

	uint32_t extract_instrument(const uint8_t *p, uint32_t size);

	uint32_t tokenize(const uint8_t *p, const uint32_t size, v_tokens_t &v_tokens);

	bool typical_pe_ids(const uint8_t *p1, const uint32_t size1, const uint8_t *p2, const uint32_t size2) const;

	void compress_lossless(uint8_t *p, uint32_t size);
	void compress_instrument(uint8_t *p, uint32_t size);

	void decompress_none(uint8_t *p, uint32_t &size);
	void decompress_lossless(uint8_t *p, uint32_t &size);
	void decompress_instrument(uint8_t *p, uint32_t &size);

	bool is_numeric(uint8_t c)
	{
		return a_numeric[c];
	}

	bool is_literal(uint8_t c)
	{
		return a_literal[c];
	}

	void init_symbol_classes();
	
	bool token_types_same(const v_tokens_t &v1, const v_tokens_t &v2)
	{
		if (v1.size() != v2.size())
			return false;

		for (size_t i = 0; i < v1.size(); ++i)
			if (get<0>(v1[i]) != get<0>(v2[i]) || get<1>(v1[i]) != get<1>(v2[i]))
				return false;

		return true;
	}

	int64_t get_int(const uint8_t *p, const uint32_t begin, const uint32_t end)
	{
		int64_t r = 0;

		for (uint32_t i = begin; i < end; ++i)
			r = r * 10 + (int64_t)(p[i] - '0');

		return r;
	}

	int32_t store_int(uint8_t *p, const uint32_t begin, int64_t val)
	{
		int n_dig = 0;

		if (val < 10ll)						n_dig = 1;
		else if (val < 100ll)				n_dig = 2;
		else if (val < 1000ll)				n_dig = 3;
		else if (val < 10000ll)				n_dig = 4;
		else if (val < 100000ll)			n_dig = 5;
		else if (val < 1000000ll)			n_dig = 6;
		else if (val < 10000000ll)			n_dig = 7;
		else if (val < 100000000ll)			n_dig = 8;
		else if (val < 1000000000ll)		n_dig = 9;
		else if (val < 10000000000ll)		n_dig = 10;
		else if (val < 100000000000ll)		n_dig = 11;
		else if (val < 1000000000000ll)		n_dig = 12;
		else if (val < 10000000000000ll)	n_dig = 13;
		else if (val < 100000000000000ll)	n_dig = 14;
		else if (val < 1000000000000000ll)	n_dig = 15;
		else
		{
			cerr << "Error in store_int(): " << val << endl;
		}
			

		for (int i = n_dig-1; i >= 0; --i)
		{
			p[begin + i] = '0' + (val % 10);
			val /= 10;
		}

		return n_dig;
	}

public:
	CIdCompressor();
	~CIdCompressor();

	bool Init();

	bool Compress(uint8_t *p, uint32_t size);
	bool CompressPE(uint8_t *p1, uint32_t size1, uint8_t *p2, uint32_t size2);

	bool Decompress(uint8_t *p, uint32_t &size);
	bool DecompressPE(uint8_t *p, uint32_t &size1, uint32_t &size2);

	void ResetReadPrev();
};

// EOF
