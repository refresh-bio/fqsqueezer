// *******************************************************************************************
// This file is a part of FQSqueezer software distributed under GNU GPL 3 licence.
// The homepage of the MSAC project is http://sun.aei.polsl.pl/REFRESH/fqsqueezer
//
// Author: Sebastian Deorowicz
// Version: 1.0
// Date   : 2019-02-22
// *******************************************************************************************

#include "quality.h"

//*****************************************************************************************************
CQualityCompressor::CQualityCompressor() : CBasicCompressor()
{
	tpl_ctx_rc = nullptr;

	for (int i = 0; i < 64; ++i)
	{
		quality_code_map_fwd[i] = 0;
		quality_code_map_rev[i] = 0;
	}
}

//*****************************************************************************************************
CQualityCompressor::~CQualityCompressor()
{
	if (tpl_ctx_rc)
		delete tpl_ctx_rc;
}

//*****************************************************************************************************
bool CQualityCompressor::Init()
{
	if (tpl_ctx_rc)
		delete tpl_ctx_rc;

	switch (params.quality_mode)
	{
	case quality_mode_t::lossless:
		tpl_ctx_rc = new CRangeCoderModel<CVectorIOStream>(rc, 64, 15, 1 << 15, nullptr, 1, params.work_mode == work_mode_t::compress);
		adjust_quality_map_lossless();
		no_bits_per_symbol = 6;
		no_ctx_symbols = 2;
		break;
	case quality_mode_t::illumina_8:
		tpl_ctx_rc = new CRangeCoderModel<CVectorIOStream>(rc, 8, 15, 1 << 15, nullptr, 1, params.work_mode == work_mode_t::compress);
		adjust_quality_map_illumina_8();
		no_bits_per_symbol = 3 + 1;
		no_ctx_symbols = 6;
		break;
	case quality_mode_t::illumina_4:
		tpl_ctx_rc = new CRangeCoderModel<CVectorIOStream>(rc, 4, 15, 1 << 15, nullptr, 1, params.work_mode == work_mode_t::compress);
		adjust_quality_map_illumina_4();
		no_bits_per_symbol = 2 + 1;
		no_ctx_symbols = 9;
		break;
	case quality_mode_t::binary:
		tpl_ctx_rc = new CRangeCoderModel<CVectorIOStream>(rc, 2, 15, 1 << 15, nullptr, 1, params.work_mode == work_mode_t::compress);
		adjust_quality_map_binary();
		no_bits_per_symbol = 1 + 1;
		no_ctx_symbols = 10;
		break;
	default:
		break;
	}

	ctx_mask = (1ull << (no_bits_per_symbol * no_ctx_symbols)) - 1ull;
	ctx_inc = 1ull << 48;

	return true;
}

//*****************************************************************************************************
void CQualityCompressor::adjust_quality_map_lossless()
{
	for (int i = 0; i < 64; ++i)
	{
		quality_code_map_fwd[i] = i;
		quality_code_map_rev[i] = i;
	}
}

//*****************************************************************************************************
void CQualityCompressor::adjust_quality_map_illumina_8()
{
	for(int i = 0; i < 2; ++i)
		quality_code_map_fwd[i] = 0;
	quality_code_map_rev[0] = 0;

	for (int i = 2; i < 10; ++i)
		quality_code_map_fwd[i] = 1;
	quality_code_map_rev[1] = 6;

	for (int i = 10; i < 20; ++i)
		quality_code_map_fwd[i] = 2;
	quality_code_map_rev[2] = 15;

	for (int i = 20; i < 25; ++i)
		quality_code_map_fwd[i] = 3;
	quality_code_map_rev[3] = 22;

	for (int i = 25; i < 30; ++i)
		quality_code_map_fwd[i] = 4;
	quality_code_map_rev[4] = 27;

	for (int i = 30; i < 35; ++i)
		quality_code_map_fwd[i] = 5;
	quality_code_map_rev[5] = 33;

	for (int i = 35; i < 40; ++i)
		quality_code_map_fwd[i] = 6;
	quality_code_map_rev[6] = 37;

	for (int i = 40; i < 64; ++i)
		quality_code_map_fwd[i] = 7;
	quality_code_map_rev[7] = 40;
}

//*****************************************************************************************************
void CQualityCompressor::adjust_quality_map_illumina_4()
{
	for (int i = 0; i < 2; ++i)
		quality_code_map_fwd[i] = 0;
	quality_code_map_rev[0] = 0;

	for (int i = 2; i < 15; ++i)
		quality_code_map_fwd[i] = 1;
	quality_code_map_rev[1] = 12;

	for (int i = 15; i < 31; ++i)
		quality_code_map_fwd[i] = 2;
	quality_code_map_rev[2] = 23;

	for (int i = 31; i < 64; ++i)
		quality_code_map_fwd[i] = 3;
	quality_code_map_rev[3] = 37;
}

//*****************************************************************************************************
void CQualityCompressor::adjust_quality_map_binary()
{
	for (uint32_t i = 0; i < params.quality_thr; ++i)
		quality_code_map_fwd[i] = 0;
	quality_code_map_rev[0] = 0;

	for (uint32_t i = params.quality_thr; i < 64u; ++i)
		quality_code_map_fwd[i] = 1;
	quality_code_map_rev[1] = params.quality_thr;
}

//*****************************************************************************************************
bool CQualityCompressor::Compress(const uint8_t *p, uint32_t size)
{
	if (params.quality_mode == quality_mode_t::none)
		return true;

	context_t context;

	reset_context(context);

	for (uint32_t i = 0; i < size; ++i)
	{
		auto p_rc = find_rc_context(context);
		auto q = quality_code_map_fwd[p[i] - 33];

		p_rc->Encode(q);

		update_context(context, q);
	}

	return true;
}

//*****************************************************************************************************
bool CQualityCompressor::Decompress(uint8_t *p, uint32_t size)
{
	if (params.quality_mode == quality_mode_t::none)
	{
		for (uint32_t i = 0; i < size; ++i)
			p[i] = 33 + params.quality_thr;

		return true;
	}

	context_t context;

	reset_context(context);

	for (uint32_t i = 0; i < size; ++i)
	{
		auto p_rc = find_rc_context(context);
		auto d = p_rc->Decode();
		auto q = quality_code_map_rev[d];
		p[i] = q + 33;

		update_context(context, d);
	}

	return true;
}

//*****************************************************************************************************
void CQualityCompressor::reset_context(context_t &ctx)
{
	ctx = ctx_mask;
}

//*****************************************************************************************************
void CQualityCompressor::update_context(context_t &ctx, uint8_t x)
{
	context_t my_ctx = ctx + ctx_inc;
	auto t = (ctx << no_bits_per_symbol) + (context_t)x;

	ctx = (my_ctx & ~ctx_mask) + (t & ctx_mask);
}

//*****************************************************************************************************
CQualityCompressor::ctx_map_t::value_type CQualityCompressor::find_rc_context(const context_t ctx)
{
	auto p = m_ctx_rc.find(ctx);

	if (p == nullptr)
		m_ctx_rc.insert(ctx, p = new CRangeCoderModel<CVectorIOStream>(*tpl_ctx_rc));

	return p;
}

// EOF
