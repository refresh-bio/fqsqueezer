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
#include "compressor.h"
#include "params.h"
#include "io.h"
#include "rc.h"
#include "sub_rc.h"

using namespace std;

class CQualityCompressor : public CBasicCompressor
{
	array<uint32_t, 96> quality_code_map_fwd;
	array<uint32_t, 96> quality_code_map_rev;

	uint32_t no_ctx_symbols;
	uint32_t no_bits_per_symbol;
	context_t ctx_mask;
	context_t ctx_inc;

	ctx_map_t m_ctx_rc;
	CRangeCoderModel<CVectorIOStream> *tpl_ctx_rc;

	void reset_context(context_t &ctx);
	void update_context(context_t &ctx, uint8_t x);

	ctx_map_t::value_type find_rc_context(const context_t ctx);

	void adjust_quality_map_lossless();
	void adjust_quality_map_illumina_8();
	void adjust_quality_map_illumina_4();
	void adjust_quality_map_binary();

public:
	CQualityCompressor();
	~CQualityCompressor();

	bool Init();

	bool Compress(const uint8_t *p, uint32_t size);
	bool Decompress(uint8_t *p, uint32_t size);
};

// EOF
