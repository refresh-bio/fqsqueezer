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
#include "params.h"
#include "compressor.h"
#include "io.h"
#include "rc.h"
#include "sub_rc.h"

class CMetaCompressor : public CBasicCompressor
{
	array<CRangeCoderModel<CVectorIOStream>*, 2> rc_meta;
	array<array<CRangeCoderModel<CVectorIOStream>*, 3>, 2> rc_meta_bytes;

	void encode_len(const uint32_t model, const uint32_t read_len);
	uint32_t decode_len(const uint32_t model);

public:
	CMetaCompressor();
	~CMetaCompressor();

	bool Init();

	bool CompressReadLen(const uint32_t read_len);
	bool CompressReadLenPE(const uint32_t read_len1, const uint32_t read_len2);

	bool DecompressReadLen(uint32_t &read_len);
	bool DecompressReadLenPE(uint32_t &read_len1, uint32_t &read_len2);
};

// EOF
