// *******************************************************************************************
// This file is a part of FQSqueezer software distributed under GNU GPL 3 licence.
// The homepage of the FQSqueezer project is http://sun.aei.polsl.pl/REFRESH/fqsqueezer
//
// Author: Sebastian Deorowicz
// Version: 1.1
// Date   : 2020-06-16
// *******************************************************************************************

#include "meta.h"

//*****************************************************************************************************
CMetaCompressor::CMetaCompressor() : CBasicCompressor()
{
	for (int i = 0; i < 2; ++i)
	{
		rc_meta[i] = nullptr;
		for (int j = 0; j < 3; ++j)
			rc_meta_bytes[i][j] = nullptr;
	}
}

//*****************************************************************************************************
CMetaCompressor::~CMetaCompressor()
{
	for (auto x : rc_meta)
		delete x;

	for (auto &x : rc_meta_bytes)
		for (auto y : x)
			delete y;
}

//*****************************************************************************************************
bool CMetaCompressor::Init()
{
	rc_meta[0] = new CRangeCoderModel<CVectorIOStream>(rc, 256, 15, 1 << 15, nullptr, 1, params.work_mode == work_mode_t::compress);
	rc_meta[1] = new CRangeCoderModel<CVectorIOStream>(rc, 256, 15, 1 << 15, nullptr, 1, params.work_mode == work_mode_t::compress);

	for(int i = 0; i < 2; ++i)
		for(int j = 0; j < 3; ++j)
			rc_meta_bytes[i][j] = new CRangeCoderModel<CVectorIOStream>(rc, 256, 15, 1 << 15, nullptr, 1, params.work_mode == work_mode_t::compress);

	return true;
}

//*****************************************************************************************************
void CMetaCompressor::encode_len(const uint32_t model, const uint32_t read_len)
{
	if (read_len < 254)
		rc_meta[model]->Encode(read_len);
	else
	{
		if (read_len < 256 * 256)
		{
			rc_meta[model]->Encode(254);
			rc_meta_bytes[model][1]->Encode(read_len >> 8);
			rc_meta_bytes[model][2]->Encode(read_len & 0xff);
		}
		else if (read_len < 256 * 256 * 256)
		{
			rc_meta[model]->Encode(255);
			rc_meta_bytes[model][0]->Encode(read_len >> 16);
			rc_meta_bytes[model][1]->Encode((read_len >> 8) & 0xff);
			rc_meta_bytes[model][2]->Encode(read_len & 0xff);
		}
		else
		{
			cerr << "Unsupported read length: " << read_len << endl;
			exit(0);
		}
	}
}

//*****************************************************************************************************
uint32_t CMetaCompressor::decode_len(const uint32_t model)
{
	uint32_t flag = rc_meta[model]->Decode();

	if (flag < 254)
		return flag;

	uint32_t read_len = 0;
	if (flag == 254)
	{
		read_len = rc_meta_bytes[model][1]->Decode() << 8;
		read_len += rc_meta_bytes[model][2]->Decode();
	}
	else
	{
		read_len = rc_meta_bytes[model][0]->Decode() << 16;
		read_len += rc_meta_bytes[model][1]->Decode() << 8;
		read_len += rc_meta_bytes[model][2]->Decode();
	}

	return read_len;
}

//*****************************************************************************************************
bool CMetaCompressor::CompressReadLen(const uint32_t read_len)
{
	encode_len(0, read_len);

	return true;
}

//*****************************************************************************************************
bool CMetaCompressor::DecompressReadLen(uint32_t &read_len)
{
	read_len = decode_len(0);

	return true;
}

//*****************************************************************************************************
bool CMetaCompressor::CompressReadLenPE(const uint32_t read_len1, const uint32_t read_len2)
{
	encode_len(0, read_len1);
	encode_len(1, read_len2);

	return true;
}

//*****************************************************************************************************
bool CMetaCompressor::DecompressReadLenPE(uint32_t &read_len1, uint32_t &read_len2)
{
	read_len1 = decode_len(0);
	read_len2 = decode_len(1);

	return true;
}

// EOF
