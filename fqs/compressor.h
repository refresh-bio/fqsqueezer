#pragma once
// *******************************************************************************************
// This file is a part of FQSqueezer software distributed under GNU GPL 3 licence.
// The homepage of the MSAC project is http://sun.aei.polsl.pl/REFRESH/fqsqueezer
//
// Author: Sebastian Deorowicz
// Version: 1.0
// Date   : 2019-02-22
// *******************************************************************************************

#include <vector>

#include "defs.h"
#include "params.h"
#include "io.h"
#include "rc.h"
#include "sub_rc.h"
#include "context_hm.h"

using namespace std;

class CBasicCompressor
{
protected:
	CParams params;
	CBasicRangeCoder<CVectorIOStream> *rc;

	typedef CContextHM<CRangeCoderModel<CVectorIOStream>> ctx_map_t;
	typedef	CContextHM<CRangeCoderModelFixedSize<CVectorIOStream, 2>> ctx_map_2_t;
	typedef	CContextHM<CRangeCoderModelFixedSize<CVectorIOStream, 4>> ctx_map_4_t;
	typedef	CContextHM<CRangeCoderModelFixedSize<CVectorIOStream, 5>> ctx_map_5_t;
	typedef	CContextHM<CRangeCoderModelFixedSize<CVectorIOStream, 16>> ctx_map_16_t;
	typedef	CContextHM<CRangeCoderModelFixedSize<CVectorIOStream, 128>> ctx_map_128_t;
	typedef	CContextHM<CRangeCoderModelFixedSize<CVectorIOStream, 256>> ctx_map_256_t;

public:
	typedef vector<pair<ctx_map_t*, CRangeCoderModel<CVectorIOStream>**>> ctx_map_desc_t;

	CBasicCompressor() : rc(nullptr) {};
	virtual ~CBasicCompressor() {};

	void SetParams(const CParams &_params);
	void SetCoder(CBasicRangeCoder<CVectorIOStream> *_rc);
};

// EOF
