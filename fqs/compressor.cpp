// *******************************************************************************************
// This file is a part of FQSqueezer software distributed under GNU GPL 3 licence.
// The homepage of the MSAC project is http://sun.aei.polsl.pl/REFRESH/fqsqueezer
//
// Author: Sebastian Deorowicz
// Version: 1.0
// Date   : 2019-02-22
// *******************************************************************************************

#include "compressor.h"

//*****************************************************************************************************
void CBasicCompressor::SetParams(const CParams &_params)
{
	params = _params;
}

//*****************************************************************************************************
void CBasicCompressor::SetCoder(CBasicRangeCoder<CVectorIOStream> *_rc)
{
	rc = _rc;
}

// EOF
