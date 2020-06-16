// *******************************************************************************************
// This file is a part of FQSqueezer software distributed under GNU GPL 3 licence.
// The homepage of the FQSqueezer project is http://sun.aei.polsl.pl/REFRESH/fqsqueezer
//
// Author: Sebastian Deorowicz
// Version: 1.1
// Date   : 2020-06-16
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
