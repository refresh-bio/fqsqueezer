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
#include <algorithm>
#include "defs.h"

using namespace std;

#pragma once

#include <vector>
#include <algorithm>
#include <array>
#include <string>
#include <map>
#include "queue.h"

using namespace std;

// *******************************************************************************************
//
// *******************************************************************************************
class CMTF
{
	vector<string> v;
	map<string, int> v_sym_pos;

	int mtf1_variant;

	int find_pos(const string x) const;

	void move_up(const int x);

public:
	CMTF();
	~CMTF();

	void Reset();
	void Insert(const string x);
	void Insert(const uint32_t pos);
	string GetString(const int x) const;
	int GetCode(const string x) const;
	bool IsPresent(const string x) const;
	int Size() const;
};

// EOF
