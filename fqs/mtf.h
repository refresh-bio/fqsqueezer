#pragma once
// *******************************************************************************************
// This file is a part of FQSqueezer software distributed under GNU GPL 3 licence.
// The homepage of the FQSqueezer project is http://sun.aei.polsl.pl/REFRESH/fqsqueezer
//
// Author: Sebastian Deorowicz
// Version: 1.1
// Date   : 2020-06-16
// *******************************************************************************************

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
	vector<std::string> v;
	map<std::string, int> v_sym_pos;

	int mtf1_variant;

	int find_pos(const std::string x) const;

	void move_up(const int x);

public:
	CMTF();
	~CMTF();

	void Reset();
	void Insert(const std::string x);
	void Insert(const uint32_t pos);
	std::string GetString(const int x) const;
	int GetCode(const std::string x) const;
	bool IsPresent(const std::string x) const;
	int Size() const;
};

// EOF
