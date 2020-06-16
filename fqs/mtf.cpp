// *******************************************************************************************
// This file is a part of FQSqueezer software distributed under GNU GPL 3 licence.
// The homepage of the FQSqueezer project is http://sun.aei.polsl.pl/REFRESH/fqsqueezer
//
// Author: Sebastian Deorowicz
// Version: 1.1
// Date   : 2020-06-16
// *******************************************************************************************

#include <iostream>
#include "mtf.h"
#include <numeric>
#include <cmath>

// *******************************************************************************************
// Constructor
CMTF::CMTF()
{
	mtf1_variant = false;
}

// *******************************************************************************************
// Destructor
CMTF::~CMTF()
{
}

// *******************************************************************************************
// Return current value for a symbol
int CMTF::find_pos(const std::string x) const
{
	auto p = v_sym_pos.find(x);

	if (p == v_sym_pos.end())
		return -1;

	return p->second;
}

// *******************************************************************************************
//
void CMTF::Reset()
{
	v.clear();
	v.shrink_to_fit();
	v_sym_pos.clear();
}

// *******************************************************************************************
// Update MTF list
void CMTF::Insert(const std::string x)
{
	int p = find_pos(x);

	if (p == -1)
	{
		v.push_back(x);
		p = (int) v.size() - 1;
		v_sym_pos[x] = p;
	}

	if (p > 0)
		move_up(p);
}

// *******************************************************************************************
void CMTF::Insert(const uint32_t pos)
{
	if (pos > 0)
		move_up(pos);
}

// *******************************************************************************************
// Move symbol to the begining of the list
inline void CMTF::move_up(const int x)
{
	auto sym = v[x];

	if (mtf1_variant)			// MTF-1 variant (only symbols from pos. 1 are moved to pos. 0; the rest are moved to pos. 1)
	{
		if (x == 1)
		{
			v[1] = v[0];
			v[0] = sym;
			v_sym_pos[v[0]] = 0;
			v_sym_pos[v[1]] = 1;
		}
		else
		{
			for (int i = x; i > 1; --i)
			{
				v[i] = v[(uint64_t) i - 1];
				v_sym_pos[v[i]] = i;
			}

			v[1] = sym;
			v_sym_pos[sym] = 1;
		}
	}
	else						// Classic MTF
	{
		for (int i = x; i > 0; --i)
		{
			v[i] = v[(uint64_t)i - 1];
			v_sym_pos[v[i]] = i;
		}

		v[0] = sym;
		v_sym_pos[sym] = 0;
	}
}

// *******************************************************************************************
// Return position of symbol in the list
int CMTF::GetCode(const std::string x) const
{
	return find_pos(x);
}

// *******************************************************************************************
// Gest symbol from the given position
std::string CMTF::GetString(const int x) const
{
	return v[x];
}

// *******************************************************************************************
// Check whether a symbol is valid
inline bool CMTF::IsPresent(const std::string x) const
{
	return find_pos(x) >= 0;
}

// *******************************************************************************************
// Return size of the MTF list
inline int CMTF::Size() const
{
	return (int) v.size();
}

// EOF
