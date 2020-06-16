// *******************************************************************************************
// This file is a part of FQSqueezer software distributed under GNU GPL 3 licence.
// The homepage of the FQSqueezer project is http://sun.aei.polsl.pl/REFRESH/fqsqueezer
//
// Author: Sebastian Deorowicz
// Version: 1.1
// Date   : 2020-06-16
// *******************************************************************************************

#include "utils.h"
#include <iostream>
#include <memory>
#include <sstream>

// ************************************************************************************
CProgress::CProgress()
{
	max_value = 0;
	show_comment = 0;
	current_value = 0;
}

// ************************************************************************************
CProgress::~CProgress()
{
}

// ************************************************************************************
void CProgress::Init(int64_t _max_value, bool _show_comment)
{
	unique_lock<mutex> lck(mtx);

	max_value = _max_value;
	current_value = 0;

	show_comment = _show_comment;
}

// ************************************************************************************
void CProgress::SetComment(std::string _comment)
{
	unique_lock<mutex> lck(mtx);

	comment = _comment;
}

// ************************************************************************************
void CProgress::Step(int64_t increment)
{
	unique_lock<mutex> lck(mtx);

	current_value += increment;

	show_progress();
}

// ************************************************************************************
void CProgress::show_progress(void)
{
	stringstream stext;

	stext.width(5);
	stext.precision(1);
	stext.setf(ios_base::fixed);

	stext << 100.0 * current_value / max_value;
	stext << "%";

	if (show_comment)
		stext << " " << comment;

	stext << "\r";

	std::string text = stext.str();

	if (text != prev_text)
	{
		cout << text;
		fflush(stdout);
		prev_text = text;
	}
}

// ************************************************************************************
std::string trim(std::string s)
{
	while (!s.empty() && (s.back() == '\n' || s.back() == '\r'))
		s.pop_back();

	return s;
}

// ************************************************************************************
uint64_t popcnt(uint64_t x)
{
	return _mm_popcnt_u64(x);
}

// ************************************************************************************
// No. of trailing 0's
uint64_t notz(uint64_t x)
{
	if (x == 0)
		return 8;
	return popcnt(~x & (x - 1ull));
}

// ************************************************************************************
void sort_copy_stats(stats_t &dest, const stats_t &src)
{
	int o0 = 0, o1 = 1, o2 = 2, o3 = 0;
	
	auto r0 = src[1] > src[0];	o0 += r0;	o1 -= r0;
	auto r1 = src[2] > src[0];	o0 += r1;	o2 -= r1;
	
	r0 = src[3] > src[0];	o0 += r0;
	r1 = src[2] > src[1];	o1 += r1;	o2 -= r1;
	r0 = src[3] > src[1];	o1 += r0;
	r1 = src[3] > src[2];	o2 += r1;
	o3 = 6 - o0 - o1 - o2;
	
	dest[o0] = src[0];
	dest[o1] = src[1];
	dest[o2] = src[2];
	dest[o3] = src[3];
}

// EOF
