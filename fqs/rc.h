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
#include "sub_rc.h"
#include <cmath>
#include <algorithm>
#include <numeric>

// *******************************************************************************************
//
// *******************************************************************************************
class CSimpleModel
{
	uint32_t n_symbols;
	uint32_t max_total;
	uint32_t *stats;
	uint32_t total;
	uint32_t adder;

	void rescale()
	{
		while (total >= max_total)
		{
			total = 0;
			for (uint32_t i = 0; i < n_symbols; ++i)
			{
				stats[i] = (stats[i] + 1) / 2;
				total += stats[i];
			}
		}
	}

public: 
	CSimpleModel(uint32_t _adder = 1) : n_symbols(0), max_total(0), stats(nullptr), total(0), adder(_adder)
	{
	};

	~CSimpleModel()
	{
		if (stats)
			delete[] stats;
	};

	CSimpleModel(const CSimpleModel &c) = delete;
	CSimpleModel& operator=(const CSimpleModel&) = delete;

	void Init(uint32_t _n_symbols, int *_init_stats, uint32_t _max_total, uint32_t _adder)
	{
		adder = _adder;

		if (stats)
		{
			if (n_symbols != _n_symbols)
			{
				delete[] stats;
				n_symbols = _n_symbols;
				stats = new uint32_t[n_symbols];
			}
		}
		else
		{
			n_symbols = _n_symbols;
			stats = new uint32_t[n_symbols];
		}

		max_total = _max_total;

		if (_init_stats)
			for (uint32_t i = 0; i < n_symbols; ++i)
				stats[i] = _init_stats[i];
		else
			fill_n(stats, n_symbols, 1);

		total = accumulate(stats, stats+n_symbols, 0u);
		rescale();
	}

	void Init(const CSimpleModel &c)
	{
		n_symbols = c.n_symbols;
		max_total = c.max_total;
		adder = c.adder;

		if (stats)
			delete[] stats;

		stats = new uint32_t[n_symbols];
		copy_n(c.stats, n_symbols, stats);
		total = accumulate(stats, stats + n_symbols, 0u);
	}

	void GetFreq(int symbol, int &sym_freq, int &left_freq, int &totf)
	{
		left_freq = 0;

		switch (symbol)
		{
			case 4: left_freq += stats[3];
			case 3: left_freq += stats[2];
			case 2: left_freq += stats[1];
			case 1: left_freq += stats[0];
			case 0: break;
			default:
				for (int i = 0; i < symbol; ++i)
					left_freq += stats[i];
		}

		sym_freq = stats[symbol];
		totf = total;
	}

	void Update(int symbol)
	{
		stats[symbol] += adder;
		total += adder;

		if (total >= max_total)
			rescale();
	}

	int GetSym(int left_freq)
	{
		int t = 0;

		for (uint32_t i = 0; i < n_symbols; ++i)
		{
			t += stats[i];
			if (t > left_freq)
				return i;
		}

		return -1;
	}

	uint32_t GetTotal()
	{
		return total;
	}

	void Merge(uint32_t *stats_to_merge)
	{
		for (uint32_t i = 0; i < n_symbols; ++i)
		{
			stats[i] += stats_to_merge[i];
			total += stats_to_merge[i];
		}
	}

	void CompleteMerge()
	{
		rescale();
	}

	uint32_t *GetStats()
	{
		return stats;
	}

	void SetStats(uint32_t *stats_to_set)
	{
		total = 0;
		for (uint32_t i = 0; i < n_symbols; ++i)
			total += stats[i] = stats_to_set[i];
	}
};

// *******************************************************************************************
//
// *******************************************************************************************
template <unsigned N_SYMBOLS> class CSimpleModelFixedSize
{
	uint32_t max_total;
	uint32_t stats[N_SYMBOLS];
	uint32_t total;
	uint32_t adder;
	size_t no_updates;

	void rescale()
	{
		while (total >= max_total)
		{
			total = 0;
			for (uint32_t i = 0; i < N_SYMBOLS; ++i)
			{
				stats[i] = (stats[i] + 1) / 2;
				total += stats[i];
			}
		}
	}

public:
	CSimpleModelFixedSize(uint32_t _adder = 1) : max_total(0), total(0), adder(_adder), no_updates(0)
	{
	};

	~CSimpleModelFixedSize()
	{
	};

	CSimpleModelFixedSize(const CSimpleModelFixedSize &c) = delete;
	CSimpleModelFixedSize& operator=(const CSimpleModelFixedSize&) = delete;

	void Init(const int *_init_stats, uint32_t _max_total, uint32_t _adder)
	{
		max_total = _max_total;
		adder = _adder;

		if (_init_stats)
			for (uint32_t i = 0; i < N_SYMBOLS; ++i)
				stats[i] = _init_stats[i];
		else
			fill_n(stats, N_SYMBOLS, 1);

		total = accumulate(stats, stats + N_SYMBOLS, 0u);
		rescale();

		no_updates = 0;
	}

	void Init(const CSimpleModelFixedSize &c)
	{
		max_total = c.max_total;
		adder = c.adder;

		copy_n(c.stats, N_SYMBOLS, stats);
		total = accumulate(stats, stats + N_SYMBOLS, 0u);

		no_updates = 0;
	}

	void GetFreq(int symbol, int &sym_freq, int &left_freq, int &totf)
	{
		left_freq = 0;

		if(N_SYMBOLS >= 5)
			for (int i = 0; i < symbol; ++i)
				left_freq += stats[i];
		else if(N_SYMBOLS > 2)
			switch (symbol)
			{
			case 4: left_freq += stats[3];
			case 3: left_freq += stats[2];
			case 2: left_freq += stats[1];
			case 1: left_freq += stats[0];
			case 0: break;
			}
		else
			switch (symbol)
			{
			case 1: left_freq += stats[0];
			case 0: break;
			}

		sym_freq = stats[symbol];
		totf = total;
	}

	void Update(int symbol)
	{
		stats[symbol] += adder;
		total += adder;

		if (total >= max_total)
			rescale();

		++no_updates;
	}

	int GetSym(int left_freq)
	{
		int t = 0;

		for (uint32_t i = 0; i < N_SYMBOLS; ++i)
		{
			t += stats[i];
			if (t > left_freq)
				return i;
		}

		return -1;
	}

	uint32_t GetTotal()
	{
		return total;
	}

	void Merge(uint32_t *stats_to_merge)
	{
		for (uint32_t i = 0; i < N_SYMBOLS; ++i)
		{
			stats[i] += stats_to_merge[i];
			total += stats_to_merge[i];
		}
	}

	void CompleteMerge()
	{
		rescale();
	}

	uint32_t *GetStats()
	{
		return stats;
	}

	void SetStats(uint32_t *stats_to_set)
	{
		total = 0;
		for (uint32_t i = 0; i < N_SYMBOLS; ++i)
			total += stats[i] = stats_to_set[i];
	}

	void GetLogStats(size_t &_no_updates, vector<float> &_v_freq)
	{
		_no_updates = no_updates;

		_v_freq.resize(N_SYMBOLS, 0u);
		auto sum = accumulate(stats, stats + N_SYMBOLS, 0u);

		if (sum)
			for (uint32_t i = 0; i < N_SYMBOLS; ++i)
				_v_freq[i] = (float)stats[i] / sum;
	}

	size_t GetNoUpdates()
	{
		return no_updates;
	}
};


// *******************************************************************************************
//
// *******************************************************************************************
template<typename T_IO_STREAM> class CRangeCoderModel
{
	CRangeEncoder<T_IO_STREAM> *rce;
	CRangeDecoder<T_IO_STREAM> *rcd;

	CSimpleModel simple_model;

	int no_symbols;
	int lg_totf;
	int totf;
	int rescale;
	uint32_t adder;
	bool compress;

public:
	CRangeCoderModel(CBasicRangeCoder<T_IO_STREAM> *rcb, int _no_symbols, int _lg_totf, int _rescale, int* _init, uint32_t _adder, bool _compress) :
		no_symbols(_no_symbols), lg_totf(_lg_totf), totf(1 << _lg_totf), rescale(_rescale), adder(_adder), compress(_compress)
	{
		simple_model.Init(no_symbols, _init, rescale, adder);

		if (compress)
		{
			rce = (CRangeEncoder<T_IO_STREAM>*) (rcb);
			rcd = nullptr;
		}
		else
		{
			rce = nullptr;
			rcd = (CRangeDecoder<T_IO_STREAM>*) (rcb);
		}
	}

	CRangeCoderModel(const CRangeCoderModel &c)
	{
		simple_model.Init(c.simple_model);
		rce = c.rce;
		rcd = c.rcd;

		no_symbols = c.no_symbols;
		lg_totf = c.lg_totf;
		totf = c.totf;
		rescale = c.rescale;
		compress = c.compress;
		adder = c.adder;
	}

	~CRangeCoderModel()
	{
	}

	void Encode(const int x)
	{
		int syfreq, ltfreq;
		simple_model.GetFreq(x, syfreq, ltfreq, totf);
		rce->EncodeFrequency(syfreq, ltfreq, totf);

		simple_model.Update(x);
	}

	int Decode()
	{
		int syfreq, ltfreq;

		totf = simple_model.GetTotal();
		ltfreq = (int) rcd->GetCumulativeFreq(totf);

		int x = simple_model.GetSym(ltfreq);

		simple_model.GetFreq(x, syfreq, ltfreq, totf);
		rcd->UpdateFrequency(syfreq, ltfreq, totf);
		simple_model.Update(x);

		return x;
	}

	CSimpleModel* GetSimpleModel()
	{
		return &simple_model;
	}

	void Init(int *init)
	{
		simple_model.Init(no_symbols, init, rescale, adder);
	}
};

// *******************************************************************************************
//
// *******************************************************************************************
template<typename T_IO_STREAM, unsigned N_SYMBOLS> class CRangeCoderModelFixedSize
{
	CRangeEncoder<T_IO_STREAM> *rce;
	CRangeDecoder<T_IO_STREAM> *rcd;

	CSimpleModelFixedSize<N_SYMBOLS> simple_model;

	bool compress;
	double est_tot_len;

public:
	CRangeCoderModelFixedSize(CBasicRangeCoder<T_IO_STREAM> *rcb, int _lg_totf, int _rescale, int* _init, uint32_t _adder, bool _compress) :
		compress(_compress)
	{
		simple_model.Init(_init, _rescale, _adder);

		if (compress)
		{
			rce = (CRangeEncoder<T_IO_STREAM>*) (rcb);
			rcd = nullptr;
		}
		else
		{
			rce = nullptr;
			rcd = (CRangeDecoder<T_IO_STREAM>*) (rcb);
		}

		est_tot_len = 0;
	}

	CRangeCoderModelFixedSize(const CRangeCoderModelFixedSize &c)
	{
		simple_model.Init(c.simple_model);
		rce = c.rce;
		rcd = c.rcd;

		compress = c.compress;

		est_tot_len = 0;
	}

	~CRangeCoderModelFixedSize()
	{
	}

	void Encode(const int x)
	{
		int syfreq, ltfreq, totf;

		simple_model.GetFreq(x, syfreq, ltfreq, totf);
		rce->EncodeFrequency(syfreq, ltfreq, totf);

		est_tot_len += rce->EstimateCodeLen(syfreq, totf);

		simple_model.Update(x);
	}

	int Decode()
	{
		int syfreq, ltfreq;
		int totf = simple_model.GetTotal();
		ltfreq = (int) rcd->GetCumulativeFreq(totf);

		int x = simple_model.GetSym(ltfreq);

		simple_model.GetFreq(x, syfreq, ltfreq, totf);
		rcd->UpdateFrequency(syfreq, ltfreq, totf);
		simple_model.Update(x);

		return x;
	}

	CSimpleModelFixedSize<N_SYMBOLS>* GetSimpleModel()
	{
		return &simple_model;
	}

	double GetEstTotLen()
	{
		return est_tot_len;
	}
};

// EOF
