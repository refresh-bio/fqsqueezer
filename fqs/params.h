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
#include <string>

using namespace std;

enum class quality_mode_t { lossless, illumina_8, illumina_4, binary, none };
enum class id_mode_t {lossless, instrument, none};
enum class dna_mode_t {se_original, se_sorted, pe_original, pe_sorted};
enum class work_mode_t {compress, decompress};

struct CParams
{
	work_mode_t work_mode;

	vector<string> v_file_names;
	string tmp_path;
	string out_path;
	string out_path2;

	size_t no_threads;

	size_t genome_size;

	dna_mode_t dna_mode;
	id_mode_t id_mode;
	quality_mode_t quality_mode;
	uint32_t quality_thr;

	bool duplicates_check;

	// internal params
	uint32_t prefix_len;
	uint32_t hmer_len;
	uint32_t bmer_len;
	uint32_t smer_len;
	uint32_t pmer_len;
	uint32_t imer_len;
	uint32_t ht_prefix_len;
	double ht_max_filling_factor;

	uint32_t verbosity;

	CParams()
	{
		work_mode = work_mode_t::compress;
		tmp_path = "./kcsd_tmp";
		out_path = "output.kcsd";
		out_path2 = "output2.kcsd";
		no_threads = 1;
		dna_mode = dna_mode_t::se_sorted;
		quality_mode = quality_mode_t::illumina_8;
		id_mode = id_mode_t::instrument;
		quality_thr = 20;
		duplicates_check = true;
		verbosity = 1;

		genome_size = 3100;

		// internal params
		prefix_len = 14;
		pmer_len = 13;
		bmer_len = 26;			
		hmer_len = 31;			
		smer_len = 15;
		imer_len = 26;
		ht_prefix_len = 11;		// Will be auto adjusted
		ht_max_filling_factor = 0.8;
	}

	void store_params(vector<uint8_t> &v_params)
	{
		v_params.push_back('K');
		v_params.push_back('C');
		v_params.push_back('S');
		v_params.push_back('D');

		v_params.push_back((uint8_t)no_threads);
		v_params.push_back((uint8_t)dna_mode);
		v_params.push_back((uint8_t)quality_mode);
		v_params.push_back((uint8_t)id_mode);
		v_params.push_back((uint8_t)quality_thr);
		v_params.push_back((uint8_t)duplicates_check);
		v_params.push_back((uint8_t)prefix_len);
		v_params.push_back((uint8_t)pmer_len);
		v_params.push_back((uint8_t)smer_len);
		v_params.push_back((uint8_t)bmer_len);
		v_params.push_back((uint8_t)imer_len);
		v_params.push_back((uint8_t)hmer_len);
		v_params.push_back((uint8_t)ht_prefix_len);
	}

	bool load_params(const vector<uint8_t> &v_params)
	{
		int i = 0;

		if (v_params.size() != 17)
			return false;

		if (v_params[i++] != 'K')	return false;
		if (v_params[i++] != 'C')	return false;
		if (v_params[i++] != 'S')	return false;
		if (v_params[i++] != 'D')	return false;

		no_threads = v_params[i++];
		dna_mode = (dna_mode_t)v_params[i++];
		quality_mode = (quality_mode_t)v_params[i++];
		id_mode = (id_mode_t)v_params[i++];
		quality_thr = v_params[i++];
		duplicates_check = (bool)v_params[i++];
		prefix_len = v_params[i++];
		pmer_len = v_params[i++];
		smer_len = v_params[i++];
		bmer_len = v_params[i++];
		imer_len = v_params[i++];
		hmer_len = v_params[i++];
		ht_prefix_len = v_params[i++];

		return true;
	}

	void adjust_kmer_sizes()
	{
		const vector<vector<uint32_t>> kmer_thr = { 
			{1, 9, 14, 17, 19, 18},
			{4, 9, 15, 18, 20, 19} ,
			{16, 10, 15, 18, 21, 20},
			{64, 11, 16, 18, 23, 21},
			{256, 12, 17, 20, 24, 22},
			{1024, 12, 17, 21, 26, 24},
			{4096, 13, 18, 21, 27, 24},
			{16384, 14, 18, 22, 27, 25},
			{65536, 15, 18, 22, 27, 25}
		};

		for(auto &x : kmer_thr)
			if (genome_size <= x[0])
			{
				prefix_len = x[1];
				pmer_len = x[2];
				smer_len = x[3];
				bmer_len = x[4];
				imer_len = x[5];
				break;
			}
	}
};

// EOF
