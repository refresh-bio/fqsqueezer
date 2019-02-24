#pragma once
// *******************************************************************************************
// This file is a part of FQSqueezer software distributed under GNU GPL 3 licence.
// The homepage of the MSAC project is http://sun.aei.polsl.pl/REFRESH/fqsqueezer
//
// Author: Sebastian Deorowicz
// Version: 1.0
// Date   : 2019-02-22
// *******************************************************************************************

#include <string>
#include <vector>
#include <utility>
#include <thread>
#include "params.h"
#include "id.h"
#include "quality.h"
#include "dna.h"
#include "io.h"
#include "meta.h"
#include "rc.h"
#include "sub_rc.h"
#include "utils.h"
#include "reads_block.h"
#include "bit_vec.h"
#include "kmer.h"
#include "ht_kmer.h"

using namespace std;

class CApplication
{
	const uint32_t BUF_READ_SIZE = 128 << 10;
	const uint64_t READS_BLOCK_SIZE = 16ull << 20;

	char dna_convert[256];
	char dna_convert_NT[256];

	CParams params;
	CProgress progress;

	CInFile *inf;
	COutFile *outf;

	vector<CMetaCompressor> v_meta_comp;
	vector<CIdCompressor> v_id_comp;
	vector<CDNACompressor> v_dna_comp;
	vector<CQualityCompressor> v_quality_comp;
	
	// k-mer structures
	TSmallIntVector<SIV_FIELD_SIZE> *siv_pmer;
	CHT_kmer<uint32_t> *ht_smer;
	CHT_kmer<uint32_t> *ht_bmer;
	CHT_pair_kmers *ht_pe_mers;

	vector<vector<vector<uint64_t>>> bmers_to_add;
	vector<vector<vector<uint64_t>>> smers_to_add;
	vector<vector<vector<uint64_t>>> pmers_to_add;
	vector<vector<vector<tuple<uint64_t, uint64_t, uint64_t>>>> pe_mers_to_add;

	CReadsBlock<CBasicFASTQFile, COutFile> *reads_block_cur, *reads_block_next;
	vector<array<vector<uint8_t>, NO_STREAMS>> vv_uint8;
	vector<array<CVectorIOStream*, NO_STREAMS>> v_vios;

	bool compress_se_original();
	bool compress_se_sorted();
	bool compress_pe_original();
	bool compress_pe_sorted();

	bool decompress();

	bool decompress_se();
	bool decompress_pe();

	bool compress_se_files(const vector<string> &v_file_names, bool original_order);
	bool decompress_se_file(const string &file_name);

	bool compress_pe_files(const vector<string> &v_file_names, bool original_order);
	bool decompress_pe_file(const string &file_name1, const string &file_name2);

	bool preprocess_se();
	bool preprocess_pe();

	bool get_read(CInFile *f, uint8_t *buffer, uint8_t *&p_end_read, uint32_t &read_len, uint8_t *&p_read, uint8_t *&p_id, uint8_t *&p_dna, uint8_t *&p_quality);
	uint64_t calc_no_synchronizations(uint64_t generation, uint64_t no_reads, uint64_t no_threads)
	{
		uint64_t r = NormalizeValue<size_t>(generation < 100ull ? 100ull - generation : 0ull, 0ull, no_reads / no_threads / 2);
		if (r)
			--r;

		return r;
	}

	string bin_name(uint32_t no) const;
	string bin_name(uint32_t no, const uint32_t id_pair) const;

	void AdjustToParams();

	void log_thread_jobs(const string &str) const;
	
public:
	CApplication();
	~CApplication();

	void SetParams(const CParams &_params);

	bool Execute();
};

// EOF
