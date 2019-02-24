// *******************************************************************************************
// This file is a part of FQSqueezer software distributed under GNU GPL 3 licence.
// The homepage of the MSAC project is http://sun.aei.polsl.pl/REFRESH/fqsqueezer
//
// Author: Sebastian Deorowicz
// Version: 1.0
// Date   : 2019-02-22
// *******************************************************************************************

#include "application.h"
#include "defs.h"
#include "io.h"
#include "queue.h"
#include <iostream>
#include <sstream>
#include <chrono>
#include <iomanip>      

using namespace std;
using namespace std::chrono;

//#define LOG_WORKER_JOBS

//*****************************************************************************************************
CApplication::CApplication()
{
	for (int i = 0; i < 256; ++i)
		dna_convert[i] = 4;
	dna_convert['A'] = 0;
	dna_convert['C'] = 1;
	dna_convert['G'] = 2;
	dna_convert['T'] = 3;

	// N is treated as T. Just for distribution into bins
	for (int i = 0; i < 256; ++i)
		dna_convert_NT[i] = 3;
	dna_convert_NT['A'] = 0;
	dna_convert_NT['C'] = 1;
	dna_convert_NT['G'] = 2;
	dna_convert_NT['T'] = 3;

	inf = new CInFile;
	outf = new COutFile;

	reads_block_cur = nullptr;
	reads_block_next = nullptr;

	siv_pmer = nullptr;
	ht_smer = nullptr;
	ht_bmer = nullptr;
	ht_pe_mers = nullptr;
}

//*****************************************************************************************************
CApplication::~CApplication()
{
	delete inf;
	delete outf;

	if (reads_block_cur)
		delete reads_block_cur;
	if (reads_block_next)
		delete reads_block_next;

	if (siv_pmer)
		delete siv_pmer;

	if (ht_smer)
		delete ht_smer;
	if (ht_bmer)
		delete ht_bmer;
	if (ht_pe_mers)
		delete ht_pe_mers;
}

//*****************************************************************************************************
void CApplication::SetParams(const CParams &_params)
{
	params = _params;

	if (params.work_mode == work_mode_t::compress)
		AdjustToParams();
}

//*****************************************************************************************************
void CApplication::AdjustToParams()
{
	siv_pmer = new TSmallIntVector<SIV_FIELD_SIZE>(2 * params.pmer_len);
	ht_smer = new CHT_kmer<uint32_t>(params.smer_len, SMER_COUNTER_BITS, params.ht_max_filling_factor);
	ht_bmer = new CHT_kmer<uint32_t>(params.bmer_len, BMER_COUNTER_BITS, params.ht_max_filling_factor);
	ht_pe_mers = new CHT_pair_kmers(params.bmer_len, params.no_threads);

	pmers_to_add.resize(params.no_threads);
	for (auto &x : pmers_to_add)
		x.resize(params.no_threads);

	smers_to_add.resize(params.no_threads);
	for (auto &x : smers_to_add)
		x.resize(params.no_threads);

	bmers_to_add.resize(params.no_threads);
	for (auto &x : bmers_to_add)
		x.resize(params.no_threads);

	pe_mers_to_add.resize(params.no_threads);
	for (auto &x : pe_mers_to_add)
		x.resize(params.no_threads);
}

//*****************************************************************************************************
bool CApplication::Execute()
{
	if (params.v_file_names.empty())
		return false;

	if (params.work_mode == work_mode_t::compress)
	{
		if (params.dna_mode == dna_mode_t::se_original)
			return compress_se_original();
		else if (params.dna_mode == dna_mode_t::se_sorted)
			return compress_se_sorted();
		else if (params.dna_mode == dna_mode_t::pe_original)
			return compress_pe_original();
		else if (params.dna_mode == dna_mode_t::pe_sorted)
			return compress_pe_sorted();
	}
	else
		return decompress();

	return false;
}

//*****************************************************************************************************
bool CApplication::decompress()
{
	bool r = inf->Open(params.v_file_names.front(), params.verbosity);

	if (!r)
	{
		cerr << "Cannot open output file: " << params.v_file_names.front() << endl;
		return false;
	}

	// Load params
	vector<uint8_t> v_params;
	v_params.resize(inf->ReadUInt(1));
	inf->Read(v_params.data(), v_params.size());

	if (!params.load_params(v_params))
	{
		cerr << "Wrong compressed format\n";
		return false;
	}

	AdjustToParams();

	if (params.dna_mode == dna_mode_t::se_original || params.dna_mode == dna_mode_t::se_sorted)
	{
		if (!decompress_se_file(params.out_path))
		{
			r = false;
			cerr << "Failed when decompressing " << params.out_path << endl;
		}
	}
	else
	{
		if (!decompress_pe_file(params.out_path, params.out_path2))
		{
			r = false;
			cerr << "Failed when decompressing " << params.out_path << " and " << params.out_path2 << endl;
		}
	}

	inf->Close();

	return r;
}

//*****************************************************************************************************
bool CApplication::compress_se_original()
{
	bool r = outf->Open(params.out_path);

	if (!r)
	{
		cerr << "Cannot open output file: " << params.out_path << endl;
		return false;
	}

	compress_se_files(params.v_file_names, true);

	outf->Close();

	return r;
}

//*****************************************************************************************************
bool CApplication::compress_se_sorted()
{
	if (!preprocess_se())
		return false;

	bool r = outf->Open(params.out_path);

	if (!r)
	{
		cerr << "Cannot open output file: " << params.out_path << endl;
		return false;
	}

	compress_se_files(params.v_file_names, false);

	outf->Close();

	for (auto &x : params.v_file_names)
		remove(x.c_str());

	return r;
}

//*****************************************************************************************************
bool CApplication::decompress_se()
{
	bool r = inf->Open(params.v_file_names.front(), params.verbosity);

	if (!r)
	{
		cerr << "Cannot open output file: " << params.v_file_names.front() << endl;
		return false;
	}

	if (!decompress_se_file(params.out_path))
	{
		r = false;
		cerr << "Failed when compressing " << params.out_path << endl;
	}

	inf->Close();

	return r;
}

//*****************************************************************************************************
bool CApplication::compress_pe_original()
{
	bool r = outf->Open(params.out_path);

	if (!r)
	{
		cerr << "Cannot open output file: " << params.out_path << endl;
		return false;
	}

	compress_pe_files(params.v_file_names, true);

	outf->Close();

	return r;
}

//*****************************************************************************************************
bool CApplication::compress_pe_sorted()
{
	if (!preprocess_pe())
		return false;

	bool r = outf->Open(params.out_path);

	if (!r)
	{
		cerr << "Cannot open output file: " << params.out_path << endl;
		return false;
	}

	compress_pe_files(params.v_file_names, false);

	outf->Close();

	return r;
}

//*****************************************************************************************************
bool CApplication::decompress_pe()
{
	bool r = inf->Open(params.v_file_names.front(), params.verbosity);

	if (!r)
	{
		cerr << "Cannot open output file: " << params.v_file_names.front() << endl;
		return false;
	}

	if (!decompress_pe_file(params.out_path, params.out_path2))
	{
		r = false;
		cerr << "Failed when compressing " << params.out_path << endl;
	}

	inf->Close();

	return r;
}

//*****************************************************************************************************
// Load a single read and place it in a buffer
bool CApplication::get_read(CInFile *f, uint8_t *buffer, uint8_t *&p_end_read, uint32_t &read_len, uint8_t *&p_read, uint8_t *&p_id, uint8_t *&p_dna, uint8_t *&p_quality)
{
	int no_eols = 0;

	read_len = 0;
	p_id = buffer;
	p_read = buffer;

	while (!f->Eof() && no_eols < 4)
	{
		int c = f->Get();
		if (c == 0xA)
		{
			++no_eols;
			if (no_eols == 1)
				p_dna = buffer + 1;
			else if (no_eols == 3)
				p_quality = buffer + 1;
		}

		if (c != 0xA && no_eols == 1)			// DNA
			++read_len;

		*buffer++ = (uint8_t)c;
	}

	// Check whether read is correct
	bool is_strange = false;
	for (auto p = p_dna; p < p_dna + read_len; ++p)
		is_strange |= (*p != 'A' && *p != 'C' && *p != 'G' && *p != 'T' && *p != 'N');

	if (is_strange)
	{
		*(p_dna + read_len) = 0;
		cerr << "Strange: " << p_dna << endl;
	}

	p_end_read = buffer;

	return no_eols == 4;
}

//*****************************************************************************************************
bool CApplication::preprocess_se()
{
	CInFile inf;
	
	vector<COutFile*> v_outf(NO_BINS, nullptr);

	uint8_t buf[BUF_SIZE];

	for (uint32_t i = 0; i < NO_BINS; ++i)
	{
		v_outf[i] = new COutFile();
		if (!v_outf[i]->Open(bin_name(i)))
		{
			cerr << "Cannot open bin file: " << bin_name(i) << endl;
			return false;
		}
	}

	for (size_t i = 0; i < params.v_file_names.size(); ++i)
	{
		if (!inf.Open(params.v_file_names[i], params.verbosity))
		{
			cerr << "Cannot open file: " << params.v_file_names[i] << endl;
			return false;
		}

		uint8_t *p_read = nullptr;
		uint8_t *p_id = nullptr;
		uint8_t *p_dna = nullptr;
		uint8_t *p_quality = nullptr;
		uint8_t *p_end_read = nullptr;
		uint32_t read_len;

		while (!inf.Eof())
		{
			get_read(&inf, buf, p_end_read, read_len, p_read, p_id, p_dna, p_quality);

			int id_file = 0;

			if(id_file == 0)
				for (uint32_t i = 0; i < BIN_PREFIX_LEN; ++i)
				{
					id_file <<= 2;
					id_file += dna_convert_NT[p_dna[i]];
				}

			v_outf[id_file]->Write(buf, p_end_read - buf);
		}

		inf.Close();
	}

	params.v_file_names.clear();

	for (uint32_t i = 0; i < NO_BINS; ++i)
	{
		v_outf[i]->Close();
		delete v_outf[i];

		params.v_file_names.push_back(bin_name(i));
	}
	
	return true;
}

//*****************************************************************************************************
bool CApplication::preprocess_pe()
{
	CInFile inf_1, inf_2;

	vector<COutFile*> v_outf_1(NO_BINS, nullptr);
	vector<COutFile*> v_outf_2(NO_BINS, nullptr);

	uint8_t buf_1[BUF_SIZE];
	uint8_t buf_2[BUF_SIZE];

	for (uint32_t i = 0; i < NO_BINS; ++i)
	{
		v_outf_1[i] = new COutFile();
		if (!v_outf_1[i]->Open(bin_name(i, 1)))
		{
			cerr << "Cannot open bin file: " << bin_name(i, 1) << endl;
			return false;
		}

		v_outf_2[i] = new COutFile();
		if (!v_outf_2[i]->Open(bin_name(i, 2)))
		{
			cerr << "Cannot open bin file: " << bin_name(i, 2) << endl;
			return false;
		}
	}

	for (size_t i = 0; i < params.v_file_names.size(); i += 2)
	{
		if (!inf_1.Open(params.v_file_names[i], params.verbosity))
		{
			cerr << "Cannot open file: " << params.v_file_names[i] << endl;
			return false;
		}
		if (!inf_2.Open(params.v_file_names[i+1], params.verbosity))
		{
			cerr << "Cannot open file: " << params.v_file_names[i+1] << endl;
			return false;
		}

		uint8_t *p_read_1 = nullptr;
		uint8_t *p_id_1 = nullptr;
		uint8_t *p_dna_1 = nullptr;
		uint8_t *p_quality_1 = nullptr;
		uint8_t *p_end_read_1 = nullptr;
		uint32_t read_len_1;

		uint8_t *p_read_2 = nullptr;
		uint8_t *p_id_2 = nullptr;
		uint8_t *p_dna_2 = nullptr;
		uint8_t *p_quality_2 = nullptr;
		uint8_t *p_end_read_2 = nullptr;
		uint32_t read_len_2;

		while (!inf_1.Eof() && !inf_2.Eof())
		{
			get_read(&inf_1, buf_1, p_end_read_1, read_len_1, p_read_1, p_id_1, p_dna_1, p_quality_1);
			get_read(&inf_2, buf_2, p_end_read_2, read_len_2, p_read_2, p_id_2, p_dna_2, p_quality_2);

			int id_file = 0;

			if (id_file == 0)
				for (uint32_t i = 0; i < BIN_PREFIX_LEN; ++i)
				{
					id_file <<= 2;
					id_file += dna_convert_NT[p_dna_1[i]];
				}

			v_outf_1[id_file]->Write(buf_1, p_end_read_1 - buf_1);
			v_outf_2[id_file]->Write(buf_2, p_end_read_2 - buf_2);
		}

		inf_1.Close();
		inf_2.Close();
	}

	params.v_file_names.clear();

	for (uint32_t i = 0; i < NO_BINS; ++i)
	{
		v_outf_1[i]->Close();
		delete v_outf_1[i];

		v_outf_2[i]->Close();
		delete v_outf_2[i];

		params.v_file_names.push_back(bin_name(i, 1));
		params.v_file_names.push_back(bin_name(i, 2));
	}

	return true;
}

//*****************************************************************************************************
bool CApplication::compress_se_files(const vector<string> &v_file_names, bool original_order)
{
	reads_block_cur  = new CReadsBlock<CBasicFASTQFile, COutFile>(READS_BLOCK_SIZE, work_mode_t::compress);
	reads_block_next = new CReadsBlock<CBasicFASTQFile, COutFile>(READS_BLOCK_SIZE, work_mode_t::compress);

	CSemaphore sem_can_read;
	CSemaphore sem_can_process(1);
	CSemaphore sem_can_compress(1);
	CSemaphore sem_compressed;

	CBarrier bar_synchro((uint32_t) params.no_threads);
	size_t no_synchronizations = 0;

	vv_uint8.resize(params.no_threads);
	v_vios.resize(params.no_threads);
	for (size_t i = 0; i < params.no_threads; ++i)
		for(size_t j = 0; j < NO_STREAMS; ++j)
		v_vios[i][j] = new CVectorIOStream(vv_uint8[i][j]);

	v_meta_comp.resize(params.no_threads);
	v_id_comp.resize(params.no_threads);
	v_dna_comp.resize(params.no_threads);
	v_quality_comp.resize(params.no_threads);

	array<bool, NO_STREAMS> a_store_stream{ true, params.id_mode != id_mode_t::none, true, params.quality_mode != quality_mode_t::none };

	// Reading thread
	thread *thr_reading = new thread([&]
	{
		for (auto &fn : v_file_names)
		{
			CBasicFASTQFile *f_b_fastq;

			if(original_order)
				f_b_fastq = (CBasicFASTQFile*) new CPlainFASTQFile();
			else
				f_b_fastq = (CBasicFASTQFile*) new CSortedFASTQFile();

			if (!f_b_fastq->Open(fn, params.verbosity))
			{
				cerr << "Cannot open file: " << fn << "*" << endl;
				continue;
			}

			while (true)
			{
				sem_can_read.WaitForZero();		// wait for reading possibility
				reads_block_next->Clear();

				if (f_b_fastq->Eof())
					break;

				reads_block_next->Read(f_b_fastq);

				sem_can_read.Inc();				// block reading possibility
				sem_can_process.Dec();			// tell that block is ready for processing
			}
			
			f_b_fastq->Close();
			delete f_b_fastq;
		}

		sem_can_process.Dec();				// tell that block is ready for processing (in fact it is empty, which tells end-of-data)
	});

	// Compressing threads
	vector<thread*> v_thr_compress(params.no_threads);
	for (size_t thread_id = 0; thread_id < params.no_threads; ++thread_id)
	{
		v_thr_compress[thread_id] = new thread([thread_id, this, &sem_can_compress, &sem_compressed, &bar_synchro, &original_order]
		{
			array<CRangeEncoder<CVectorIOStream>, NO_STREAMS> c_rc_enc{
				CRangeEncoder<CVectorIOStream>(*v_vios[thread_id][0]),
				CRangeEncoder<CVectorIOStream>(*v_vios[thread_id][1]),
				CRangeEncoder<CVectorIOStream>(*v_vios[thread_id][2]),
				CRangeEncoder<CVectorIOStream>(*v_vios[thread_id][3]) };

			CMetaCompressor &meta_comp = v_meta_comp[thread_id];
			CIdCompressor &id_comp = v_id_comp[thread_id];
			CDNACompressor &dna_comp = v_dna_comp[thread_id];
			CQualityCompressor &quality_comp = v_quality_comp[thread_id];
			int generation = 0;

			meta_comp.SetParams(params);
			id_comp.SetParams(params);
			dna_comp.SetParams(params);
			quality_comp.SetParams(params);

			meta_comp.SetCoder(&c_rc_enc[STREAM_META]);
			meta_comp.Init();

			id_comp.SetCoder(&c_rc_enc[STREAM_ID]);
			id_comp.Init();

			dna_comp.SetCoder(&c_rc_enc[STREAM_DNA]);
			dna_comp.Init();
			dna_comp.SetKmerDS(siv_pmer, ht_smer, ht_bmer, ht_pe_mers, &pmers_to_add, &smers_to_add, &bmers_to_add, &pe_mers_to_add, thread_id);

			quality_comp.SetCoder(&c_rc_enc[STREAM_QUALITY]);
			quality_comp.Init();

			while (true)
			{
				sem_can_compress.WaitForZero(generation);

				if (reads_block_cur->Empty())
					break;

				uint64_t my_first = reads_block_cur->v_ranges[thread_id].first;
				uint64_t my_last = reads_block_cur->v_ranges[thread_id].second;
				uint64_t no_synchronizations = reads_block_cur->no_synchronizations;

				int synchro_generation = 0;
				uint64_t next_synchro = ((uint64_t) synchro_generation + 1) * (my_last - my_first) / (no_synchronizations + 1) + my_first;

				dna_comp.ResetReadPrev();
				id_comp.ResetReadPrev();

				for(auto &x : c_rc_enc)
					x.Start();

				for (uint64_t i = my_first; i < my_last; ++i)
				{
					auto &cur_read = reads_block_cur->v_reads[i];
					meta_comp.CompressReadLen(cur_read.read_len());
					
					id_comp.Compress(cur_read.id, (uint32_t) cur_read.id_len());

					if (original_order)
						dna_comp.CompressDirect(cur_read.dna, cur_read.read_len(), cur_read.quality, true);
					else
						dna_comp.CompressSorted(cur_read.dna, cur_read.read_len(), true);
					quality_comp.Compress(cur_read.quality, cur_read.read_len());

					if (i == next_synchro)
					{
						bar_synchro.count_down_and_wait();

						dna_comp.InsertKmersToHT();
						bar_synchro.count_down_and_wait();
						dna_comp.ClearKmersToHT();

						synchro_generation++;
						next_synchro = ((uint64_t) synchro_generation + 1) * (my_last - my_first) / (no_synchronizations + 1) + my_first;

						bar_synchro.count_down_and_wait();
					}
				}

				bar_synchro.count_down_and_wait();
				dna_comp.InsertKmersToHT();
				bar_synchro.count_down_and_wait();
				dna_comp.ClearKmersToHT();
				bar_synchro.count_down_and_wait();

				for(auto &x : c_rc_enc)
					x.End();

				sem_compressed.Dec(generation);
				++generation;
			}
		});
	}

	// Store params
	vector<uint8_t> v_params;
	params.store_params(v_params);
	outf->WriteUInt(v_params.size(), 1);
	outf->Write(v_params.data(), v_params.size());

	// Management
	int generation = 0;
	size_t total_stored = 0;
	size_t total_input = 0;
	array<size_t, NO_STREAMS> stream_stored = { 0, 0, 0, 0 };

	while (true)
	{
		sem_can_process.WaitForZero();
		swap(reads_block_cur, reads_block_next);
		sem_can_process.Inc();
		sem_can_read.Dec();

		if (reads_block_cur->Empty())
		{
			sem_can_compress.Dec_notify_all(generation);
			break;
		}

		reads_block_cur->PartitionForWorkers(params.no_threads);

		no_synchronizations = calc_no_synchronizations(generation, reads_block_cur->v_reads.size(), params.no_threads);
		reads_block_cur->no_synchronizations = no_synchronizations;
		
		sem_compressed.IncNum((int) params.no_threads, (int) generation);
		sem_can_compress.Dec_notify_all(generation);

		sem_compressed.WaitForZero(generation);

		if (reads_block_cur->filled_size)
		{
			outf->WriteUIntVar(reads_block_cur->v_reads.size());
			total_input += reads_block_cur->filled_size;
		}

		for (size_t i = 0; i < params.no_threads; ++i)
		{
			outf->WriteUIntVar(reads_block_cur->v_reads[reads_block_cur->v_ranges[i].first].id - reads_block_cur->v_reads[0].id);
			for (size_t j = 0; j < NO_STREAMS; ++j)
			{
				if (a_store_stream[j])
				{
					outf->WriteUIntVar(vv_uint8[i][j].size());
					outf->Write(vv_uint8[i][j].data(), vv_uint8[i][j].size());
					total_stored += vv_uint8[i][j].size();
					stream_stored[j] += vv_uint8[i][j].size();
				}
				vv_uint8[i][j].clear();
			}
		}

		++generation;
		sem_can_compress.Inc(generation);

		if (params.verbosity == 2)
		{
			cout << "Block no. " << generation << "  *  in_size: " << total_input << "   out_size: " << total_stored <<
				"   SIV aff: " << siv_pmer->avg_filling_factor() << "   "
				"   ratio: " << (1.0 * total_input / total_stored) << endl;
			cout << "      p_items: " << siv_pmer->get_no_pmers() << "   s_items: " << ht_smer->get_no_kmers() <<
				"   b_items: " << ht_bmer->get_no_kmers() << "   pe_items: " << ht_pe_mers->get_no_pairs() << endl;
			fflush(stdout);
		}
		else if (params.verbosity == 1)
		{
			cout << "Processed " << total_input << "      compression ratio: " << (1.0 * total_input / total_stored) << "\r";
			fflush(stdout);
		}
	}

	if (params.verbosity > 0)
	{
		if (params.verbosity == 1)
			cout << endl;
		cout << "Processed " << total_input << "      compression ratio: " << (1.0 * total_input / total_stored) << "\n";
		cout << "***** Streams: " << endl;
		cout << "  Meta    : " << stream_stored[STREAM_META] << endl;
		cout << "  Id      : " << stream_stored[STREAM_ID] << endl;
		cout << "  DNA     : " << stream_stored[STREAM_DNA] << endl;
		cout << "  Quality : " << stream_stored[STREAM_QUALITY] << endl;
		fflush(stdout);
	}

	thr_reading->join();
	delete thr_reading;

	for (auto x : v_thr_compress)
	{
		x->join();
		delete x;
	}
	v_thr_compress.clear();

	for (auto x : v_vios)
		for(auto y : x)
			delete y;
	v_vios.clear();
	vv_uint8.clear();

	v_meta_comp.clear();
	v_id_comp.clear();
	v_dna_comp.clear();
	v_quality_comp.clear();

	return true;
}

//*****************************************************************************************************
bool CApplication::decompress_se_file(const string &file_name)
{
	COutFile f_fastq;
	bool r = true;

	if (!f_fastq.Open(file_name))
		return false;

	reads_block_cur = new CReadsBlock<CBasicFASTQFile, COutFile>(READS_BLOCK_SIZE, work_mode_t::decompress);

	CSemaphore sem_can_decompress(1);
	CSemaphore sem_decompressed((int) params.no_threads);

	CBarrier bar_synchro((int) params.no_threads);
	size_t no_synchronizations = 0;

	vv_uint8.resize(params.no_threads);
	v_vios.resize(params.no_threads);
	for (size_t i = 0; i < params.no_threads; ++i)
		for(size_t j = 0; j < NO_STREAMS; ++j)
			v_vios[i][j] = new CVectorIOStream(vv_uint8[i][j]);

	v_meta_comp.resize(params.no_threads);
	v_id_comp.resize(params.no_threads);
	v_dna_comp.resize(params.no_threads);
	v_quality_comp.resize(params.no_threads);

	array<bool, NO_STREAMS> a_store_stream{ true, params.id_mode != id_mode_t::none, true, params.quality_mode != quality_mode_t::none };

	// Decompressing threads
	vector<thread*> v_thr_decompress(params.no_threads);
	for (size_t thread_id = 0; thread_id < params.no_threads; ++thread_id)
	{
		v_thr_decompress[thread_id] = new thread([thread_id, this, &sem_can_decompress, &sem_decompressed, &bar_synchro]
		{
			array<CRangeDecoder<CVectorIOStream>, NO_STREAMS> c_rc_dec{
				CRangeDecoder<CVectorIOStream>(*v_vios[thread_id][0]),
				CRangeDecoder<CVectorIOStream>(*v_vios[thread_id][1]),
				CRangeDecoder<CVectorIOStream>(*v_vios[thread_id][2]),
				CRangeDecoder<CVectorIOStream>(*v_vios[thread_id][3])};

			CMetaCompressor &meta_comp = v_meta_comp[thread_id];
			CIdCompressor &id_comp = v_id_comp[thread_id];
			CDNACompressor &dna_comp = v_dna_comp[thread_id];
			CQualityCompressor &quality_comp = v_quality_comp[thread_id];

			int generation = 0;

			meta_comp.SetParams(params);
			id_comp.SetParams(params);
			dna_comp.SetParams(params);
			quality_comp.SetParams(params);

			meta_comp.SetCoder(&c_rc_dec[STREAM_META]);
			meta_comp.Init();
			id_comp.SetCoder(&c_rc_dec[STREAM_ID]);
			id_comp.Init();
			dna_comp.SetCoder(&c_rc_dec[STREAM_DNA]);
			dna_comp.Init();
			dna_comp.SetKmerDS(siv_pmer, ht_smer, ht_bmer, ht_pe_mers, &pmers_to_add, &smers_to_add, &bmers_to_add, &pe_mers_to_add, thread_id);
			quality_comp.SetCoder(&c_rc_dec[STREAM_QUALITY]);
			quality_comp.Init();

			while (true)
			{
				sem_can_decompress.WaitForZero(generation);

				if (vv_uint8[thread_id][STREAM_META].empty())
					break;

				uint64_t my_first = reads_block_cur->v_ranges[thread_id].first;
				uint64_t my_last = reads_block_cur->v_ranges[thread_id].second;
				uint64_t my_offset = reads_block_cur->v_offsets[thread_id].first;
				uint64_t no_synchronizations = reads_block_cur->no_synchronizations;

				int synchro_generation = 0;
				uint64_t next_synchro = ((uint64_t) synchro_generation + 1) * (my_last - my_first) / (no_synchronizations + 1) + my_first;

				dna_comp.ResetReadPrev();
				id_comp.ResetReadPrev();

				for(auto &x : c_rc_dec)
					x.Start();

				uint8_t *p = reads_block_cur->output_FASTQ + my_offset;
				uint8_t *p_start = p;

				for (uint64_t i = my_first; i < my_last; ++i)
				{
					uint32_t read_len;
					uint32_t id_size = 0;
					meta_comp.DecompressReadLen(read_len);
					id_comp.Decompress(p, id_size);
					p += id_size;
					dna_comp.DecompressSE(p, read_len);
					p += read_len;
					*p++ = 0xA;
					*p++ = '+';
					*p++ = 0xA;
					quality_comp.Decompress(p, read_len);
					p += read_len;
					*p++ = 0xA;

					if (i >= next_synchro)
					{
						bar_synchro.count_down_and_wait();

						dna_comp.InsertKmersToHT();
						bar_synchro.count_down_and_wait();
						// Sprawdza czy bêdzie potrzebna restrukturyzacja pe-merów w którejœ czêœci
						bool need_restruct_pe = dna_comp.NeedRestructHT();

						bar_synchro.count_down_and_wait();
						if (need_restruct_pe)
							dna_comp.DoRestruct();
						dna_comp.ClearKmersToHT();

						synchro_generation++;
						next_synchro = ((uint64_t) synchro_generation + 1) * (my_last - my_first) / (no_synchronizations + 1) + my_first;

						bar_synchro.count_down_and_wait();
					}
				}

				reads_block_cur->v_offsets[thread_id].second = my_offset + (p - p_start);

				bar_synchro.count_down_and_wait();
				dna_comp.InsertKmersToHT();
				bar_synchro.count_down_and_wait();
				dna_comp.ClearKmersToHT();
				bar_synchro.count_down_and_wait();

				for(auto &x : c_rc_dec)
					x.End();

				sem_decompressed.Dec(generation);
				++generation;
			}
		});
	}

	// Management
	int generation = 0;
	while (!inf->Eof())
	{
		if (params.verbosity == 2)
		{
			cout << "Fpos: " << inf->GetPos() << " / " << inf->FileSize() << endl;
			fflush(stdout);
		}
		else if (params.verbosity == 1)
		{
			cout << "Decompressed: " << inf->GetPos() << " of " << inf->FileSize() << "\r";
			fflush(stdout);
		}
		
		uint64_t no_reads_in_block = inf->ReadUIntVar();

		reads_block_cur->v_offsets.resize(params.no_threads);

		for (size_t i = 0; i < params.no_threads; ++i)
		{
			reads_block_cur->v_offsets[i].first = inf->ReadUIntVar();
			for (size_t j = 0; j < NO_STREAMS; ++j)
			{
				vv_uint8[i][j].clear();
				if (a_store_stream[j])
				{
					v_vios[i][j]->RestartRead();
					vv_uint8[i][j].resize(inf->ReadUIntVar());
					inf->Read(vv_uint8[i][j].data(), vv_uint8[i][j].size());
				}
			}
		}

		if (params.verbosity == 2)
		{
			cout << "Block: " << generation << endl;
			fflush(stdout);
		}

		reads_block_cur->PartitionForWorkers(params.no_threads, no_reads_in_block);
		no_synchronizations = calc_no_synchronizations(generation, no_reads_in_block, params.no_threads);
		reads_block_cur->no_synchronizations = no_synchronizations;

		sem_can_decompress.Dec_notify_all(generation);

		sem_decompressed.WaitForZero(generation);
		++generation;

		sem_decompressed.IncNum((int) params.no_threads, (int) generation);
		sem_can_decompress.Inc(generation);

		for(size_t thr = 0; thr < params.no_threads; ++thr)
			f_fastq.Write(reads_block_cur->output_FASTQ + reads_block_cur->v_offsets[thr].first, 
				reads_block_cur->v_offsets[thr].second - reads_block_cur->v_offsets[thr].first);
	}
	
	if (params.verbosity == 1)
	{
		cout << "Decompressed: " << inf->GetPos() << " of " << inf->FileSize() << "\r";
		fflush(stdout);
	}

	for (size_t i = 0; i < params.no_threads; ++i)
		for(size_t j = 0; j < NO_STREAMS; ++j)
			vv_uint8[i][j].clear();
	sem_can_decompress.Dec_notify_all(generation);

	for (auto x : v_thr_decompress)
	{
		x->join();
		delete x;
	}
	v_thr_decompress.clear();

	v_meta_comp.clear();
	v_id_comp.clear();
	v_dna_comp.clear();
	v_quality_comp.clear();

	f_fastq.Close();

	if (params.verbosity == 1)
	{
		cout << endl;
		fflush(stdout);
	}

	return r;
}

//*****************************************************************************************************
bool CApplication::compress_pe_files(const vector<string> &v_file_names, bool original_order)
{
	reads_block_cur = new CReadsBlock<CBasicFASTQFile, COutFile>(READS_BLOCK_SIZE, work_mode_t::compress);
	reads_block_next = new CReadsBlock<CBasicFASTQFile, COutFile>(READS_BLOCK_SIZE, work_mode_t::compress);

	CSemaphore sem_can_read;
	CSemaphore sem_can_process(1);
	CSemaphore sem_can_compress(1);
	CSemaphore sem_compressed;

	CBarrier bar_synchro((uint32_t)params.no_threads);
	size_t no_synchronizations = 0;

	vv_uint8.resize(params.no_threads);
	v_vios.resize(params.no_threads);
	for (size_t i = 0; i < params.no_threads; ++i)
		for(size_t j = 0; j < NO_STREAMS; ++j)
			v_vios[i][j] = new CVectorIOStream(vv_uint8[i][j]);

	v_meta_comp.resize(params.no_threads);
	v_id_comp.resize(params.no_threads);
	v_dna_comp.resize(params.no_threads);
	v_quality_comp.resize(params.no_threads);

	array<bool, NO_STREAMS> a_store_stream{ true, params.id_mode != id_mode_t::none, true, params.quality_mode != quality_mode_t::none };

	// Reading thread
	thread *thr_reading = new thread([&]
	{
		for(size_t i_fn = 0; i_fn < v_file_names.size(); i_fn += 2)
		{
			auto fn_1 = v_file_names[i_fn];
			auto fn_2 = v_file_names[i_fn+1];

			CBasicFASTQFile *f_b_fastq_1;
			CBasicFASTQFile *f_b_fastq_2;

			if (original_order)
			{
				f_b_fastq_1 = (CBasicFASTQFile*) new CPlainFASTQFile();
				f_b_fastq_2 = (CBasicFASTQFile*) new CPlainFASTQFile();
			}
			else
			{
				f_b_fastq_1 = (CBasicFASTQFile*) new CSortedFASTQFile();
				f_b_fastq_2 = (CBasicFASTQFile*) new CSortedFASTQFile((CSortedFASTQFile*) f_b_fastq_1);
			}

			if (!f_b_fastq_1->Open(fn_1, params.verbosity))
			{
				cerr << "Cannot open file: " << fn_1 << "*" << endl;
				continue;
			}

			if (!f_b_fastq_2->Open(fn_2, params.verbosity))
			{
				f_b_fastq_1->Close();
				cerr << "Cannot open file: " << fn_2 << "*" << endl;
				continue;
			}

			while (true)
			{
				sem_can_read.WaitForZero();		// wait for reading possibility
				reads_block_next->Clear();

				if (f_b_fastq_1->Eof() || f_b_fastq_2->Eof())
					break;

				reads_block_next->Read(f_b_fastq_1, f_b_fastq_2);

				sem_can_read.Inc();				// block reading possibility
				sem_can_process.Dec();			// tell that block is ready for processing
			}

			f_b_fastq_1->Close();
			f_b_fastq_2->Close();
			delete f_b_fastq_1;
			delete f_b_fastq_2;
		}

		sem_can_process.Dec();				// tell that block is ready for processing (in fact it is empty, which tells end-of-data)
	});

	// Compressing threads
	vector<thread*> v_thr_compress(params.no_threads);
	for (size_t thread_id = 0; thread_id < params.no_threads; ++thread_id)
	{
		v_thr_compress[thread_id] = new thread([thread_id, this, &sem_can_compress, &sem_compressed, &bar_synchro, &original_order]
		{
			array<CRangeEncoder<CVectorIOStream>, 4> c_rc_enc{
				CRangeEncoder<CVectorIOStream>(*v_vios[thread_id][0]),
				CRangeEncoder<CVectorIOStream>(*v_vios[thread_id][1]),
				CRangeEncoder<CVectorIOStream>(*v_vios[thread_id][2]),
				CRangeEncoder<CVectorIOStream>(*v_vios[thread_id][3]) };

			CMetaCompressor &meta_comp = v_meta_comp[thread_id];
			CIdCompressor &id_comp = v_id_comp[thread_id];
			CDNACompressor &dna_comp = v_dna_comp[thread_id];
			CQualityCompressor &quality_comp = v_quality_comp[thread_id];
			int generation = 0;

			meta_comp.SetParams(params);
			id_comp.SetParams(params);
			dna_comp.SetParams(params);
			quality_comp.SetParams(params);

			meta_comp.SetCoder(&c_rc_enc[STREAM_META]);
			meta_comp.Init();
			id_comp.SetCoder(&c_rc_enc[STREAM_ID]);
			id_comp.Init();
			dna_comp.SetCoder(&c_rc_enc[STREAM_DNA]);
			dna_comp.Init();
			dna_comp.SetKmerDS(siv_pmer, ht_smer, ht_bmer, ht_pe_mers, &pmers_to_add, &smers_to_add, &bmers_to_add, &pe_mers_to_add, thread_id);

			quality_comp.SetCoder(&c_rc_enc[STREAM_QUALITY]);
			quality_comp.Init();

			while (true)
			{
				sem_can_compress.WaitForZero(generation);

				if (reads_block_cur->Empty())
					break;

				uint64_t my_first = reads_block_cur->v_ranges[thread_id].first;
				uint64_t my_last = reads_block_cur->v_ranges[thread_id].second;
				uint64_t no_synchronizations = reads_block_cur->no_synchronizations;

				int synchro_generation = 0;
				uint64_t next_synchro = ((uint64_t) synchro_generation + 1) * (my_last - my_first) / (no_synchronizations + 1) + my_first;

				dna_comp.ResetReadPrev();
				id_comp.ResetReadPrev();

				for(auto &x : c_rc_enc)
					x.Start();

				log_thread_jobs("Thread " + to_string(thread_id) + " : before block processing");
				for (uint64_t i = my_first; i < my_last; i += 2)
				{
					auto &cur_read_1 = reads_block_cur->v_reads[i];
					auto &cur_read_2 = reads_block_cur->v_reads[i+1];

					meta_comp.CompressReadLenPE(cur_read_1.read_len(), cur_read_2.read_len());
					id_comp.CompressPE(cur_read_1.id, cur_read_1.id_len(), cur_read_2.id, cur_read_2.id_len());
					dna_comp.CompressPE(cur_read_1.dna, cur_read_1.read_len(), cur_read_2.dna, cur_read_2.read_len());
					quality_comp.Compress(cur_read_1.quality, cur_read_1.read_len());
					quality_comp.Compress(cur_read_2.quality, cur_read_2.read_len());

					if (i >= next_synchro)
					{
						bar_synchro.count_down_and_wait();

						dna_comp.InsertKmersToHT();
						bar_synchro.count_down_and_wait();
						dna_comp.ClearKmersToHT();

						synchro_generation++;
						next_synchro = ((uint64_t) synchro_generation + 1) * (my_last - my_first) / (no_synchronizations + 1) + my_first;

						bar_synchro.count_down_and_wait();
					}
				}
				log_thread_jobs("Thread " + to_string(thread_id) + " : after block processing");

				bar_synchro.count_down_and_wait();
				log_thread_jobs("Thread " + to_string(thread_id) + " : before kmer inserting");
				dna_comp.InsertKmersToHT();
				log_thread_jobs("Thread " + to_string(thread_id) + " : after kmer inserting");
				bar_synchro.count_down_and_wait();
				log_thread_jobs("Thread " + to_string(thread_id) + " : before kmer clearing");
				dna_comp.ClearKmersToHT();
				log_thread_jobs("Thread " + to_string(thread_id) + " : after kmer inserting");
				bar_synchro.count_down_and_wait();

				for(auto &x : c_rc_enc)
					x.End();

				log_thread_jobs("Thread " + to_string(thread_id) + " : end of generation");
				sem_compressed.Dec(generation);
				++generation;
			}
		});
	}

	// Store params
	vector<uint8_t> v_params;
	params.store_params(v_params);
	outf->WriteUInt(v_params.size(), 1);
	outf->Write(v_params.data(), v_params.size());

	// Management
	//	int block_no = 0;
	int generation = 0;
	size_t total_stored = 0;
	array<size_t, 4> stream_stored = { 0, 0, 0, 0 };
	size_t total_input = 0;

	while (true)
	{
		sem_can_process.WaitForZero();
		log_thread_jobs("***** Main thread: before swap read block");
		swap(reads_block_cur, reads_block_next);
		sem_can_process.Inc();
		sem_can_read.Dec();

		if (reads_block_cur->Empty())
		{
			sem_can_compress.Dec_notify_all(generation);
			break;
		}

		log_thread_jobs("***** Main thread: before partition for workers");
		reads_block_cur->PartitionForWorkers(params.no_threads);

		no_synchronizations = calc_no_synchronizations(generation, reads_block_cur->v_reads.size(), params.no_threads);
		reads_block_cur->no_synchronizations = no_synchronizations;

		sem_compressed.IncNum((int)params.no_threads, (int)generation);
		sem_can_compress.Dec_notify_all(generation);

		sem_compressed.WaitForZero(generation);

		if (reads_block_cur->filled_size)
		{
			outf->WriteUIntVar(reads_block_cur->v_reads.size());
			total_input += reads_block_cur->filled_size;
		}

		log_thread_jobs("***** Main thread: before saving block data");
		for (size_t i = 0; i < params.no_threads; ++i)
		{
			outf->WriteUIntVar(reads_block_cur->v_reads[reads_block_cur->v_ranges[i].first].id - reads_block_cur->v_reads[0].id);
			for (size_t j = 0; j < NO_STREAMS; ++j)
			{
				if (a_store_stream[j])
				{
					outf->WriteUIntVar(vv_uint8[i][j].size());
					outf->Write(vv_uint8[i][j].data(), vv_uint8[i][j].size());
					total_stored += vv_uint8[i][j].size();
					stream_stored[j] += vv_uint8[i][j].size();
				}
				vv_uint8[i][j].clear();
			}
		}
		log_thread_jobs("***** Main thread: after saving block data");

		++generation;
		sem_can_compress.Inc(generation);

		if (params.verbosity == 2)
		{
			cout << "Block no. " << generation << "  *  in_size: " << total_input << "   out_size: " << total_stored <<
				"   SIV aff: " << siv_pmer->avg_filling_factor() << "   "
				"   ratio: " << (1.0 * total_input / total_stored) << endl;
			cout << "      p_items: " << siv_pmer->get_no_pmers() << "   s_items: " << ht_smer->get_no_kmers() <<
				"   b_items: " << ht_bmer->get_no_kmers() << "   pe_items: " << ht_pe_mers->get_no_pairs() << endl;
			fflush(stdout);
		}
		else if (params.verbosity == 1)
		{
			cout << "Processed " << total_input << "      compression ratio: " << (1.0 * total_input / total_stored) << "\r";
			fflush(stdout);
		}
	}

	log_thread_jobs("***** Main thread: before thr_reading->join()");

	thr_reading->join();
	delete thr_reading;

	log_thread_jobs("***** Main thread: before thr_compress->join()");
	for (auto x : v_thr_compress)
	{
		x->join();
		delete x;
	}
	v_thr_compress.clear();

	log_thread_jobs("***** Main thread: before v_vios delete");
	for (auto x : v_vios)
		for(auto y : x)
			delete y;
	v_vios.clear();
	vv_uint8.clear();

	log_thread_jobs("***** Main thread: before v_meta_comp clear()");
	v_meta_comp.clear();
	log_thread_jobs("***** Main thread: before v_id_comp clear()");
	v_id_comp.clear();
	log_thread_jobs("***** Main thread: before v_dna_comp clear()");
	v_dna_comp.clear();
	log_thread_jobs("***** Main thread: before v_quality_comp clear()");
	v_quality_comp.clear();
	log_thread_jobs("***** Main thread: after v_* clear()");

	if (params.verbosity > 0)
	{
		if (params.verbosity == 1)
			cout << endl;
		cout << "Processed " << total_input << "      compression ratio: " << (1.0 * total_input / total_stored) << "\n";
		cout << "***** Streams: " << endl;
		cout << "  Meta    : " << stream_stored[STREAM_META] << endl;
		cout << "  Id      : " << stream_stored[STREAM_ID] << endl;
		cout << "  DNA     : " << stream_stored[STREAM_DNA] << endl;
		cout << "  Quality : " << stream_stored[STREAM_QUALITY] << endl;
		fflush(stdout);
	}

	return true;
}

//*****************************************************************************************************
bool CApplication::decompress_pe_file(const string &file_name1, const string &file_name2)
{
	COutFile f_fastqs[2];
	bool r = true;

	if (!f_fastqs[0].Open(file_name1))
	{
		cerr << "Cannot open " << file_name1 << endl;
		return false;
	}
	if (!f_fastqs[1].Open(file_name2))
	{
		cerr << "Cannot open " << file_name2 << endl;
		f_fastqs[0].Close();
		return false;
	}

	reads_block_cur = new CReadsBlock<CBasicFASTQFile, COutFile>(READS_BLOCK_SIZE, work_mode_t::decompress);

	CSemaphore sem_can_decompress(1);
	CSemaphore sem_decompressed((int)params.no_threads);

	CBarrier bar_synchro((int)params.no_threads);
	size_t no_synchronizations = 0;

	vv_uint8.resize(params.no_threads);
	v_vios.resize(params.no_threads);
	for (size_t i = 0; i < params.no_threads; ++i)
		for (size_t j = 0; j < NO_STREAMS; ++j)
			v_vios[i][j] = new CVectorIOStream(vv_uint8[i][j]);

	v_meta_comp.resize(params.no_threads);
	v_id_comp.resize(params.no_threads);
	v_dna_comp.resize(params.no_threads);
	v_quality_comp.resize(params.no_threads);

	array<bool, NO_STREAMS> a_store_stream{ true, params.id_mode != id_mode_t::none, true, params.quality_mode != quality_mode_t::none };

	// Decompressing threads
	vector<thread*> v_thr_decompress(params.no_threads);
	for (size_t thread_id = 0; thread_id < params.no_threads; ++thread_id)
	{
		v_thr_decompress[thread_id] = new thread([thread_id, this, &sem_can_decompress, &sem_decompressed, &bar_synchro]
		{
			array<CRangeDecoder<CVectorIOStream>, NO_STREAMS> c_rc_dec{
				CRangeDecoder<CVectorIOStream>(*v_vios[thread_id][0]),
				CRangeDecoder<CVectorIOStream>(*v_vios[thread_id][1]),
				CRangeDecoder<CVectorIOStream>(*v_vios[thread_id][2]),
				CRangeDecoder<CVectorIOStream>(*v_vios[thread_id][3]) };

			CMetaCompressor &meta_comp = v_meta_comp[thread_id];
			CIdCompressor &id_comp = v_id_comp[thread_id];
			CDNACompressor &dna_comp = v_dna_comp[thread_id];
			CQualityCompressor &quality_comp = v_quality_comp[thread_id];

			int generation = 0;

			meta_comp.SetParams(params);
			id_comp.SetParams(params);
			dna_comp.SetParams(params);
			quality_comp.SetParams(params);

			meta_comp.SetCoder(&c_rc_dec[STREAM_META]);
			meta_comp.Init();
			id_comp.SetCoder(&c_rc_dec[STREAM_ID]);
			id_comp.Init();
			dna_comp.SetCoder(&c_rc_dec[STREAM_DNA]);
			dna_comp.Init();
			dna_comp.SetKmerDS(siv_pmer, ht_smer, ht_bmer, ht_pe_mers, &pmers_to_add, &smers_to_add, &bmers_to_add, &pe_mers_to_add, thread_id);
			quality_comp.SetCoder(&c_rc_dec[STREAM_QUALITY]);
			quality_comp.Init();

			while (true)
			{
				sem_can_decompress.WaitForZero(generation);

				if (vv_uint8[thread_id][STREAM_META].empty())
					break;

				uint64_t my_first = reads_block_cur->v_ranges[thread_id].first;
				uint64_t my_last = reads_block_cur->v_ranges[thread_id].second;
				uint64_t my_offset = reads_block_cur->v_offsets[thread_id].first;
				uint64_t no_synchronizations = reads_block_cur->no_synchronizations;

				int synchro_generation = 0;
				uint64_t next_synchro = ((uint64_t) synchro_generation + 1) * (my_last - my_first) / (no_synchronizations + 1) + my_first;

				dna_comp.ResetReadPrev();
				id_comp.ResetReadPrev();

				for (auto &x : c_rc_dec)
					x.Start();

				uint8_t *p = reads_block_cur->output_FASTQ + my_offset;
				uint8_t *p_start = p;

				for (uint64_t i = my_first; i < my_last; i += 2)
				{
					uint32_t read_len1, read_len2;
					uint32_t id_size1 = 0;
					uint32_t id_size2 = 0;

					meta_comp.DecompressReadLenPE(read_len1, read_len2);
					id_comp.DecompressPE(p, id_size1, id_size2);
					p += id_size1;
					p += id_size2;
					dna_comp.DecompressPE(p, read_len1, p + read_len1 + 1, read_len2);
					p += read_len1;
					*p++ = 0xA;
					p += read_len2;
					*p++ = 0xA;
					*p++ = '+';
					*p++ = 0xA;
					*p++ = '+';
					*p++ = 0xA;
					quality_comp.Decompress(p, read_len1);
					p += read_len1;
					*p++ = 0xA;
					quality_comp.Decompress(p, read_len2);
					p += read_len2;
					*p++ = 0xA;

					if (i >= next_synchro)
					{
						bar_synchro.count_down_and_wait();

						dna_comp.InsertKmersToHT();
						bar_synchro.count_down_and_wait();
						dna_comp.ClearKmersToHT();

						synchro_generation++;
						next_synchro = ((uint64_t) synchro_generation + 1) * (my_last - my_first) / (no_synchronizations + 1) + my_first;

						bar_synchro.count_down_and_wait();
					}
				}

				reads_block_cur->v_offsets[thread_id].second = my_offset + (p - p_start);

				bar_synchro.count_down_and_wait();
				dna_comp.InsertKmersToHT();
				bar_synchro.count_down_and_wait();
				dna_comp.ClearKmersToHT();
				bar_synchro.count_down_and_wait();

				for (auto &x : c_rc_dec)
					x.End();

				sem_decompressed.Dec(generation);
				++generation;
			}
		});
	}

	// Management
	int generation = 0;
	while (!inf->Eof())
	{
		if (params.verbosity == 2)
		{
			cout << "Fpos: " << inf->GetPos() << " / " << inf->FileSize() << endl;
			fflush(stdout);
		}
		else if (params.verbosity == 1)
		{
			cout << "Decompressed: " << inf->GetPos() << " of " << inf->FileSize() << "\r";
			fflush(stdout);
		}

		uint64_t no_reads_in_block = inf->ReadUIntVar();

		reads_block_cur->v_offsets.resize(params.no_threads);

		for (size_t i = 0; i < params.no_threads; ++i)
		{
			reads_block_cur->v_offsets[i].first = inf->ReadUIntVar();
			for (size_t j = 0; j < NO_STREAMS; ++j)
			{
				vv_uint8[i][j].clear();
				if (a_store_stream[j])
				{
					v_vios[i][j]->RestartRead();
					vv_uint8[i][j].resize(inf->ReadUIntVar());
					inf->Read(vv_uint8[i][j].data(), vv_uint8[i][j].size());
				}
			}
		}

		if (params.verbosity == 2)
		{
			cout << "Block: " << generation << endl;
			fflush(stdout);
		}

		reads_block_cur->PartitionForWorkers(params.no_threads, no_reads_in_block);
		no_synchronizations = calc_no_synchronizations(generation, no_reads_in_block, params.no_threads);
		reads_block_cur->no_synchronizations = no_synchronizations;

		sem_can_decompress.Dec_notify_all(generation);

		sem_decompressed.WaitForZero(generation);
		++generation;

		sem_decompressed.IncNum((int)params.no_threads, (int)generation);
		sem_can_decompress.Inc(generation);

		// Split and write
		int fastq_idx = 0;
		auto p = reads_block_cur->output_FASTQ;

		for (size_t thr = 0; thr < params.no_threads; ++thr)
		{
			size_t line_start_idx = reads_block_cur->v_offsets[thr].first;
			for(size_t i = reads_block_cur->v_offsets[thr].first; i < reads_block_cur->v_offsets[thr].second; ++i)
			{
				if (p[i] == 0xA)
				{
					f_fastqs[fastq_idx].Write(p + line_start_idx, i + 1 - line_start_idx);
					line_start_idx = i + 1;
					fastq_idx = !fastq_idx;
				}
			}
		}
	}

	if (params.verbosity == 1)
	{
		cout << "Decompressed: " << inf->GetPos() << " of " << inf->FileSize() << "\r";
		fflush(stdout);
	}
	
	for (size_t i = 0; i < params.no_threads; ++i)
		for (size_t j = 0; j < NO_STREAMS; ++j)
			vv_uint8[i][j].clear();
	sem_can_decompress.Dec_notify_all(generation);

	for (auto x : v_thr_decompress)
	{
		x->join();
		delete x;
	}
	v_thr_decompress.clear();

	v_meta_comp.clear();
	v_id_comp.clear();
	v_dna_comp.clear();
	v_quality_comp.clear();

	f_fastqs[0].Close();
	f_fastqs[1].Close();

	if (params.verbosity == 1)
	{
		cout << endl;
		fflush(stdout);
	}

	return r;
}

//*****************************************************************************************************
string CApplication::bin_name(uint32_t no) const
{
	char s[10];

	int log4_n_bins = (int) ilog2(NO_BINS - 1) / 2;

	s[log4_n_bins] = 0;
	for (int i = log4_n_bins - 1; i >= 0; --i)
	{
		s[i] = dna_alphabet[no % 4];
		no /= 4;
	}

	return params.tmp_path + "kcsd_" + string(s) + ".bin";
}

//*****************************************************************************************************
string CApplication::bin_name(uint32_t no, const uint32_t id_pair) const
{
	char s[10];

	int log4_n_bins = (int)ilog2(NO_BINS - 1) / 2;

	s[log4_n_bins] = 0;
	for (int i = log4_n_bins - 1; i >= 0; --i)
	{
		s[i] = dna_alphabet[no % 4];
		no /= 4;
	}

	return params.tmp_path + "kcsd_" + string(s) + "_" + to_string(id_pair) + ".bin";
}

//*****************************************************************************************************
void CApplication::log_thread_jobs(const string &str) const
{
#ifndef LOG_WORKER_JOBS
	return;
#endif
	stringstream tmp_str;

	auto timeASMs = TimeHelpers::TimeFromEpochInMilliSeconds<uint64_t>();

	auto tt = system_clock::to_time_t(chrono::system_clock::now());
	tmp_str << ctime(&tt) << "   " << std::fixed << std::setprecision(3) << timeASMs / 1000.0 << " : " << str << endl;

	auto ts = tmp_str.str();
	cerr << ts;
}

// EOF
