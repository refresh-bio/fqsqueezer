#pragma once
// *******************************************************************************************
// This file is a part of FQSqueezer software distributed under GNU GPL 3 licence.
// The homepage of the MSAC project is http://sun.aei.polsl.pl/REFRESH/fqsqueezer
//
// Author: Sebastian Deorowicz
// Version: 1.0
// Date   : 2019-02-22
// *******************************************************************************************

#include "defs.h"
#include "params.h"
#include "io.h"
#include <vector>

using namespace std;

template<typename IN_FILE_TYPENAME, typename OUT_FILE_TYPENAME> class CReadsBlock
{
public:
	uint64_t block_size;
	work_mode_t work_mode;
	uint64_t filled_size;

	const uint64_t block_size_margin = 10240;

	uint8_t *input_FASTQ;
	uint8_t *output_FASTQ;

	vector<read_desc_t> v_reads;
	vector<pair<uint64_t, uint64_t>> v_ranges;
	vector<pair<uint64_t, uint64_t>> v_offsets;
	uint64_t no_synchronizations;

	bool get_read(IN_FILE_TYPENAME *f, uint8_t *buffer, read_desc_t &read_desc)
	{
		int no_eols = 0;

		int read_len = 0;

		uint8_t *p = buffer;

		read_desc.id = p;

		while (!f->Eof() && no_eols < 4)
		{
			int c = f->Get();
			if (c == 0xA)
			{
				++no_eols;
				if (no_eols == 1)
					read_desc.dna = p + 1;
				else if (no_eols == 2)
					read_desc.plus = p + 1;
				else if (no_eols == 3)
					read_desc.quality = p + 1;
			}

			if (c != 0xA && no_eols == 1)			// DNA
				++read_len;

			*p++ = c;
		}

		// Verify read validity
		bool is_strange = false;
		for (auto p = read_desc.dna; p + 1 < read_desc.plus; ++p)
			is_strange |= (*p != 'A' && *p != 'C' && *p != 'G' && *p != 'T' && *p != 'N');

		if (is_strange)
		{
			read_desc.plus = 0;
			cout << "Strange: " << read_desc.dna << endl;
		}
		return no_eols == 4;
	}

public:
	CReadsBlock(uint64_t _block_size, work_mode_t _work_mode)
	{
		block_size = _block_size;
		work_mode = _work_mode;

		no_synchronizations = 0;

		if (work_mode == work_mode_t::compress)
		{
			input_FASTQ = new uint8_t[block_size];
			output_FASTQ = nullptr;
		}
		else
		{
			input_FASTQ = nullptr;
			output_FASTQ = new uint8_t[block_size];
#ifdef DECOMPRESS_ONLY_DNA
			fill_n(output_FASTQ, block_size, 0);
#endif
		}

		filled_size = 0;
	}

	~CReadsBlock()
	{
		if (input_FASTQ)
			delete[] input_FASTQ;
		if (output_FASTQ)
			delete[] output_FASTQ;
	}

	void Clear()
	{
		filled_size = 0;
	}

	bool Empty()
	{
		return filled_size == 0;
	}

	bool Read(IN_FILE_TYPENAME *in_f)
	{
		read_desc_t read_desc;

		v_reads.clear();
		filled_size = 0;

		while (!in_f->Eof())
		{
			if (!get_read(in_f, input_FASTQ + filled_size, read_desc))
				break;

			filled_size += read_desc.read_size();

			v_reads.push_back(read_desc);

			if (block_size - filled_size < block_size_margin)
				break;
		}

		return true;
	}

	bool Read(IN_FILE_TYPENAME *in_f_1, IN_FILE_TYPENAME *in_f_2)
	{
		read_desc_t read_desc_1;
		read_desc_t read_desc_2;

		v_reads.clear();
		filled_size = 0;

		while (!in_f_1->Eof() && !in_f_2->Eof())
		{
			if (!get_read(in_f_1, input_FASTQ + filled_size, read_desc_1))
				break;
			filled_size += read_desc_1.read_size();
			v_reads.push_back(read_desc_1);

			if (!get_read(in_f_2, input_FASTQ + filled_size, read_desc_2))
				break;
			filled_size += read_desc_2.read_size();
			v_reads.push_back(read_desc_2);
		
			if (block_size - filled_size < 2 * block_size_margin)
				break;
		}

		return true;
	}

	bool Write(OUT_FILE_TYPENAME *out_f)
	{
		if (!out_f)
			return false;

		out_f->Write(output_FASTQ, filled_size);
		filled_size = 0;

		return true;
	}

	bool Write(OUT_FILE_TYPENAME *out_f_1, OUT_FILE_TYPENAME *out_f_2)
	{
		return false;
	}
	
	uint64_t Size()
	{
		return filled_size;
	}

	void PartitionForWorkers(uint64_t no_workers)
	{
		PartitionForWorkers(no_workers, v_reads.size());
	}

	void PartitionForWorkers(uint64_t no_workers, uint64_t no_reads)
	{
		v_ranges.clear();

		uint64_t lower = 0;
		uint64_t upper;

		for (uint64_t i = 0; i < no_workers; ++i)
		{
			upper = (i + 1) * no_reads / no_workers;
			if(i < no_workers - 1)
				upper &= ~1ull;		
			
			v_ranges.push_back(make_pair(lower, upper));

			lower = upper;
		}
	}
};

// EOF
