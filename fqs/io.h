#pragma once
// *******************************************************************************************
// This file is a part of FQSqueezer software distributed under GNU GPL 3 licence.
// The homepage of the FQSqueezer project is http://sun.aei.polsl.pl/REFRESH/fqsqueezer
//
// Author: Sebastian Deorowicz
// Version: 1.1
// Date   : 2020-06-16
// *******************************************************************************************

#include <algorithm>
#include <vector>
#include <string>
#include <cstdint>
#include <cstring>

#include <iostream>

using namespace std;

#ifndef _WIN32
#define my_fseek	fseek
#define my_ftell	ftell
#else
#define my_fseek	_fseeki64
#define my_ftell	_ftelli64
#endif

// *******************************************************************************************
// Buffered input file
class CInFile
{
	const size_t BUFFER_SIZE = 128 << 20;

	FILE *f;
	uint8_t *buffer;
	size_t buffer_pos;
	size_t buffer_filled;

	size_t file_size;
	size_t before_buffer_bytes;

public:
	CInFile() : f(nullptr), buffer(nullptr), buffer_pos(0), buffer_filled(0), file_size(0), before_buffer_bytes(0)
	{};

	~CInFile()
	{
		if (f)
			fclose(f);
		if (buffer)
			delete[] buffer;
	}

	bool Open(std::string file_name, uint32_t verbosity)
	{
		if (f)
			return false;

		f = fopen(file_name.c_str(), "rb");
		if (!f)
			return false;

		my_fseek(f, 0, SEEK_END);
		file_size = my_ftell(f);
		my_fseek(f, 0, SEEK_SET);
		before_buffer_bytes = 0;

		buffer = new uint8_t[BUFFER_SIZE];
		buffer_pos = 0;
		buffer_filled = 0;

		if (verbosity == 2)
		{
			cout << std::string("Opening file of size: ") + to_string(file_size) + "\n";
			fflush(stdout);
		}

		return true;
	}

	bool Close()
	{
		if (f)
		{
			fclose(f);
			f = nullptr;
		}
		if (buffer)
		{
			delete[] buffer;
			buffer = nullptr;
		}

		return true;
	}

	int Get()
	{
		if (buffer_pos < buffer_filled)
			return buffer[buffer_pos++];

		if (feof(f))
			return EOF;

		before_buffer_bytes += buffer_filled;

		buffer_filled = fread(buffer, 1, BUFFER_SIZE, f);
		if (buffer_filled == 0)
			return EOF;

		buffer_pos = 0;
		return buffer[buffer_pos++];
	}

	uint64_t ReadUInt(const int no_bytes)
	{
		uint64_t x = 0;
		uint64_t shift = 0;

		for (int i = 0; i < no_bytes; ++i)
		{
			uint64_t c = Get();
			x += c << shift;
			shift += 8;
		}

		return x;
	}

	uint64_t ReadUIntVar()
	{
		uint64_t x = 0;
		uint64_t c = Get();

		if ((c >> 7) == 0)		// [0, 0x8000)
		{
			x = c << 8;
			x += (uint64_t)Get();
		}
		else if ((c >> 6) == 0b10)	// [0x8000, 0x400000)
		{
			x = (c & 0x3f) << 16;
			x += (uint64_t) Get() << 8;
			x += (uint64_t)Get();
		}
		else if ((c >> 6) == 0b11)	// [0x80000, 0x4000000)
		{
			x = (c & 0x3f) << 24;
			x += (uint64_t)Get() << 16;
			x += (uint64_t)Get() << 8;
			x += (uint64_t)Get();
		}

		return x;
	}
		
	void Read(uint8_t *ptr, uint64_t size)
	{
		if (before_buffer_bytes + buffer_pos + size > file_size)
			size = file_size - (before_buffer_bytes + buffer_pos);

		uint64_t to_read = size;

		while (buffer_pos + to_read > BUFFER_SIZE)
		{
			memcpy(ptr, buffer + buffer_pos, BUFFER_SIZE - buffer_pos);
			ptr += BUFFER_SIZE - buffer_pos;
			to_read -= BUFFER_SIZE - buffer_pos;

			before_buffer_bytes += buffer_filled;
			buffer_filled = fread(buffer, 1, BUFFER_SIZE, f);
			buffer_pos = 0;
		}

		memcpy(ptr, buffer + buffer_pos, to_read);
		buffer_pos += to_read;
	}

	bool Eof() const
	{
		return before_buffer_bytes + buffer_pos >= file_size;
	}

	size_t FileSize() const
	{
		if (f)
			return file_size;
		else
			return 0;
	}

	size_t GetPos() const
	{
		return before_buffer_bytes + buffer_pos;
	}
};

// *******************************************************************************************
// Buffered output file
class COutFile
{
	const size_t BUFFER_SIZE = 8 << 20;

	FILE *f;
	uint8_t *buffer;
	size_t buffer_pos;
	bool success;

public:
	COutFile() : f(nullptr), buffer(nullptr), buffer_pos(0), success(false)
	{};

	~COutFile()
	{
		if (f)
			Close();
		if (buffer)
			delete[] buffer;
	}

	bool Open(const std::string &file_name)
	{
		if (f)
			return false;

		f = fopen(file_name.c_str(), "wb");
		if (!f)
			return false;

		buffer = new uint8_t[BUFFER_SIZE];
		buffer_pos = 0;
		success = true;

		return true;
	}

	bool Close()
	{
		if (buffer_pos)
			success &= fwrite(buffer, 1, buffer_pos, f) == buffer_pos;

		if (f)
		{
			fclose(f);
			f = nullptr;
		}
		if (buffer)
		{
			delete[] buffer;
			buffer = nullptr;
		}

		return success;
	}

	void Put(const uint8_t c)
	{
		if (buffer_pos == BUFFER_SIZE)
		{
			success &= fwrite(buffer, 1, BUFFER_SIZE, f) == BUFFER_SIZE;
			buffer_pos = 0;
		}

		buffer[buffer_pos++] = c;
	}

	void Write(const uint8_t *p, size_t n)
	{
		uint8_t *q = (uint8_t *)p;

		while (buffer_pos + n > BUFFER_SIZE)
		{
			size_t small_n = BUFFER_SIZE - buffer_pos;
			memcpy(buffer + buffer_pos, q, small_n);
			success &= fwrite(buffer, 1, BUFFER_SIZE, f) == BUFFER_SIZE;

			buffer_pos = 0;
			n -= small_n;
			q += small_n;
		}

		memcpy(buffer+buffer_pos, q, n);
		buffer_pos += n;
	}

	void WriteUInt(uint64_t x, int no_bytes)
	{
		for (int i = 0; i < no_bytes; ++i)
		{
			Put(x & 0xff);
			x >>= 8;
		}
	}

	// [0, 0x8000) -> 0 [wart. 15-bit]
	// [0x8000, 0x400000) -> 10 [wart. 22-bit]
	// [0x400000, 0x40000000) -> 11 [wart. 30-bit]
	// wieksze wart. nie sa wspierane
	void WriteUIntVar(uint64_t x)
	{
		if (x < 0x8000ull)		
		{
			Put((uint8_t) (x >> 8));
			Put((uint8_t) (x & 0xff));
		}
		else if (x < 0x400000ull)
		{
			Put((uint8_t) (0x80 + (x >> 16)));
			Put((uint8_t) ((x >> 8) & 0xff));
			Put((uint8_t) (x & 0xff));
		}
		else if (x < 0x40000000ull)
		{
			Put((uint8_t) (0xc0 + (x >> 24)));
			Put((uint8_t) ((x >> 16) & 0xff));
			Put((uint8_t) ((x >> 8) & 0xff));
			Put((uint8_t) (x & 0xff));
		}
		else
			cerr << "To large value\n";
	}

	void Write(const std::string &s)
	{
		Write((uint8_t*)s.c_str(), s.size());
	}

	void Write(const std::string &s, const size_t start_pos, const size_t len)
	{
		Write((uint8_t*)s.c_str() + start_pos, len);
	}
};

// *******************************************************************************************
// Class for storage of range coder compressed data
class CVectorIOStream
{
	vector<uint8_t> &v;
	size_t read_pos;

public:
	CVectorIOStream(vector<uint8_t> &_v) : v(_v), read_pos(0)
	{}

	void RestartRead()
	{
		read_pos = 0;
	}

	bool Eof() const
	{
		return read_pos >= v.size();
	}

	uint8_t GetByte() 
	{
		return v[read_pos++];
	}

	void PutByte(uint8_t x)
	{
		v.push_back(x);
	}

	size_t Size()
	{
		return v.size();
	}
};

// *******************************************************************************************
class CBasicFASTQFile
{
public:
	CBasicFASTQFile()
	{};

	virtual ~CBasicFASTQFile() {}
	virtual bool Open(const std::string &file_name, uint32_t verbosity) = 0;
	virtual bool Close() = 0;
	virtual int Get() = 0;
	virtual bool Eof() = 0;
	virtual size_t FileSize() = 0;
	virtual size_t GetPos() = 0;
};

// *******************************************************************************************
class CPlainFASTQFile : public CBasicFASTQFile
{
	CInFile in_file;

public:
	CPlainFASTQFile()
	{};

	virtual ~CPlainFASTQFile()
	{}

	virtual bool Open(const std::string &file_name, uint32_t verbosity)
	{
		return in_file.Open(file_name, verbosity);
	}

	virtual bool Close()
	{
		return in_file.Close();
	}

	virtual int Get()
	{
		return in_file.Get();
	}

	virtual bool Eof()
	{
		return in_file.Eof();
	}

	virtual size_t FileSize()
	{
		return in_file.FileSize();
	}

	virtual size_t GetPos()
	{
		return in_file.GetPos();
	}
};

// *******************************************************************************************
// Sorted input FASTQ file
class CSortedFASTQFile : public CBasicFASTQFile
{
	char dna_convert_NT[256];

	FILE *f;
	std::string fn;
	CSortedFASTQFile *f_sorting_order;
	uint8_t *buffer;
	size_t buffer_pos;
	size_t buffer_filled;
	size_t file_size;

	vector<pair<read_desc_t, uint32_t>> v_reads;
	size_t i_read;
	size_t n_reads;
	size_t cur_read_pos;
	size_t cur_read_size;

	bool get_read(read_desc_t &read_desc)
	{
		int no_eols = 0;
		int read_len = 0;

		uint8_t *p = buffer + buffer_pos;

		read_desc.id = p;

		while (buffer_pos < buffer_filled && no_eols < 4)
		{
			int c = buffer[buffer_pos];
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

			++buffer_pos;
			++p;
		}

		return no_eols == 4;
	}

	// Find read
	void find_reads()
	{
		read_desc_t read_desc;
		uint32_t read_id = 0;

		v_reads.clear();

		for (; get_read(read_desc); ++read_id)
			v_reads.push_back(make_pair(read_desc, read_id));

		v_reads.shrink_to_fit();
	}

	// Sort reads
	void sort_reads()
	{
		sort(v_reads.begin(), v_reads.end(), [&](pair<read_desc_t, uint32_t> &x, pair<read_desc_t, uint32_t> &y) {
			auto x_len = x.first.read_len();
			auto y_len = y.first.read_len();
			auto min_len = min(x_len, y_len);

			for (uint32_t i = 0; i < min_len; ++i)
			{
				if (dna_convert_NT[x.first.dna[i]] == dna_convert_NT[y.first.dna[i]])
					continue;

				return dna_convert_NT[x.first.dna[i]] < dna_convert_NT[y.first.dna[i]];
			}

			if (x_len != y_len)
				return x_len < y_len;

			// If the reads seem to be equal in case of ignoring Ns then try lexicographical comparison
			for (uint32_t i = 0; i < min_len; ++i)
			{
				if (x.first.dna[i] == y.first.dna[i])
					continue;

				return x.first.dna[i] < y.first.dna[i];
			}

			return false;
		});
	}

	//Return ordering of the sorted collection of reads
	void get_ordering(vector<uint32_t> &v_ordering)
	{
		v_ordering.clear();
		v_ordering.resize(v_reads.size());

		for (size_t i = 0; i < v_reads.size(); ++i)
			v_ordering[i] = v_reads[i].second;
	}
	
	// Adjust sorting order according to other file
	void adjust_ordering()
	{
		vector<uint32_t> v_ordering;
		f_sorting_order->get_ordering(v_ordering);

		vector<pair<read_desc_t, uint32_t>> v_old_reads(v_reads);

		for (size_t i = 0; i < v_reads.size(); ++i)
			v_reads[i] = v_old_reads[v_ordering[i]];
	}

	void prepare_dna_convert()
	{
		// N is treated as T. Just for distribution into bins
		for (int i = 0; i < 256; ++i)
			dna_convert_NT[i] = 3;
		dna_convert_NT['A'] = 0;
		dna_convert_NT['C'] = 1;
		dna_convert_NT['G'] = 2;
		dna_convert_NT['T'] = 3;
	}

public:
	CSortedFASTQFile(CSortedFASTQFile *_f_sorting_order = nullptr) : f(nullptr), f_sorting_order(_f_sorting_order), buffer(nullptr), i_read(0)
	{
		prepare_dna_convert();
	}
		
	virtual ~CSortedFASTQFile()
	{
		if (f)
		{
			fclose(f);
			remove(fn.c_str());
		}
		if (buffer)
			delete[] buffer;
	}

	virtual bool Open(const std::string &file_name, uint32_t verbosity)
	{
		if (f)
			return false;

		f = fopen(file_name.c_str(), "rb");
		if (!f)
			return false;

		fn = file_name;

		my_fseek(f, 0, SEEK_END);
		file_size = my_ftell(f);
		my_fseek(f, 0, SEEK_SET);

		buffer = new uint8_t[file_size];

		buffer_pos = 0;

		if (verbosity == 2)
		{
			cout << std::string("Opening file of size: ") + to_string(file_size) + "\n";
			fflush(stdout);
		}

		buffer_filled = fread(buffer, 1, file_size, f);
		
		find_reads();
		if (!f_sorting_order)
			sort_reads();
		else
			adjust_ordering();

		n_reads = v_reads.size();
		i_read = 0;
		cur_read_pos = 0;
		
		if (v_reads.empty())
			cur_read_size = 0;
		else
			cur_read_size = v_reads[i_read].first.read_size();

		return true;
	}

	virtual bool Close()
	{
		if (f)
		{
			fclose(f);
			remove(fn.c_str());
			f = nullptr;
		}
		if (buffer)
		{
			delete[] buffer;
			buffer = nullptr;
		}

		return true;
	}

	virtual int Get()
	{
		if (i_read >= n_reads)
			return EOF;

		// .id points to the beginning of the read, so it is used here in this way
		int r = v_reads[i_read].first.id[cur_read_pos];

		if (++cur_read_pos == cur_read_size)
		{
			++i_read;
			cur_read_pos = 0;
			if(i_read < n_reads)
				cur_read_size = v_reads[i_read].first.read_size();
		}

		return r;
	}

	virtual bool Eof()
	{
		return i_read >= n_reads;
	}

	virtual size_t FileSize()
	{
		if (f)
			return file_size;
		else
			return 0;
	}

	virtual size_t GetPos()
	{
		return buffer_pos;
	}
};

// EOF
