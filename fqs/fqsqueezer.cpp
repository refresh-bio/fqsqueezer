// *******************************************************************************************
// This file is a part of FQSqueezer software distributed under GNU GPL 3 licence.
// The homepage of the FQSqueezer project is http://sun.aei.polsl.pl/REFRESH/fqsqueezer
//
// Author: Sebastian Deorowicz
// Version: 1.1
// Date   : 2020-06-16
// *******************************************************************************************

#include "defs.h"
#include "application.h"
#include "io.h"
#include "params.h"
#include <iostream>
#include <string>
#include <chrono>

CInFile infile;
COutFile outfiles[257];

CParams params;
CApplication *application;

bool parse_params(int argc, char **argv);
void usage();

using namespace std::chrono;

//*****************************************************************************************************
void usage()
{
	cerr << FQS_VERSION << endl;
	cerr << "Usage for compression:\n";
	cerr << "  fqs-0.1 e [compression-options] <input.fastq>\n";
	cerr << "  fqs-0.1 e [compression-options] <input1.fastq> <input2.fastq>\n";
	cerr << "  fqs-0.1 e [compression-options] @<file_list>\n";
	cerr << "Usage for decompression:\n";
	cerr << "  fqs-0.1 d [decompression-options] <input.fastq>\n";
	cerr << "compression-options:\n";
	cerr << "  -s          - single-end data\n";
	cerr << "  -p          - paired-end data\n";
	cerr << "  -t <num>    - number of threads (default: " << params.no_threads << ")\n";
	cerr << "  -gs <num>   - approx. length of genome in Mbp (default: " << params.genome_size << ")\n";
	cerr << "  -tmp <path> - path to temporary files (default: " << params.tmp_path << ")\n";
	cerr << "  -out <path> - path to output file (default: " << params.out_path << ")\n";
	cerr << "  -om <s|o>   - order of reads: s (sorted), o (original) (default: " <<
		((params.dna_mode == dna_mode_t::pe_sorted || params.dna_mode == dna_mode_t::se_sorted) ? "s" : "o") << ")\n";
	cerr << "  -qm <o|8|4|2|n> - quality mode: o (original), 8 (Illumina 8-lev.), 4 (Illumina 4-lev.), 2 (binary thr), n (none) (default: " <<
		(params.quality_mode == quality_mode_t::lossless ? "o" :
			params.quality_mode == quality_mode_t::illumina_8 ? "8" :
			params.quality_mode == quality_mode_t::illumina_4 ? "4" : 
			params.quality_mode == quality_mode_t::binary ? "2" :
			"none") << ")\n";
	cerr << "  -qt <num>   - threshold of quality for \"trusted\" base (default: " << params.quality_thr << ")\n";
	cerr << "  -im <o|i|n> - id mode: o (original), i (instrument only), n (none) (default: " <<
		(params.id_mode == id_mode_t::lossless ? "l" : 
			params.id_mode == id_mode_t::instrument ? "i" : "n") << ")\n";
#ifdef DEVELOPMENT_MODE
	cerr << "  -d          - enable checking of duplicate reads (default: " << (bool) params.duplicates_check << ")\n";
	cerr << "  --prefix-len <num> - prefix length (default: " << params.prefix_len << ")\n";
	cerr << "  --pmer-len <num> - prefix length (default: " << params.pmer_len << ")\n";
	cerr << "  --smer-len <num> - k-mer length (default: " << params.smer_len << ")\n";
	cerr << "  --bmer-len <num> - short length (default: " << params.bmer_len << ")\n";
	cerr << "  --imer-len <num> - i-mer length (default: " << params.imer_len << ")\n";
	cerr << "  --hmer-len <num> - short length (default: " << params.hmer_len << ")\n";
#endif
	cerr << "  -v <num>    - verbosity from 0 to 2 (default: " << (bool)params.verbosity << ")\n";
	cerr << "decompression-options:\n";
	cerr << "  -out <path> - path to 1st (or only in SE mode) output file (default: " << params.out_path << ")\n";
	cerr << "  -out2 <path> - path to 2nd output file in PE mode (default: " << params.out_path2 << ")\n";
}

//*****************************************************************************************************
bool parse_params(int argc, char **argv)
{
	if (argc < 3)
		return false;

	if (std::string(argv[1]) == "e")
		params.work_mode = work_mode_t::compress;
	else if (std::string(argv[1]) == "d")
		params.work_mode = work_mode_t::decompress;
	else
	{
		cerr << "Invalid mode: " << std::string(argv[1]) << endl;
		return false;
	}

	int i;

	for (i = 2; i < argc; ++i)
	{
		std::string par = std::string(argv[i]);

		if (par[0] != '-')
			break;
		
		if (par == "-s")
		{
			if (params.dna_mode == dna_mode_t::pe_original)
				params.dna_mode = dna_mode_t::se_original;
			else if (params.dna_mode == dna_mode_t::pe_sorted)
				params.dna_mode = dna_mode_t::se_sorted;
		}
		else if (par == "-p")
		{
			if (params.dna_mode == dna_mode_t::se_original)
				params.dna_mode = dna_mode_t::pe_original;
			else if (params.dna_mode == dna_mode_t::se_sorted)
				params.dna_mode = dna_mode_t::pe_sorted;
		}
		else if (par == "-t" && i + 1 < argc)
		{
			params.no_threads = NormalizeValue(atoi(argv[i + 1]), 1, 64);
			++i;
		}
		else if (par == "-gs" && i + 1 < argc)
		{
			params.genome_size = NormalizeValue(atoi(argv[i + 1]), 1, 32768);
			++i;
		}
		else if (par == "-tmp" && i + 1 < argc)
		{
			params.tmp_path = std::string(argv[i + 1]);
			++i;
		}
		else if (par == "-out" && i + 1 < argc)
		{
			params.out_path = std::string(argv[i + 1]);
			++i;
		}
		else if (par == "-out2" && i + 1 < argc)
		{
			params.out_path2 = std::string(argv[i + 1]);
			++i;
		}
		else if (par == "-om" && i + 1 < argc)
		{
			if (std::string(argv[i + 1]) == "o")
			{
				if (params.dna_mode == dna_mode_t::pe_sorted)
					params.dna_mode = dna_mode_t::pe_original;
				else if (params.dna_mode == dna_mode_t::se_sorted)
					params.dna_mode = dna_mode_t::se_original;
			}
			else if (std::string(argv[i + 1]) == "s")
			{
				if (params.dna_mode == dna_mode_t::pe_original)
					params.dna_mode = dna_mode_t::pe_sorted;
				else if (params.dna_mode == dna_mode_t::se_original)
					params.dna_mode = dna_mode_t::se_sorted;
			}
			else
				cerr << "Unknown order option: " << par << " " << std::string(argv[i + 1]) << endl;
			++i;
		}
		else if (par == "-qm" && i + 1 < argc)
		{
			std::string par2 = std::string(argv[i + 1]);
			if (par2 == "o")
				params.quality_mode = quality_mode_t::lossless;
			else if (par2 == "8")
				params.quality_mode = quality_mode_t::illumina_8;
			else if (par2 == "4")
				params.quality_mode = quality_mode_t::illumina_4;
			else if (par2 == "2")
				params.quality_mode = quality_mode_t::binary;
			else if (par2 == "n")
				params.quality_mode = quality_mode_t::none;
			else
				cerr << "Unknown quality mode option: " << par << " " << par2 << endl;
			++i;
		}
		else if (par == "-qt" && i + 1 < argc)
		{
			params.quality_thr = atoi(argv[i + 1]);
			++i;
		}
		else if (par == "-im" && i + 1 < argc)
		{
			std::string par2 = std::string(argv[i + 1]);
			if (par2 == "o")
				params.id_mode = id_mode_t::lossless;
			else if (par2 == "i")
				params.id_mode = id_mode_t::instrument;
			else if (par2 == "n")
				params.id_mode = id_mode_t::none;
			else
				cerr << "Unknown instrument mode option: " << par << " " << par2 << endl;
			++i;
		}
#ifdef DEVELOPMENT_MODE
		else if (par == "-d")
		{
			params.duplicates_check = true;
		}
		else if (par == "--prefix-len" && i + 1 < argc)
		{
			params.prefix_len = NormalizeValue(atoi(argv[i + 1]), 1, 64);
			++i;
		}
		else if (par == "--pmer-len" && i + 1 < argc)
		{
			params.pmer_len = NormalizeValue(atoi(argv[i + 1]), 1, 64);
			++i;
		}
		else if (par == "--smer-len" && i + 1 < argc)
		{
			params.smer_len = NormalizeValue(atoi(argv[i + 1]), 1, 64);
			++i;
		}
		else if (par == "--bmer-len" && i + 1 < argc)
		{
			params.bmer_len = NormalizeValue(atoi(argv[i + 1]), 1, 64);
			++i;
		}
		else if (par == "--imer-len" && i + 1 < argc)
		{
			params.imer_len = NormalizeValue(atoi(argv[i + 1]), 1, 64);
			++i;
		}
		else if (par == "--hmer-len" && i + 1 < argc)
		{
			params.hmer_len = NormalizeValue(atoi(argv[i + 1]), 1, 64);
			++i;
		}
#endif
		else if (par == "-v" && i + 1 < argc)
		{
			params.verbosity = NormalizeValue(atoi(argv[i + 1]), 0, 2);
			++i;
		}
		else
			cerr << "unknown option: " << par << endl;
	}

	if (i >= argc)
	{
		cerr << "No input files given\n";
		return false;
	}

	if (argv[i][0] == '@')
	{
		std::string fn = std::string(argv[i] + 1);
		FILE *f_list = fopen(argv[i] + 1, "rt");
		if (!f_list)
		{
			cerr << "Cannot read file list file\n";
			return false;
		}

		char s[1024];
		while (true)
		{
			if (fgets(s, 1024, f_list) == nullptr)
				break;
			if (feof(f_list))
				break;

			params.v_file_names.push_back(trim(std::string(s)));
		}
		fclose(f_list);

		if (params.v_file_names.empty())
		{
			cerr << "Empty file list\n";
			return false;
		}
	}
	else if (params.dna_mode == dna_mode_t::se_original || params.dna_mode == dna_mode_t::se_sorted)
		params.v_file_names.push_back(std::string(argv[i]));
	else if (i + 1 < argc)
	{
		params.v_file_names.push_back(std::string(argv[i]));
		params.v_file_names.push_back(std::string(argv[i + 1]));
	}
	else
	{
		cerr << "Two few file names in paired-end mode\n";
		return false;
	}

	params.adjust_kmer_sizes();

	return true;
}

//*****************************************************************************************************
int main(int argc, char **argv)
{
	if (!parse_params(argc, argv))
	{
		usage();
		return 0;
	}

	auto t1 = high_resolution_clock::now();

#ifdef _WIN32
	_setmaxstdio(2 * (1u << (2 * BIN_PREFIX_LEN)) + 5);
#else
#endif

	application = new CApplication();
	application->SetParams(params);
	application->Execute();
	delete application;

	auto t2 = high_resolution_clock::now();

	std::chrono::duration<double> diff = t2 - t1;
	std::cerr << "Processing time: " << diff.count() << " s\n";

	return 0;
}

// EOF
