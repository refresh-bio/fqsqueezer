// *******************************************************************************************
// This file is a part of FQSqueezer software distributed under GNU GPL 3 licence.
// The homepage of the MSAC project is http://sun.aei.polsl.pl/REFRESH/fqsqueezer
//
// Author: Sebastian Deorowicz
// Version: 1.0
// Date   : 2019-02-22
// *******************************************************************************************

#include "id.h"
#include "utils.h"

//*****************************************************************************************************
CIdCompressor::CIdCompressor() : CBasicCompressor()
{
	rc_id_mtf_flag = nullptr;
	for(int i = 0; i < 7; ++i)
		rc_id_mtf_code[i] = nullptr;

	for (int i = 0; i < 4; ++i)
		rc_id_mtf_byte[i] = nullptr;

	init_symbol_classes();
}

//*****************************************************************************************************
CIdCompressor::~CIdCompressor()
{
	if (rc_id_mtf_flag)
		delete rc_id_mtf_flag;
		
	for (int i = 0; i < 7; ++i)
		if(rc_id_mtf_code[i])
			delete rc_id_mtf_code[i];

	for (int i = 0; i < 4; ++i)
		if (rc_id_mtf_byte[i])
			delete rc_id_mtf_byte[i];

	delete tpl_ctx_rc_flags;
	delete tpl_ctx_rc_pe_flags;

	delete tpl_ctx_rc_numeric_size;

	delete tpl_ctx_rc_numeric;
	delete tpl_ctx_rc_numeric_small;
	delete tpl_ctx_rc_literal;

	delete tpl_ctx_rc_literal_same;
	delete tpl_ctx_rc_literal_same_length;

	delete tpl_ctx_rc_plain;
}

//*****************************************************************************************************
void CIdCompressor::init_symbol_classes()
{
	fill_n(a_literal.begin(), 256, false);
	fill_n(a_numeric.begin(), 256, false);

	for (int i = '0'; i <= '9'; ++i)
		a_numeric[i] = a_literal[i] = true;

	for (int i = 'A'; i <= 'Z'; ++i)
		a_literal[i] = true;
	for (int i = 'a'; i <= 'z'; ++i)
		a_literal[i] = true;
	a_literal['@'] = true;
}

//*****************************************************************************************************
bool CIdCompressor::Init()
{
	if (rc_id_mtf_flag)
		delete rc_id_mtf_flag;

	for (int i = 0; i < 7; ++i)
		if (rc_id_mtf_code[i])
			delete rc_id_mtf_code[i];

	for (int i = 0; i < 4; ++i)
		if (rc_id_mtf_byte[i])
			delete rc_id_mtf_byte[i];

	rc_id_mtf_flag = new CRangeCoderModel<CVectorIOStream>(rc, 11, 15, 1 << 15, nullptr, 1, params.work_mode == work_mode_t::compress);

	for(int i = 0; i < 7; ++i)
		rc_id_mtf_code[i] = new CRangeCoderModel<CVectorIOStream>(rc, 2 << i, 15, 1 << 15, nullptr, 1, params.work_mode == work_mode_t::compress);

	for (int i = 0; i < 4; ++i)
		rc_id_mtf_byte[i] = new CRangeCoderModel<CVectorIOStream>(rc, 256, 15, 1 << 15, nullptr, 1, params.work_mode == work_mode_t::compress);

	tpl_ctx_rc_flags = new CRangeCoderModelFixedSize<CVectorIOStream, 2>(rc, 15, 1 << 15, nullptr, 1, params.work_mode == work_mode_t::compress);
	tpl_ctx_rc_pe_flags = new CRangeCoderModelFixedSize<CVectorIOStream, 2>(rc, 15, 1 << 15, nullptr, 1, params.work_mode == work_mode_t::compress);

	tpl_ctx_rc_numeric_size = new CRangeCoderModelFixedSize<CVectorIOStream, 256>(rc, 15, 1 << 15, nullptr, 1, params.work_mode == work_mode_t::compress);
	
	tpl_ctx_rc_numeric = new CRangeCoderModel<CVectorIOStream>(rc, 256, 15, 1 << 15, nullptr, 1, params.work_mode == work_mode_t::compress);
	tpl_ctx_rc_numeric_small = new CRangeCoderModelFixedSize<CVectorIOStream, 4>(rc, 15, 1 << 15, nullptr, 1, params.work_mode == work_mode_t::compress);
	tpl_ctx_rc_literal = new CRangeCoderModel<CVectorIOStream>(rc, 128, 15, 1 << 15, nullptr, 1, params.work_mode == work_mode_t::compress);

	tpl_ctx_rc_literal_same = new CRangeCoderModelFixedSize<CVectorIOStream, 2>(rc, 15, 1 << 15, nullptr, 1, params.work_mode == work_mode_t::compress);
	tpl_ctx_rc_literal_same_length = new CRangeCoderModelFixedSize<CVectorIOStream, 2>(rc, 15, 1 << 15, nullptr, 1, params.work_mode == work_mode_t::compress);

	tpl_ctx_rc_plain = new CRangeCoderModelFixedSize<CVectorIOStream, 128>(rc, 15, 1 << 15, nullptr, 1, params.work_mode == work_mode_t::compress);

	mtf.Reset();

	ResetReadPrev();

	id_prev = nullptr;
	size_prev = 0;
	v_tokens_cur.clear();
	v_tokens_prev.clear();
	v_deltas.clear();

	ctx_flags = 0;
	ctx_pe_flags = 0;

	return true;
}

//*****************************************************************************************************
void CIdCompressor::ResetReadPrev()
{
	id_prev = nullptr;
	size_prev = 0;
	v_tokens_cur.clear();
	v_tokens_prev.clear();
	v_deltas.clear();

	ctx_flags = 0;
	ctx_pe_flags = 0;
}

//*****************************************************************************************************
bool CIdCompressor::Compress(uint8_t *p, uint32_t size)
{
	if (params.id_mode == id_mode_t::none)
		return true;
	else if (params.id_mode == id_mode_t::lossless)
		compress_lossless(p, size);
	else
		compress_instrument(p, size);

	return true;
}

//*****************************************************************************************************
bool CIdCompressor::CompressPE(uint8_t *p1, uint32_t size1, uint8_t *p2, uint32_t size2)
{
	if (params.id_mode == id_mode_t::none)
		return true;
	else if (params.id_mode == id_mode_t::lossless)
	{
		auto rc_pe = find_rc_context(m_ctx_rc_pe_flags, ctx_pe_flags, tpl_ctx_rc_pe_flags);

		if (typical_pe_ids(p1, size1, p2, size2))
		{
			ctx_pe_flags = ((ctx_pe_flags << 1) + 1) & 0xff;
			rc_pe->Encode(1);
			compress_lossless(p1, size1);
		}
		else
		{
			ctx_pe_flags = ((ctx_pe_flags << 1) + 0) & 0xff;
			rc_pe->Encode(0);
			compress_lossless(p1, size1);
			compress_lossless(p2, size2);
		}
	}
	else
	{
		compress_instrument(p1, size1);
		compress_instrument(p2, size2);
	}

	return true;
}

//*****************************************************************************************************
bool CIdCompressor::Decompress(uint8_t *p, uint32_t &size)
{
	if (params.id_mode == id_mode_t::none)
		decompress_none(p, size);
	else if (params.id_mode == id_mode_t::lossless)
		decompress_lossless(p, size);
	else
		decompress_instrument(p, size);

	return size > 1;
}

//*****************************************************************************************************
bool CIdCompressor::DecompressPE(uint8_t *p, uint32_t &size1, uint32_t &size2)
{
	if (params.id_mode == id_mode_t::none)
	{
		decompress_none(p, size1);
		decompress_none(p+size1, size2);
	}
	else if (params.id_mode == id_mode_t::lossless)
	{
		auto rc_pe = find_rc_context(m_ctx_rc_pe_flags, ctx_pe_flags, tpl_ctx_rc_pe_flags);

		if (rc_pe->Decode() == 1)			// typical PE ids
		{
			ctx_pe_flags = ((ctx_pe_flags << 1) + 1) & 0xff;
			decompress_lossless(p, size1);
			copy_n(p, size1, p + size1);
			size2 = size1;
			p[size1 + size2 - 2] = '2';
		}
		else
		{
			ctx_pe_flags = ((ctx_pe_flags << 1) + 0) & 0xff;
			decompress_lossless(p, size1);
			decompress_lossless(p+size1, size2);
		}
	}
	else
	{
		decompress_instrument(p, size1);
		decompress_instrument(p+size1, size2);
	}

	return size1 > 1 && size2 > 1;
}

//*****************************************************************************************************
uint32_t CIdCompressor::extract_instrument(const uint8_t *p, uint32_t size)
{
	for (uint32_t i = 0; i < size; ++i)
		if (p[i] == '.' || p[i] == ' ' || p[i] == ':')
			return i;
	
	return size;
}

//*****************************************************************************************************
// Check whether ids od PE reads are typical, i.e., differ only on the last symbol (1 in 1st read and 2 in 2nd read)
bool CIdCompressor::typical_pe_ids(const uint8_t *p1, const uint32_t size1, const uint8_t *p2, const uint32_t size2) const
{
	if (size1 != size2 || size1 < 3)
		return false;

	if (!equal(p1, p1 + size1 - 2, p2))
		return false;

	if (p1[size1 - 2] != '1' || p2[size2 - 2] != '2')
		return false;

	return true;
}

//*****************************************************************************************************
void CIdCompressor::compress_lossless(uint8_t *p, uint32_t size)
{
	auto n_tokens = tokenize(p, size, v_tokens_cur);

	auto rc_flags = find_rc_context(m_ctx_rc_flags, ctx_flags, tpl_ctx_rc_flags);
	if (token_types_same(v_tokens_cur, v_tokens_prev))
	{
		rc_flags->Encode(1);
		ctx_flags = ((ctx_flags << 1) + 1) & 0xff;

		for (uint32_t i = 0; i < n_tokens; ++i)
		{
			if (get<0>(v_tokens_cur[i]) == token_type_t::literal)
			{
				bool same_length = (get<3>(v_tokens_cur[i]) - get<2>(v_tokens_cur[i])) == (get<3>(v_tokens_prev[i]) - get<2>(v_tokens_prev[i]));
				bool same = false;
				if (same_length)
					same = equal(p + get<2>(v_tokens_cur[i]), p + get<3>(v_tokens_cur[i]), id_prev + get<2>(v_tokens_prev[i]));

				auto rc_same = find_rc_context(m_ctx_rc_literal_same, i, tpl_ctx_rc_literal_same);
				
				if (same)
					rc_same->Encode(1);
				else
				{
					rc_same->Encode(0);
					auto rc_same_length = find_rc_context(m_ctx_rc_literal_same_length, i, tpl_ctx_rc_literal_same_length);
					if (same_length)
					{
						rc_same_length->Encode(1);
						for (uint32_t j = 0; j < get<3>(v_tokens_cur[i]) - get<2>(v_tokens_cur[i]); ++j)
						{
							auto rc = find_rc_context(m_ctx_rc_literal, ctx_flags + (1ll << 32) + j, tpl_ctx_rc_literal);
							if (p[j + get<2>(v_tokens_cur[i])] == id_prev[j + get<2>(v_tokens_prev[i])])
								rc->Encode(0);
							else
								rc->Encode(p[j + get<2>(v_tokens_cur[i])]);
						}
					}
					else
					{
						rc_same_length->Encode(0);
						for (uint32_t j = get<2>(v_tokens_cur[i]); j < get<3>(v_tokens_cur[i]); ++j)
						{
							auto rc = find_rc_context(m_ctx_rc_literal, ctx_flags + j - get<2>(v_tokens_cur[i]), tpl_ctx_rc_literal);
							rc->Encode(p[j]);
						}

						auto rc = find_rc_context(m_ctx_rc_literal, ctx_flags + (get<3>(v_tokens_cur[i]) - get<2>(v_tokens_cur[i])), tpl_ctx_rc_literal);
						rc->Encode(0);
					}
				}
			}
			else
			{
				int64_t v_prev = get_int(id_prev, get<2>(v_tokens_prev[i]), get<3>(v_tokens_prev[i]));
				int64_t v_cur = get_int(p, get<2>(v_tokens_cur[i]), get<3>(v_tokens_cur[i]));
				int64_t delta = v_cur - v_prev;

				context_t ctx_type = (uint64_t)i << 40;
				ctx_type += (uint64_t)ilog2(abs(v_deltas[i])) << 31;
				ctx_type += (uint64_t) (v_deltas[i] < 0) << 30;
				auto rc_numeric_size = find_rc_context(m_ctx_rc_numeric_size, ctx_type, tpl_ctx_rc_numeric_size);
				auto rc_numeric_small = find_rc_context(m_ctx_rc_numeric_small, ctx_type, tpl_ctx_rc_numeric_small);

				v_deltas[i] = delta;

				if (delta >= -1 && delta <= 1)
				{
					rc_numeric_small->Encode((int) (delta + 1));
				}
				else
				{
					rc_numeric_small->Encode(3);
					int n_bytes = 0;
					if (delta >= -123 && delta <= 123)
						rc_numeric_size->Encode(((uint32_t)(delta + 123)) & 0xff);
					else if (delta > 0 && delta < 0x10000ll)
					{
						rc_numeric_size->Encode(247);
						n_bytes = 2;
						ctx_type += 0x10;
					}
					else if (delta > 0 && delta < 0x1000000ll)
					{
						rc_numeric_size->Encode(248);
						n_bytes = 3;
						ctx_type += 0x20;
					}
					else if (delta > 0 && delta < 0x100000000)
					{
						rc_numeric_size->Encode(249);
						n_bytes = 4;
						ctx_type += 0x30;
					}
					else if (delta > 0)
					{
						rc_numeric_size->Encode(250);
						n_bytes = 8;
						ctx_type += 0x40;
					}
					else if (delta < 0 && delta > -0x10000ll)
					{
						rc_numeric_size->Encode(251);
						delta = -delta;
						n_bytes = 2;
						ctx_type += 0x50;
					}
					else if (delta < 0 && delta > -0x1000000ll)
					{
						rc_numeric_size->Encode(252);
						delta = -delta;
						n_bytes = 3;
						ctx_type += 0x60;
					}
					else if (delta < 0 && delta > -0x100000000ll)
					{
						rc_numeric_size->Encode(253);
						delta = -delta;
						n_bytes = 4;
						ctx_type += 0x70;
					}
					else if (delta < 0)
					{
						rc_numeric_size->Encode(254);
						delta = -delta;
						n_bytes = 8;
						ctx_type += 0x80;
					}

					for (int j = 0; j < n_bytes; ++j)
					{
						auto rc = find_rc_context(m_ctx_rc_numeric_size, ctx_type + j, tpl_ctx_rc_numeric_size);
						rc->Encode(((uint64_t) delta >> (8 * j)) & 0xff);
					}
				}
			}
		}
	}
	else
	{
		rc_flags->Encode(0);
		ctx_flags = ((ctx_flags << 1) + 0) & 0xff;
		
		// Encode id plain
		for (uint32_t i = 0; i < size; ++i)
		{
			auto rc = find_rc_context(m_ctx_rc_plain, i, tpl_ctx_rc_plain);
			rc->Encode(p[i]);
		}

		v_deltas.clear();
		v_deltas.resize(v_tokens_cur.size(), 0u);
	}
	
	v_tokens_prev.swap(v_tokens_cur);
	id_prev = p;
	size_prev = size;
}

//*****************************************************************************************************
void CIdCompressor::compress_instrument(uint8_t *p, uint32_t size)
{
	uint32_t instr_name_size = extract_instrument(p, size);
	p[instr_name_size] = 0;

	string instr_name((char *)p);

	int mtf_code = mtf.GetCode(instr_name);

	if (mtf_code < 0)
	{
		rc_id_mtf_flag->Encode(0);					// new instrument name
		compress_lossless(p, instr_name_size+1);
	}
	else if (mtf_code < 2)
	{
		rc_id_mtf_flag->Encode(mtf_code + 1);
	}
	else if (mtf_code < 4)
	{
		rc_id_mtf_flag->Encode(3);
		rc_id_mtf_code[0]->Encode(mtf_code - 2);
	}
	else if (mtf_code < 8)
	{
		rc_id_mtf_flag->Encode(4);
		rc_id_mtf_code[1]->Encode(mtf_code - 4);
	}
	else if (mtf_code < 16)
	{
		rc_id_mtf_flag->Encode(5);
		rc_id_mtf_code[2]->Encode(mtf_code - 8);
	}
	else if (mtf_code < 32)
	{
		rc_id_mtf_flag->Encode(6);
		rc_id_mtf_code[3]->Encode(mtf_code - 16);
	}
	else if (mtf_code < 64)
	{
		rc_id_mtf_flag->Encode(7);
		rc_id_mtf_code[4]->Encode(mtf_code - 32);
	}
	else if (mtf_code < 128)
	{
		rc_id_mtf_flag->Encode(8);
		rc_id_mtf_code[5]->Encode(mtf_code - 64);
	}
	else if (mtf_code < 256)
	{
		rc_id_mtf_flag->Encode(9);
		rc_id_mtf_code[6]->Encode(mtf_code - 128);
	}
	else
	{
		rc_id_mtf_flag->Encode(10);					// large MTF code
		
		for (int i = 0; i < 4; ++i)
		{
			rc_id_mtf_byte[i]->Encode(mtf_code & 0xff);
			mtf_code >>= 8;
		}
	}

	mtf.Insert(instr_name);
}

//*****************************************************************************************************
void CIdCompressor::decompress_none(uint8_t *p, uint32_t &size)
{
	size = 0;

	p[size++] = '@';
	p[size++] = 0xA;
}

//*****************************************************************************************************
void CIdCompressor::decompress_lossless(uint8_t *p, uint32_t &size)
{
	size = 0;

	auto rc_flags = find_rc_context(m_ctx_rc_flags, ctx_flags, tpl_ctx_rc_flags);
	
	if (rc_flags->Decode() == 1)	// tokens of the same type
	{
		ctx_flags = ((ctx_flags << 1) + 1) & 0xff;
		auto n_tokens = v_tokens_prev.size();

		for (uint32_t i = 0; i < n_tokens; ++i)
		{
			if (get<0>(v_tokens_prev[i]) == token_type_t::literal)
			{
				auto rc_same = find_rc_context(m_ctx_rc_literal_same, i, tpl_ctx_rc_literal_same);

				if (rc_same->Decode() == 1)		// same token
				{
					copy_n(id_prev + get<2>(v_tokens_prev[i]), get<3>(v_tokens_prev[i]) - get<2>(v_tokens_prev[i]), p+size);
					size += get<3>(v_tokens_prev[i]) - get<2>(v_tokens_prev[i]);
				}
				else
				{
					auto rc_same_length = find_rc_context(m_ctx_rc_literal_same_length, i, tpl_ctx_rc_literal_same_length);
					
					if (rc_same_length->Decode() == 1)	// same lengths
					{
						for (uint32_t j = 0; j < get<3>(v_tokens_prev[i]) - get<2>(v_tokens_prev[i]); ++j)
						{
							// !!! Fix mixing literal contexts
							auto rc = find_rc_context(m_ctx_rc_literal, ctx_flags + (1ll << 32) + j, tpl_ctx_rc_literal);
							auto d = rc->Decode();

							if (d == 0)
								p[size++] = id_prev[j + get<2>(v_tokens_prev[i])];
							else
								p[size++] = d;
						}
					}
					else
					{
						for(uint32_t j = 0; ; ++j)
						{
							auto rc = find_rc_context(m_ctx_rc_literal, ctx_flags + j, tpl_ctx_rc_literal);
							auto d = rc->Decode();
							if (d)
								p[size++] = d;
							else
								break;
						}
					}
				}
			}
			else
			{
				int64_t v_prev = get_int(id_prev, get<2>(v_tokens_prev[i]), get<3>(v_tokens_prev[i]));
				int64_t v_cur;
				int64_t delta;

				context_t ctx_type = (uint64_t)i << 40;
				ctx_type += (uint64_t)ilog2(abs(v_deltas[i])) << 31;
				ctx_type += (uint64_t)(v_deltas[i] < 0) << 30;
				auto rc_numeric_size = find_rc_context(m_ctx_rc_numeric_size, ctx_type, tpl_ctx_rc_numeric_size);
				auto rc_numeric_small = find_rc_context(m_ctx_rc_numeric_small, ctx_type, tpl_ctx_rc_numeric_small);
				
				int64_t d = rc_numeric_small->Decode();

				if (d < 3)
				{
					delta = d - 1ll;
				}
				else
				{
					int n_bytes = 0;
					delta = 0;

					d = rc_numeric_size->Decode();
					if (d <= 246)
					{
						delta = d - 123ll;
					}
					else if (d == 247)
					{
						n_bytes = 2;
						ctx_type += 0x10;
					}
					else if (d == 248)
					{
						n_bytes = 3;
						ctx_type += 0x20;
					}
					else if (d == 249)
					{
						n_bytes = 4;
						ctx_type += 0x30;
					}
					else if (d == 250)
					{
						n_bytes = 8;
						ctx_type += 0x40;
					}
					else if (d == 251)
					{
						n_bytes = 2;
						ctx_type += 0x50;
					}
					else if (d == 252)
					{
						n_bytes = 3;
						ctx_type += 0x60;
					}
					else if (d == 253)
					{
						n_bytes = 4;
						ctx_type += 0x70;
					}
					else if (d == 254)
					{
						n_bytes = 8;
						ctx_type += 0x80;
					}

					for (int j = 0; j < n_bytes; ++j)
					{
						auto rc = find_rc_context(m_ctx_rc_numeric_size, ctx_type + j, tpl_ctx_rc_numeric_size);
						delta += (int64_t) rc->Decode() << (8 * j);
					}
					
					if (d >= 251)
						delta = -delta;
				}

				v_deltas[i] = delta;

				v_cur = delta + v_prev;

				size += store_int(p, size, v_cur);
			}

			p[size++] = get<1>(v_tokens_prev[i]);
		}

		tokenize(p, size, v_tokens_cur);
	}
	else
	{
		ctx_flags = ((ctx_flags << 1) + 0) & 0xff;

		// Decode id plain
		for (uint32_t i = 0; ; ++i)
		{
			auto rc = find_rc_context(m_ctx_rc_plain, i, tpl_ctx_rc_plain);
			p[i] = rc->Decode();

			if (p[i] == 0 || p[i] == 0xA)
			{
				size = i+1;
				break;
			}
		}

		tokenize(p, size, v_tokens_cur);

		v_deltas.clear();
		v_deltas.resize(v_tokens_cur.size(), 0u);
	}

	v_tokens_prev.swap(v_tokens_cur);
	id_prev = p;
	size_prev = size;
}

//*****************************************************************************************************
void CIdCompressor::decompress_instrument(uint8_t *p, uint32_t &size)
{
	uint32_t flag = rc_id_mtf_flag->Decode();
	int mtf_code = 0;
	string instr_name;

	if (flag == 0)			// new instrument name
	{
		decompress_lossless(p, size);
		instr_name = string((char*) p);
		p[size-1] = '\n';

		mtf.Insert(instr_name);
	}
	else if (flag < 3)
	{
		mtf_code = flag - 1;
	}
	else if (flag == 3)
	{
		mtf_code = rc_id_mtf_code[0]->Decode() + 2;
	}
	else if (flag == 4)
	{
		mtf_code = rc_id_mtf_code[1]->Decode() + 4;
	}
	else if (flag == 5)
	{
		mtf_code = rc_id_mtf_code[2]->Decode() + 8;
	}
	else if (flag == 6)
	{
		mtf_code = rc_id_mtf_code[3]->Decode() + 16;
	}
	else if (flag == 7)
	{
		mtf_code = rc_id_mtf_code[4]->Decode() + 32;
	}
	else if (flag == 8)
	{
		mtf_code = rc_id_mtf_code[5]->Decode() + 64;
	}
	else if (flag == 9)
	{
		mtf_code = rc_id_mtf_code[6]->Decode() + 128;
	}
	else if (flag == 10)
	{
		mtf_code = 0;

		for (int i = 0; i < 4; ++i)
			mtf_code += ((uint32_t)rc_id_mtf_byte[i]->Decode()) << (8 * i);
	}

	if (flag > 0)
	{
		mtf.Insert(mtf_code);
		auto s = mtf.GetString(0);
		size = (uint32_t) s.size();
		copy_n(s.begin(), size, p);
		p[size++] = '\n';
	}
}

//*****************************************************************************************************
uint32_t CIdCompressor::tokenize(const uint8_t *p, const uint32_t size, v_tokens_t &v_tokens)
{
	uint32_t token_start_pos = 0;
	token_type_t token_type = token_type_t::numeric;

	v_tokens.clear();

	for (uint32_t i = 0; i < size; ++i)
	{
		if (!is_literal(p[i]))
		{
			if (token_type == token_type_t::numeric && ((i - token_start_pos >= 11) || (i == token_start_pos)))
				token_type = token_type_t::literal;
			v_tokens.emplace_back(make_tuple(token_type, p[i], token_start_pos, i));

			token_type = token_type_t::numeric;
			token_start_pos = i + 1;
		}
		else if (!is_numeric(p[i]))
			token_type = token_type_t::literal;
	}

	return (uint32_t) v_tokens.size();
}


//*****************************************************************************************************
CIdCompressor::ctx_map_2_t::value_type CIdCompressor::find_rc_context(ctx_map_2_t &m_ctx_rc, context_t ctx, CRangeCoderModelFixedSize<CVectorIOStream, 2> *tpl, const stats_t *counts)
{
	auto p = m_ctx_rc.find(ctx);

	if (p == nullptr)
		m_ctx_rc.insert(ctx, p = new CRangeCoderModelFixedSize<CVectorIOStream, 2>(*tpl));

	return p;
}

//*****************************************************************************************************
CIdCompressor::ctx_map_4_t::value_type CIdCompressor::find_rc_context(ctx_map_4_t &m_ctx_rc, context_t ctx, CRangeCoderModelFixedSize<CVectorIOStream, 4> *tpl, const stats_t *counts)
{
	auto p = m_ctx_rc.find(ctx);

	if (p == nullptr)
		m_ctx_rc.insert(ctx, p = new CRangeCoderModelFixedSize<CVectorIOStream, 4>(*tpl));

	return p;
}
//*****************************************************************************************************
CIdCompressor::ctx_map_128_t::value_type CIdCompressor::find_rc_context(ctx_map_128_t &m_ctx_rc, context_t ctx, CRangeCoderModelFixedSize<CVectorIOStream, 128> *tpl, const stats_t *counts)
{
	auto p = m_ctx_rc.find(ctx);

	if (p == nullptr)
		m_ctx_rc.insert(ctx, p = new CRangeCoderModelFixedSize<CVectorIOStream, 128>(*tpl));

	return p;
}

//*****************************************************************************************************
CIdCompressor::ctx_map_256_t::value_type CIdCompressor::find_rc_context(ctx_map_256_t &m_ctx_rc, context_t ctx, CRangeCoderModelFixedSize<CVectorIOStream, 256> *tpl, const stats_t *counts)
{
	auto p = m_ctx_rc.find(ctx);

	if (p == nullptr)
		m_ctx_rc.insert(ctx, p = new CRangeCoderModelFixedSize<CVectorIOStream, 256>(*tpl));

	return p;
}

//*****************************************************************************************************
CIdCompressor::ctx_map_t::value_type CIdCompressor::find_rc_context(ctx_map_t &m_ctx_rc, context_t ctx, CRangeCoderModel<CVectorIOStream> *tpl, const stats_t *counts)
{
	auto p = m_ctx_rc.find(ctx);

	if (p == nullptr)
		m_ctx_rc.insert(ctx, p = new CRangeCoderModel<CVectorIOStream>(*tpl));

	return p;
}

// EOF