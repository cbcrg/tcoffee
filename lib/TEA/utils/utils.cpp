/*
 * utils.cpp
 *
 *  Created on: Jul 28, 2012
 *      Author: Carsten Kemena
 *
 *   Copyright 2011 Carsten Kemena
 */

/*
 *
 * This file is part of BioTools++.
 *
 * BioTools++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BioTools++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BioTools++.  If not, see <http://www.gnu.org/licenses/>.
 *
 */


#include "utils.h"

namespace BioTools {
namespace Utils {

using namespace std;

StrTok::StrTok(char *str):_str(str)
{}


StrTok::StrTok()
{}


char *
StrTok::next(const char *sep)
{
	if (*_str == '\0')
		return NULL;

	int go_on = 1;
	const char *tmp_sep;
	while (*_str != '\0' && go_on)
	{

		go_on = 0;
		tmp_sep = sep;
		while (*tmp_sep != '\0')
		{
			if (*tmp_sep == *_str)
				go_on = 1;
			++tmp_sep;
		}
		if (go_on)
			++_str;
	}

	if (*_str == '\0')
		return NULL;
	char *ret=_str;
	go_on = 1;
	while (*_str != '\0' && go_on)
	{
		tmp_sep = sep;
		while (*tmp_sep != '\0')
		{
			if (*tmp_sep == *_str)
				go_on = 0;
			++tmp_sep;
		}
		if (go_on)
			++_str;
	}
	if (*_str != '\0')
	{
		*_str = '\0';
		++_str;
	}
	return ret;
}


void
str_upper(string &str)
{
	size_t str_len = str.length();
	for (unsigned int i = 0; i < str_len; ++i)
	{
	    str[i]=toupper(str[i]);
	}
}

void
str_lower(string &str)
{
	char c;
	size_t str_len = str.length();
	for (unsigned int i = 0; i < str_len; ++i)
	{
		c=str[i];
	    str[i]=tolower(c);
	}
}


string
get_alphabet(std::string alphabet)
{
	if (alphabet=="full")
		return string("A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y");
	if (alphabet=="SE-B-14")
		return string("A,C,D,EQ,FY,G,H,IV,KR,LM,N,P,ST,W");
	if (alphabet=="SE-B-10")
		return string("AST,C,DN,EQ,FY,G,HW,ILMV,KR,P");
	if (alphabet=="SE-V-10")
		return string("AST,C,DEN,FY,G,H,ILMV,KQR,P,W");
	if (alphabet=="Li-A-10")
		return string("AC,DE,FWY,G,HN,IV,KQR,LM,P,ST");
	if (alphabet=="Li-B-10")
		return string("AST,C,DEQ,FWY,G,HN,IV,KR,LM,P");
	if (alphabet=="Solis-D-10")
		return string("AM,C,DNS,EKQR,F,GP,HT,IV,LY,W");
	if (alphabet=="Solis-G-10")
		return string("AEFIKLMQRVW,C,D,G,H,N,P,S,T,Y");
	if (alphabet=="Murphy-10")
		return string("A,C,DENQ,FWY,G,H,ILMV,KR,P,ST");
	if (alphabet=="SE-B-8")
		return string("AST,C,DHN,EKQR,FWY,G,ILMV,P");
	if (alphabet=="SE-B-6")
		return string("AST,CP,DEHKNQR,FWY,G,ILMV");
	if (alphabet=="Dayhoff-6")
		return string("AGPST,C,DENQ,FWY,HKR,ILMV");
	return "";
}



short *
encode(std::string alphabet)
{
	string alpha = get_alphabet(alphabet);
	if (alpha.empty())
		return NULL;

	short *encoded = new short[256];
	size_t len=alpha.size();
	size_t i,x=0;
	for (i=0; i<len; ++i)
	{
		if (alpha[i]==',')
			++x;
	}
	++x;
	for (i=0; i<256; ++i)
		encoded[i] = x;
	x=0;
	for (i=0; i<len; ++i)
	{
		if (alpha[i] !=',')
		{
			encoded[static_cast<short>(alpha[i])] = x;
			encoded[static_cast<short>(tolower(alpha[i]))] = x;
		}
		else
			++x;
	}
	return encoded;
}


} /* namespace Utils */
} /* namespace BioTools */


