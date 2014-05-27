/*
 * SkipList_test.h
 *
 *  Created on: Aug 8, 2012
 *      Author: Carsten Kemena
 *
 *   Copyright 2011 Carsten Kemena
 *
 *
 * This file is part of BioTools.
 *
 * BioTools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BioTools is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BioTools.  If not, see <http://www.gnu.org/licenses/>.
 *
 */


// C header
#include <cstdlib>

// CxxTest header
#include <cxxtest/TestSuite.h>

#include "../src/lib/Container/SkipList.h"

using namespace std;
using namespace BioTools;
using namespace BioTools::Container;

class SkipList_Test : public CxxTest::TestSuite
{
public:

	void test_insert_find()
	{
		SkipList<int, int> list(9,0.5);
		list.insert(5, 6);
		TS_ASSERT_EQUALS(list.find(5)->value, 6);
		list.insert(4, 5);
		TS_ASSERT_EQUALS(list.find(5)->value, 6);
		TS_ASSERT_EQUALS(list.find(4)->value, 5);
		list.insert(7, 8);
		TS_ASSERT_EQUALS(list.find(7)->value, 8);
		list.insert(2, 3);
		TS_ASSERT_EQUALS(list.find(2)->value, 3);
	}

	void test_find_smaller()
	{
		SkipList<int, int> list(9,0.5);
		list.insert(5, 6);
		list.insert(4, 5);
		TS_ASSERT_EQUALS(list.find_smaller(4), list.end());
		list.insert(7, 8);
		TS_ASSERT_EQUALS(list.find_smaller(10)->value, 8);
		list.insert(2, 3);
	}

};
