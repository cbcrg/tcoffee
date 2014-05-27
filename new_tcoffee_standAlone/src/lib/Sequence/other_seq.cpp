/*
 * other_seq.cpp
 *
 *  Created on: Aug 4, 2012
 *      Author: ck
 */

#include "other_seq.h"


namespace BioTools {
namespace Utils {

using namespace BioTools::Container;

struct SkipList_Point_Fragment
{
	int x;
	int y;
	SkipList_Point_Fragment *pre;
	double score;
};


}
}
