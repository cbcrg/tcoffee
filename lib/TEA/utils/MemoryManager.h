/*
 * MemoryManager.h
 *
 *  Created on: May 1, 2012
 *      Author: ck
 */

#ifndef MEMORYMANAGER_H_
#define MEMORYMANAGER_H_

namespace BioTools {
namespace Utils {


template<int object_size, int chunk_size>
class MemoryManager {

public:
	MemoryManager();
	virtual ~MemoryManager();
};



} /* namespace Utils */
} /* namespace BioTools */
#endif /* MEMORYMANAGER_H_ */
