find_package(CxxTest)

FIND_PACKAGE(OpenMP)
 find_package(OpenMP)
 if (OPENMP_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
 endif (OPENMP_FOUND)
 
 
if(CXXTEST_FOUND)
include_directories(${CXXTEST_INCLUDE_DIR})
enable_testing()
CXXTEST_ADD_TEST(seq_module_test seq_module_test.cpp ${CMAKE_CURRENT_SOURCE_DIR}/seq_module_test.h)
CXXTEST_ADD_TEST(SkipList_test SkipList_test.cpp ${CMAKE_CURRENT_SOURCE_DIR}/SkipList_test.h)
target_link_libraries(seq_module_test biotools_aln)
target_link_libraries(SkipList_test biotools_aln)
endif()