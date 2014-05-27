cmake_minimum_required(VERSION 2.6)


project (BioTools C CXX)

set(CMAKE_CXX_FLAGS "-O3 -fopenmp -std=c++0x -Wall -pedantic -Wextra -Wunused ")
#-std=c++11
add_subdirectory(doc)
add_subdirectory(testsuit)



################################################
#         find necessary libraries             #
################################################

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)

FIND_PACKAGE(Boost COMPONENTS system filesystem program_options REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIRS})
ADD_DEFINITIONS( "-DHAS_BOOST" )

FIND_PACKAGE(OpenMP)
 find_package(OpenMP)
 if (OPENMP_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
 endif (OPENMP_FOUND)


IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)
  FIND_PACKAGE(Git)
  IF(GIT_FOUND)
    EXECUTE_PROCESS(
      COMMAND ${GIT_EXECUTABLE} describe --tags --always
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      OUTPUT_VARIABLE "MY_BUILD_VERSION"
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    MESSAGE( STATUS "Git version: ${MY_BUILD_VERSION}" )
  ELSE(GIT_FOUND)
    SET(MY_BUILD_VERSION 0)
  ENDIF(GIT_FOUND)
ENDIF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/Version.hpp.in ${CMAKE_CURRENT_BINARY_DIR}/Version.hpp @ONLY)
include_directories(${CMAKE_CURRENT_BINARY_DIR})



################################################
#                    MSAA                      #
################################################

SET(msaa_src ./src/programs/msaa.cpp ./src/lib/utils/filesystem.cpp ./src/lib/utils/utils.cpp ./src/lib/utils/graph_algorithms.cpp ./src/lib/utils/ScoringMatrix.cpp ./src/lib/Sequence/Sequence.cpp ./src/lib/Sequence/SequenceSet.cpp ./src/lib/Sequence/Alignment.cpp ./src/lib/Sequence/SequenceSet_io.cpp ./src/lib/Sequence/Alignment_io.cpp ./src/lib/Sequence/aln_analysis.cpp)
SET(msaa_exe msaa)
ADD_EXECUTABLE(${msaa_exe} ${msaa_src})
target_link_libraries(${msaa_exe} 
	${Boost_LIBRARIES}
)

################################################
#                    aln_trimmer               #
################################################

SET(aln_trimmer_src ./src/programs/aln_trimmer.cpp ./src/lib/utils/filesystem.cpp ./src/lib/utils/utils.cpp ./src/lib/utils/graph_algorithms.cpp ./src/lib/utils/ScoringMatrix.cpp ./src/lib/Sequence/Sequence.cpp ./src/lib/Sequence/SequenceSet.cpp ./src/lib/Sequence/Alignment.cpp ./src/lib/Sequence/SequenceSet_io.cpp ./src/lib/Sequence/Alignment_io.cpp ./src/lib/Sequence/aln_analysis.cpp)
SET(aln_trimmer_exe aln_trimmer)
ADD_EXECUTABLE(${aln_trimmer_exe} ${aln_trimmer_src})
target_link_libraries(${aln_trimmer_exe} 
	${Boost_LIBRARIES}
)




################################################
#                  setAnalyzer                 #
################################################

SET(setanalyser_src ./src/programs/setAnalyser.cpp ./src/lib/utils/filesystem.cpp ./src/lib/utils/utils.cpp ./src/lib/Sequence/Sequence.cpp ./src/lib/Sequence/SequenceSet.cpp ./src/lib/Sequence/Alignment.cpp ./src/lib/Sequence/SequenceSet_io.cpp ./src/lib/Sequence/Alignment_io.cpp ./src/lib/Sequence/aln_analysis.cpp)
SET(setanalyser_exe setAnalyser)
ADD_EXECUTABLE(${setanalyser_exe} ${setanalyser_src})
target_link_libraries(${setanalyser_exe} 
	${Boost_LIBRARIES}
)



################################################
#                  maskHelper                 #
################################################

SET(maskHelper_src ./src/programs/maskHelper.cpp ./src/lib/utils/ ./src/lib/utils/filesystem.cpp ./src/lib/Sequence/Sequence.cpp ./src/lib/Sequence/SequenceSet.cpp ./src/lib/utils/utils.cpp ./src/lib/Sequence/SequenceSet_io.cpp )
SET(maskHelper_exe maskHelper)
ADD_EXECUTABLE(${maskHelper_exe} ${maskHelper_src})
target_link_libraries(${maskHelper_exe} 
	${Boost_LIBRARIES}
)



################################################
#                     tea                      #
################################################

SET(tea_src ./src/programs/tea/tea.cpp 
             ./src/lib/utils/utils.cpp
             ./src/lib/Clustering/kmeans.cpp 
             ./src/lib/Clustering/clustering.cpp ./src/lib/Clustering/Vector.cpp ./src/lib/Clustering/Tree.cpp  ./src/lib/utils/Matrix.cpp ./src/lib/utils/fast_math.cpp ./src/lib/utils/Library.cpp ./src/lib/utils/ScoringMatrix.cpp ./src/lib/Sequence/Sequence.cpp ./src/lib/Sequence/SequenceSet.cpp ./src/lib/Sequence/aln_analysis.cpp ./src/lib/Sequence/Alignment.cpp ./src/lib/Sequence/SequenceSet_io.cpp ./src/lib/Sequence/Alignment_io.cpp ./src/lib/HMM/HMM.cpp ./src/lib/HMM/fw_bw.cpp ./src/lib/Sequence/gotoh.cpp ./src/lib/Sequence/other_seq.cpp)
SET(tea_exe tea)
ADD_EXECUTABLE(${tea_exe} ${tea_src})
target_link_libraries(${tea_exe} 
	${Boost_LIBRARIES}
)

install(TARGETS ${msaa_exe} ${maskHelper_exe} ${setanalyser_exe} ${tea_exe} DESTINATION .)



################################################
#                Alignment library             #
################################################

FILE (GLOB_RECURSE SOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/src/lib/*.cpp)

SET(aln_lib_src ${SOURCE_FILES})
add_library (biotools_aln ${aln_lib_src})
 




