cmake_minimum_required(VERSION 3.16)
project(SpeedPPR)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-march=native -msse2 -msse -DHAVE_SSE2 -ffast-math -mfpmath=sse -use_fast_math -pthread -std=c++11 -O3 -DNDEBUG")
#set(CMAKE_CXX_FLAGS "-march=native -ffast-math -use_fast_math -pthread -std=c++11 -O3 -DNDEBUG")
#set(CMAKE_CXX_FLAGS "-march=native -ffast-math -use_fast_math -pthread -std=c++11")

set(SFMT_INCLUDE_DIR "SFMT-src-1.5.1/")
INCLUDE_DIRECTORIES(${SFMT_INCLUDE_DIR})
set(FAST_DOUBLE_PARSER_INCLUDE_DIR "fast_double_parser-master/include")
INCLUDE_DIRECTORIES(${FAST_DOUBLE_PARSER_INCLUDE_DIR})

add_executable(SpeedPPR main.cpp
        Graph.h BasicDefinition.h
        HelperFunctions.h IteratedMethods.h HelperFunctions.cpp
        QueryLoader.h
        BatchRandomWalk.cpp BatchRandomWalk.h
        MyQueue.h MyRandom.cpp MyRandom.h
        SpeedPPR.cpp SpeedPPR.h
        SFMT-src-1.5.1/SFMT.c SFMT-src-1.5.1/SFMT.h SFMT-src-1.5.1/SFMT-sse2.h
        CleanGraph.h QueryGenerator.h)
target_compile_definitions(SpeedPPR PRIVATE SFMT_MEXP=19937)