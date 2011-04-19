cmake_minimum_required(VERSION 2.8)

project (oifits-sim)

set(CMAKE_BUILD_TYPE DEBUG)

#set(CMAKE_CXX_FLAGS "-Wall -Wcast-align -Wcast-qual -Wextra -Wfloat-equal -Wredundant-decls -Wshadow -Wwrite-strings -pedantic -Wno-long-long")

#set(CMAKE_CXX_FLAGS "-Wall -pedantic")

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

find_package(GLIB2 REQUIRED)
include_directories(${GLIB2_INCLUDE_DIRS})
set(LIBS ${LIBS} ${GLIB2_LIBRARIES})

find_package(CFITSIO REQUIRED)
include_directories(${CFITSIO_INCLUDE_DIRS})
set(LIBS ${LIBS} ${CFITSIO_LIBRARIES})

find_package(FFTW REQUIRED)
include_directories(${FFTW_INCLUDE_DIRS})
set(LIBS ${LIBS} ${FFTW_LIBRARIES})

include_directories(../lib/textio)
link_directories(. ../lib/textio)
set(LIBS ${LIBS} textio)

include_directories(../lib/oifitslib)
link_directories(. ../lib/oifitslib)
set(LIBS ${LIBS} oifits)

file(GLOB_RECURSE CPP_SOURCE . *.cpp)
file(GLOB_RECURSE C_SOURCE . *.c)

add_executable(oifits-sim ${CPP_SOURCE} ${C_SOURCE})

#SET_TARGET_PROPERTIES(oifits-sim PROPERTIES LINKER_LANGUAGE CPP)

target_link_libraries(oifits-sim ${LIBS})


