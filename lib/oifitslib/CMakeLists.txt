cmake_minimum_required(VERSION 2.6)
project(oifits)

# Include the required libraires for oifitslib:
find_package(GLIB2 REQUIRED)
include_directories(${GLIB2_INCLUDE_DIRS})
set(LIBS ${LIBS} ${GLIB2_LIBRARIES})

find_package(CFITSIO REQUIRED)
include_directories(${CFITSIO_INCLUDE_DIRS})
set(LIBS ${LIBS} ${CFITSIO_LIBRARIES})

# Now build the library:
add_library(oifits read_fits write_fits free_fits oifile oifilter oicheck oimerge)
