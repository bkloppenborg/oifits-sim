cmake_minimum_required(VERSION 2.6)
project(oifits)

# Include the required libraires for oifitslib:
find_package(GLIB2 REQUIRED)
include_directories(${GLIB2_INCLUDE_DIRS})

find_package(CFITSIO REQUIRED)
include_directories(${CFITSIO_INCLUDE_DIRS})

# Now build the library:
add_library(oifits read_fits write_fits free_fits oifile oifilter oicheck oimerge)
target_link_libraries(oifits ${GLIB2_LIBRARIES} ${CFITSIO_LIBRARIES})
