# Define the header and source files (and directories)
set(HEADERS hdf5.h io_fwd.h)

set(SOURCES hdf5.cpp)

add_library(alps_io ${HEADERS} ${SOURCES})
target_link_libraries(
  alps_io
  PUBLIC alps_base alps_decomp Kokkos::kokkoscore HighFive
  INTERFACE alps_src_root
  PRIVATE project_warnings
)
if(NOT ALPS_BUILD_HDF5)
  target_include_directories(alps_io PUBLIC ${HDF5_INCLUDE_DIRS})
endif()
target_formatter_setup(alps_io)
target_static_analysis_setup(alps_io)
