set(HEADERS
    collectives.h
    comm.h
    config.h
    constants.h
    datatype.h
    environment.h
    operator.h
    point2point.h
    request.h
    rma_window.h
    shared_window.h
    status.h
    utility.h
)
set(SOURCES
    environment.cpp
    comm.cpp
    status.cpp
    request.cpp
    rma_window.cpp
    shared_window.cpp
    utility.cpp
)

add_library(mpipp ${HEADERS} ${SOURCES})
target_link_libraries(
  mpipp
  PUBLIC MPI::MPI_CXX
         spdlog
         fmt::fmt
         nonstd::span-lite
         alps_src_root
  PRIVATE enum.hpp
)
target_compile_host_only(mpipp)
target_formatter_setup(mpipp)
target_static_analysis_setup(mpipp)
