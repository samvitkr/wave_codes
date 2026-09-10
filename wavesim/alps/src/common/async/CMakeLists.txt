set(HEADERS
    thread_pool.h
    event_pool.h
    event_openmp.h
    event_cuda.h
    event.h
    streams.h
)
set(SOURCES thread_pool.cpp event_pool.cpp event_openmp.cpp streams.cpp)
if(${PROJECT_NAME_UPPERCASE}_ENABLE_CUDA OR ${PROJECT_NAME_UPPERCASE}_ENABLE_HIP)
  list(APPEND SOURCES event_cuda.cpp)
endif()
add_library(alps_async STATIC ${HEADERS} ${SOURCES})
target_link_libraries(
  alps_async
  PUBLIC Kokkos::kokkoscore alps_base
  INTERFACE alps_src_root
  PRIVATE project_warnings
)
target_formatter_setup(alps_async)
