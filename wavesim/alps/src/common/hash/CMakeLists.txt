add_library(alps_hash INTERFACE)
target_sources(alps_hash PUBLIC hash.h)
target_link_libraries(alps_hash INTERFACE alps_src_root)
target_formatter_setup(alps_hash)
