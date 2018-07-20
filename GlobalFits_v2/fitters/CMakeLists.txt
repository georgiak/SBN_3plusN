add_executable (globalFit globalFit_v2.cc)

include_directories (${CMAKE_SOURCE_DIR}/inc)

target_link_libraries (globalFit DatasetLib)
