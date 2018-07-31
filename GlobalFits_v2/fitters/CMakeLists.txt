add_executable (globalFit globalFit_v2.cc)
add_executable (process_3p1 process_threeplusone.cc)

include_directories (${CMAKE_SOURCE_DIR}/inc)
include_directories (${CMAKE_SOURCE_DIR}/tinyxml)

target_link_libraries (globalFit DatasetLib)
target_link_libraries (globalFit TinyXMLLib)

target_link_libraries (process_3p1 DatasetLib)
target_link_libraries (process_3p1 TinyXMLLib)
