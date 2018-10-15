add_executable (bruteforce3p1 bruteforce3p1.cc)
add_executable (process_3p1 process_threeplusone.cc)
add_executable (ntupleGridder ntupleGridder_old.cc)
add_executable (globalFit globalFit.cc)
add_executable (singlepoint singlepoint.cc)

include_directories (${CMAKE_SOURCE_DIR}/inc)
include_directories (${CMAKE_SOURCE_DIR}/tinyxml)

target_link_libraries (bruteforce3p1 DatasetLib)
target_link_libraries (bruteforce3p1 TinyXMLLib)

target_link_libraries (process_3p1 DatasetLib)
target_link_libraries (process_3p1 TinyXMLLib)

target_link_libraries (ntupleGridder DatasetLib)
target_link_libraries (ntupleGridder TinyXMLLib)

target_link_libraries (globalFit DatasetLib)
target_link_libraries (globalFit TinyXMLLib)

target_link_libraries (singlepoint DatasetLib)
target_link_libraries (singlepoint TinyXMLLib)
