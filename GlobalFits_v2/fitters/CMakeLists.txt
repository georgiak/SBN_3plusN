add_executable (ntuple_metrics ntuple_metrics.cc)
add_executable (bruteforce3p1 bruteforce3p1.cc)
add_executable (process_3p1 process_threeplusone.cc)
add_executable (process_3p3 process_3p3.cc)
add_executable (ntupleGridder ntupleGridder_old.cc)
add_executable (globalFit globalFit.cc)
add_executable (singlepoint singlepoint.cc)
add_executable (process_global process_global.cc)
add_executable (process_ntuple  process_ntuple.cc)

include_directories (${CMAKE_SOURCE_DIR}/inc)
include_directories (${CMAKE_SOURCE_DIR}/tinyxml)

target_link_libraries (bruteforce3p1 DatasetLib)
target_link_libraries (bruteforce3p1 TinyXMLLib)

target_link_libraries (process_3p1 DatasetLib)
target_link_libraries (process_3p1 TinyXMLLib)

target_link_libraries (process_3p3 DatasetLib)
target_link_libraries (process_3p3 TinyXMLLib)

target_link_libraries (ntupleGridder DatasetLib)
target_link_libraries (ntupleGridder TinyXMLLib)

target_link_libraries (globalFit DatasetLib)
target_link_libraries (globalFit TinyXMLLib)

target_link_libraries (singlepoint DatasetLib)
target_link_libraries (singlepoint TinyXMLLib)

target_link_libraries (process_global DatasetLib)
target_link_libraries (process_global TinyXMLLib)

target_link_libraries (process_ntuple DatasetLib)
target_link_libraries (process_ntuple TinyXMLLib)

target_link_libraries (ntuple_metrics DatasetLib)
target_link_libraries (ntuple_metrics TinyXMLLib)
