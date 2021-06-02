add_executable (bruteforce3p1_angles bruteforce3p1_angles.cc)
#add_executable (noosc noosc.cc)
add_executable (HV_3p1 hypervolume_threeplusone.cc)
add_executable (PL_v1 plotLooper_v1.cc)
add_executable (uBFitter uBfitter.cc)
add_executable (uBFitter_frequensens uBfitter_frequensens.cc)

include_directories (${CMAKE_SOURCE_DIR}/inc)
include_directories (${CMAKE_SOURCE_DIR}/tinyxml)

target_link_libraries (bruteforce3p1_angles DatasetLib)
target_link_libraries (bruteforce3p1_angles TinyXMLLib)

#target_link_libraries (noosc DatasetLib)
#target_link_libraries (noosc TinyXMLLib)

target_link_libraries (PL_v1 DatasetLib)
target_link_libraries (PL_v1 TinyXMLLib)

target_link_libraries (HV_3p1 DatasetLib)
target_link_libraries (HV_3p1 TinyXMLLib)

target_link_libraries (uBFitter DatasetLib)
target_link_libraries (uBFitter TinyXMLLib)
target_link_libraries (uBFitter ubLib)

target_link_libraries (uBFitter_frequensens DatasetLib)
target_link_libraries (uBFitter_frequensens TinyXMLLib)
target_link_libraries (uBFitter_frequensens ubLib)
