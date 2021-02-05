#add_executable (ntuple_metrics ntuple_metrics.cc)
#add_executable (bruteforce3p1 bruteforce3p1.cc)
#add_executable (process_3p1 process_threeplusone.cc)
#add_executable (process_3p3 process_3p3.cc)
add_executable (bruteforce3p1_angles bruteforce3p1_angles.cc)
add_executable (noosc noosc.cc)
add_executable (PL_appearance plotLooper_appearance.cc)
add_executable (PL_nuedis plotLooper_nuedis.cc)
add_executable (PL_numudis plotLooper_numudis.cc)
add_executable (PL_numudis_exclusion plotLooper_numudis_exclusion.cc)
add_executable (HV_3p1 hypervolume_threeplusone.cc)
add_executable (ICtoTH3D ICtoTH3D.cc)
add_executable (PL_v1 plotLooper_v1.cc)

#add_executable (ntupleGridder ntupleGridder_old.cc)
#add_executable (globalFit globalFit.cc)
#add_executable (singlepoint singlepoint.cc)
#add_executable (process_global process_global.cc)
#add_executable (process_ntuple  process_ntuple.cc)

include_directories (${CMAKE_SOURCE_DIR}/inc)
include_directories (${CMAKE_SOURCE_DIR}/tinyxml)

target_link_libraries (bruteforce3p1_angles DatasetLib)
target_link_libraries (bruteforce3p1_angles TinyXMLLib)

target_link_libraries (noosc DatasetLib)
target_link_libraries (noosc TinyXMLLib)

target_link_libraries (PL_v1 DatasetLib)
target_link_libraries (PL_v1 TinyXMLLib)

target_link_libraries (PL_appearance DatasetLib)
target_link_libraries (PL_appearance TinyXMLLib)

target_link_libraries (PL_nuedis DatasetLib)
target_link_libraries (PL_nuedis TinyXMLLib)

target_link_libraries (PL_numudis DatasetLib)
target_link_libraries (PL_numudis TinyXMLLib)

target_link_libraries (PL_numudis_exclusion DatasetLib)
target_link_libraries (PL_numudis_exclusion TinyXMLLib)

target_link_libraries (HV_3p1 DatasetLib)
target_link_libraries (HV_3p1 TinyXMLLib)

target_link_libraries (ICtoTH3D DatasetLib)
target_link_libraries (ICtoTH3D TinyXMLLib)

#target_link_libraries (process_3p1 DatasetLib)
#target_link_libraries (process_3p1 TinyXMLLib)

#target_link_libraries (process_3p3 DatasetLib)
#target_link_libraries (process_3p3 TinyXMLLib)

#target_link_libraries (ntupleGridder DatasetLib)
#target_link_libraries (ntupleGridder TinyXMLLib)

#target_link_libraries (globalFit DatasetLib)
#target_link_libraries (globalFit TinyXMLLib)

#target_link_libraries (singlepoint DatasetLib)
#target_link_libraries (singlepoint TinyXMLLib)

#target_link_libraries (process_global DatasetLib)
#target_link_libraries (process_global TinyXMLLib)

#target_link_libraries (process_ntuple DatasetLib)
#target_link_libraries (process_ntuple TinyXMLLib)

#target_link_libraries (ntuple_metrics DatasetLib)
#target_link_libraries (ntuple_metrics TinyXMLLib)
