add_library (DatasetLib MiniBooNE.cc OscTools.cc)

include_directories ( ${CMAKE_SOURCE_DIR}/inc)

target_link_libraries(DatasetLib ${ROOT_LIBRARIES})
target_link_libraries(DatasetLib ${GSL_LIBRARIES})
