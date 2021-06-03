add_executable (ICtoTH3D ICtoTH3D.cc)

include_directories (${CMAKE_SOURCE_DIR}/inc)
include_directories (${CMAKE_SOURCE_DIR}/tinyxml)

target_link_libraries (ICtoTH3D DatasetLib)
target_link_libraries (ICtoTH3D TinyXMLLib)
