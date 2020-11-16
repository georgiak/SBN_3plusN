add_library (DatasetLib NuMI.cc NOMAD.cc XSec.cc KARMEN.cc LSND_loglikelihood.cc MINOS.cc MiniBooNE_combined.cc CDHS.cc DANSS.cc NEOS.cc NuMI.cc CCFR.cc Bugey.cc NOMAD.cc KARMEN.cc MiniBooNE_dis.cc LSND_loglikelihood.cc Gallium.cc OscTools.cc Atm.cc XMLTools.cc)

# on hold:  MiniBooNE.cc, PROSPECT.cc, FromChi2Surf.cc,

include_directories ( ${CMAKE_SOURCE_DIR}/inc)
include_directories ( ${CMAKE_SOURCE_DIR}/tinyxml)

target_link_libraries(DatasetLib ${ROOT_LIBRARIES})
target_link_libraries(DatasetLib ${GSL_LIBRARIES})
