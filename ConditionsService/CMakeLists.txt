cet_make_library(
    SOURCE
src/AcceleratorParams.cc
src/CalorimeterCalibrations.cc
src/ConditionsService_service.cc
src/ExtMonFNALConditions.cc
LIBRARIES PUBLIC
cetlib::cetlib
)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
