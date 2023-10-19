cet_make_library(
    SOURCE
src/CosmicCORSIKA.cc
src/ExtMonFNALMARSUtils.cc
src/FromCorsikaBinary_source.cc
src/FromEMFMARSFileWeighted_source.cc
src/FromExtMonFNALMARSFile_source.cc
src/FromSTMTestBeamData_source.cc
src/FromTrackerPrototypeData_source.cc
src/PBISequence_source.cc
src/STMTestBeamFileNameDecoder.cc
src/STMTestBeamHeaders.cc
LIBRARIES PUBLIC
cetlib::cetlib
)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
