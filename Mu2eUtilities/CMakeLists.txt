cet_make_library(
    SOURCE
      src/BinnedSpectrum.cc
      src/BuildLinearFitMatrixSums.cc
      src/CaloPulseShape.cc
      src/checkSimParticleCollection.cc
      src/compressPdgId.cc
      src/ConversionSpectrum.cc
      src/CoordinateString.cc
      src/CosmicTrackUtils.cc
      src/CzarneckiSpectrum.cc
      src/EjectedProtonSpectrum.cc
      src/EventWeightHelper.cc
      src/fromStrings.cc
      src/GammaPairConversionSpectrum.cc
      src/GeneratorSummaryHistograms.cc
      src/HelixTool.cc
      src/LsqSums2.cc
      src/LsqSums4.cc
      src/McUtilsToolBase.cc
      src/MedianCalculator.cc
      src/ModuleHistToolBase.cc
      src/MuonCaptureSpectrum.cc
      src/MVATools.cc
      src/ParametricFit.cc
      src/particleEnteringG4Volume.cc
      src/PhysicalVolumeMultiHelper.cc
      src/PionCaptureSpectrum.cc
      src/PointLinePCA_XYZ.cc
      src/PoissonHistogramBinning.cc
      src/polar3Vector.cc
      src/ProtonPulseRandPDF.cc
      src/RandomLimitedExpo.cc
      src/RandomUnitSphere.cc
      src/ReSeedByEventID.cc
      src/rm48.cc
      src/ShankerWatanabeSpectrum.cc
      src/SimParticleCollectionPrinter.cc
      src/SimParticleGetTau.cc
      src/simParticleList.cc
      src/SimParticleParentGetter.cc
      src/SimpleSpectrum.cc
      src/SortedStepPoints.cc
      src/STMUtils.cc
      src/StopWatch.cc
      src/Table.cc
      src/TrackerBFieldInfo.cc
      src/TrackPatRecType.cc
      src/TrackTool.cc
      src/TriggerResultsNavigator.cc
      src/TwoLinePCA.cc
      src/TwoLinePCA_XYZ.cc
      src/VectorVolume.cc
    LIBRARIES PUBLIC
      art_root_io::tfile_support
      Offline::BFieldGeom
      Offline::ConditionsService
      Offline::ConfigTools
      Offline::DataProducts
      Offline::GeneralUtilities
      Offline::GeometryService
      Offline::GeomPrimitives
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::RecoDataProducts
      Offline::SeedService
      Offline::STMConditions
      Offline::StoppingTargetGeom
      Offline::TrackerGeom
      ROOT::Hist
      ROOT::Matrix
      XercesC::XercesC
)


configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/acDipoleTransmissionFunction_20160511.txt   ${CURRENT_BINARY_DIR} data/acDipoleTransmissionFunction_20160511.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/potTimingDistribution_20160511.txt   ${CURRENT_BINARY_DIR} data/potTimingDistribution_20160511.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/potTimingDistribution_38ms.txt   ${CURRENT_BINARY_DIR} data/potTimingDistribution_38ms.txt   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/potTimingDistribution_avg.txt   ${CURRENT_BINARY_DIR} data/potTimingDistribution_avg.txt   COPYONLY)

install(DIRECTORY data DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Offline/Mu2eUtilities)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
