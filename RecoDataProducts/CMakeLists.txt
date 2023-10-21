cet_make_library(
    SOURCE
      src/AlgorithmID.cc
      src/BkgClusterFlag.cc
      src/BkgQual.cc
      src/ComboHit.cc
      src/CosmicTrack.cc
      src/CrvDigi.cc
      src/CrvRecoPulse.cc
      src/ExtMonFNALRawCluster.cc
      src/ExtMonFNALRawHit.cc
      src/ExtMonFNALRecoCluster.cc
      src/ExtMonFNALRecoClusterCollection.cc
      src/ExtMonFNALTrkFitQuality.cc
      src/ExtMonFNALTrkParam.cc
      src/ExtMonUCITofHit.cc
      src/HelixVal.cc
      src/KalSeed.cc
      src/KalSegment.cc
      src/PIDProduct.cc
      src/RobustHelix.cc
      src/StrawDigi.cc
      src/StrawDigiFlag.cc
      src/StrawHit.cc
      src/StrawHitFlag.cc
      src/TrackClusterMatch.cc
      src/TrackSummary.cc
      src/TrkCaloHitPID.cc
      src/TrkExtTraj.cc
      src/TrkExtTrajPoint.cc
      src/TrkFitDirection.cc
      src/TrkFitFlag.cc
      src/TrkQual.cc
      src/TrkToCaloExtrapol.cc
    LIBRARIES PUBLIC
      BTrk_BbrGeom
      BTrk_KalmanTrack
      BTrk_ProbTools
      KinKal_Trajectory
      
      offline::DataProducts
      offline::GeneralUtilities
      offline::KinKalGeom
      offline::Mu2eKinKal
      offline::TrackerConditions
)

art_dictionary( NO_CHECK_CLASS_VERSION # For some reason this segfaults
    CLASSES_DEF_XML ${CMAKE_CURRENT_SOURCE_DIR}/src/classes_def.xml
    CLASSES_H ${CMAKE_CURRENT_SOURCE_DIR}/src/classes.h
     DICTIONARY_LIBRARIES
      offline::RecoDataProducts
)
install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
