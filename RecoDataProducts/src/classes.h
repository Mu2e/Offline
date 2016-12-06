//
// Build a dictionary.
//
// $Id: classes.h,v 1.29 2014/05/27 20:09:35 gandr Exp $
// $Author: gandr $
// $Date: 2014/05/27 20:09:35 $
//
// Original author Rob Kutschke
//

#define ENABLE_MU2E_GENREFLEX_HACKS

#include <vector>
//#include <array>

#include "art/Persistency/Common/Wrapper.h"
#include "art/Persistency/Common/Assns.h"


#include "RecoDataProducts/inc/CaloDigiPacked.hh"
#include "RecoDataProducts/inc/CaloDigiCollection.hh"
#include "RecoDataProducts/inc/RecoCaloDigiCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StereoHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/HelixHit.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StrawClusterCollection.hh"
#include "RecoDataProducts/inc/CrvRecoPulsesCollection.hh"
#include "RecoDataProducts/inc/CrvCoincidenceCheckResult.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/KalRepPayloadCollection.hh"
#include "RecoDataProducts/inc/KalRepExtensionPayloadCollection.hh"
#include "RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "RecoDataProducts/inc/TrackSummaryRecoMap.hh"
#include "RecoDataProducts/inc/ExtMonUCITofHitCollection.hh"
#include "RecoDataProducts/inc/SubEventCollection.hh"
#include "RecoDataProducts/inc/TimeClusterCollection.hh"
#include "RecoDataProducts/inc/HelixVal.hh"
#include "RecoDataProducts/inc/RobustHelix.hh"
#include "RecoDataProducts/inc/TrkFitFlag.hh"
#include "RecoDataProducts/inc/ZRotStrawHitMapCollection.hh"
#include "RecoDataProducts/inc/TrackerHitByID.hh"
#include "RecoDataProducts/inc/CaloProtoClusterCollection.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "RecoDataProducts/inc/TrackSeedCollection.hh"
#include "RecoDataProducts/inc/HelixSeedCollection.hh"
#include "RecoDataProducts/inc/TrkExtTrajCollection.hh"
#include "RecoDataProducts/inc/TrkCaloIntersectCollection.hh"
#include "RecoDataProducts/inc/TrkCaloMatchCollection.hh"
#include "RecoDataProducts/inc/AvikPIDProductCollection.hh"
#include "RecoDataProducts/inc/AvikPIDNewProductCollection.hh"
#include "RecoDataProducts/inc/PIDProductCollection.hh"
#include "RecoDataProducts/inc/ExtMonFNALRawHit.hh"
#include "RecoDataProducts/inc/ExtMonFNALRawHitCollection.hh"
#include "RecoDataProducts/inc/ExtMonFNALRawCluster.hh"
#include "RecoDataProducts/inc/ExtMonFNALRawClusterCollection.hh"
#include "RecoDataProducts/inc/ExtMonFNALRecoCluster.hh"
#include "RecoDataProducts/inc/ExtMonFNALRecoClusterCollection.hh"

#include "RecoDataProducts/inc/ExtMonFNALTrkParam.hh"
#include "RecoDataProducts/inc/ExtMonFNALTrkFitQuality.hh"
#include "RecoDataProducts/inc/ExtMonFNALTrkClusterResiduals.hh"
#include "RecoDataProducts/inc/ExtMonFNALTrkFit.hh"
#include "RecoDataProducts/inc/ExtMonFNALTrkFitCollection.hh"

#include "RecoDataProducts/inc/StrawDigiCollection.hh"
#include "RecoDataProducts/inc/TrackSummary.hh"

#include "DataProducts/inc/CRSScintillatorBarIndex.hh"

#include "RecoDataProducts/inc/TrackCaloAssns.hh" 
// BTrk Kalman fit persistence
#include "RecoDataProducts/inc/KalSeedCollection.hh"

// Cannot use the typedefs in here - not sure why.
template class art::Ptr<mu2e::CaloHit>;
template class std::vector<art::Ptr<mu2e::CaloHit> >;
template class art::Ptr<mu2e::ExtMonUCITofHit>;
template class std::vector<art::Ptr<mu2e::ExtMonUCITofHit> >;
template class art::Ptr<mu2e::StrawHit>;
template class std::vector<art::Ptr<mu2e::StrawHit> >;
template class art::Ptr<mu2e::CaloDigiPacked>;
template class std::vector<art::Ptr<mu2e::CaloDigiPacked> >;
template class art::Ptr<mu2e::CaloDigi>;
template class std::vector<art::Ptr<mu2e::CaloDigi> >;
template class art::Ptr<mu2e::RecoCaloDigi>;
template class std::vector<art::Ptr<mu2e::RecoCaloDigi> >;
template class art::Ptr<mu2e::StereoHit>;
template class std::vector<art::Ptr<mu2e::StereoHit> >;
template class art::Ptr<mu2e::StrawHitPosition>;
template class std::vector<art::Ptr<mu2e::StrawHitPosition> >;
template class art::Ptr<mu2e::HelixHit>;
template class std::vector<art::Ptr<mu2e::HelixHit> >;
template class art::Ptr<mu2e::TrkFitFlag>;
template class art::Ptr<mu2e::StrawHitFlag>;
template class std::vector<art::Ptr<mu2e::StrawHitFlag> >;
template class art::Ptr<mu2e::TimeCluster>;
template class std::vector<art::Ptr<mu2e::TimeCluster> >;
template class art::Ptr<mu2e::RobustHelix>;
template class art::Ptr<mu2e::HelixVal>;
template class std::vector<art::Ptr<mu2e::HelixVal> >;
template class art::Ptr<mu2e::HelixCov>;
template class std::vector<art::Ptr<mu2e::HelixCov> >;
template class std::multimap<unsigned long int, art::Ptr<mu2e::StrawHit> >;
template class art::Ptr<mu2e::CaloCrystalHit>;
template class std::vector<art::Ptr<mu2e::CaloCrystalHit> >;
template class art::Ptr<mu2e::CaloProtoCluster>;
template class std::vector<art::Ptr<mu2e::CaloProtoCluster> >;
template class art::Ptr<mu2e::CaloCluster>;
template class std::vector<art::Ptr<mu2e::CaloCluster> >;
template class art::Ptr<mu2e::TrackSeed>;
template class std::vector<art::Ptr<mu2e::TrackSeed> >;
template class art::Ptr<mu2e::TrkExtTrajPoint>;
template class std::vector<mu2e::TrkExtTrajPoint>;
template class art::Ptr<mu2e::TrkExtTraj>;
template class std::vector<mu2e::TrkExtTraj>;
template class art::Ptr<mu2e::TrkCaloIntersect>;
template class std::vector<mu2e::TrkCaloIntersect>;
template class art::Ptr<mu2e::TrkCaloMatch>;
template class std::vector<mu2e::TrkCaloMatch>;
template class art::Ptr<mu2e::PIDProduct>;
template class std::vector<mu2e::PIDProduct>;
template class art::Ptr<mu2e::AvikPIDProduct>;
template class std::vector<mu2e::AvikPIDProduct>;
template class art::Ptr<mu2e::AvikPIDNewProduct>;
template class std::vector<mu2e::AvikPIDNewProduct>;
template class std::vector<std::pair<unsigned int, unsigned int> >;
template class std::vector<mu2e::CrvCoincidenceCheckResult::CoincidenceCombination>;
template class std::vector<mu2e::CrvRecoPulses::CrvSingleRecoPulse>;
template class std::pair<mu2e::CRSScintillatorBarIndex,mu2e::CrvRecoPulses>;
template class std::map<mu2e::CRSScintillatorBarIndex,mu2e::CrvRecoPulses>;
template class art::Ptr<KalRep>;
template class std::vector<art::Ptr<KalRep> >;

template class art::Wrapper<mu2e::StrawHitCollection>;
template class art::Wrapper<mu2e::CaloDigiPacked>;
template class art::Wrapper<mu2e::CaloDigi>;
template class art::Wrapper<mu2e::CaloDigiCollection>;
template class art::Wrapper<mu2e::RecoCaloDigi>;
template class art::Wrapper<mu2e::RecoCaloDigiCollection>;
template class art::Wrapper<mu2e::StereoHitCollection>;
template class art::Wrapper<mu2e::StrawHitPositionCollection>;
template class art::Wrapper<mu2e::HelixHitCollection>;
template class art::Wrapper<mu2e::StrawHitFlagCollection>;
template class art::Wrapper<mu2e::StrawClusterCollection>;
template class art::Wrapper<mu2e::CaloHitCollection>;
template class art::Wrapper<mu2e::CaloCrystalHitCollection>;
template class art::Wrapper<mu2e::CaloProtoClusterCollection>;
template class art::Wrapper<mu2e::CaloClusterCollection>;
template class art::Wrapper<mu2e::TrkCaloIntersectCollection>;
template class art::Wrapper<mu2e::TrkCaloMatchCollection>;
template class art::Wrapper<mu2e::KalRepPayloadCollection>;
template class art::Wrapper<mu2e::KalRepExtensionPayloadCollection>;
template class art::Wrapper<mu2e::KalRepPtrCollection>;
template class art::Wrapper<mu2e::TimeCluster>;
template class art::Wrapper<mu2e::TimeClusterCollection>;
template class art::Ptr<mu2e::HelixSeed>;
template class art::Wrapper<mu2e::HelixSeed>;
template class art::Wrapper<mu2e::HelixSeedCollection>;

template class std::pair<art::Ptr<art::Ptr<KalRep> >, art::Ptr<mu2e::TrackSummary> >;
template class std::pair<art::Ptr<mu2e::TrackSummary>, art::Ptr<art::Ptr<KalRep> > >;
template class art::Assns<art::Ptr<KalRep>, mu2e::TrackSummary>;
template class art::Assns<mu2e::TrackSummary, art::Ptr<KalRep> >;
template class art::Wrapper<art::Assns<art::Ptr<KalRep>, mu2e::TrackSummary> >;
template class art::Wrapper<art::Assns<mu2e::TrackSummary, art::Ptr<KalRep> > >;

template class art::Wrapper<mu2e::ExtMonUCITofHitCollection>;
template class art::Wrapper<mu2e::SubEventCollection>;
template class art::Wrapper<mu2e::ZRotStrawHitMapCollection>;
template class art::Wrapper<mu2e::TrackerHitByID>;
template class art::Wrapper<mu2e::TrackSeedCollection>;
template class art::Wrapper<mu2e::TrkExtTrajCollection>;
template class art::Wrapper<mu2e::PIDProductCollection>;
template class art::Wrapper<mu2e::AvikPIDProductCollection>;
template class art::Wrapper<mu2e::AvikPIDNewProductCollection>;
template class art::Wrapper<std::vector <mu2e::CrvCoincidenceCheckResult::CoincidenceCombination> >;
template class art::Wrapper<std::vector <mu2e::CrvRecoPulses::CrvSingleRecoPulse> >;
template class art::Wrapper<mu2e::CrvRecoPulsesCollection>;
template class art::Wrapper<mu2e::CrvCoincidenceCheckResult>;

template class std::vector<mu2e::ExtMonFNALRawHit>;
template class art::Wrapper<mu2e::ExtMonFNALRawHitCollection>;

template class art::Ptr<mu2e::ExtMonFNALRawHit>;
template class art::PtrVector<mu2e::ExtMonFNALRawHit>;

template class art::Ptr<mu2e::ExtMonFNALRawCluster>;
template class std::vector<mu2e::ExtMonFNALRawCluster>;
template class art::Wrapper<mu2e::ExtMonFNALRawClusterCollection>;

template class std::vector<mu2e::ExtMonFNALRecoCluster>;
template class std::vector<std::vector<mu2e::ExtMonFNALRecoCluster> >;
template class art::Wrapper<mu2e::ExtMonFNALRecoClusterCollection>;

template class std::vector<mu2e::ExtMonFNALTrkParam>;

template class art::Ptr<mu2e::ExtMonFNALRecoCluster>;
template class art::Ptr<mu2e::ExtMonFNALTrkParam>;
template class std::pair<art::Ptr<mu2e::ExtMonFNALRecoCluster>, art::Ptr<mu2e::ExtMonFNALTrkParam> >;
template class std::pair<art::Ptr<mu2e::ExtMonFNALTrkParam>, art::Ptr<mu2e::ExtMonFNALRecoCluster> >;

template class std::vector<mu2e::ExtMonFNALTrkClusterResiduals>;
template class std::vector<mu2e::ExtMonFNALTrkFit>;
template class art::Wrapper<mu2e::ExtMonFNALTrkFitCollection>;

//template class std::array<unsigned short,10>; // used in StrawDigi
//template class std:array<unsigned long,2>; // used in StrawDigi
template class art::Ptr<mu2e::StrawDigi>;
template class std::vector<art::Ptr<mu2e::StrawDigi> >;
template class art::Wrapper<mu2e::StrawDigiCollection>;

template class std::vector<mu2e::TrackSummary>;
template class art::Wrapper<std::vector<mu2e::TrackSummary> >;
template class art::Ptr<mu2e::TrackSummary>;

template class std::pair<art::Ptr<mu2e::KalRepPtr>,art::Ptr<mu2e::CaloCluster> >;
template class std::pair<art::Ptr<mu2e::CaloCluster>,art::Ptr<mu2e::KalRepPtr> >;
template class std::vector<mu2e::TrackCaloMatchInfo>;
template class art::Assns<mu2e::KalRepPtr,mu2e::CaloCluster>;
template class art::Assns<mu2e::CaloCluster,mu2e::KalRepPtr>;
template class art::Assns<mu2e::KalRepPtr,mu2e::CaloCluster, mu2e::TrackCaloMatchInfo>;
template class art::Assns<mu2e::CaloCluster,mu2e::KalRepPtr, mu2e::TrackCaloMatchInfo>;
template class art::Wrapper<art::Assns<mu2e::KalRepPtr,mu2e::CaloCluster, mu2e::TrackCaloMatchInfo> >;
template class art::Wrapper<art::Assns<mu2e::CaloCluster,mu2e::KalRepPtr, mu2e::TrackCaloMatchInfo> >;

// BTrk Kalman fit persistence
template class std::vector<mu2e::TrkStrawHitSeed>;
template class std::vector<mu2e::KalSegment>;
template class art::Wrapper<mu2e::KalSeedCollection>;
template class art::Ptr<mu2e::KalSeed>;
template class art::Wrapper<mu2e::TrkStrawHitSeed>;

#undef ENABLE_MU2E_GENREFLEX_HACKS
