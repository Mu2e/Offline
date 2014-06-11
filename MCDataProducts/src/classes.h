//
// Build a dictionary.
//
// $Id: classes.h,v 1.33 2014/06/11 00:20:27 gandr Exp $
// $Author: gandr $
// $Date: 2014/06/11 00:20:27 $
//
// Original author Rob Kutschke
//

//
// For every type that is described in classes_def.xml, the header must
// be included in this file. It is OK if the header comes in as the
// result of including another header.
//
// For every type described in classes_def.xml that is a templated type
// there must be an instantiation of that template in this file.
// Again, it is OK if a template is instantiated as the result of
// instantiating another template.
//
#include <vector>
#include <map>

#include "art/Persistency/Common/Wrapper.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/Assns.h"
#include "art/Persistency/Common/RNGsnapshot.h"
#include "cetlib/map_vector.h"

#include "MCDataProducts/inc/CaloCrystalOnlyHitCollection.hh"
#include "MCDataProducts/inc/CaloHitMCTruthCollection.hh"
#include "MCDataProducts/inc/CaloHitSimPartMCCollection.hh"
#include "MCDataProducts/inc/ExtMonUCITofHitMCTruthCollection.hh"
#include "MCDataProducts/inc/G4BeamlineInfo.hh"
#include "MCDataProducts/inc/G4BeamlineInfoCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfo.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoMultiCollection.hh"
#include "MCDataProducts/inc/MCTrajectoryCollection.hh"
#include "MCDataProducts/inc/PointTrajectoryCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/SimParticlePtrCollection.hh"
#include "MCDataProducts/inc/StatusG4.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/MixingSummary.hh"
#include "MCDataProducts/inc/VisibleGenElTrackCollection.hh"
#include "MCDataProducts/inc/GenParticleSPMHistory.hh"
#include "MCDataProducts/inc/GenSimParticleLink.hh"
#include "MCDataProducts/inc/ExtMonFNALSimHit.hh"
#include "MCDataProducts/inc/ExtMonFNALSimHitCollection.hh"
#include "MCDataProducts/inc/ExtMonFNALHitTruthAssn.hh"
#include "MCDataProducts/inc/ExtMonFNALRecoClusterTruthAssn.hh"
#include "MCDataProducts/inc/ExtMonFNALPatRecTruthAssns.hh"
#include "MCDataProducts/inc/MARSInfo.hh"
#include "MCDataProducts/inc/MARSInfoCollection.hh"
#include "MCDataProducts/inc/GenParticleMARSAssns.hh"
#include "MCDataProducts/inc/SimParticleMARSAssns.hh"
#include "MCDataProducts/inc/SimParticleTimeMap.hh"

#include "MCDataProducts/inc/StepFilterMode.hh"
#include "MCDataProducts/inc/GenEventCount.hh"
#include "MCDataProducts/inc/EventWeight.hh"

#include "MCDataProducts/inc/TrackSummaryTruthAssns.hh"

// For cet::map_vector<T> instantiate the component pair<> and vector<pair<>> templates.
template class std::pair<cet::map_vector_key,mu2e::SimParticle>;
template class std::pair<cet::map_vector_key,mu2e::PointTrajectory>;
template class std::pair<art::Ptr<mu2e::SimParticle>, mu2e::MCTrajectory>;
template class std::vector<std::pair<cet::map_vector_key,mu2e::SimParticle> >;
template class std::vector<std::pair<cet::map_vector_key,mu2e::PointTrajectory> >;
template class std::map<art::Ptr<mu2e::SimParticle>, mu2e::MCTrajectory>;

// Cannot use the typedefs in here - not sure why.
template class art::Ptr<mu2e::StepPointMC>;
template class std::vector<art::Ptr<mu2e::StepPointMC> >;
template class std::vector<std::vector<art::Ptr<mu2e::StepPointMC> > >;
template class std::pair<CLHEP::Hep3Vector,CLHEP::HepLorentzVector>;
template class std::map<art::Ptr<mu2e::SimParticle>::key_type,mu2e::GenElHitData>;// GenElHitDataCollection;
template class std::vector<mu2e::ExtMonFNALSimHit>;

template class art::Wrapper<mu2e::GenParticleCollection>;
template class art::Wrapper<mu2e::StepPointMCCollection>;
template class art::Wrapper<mu2e::PhysicalVolumeInfoCollection>;
template class art::Wrapper<mu2e::PhysicalVolumeInfoMultiCollection>;
template class art::Wrapper<mu2e::SimParticleCollection>;
template class art::Wrapper<mu2e::StrawHitMCTruthCollection>;
template class art::Wrapper<mu2e::StrawDigiMCCollection>;
template class art::Wrapper<mu2e::G4BeamlineInfo>;
template class art::Wrapper<mu2e::G4BeamlineInfoCollection>;
template class art::Wrapper<mu2e::CaloHitMCTruthCollection>;
template class art::Wrapper<mu2e::CaloHitSimPartMCCollection>;
template class art::Wrapper<mu2e::CaloCrystalOnlyHitCollection>;
template class art::Wrapper<mu2e::ExtMonUCITofHitMCTruthCollection>;
template class art::Wrapper<mu2e::PointTrajectoryCollection>;
template class art::Wrapper<mu2e::MCTrajectoryCollection>;
template class art::Wrapper<mu2e::StatusG4>;
template class art::Wrapper<mu2e::PtrStepPointMCVectorCollection>;
template class art::Wrapper<mu2e::MixingSummary>;
template class art::Wrapper<std::vector<art::RNGsnapshot> >;
template class art::Wrapper<mu2e::ExtMonFNALSimHitCollection>;

template class art::Wrapper<mu2e::GenElHitData>;
template class art::Wrapper<mu2e::VisibleGenElTrackCollection>;

// A way to instantiate a typedef without keeping sync with its definition
namespace {
  struct Instantiations {
    mu2e::GenParticleSPMHistory gpspmh;
    mu2e::GenSimParticleLink gspl;
    mu2e::GenParticleMARSAssns gpwa;
    mu2e::SimParticleMARSAssns spwa;
    mu2e::SimParticleTimeMap sppptm;
  };
}

template class std::pair<cet::map_vector_key,mu2e::PhysicalVolumeInfo>;
template class std::vector<std::pair<cet::map_vector_key,mu2e::PhysicalVolumeInfo> >;
template class std::pair<unsigned int,cet::map_vector<mu2e::PhysicalVolumeInfo> >;
template class std::vector<std::pair<unsigned int,cet::map_vector<mu2e::PhysicalVolumeInfo> > >;

template class art::Wrapper<mu2e::GenParticleSPMHistory>;
template class art::Wrapper<mu2e::GenSimParticleLink>;

template class std::pair<art::Ptr<mu2e::GenParticle>, art::Ptr<mu2e::MARSInfo> >;
template class std::pair<art::Ptr<mu2e::MARSInfo>, art::Ptr<mu2e::GenParticle> >;
template class art::Wrapper<mu2e::MARSInfoCollection>;
template class art::Wrapper<mu2e::GenParticleMARSAssns>;

template class std::pair<art::Ptr<mu2e::SimParticle>, art::Ptr<mu2e::MARSInfo> >;
template class std::pair<art::Ptr<mu2e::MARSInfo>, art::Ptr<mu2e::SimParticle> >;
template class art::Wrapper<mu2e::SimParticleMARSAssns>;

template class art::Wrapper<mu2e::SimParticleTimeMap>;
template class std::pair<art::Ptr<mu2e::SimParticle>,double>;

template class std::vector<mu2e::ExtMonFNALHitTruthBits>;
template class std::pair<art::Ptr<mu2e::SimParticle>,art::Ptr<mu2e::ExtMonFNALRawHit> >;
template class std::pair<art::Ptr<mu2e::ExtMonFNALRawHit>,art::Ptr<mu2e::SimParticle> >;
template class art::Wrapper<art::Assns<mu2e::SimParticle,mu2e::ExtMonFNALRawHit,mu2e::ExtMonFNALHitTruthBits> >;
template class art::Wrapper<art::Assns<mu2e::ExtMonFNALRawHit,mu2e::SimParticle,mu2e::ExtMonFNALHitTruthBits> >;

template class std::vector<mu2e::ExtMonFNALRecoClusterTruthBits>;
template class std::pair<art::Ptr<mu2e::SimParticle>,art::Ptr<mu2e::ExtMonFNALRecoCluster> >;
template class std::pair<art::Ptr<mu2e::ExtMonFNALRecoCluster>,art::Ptr<mu2e::SimParticle> >;
template class art::Wrapper<art::Assns<mu2e::SimParticle,mu2e::ExtMonFNALRecoCluster,mu2e::ExtMonFNALRecoClusterTruthBits> >;
template class art::Wrapper<art::Assns<mu2e::ExtMonFNALRecoCluster,mu2e::SimParticle,mu2e::ExtMonFNALRecoClusterTruthBits> >;

template class std::vector<mu2e::ExtMonFNALTrkMatchInfo>;
template class std::pair<art::Ptr<mu2e::SimParticle>,art::Ptr<mu2e::ExtMonFNALTrkFit> >;
template class std::pair<art::Ptr<mu2e::ExtMonFNALTrkFit>,art::Ptr<mu2e::SimParticle> >;
template class art::Wrapper<art::Assns<mu2e::SimParticle,mu2e::ExtMonFNALTrkFit,mu2e::ExtMonFNALTrkMatchInfo> >;
template class art::Wrapper<art::Assns<mu2e::ExtMonFNALTrkFit,mu2e::SimParticle,mu2e::ExtMonFNALTrkMatchInfo> >;

template class std::vector<art::Ptr<mu2e::SimParticle> >;
template class art::Wrapper<std::vector<art::Ptr<mu2e::SimParticle> > >;

template class art::Wrapper<mu2e::GenEventCount>;
template class art::Wrapper<mu2e::EventWeight>;

template class std::vector<mu2e::TrackSummaryMatchInfo>;
template class std::pair<art::Ptr<mu2e::SimParticle>,art::Ptr<mu2e::TrackSummary> >;
template class std::pair<art::Ptr<mu2e::TrackSummary>,art::Ptr<mu2e::SimParticle> >;
template class art::Assns<mu2e::SimParticle,mu2e::TrackSummary,void>;
template class art::Assns<mu2e::TrackSummary,mu2e::SimParticle,void>;
template class art::Assns<mu2e::SimParticle,mu2e::TrackSummary,mu2e::TrackSummaryMatchInfo>;
template class art::Assns<mu2e::TrackSummary,mu2e::SimParticle,mu2e::TrackSummaryMatchInfo>;
template class art::Wrapper<art::Assns<mu2e::TrackSummary,mu2e::SimParticle,mu2e::TrackSummaryMatchInfo> >;
template class art::Wrapper<art::Assns<mu2e::SimParticle,mu2e::TrackSummary,mu2e::TrackSummaryMatchInfo> >;
