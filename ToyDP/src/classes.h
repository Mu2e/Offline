//
// Build a dictionary.
//
// $Id: classes.h,v 1.25 2011/05/18 19:46:19 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/18 19:46:19 $
//
// Original author Rob Kutschke
//
// Notes:
// 1) The system is not able to deal with
//    art::Wrapper<std::vector<std::string> >;
//    The problem is somewhere inside root's reflex mechanism
//    and Philippe Canal says that it is ( as of March 2010) a
//    known problem.  He also says that they do not have any
//    plans to fix it soon.  We can always work around it
//    by putting the string inside another object.

#include <vector>

#include "art/Persistency/Common/Wrapper.h"


#include "ToyDP/inc/CrudeStrawHitPData.hh"
#include "ToyDP/inc/HoughCircleCollection.hh"
#include "ToyDP/inc/PhysicalVolumeInfoCollection.hh"
#include "ToyDP/inc/SimParticleCollection.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "ToyDP/inc/StrawHitCollection.hh"
#include "ToyDP/inc/StrawHitMCTruthCollection.hh"
#include "ToyDP/inc/StrawClusterCollection.hh"
#include "ToyDP/inc/DPIndexVectorCollection.hh"
#include "ToyDP/inc/ToyGenParticleCollection.hh"
#include "ToyDP/inc/G4BeamlineInfo.hh"
#include "ToyDP/inc/G4BeamlineInfoCollection.hh"
#include "ToyDP/inc/CaloHitCollection.hh"
#include "ToyDP/inc/CaloHitMCTruthCollection.hh"
#include "ToyDP/inc/CaloCrystalHitCollection.hh"
#include "ToyDP/inc/CaloCrystalOnlyHitCollection.hh"
#include "ToyDP/inc/PointTrajectoryCollection.hh"
#include "ToyDP/inc/StatusG4.hh"

//
// I am not 100% clear what needs to be here.  I do know:
//
// 1) There must be a line for every Wrapper that appears
//    in classes_def.xml
// 2) For any map that appears in classes_def.xml, there
//    must be two lines in this file: one for the map
//    and one for the underlying pair type.
//

template class std::pair<MapVectorKey,mu2e::SimParticle>;
template class std::map<MapVectorKey,mu2e::SimParticle>;
template class std::pair<MapVectorKey,mu2e::PointTrajectory>;
template class std::map<MapVectorKey,mu2e::PointTrajectory>;
template class std::vector<uint32_t>;

template class art::Wrapper<mu2e::ToyGenParticleCollection>;
template class art::Wrapper<mu2e::StepPointMCCollection>;
template class art::Wrapper<mu2e::PhysicalVolumeInfoCollection>;
template class art::Wrapper<mu2e::CrudeStrawHitPData>;
template class art::Wrapper<mu2e::SimParticleCollection>;
template class art::Wrapper<mu2e::HoughCircleCollection>;
template class art::Wrapper<mu2e::StrawHitCollection>;
template class art::Wrapper<mu2e::StrawClusterCollection>;
template class art::Wrapper<mu2e::StrawHitMCTruthCollection>;
template class art::Wrapper<mu2e::G4BeamlineInfo>;
template class art::Wrapper<mu2e::G4BeamlineInfoCollection>;
template class art::Wrapper<mu2e::CaloHitCollection>;
template class art::Wrapper<mu2e::CaloHitMCTruthCollection>;
template class art::Wrapper<mu2e::CaloCrystalHitCollection>;
template class art::Wrapper<mu2e::CaloCrystalOnlyHitCollection>;
template class art::Wrapper<mu2e::PointTrajectoryCollection>;
template class art::Wrapper<mu2e::DPIndexVectorCollection>;
template class art::Wrapper<mu2e::StatusG4>;
