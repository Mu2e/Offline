//
// Build a dictionary.
//
// $Id: classes.h,v 1.1 2011/05/24 17:16:44 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/05/24 17:16:44 $
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

#include "MCDataProducts/inc/CaloCrystalOnlyHitCollection.hh"
#include "MCDataProducts/inc/CaloHitMCTruthCollection.hh"
#include "MCDataProducts/inc/G4BeamlineInfo.hh"
#include "MCDataProducts/inc/G4BeamlineInfoCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/PointTrajectoryCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StatusG4.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"

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

template class art::Wrapper<mu2e::GenParticleCollection>;
template class art::Wrapper<mu2e::StepPointMCCollection>;
template class art::Wrapper<mu2e::PhysicalVolumeInfoCollection>;
template class art::Wrapper<mu2e::SimParticleCollection>;
template class art::Wrapper<mu2e::StrawHitMCTruthCollection>;
template class art::Wrapper<mu2e::G4BeamlineInfo>;
template class art::Wrapper<mu2e::G4BeamlineInfoCollection>;
template class art::Wrapper<mu2e::CaloHitMCTruthCollection>;
template class art::Wrapper<mu2e::CaloCrystalOnlyHitCollection>;
template class art::Wrapper<mu2e::PointTrajectoryCollection>;
template class art::Wrapper<mu2e::StatusG4>;
