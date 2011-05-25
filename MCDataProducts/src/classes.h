//
// Build a dictionary.
//
// $Id: classes.h,v 1.3 2011/05/25 01:23:14 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/05/25 01:23:14 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "art/Persistency/Common/Wrapper.h"
#include "cetlib/map_vector.h"

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

template class std::pair<cet::map_vector_key,mu2e::SimParticle>;
template class std::pair<cet::map_vector_key,mu2e::PointTrajectory>;

template class std::vector<std::pair<cet::map_vector_key,mu2e::SimParticle> >;
template class std::vector<std::pair<cet::map_vector_key,mu2e::PointTrajectory> >;

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
