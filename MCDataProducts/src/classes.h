//
// Build a dictionary.
//
// $Id: classes.h,v 1.8 2011/12/30 20:31:46 youzy Exp $
// $Author: youzy $
// $Date: 2011/12/30 20:31:46 $
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

#include "art/Persistency/Common/Wrapper.h"
#include "cetlib/map_vector.h"

#include "MCDataProducts/inc/CaloCrystalOnlyHitCollection.hh"
#include "MCDataProducts/inc/CaloHitMCTruthCollection.hh"
#include "MCDataProducts/inc/ExtMonUCITofHitMCTruthCollection.hh"
#include "MCDataProducts/inc/G4BeamlineInfo.hh"
#include "MCDataProducts/inc/G4BeamlineInfoCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/PointTrajectoryCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StatusG4.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/MixingSummary.hh"
#include "art/Persistency/Common/RNGsnapshot.h"
#include "MCDataProducts/inc/VisibleGenElTrackCollection.hh"

// For cet::map_vector<T> instantiate the component pair<> and vector<pair<>> templates.
template class std::pair<cet::map_vector_key,mu2e::SimParticle>;
template class std::pair<cet::map_vector_key,mu2e::PointTrajectory>;

template class std::vector<std::pair<cet::map_vector_key,mu2e::SimParticle> >;
template class std::vector<std::pair<cet::map_vector_key,mu2e::PointTrajectory> >;

// Cannot use the typedefs in here - not sure why.
template class art::Ptr<mu2e::StepPointMC>;
template class std::vector<art::Ptr<mu2e::StepPointMC> >;
template class std::vector<std::vector<art::Ptr<mu2e::StepPointMC> > >;

template class art::Wrapper<mu2e::GenParticleCollection>;
template class art::Wrapper<mu2e::StepPointMCCollection>;
template class art::Wrapper<mu2e::PhysicalVolumeInfoCollection>;
template class art::Wrapper<mu2e::SimParticleCollection>;
template class art::Wrapper<mu2e::StrawHitMCTruthCollection>;
template class art::Wrapper<mu2e::G4BeamlineInfo>;
template class art::Wrapper<mu2e::G4BeamlineInfoCollection>;
template class art::Wrapper<mu2e::CaloHitMCTruthCollection>;
template class art::Wrapper<mu2e::CaloCrystalOnlyHitCollection>;
template class art::Wrapper<mu2e::ExtMonUCITofHitMCTruthCollection>;
template class art::Wrapper<mu2e::PointTrajectoryCollection>;
template class art::Wrapper<mu2e::StatusG4>;
template class art::Wrapper<mu2e::PtrStepPointMCVectorCollection>;
template class art::Wrapper<mu2e::MixingSummary>;
template class art::Wrapper<std::vector<art::RNGsnapshot> >;
template class art::Wrapper<mu2e::VisibleGenElTrackCollection>;

