//
// Build a dictionary.
//
// $Id: classes.h,v 1.13 2010/09/29 19:37:58 logash Exp $
// $Author: logash $
// $Date: 2010/09/29 19:37:58 $
//
// Original author Rob Kutschke
//
// Notes:
// 1) The system is not able to deal with
//    edm::Wrapper<std::vector<std::string> >;
//    The problem is somewhere inside root's reflex mechanism
//    and Philippe Canal says that it is ( as of March 2010) a
//    known problem.  He also says that they do not have any
//    plans to fix it soon.  We can always work around it
//    by putting the string inside another object.

#include <vector>

#include "DataFormats/Common/interface/AssociationVector.h"
#include "DataFormats/Common/interface/OwnVector.h"
#include "DataFormats/Common/interface/SortedCollection.h"
#include "DataFormats/Common/interface/Wrapper.h"


#include "ToyDP/inc/CrudeStrawHitPData.hh"
#include "ToyDP/inc/HoughCircleCollection.hh"
#include "ToyDP/inc/PhysicalVolumeInfoCollection.hh"
#include "ToyDP/inc/SimParticleCollection.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "ToyDP/inc/StrawHitCollection.hh"
#include "ToyDP/inc/StrawHitMCTruthCollection.hh"
#include "ToyDP/inc/StrawHitMCPtrCollection.hh"
#include "ToyDP/inc/ToyGenParticleCollection.hh"
#include "ToyDP/inc/ToyHitCollection.hh"
#include "ToyDP/inc/G4BeamlineInfo.hh"
#include "ToyDP/inc/G4BeamlineInfoCollection.hh"
#include "ToyDP/inc/CaloHitCollection.hh"
#include "ToyDP/inc/CaloHitMCTruthCollection.hh"

//
// Only include objects that we would like to be able to put into the event.
// Do not include the objects they contain internally.
//

template class edm::Wrapper<mu2e::ToyHitCollection>;
template class edm::Wrapper<mu2e::ToyGenParticleCollection>;
template class edm::Wrapper<mu2e::StepPointMCCollection>;
template class edm::Wrapper<mu2e::PhysicalVolumeInfoCollection>;
template class edm::Wrapper<mu2e::CrudeStrawHitPData>;
template class edm::Wrapper<mu2e::SimParticleCollection>;
template class edm::Wrapper<mu2e::HoughCircleCollection>;
template class edm::Wrapper<mu2e::StrawHitCollection>;
template class edm::Wrapper<mu2e::StrawHitMCTruthCollection>;
template class edm::Wrapper<mu2e::StrawHitMCPtrCollection>;
template class edm::Wrapper<mu2e::G4BeamlineInfo>;
template class edm::Wrapper<mu2e::G4BeamlineInfoCollection>;
template class edm::Wrapper<mu2e::CaloHitCollection>;
template class edm::Wrapper<mu2e::CaloHitMCTruthCollection>;
