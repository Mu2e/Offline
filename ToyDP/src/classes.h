//
// Build a dictionary.
//
// $Id: classes.h,v 1.7 2010/03/23 20:37:00 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/03/23 20:37:00 $
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

#include "DataFormats/Common/interface/SortedCollection.h"
#include "DataFormats/Common/interface/OwnVector.h"
#include "DataFormats/Common/interface/AssociationVector.h"
#include "DataFormats/Common/interface/Wrapper.h"


#include "ToyDP/inc/ToyHitCollection.hh"
#include "ToyDP/inc/ToyGenParticleCollection.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "ToyDP/inc/CrudeStrawHitPData.hh"
#include "ToyDP/inc/RandomEngineState.hh"
#include "ToyDP/inc/PhysicalVolumeInfoCollection.hh"
#include "ToyDP/inc/SimParticleCollection.hh"

//
// Only include objects that we would like to be able to put into the event.
// Do not include the objects they contain internally.
//

template class edm::Wrapper<mu2e::ToyHitCollection>;
template class edm::Wrapper<mu2e::ToyGenParticleCollection>;
template class edm::Wrapper<mu2e::StepPointMCCollection>;
template class edm::Wrapper<mu2e::PhysicalVolumeInfoCollection>;
template class edm::Wrapper<mu2e::CrudeStrawHitPData>;
template class edm::Wrapper<std::vector<mu2e::RandomEngineState> >;
template class edm::Wrapper<mu2e::SimParticleCollection>;

