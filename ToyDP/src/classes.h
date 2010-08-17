//
// Build a dictionary.
//
// $Id: classes.h,v 1.10 2010/08/17 15:18:39 wb Exp $
// $Author: wb $
// $Date: 2010/08/17 15:18:39 $
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
#include "ToyDP/inc/ToyGenParticleCollection.hh"
#include "ToyDP/inc/ToyHitCollection.hh"

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
