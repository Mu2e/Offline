//
// Build a dictionary.
//
// $Id: classes.h,v 1.2 2009/10/06 23:20:00 kutschke Exp $
// $Author: kutschke $
// $Date: 2009/10/06 23:20:00 $
// 
// Original author Rob Kutschke

#include "DataFormats/Common/interface/SortedCollection.h"
#include "DataFormats/Common/interface/OwnVector.h"
#include "DataFormats/Common/interface/AssociationVector.h"
#include "DataFormats/Common/interface/Wrapper.h"


#include "ToyDP/inc/ToyHitCollection.hh"
#include "ToyDP/inc/ToyGenParticleCollection.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"

//
// Only include objects that we would like to be able to put into the event.
//

namespace {
struct dictionary {
  edm::Wrapper<mu2e::ToyHitCollection>         dummy301;
  edm::Wrapper<mu2e::ToyGenParticleCollection> dummy302;
  edm::Wrapper<mu2e::StepPointMCCollection>    dummy303;
};
}
