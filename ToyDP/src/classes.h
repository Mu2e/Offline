//
// Build a dictionary.
//
// $Id: classes.h,v 1.5 2010/03/05 16:07:38 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/03/05 16:07:38 $
// 
// Original author Rob Kutschke

#include <vector>
#include <string>

#include "DataFormats/Common/interface/SortedCollection.h"
#include "DataFormats/Common/interface/OwnVector.h"
#include "DataFormats/Common/interface/AssociationVector.h"
#include "DataFormats/Common/interface/Wrapper.h"


#include "ToyDP/inc/ToyHitCollection.hh"
#include "ToyDP/inc/ToyGenParticleCollection.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "ToyDP/inc/CrudeStrawHitPData.hh"

//
// Only include objects that we would like to be able to put into the event.
// Do not include the objects they contain internally.
//

//template <> class vector<std::vector<uint32_t> >;

namespace {
struct dictionary {
  edm::Wrapper<mu2e::ToyHitCollection>               dummy301;
  edm::Wrapper<mu2e::ToyGenParticleCollection>       dummy302;
  edm::Wrapper<mu2e::StepPointMCCollection>          dummy303;
  edm::Wrapper<mu2e::CrudeStrawHitPData>             dummy304;
  edm::Wrapper<std::vector<std::string> >            dummy305;

};
}
