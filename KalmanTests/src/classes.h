#include "DataFormats/Common/interface/SortedCollection.h"
#include "DataFormats/Common/interface/OwnVector.h"
#include "DataFormats/Common/interface/AssociationVector.h"
#include "DataFormats/Common/interface/Wrapper.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"
using namespace CLHEP;
#include "TrkBase/TrkExchangePar.hh"
#include <vector>
template class edm::Wrapper<CLHEP::HepVector>;
//template class edm::Wrapper<CLHEP::Hep3Vector>;
template class edm::Wrapper<CLHEP::HepMatrix>;

//template class edm::Wrapper<TrkExchangePar>;

