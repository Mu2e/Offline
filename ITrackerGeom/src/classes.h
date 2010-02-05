#include "DataFormats/Common/interface/SortedCollection.h"
#include "DataFormats/Common/interface/OwnVector.h"
#include "DataFormats/Common/interface/AssociationVector.h"
#include "DataFormats/Common/interface/Wrapper.h"

#include "ITrackerGeom/inc/ITrackerWiredata.hh"

namespace {
struct dictionary {
	ITrackerWiredata		dummy0;
	edm::Wrapper<ITrackerWiredata>			dummy1;
	IlcDCHwiredata			dummy2;
	edm::Wrapper<IlcDCHwiredata>			dummy3;
};
}
