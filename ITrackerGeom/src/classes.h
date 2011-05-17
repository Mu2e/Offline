#include "art/Persistency/Common/SortedCollection.h"
#include "art/Persistency/Common/OwnVector.h"
#include "art/Persistency/Common/AssociationVector.h"
#include "art/Persistency/Common/Wrapper.h"

#include "ITrackerGeom/inc/ITrackerWiredata.hh"

namespace {
struct dictionary {
	ITrackerWiredata		dummy0;
	art::Wrapper<ITrackerWiredata>			dummy1;
	IlcDCHwiredata			dummy2;
	art::Wrapper<IlcDCHwiredata>			dummy3;
};
}
