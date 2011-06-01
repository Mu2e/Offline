#include "art/Persistency/Common/Wrapper.h"

#include "ITrackerGeom/inc/ITrackerWiredata.hh"

namespace {
struct dictionary {
        art::Wrapper<ITrackerWiredata>                  dummy1;
        art::Wrapper<IlcDCHwiredata>                    dummy3;
};
}
