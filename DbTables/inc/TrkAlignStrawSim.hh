#ifndef DbTables_TrkAlignStrawSim_hh
#define DbTables_TrkAlignStrawSim_hh

#include "Offline/DbTables/inc/DbTable.hh"
#include "Offline/DbTables/inc/TrkStrawEndAlign.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include <iomanip>
#include <map>
#include <sstream>
#include <string>

namespace mu2e {

class TrkAlignStrawSim : public TrkAlignStraw {
 public:
  constexpr static const char* cxname = "TrkAlignStrawSim";
  TrkAlignStrawSim() : TrkAlignStraw(cxname,"trk.alignstrawsim") {}
};

}  // namespace mu2e
#endif
