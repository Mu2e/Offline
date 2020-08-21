#ifndef TrkExt_TrkExtTrajCollection_hh
#define TrkExt_TrkExtTrajCollection_hh

//
// Define a type for a collection of TrkExtTraj objects.
//
//
// Original author MyeongJae Lee
//

#include <vector>

#include "RecoDataProducts/inc/TrkExtTraj.hh"

namespace mu2e {
   typedef std::vector<mu2e::TrkExtTraj> TrkExtTrajCollection;
}

#endif 
