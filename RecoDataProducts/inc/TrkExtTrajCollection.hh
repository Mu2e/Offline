#ifndef TrkExt_TrkExtTrajCollection_hh
#define TrkExt_TrkExtTrajCollection_hh

//
// Define a type for a collection of TrkExtTraj objects.
//
// $Id: TrkExtTrajCollection.hh,v 1.1 2012/08/04 00:16:02 mjlee Exp $
// $Author: mjlee $
// $Date: 2012/08/04 00:16:02 $
//
// Original author MyeongJae Lee
//

#include <vector>

#include "RecoDataProducts/inc/TrkExtTraj.hh"

namespace mu2e {
   typedef std::vector<mu2e::TrkExtTraj> TrkExtTrajCollection;
}

#endif 
