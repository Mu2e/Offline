//
// Build a dictionary.
//
// $Id: classes.h,v 1.2 2011/06/11 03:17:48 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/06/11 03:17:48 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "art/Persistency/Common/Wrapper.h"

using namespace CLHEP;
#include "TrkBase/TrkRecoTrk.hh"
#include "KalmanTests/inc/TrkRecoTrkCollection.hh" 


template class art::Wrapper<mu2e::TrkRecoTrkCollection>;
