#ifndef KalmanTests_TrkRecoTrkCollection_hh
#define KalmanTests_TrkRecoTrkCollection_hh

//
// Define a type for a collection of TrkRecoTrk objects.
//
// $Id: TrkRecoTrkCollection.hh,v 1.2 2012/07/10 19:31:44 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/07/10 19:31:44 $
//
// Original author Rob Kutschke
//


// Any class that includes this header needs the following
// using declaration before including this header.  This will
// be needed until the BaBar code is modified.
// using namespace CLHEP:

#include "BaBar/BaBar.hh"
#include "TrkBase/TrkRecoTrk.hh"
#include "GeneralUtilities/inc/OwningPointerCollection.hh"

namespace mu2e {
  typedef mu2e::OwningPointerCollection<TrkRecoTrk> TrkRecoTrkCollection;
}

#endif /* KalmanTests_TrkRecoTrkCollection_hh */
