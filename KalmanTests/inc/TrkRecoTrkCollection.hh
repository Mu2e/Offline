#ifndef KalmanTests_TrkRecoTrkCollection_hh
#define KalmanTests_TrkRecoTrkCollection_hh

//
// Define a type for a collection of TrkRecoTrk objects.
//
// $Id: TrkRecoTrkCollection.hh,v 1.1 2011/06/11 03:17:48 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/06/11 03:17:48 $
//
// Original author Rob Kutschke
//


// Any class that includes this header needs the following
// using declaration before including this header.  This will
// be needed until the BaBar code is modified.
// using namespace CLHEP:

#include "TrkBase/TrkRecoTrk.hh"
#include "GeneralUtilities/inc/OwningPointerCollection.hh"

namespace mu2e {
  typedef mu2e::OwningPointerCollection<TrkRecoTrk> TrkRecoTrkCollection;
}

#endif /* KalmanTests_TrkRecoTrkCollection_hh */
