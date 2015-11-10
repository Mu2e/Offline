#ifndef HitMakers_formStepPointMCStrawHit_hh
#define HitMakers_formStepPointMCStrawHit_hh
//
// Free function to create the hall walls and hall interior inside the earthen overburden.
//
// $Id: formStepPointMCStrawHit.hh,v 1.2 2013/03/15 15:52:04 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/15 15:52:04 $
//
// Original author KLG
//

// Mu2e includes.

#include "ConditionsService/inc/ConditionsHandle.hh"
#include "GlobalConstantsService/inc/MassCache.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
#include "DataProducts/inc/StrawIndex.hh"
#include "MCDataProducts/inc/StepPointMCStrawHit.hh"
#include "TrackerGeom/inc/Tracker.hh"

// Art includes.

#include "art/Persistency/Common/Ptr.h"

// Other includes.

#include "CLHEP/Random/RandGaussQ.h"

// C++ includes.

#include <memory>

namespace mu2e {

  std::unique_ptr<StepPointMCStrawHit> formStepPointMCStrawHit(
                     art::Ptr<StepPointMC> const spmcp,
                     StrawIndex const & straw_id,
                     double _minimumLength,
                     bool   _enableFlightTimeCorrection,
                     MassCache & cache,
                     CLHEP::RandGaussQ & _gaussian,
                     Tracker const & tracker,
                     ConditionsHandle<TrackerCalibrations> const & trackerCalibrations
                     );
}

#endif /* HitMakers_formStepPointMCStrawHit_hh */
