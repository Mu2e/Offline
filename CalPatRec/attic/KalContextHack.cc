// Fixme:
// Temporary hack needed to make closed links work.
//
// The underlying issue is that BTrk/KalmanTrack/KalContext.hh
// is not copyable; therefore it has two functions that are declared
// but unimplemented: the copy c'tor and copy assignment.
//
// This works fine until we try to make a dictionary with a
// closed link. The hack is to provide implementations that throw
// if they are ever called.
//

#include "BTrk/KalmanTrack/KalContext.hh"
#include "cetlib_except/exception.h"

KalContext::KalContext(const KalContext& ){
  throw cet::exception("RECO")
    <<"mu2e::KalFitHack: calling copy c'tor of KalContex\nl";
}

KalContext& KalContext::operator = (const KalContext& ){
  throw cet::exception("RECO")
    <<"mu2e::KalFitHack: calling copy assignment of KalContext\nl";
}
