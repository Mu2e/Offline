//
//  'global' information used by all TrkHit subclasses controlling access to the fit, residuals,
//  simulated annealing parameters, etc.
//
//  Original author David Brown LBNL
//
#ifndef BtrkData_TrkHitContext_hh
#define BtrkData_TrkHitContext_hh

#include "BTrk/TrkBase/TrkHit.hh"
namespace mu2e {
  struct TrkHitContext {
    double _exterr; // error assocated with temperature in simulated annealing
    double _maxdriftpull; // parameter to define physical drift
  };
}
#endif
