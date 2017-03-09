//
//  Base clase for physics interface to track fit event objects
//  Dave Brown (LBNL) 2/12/17
//
#ifndef TrkAdapter_TrkInter_hh
#define TrkAdapter_TrkInter_hh
#include "DataProducts/inc/PDGCode.hh"
// data
#include "RecoDataProducts/inc/TrkFitFlag.hh"
// CLHEP
#include "CLHEP/Vector/ThreeVector.h"
// BTrk
#include "BTrk/TrkBase/TrkParticle.hh"

namespace mu2e {
  class TrkInter {
    TrkInter() {};
    virtual ~TrkInter() = 0;
// basic description of the particle interpretation of this track
// particle type used in the physics interpretation.  For now, include
// both sim and reco interfaces, we should standardize FIXME!
    PDGCode::type pdgId();
    TrkParticle particle();

// translate some standard positions into times.  The returned vector will
// give the time-ordered times where the condition is satisfied
    void crossesZ(std::vector<double>& times,double z);

// position wrt Tracker center of the particle as a function of time in the microbunch (mm)
    virtual CLHEP::Hep3Vector position(double time) const = 0;
// direction (unit vector) of the particle as a function of time in the microbunch
    virtual CLHEP::Hep3Vector direction(double time) const = 0;
// momentum of the particle as a function of time in the microbunch (MeV/c)
    virtual CLHEP::Hep3Vector momentum(double time) const = 0;
// information about the fit itself.  First, global information
// fit status
    virtual TrkFitFlag status() const = 0;
// fit consistency value
    virtual double fitcon() const = 0;
// it's not clear every fit will have a valid 'chisquared', but for now keep these FIXME!
    virtual double chisquared() const = 0;
    virtual unsigned nDOF() const = 0;

  };
}
#endif
