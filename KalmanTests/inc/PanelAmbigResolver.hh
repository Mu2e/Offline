//
// class to resolve hit ambiguities by panel, assuming a reasonable track
// fit as input
//
// $Id: PanelAmbigResolver.hh,v 1.5 2014/08/01 18:56:10 gandr Exp $
// $Author: gandr $ 
// $Date: 2014/08/01 18:56:10 $
//
#ifndef PanelAmbigResolver_HH
#define PanelAmbigResolver_HH
#include "BaBar/BaBar.hh"
#include "KalmanTests/inc/TrkStrawHitState.hh"
#include "KalmanTests/inc/AmbigResolver.hh"
#include "TrackerGeom/inc/LayerId.hh"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Vector/ThreeVector.h"
#ifndef __GCCXML__
#include "fhiclcpp/ParameterSet.h"
#endif/*__GCCXML__*/
#include <cstddef>
#include <vector>

class KalRep;
class TrkT0;

namespace mu2e {

// class to store the 1-d projection of a hit.
// this direction (u) is perpendicular to the wire and to the track as it hits the panel
  struct TSHUInfo {
    // construct from a TrkStrawHit, the U direction and the U origin in mu2e coordinates
    TSHUInfo(const TrkStrawHit* tsh,CLHEP::Hep3Vector const& udir, CLHEP::Hep3Vector const& uorigin);
    const TrkStrawHit* _tsh; // reference the straw hit, to allow consistency checking
    double _upos; // hit position in u relative to the first 
    double _uwt; // cache of the hit weight = 1/error^2
  };

// struct to hold hit and other information for a single panel.
// It also holds the result of the optimization
  struct PanelInfo {
    PanelInfo() : _utpos(0.0), _utwt(0.0) {}
// information from the track fit, excluding the hits from panel, projected into this panel's measurement direction
    double _utpos;
    double _utwt; // weight = 1/error^2
// information from all the hits
    std::vector<TSHUInfo> _uinfo;
  };

// struct to hold the result of optimizing the panel parameters.
  struct PanelResult {
    PanelResult(TrkStrawHitStateVector const& state) : _state(state),
    _chisq(0.0), _status(-1) {}
    TrkStrawHitStateVector _state; // ambiguity/activity state of each hit in this panel
    double _chisq; // chisquared of this state, including all penalties
    CLHEP::HepVector _delta; // change in parameters at optimum; U is parameter 0, deltaT is parameter 1
    CLHEP::HepSymMatrix _dcov; // covariance of parameter changes
    int _status; // status of matrix inversion
  };

  class PanelAmbigResolver : public AmbigResolver {
    public:
// construct from parameter set
#ifndef __GCCXML__
      explicit PanelAmbigResolver(fhicl::ParameterSet const&);
#endif/*__GCCXML__*/
      virtual ~PanelAmbigResolver();
// resolve a track.  Depending on the configuration, this might
// update the hit state and the t0 value.
      virtual void resolveTrk(KalFitResult& kfit) const;
    private:
// resolve the ambiguity on a single panel
      void resolvePanel(TSHV& phits, KalFitResult& kfit) const;
// fill information about a given panel's track and hits
      bool fillPanelInfo(TSHV const& phits, const KalRep* krep, PanelInfo& pinfo) const;
// compute the panel result for a given ambiguity/activity state and the ionput t0
      void fillResult(PanelInfo const& pinfo,TrkT0 const& t0, PanelResult& result) const;
// parameters
      double _minsep; // minimum chisquared separation between best solution and the rest to consider a panel resolved
      double _inactivepenalty; // chisquared penalty for an inactive hit
      double _penaltyres; // resolution term to add to hits if ambiguity/activity can't be resolved
      double _nullerr2; // additional error (squared) for hits with null ambiguity
      TSHSSV _allowed; // allowed states of a TrkStrawHit
  };
}

#endif

