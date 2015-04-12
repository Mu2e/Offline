///////////////////////////////////////////////////////////////////////////////
// 2015-04-08 P.Murat: 
// -------------------
// resolve hit ambiguities by panel, assuming a reasonable track fit as input
///////////////////////////////////////////////////////////////////////////////

#ifndef DoubletAmbigResolver_HH
#define DoubletAmbigResolver_HH

#include "BaBar/BaBar.hh"
#include "KalmanTests/inc/AmbigResolver.hh"

#ifndef __GCCXML__
#include "fhiclcpp/ParameterSet.h"
#endif/*__GCCXML__*/
#include <cstddef>
#include <vector>

class TrkDifTraj;
class KalRep;

#include "KalmanTests/inc/Doublet.hh"

namespace mu2e {

  class DoubletAmbigResolver : public AmbigResolver {
    public:
      enum trajtype {reftraj=0};

  protected:
    int    _iherr;
    int    _debugLevel;

    double _mindrift;                   // minimum drift to assign an ambiguity.  Below this, an ambiguity of '0' is defined
    double _zeropenalty;                // special penalty for drifts below the minimum
    bool   _penalty;			// apply penalty or notA
                                        // exponential + linear fit to ambiguity mis-assignment
    double _expnorm;
    double _lambda;
    double _offset;
    double _slope;
					// need to define !!!!
    double _sigmaSlope;
    double _maxDoubletChi2;
    double _scaleErrDoublet;
    double _minDriftDoublet;
    double _deltaDriftDoublet;

    double _decisionMode;
    int    _sign[4][2];
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
  public:
#ifndef __GCCXML__
    explicit DoubletAmbigResolver(fhicl::ParameterSet const&, int Iherr);
#endif/*__GCCXML__*/
    virtual ~DoubletAmbigResolver();

    void findLines       (Hep3Vector* Pos, double* R, double* Slopes) const ;
    void findDoublets    (KalFitResult& KRes) const ;

    void markDoublet     (KalFitResult& KRes, Doublet *doublet, int index0, int index1) const;
    void markMultiplets  (KalFitResult& Kres) const;
    void resolveSingleHit(KalFitResult& Kres, mu2e::TrkStrawHit* Hit) const ;
    
					// resolve a track.  Depending on the configuration, this might
					// penalty function depends on the drift radius

    double penaltyError(double rdrift) const;
					// update the hit state and the t0 value.
    virtual void resolveTrk(KalFitResult& kfit) const;
  };
}
#endif
