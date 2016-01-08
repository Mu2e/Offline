///////////////////////////////////////////////////////////////////////////////
// 2015-04-08 P.Murat: 
// -------------------
// resolve hit ambiguities by panel, assuming a reasonable track fit as input
///////////////////////////////////////////////////////////////////////////////

#ifndef DoubletAmbigResolver_HH
#define DoubletAmbigResolver_HH

#include "BTrk/BaBar/BaBar.hh"
#include "TrkReco/inc/AmbigResolver.hh"

#ifndef __GCCXML__
#include "fhiclcpp/ParameterSet.h"
#endif/*__GCCXML__*/
#include <cstddef>
#include <vector>
#include "CLHEP/Vector/ThreeVector.h"

class TrkDifTraj;
class KalRep;

#include "RecoDataProducts/inc/Doublet.hh"

namespace mu2e {

  class Straw;

  class DoubletAmbigResolver : public AmbigResolver {
    public:
      enum trajtype {reftraj=0};

    struct Data_t {
      int               index[2];
      int               ibest;
      int               inext;
      
      double            chi2min;
      double            chi2next;

      const mu2e::Straw *straw[2];
      double            rdrift[2];
      double            doca[4][2];
      double            chi2[4];
      double            trkslope;
      double            lineSlopes[4];
      CLHEP::Hep3Vector spos [2];
      CLHEP::Hep3Vector sposr[2];
      CLHEP::Hep3Vector tpos [2];
      CLHEP::Hep3Vector tposr[2];
    };

  protected:
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
    int    _excludeBothHits;            // when calculating residuals to choose the drift signs 

    int    _sign[4][2];
    int    _iter; // iteration
    int	   _Final; // final iteration
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
  public:
#ifndef __GCCXML__
    explicit DoubletAmbigResolver(fhicl::ParameterSet const&, double ExtErr, int Iter, int Final);
#endif/*__GCCXML__*/
    virtual ~DoubletAmbigResolver();
    
    int calculateDoubletParameters(const KalRep* KRep, Doublet* HitDoublet, Data_t* R) const ;

    void findLines       (CLHEP::Hep3Vector* Pos, double* R, double* Slopes) const ;
      
    void findDoublets    (const KalRep* KRep, vector<Doublet>* DCol) const ;
      
    void markDoublet     (KalRep* KRes, Doublet *HitDoublet, int Index0, int Index1) const;
    void markMultiplets  (KalRep* Kres, vector<Doublet>* dcol) const;
    void resolveSingleHit(KalRep* Kres, mu2e::TrkStrawHit* Hit) const ;
    
					// resolve a track.  Depending on the configuration, this might
					// penalty function depends on the drift radius

    double penaltyError(double rdrift) const;
					// update the hit state and the t0 value.
    virtual bool resolveTrk(KalRep* KRes) const;
  };
}
#endif
