///////////////////////////////////////////////////////////////////////////////
// 2015-04-08 P.Murat:
// -------------------
// resolve hit ambiguities by panel, assuming a reasonable track fit as input
///////////////////////////////////////////////////////////////////////////////

#ifndef DoubletAmbigResolver_HH
#define DoubletAmbigResolver_HH

// #include "BTrk/BaBar/BaBar.hh"
   #include "Offline/TrkReco/inc/AmbigResolver.hh"

#ifndef __GCCXML__
#include "fhiclcpp/ParameterSet.h"
#endif/*__GCCXML__*/
#include <cstddef>
#include <vector>
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

class TrkDifTraj;
class KalRep;

#include "Offline/BTrkData/inc/Doublet.hh"

namespace mu2e {

  class Straw;

  class DoubletAmbigResolver : public AmbigResolver {
    public:
      enum trajtype {reftraj=0};

      struct Data_t {
        int                index[2];
        int                ibest;
        int                inext;

        double             chi2min;
        double             chi2next;

        const mu2e::Straw* straw     [2];
        double             rdrift    [2];
        double             doca      [4][2];
        double             chi2Slope [4]; // slope contribution to the total chi2
        double             chi2Coord [4]; // coordinate contribution to the total chi2
        double             chi2      [4]; // total chi2
        double             trkslope;
        double             lineSlopes[4];
        CLHEP::Hep3Vector  spos      [2];
        CLHEP::Hep3Vector  sposr     [2];
        CLHEP::Hep3Vector  tpos      [2];
        CLHEP::Hep3Vector  tposr     [2];
      };

    protected:
      double _extErr;
      int    _debugLevel;

      double _mindrift;                   // minimum drift to assign an ambiguity.  Below this, an ambiguity of '0' is defined
      double _zeropenalty;                // special penalty for drifts below the minimum
      bool   _penalty;      // apply penalty or notA
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
      double _minChi2Ratio;               // if chi2(best)/chi2(next) < _minChi2Ratio, doublet is "well measured"
      double _tempScale;                  //
      double _penaltyScale;               //
      int    _useMeanResidual;            // as a penalty, if _meanResidual is defined
      double _maxMeanResidual;            // max mean residual - can't make it very large
      double _maxHitChi;

      int    _sign[4][2];
      int    _iter;                       // iteration
      int    _Final;                      // final iteration
      double _meanResidual;               // mean residual (est) for active hits on a track
      //-----------------------------------------------------------------------------
      // constructors and destructor
      //-----------------------------------------------------------------------------
    public:
#ifndef __GCCXML__
      explicit DoubletAmbigResolver(fhicl::ParameterSet const&, double ExtErr, int Iter, int Final);
#endif/*__GCCXML__*/
      virtual ~DoubletAmbigResolver();

      int  calculateDoubletParameters(const KalRep* KRep, Doublet* HitDoublet, Data_t* R) const ;
      void defineHitDriftSign        (mu2e::TrkStrawHit* Hit, int I, Data_t* R) const ;

      void findLines                 (CLHEP::Hep3Vector* Pos, double* R, double* Slopes) const ;
      void findDoublets              (const KalRep* KRep, vector<Doublet>* ListOfDoublets) const ;
      //-----------------------------------------------------------------------------
      // three functions below modify KRes - update hit drift directions and assign penalty errors
      //-----------------------------------------------------------------------------
      void markDoublet               (KalRep* KRes, Doublet *HitDoublet, int Index0, int Index1) const;
      void markMultiplets            (KalRep* Kres, vector<Doublet>* ListOfDoublets) const ;
      void resolveSingleHit          (KalRep* Kres, mu2e::TrkStrawHit* Hit) const ;
      //-----------------------------------------------------------------------------
      // overloaded functions of the base class
      // resolve a track.  Depending on the configuration, this might
      // make penalty error virtual ? penalty function depends on the drift radius
      //-----------------------------------------------------------------------------
      double penaltyError(double rdrift) const;

      // update the hit state and the t0 value.
      virtual bool resolveTrk   (KalRep* KRes) const;

      virtual void initHitErrors(KalRep* krep) const ;

  };
}
#endif
