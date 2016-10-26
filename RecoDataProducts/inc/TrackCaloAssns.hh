#ifndef RecoDataProducts_inc_TrackCaloAssns_hh
#define RecoDataProducts_inc_TrackCaloAssns_hh

#include "RecoDataProducts/inc/KalRepCollection.hh"
#include "RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/Assns.h"

namespace mu2e {
  class TrackCaloMatchInfo {
  public:
    TrackCaloMatchInfo() {};
    TrackCaloMatchInfo(
		         double chi2
		       , double chi2Pos
		       , double chi2Time
		       , double kFrac
		       , double kinetic
		       , double eOverP)
    { 
      _chi2 = chi2;
      _chi2Pos = chi2Pos;
      _chi2Time = chi2Time;
      _kFrac = kFrac;
      _kinetic = kinetic;
      _eOverP = eOverP;
    }


    double                         getChi2()         const {return _chi2;}
    double                         getChi2Pos()      const {return _chi2Pos;}
    double                         getChi2Time()     const {return _chi2Time;}
    double                         getKineticFrac()  const {return _kFrac;}
    double                         getKinetic()      const {return _kinetic;}
    double                         getEOverP()       const {return _eOverP;}


  private:
    double _chi2;
    double _chi2Pos;
    double _chi2Time;
    double _kFrac;
    double _kinetic;
    double _eOverP;

  };

 
  // use as many-to-many Assns (want track, cluster, payload not Ptr's)
     typedef art::Assns<KalRepPtr,CaloCluster,TrackCaloMatchInfo> TrackCaloMatchAssns;
  //   typedef art::Assns<KalRep,CaloCluster,TrackCaloMatchInfo> TrackCaloMatchAssns;
}
#endif /*RecoDataProducts_inc_TrackCaloAssns_hh*/

