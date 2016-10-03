//
// Container for the info of the extrapolated trajectory on the calorimeter
//
// $Id: TrkCaloIntersect.hh,v 1.9 2014/08/20 14:23:09 murat Exp $
// $Author: murat $
// $Date: 2014/08/20 14:23:09 $
//
// Original author G. Pezzullo / amended by B. Echenard
//


#ifndef RecoDataProducts_TrkCaloMatch_hh
#define RecoDataProducts_TrkCaloMatch_hh


// Mu2e includes:
#include "canvas/Persistency/Common/Ptr.h"


//tracker includes
#include "BTrk/KalmanTrack/KalRep.hh"
#include "RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/TrkCaloIntersect.hh"




namespace mu2e {


  class TrkCaloMatch {


       public:

	   TrkCaloMatch() {}

	   TrkCaloMatch(art::Ptr<TrkCaloIntersect> const& intersect, art::Ptr<CaloCluster> const& cluster, int cluId, double chi2, double chi2Pos, double chi2Time) :
	      _intersect(intersect), _cluster(cluster), _cluId(cluId), _chi2(chi2), _chi2Pos(chi2Pos), _chi2Time(chi2Time)
	   {}

	   ~TrkCaloMatch(){}

	   art::Ptr<TrkCaloIntersect>  const&  intersect()  const {return _intersect;}
	   KalRepPtr                   const&  trk()        const {return _intersect->trk();}
	   int                                 trkId()      const {return _intersect->trkId();}
	   
	   art::Ptr<CaloCluster>       const&  cluster()    const {return _cluster;}
	   int                                 cluId()      const {return _cluId;}
           
	   double                              chi2()       const {return _chi2;} 
           double                              chi2Pos()    const {return _chi2Pos;} 
           double                              chi2Time()   const {return _chi2Time;} 



       private:

	   art::Ptr<TrkCaloIntersect>   _intersect;
	   art::Ptr<CaloCluster>        _cluster;
	   int                          _cluId;	   	   
	   double                       _chi2;	   
	   double                       _chi2Pos;	   
	   double                       _chi2Time;	   

   };


} 

#endif 

