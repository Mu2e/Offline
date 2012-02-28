//
// Calorimeter's data cluster container
//
// $Id: CaloCluster.hh,v 1.1 2012/02/28 22:26:01 gianipez Exp $
// $Author: gianipez $
// $Date: 2012/02/28 22:26:01 $
//
// Original author G. Pezzullo & G. Tassielli
//


#ifndef RecoDataProducts_CaloCluster_hh
#define RecoDataProducts_CaloCluster_hh

// Mu2e includes:
#include "art/Persistency/Common/Ptr.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"


// C++ includes
#include <vector>

#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHit.hh"


namespace mu2e {
  
  
  typedef art::Ptr< CaloCrystalHit>                  CaloCrystalHitPtr;
  typedef std::vector<CaloCrystalHitPtr>      CaloCrystalHitPtrVector;

  
  struct CaloCluster{

  private:


  public:
    int                                    _iVane;    //number of vane
    size_t                              _nCrystal;
    float                                   _time;    //(ns)
    //float                                     _dt;    //(ns)
    float                                 _energy;    //(MeV)
    CLHEP::Hep3Vector                _impactPoint;    //center of gravity of the cluster in the mu2e frame
    CLHEP::Hep3Vector           _impactPointError;    //
    CaloCrystalHitPtrVector      _caloCrystalHits;
    
 // public:
    
    CaloCluster():
      _iVane(0),
      _nCrystal(0),
      _time(0.),
      //_dt(0.),
      _energy(0.){
    }

    CaloCluster(int                                    iVane):
      _iVane(iVane),
      _nCrystal(0),
      _time(0.),
      //_dt(0.),
      _energy(0.){
    }

    CaloCluster(int                                    iVane,
		float                                   time,
		//float                                     dt,
		float                                 energy,
		CaloCrystalHitPtrVector              caloCrystalHits):
      _iVane(iVane),
      _time(time),
      //_dt(dt),
      _energy(energy){

            for(CaloCrystalHitPtrVector::iterator it = caloCrystalHits.begin(); it!= caloCrystalHits.end(); ++it){
                    _caloCrystalHits.push_back(*it);
            }
            _nCrystal=caloCrystalHits.size();
    }

    void AddHit (CaloCrystalHitPtr &a){
            _caloCrystalHits.push_back(a);
            _time*=(float)_nCrystal;
            ++_nCrystal;
            _time += a->time();
            _time /= (float)_nCrystal;
            _energy += a->energyDep();
    }

  };
  
} // namespace mu2e

#endif /* RecoDataProducts_CaloCluster_hh */
