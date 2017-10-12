#ifndef RecoDataProducts_CaloRecoDigiFast_hh
#define RecoDataProducts_CaloRecoDigiFast_hh

#include "CLHEP/Vector/TwoVector.h"

namespace mu2e
{

  class CaloRecoDigiFast {

     public:

       CaloRecoDigiFast(): _time(0),_energy(0),_pos(0,0)
       {} 

       CaloRecoDigiFast(double time, double energy, double x, double y) :
	   _time(time),_energy(energy),_pos(x,y)
       {}

       double time()                  const { return _time;} 
       double energy()                const { return _energy;} 
       const CLHEP::Hep2Vector& pos() const { return _pos;} 


     private:

       int                _time;
       double             _energy; 
       CLHEP::Hep2Vector  _pos; 
  };


}
#endif 
