#ifndef MCDataProducts_CaloHitSimPartMC_hh
#define MCDataProducts_CaloHitSimPartMC_hh

//
// Record information about SimParticles and StepPonintMC generating hits in the calorimeter. 
// The posiiton refers to the SimParticle's position when it enters the crystal
//
// Original author Bertrand Echenard
//

// Mu2e includes
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/StepPointMC.hh"

//other include
#include "CLHEP/Vector/ThreeVector.h"
#include <iostream>
#include <vector>




namespace mu2e {


   class CaloHitSimPartMC {

       public:

	  CaloHitSimPartMC(): _simPart(),_edep(),_time(),_mom(),_position() {}


	  std::vector<art::Ptr<SimParticle> > const& simParticles() const   {return _simPart;}
	  std::vector<double>                 const& eDep() const           {return _edep;}
	  std::vector<double>                 const& time() const           {return _time;}
	  std::vector<double>                 const& momentum() const       {return _mom;}
	  std::vector<CLHEP::Hep3Vector>      const& position() const       {return _position;}
	  
	  void add(art::Ptr<SimParticle> const& simPtr, double edep, double time, double momentum, CLHEP::Hep3Vector const& position );
	  void print(std::ostream& ost = std::cout) const;


       private:

	    std::vector<art::Ptr<SimParticle> >  _simPart;
	    std::vector<double>                  _edep;
	    std::vector<double>                  _time;
	    std::vector<double>                  _mom;
	    std::vector<CLHEP::Hep3Vector>       _position;
   };

} 

#endif
