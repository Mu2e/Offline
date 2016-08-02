#ifndef MCDataProducts_CaloClusterContentMC_hh
#define MCDataProducts_CaloClusterContentMC_hh

//
// Record information about cluster content for each simParticle, patlod for the art::Assn<CaloCluster, SimParticle>
//
// Original author Bertrand Echenard
//


#include "MCDataProducts/inc/SimParticle.hh"
#include "CLHEP/Vector/ThreeVector.h"




namespace mu2e {


   class CaloClusterContentMC {


       public:

	  CaloClusterContentMC(): _simPart(),_eDep(0),_time(0),_mom(),_position() {}
	  
	  CaloClusterContentMC(art::Ptr<SimParticle> sim, double eDep, double time, double mom, const CLHEP::Hep3Vector& position) : 
	    _simPart(sim),_eDep(eDep),_time(time),_mom(mom),_position(position) 
	  {}


	  const art::Ptr<SimParticle>& simParticles() const {return _simPart;}
	  double                       eDep()         const {return _eDep;}
	  double                       time()         const {return _time;}
	  double                       momentum()     const {return _mom;}
	  const CLHEP::Hep3Vector&     position()     const {return _position;}
	  



       private:

	    art::Ptr<SimParticle>   _simPart;
	    double                  _eDep;
	    double                  _time;
	    double                  _mom;
	    CLHEP::Hep3Vector       _position;
   };

} 

#endif
