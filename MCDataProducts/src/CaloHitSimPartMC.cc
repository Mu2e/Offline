//
// Original author Bertrand Echenard
//

// C++ includes
#include <ostream>
#include "cetlib_except/exception.h"

//Mu2e and CLHEP includes
#include "MCDataProducts/inc/CaloHitSimPartMC.hh"
#include "CLHEP/Vector/ThreeVector.h"



namespace mu2e {

   
     void CaloHitSimPartMC::print( std::ostream& ost ) const 
     {
	for (unsigned int i=0;i<_simPart.size();++i) 
           ost<<"Hit SimParticle content "<<_simPart[i]->id().asInt()<<" edep="<<_edep[i]<<" time="<<_time[i]<<" mom="<<_mom[i]<<" x,y,z="<<_position[i]<<std::endl;
     }


     void CaloHitSimPartMC::add(art::Ptr<SimParticle> const& simPtr, double edep, double time, double momentum, CLHEP::Hep3Vector const& position)
     {
	_simPart.push_back(simPtr);
	_edep.push_back(edep);   
	_time.push_back(time);   
	_mom.push_back(momentum);   
	_position.push_back(position);   
     }

} 
