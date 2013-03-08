//
// Original author Bertrand Echenard
//

// C++ includes
#include <ostream>

// Framework includes.
#include "cetlib/exception.h"

// Mu2e includes
#include "MCDataProducts/inc/CaloHitSimPartMC.hh"


namespace mu2e {

   // Print the information found in this hit.
   void CaloHitSimPartMC::print( std::ostream& ost ) const 
   {
      for (unsigned int i=0;i<_simPart.size();++i) 
         ost<<"Hit SimParticle content "<<_simPart[i]->id().asInt()<<" "<<_step[i]->time()<<" "<<_step[i]->momentum().mag()<<"  "<<_edep[i]<<std::endl;
   }

   void CaloHitSimPartMC::add(art::Ptr<SimParticle> const& simPtr, art::Ptr<StepPointMC> const& stepPtr, double edep)
   {
      _simPart.push_back(simPtr);
      _step.push_back(stepPtr);
      _edep.push_back(edep);   
   }

} // namespace mu2e
