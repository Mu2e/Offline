#ifndef MCDataProducts_CaloHitSimPartMC_hh
#define MCDataProducts_CaloHitSimPartMC_hh

//
// Record information about SimParticles and StepPointMC generating hits in the calorimeter
//
// Original author Bertrand Echenard
//

// Mu2e includes
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/StepPointMC.hh"


// c++ includes
#include <iostream>
#include <vector>




namespace mu2e {

   class CaloHitSimPartMC {

       public:

	  CaloHitSimPartMC(): _simPart(),_step(),_edep() {}

	  CaloHitSimPartMC(std::vector<art::Ptr<SimParticle> > const& simPart, 
                	   std::vector<art::Ptr<StepPointMC> > const& step, 
			   std::vector<double> edep) : 
			   _simPart(simPart),_step(step),_edep(edep) {}

	  // Accept compiler generated versions of d'tor, copy c'tor, assignment operator.

	  std::vector<art::Ptr<SimParticle> > const& simParticles() const   {return _simPart;}
	  std::vector<art::Ptr<StepPointMC> > const& stepPoints() const     {return _step;}
	  std::vector<double>                 const& eDep() const           {return _edep;}

	  
	  void add(art::Ptr<SimParticle> const& simPtr, art::Ptr<StepPointMC> const& stepPtr, double edep);
	  void print(std::ostream& ost = std::cout) const;


       private:

	    std::vector<art::Ptr<SimParticle> >  _simPart;
	    std::vector<art::Ptr<StepPointMC> >  _step;
	    std::vector<double>                  _edep;
   };


} 

#endif /* MCDataProducts_CaloHitSimPartMC_hh */
