#ifndef CaloMC_CrystalMCUtil_hh
#define CaloMC_CrystalMCUtil_hh
//
// CaloReadout Utility to study the MC content of a calo hit:
//   extract the SimParticle(s) originating from outside the Calo producing the hit
//   extract the first StepPointMC in the hit from this SimParticle (to estimate the energy/time)
//
// Original author B. Echenard
//

// C++ includes
#include <vector>
#include <map>
#include <unordered_map>

// Mu2e includes.
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/CaloHitSimPartMCCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"



namespace mu2e {


     class CaloHitMCUtil {

	   public:

	       CaloHitMCUtil(art::Ptr<StepPointMC> const& step, double edep) : _step(step),_edepTot(edep) {};

	       art::Ptr<StepPointMC>  const& step()    const {return _step;}
	       double                        edepTot() const {return _edepTot;}


	       void update(art::Ptr<StepPointMC> const& step, double edep)
	       {
        	  _edepTot += edep;	
		  if ( step->momentum().mag() > _step->momentum().mag() ) _step = step;
	       }


	   private:

	       art::Ptr<StepPointMC> _step;
	       double                _edepTot;
     };




     class CrystalMCUtil {

	 public:

	   CrystalMCUtil() {}
	   void fillSimMother(Calorimeter const& cal, PtrStepPointMCVector const& mcptr,CaloHitSimPartMC& caloHitSimPartMC, 
	                      const std::unordered_map<const StepPointMC*, double> &timeMap);


	 private:

	   typedef art::Ptr<StepPointMC>   StepPtr;
	   typedef art::Ptr<SimParticle>   SimPtr;
	   typedef std::vector<StepPtr >   StepPtrs;
	   typedef std::vector<SimPtr >    SimPtrs;
           typedef std::map<SimPtr, CaloHitMCUtil> SimStepMap;

     };

} 

#endif
