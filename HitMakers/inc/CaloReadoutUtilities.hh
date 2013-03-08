#ifndef HitMakers_CaloReadoutUtilities_hh
#define HitMakers_CaloReadoutUtilities_hh
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

// Mu2e includes.
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "HitMakers/inc/CaloHitSimUtil.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/CaloHitSimPartMCCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"





namespace mu2e {


     class CaloReadoutUtilities{

	private:

	  typedef art::Ptr<StepPointMC>   StepPtr;
	  typedef std::vector<StepPtr >   StepPtrs;
	  typedef art::Ptr<SimParticle>   SimPtr;
	  typedef std::vector<SimPtr >    SimPtrs;
          typedef std::map<SimPtr, CaloHitSimUtil> SimStepMap;


	public:

	  CaloReadoutUtilities() {}
	  void fillSimMother(Calorimeter const& cal, 
                             PtrStepPointMCVector const& mcptr,
			     CaloHitSimPartMC& caloHitSimPartMC);
			     

     };

} 

#endif /* HitMakers_CaloReadoutUtilities_hh */
