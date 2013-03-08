#ifndef HitMakers_CaloHitSimUtil_hh
#define HitMakers_CaloHitSimUtil_hh
//
// Helper class for CaloReadoutUtilities
//
// Original author B. Echenard
//


// Mu2e includes.
#include "MCDataProducts/inc/StepPointMCCollection.hh"



namespace mu2e {


	 class CaloHitSimUtil {

	       public:

		   CaloHitSimUtil(art::Ptr<StepPointMC> const& step, double edep) : _step(step),_edepTot(edep) {};

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

} // end namespace mu2e

#endif

