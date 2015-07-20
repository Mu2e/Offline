#ifndef Mu2eUtilities_CaloHitMCNavigator_hh
#define Mu2eUtilities_CaloHitMCNavigator_hh
//
// Helper class to navigate the MC Truth information associated with a CaloHit.
//
// $Id: CaloHitMCNavigator.hh,v 1.1 2013/03/08 01:22:32 echenard Exp $
// $Author: echenard $
// $Date: 2013/03/08 01:22:32 $
//
// Original author Rob Kutschke
//

#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "MCDataProducts/inc/CaloHitMCTruthCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/CaloHitSimPartMCCollection.hh"

namespace mu2e {

  class CaloHitMCNavigator{

      public:

	   CaloHitMCNavigator( CaloHitCollection const&               hits,
                               CaloHitMCTruthCollection const&        truth,
                               CaloHitSimPartMCCollection const&      sims);

	   // Accept compiler supplied d'tor, copy c'tor and assignment operator.


	   // Access the underlying collections.
	   CaloHitCollection              const& hits()   const {return *_hits;}
	   CaloHitMCTruthCollection       const& truths() const {return *_truth;}
	   CaloHitSimPartMCCollection     const& sims()   const {return *_sims;}


	   // Find index in the correpsonding vector
	   size_t index(CaloHit const& hit) const                {return &hit-&_hits->at(0);}
	   size_t index(CaloHitMCTruth const& truth) const       {return &truth-&_truth->at(0);}
	   size_t index(CaloHitSimPartMC const& sim) const       {return &sim-&_sims->at(0);}

	   template <class T>  CaloHit              const& hit(T const& element) const   {return _hits->at(index(element));}
	   template <class T>  CaloHitMCTruth       const& truth(T const& element) const {return _truth->at(index(element));}
	   template <class T>  CaloHitSimPartMC     const& sim(T const& element)  const  {return _sims->at(index(element));}



	 private:

	   CaloHitCollection const*              _hits;
	   CaloHitMCTruthCollection const*       _truth;
	   CaloHitSimPartMCCollection const*     _sims;

       };

}

#endif

