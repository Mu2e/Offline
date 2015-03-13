#ifndef RecoDataProducts_CaloProtoCluster_hh
#define RecoDataProducts_CaloProtoCluster_hh
//
// Calorimeter's data proto cluster container
//
// $Id: CaloProtoCluster.hh,v 1.6 2013/03/08 01:22:32 echenard Exp $
// $Author: echenard $
// $Date: 2013/03/08 01:22:32 $
//
// Original author B. Echenard
//

// Mu2e includes:
#include "art/Persistency/Common/Ptr.h"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHit.hh"

// C++ includes
#include <vector>

namespace mu2e {


   class CaloProtoCluster {
       

       public:
    
            typedef art::Ptr< CaloCrystalHit>      CaloCrystalHitPtr;
            typedef std::vector<CaloCrystalHitPtr> CaloCrystalHitPtrVector;

	    CaloProtoCluster() : 
	       _time(-1),_energyDep(-1),_caloCrystalHitsPtrVector(),_isSplit(false)
	    {}

	    CaloProtoCluster(double time, double energy, CaloCrystalHitPtrVector caloCrystalHits, bool isSplit) : 
	       _time(time),_energyDep(energy),_caloCrystalHitsPtrVector(caloCrystalHits),_isSplit(isSplit)
	    {}

	    double                                             time() const{return _time;}            
	    double                                        energyDep() const{return _energyDep;}       
	    CaloCrystalHitPtrVector const& caloCrystalHitsPtrVector() const{return _caloCrystalHitsPtrVector;}
	    bool                                            isSplit() const{return _isSplit;} 


	 private:
	 
	    double                   _time;       
	    double                   _energyDep;  
	    CaloCrystalHitPtrVector  _caloCrystalHitsPtrVector;
	    bool                     _isSplit;    

   };

} 

#endif /* RecoDataProducts_CaloProtoCluster_hh */
