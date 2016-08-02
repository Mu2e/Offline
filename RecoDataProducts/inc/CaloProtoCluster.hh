#ifndef RecoDataProducts_CaloProtoCluster_hh
#define RecoDataProducts_CaloProtoCluster_hh

#include "art/Persistency/Common/Ptr.h"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include <vector>

namespace mu2e {


   class CaloProtoCluster {
       
       public:
    
            typedef art::Ptr< CaloCrystalHit>      CaloCrystalHitPtr;
            typedef std::vector<CaloCrystalHitPtr> CaloCrystalHitPtrVector;

	    CaloProtoCluster() : 
	       _time(-1),_energyDep(-1),_caloCrystalHitPtrVector(),_isSplit(false)
	    {}

	    CaloProtoCluster(double time, double energy, CaloCrystalHitPtrVector CaloCrystalHit, bool isSplit) : 
	       _time(time),_energyDep(energy),_caloCrystalHitPtrVector(CaloCrystalHit),_isSplit(isSplit)
	    {}


	    double                         time()                    const {return _time;}            
	    double                         energyDep()               const {return _energyDep;}       
	    const CaloCrystalHitPtrVector& caloCrystalHitPtrVector() const {return _caloCrystalHitPtrVector;}
	    bool                           isSplit()                 const {return _isSplit;} 


        private:
	 
	    double                   _time;       
	    double                   _energyDep;  
	    CaloCrystalHitPtrVector  _caloCrystalHitPtrVector;
	    bool                     _isSplit;    

   };

} 

#endif 
