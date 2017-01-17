#ifndef RecoDataProducts_CaloProtoCluster_hh
#define RecoDataProducts_CaloProtoCluster_hh

#include "canvas/Persistency/Common/Ptr.h"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include <vector>

namespace mu2e {


   class CaloProtoCluster {
       
       public:
    
            typedef art::Ptr< CaloCrystalHit>      CaloCrystalHitPtr;
            typedef std::vector<CaloCrystalHitPtr> CaloCrystalHitPtrVector;

	    CaloProtoCluster() : 
	       _time(0.),_timeErr(0.0),_energyDep(0.),_energyDepErr(0.0),_caloCrystalHitsPtrVector(),_isSplit(false)
	    {}

	    CaloProtoCluster(double time, double timeErr, double energy, double energyErr, 
	                     CaloCrystalHitPtrVector CaloCrystalHit, bool isSplit) : 
	       _time(time),_timeErr(timeErr),_energyDep(energy),_energyDepErr(energyErr),
	       _caloCrystalHitsPtrVector(CaloCrystalHit),_isSplit(isSplit)
	    {}


	    double                         time()                     const {return _time;}            
	    double                         timeErr()                  const {return _timeErr;}            
	    double                         energyDep()                const {return _energyDep;}       
	    double                         energyDepErr()             const {return _energyDepErr;}       
	    const CaloCrystalHitPtrVector& caloCrystalHitsPtrVector() const {return _caloCrystalHitsPtrVector;}
	    bool                           isSplit()                  const {return _isSplit;} 



        private:
	 
	    double                   _time;       
	    double                   _timeErr;       
	    double                   _energyDep;  
	    double                   _energyDepErr;  
	    CaloCrystalHitPtrVector  _caloCrystalHitsPtrVector;
	    bool                     _isSplit;    

   };

} 

#endif 
