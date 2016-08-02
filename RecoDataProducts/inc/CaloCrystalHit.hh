#ifndef RecoDataProducts_CaloCrystalHit_hh
#define RecoDataProducts_CaloCrystalHit_hh

#include "RecoDataProducts/inc/CaloRecoDigi.hh"
#include "RecoDataProducts/inc/CaloRecoDigiCollection.hh"
#include "art/Persistency/Common/Ptr.h"

#include <vector>


namespace mu2e {

  class CaloCrystalHit
  {

      public:

	 CaloCrystalHit(): _crystalId(-1),_time(0.),_eDep(0.),_eDepTot(0.),_nROId(0),_recoCaloDigis()
	 {}

	 CaloCrystalHit(int crystalId, double time, double eDep, double eDepTot, int nRoid, 
	                std::vector<art::Ptr<CaloRecoDigi> > &CaloRecoDigi)  :
              _crystalId(crystalId),_time(time),_eDep(eDep),_eDepTot(eDepTot),_nROId(nRoid),
	     _recoCaloDigis(CaloRecoDigi)
	 {}


	 int    id()                                                 const { return _crystalId; }
	 double time()                                               const { return _time;}
         double energyDep()                                          const { return _eDep;} 
	 double energyDepTot()                                       const { return _eDepTot;}
	 int    nROId()                                              const { return _nROId;}
	 const std::vector<art::Ptr<CaloRecoDigi> >& recoCaloDigis() const {return _recoCaloDigis;}


       private:

	 int        _crystalId;
	 double     _time;             
	 double     _eDep;        
	 double     _eDepTot;   
	 int        _nROId;
	 std::vector<art::Ptr<CaloRecoDigi> > _recoCaloDigis;

  };


}

#endif
