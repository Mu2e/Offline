#ifndef RecoDataProducts_CaloHit_hh
#define RecoDataProducts_CaloHit_hh

#include "RecoDataProducts/inc/CaloRecoDigi.hh"
#include "RecoDataProducts/inc/CaloRecoDigi.hh"
#include "canvas/Persistency/Common/Ptr.h"

#include <vector>
#include <map>

namespace mu2e {

  class CaloHit
  {
      public:
	 CaloHit(): _crystalId(-1),_nROId(0),_time(0.),_timeErr(0.),_eDep(0.),_eDepErr(0.),_recoCaloDigis()
	 {}

	 CaloHit(int crystalId, int nRoid, double time, double timeErr, double eDep, double eDepErr,  
	                std::vector<art::Ptr<CaloRecoDigi> > &CaloRecoDigi)  :
                        _crystalId(crystalId),_nROId(nRoid),_time(time),_timeErr(timeErr),_eDep(eDep),_eDepErr(eDepErr),
			_recoCaloDigis(CaloRecoDigi)
         {}

	       int                                    id             () const { return _crystalId; }
	       int                                    nROId          () const { return _nROId;}
	       double                                 time           () const { return _time;}
	       double                                 timeErr        () const { return _timeErr;}
               double                                 energyDep      () const { return _eDep;} 
               double                                 energyDepErr   () const { return _eDepErr;} 
	       double                                 energyDepTot   () const { return _eDep*_nROId;}
	       double                                 energyDepTotErr() const { return _eDepErr*_nROId;}
	 const std::vector<art::Ptr<CaloRecoDigi>>&   recoCaloDigis  () const { return _recoCaloDigis;}


       private:
	 int        _crystalId;
	 int        _nROId;
	 double     _time;             
	 double     _timeErr;             
	 double     _eDep;        
	 double     _eDepErr;        
	 std::vector<art::Ptr<CaloRecoDigi> > _recoCaloDigis;

  };


  typedef std::vector<mu2e::CaloHit> CaloHitCollection;
  typedef std::map<art::Ptr<CaloHit>,art::Ptr<CaloHit>> CaloHitRemapping;
}

#endif
