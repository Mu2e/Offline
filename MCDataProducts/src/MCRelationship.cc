//
// MC relationships of some objects.  Extracted from KalDiag
// Dave Brown, LBNL, 20 Jun 2016
//
#include "MCDataProducts/inc/MCRelationship.hh"
#include <vector>
namespace mu2e
{
  using std::vector;

  MCRelationship::MCRelationship(StrawDigiMC const& mcd1, StrawDigiMC const& mcd2) :
    MCRelationship(mcd1.earlyStepPointMC()->simParticle(),mcd2.earlyStepPointMC()->simParticle())
  {}

  MCRelationship::MCRelationship(StrawDigiMC const& mcd, SPPtr const& spp) :
    MCRelationship(mcd.earlyStepPointMC()->simParticle(),spp)
  {}

  MCRelationship::MCRelationship(SPPtr const& sppi,SPPtr const& sppj) : _rel(none) {
    if(sppi.isNonnull() && sppj.isNonnull()){
      if(sppi == sppj){
	_rel= same;
      } else {
	SPPtr pi = sppi->originParticle().parent();
	SPPtr pj = sppj->originParticle().parent();
	if(pi.isNonnull() && pi == sppj){
	  _rel = daughter;
	} else if(pj.isNonnull() && pj == sppi) {
	  _rel = mother;
	} else if(pi.isNonnull() && pj.isNonnull()){
	  if( pi == pj){
	    _rel = sibling;
	  } else {
	    vector<SPPtr > pvi, pvj;
	    pvi.push_back(sppi);
	    pvj.push_back(sppj);
	    while(pi.isNonnull()){
	      pvi.push_back(pi);
	      pi = pi->originParticle().parent();
	    }
	    while(pj.isNonnull()){
	      pvj.push_back(pj);
	      pj = pj->originParticle().parent();
	    }
	    if(find(pvi.begin(),pvi.end(),sppj) != pvi.end()){
	      _rel = udaughter;
	    } else if(find(pvj.begin(),pvj.end(),sppi) != pvj.end()){
	      _rel = umother;
	    } else {
	      for(size_t ii=0;ii<pvj.size();++ii){
		if(find(pvi.begin(),pvi.end(),pvj[ii]) != pvi.end()){
		  _rel = usibling;
		  break;
		} 
	      }
	      for(size_t ii=0;ii<pvi.size();++ii){
		if( find(pvj.begin(),pvj.end(),pvi[ii]) != pvj.end()){
		  _rel = usibling;
		  break;
		}
	      }
	    }
	  }
	}
      }
    }
  }
}
 
