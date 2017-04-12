//
// MC relationships of some objects.  Extracted from KalDiag
// Dave Brown, LBNL, 20 Jun 2016
//
#include "MCDataProducts/inc/MCRelationship.hh"
#include <vector>
namespace mu2e
{
  using std::vector;

  MCRelationship::relation MCRelationship::relationship(StrawDigiMC const& mcd1, StrawDigiMC const& mcd2) {
    SPPtr ptr1, ptr2;
    if(mcd1.stepPointMC(StrawDigi::zero).isNonnull())
      ptr1 = mcd1.stepPointMC(StrawDigi::zero)->simParticle();
    else if(mcd1.stepPointMC(StrawDigi::one).isNonnull())
      ptr1 = mcd1.stepPointMC(StrawDigi::one)->simParticle();
    if(mcd2.stepPointMC(StrawDigi::zero).isNonnull())
      ptr2 = mcd2.stepPointMC(StrawDigi::zero)->simParticle();
    else if(mcd2.stepPointMC(StrawDigi::one).isNonnull())
      ptr2 = mcd2.stepPointMC(StrawDigi::one)->simParticle();
    return relationship(ptr1,ptr2);
  }

  MCRelationship::relation MCRelationship::relationship(SPPtr const& sppi,SPPtr const& sppj) {
    if(sppi.isNull() || sppj.isNull()) return none;
    if(sppi == sppj)return same;
    SPPtr pi = sppi->originParticle().parent();
    SPPtr pj = sppj->originParticle().parent();

    if(pi.isNonnull() && pi == sppj)return daughter;
    if(pj.isNonnull() && pj == sppi)return mother;
    if(pi.isNonnull() && pj.isNonnull()){
      if( pi == pj)return sibling;
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
      vector<art::Ptr<SimParticle> >::iterator ifnd;
      ifnd = find(pvi.begin(),pvi.end(),sppj);
      if(ifnd != pvi.end())return udaughter;
      ifnd = find(pvj.begin(),pvj.end(),sppi);
      if(ifnd != pvj.end())return umother;
      for(size_t ii=0;ii<pvj.size();++ii){
	ifnd = find(pvi.begin(),pvi.end(),pvj[ii]);
	if(ifnd != pvi.end())return usibling;
      }
      for(size_t ii=0;ii<pvi.size();++ii){
	ifnd = find(pvj.begin(),pvj.end(),pvi[ii]);
	if(ifnd != pvj.end())return usibling;
      }
    }
    return none;
  }
}
 
