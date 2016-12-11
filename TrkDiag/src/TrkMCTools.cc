//
// Namespace for collecting tools used in MC truth evaluation
// Original author: Dave Brown (LBNL) 8/10/2016
//
#include "TrkDiag/inc/TrkMCTools.hh"
#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include <map>

namespace mu2e {
  namespace TrkMCTools {

    int stepPoint(art::Ptr<StepPointMC>& sp, StrawDigiMC const& mcdigi) {
      int retval(-1);
      if(mcdigi.hasTDC(StrawDigi::zero) && mcdigi.hasTDC(StrawDigi::one)) {
	if(mcdigi.wireEndTime(StrawDigi::zero) < mcdigi.wireEndTime(StrawDigi::one)) {
	  sp = mcdigi.stepPointMC(StrawDigi::zero);
	  retval = StrawDigi::zero;
	} else {
	  sp = mcdigi.stepPointMC(StrawDigi::one);
	  retval = StrawDigi::one;
	};
      } else if (mcdigi.hasTDC(StrawDigi::zero)) {
	sp = mcdigi.stepPointMC(StrawDigi::zero);
	retval = StrawDigi::zero;
      } else if (mcdigi.hasTDC(StrawDigi::one)) {
	sp = mcdigi.stepPointMC(StrawDigi::one);
	retval = StrawDigi::one;
      }
      return retval;
    }

    bool CEDigi(StrawDigiMC const& mcdigi) {
      bool conversion(false);
      art::Ptr<StepPointMC> spmcp;
      if(stepPoint(spmcp,mcdigi) >= 0){
	art::Ptr<SimParticle> const& spp = spmcp->simParticle();
	int gid(-1);
	if(spp->genParticle().isNonnull())
	  gid = spp->genParticle()->generatorId().id();
      // a conversion electron is an electron from the CE generator.  The momentum requirement
      // removes cases where the CE loses a catastrophic amount of energy (ie albedo backsplash
      // from the calorimeter).
	conversion = (spp->pdgId() == 11 && gid == 2 && spmcp->momentum().mag()>90.0);
      }
      return conversion;
    }

    unsigned primaryParticle(art::Ptr<SimParticle>& pspp, std::vector<StrawHitIndex> const& hits, const StrawDigiMCCollection* mcdigis) {
      unsigned retval(0);
      std::map<art::Ptr<SimParticle>, unsigned> spmap;

      for( auto hi : hits ) {
	StrawDigiMC const& mcdigi = mcdigis->at(hi);
	art::Ptr<SimParticle> spp;
	if(simParticle(spp,mcdigi) > 0 ){
	  auto ispp = spmap.find(spp);
	  if(ispp != spmap.end())
	    ++(ispp->second);
	  else
	    spmap[spp] = 1;
	}
      }
      for( auto imap : spmap ) {
	if(imap.second > retval){
	  retval = imap.second;
	  pspp = imap.first;
	}
     }
     return retval;
    }

    unsigned countDigis(art::Ptr<SimParticle> const& pspp, const StrawDigiMCCollection* mcdigis) {
      unsigned retval(0);
      for( auto mcdigi : *mcdigis ) {
	art::Ptr<SimParticle> spp;
	if(simParticle(spp,mcdigi) > 0 && spp == pspp) {
	  ++retval;
	}
      }
      return retval;
    }

    unsigned simParticle(art::Ptr<SimParticle>& spp, StrawDigiMC const& mcdigi) {
      unsigned retval(0);
      if(mcdigi.hasTDC(StrawDigi::zero) &&
	  mcdigi.hasTDC(StrawDigi::one) &&
	  mcdigi.stepPointMC(StrawDigi::zero)->simParticle() ==
	  mcdigi.stepPointMC(StrawDigi::one)->simParticle() ) {
	spp = mcdigi.stepPointMC(StrawDigi::zero)->simParticle();
	retval = 2;
      } else if(mcdigi.hasTDC(StrawDigi::zero)) {
	spp = mcdigi.stepPointMC(StrawDigi::zero)->simParticle();
	retval = 1;
      } else if(mcdigi.hasTDC(StrawDigi::one)) {
	spp = mcdigi.stepPointMC(StrawDigi::one)->simParticle();
	retval = 1;
      }
      return retval;
    }
  }
}
