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
      if(mcdigi.wireEndTime(StrawEnd::cal) < mcdigi.wireEndTime(StrawEnd::hv)) {
	sp = mcdigi.stepPointMC(StrawEnd::cal);
	retval = StrawEnd::cal;
      } else {
	sp = mcdigi.stepPointMC(StrawEnd::hv);
	retval = StrawEnd::hv;
      };
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
      if( mcdigi.stepPointMC(StrawEnd::cal)->simParticle() ==
	  mcdigi.stepPointMC(StrawEnd::hv)->simParticle() ) {
	spp = mcdigi.stepPointMC(StrawEnd::cal)->simParticle();
	retval = 2;
      }
      return retval;
    }


    void findMCTrk(const KalSeed& kseed,art::Ptr<SimParticle>& spp, const StrawDigiMCCollection& mcdigis) {
      static art::Ptr<SimParticle> nullvec;
      spp = nullvec;
      std::vector<spcount> sct;
      findMCTrk(kseed,sct, mcdigis);
      if(sct.size()>0)
	spp = sct[0]._spp;
    }

    void findMCTrk(const KalSeed& kseed,std::vector<spcount>& sct, const StrawDigiMCCollection& mcdigis) {
      sct.clear();
      // find the SimParticles which contributed hits.
      // loop through the straw hits from the track
      static StrawHitFlag active(StrawHitFlag::active);
      for(const auto& i_hit : kseed.hits()) {
	// loop over the hits and find the associated steppoints
	if(i_hit.flag().hasAllProperties(active)){
	  StrawDigiMC const& mcdigi = mcdigis.at(i_hit.index());
	  StrawEnd itdc;
	  art::Ptr<SimParticle> spp = mcdigi.stepPointMC(itdc)->simParticle();
	  // see if this particle has already been found; if so, increment, if not, add it
	  bool found(false);
	  for(size_t isp=0;isp<sct.size();++isp){
	    // count direct daughter/parent as part the same particle
	    if(sct[isp]._spp == spp ){
	      found = true;
	      sct[isp].append(spp);
	      break;
	    }
	  }
	  if(!found)sct.push_back(spp);
	}
      }
      // sort by # of contributions
      sort(sct.begin(),sct.end(),spcountcomp());
    }

    void findMCSteps(const StepPointMCCollection& mcsteps, cet::map_vector_key const& trkid, std::vector<int> const& vids, std::vector<MCStepItr>& steps) {
      steps.clear();
      // Loop over the step points, and find the one corresponding to the given detector
      for( MCStepItr imcs =mcsteps.begin();imcs!= mcsteps.end();imcs++){
	if(vids.size() == 0 ||  (imcs->trackId() == trkid && find(vids.begin(),vids.end(),imcs->volumeId()) != vids.end())){
	  steps.push_back(imcs);
	}
      }
      // sort these in time
      sort(steps.begin(),steps.end(),timecomp());
    }

    void countHits(const KalSeed& kseed, art::Ptr<SimParticle>& spp, const StrawDigiMCCollection& mcdigis, int& nactive, int& nhits, int& ngood, int& nambig) {
      nactive = 0; nhits = 0; ngood = 0; nambig = 0;
      static StrawHitFlag active(StrawHitFlag::active);
      //      Hep3Vector mcmomvec = spp->startMomentum();
      //      double mcmom = mcmomvec.mag();
      for(const auto& ihit : kseed.hits()) {
	StrawDigiMC const& mcdigi = mcdigis.at(ihit.index());
	StrawEnd itdc;
	art::Ptr<StepPointMC> const& spmcp = mcdigi.stepPointMC(itdc);
	if(spp == spmcp->simParticle()){
	  ++nhits;
	  // easiest way to get MC ambiguity is through info object
	  /*	  TrkStrawHitInfoMC tshinfomc;
	  fillHitInfoMC(spp,mcdigi,tsh->straw(),tshinfomc);
	  // count hits with at least givien fraction of the original momentum as 'good'
	  if(tshinfomc._mom/mcmom > _mingood )++trkinfomc._ngood;
	  */
	  if(ihit.flag().hasAllProperties(active)){
	    ++nactive;
	    // count hits with correct left-right iguity
	    //	    if(tsh->ambig()*tshinfomc._ambig > 0)++trkinfomc._nambig;
	  }
	}
      }      
    }

    void countDigis(art::Ptr<SimParticle>& spp, const StrawDigiMCCollection& mcdigis, int& ndigi, int& ndigigood) {
      ndigi = 0; ndigigood = 0;
      for(auto imcd = mcdigis.begin(); imcd !=mcdigis.end();++imcd){
	if( imcd->stepPointMC(StrawEnd::cal)->simParticle() == spp){
	  ndigi++;
	  /*	  if(imcd->stepPointMC(StrawEnd::cal)->momentum().mag()/spp->startMomentum().mag() > _mingood) {
		  ndigigood++;
	  }
	  */
	}
      }
    }
  }
}
