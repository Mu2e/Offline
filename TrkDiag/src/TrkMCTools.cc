//
// Namespace for collecting tools used in MC truth evaluation
// Original author: Dave Brown (LBNL) 8/10/2016
//
#include "TrkDiag/inc/TrkMCTools.hh"
#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/MCRelationship.hh"

#include "TrackerGeom/inc/Tracker.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"
#include "BTrk/TrkBase/TrkHelixUtils.hh"

#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "BFieldGeom/inc/BFieldManager.hh"

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


    void findMCTrk(const KalSeed& kseed,art::Ptr<SimParticle>& spp, StrawDigiMCCollection const& mcdigis) {
      static art::Ptr<SimParticle> nullvec;
      spp = nullvec;
      std::vector<spcount> sct;
      findMCTrk(kseed,sct, mcdigis);
      if(sct.size()>0)
	spp = sct[0]._spp;
    }

    void findMCTrk(const KalSeed& kseed,std::vector<spcount>& sct, StrawDigiMCCollection const& mcdigis) {
      sct.clear();
      // find the SimParticles which contributed hits.
      // loop through the straw hits from the track
      static StrawHitFlag active(StrawHitFlag::active);
      for(const auto& i_hit : kseed.hits()) {
	// loop over the hits and find the associated steppoints
	bool isactive = i_hit.flag().hasAllProperties(active);
	StrawDigiMC const& mcdigi = mcdigis.at(i_hit.index());
	art::Ptr<SimParticle> spp = mcdigi.earlyStepPointMC()->simParticle();
	// see if this particle has already been found; if so, increment, if not, add it
	bool found(false);
	for(size_t isp=0;isp<sct.size();++isp){
	  // count direct daughter/parent as part the same particle
	  if(sct[isp]._spp == spp ){
	    found = true;
	    sct[isp].append(spp,isactive);
	    break;
	  }
	}
	if(!found)sct.push_back(spcount(spp,isactive));
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

    void countDigis(const KalSeedMC& kseedmc, const KalSeed& kseed, int& ndigi, int& ndigigood, int& nambig) {
      ndigi = 0; ndigigood = 0, nambig = 0;
      
      for(size_t i_digi = 0; i_digi < kseedmc._tshmcs.size(); ++i_digi) {
	const auto& i_tshmc = kseedmc._tshmcs.at(i_digi);

	if (kseedmc.simParticle(i_tshmc._spindex)._rel == MCRelationship::same) {
	  ++ndigi;

	  ++ndigigood;

	  //	  if(spmcp->momentum().mag()/spp->startMomentum().mag() > mingood) {
	  //	  ++ngood;
	  // }

	  // easiest way to get MC ambiguity is through info object
	  TrkStrawHitInfoMC tshinfomc;
	  fillHitInfoMC(kseedmc,tshinfomc,i_tshmc);  

	  const auto& ihit = kseed.hits().at(i_digi);
	  if(ihit.ambig()*tshinfomc._ambig > 0) {
	    ++nambig;
	  }
	}
      }
    }

    void fillTrkInfoMC(const KalSeedMC& kseedmc, const KalSeed& kseed, 
    TrkInfoMC& trkinfomc) {
      // use the primary match of the track
      if(kseedmc.simParticles().size() > 0){
	auto const& simp = kseedmc.simParticles().front();
	trkinfomc._gen = simp._gid.id();
	trkinfomc._pdg = simp._pdg;
	trkinfomc._proc = simp._proc;
	trkinfomc._nhits = simp._nhits; // number of hits from the primary particle
	trkinfomc._nactive = simp._nactive; // number of active hits from the primary particle
	trkinfomc._prel = simp._rel.relationship(); // relationship of the track primary to the event primary
      }

      int ndigi = -1, ndigigood = -1, nambig = -1;
      TrkMCTools::countDigis(kseedmc, kseed, ndigi, ndigigood, nambig);
      trkinfomc._ndigi = ndigi;
      trkinfomc._ndigigood = ndigigood;
      trkinfomc._nambig = nambig;
    }

    void fillTrkInfoMCStep(const KalSeedMC& kseedmc, TrkInfoMCStep& trkinfomcstep, const PrimaryParticle& primary) {

      const auto& genParticle = primary.primary();

      trkinfomcstep._time = genParticle.time();
      trkinfomcstep._mom = std::sqrt(genParticle.momentum().px()*genParticle.momentum().px() + 
				     genParticle.momentum().py()*genParticle.momentum().py() + 
				     genParticle.momentum().pz()*genParticle.momentum().pz());
      GeomHandle<DetectorSystem> det;
      trkinfomcstep._pos = det->toDetector(genParticle.position());

      GeomHandle<BFieldManager> bfmgr;
      GlobalConstantsHandle<ParticleDataTable> pdt;      
      static CLHEP::Hep3Vector vpoint_mu2e = det->toMu2e(CLHEP::Hep3Vector(0.0,0.0,0.0));
      static double bz = bfmgr->getBField(vpoint_mu2e).z();

      CLHEP::HepVector parvec(5,0);
      double hflt(0.0);
      HepPoint ppos(trkinfomcstep._pos._x, trkinfomcstep._pos._y, trkinfomcstep._pos._z);
      CLHEP::Hep3Vector mom = genParticle.momentum().vect();
      double charge = pdt->particle(genParticle.pdgId()).ref().charge();
      TrkHelixUtils::helixFromMom( parvec, hflt,ppos, mom,charge,bz);
      trkinfomcstep._hpar = helixpar(parvec);
    }

    void fillTrkInfoMCStep(const KalSeedMC& kseedmc, TrkInfoMCStep& trkinfomcstep,
      std::vector<int> const& vids) {

      GeomHandle<BFieldManager> bfmgr;
      GeomHandle<DetectorSystem> det;
      GlobalConstantsHandle<ParticleDataTable> pdt;
      static CLHEP::Hep3Vector vpoint_mu2e = det->toMu2e(CLHEP::Hep3Vector(0.0,0.0,0.0));
      static double bz = bfmgr->getBField(vpoint_mu2e).z();

      const auto& mcsteps = kseedmc._vdsteps;
      for (const auto& i_mcstep : mcsteps) {
	for(auto vid : vids) {
	  if (i_mcstep._vdid == vid) {
	    trkinfomcstep._time = i_mcstep._time;
	    trkinfomcstep._mom = std::sqrt(i_mcstep._mom.mag2());
	    trkinfomcstep._pos = Geom::Hep3Vec(i_mcstep._pos);

	    CLHEP::HepVector parvec(5,0);
	    double hflt(0.0);
	    HepPoint ppos(trkinfomcstep._pos._x, trkinfomcstep._pos._y, trkinfomcstep._pos._z);
	    CLHEP::Hep3Vector mom = Geom::Hep3Vec(i_mcstep._mom);
	    double charge = pdt->particle(kseedmc.simParticle()._pdg).ref().charge();
	    TrkHelixUtils::helixFromMom( parvec, hflt,ppos, mom,charge,bz);
	    trkinfomcstep._hpar = helixpar(parvec);

	    break; // only do one step, don't want to keep overwriting
	  }
	}
      }
    }

    void fillHitInfoMCs(const KalSeedMC& kseedmc, std::vector<TrkStrawHitInfoMC>& tshinfomcs) {
      tshinfomcs.clear();

      for(const auto& i_tshmc : kseedmc._tshmcs) {
      	TrkStrawHitInfoMC tshinfomc;
	fillHitInfoMC(kseedmc, tshinfomc, i_tshmc);
      	tshinfomcs.push_back(tshinfomc);
      }
    }

    void fillHitInfoMC(const KalSeedMC& kseedmc, TrkStrawHitInfoMC& tshinfomc, const TrkStrawHitMC& tshmc) {
      const Tracker& tracker = *GeomHandle<Tracker>();

      const SimPartStub& simPart = kseedmc.simParticle(tshmc._spindex);
      tshinfomc._pdg = simPart._pdg;
      tshinfomc._proc = simPart._proc;
      tshinfomc._gen = simPart._gid.id();
      tshinfomc._rel = simPart._rel.relationship();
      tshinfomc._t0 = tshmc._time;
      tshinfomc._edep = tshmc._energySum;
      tshinfomc._mom = std::sqrt(tshmc._mom.mag2());
      tshinfomc._cpos  = tshmc._cpos; 
	
      // find the step midpoint
      const Straw& straw = tracker.getStraw(tshmc._strawId);
      CLHEP::Hep3Vector mcsep = Geom::Hep3Vec(tshmc._cpos)-straw.getMidPoint();
      tshinfomc._len = mcsep.dot(straw.getDirection());
      CLHEP::Hep3Vector mdir = Geom::Hep3Vec(tshmc._mom.unit());
      CLHEP::Hep3Vector mcperp = (mdir.cross(straw.getDirection())).unit();
      double dperp = mcperp.dot(mcsep);
      tshinfomc._twdot = mdir.dot(straw.getDirection());
      tshinfomc._dist = fabs(dperp);
      tshinfomc._ambig = dperp > 0 ? -1 : 1; // follow TrkPoca convention
      // use 2-line POCA here
      TwoLinePCA pca(Geom::Hep3Vec(tshmc._cpos),mdir,straw.getMidPoint(),straw.getDirection());
      tshinfomc._doca = pca.dca();
    }

    void fillCaloClusterInfoMC(CaloClusterMC const& ccmc, CaloClusterInfoMC& ccimc) {
      ccimc._nsim = ccmc.energyDeposits().size();
      ccimc._etot = ccmc.totalEnergyDeposit();
      ccimc._tavg = ccmc.averageTime();
      if(ccmc.energyDeposits().size() > 0){
	auto const& primary = ccmc.energyDeposits().front();
	ccimc._eprimary = primary.energyDeposit();
	ccimc._tprimary = primary.time();
	ccimc._prel = primary._rel.relationship();
      }
    }
  }
}
