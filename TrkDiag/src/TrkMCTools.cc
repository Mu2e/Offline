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

    void findMCTrk(const KalSeed& kseed,std::vector<spcount>& sct, StrawDigiMCCollection const& mcdigis, bool saveall) {
      sct.clear();
      // find the SimParticles which contributed hits.
      // loop through the straw hits from the track
      static StrawHitFlag active(StrawHitFlag::active);
      for(const auto& tshs : kseed.hits()) {
	// loop over the hits and find the associated steppoints
	bool isactive = tshs.flag().hasAllProperties(active);
	StrawDigiMC const& mcdigi = mcdigis.at(tshs.index());
	art::Ptr<SimParticle> spp = mcdigi.earlyStepPointMC()->simParticle();
	// see if this particle has already been found; if so, increment, if not, add it
	bool found(false);
	for(auto& spc : sct ) {
	  if(spc._spp == spp ){
	    found = true;
	    spc.append(spp,isactive);
	    break;
	  }
	}
	if(!found)sct.push_back(spcount(spp,isactive));
      }
      if(saveall){
	// add the SimParticles that contributed non-trigger energy.  These have 0 count
	for(const auto& tshs : kseed.hits()) {
	  StrawDigiMC const& mcdigi = mcdigis.at(tshs.index());
	  for(auto const& spmc : mcdigi.stepPointMCs()){
	    bool found(false);
	    for(auto& spc : sct ) {
	      if(spc._spp == spmc->simParticle() ){
		found = true;
		break;
	      }
	    }
	    if(!found)sct.push_back(spcount(spmc->simParticle()));
	  }
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

    void countDigis(const KalSeedMC& kseedmc, const KalSeed& kseed, int& ndigi, int& ndigigood, int& nambig) {
      static double mingood = 0.9; // this is clumsy!
      ndigi = 0; ndigigood = 0, nambig = 0;
      // find the first segment momentum as reference
      double simmom = 1.0;
      if(kseedmc.simParticles().size()>0)
	simmom = sqrt(kseedmc.simParticles().front()._mom.mag2());
      for(size_t i_digi = 0; i_digi < kseedmc._tshmcs.size(); ++i_digi) {
	const auto& tshmc = kseedmc._tshmcs.at(i_digi);

	if (kseedmc.simParticle(tshmc._spindex)._rel == MCRelationship::same) {
	  ++ndigi;
	  if(sqrt(tshmc.particleMomentum().mag2())/simmom > mingood)	  ++ndigigood;

	  // easiest way to get MC ambiguity is through info object
	  TrkStrawHitInfoMC tshinfomc;
	  fillHitInfoMC(kseedmc,tshinfomc,tshmc);  
	  // the MCDigi list can be longer than the # of TrkStrawHits in the seed:
	  if(i_digi < kseed.hits().size()){ 
	    const auto& ihit = kseed.hits().at(i_digi);
	    if(ihit.ambig()*tshinfomc._ambig > 0) {
	      ++nambig;
	    }
	  }
	}
      }
    }

    void primaryRelation(PrimaryParticle const& primary,
	StrawDigiMCCollection const& sdmccol, std::vector<StrawDigiIndex> const& indices,
	art::Ptr<SimParticle>& primarysim, unsigned& nprimary, MCRelationship& mcrel) {
      // reset
      primarysim = art::Ptr<SimParticle>();
      nprimary = 0;
      mcrel = MCRelationship();
      // loop over primary sim particles
      for( auto spp : primary.primarySimParticles()){
	unsigned count(0);
	art::Ptr<SimParticle> sp;
	for(auto sdi : indices) {
	// find relation of this digi to the primary
	  MCRelationship prel(spp,sdmccol.at(sdi).earlyStepPointMC()->simParticle());
	  // count the highest relationship for these digis
	  if(prel == mcrel)
	    count++;
	  else if(prel > mcrel){
	    mcrel = prel;
	    count = 1;
	    sp = sdmccol.at(sdi).earlyStepPointMC()->simParticle();
	  }
	}
	if(count > nprimary){
	  nprimary = count;
	  primarysim = sp;
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
      tshinfomc._rel = simPart._rel;
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
	ccimc._prel = primary._rel;
      }
    }
  }

  void TrkMCHelper::fillTrkInfoMC(const KalSeedMC& kseedmc, TrkInfoMC& trkinfomc) {
    // use the primary match of the track
    // primary associated SimParticle
    auto trkprimary = kseedmc.simParticle().simParticle(_spcH);
    if(kseedmc.simParticles().size() > 0){
      auto const& simp = kseedmc.simParticles().front();
      trkinfomc._gen = simp._gid.id();
      trkinfomc._pdg = simp._pdg;
      trkinfomc._proc = simp._proc;
      trkinfomc._nhits = simp._nhits; // number of hits from the primary particle
      trkinfomc._nactive = simp._nactive; // number of active hits from the primary particle
      trkinfomc._prel = simp._rel; // relationship of the track primary to the event primary
    }

    fillTrkInfoMCDigis(kseedmc, trkinfomc);

    // fill the origin information of this SimParticle
    GeomHandle<DetectorSystem> det;
    trkinfomc._otime = trkprimary->startGlobalTime() + _toff.totalTimeOffset(trkprimary);
    trkinfomc._opos = Geom::toXYZVec(det->toDetector(trkprimary->startPosition()));
    trkinfomc._omom = Geom::toXYZVec(trkprimary->startMomentum());
  }

  void TrkMCHelper::fillTrkInfoMCDigis(const KalSeedMC& kseedmc, TrkInfoMC& trkinfomc) {
    trkinfomc._ndigi = 0; trkinfomc._ndigigood = 0, trkinfomc._nambig = 0;
    // find the first segment momentum as reference
    double simmom = 1.0;
    if(kseedmc.simParticles().size()>0)
      simmom = sqrt(kseedmc.simParticles().front()._mom.mag2());
    for(size_t i_digi = 0; i_digi < kseedmc._tshmcs.size(); ++i_digi) {
      const auto& tshmc = kseedmc._tshmcs.at(i_digi);

      if (kseedmc.simParticle(tshmc._spindex)._rel == MCRelationship::same) {
	++trkinfomc._ndigi;
	if(sqrt(tshmc.particleMomentum().mag2())/simmom > _mingood) {
	  ++trkinfomc._ndigigood;
	}

	// easiest way to get MC ambiguity is through info object
	TrkStrawHitInfoMC tshinfomc;
	fillHitInfoMC(kseedmc,tshinfomc,tshmc);  
	// the MCDigi list can be longer than the # of TrkStrawHits in the seed:
	/*	if(i_digi < kseed.hits().size()){ 
	  const auto& ihit = kseed.hits().at(i_digi);
	  if(ihit.ambig()*tshinfomc._ambig > 0) {
	    ++trkinfomc._nambig; // TODO
	  }
	}
	*/
      }
    }
  }

  void TrkMCHelper::fillHitInfoMC(const KalSeedMC& kseedmc, TrkStrawHitInfoMC& tshinfomc, const TrkStrawHitMC& tshmc) {
    const Tracker& tracker = *GeomHandle<Tracker>();

    const SimPartStub& simPart = kseedmc.simParticle(tshmc._spindex);
    tshinfomc._pdg = simPart._pdg;
    tshinfomc._proc = simPart._proc;
    tshinfomc._gen = simPart._gid.id();
    tshinfomc._rel = simPart._rel;
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

  void TrkMCHelper::fillGenAndPriInfo(const KalSeedMC& kseedmc, const PrimaryParticle& primary, GenInfo& priinfo, GenInfo& geninfo) {
    auto trkprimary = kseedmc.simParticle().simParticle(_spcH);

    fillPriInfo(primary, priinfo);

    // go through the SimParticles of this primary, and find the one most related to the
    // downstream fit (KalSeedMC)
    auto bestprimarysp = primary.primarySimParticles().front();
    MCRelationship bestrel;
    for(auto const& spp : primary.primarySimParticles()){
      MCRelationship mcrel(spp,trkprimary);
      if(mcrel > bestrel){
	bestrel = mcrel;
	bestprimarysp = spp;
      }
    }
    priinfo._time = primary.primary().time() + _toff.totalTimeOffset(bestprimarysp);

    const auto& gp = bestprimarysp->genParticle();
    fillGenInfo(gp, geninfo);
    geninfo._time = gp->time() + _toff.totalTimeOffset(bestprimarysp);
  }

  void TrkMCHelper::fillPriInfo(const PrimaryParticle& primary, GenInfo& priinfo) {

    GeomHandle<DetectorSystem> det;
    // fill primary info from the primary GenParticle
    const auto& genParticle = primary.primary();
    priinfo._pdg = genParticle.pdgId();
    priinfo._gen = genParticle.generatorId().id();
    priinfo._mom = Geom::toXYZVec(genParticle.momentum());
    priinfo._pos = Geom::toXYZVec(det->toDetector(genParticle.position()));
    priinfo._time = genParticle.time(); // NB doesn't have time offsets applied
  }

  void TrkMCHelper::fillGenInfo(const art::Ptr<GenParticle>& gp, GenInfo& geninfo) {

    GeomHandle<DetectorSystem> det;

    if(gp.isNonnull()){
      geninfo._pdg = gp->pdgId();
      geninfo._gen = gp->generatorId().id();
      geninfo._mom = Geom::toXYZVec(gp->momentum());
      geninfo._pos = Geom::toXYZVec(det->toDetector(gp->position()));
      geninfo._time = gp->time(); // NB doesn't have time offsets applied
    }
  }

  void TrkMCHelper::fillTrkInfoMCStep(const KalSeedMC& kseedmc, TrkInfoMCStep& trkinfomcstep,
				      std::vector<int> const& vids) {

    GeomHandle<BFieldManager> bfmgr;
    GeomHandle<DetectorSystem> det;
    GlobalConstantsHandle<ParticleDataTable> pdt;
    static CLHEP::Hep3Vector vpoint_mu2e = det->toMu2e(CLHEP::Hep3Vector(0.0,0.0,0.0));
    static double bz = bfmgr->getBField(vpoint_mu2e).z();

    const auto& mcsteps = kseedmc._vdsteps;
    double dmin = FLT_MAX;
    for (const auto& i_mcstep : mcsteps) {
      for(auto vid : vids) {
	if (i_mcstep._vdid == vid) {
	  // take the earliest time if there are >1; this avoids picking up
	  // the albedo track from the calorimeter which can re-enter the tracker
	  if(i_mcstep._time < dmin){
	    dmin = i_mcstep._time;
	    trkinfomcstep._time = i_mcstep._time;
	    trkinfomcstep._mom = Geom::Hep3Vec(i_mcstep._mom);
	    trkinfomcstep._pos = Geom::Hep3Vec(i_mcstep._pos);

	    CLHEP::HepVector parvec(5,0);
	    double hflt(0.0);
	    HepPoint ppos(trkinfomcstep._pos.x(), trkinfomcstep._pos.y(), trkinfomcstep._pos.z());
	    CLHEP::Hep3Vector mom = Geom::Hep3Vec(i_mcstep._mom);
	    double charge = pdt->particle(kseedmc.simParticle()._pdg).ref().charge();
	    TrkHelixUtils::helixFromMom( parvec, hflt,ppos, mom,charge,bz);
	    trkinfomcstep._hpar = helixpar(parvec);
	  }
	}
      }
    }
  }
}
