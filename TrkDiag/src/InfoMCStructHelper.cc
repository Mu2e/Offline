//
// Namespace for collecting tools used in MC truth evaluation
// Original author: Dave Brown (LBNL) 8/10/2016
//
#include "TrkDiag/inc/InfoMCStructHelper.hh"
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
  void InfoMCStructHelper::fillTrkInfoMC(const KalSeedMC& kseedmc, TrkInfoMC& trkinfomc) {
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

  void InfoMCStructHelper::fillTrkInfoMCDigis(const KalSeedMC& kseedmc, TrkInfoMC& trkinfomc) {
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

  void InfoMCStructHelper::fillHitInfoMC(const KalSeedMC& kseedmc, TrkStrawHitInfoMC& tshinfomc, const TrkStrawHitMC& tshmc) {
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

  void InfoMCStructHelper::fillGenAndPriInfo(const KalSeedMC& kseedmc, const PrimaryParticle& primary, GenInfo& priinfo, GenInfo& geninfo) {
    auto trkprimary = kseedmc.simParticle().simParticle(_spcH);

    fillPriInfo(primary, priinfo);

    // go through the SimParticles of this primary, and find the one most related to the
    // downstream fit (KalSeedMC)

    if (primary.primarySimParticles().empty()) {
      throw cet::exception("Simulation")<<"InfoMCStructHelper: No Primary Particle found" << std::endl;
    }
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

  void InfoMCStructHelper::fillPriInfo(const PrimaryParticle& primary, GenInfo& priinfo) {

    GeomHandle<DetectorSystem> det;
    // fill primary info from the primary GenParticle
    const auto& genParticle = primary.primary();
    priinfo._pdg = genParticle.pdgId();
    priinfo._gen = genParticle.generatorId().id();
    priinfo._mom = Geom::toXYZVec(genParticle.momentum());
    priinfo._pos = Geom::toXYZVec(det->toDetector(genParticle.position()));
    priinfo._time = genParticle.time(); // NB doesn't have time offsets applied
  }

  void InfoMCStructHelper::fillGenInfo(const art::Ptr<GenParticle>& gp, GenInfo& geninfo) {

    GeomHandle<DetectorSystem> det;

    if(gp.isNonnull()){
      geninfo._pdg = gp->pdgId();
      geninfo._gen = gp->generatorId().id();
      geninfo._mom = Geom::toXYZVec(gp->momentum());
      geninfo._pos = Geom::toXYZVec(det->toDetector(gp->position()));
      geninfo._time = gp->time(); // NB doesn't have time offsets applied
    }
  }

  void InfoMCStructHelper::fillTrkInfoMCStep(const KalSeedMC& kseedmc, TrkInfoMCStep& trkinfomcstep,
					     std::vector<int> const& vids, double target_time) {

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
	  // take the VD step with the time closest to target_time
	  // this is so that we can take the correct step when looking at upstream/downstream trachs
	  double corrected_time = fmod(i_mcstep._time, 1695); // VDStep is created with the time offsets included
	  if(fabs(target_time - corrected_time) < dmin){
	    dmin = fabs(target_time - corrected_time);//i_mcstep._time;
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

  void InfoMCStructHelper::fillHitInfoMCs(const KalSeedMC& kseedmc, std::vector<TrkStrawHitInfoMC>& tshinfomcs) {
    tshinfomcs.clear();

    for(const auto& i_tshmc : kseedmc._tshmcs) {
      TrkStrawHitInfoMC tshinfomc;
      fillHitInfoMC(kseedmc, tshinfomc, i_tshmc);
      tshinfomcs.push_back(tshinfomc);
    }
  }

  void InfoMCStructHelper::fillCaloClusterInfoMC(CaloClusterMC const& ccmc, CaloClusterInfoMC& ccimc) {
    auto const& edeps = ccmc.energyDeposits();
    ccimc._nsim = edeps.size();
    ccimc._etot = ccmc.totalEnergyDep();
    if(ccmc.energyDeposits().size() > 0){
      auto const& primary = edeps.front();
      ccimc._eprimary = primary.energyDep();
      ccimc._tavg = primary.time(); // this is unnecessary FIXMI
      ccimc._tprimary = primary.time();
      ccimc._prel = primary.rel();
    }
  }
}
