///////////////////////////////////////////////////////////////////////////////
#include "CalPatRec/inc/ObjectDumpUtils.hh"

#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Selector.h"

#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

#include "TrackerGeom/inc/Tracker.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"

#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "RecoDataProducts/inc/CaloProtoCluster.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/TrackSeed.hh"
#include "RecoDataProducts/inc/TrackSeedCollection.hh"

#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StrawGasStep.hh"
#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfo.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"

#include "BTrkData/inc/TrkStrawHit.hh"

#include "RecoDataProducts/inc/KalRepPtrCollection.hh"

#include "RecoDataProducts/inc/TrkCaloIntersectCollection.hh"

// #include "CalPatRec/inc/CalTimePeak.hh"


//BaBar includes
#include "BTrk/BbrGeom/TrkLineTraj.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/KalmanTrack/KalHit.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/TrkBase/HelixParams.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"

namespace mu2e {

  const SimParticleTimeOffset           *ObjectDumpUtils::_TimeOffsets(NULL);
  const StrawDigiMCCollection           *ObjectDumpUtils::_ListOfMCStrawHits(NULL);
  std::string                            ObjectDumpUtils::_FlagBgrHitsModuleLabel;

// //______________________________________________________________________________
// ObjectDumpUtils::ObjectDumpUtils() {

//   _FlagBgrHitsModuleLabel = "FlagBkgHits";

//   std::vector<std::string> VS;
//   VS.push_back(std::string("protonTimeMap"));
//   VS.push_back(std::string("muonTimeMap"));
  
//   fhicl::ParameterSet  pset;
//   pset.put("inputs", VS);
//   fgTimeOffsets = new mu2e::SimParticleTimeOffset(pset);

// }

// //______________________________________________________________________________
// ObjectDumpUtils::~ObjectDumpUtils() {
// }


//-----------------------------------------------------------------------------
void ObjectDumpUtils::printEventHeader(const art::Event* Event, const char* Message) {

  printf(" Run / Subrun / Event : %10i / %10i / %10i : %s\n",
	 Event->run(),
	 Event->subRun(),
	 Event->event(),
	 Message);
}


//-----------------------------------------------------------------------------
void ObjectDumpUtils::printCaloProtoCluster(const mu2e::CaloProtoCluster* Cluster, const char* Opt) {

  TString opt = Opt;

  int section_id(-1), iz, ir;

  const mu2e::Calorimeter       * cal(NULL);
  const mu2e::Crystal           *cr;
  const CLHEP::Hep3Vector       *pos;

  art::ServiceHandle<mu2e::GeometryService> geom;
  mu2e::GeomHandle  <mu2e::Calorimeter>     cg;

  cal = cg.get();

  if ((opt == "") || (opt == "banner")) {
    printf("-----------------------------------------------------------------------------------------------------\n");
    printf("       Address  SectionID  IsSplit  NC    Time    Energy      \n");
    printf("-----------------------------------------------------------------------------------------------------\n");
  }

  const mu2e::CaloProtoCluster::CaloCrystalHitPtrVector caloClusterHits = Cluster->caloCrystalHitsPtrVector();
  int nh = caloClusterHits.size();

  if ((opt == "") || (opt.Index("data") >= 0)) {

    printf("%16p  %3i %5i %5i %10.3f %10.3f\n",
	   Cluster,
	   section_id,
	   nh,
	   Cluster->isSplit(),
	   Cluster->time(),
	   Cluster->energyDep()
	   );
  }

  if (opt.Index("hits") >= 0) {
//-----------------------------------------------------------------------------
// print individual crystals in local disk coordinate system
//-----------------------------------------------------------------------------
    for (int i=0; i<nh; i++) {
      const mu2e::CaloCrystalHit* hit = &(*caloClusterHits.at(i));
      int id = hit->id();

      //      pos = cg->crystalOriginInSection(id);

      cr  = &cal->crystal(id);
      pos = &cr->localPosition();

      iz = -1;
      ir = -1;

      printf("%6i     %10.3f %5i %5i %8.3f %10.3f %10.3f %10.3f %10.3f\n",
	     id,
	     hit->time(),
	     iz,ir,
	     hit->energyDep(),
	     pos->x(),
	     pos->y(),
	     pos->z(),
	     hit->energyDepTot()
	     );
    }
  }
}


//-----------------------------------------------------------------------------
  void ObjectDumpUtils::printCaloProtoClusterCollection(const mu2e::CaloProtoClusterCollection* Coll) {

  const mu2e::CaloProtoCluster                 *cluster;

  int banner_printed(0), nclusters;

  nclusters = Coll->size();

  for (int i=0; i<nclusters; i++) {
    cluster = &Coll->at(i);
    if (banner_printed == 0) {
      printCaloProtoCluster(cluster, "banner");
      banner_printed = 1;
    }
    printCaloProtoCluster(cluster,"data");
  }
}

//-----------------------------------------------------------------------------
void ObjectDumpUtils::printKalRep(const KalRep* Krep, const char* Opt, const char* Prefix) {

  std::string opt = Opt;

  if ((opt == "") || (opt == "banner")) {
    printf("-----------------------------------------------------------------------------------------------");
    printf("-----------------------------------------------------\n");
    printf("%s  Address        TrkID     N  NA      P      T0      MomErr   T0Err    pT      costh    Omega",Prefix);
    printf("      D0      Z0      Phi0   TanDip    Chi2    FCons\n");
    printf("-----------------------------------------------------------------------------------------------");
    printf("-----------------------------------------------------\n");
  }

  if ((opt == "") || (opt.find("data") >= 0)) {
    double chi2   = Krep->chisq();

    int    nhits(0);

    const TrkHitVector* hits = &Krep->hitVector();
    for (auto ih=hits->begin(); ih != hits->end(); ++ih) {
      nhits++;
    }

    int    nact   = Krep->nActive();
    double t0     = Krep->t0().t0();
    double t0err  = Krep->t0().t0Err();
//-----------------------------------------------------------------------------
// in all cases define momentum at lowest Z - ideally, at the tracker entrance
//-----------------------------------------------------------------------------
    double s1     = Krep->firstHit()->kalHit()->hit()->fltLen();
    double s2     = Krep->lastHit ()->kalHit()->hit()->fltLen();
    double s      = std::min(s1,s2);

    double d0     = Krep->helix(s).d0();
    double z0     = Krep->helix(s).z0();
    double phi0   = Krep->helix(s).phi0();
    double omega  = Krep->helix(s).omega();
    double tandip = Krep->helix(s).tanDip();


    CLHEP::Hep3Vector fitmom = Krep->momentum(s);
    CLHEP::Hep3Vector momdir = fitmom.unit();
    BbrVectorErr      momerr = Krep->momentumErr(s);

    HepVector momvec(3);
    for (int i=0; i<3; i++) momvec[i] = momdir[i];

    double sigp = sqrt(momerr.covMatrix().similarity(momvec));

    double fit_consistency = Krep->chisqConsistency().consistency();
    int q         = Krep->charge();

    Hep3Vector trk_mom;

    trk_mom       = Krep->momentum(s);
    double mom    = trk_mom.mag();
    double pt     = trk_mom.perp();
    double costh  = trk_mom.cosTheta();

    char form[100];
    int nc = strlen(Prefix);
    for (int i=0; i<nc; i++) form[i]=' ';
    form[nc] = 0;
    printf("%s",form);

    printf("  %-16p %3i   %3i %3i %8.3f %8.3f %8.4f %7.4f %7.3f %8.4f",
	   Krep,
	   -1,
	   nhits,
	   nact,
	   q*mom,t0,sigp,t0err,pt,costh
	   );

    printf(" %8.5f %8.3f %8.3f %8.4f %7.4f",
	   omega,d0,z0,phi0,tandip
	   );
    printf(" %8.3f %10.3e\n",
	   chi2,
	   fit_consistency);
  }

  if (opt.find("hits") >= 0) {
//-----------------------------------------------------------------------------
// print detailed information about the track hits
//-----------------------------------------------------------------------------
    const TrkHitVector* hits = &Krep->hitVector();

    printf("--------------------------------------------------------------------");
    printf("---------------------------------------------------------------");
    printf("------------------------------------------------------\n");
    //    printf(" ih  SInd U A     len         x        y        z      HitT    HitDt");
    printf(" ih  SInd A     len         x        y        z      HitT    HitDt");
    printf(" Ch Pl  L  W     T0       Xs      Ys        Zs     resid sigres");
    printf("    Rdrift   mcdoca totErr hitErr  t0Err penErr extErr\n");
    printf("--------------------------------------------------------------------");
    printf("---------------------------------------------------------------");
    printf("------------------------------------------------------\n");

    mu2e::TrkStrawHit* hit;
    CLHEP::Hep3Vector  pos;
    int i = 0;

    for (auto it=hits->begin(); it!=hits->end(); it++) {
      // TrkStrawHit inherits from TrkHitOnTrk

      hit = (mu2e::TrkStrawHit*) (*it);

      const mu2e::ComboHit* sh = &hit->comboHit();
      mu2e::Straw*   straw = (mu2e::Straw*) &hit->straw();

      hit->hitPosition(pos);

      double    len  = hit->fltLen();
      HepPoint  plen = Krep->position(len);
//-----------------------------------------------------------------------------
// find MC truth DOCA in a given straw
// start from finding the right vector of StepPointMC's
//-----------------------------------------------------------------------------
      int vol_id, nstraws(0);

      if (_ListOfMCStrawHits) nstraws = _ListOfMCStrawHits->size();

      const mu2e::StrawGasStep* step(0);

      for (int i=0; i<nstraws; i++) {

	const mu2e::StrawDigiMC* mcdigi = &_ListOfMCStrawHits->at(i);

	const mu2e::StrawGasStep   *step;
	if (mcdigi->wireEndTime(mu2e::StrawEnd::cal) < mcdigi->wireEndTime(mu2e::StrawEnd::hv)) {
	  step = mcdigi->strawGasStep(mu2e::StrawEnd::cal).get();
	}
	else {
	  step = mcdigi->strawGasStep(mu2e::StrawEnd::hv ).get();
	}

//	vol_id = step->volumeId();
	vol_id = step->strawId().asUint16();
 	if (vol_id == straw->id().asUint16()) {
 					// step found - use the first one in the straw
 	  break;
 	}
      }

      double mcdoca = -99.0;

      if (step) {
	const Hep3Vector* v1 = &straw->getMidPoint();
	HepPoint p1(v1->x(),v1->y(),v1->z());

	Hep3Vector v2 = step->position();
	HepPoint    p2(v2.x(),v2.y(),v2.z());

	TrkLineTraj trstraw(p1,straw->getDirection()  ,0.,0.);
	TrkLineTraj trstep (p2,Geom::Hep3Vec(step->momentum().unit()),0.,0.);

	TrkPoca poca(trstep, 0., trstraw, 0.);

	mcdoca = poca.doca();
      }

      //      printf("%3i %5i %1i %1i %9.3f %8.3f %8.3f %9.3f %8.3f %7.3f",
      printf("%3i %5i %1i %9.3f %8.3f %8.3f %9.3f %8.3f %7.3f",
	     ++i,
	     straw->id().asUint16(),
	     //	     hit->isUsable(),
	     hit->isActive(),
	     len,
	     //	     hit->hitRms(),
	     plen.x(),plen.y(),plen.z(),
	     sh->time(), 0.//sh->dt()//FIXME!
	     );

      printf(" %2i %2i %2i %2i",
	     straw->id().getPlane(),
	     straw->id().getPanel(),
	     straw->id().getLayer(),
	     straw->id().getStraw()
	     );

      printf(" %8.3f",hit->hitT0().t0());

      double res, sigres;
      hit->resid(res, sigres, true);

      printf("%8.3f %8.3f %9.3f %7.3f %7.3f",
	     pos.x(),
	     pos.y(),
	     pos.z(),
	     res,
	     sigres
	     );

      if (hit->isActive()) {
	if      (hit->ambig()       == 0) printf(" * %6.3f",hit->driftRadius());
	else if (hit->ambig()*mcdoca > 0) printf("   %6.3f",hit->driftRadius()*hit->ambig());
	else                              printf(" ? %6.3f",hit->driftRadius()*hit->ambig());
      }
      else {
//-----------------------------------------------------------------------------
// do not analyze correctness of the drift sign determination for hits not
// marked as 'active'
//-----------------------------------------------------------------------------
	printf("   %6.3f",hit->driftRadius());
      }


      printf("  %7.3f",mcdoca);
      printf(" %6.3f %6.3f %6.3f %6.3f %6.3f",
	     hit->totalErr(),
	     hit->hitErr(),
	     hit->t0Err(),
	     hit->penaltyErr(),
	     hit->temperature()
	     );
//-----------------------------------------------------------------------------
// test: calculated residual in fTmp[0]
//-----------------------------------------------------------------------------
//       Test_000(Krep,hit);
//       printf(" %7.3f",fTmp[0]);

      printf("\n");
    }
  }
}


//-----------------------------------------------------------------------------
void ObjectDumpUtils::printKalRepCollection(const art::Event* Event        ,
					    const KalRepPtrCollection* Coll,
					    int               PrintHits    ) {

  art::Handle<mu2e::StrawDigiMCCollection> mcdigiH;

  Event->getByLabel("makeSD","",mcdigiH);
  if (mcdigiH.isValid()) {
    _ListOfMCStrawHits = (mu2e::StrawDigiMCCollection*) mcdigiH.product();
  }
  else {
    _ListOfMCStrawHits = NULL;
    printf(">>> ERROR in ObjectDumpUtils::printKalRepCollection: failed to locate StrawDigiMCCollection by makeSD\n");
  }

  int ntrk = Coll->size();

  const KalRep *trk;

  //  int banner_printed = 0;
  for (int i=0; i<ntrk; i++) {
    art::Ptr<KalRep> kptr = Coll->at(i);
//     Event->get(kptr.id(), krepsHandle);
//     fhicl::ParameterSet const& pset = krepsHandle.provenance()->parameterSet();
//     string module_type = pset.get<std::string>("module_type");

    trk = kptr.get();
    printKalRep(trk,"banner+data+hits",""); // module_type.data());
  }

}


}
