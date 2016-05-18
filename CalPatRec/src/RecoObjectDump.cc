///////////////////////////////////////////////////////////////////////////////
#include "CalPatRec/inc/RecoObjectDump.hh"

#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Selector.h"

#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

#include "TTrackerGeom/inc/TTracker.hh"
#include "CalorimeterGeom/inc/VaneCalorimeter.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"

#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "RecoDataProducts/inc/CaloProtoCluster.hh"
#include "RecoDataProducts/inc/CaloProtoClusterCollection.hh"
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
#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruth.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfo.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"

#include "BTrkData/inc/TrkStrawHit.hh"
#include "RecoDataProducts/inc/KalRepPtrCollection.hh"

#include "TrackCaloMatching/inc/TrkToCaloExtrapolCollection.hh"
#include "RecoDataProducts/inc/TrkCaloIntersectCollection.hh"
#include "TrackCaloMatching/inc/TrackClusterMatch.hh"

#include "CalPatRec/inc/CalTimePeak.hh"

#include "Stntuple/base/TNamedHandle.hh"


//BaBar includes
#include "BTrk/BbrGeom/TrkLineTraj.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/KalmanTrack/KalHit.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/TrkBase/HelixParams.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"

namespace mu2e {

  const SimParticleTimeOffset           *RecoObjectDump::_TimeOffsets(NULL);
  const PtrStepPointMCVectorCollection  *RecoObjectDump::_ListOfMCStrawHits(NULL);
  std::string                            RecoObjectDump::_FlagBgrHitsModuleLabel;

// //______________________________________________________________________________
// RecoObjectDump::RecoObjectDump() {

//   _FlagBgrHitsModuleLabel = "FlagBkgHits";

//   std::vector<std::string> VS;
//   VS.push_back(std::string("protonTimeMap"));
//   VS.push_back(std::string("muonTimeMap"));
  
//   fhicl::ParameterSet  pset;
//   pset.put("inputs", VS);
//   fgTimeOffsets = new mu2e::SimParticleTimeOffset(pset);

// }

// //______________________________________________________________________________
// RecoObjectDump::~RecoObjectDump() {
// }


//-----------------------------------------------------------------------------
  void RecoObjectDump::printEventHeader(const art::Event* Event, const char* Message) {

  printf(" Run / Subrun / Event : %10i / %10i / %10i : %s\n",
	 Event->run(),
	 Event->subRun(),
	 Event->event(),
	 Message);
}

//-----------------------------------------------------------------------------
void RecoObjectDump::printKalRep(const KalRep* Krep, const char* Opt, const char* Prefix) {

  string opt = Opt;
  
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

      const mu2e::StrawHit* sh = &hit->strawHit();
      mu2e::Straw*   straw = (mu2e::Straw*) &hit->straw();

      hit->hitPosition(pos);

      double    len  = hit->fltLen();
      HepPoint  plen = Krep->position(len);
//-----------------------------------------------------------------------------
// find MC truth DOCA in a given straw
// start from finding the right vector of StepPointMC's
//-----------------------------------------------------------------------------
      int vol_id;
      int nstraws = _ListOfMCStrawHits->size();

      const mu2e::StepPointMC* step(0);

      for (int i=0; i<nstraws; i++) {
	mu2e::PtrStepPointMCVector  const& mcptr(_ListOfMCStrawHits->at(i));
	step = &(*mcptr.at(0));
	vol_id = step->volumeId();
 	if (vol_id == straw->index().asInt()) {
 					// step found - use the first one in the straw
 	  break;
 	}
      }

      double mcdoca = -99.0;

      if (step) {
	const Hep3Vector* v1 = &straw->getMidPoint();
	HepPoint p1(v1->x(),v1->y(),v1->z());

	const Hep3Vector* v2 = &step->position();
	HepPoint    p2(v2->x(),v2->y(),v2->z());

	TrkLineTraj trstraw(p1,straw->getDirection()  ,0.,0.);
	TrkLineTraj trstep (p2,step->momentum().unit(),0.,0.);

	TrkPoca poca(trstep, 0., trstraw, 0.);
    
	mcdoca = poca.doca();
      }

      //      printf("%3i %5i %1i %1i %9.3f %8.3f %8.3f %9.3f %8.3f %7.3f",
      printf("%3i %5i %1i %9.3f %8.3f %8.3f %9.3f %8.3f %7.3f",
	     ++i,
	     straw->index().asInt(), 
	     //	     hit->isUsable(),
	     hit->isActive(),
	     len,
	     //	     hit->hitRms(),
	     plen.x(),plen.y(),plen.z(),
	     sh->time(), sh->dt()
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
	     hit->extErr()
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
void RecoObjectDump::printKalRepCollection(const art::Event* Event         , 
					   const KalRepPtrCollection* Coll,
					   int               PrintHits     ) {

  art::Handle<mu2e::PtrStepPointMCVectorCollection> mcptrHandle;

  Event->getByLabel("makeSH","StrawHitMCPtr",mcptrHandle);
  if (mcptrHandle.isValid()) {
    _ListOfMCStrawHits = (mu2e::PtrStepPointMCVectorCollection*) mcptrHandle.product();
  }
  else {
    _ListOfMCStrawHits = NULL;
    printf(">>> ERROR in RecoObjectDump::printKalRepCollection: failed to locate StepPointMCCollection makeSH:StrawHitMCPtr\n");
  }

  int ntrk = Coll->size();

  const KalRep *trk;

  int banner_printed = 0;
  for (int i=0; i<ntrk; i++) {
    art::Ptr<KalRep> kptr = Coll->at(i);
//     Event->get(kptr.id(), krepsHandle);
//     fhicl::ParameterSet const& pset = krepsHandle.provenance()->parameterSet();
//     string module_type = pset.get<std::string>("module_type");
 
    trk = kptr.get();
    if (banner_printed == 0) {
      printKalRep(trk,"banner",""); // module_type.data());
      if (PrintHits == 0) banner_printed = 1;
    }
    printKalRep(trk,"data",""); // module_type.data());
    if (PrintHits > 0) printKalRep(trk,"hits");
  }
 
}


}

