//

#include "mod/TAnaDump.hh"
#include "TROOT.h"

#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
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

#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruth.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfo.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"

#include "KalmanTests/inc/TrkStrawHit.hh"
#include "KalmanTests/inc/KalRepPtrCollection.hh"

#include "TrackCaloMatching/inc/TrkToCaloExtrapolCollection.hh"
#include "TrackCaloMatching/inc/TrackClusterMatch.hh"

#include "CalPatRec/inc/CalTimePeak.hh"

#include "Stntuple/base/TNamedHandle.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"


//BaBar includes
#include "TrajGeom/TrkLineTraj.hh"
#include "TrkBase/TrkPoca.hh"

ClassImp(TAnaDump)

TAnaDump* TAnaDump::fgInstance = 0;

namespace {
  mu2e::SimParticleTimeOffset           *fgTimeOffsets(NULL);
  mu2e::PtrStepPointMCVectorCollection  *fgListOfMCStrawHits;
}
//______________________________________________________________________________
TAnaDump::TAnaDump() {
//   if (! TROOT::Initialized()) {
//     static TROOT a ("ROOT@Mu2e","hmm",initfuncs);
//   }
  fEvent = 0;
  fListOfObjects          = new TObjArray();
  fFlagBgrHitsModuleLabel = "FlagBkgHits";

  std::vector<std::string> VS;
  VS.push_back(std::string("protonTimeMap"));
  VS.push_back(std::string("muonTimeMap"));
  
  fhicl::ParameterSet  pset;
  pset.put("inputs", VS);
  fgTimeOffsets = new mu2e::SimParticleTimeOffset(pset);

}

//------------------------------------------------------------------------------
TAnaDump* TAnaDump::Instance() {
  static TAnaDump::Cleaner cleaner;

  if  (!  fgInstance) fgInstance  = new TAnaDump();
  return fgInstance;
}



//______________________________________________________________________________
TAnaDump::~TAnaDump() {
  fListOfObjects->Delete();
  delete fListOfObjects;
}

//------------------------------------------------------------------------------
TAnaDump::Cleaner::Cleaner() {
}

//------------------------------------------------------------------------------
  TAnaDump::Cleaner::~Cleaner() {
    if (TAnaDump::fgInstance) {
      delete TAnaDump::fgInstance;
      TAnaDump::fgInstance = 0;
    }
  }


//-----------------------------------------------------------------------------
void TAnaDump::AddObject(const char* Name, void* Object) {
  TNamedHandle* h = new TNamedHandle(Name,Object);
  fListOfObjects->Add(h);
}

//-----------------------------------------------------------------------------
void* TAnaDump::FindObject(const char* Name) {
  void* o(NULL);

  TNamedHandle* h = (TNamedHandle*) fListOfObjects->FindObject(Name);
  if (h != NULL) {
    o = h->Object();
  }
  return o;
}

//-----------------------------------------------------------------------------
// print position of the cluster in the tracker system
//-----------------------------------------------------------------------------
void TAnaDump::printCaloCluster(const mu2e::CaloCluster* Cl, const char* Opt) {
  int row, col;
  TString opt = Opt;

  art::ServiceHandle<mu2e::GeometryService> geom;
  mu2e::GeomHandle  <mu2e::Calorimeter>     cal;
  Hep3Vector        gpos, tpos;

  if ((opt == "") || (opt == "banner")) {
    printf("-----------------------------------------------------------------------------------------------");
    printf("-------------------------------\n");
    printf("       Address  VaneID  Parent  NC       Time    Row   Col   Energy  X(loc)     Y(loc)   Z(loc)");
    printf("        X          Y          Z\n");
    printf("-----------------------------------------------------------------------------------------------");
    printf("-------------------------------\n");
  }
 
  if ((opt == "") || (opt.Index("data") >= 0)) {
    row = -1; // Cl->cogRow();
    col = -1; // Cl->cogColumn();
    
    const mu2e::CaloCluster::CaloCrystalHitPtrVector caloClusterHits = Cl->caloCrystalHitsPtrVector();
    int nh = caloClusterHits.size();

    if ((row < 0) || (row > 9999)) row = -9999;
    if ((col < 0) || (col > 9999)) col = -9999;
//-----------------------------------------------------------------------------
// transform cluster coordinates to the tracker coordiante system
//-----------------------------------------------------------------------------
    gpos = cal->fromSectionFrameFF(Cl->sectionId(),Cl->cog3Vector());
    tpos = cal->toTrackerFrame(gpos);

    printf("%16p %5i %7i %3i  %10.3f %5i %5i %8.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n",
	   Cl,
	   Cl->sectionId(),
	   -999, 
	   nh,
	   Cl->time(),
	   row, col,
	   Cl->energyDep(),
	   Cl->cog3Vector().x(),
	   Cl->cog3Vector().y(),
	   Cl->cog3Vector().z(),
	   tpos.x(),
	   tpos.y(),
	   tpos.z() 
	   );
  }
  
  if (opt.Index("hits") >= 0) {

    Hep3Vector pos;
    int  iz, ir;
//-----------------------------------------------------------------------------
// print individual crystals in local vane coordinate system
//-----------------------------------------------------------------------------
    const mu2e::CaloCluster::CaloCrystalHitPtrVector caloClusterHits = Cl->caloCrystalHitsPtrVector();
    int nh = caloClusterHits.size();

    for (int i=0; i<nh; i++) {
      const mu2e::CaloCrystalHit* hit = &(*caloClusterHits.at(i));
      int id = hit->id();
      
      pos = cal->crystalOriginInSection(id);

      if(geom->hasElement<mu2e::VaneCalorimeter>() ){
	mu2e::GeomHandle<mu2e::VaneCalorimeter> cgvane;
	iz  = cgvane->nCrystalX();
	ir  = cgvane->nCrystalY();
      }
      else {
	iz = -1;
	ir = -1;
      }

      printf("%6i     %10.3f %5i %5i %8.3f %10.3f %10.3f %10.3f %10.3f\n",
	     id,
	     hit->time(),
	     iz,ir,
	     hit->energyDep(),
	     pos.x(),
	     pos.y(),
	     pos.z(),
	     hit->energyDepTotal()
	     );
    }
  }
}


//-----------------------------------------------------------------------------
void TAnaDump::printCaloClusterCollection(const char* ModuleLabel, 
					  const char* ProductName,
					  const char* ProcessName) {

  printf(">>>> ModuleLabel = %s\n",ModuleLabel);

  //data about hits in the calorimeter crystals

  art::Handle<mu2e::CaloClusterCollection> handle;
  const mu2e::CaloClusterCollection* caloCluster;

  if (ProductName[0] != 0) {
    art::Selector  selector(art::ProductInstanceNameSelector(ProductName) &&
			    art::ProcessNameSelector(ProcessName)         && 
			    art::ModuleLabelSelector(ModuleLabel)            );
    fEvent->get(selector,handle);
  }
  else {
    art::Selector  selector(art::ProcessNameSelector(ProcessName)         && 
			    art::ModuleLabelSelector(ModuleLabel)            );
    fEvent->get(selector,handle);
  }
//-----------------------------------------------------------------------------
// make sure collection exists
//-----------------------------------------------------------------------------
  if (! handle.isValid()) {
    printf("TAnaDump::printCaloClusterCollection: no CaloClusterCollection ");
    printf("for module %s and ProductName=%s found, BAIL OUT\n",
	   ModuleLabel,ProductName);
    return;
  }

  caloCluster = handle.product();

  int nhits = caloCluster->size();

  const mu2e::CaloCluster* hit;


  int banner_printed = 0;
  for (int i=0; i<nhits; i++) {
    hit = &caloCluster->at(i);
    if (banner_printed == 0) {
      printCaloCluster(hit, "banner");
      banner_printed = 1;
    }
    printCaloCluster(hit,"data");
  }
 
}


//-----------------------------------------------------------------------------
void TAnaDump::printCaloProtoCluster(const mu2e::CaloProtoCluster* Cluster, const char* Opt) {

  TString opt = Opt;

  int section_id(-1), iz, ir;
  Hep3Vector pos;

  art::ServiceHandle<mu2e::GeometryService> geom;
  mu2e::GeomHandle  <mu2e::Calorimeter>     cg;

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
// print individual crystals in local vane coordinate system
//-----------------------------------------------------------------------------
    for (int i=0; i<nh; i++) {
      const mu2e::CaloCrystalHit* hit = &(*caloClusterHits.at(i));
      int id = hit->id();
      
      pos = cg->crystalOriginInSection(id);

      if (geom->hasElement<mu2e::VaneCalorimeter>()) {
	mu2e::GeomHandle<mu2e::VaneCalorimeter> cgvane;
	iz  = cgvane->nCrystalX();
	ir  = cgvane->nCrystalY();
      }
      else {
	iz = -1;
	ir = -1;
      }
      
      printf("%6i     %10.3f %5i %5i %8.3f %10.3f %10.3f %10.3f %10.3f\n",
	     id,
	     hit->time(),
	     iz,ir,
	     hit->energyDep(),
	     pos.x(),
	     pos.y(),
	     pos.z(),
	     hit->energyDepTotal()
	     );
    }
  }
}


//-----------------------------------------------------------------------------
void TAnaDump::printCaloProtoClusterCollection(const char* ModuleLabel, 
					       const char* ProductName,
					       const char* ProcessName) {

  art::Handle<mu2e::CaloProtoClusterCollection> handle;
  const mu2e::CaloProtoClusterCollection       *coll;
  const mu2e::CaloProtoCluster                 *cluster;

  int banner_printed(0), nclusters;

  if (ProductName[0] != 0) {
    art::Selector  selector(art::ProductInstanceNameSelector(ProductName) &&
			    art::ProcessNameSelector(ProcessName)         && 
			    art::ModuleLabelSelector(ModuleLabel)            );
    fEvent->get(selector,handle);
  }
  else {
    art::Selector  selector(art::ProcessNameSelector(ProcessName)         && 
			    art::ModuleLabelSelector(ModuleLabel)            );
    fEvent->get(selector,handle);
  }
//-----------------------------------------------------------------------------
// make sure collection exists
//-----------------------------------------------------------------------------
  if (! handle.isValid()) {
    printf("TAnaDump::printCaloProtoClusterCollection: no CaloProtoClusterCollection ");
    printf("for module %s and ProductName=%s found, BAIL OUT\n",
	   ModuleLabel,ProductName);
    return;
  }

  coll      = handle.product();
  nclusters = coll->size();

  for (int i=0; i<nclusters; i++) {
    cluster = &coll->at(i);
    if (banner_printed == 0) {
      printCaloProtoCluster(cluster, "banner");
      banner_printed = 1;
    }
    printCaloProtoCluster(cluster,"data");
  }
}

//-----------------------------------------------------------------------------
void TAnaDump::printCalTimePeak(const mu2e::CalTimePeak* TPeak, const char* Opt) {

  const mu2e::StrawHit*      hit;
  int                        flags;
  const mu2e::StepPointMC*   step(NULL);

  TString opt = Opt;
  opt.ToLower();

  if ((opt == "") || (opt.Index("banner") >= 0)) {
    printf("------------------------------------------------\n");
    printf("    Energy CalPatRec      Z         T0    Nhits \n");
    printf("------------------------------------------------\n");
  }
 
  if ((opt == "") || (opt.Index("data") >= 0)) {
    printf("%10.3f %8i %10.3f %10.3f %5i \n",
	   TPeak->Cluster()->energyDep(), 
	   TPeak->CprIndex (),
	   TPeak->ClusterZ (),
	   TPeak->ClusterT0(),
	   TPeak->NHits    ());

    if (opt.Index("hits") >= 0) {
//-----------------------------------------------------------------------------
// print straw hits in the list
//-----------------------------------------------------------------------------
      printStrawHit(0,0,"banner",0,0);

      int  nhits, loc;
      nhits = TPeak->NHits();
      for (int i=0; i<nhits; i++) {
	loc   = TPeak->_index[i]._index;
	hit   = &TPeak->_shcol->at(loc);
	flags = *((int*) &TPeak->_shfcol->at(loc));

	printStrawHit(hit,step,"data",i,flags);
      }
    }
  }
}



//-----------------------------------------------------------------------------
void TAnaDump::printCalTimePeakCollection(const char* ModuleLabel, 
					  const char* ProductName, 
					  const char* ProcessName) {

  art::Handle<mu2e::CalTimePeakCollection>  handle;
  const mu2e::CalTimePeakCollection*        coll(0);
  const mu2e::CalTimePeak*                  tp(0);

  if (ProductName[0] != 0) {
    art::Selector  selector(art::ProductInstanceNameSelector(ProductName) &&
			    art::ProcessNameSelector        (ProcessName) && 
			    art::ModuleLabelSelector        (ModuleLabel)    );
    fEvent->get(selector, handle);
  }
  else {
    art::Selector  selector(art::ProcessNameSelector(ProcessName)         && 
			    art::ModuleLabelSelector(ModuleLabel)            );
    fEvent->get(selector, handle);
  }

  if (handle.isValid()) coll = handle.product();
  else {
    printf(">>> ERROR in TAnaDump::printSimParticleCollection: failed to locate collection");
    printf(". BAIL OUT. \n");
    return;
  }

  int banner_printed(0);

  for ( mu2e::CalTimePeakCollection::const_iterator j=coll->begin(); j != coll->end(); ++j ){
    tp = &(*j);

    if (banner_printed == 0) {
      printCalTimePeak(tp,"banner");
      banner_printed = 1;
    }
    printCalTimePeak(tp,"data+hits");
  }
}

//-----------------------------------------------------------------------------
void TAnaDump::printEventHeader() {

  printf(" Run / Subrun / Event : %10i / %10i / %10i\n",
	 fEvent->run(),
	 fEvent->subRun(),
	 fEvent->event());
}

//-----------------------------------------------------------------------------
void TAnaDump::printKalRep(const KalRep* Trk, const char* Opt, const char* Prefix) {

  TString opt = Opt;
  
  if ((opt == "") || (opt == "banner")) {
    printf("---------------------------------------------------------------------------------");
    printf("-----------------------------------------------------\n");
    printf("%s",Prefix);
    printf("  TrkID       Address    N  NA      P    sigma(p)    pT     costh    T0      T0Err   Omega");
    printf("      D0       Z0      Phi0   TanDip    Chi2    FCons\n");
    printf("---------------------------------------------------------------------------------");
    printf("-----------------------------------------------------\n");
  }
 
  if ((opt == "") || (opt.Index("data") >= 0)) {
    double chi2   = Trk->chisq();

    int    nhits(0);

    const TrkHotList* hots = Trk->hotList();
    for (TrkHotList::hot_iterator ihot=hots->begin(); ihot != hots->end(); ++ihot) {
      nhits++;
    }

    int    nact   = Trk->nActive();
    double t0     = Trk->t0().t0();
    double t0err  = Trk->t0().t0Err();
//-----------------------------------------------------------------------------
// in all cases define momentum at lowest Z - ideally, at the tracker entrance
//-----------------------------------------------------------------------------
    double s1     = Trk->firstHit()->kalHit()->hitOnTrack()->fltLen();
    double s2     = Trk->lastHit ()->kalHit()->hitOnTrack()->fltLen();
    double s      = std::min(s1,s2);

    double d0     = Trk->helix(s).d0();
    double z0     = Trk->helix(s).z0();
    double phi0   = Trk->helix(s).phi0();
    double omega  = Trk->helix(s).omega();
    double tandip = Trk->helix(s).tanDip();


    CLHEP::Hep3Vector fitmom = Trk->momentum(s);
    CLHEP::Hep3Vector momdir = fitmom.unit();
    BbrVectorErr      momerr = Trk->momentumErr(s);
    
    HepVector momvec(3);
    for (int i=0; i<3; i++) momvec[i] = momdir[i];
    
    double sigp = sqrt(momerr.covMatrix().similarity(momvec));
  
    double fit_consistency = Trk->chisqConsistency().consistency();
    int q         = Trk->charge();

    Hep3Vector trk_mom;

    trk_mom       = Trk->momentum(s);
    double mom    = trk_mom.mag();
    double pt     = trk_mom.perp();
    double costh  = trk_mom.cosTheta();

    printf("%5i %16p %3i %3i %8.3f %8.5f %7.3f %8.4f %7.3f %7.4f",
	   -1,
	   Trk,
	   nhits,
	   nact,
	   q*mom,sigp,pt,costh,t0,t0err
	   );

    printf(" %8.5f %8.3f %8.3f %8.4f %7.4f",
	   omega,d0,z0,phi0,tandip
	   );
    printf(" %8.3f %10.3e\n",
	   chi2,
	   fit_consistency);
  }

  if (opt.Index("hits") >= 0) {
  //-----------------------------------------------------------------------------
// print detailed information about the track hits
//-----------------------------------------------------------------------------
    const TrkHotList* hot_list = Trk->hotList();

    printf("--------------------------------------------------------------------");
    printf("---------------------------------------------------------------");
    printf("------------------------------------------------------\n");
    printf(" ih  SInd U A     len         x        y        z      HitT    HitDt");
    printf(" Ch Pl  L  W     T0       Xs      Ys        Zs     resid sigres");
    printf("    Rdrift   mcdoca totErr hitErr  t0Err penErr extErr\n");
    printf("--------------------------------------------------------------------");
    printf("---------------------------------------------------------------");
    printf("------------------------------------------------------\n");

    Hep3Vector pos;
    int i = 0;

    for (auto it=hot_list->begin(); it<hot_list->end(); it++) {
      // TrkStrawHit inherits from TrkHitOnTrk

      const mu2e::TrkStrawHit* hit = (const mu2e::TrkStrawHit*) &(*it);

      const mu2e::StrawHit* sh = &hit->strawHit();
      mu2e::Straw*   straw = (mu2e::Straw*) &hit->straw();

      hit->hitPosition(pos);

      double    len  = hit->fltLen();
      HepPoint  plen = Trk->position(len);
//-----------------------------------------------------------------------------
// find MC truth DOCA in a given straw
// start from finding the right vector of StepPointMC's
//-----------------------------------------------------------------------------
      int vol_id;
      int nstraws = fgListOfMCStrawHits->size();

      const mu2e::StepPointMC* step(0);

      for (int i=0; i<nstraws; i++) {
	mu2e::PtrStepPointMCVector  const& mcptr(fgListOfMCStrawHits->at(i));
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

      printf("%3i %5i %1i %1i %9.3f %8.3f %8.3f %9.3f %8.3f %7.3f",
	     ++i,
	     straw->index().asInt(), 
	     hit->isUsable(),
	     hit->isActive(),
	     len,
	     //	     hit->hitRms(),
	     plen.x(),plen.y(),plen.z(),
	     sh->time(), sh->dt()
	     );

      printf(" %2i %2i %2i %2i",
	     straw->id().getDevice(),
	     straw->id().getSector(),
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
      
      if      (hit->ambig()       == 0) printf(" * %6.3f",hit->driftRadius());
      else if (hit->ambig()*mcdoca > 0) printf("   %6.3f",hit->driftRadius()*hit->ambig());
      else                              printf(" ? %6.3f",hit->driftRadius()*hit->ambig());

      printf("  %7.3f",mcdoca);
      printf(" %6.3f %6.3f %6.3f %6.3f %6.3f\n",		 
	     hit->totalErr(),
	     hit->hitErr(),
	     hit->t0Err(),
	     hit->penaltyErr(),
	     hit->extErr()
	     );
    }
  }
}

//-----------------------------------------------------------------------------
void TAnaDump::printKalRepCollection(const char* ModuleLabel, 
				     const char* ProductName,
				     const char* ProcessName,
				     int hitOpt) {

  art::Handle<mu2e::KalRepPtrCollection> krepsHandle;

  if (ProductName[0] != 0) {
    art::Selector  selector(art::ProductInstanceNameSelector(ProductName) &&
			    art::ProcessNameSelector(ProcessName)         && 
			    art::ModuleLabelSelector(ModuleLabel)            );
    fEvent->get(selector,krepsHandle);
  }
  else {
    art::Selector  selector(art::ProcessNameSelector(ProcessName)         && 
			    art::ModuleLabelSelector(ModuleLabel)            );
    fEvent->get(selector,krepsHandle);
  }
//-----------------------------------------------------------------------------
// make sure collection exists
//-----------------------------------------------------------------------------
  if (! krepsHandle.isValid()) {
    printf("TAnaDump::printKalRepCollection: no KalRepPtrCollection for module %s, BAIL OUT\n",
	   ModuleLabel);
    return;
  }

  art::Handle<mu2e::PtrStepPointMCVectorCollection> mcptrHandle;
  if (ModuleLabel[0] != 0) {
    fEvent->getByLabel("makeSH","StrawHitMCPtr",mcptrHandle);
    if (mcptrHandle.isValid()) {
      fgListOfMCStrawHits = (mu2e::PtrStepPointMCVectorCollection*) mcptrHandle.product();
    }
    else {
      printf(">>> ERROR in TAnaDump::printKalRepCollection: failed to locate StepPointMCCollection makeSH:StrawHitMCPtr\n");
      fgListOfMCStrawHits = NULL;
    }
  }

  int ntrk = krepsHandle->size();

  const KalRep *trk;

  int banner_printed = 0;
  for (int i=0; i<ntrk; i++) {
    trk = /*(KalRep*)*/ krepsHandle->at(i).get();
    if (banner_printed == 0) {
      printKalRep(trk,"banner");
      banner_printed = 1;
    }
    printKalRep(trk,"data");
    if(hitOpt>0) printKalRep(trk,"hits");
  }
 
}


//-----------------------------------------------------------------------------
void TAnaDump::printGenParticle(const mu2e::GenParticle* P, const char* Opt) {

  TString opt = Opt;
  
  if ((opt == "") || (opt == "banner")) {
    printf("------------------------------------------------------------------------------------\n");
    printf("Index                 generator     PDG      Time      Momentum       Pt       CosTh\n");
    printf("------------------------------------------------------------------------------------\n");
  }
  
  if ((opt == "") || (opt == "data")) {
    int    gen_code   = P->generatorId().id();
    std::string gen_name   = P->generatorId().name();
    int    pdg_code   = P->pdgId();
    double time       = P->time();
    
    double mom   = P->momentum().vect().mag();
    double pt    = P->momentum().vect().perp();
    double costh = P->momentum().vect().cosTheta();
    
    printf("%5i %2i:%-26s %3i %10.3f %10.3f %10.3f %10.3f\n",
	   -1,gen_code,gen_name.data(),pdg_code,time,mom,pt,costh);
  }
}

//-----------------------------------------------------------------------------
// there could be multiple collections in the event
//-----------------------------------------------------------------------------
void TAnaDump::printGenParticleCollections() {
  
  std::vector<art::Handle<mu2e::GenParticleCollection>> list_of_gp;

  const mu2e::GenParticleCollection*        coll(0);
  const mu2e::GenParticle*        genp(0);

  const art::Provenance* prov;

  //  art::Selector  selector(art::ProductInstanceNameSelector("mu2e::GenParticleCollection"));
  art::Selector  selector(art::ProductInstanceNameSelector(""));

  fEvent->getMany(selector, list_of_gp);

  const art::Handle<mu2e::GenParticleCollection>* handle;

  int banner_printed;
  for (std::vector<art::Handle<mu2e::GenParticleCollection>> ::const_iterator it = list_of_gp.begin();
       it != list_of_gp.end(); it++) {
    handle = it.operator -> ();

    if (handle->isValid()) {
      coll = handle->product();
      prov = handle->provenance();

      printf("moduleLabel = %-20s, producedClassname = %-30s, productInstanceName = %-20s\n",
	     prov->moduleLabel().data(),
	     prov->producedClassName().data(),
	     prov->productInstanceName().data());

      banner_printed = 0;
      for (std::vector<mu2e::GenParticle>::const_iterator ip = coll->begin();
	   ip != coll->end(); ip++) {
	genp = ip.operator -> ();
	if (banner_printed == 0) {
	  printGenParticle(genp,"banner");
	  banner_printed = 1;
	}
	printGenParticle(genp,"data");
      }

      
    }
    else {
      printf(">>> ERROR in TAnaDump::printStepPointMCCollection: failed to locate collection");
      printf(". BAIL OUT. \n");
      return;
    }
  }
}


// //-----------------------------------------------------------------------------
//   void TAnaDump::printCaloHit(const CaloHit* Hit, const char* Opt) {
//     //    int row, col;
//     TString opt = Opt;

//     if ((opt == "") || (opt == "banner")) {
//       printf("--------------------------------------\n");
//       printf("RID      Time   Energy                \n");
//       printf("--------------------------------------\n");
//     }
    
//     if ((opt == "") || (opt == "data")) {
//       printf("%7i  %10.3f %10.3f \n",
// 	     Hit->id(),
// 	     Hit->time(),
// 	     Hit->energyDep()); 
//     }
//   }


//-----------------------------------------------------------------------------
void TAnaDump::printCaloHits(const char* ModuleLabel, 
			     const char* ProductName, 
			     const char* ProcessName) {

  printf(">>>> ModuleLabel = %s\n",ModuleLabel);

  //data about hits in the calorimeter crystals

  art::Selector  selector(art::ProductInstanceNameSelector(ProductName) &&
			  art::ProcessNameSelector(ProcessName)         && 
			  art::ModuleLabelSelector(ModuleLabel)            );

  art::Handle<mu2e::CaloHitCollection> caloHitsHandle;

  fEvent->get(selector,caloHitsHandle);

  const mu2e::CaloHitCollection* caloHits;

  caloHits = caloHitsHandle.operator ->();

  int nhits = caloHits->size();

  const mu2e::CaloHit* hit;

  printf("--------------------------------------\n");
  printf("RID      Time   Energy                \n");
  printf("--------------------------------------\n");

  for (int ic=0; ic<nhits; ic++) {
    hit  = &caloHits->at(ic);
    printf("%7i  %10.3f %10.3f \n",
	   hit->id(),
	   hit->time(),
	   hit->energyDep()); 
  }
}


//-----------------------------------------------------------------------------
void TAnaDump::printDiskCalorimeter() {
  const mu2e::DiskCalorimeter* cal;
  const mu2e::Disk* disk;
  
  art::ServiceHandle<mu2e::GeometryService> geom;
    
  if (geom->hasElement<mu2e::DiskCalorimeter>() ) {
    mu2e::GeomHandle<mu2e::DiskCalorimeter> dc;
    cal = dc.operator->();
  }
  else {
    printf(">>> ERROR: disk calorimeter not found.\n");
    return;
  }

  int nd = cal->nDisk();
  printf(" ndisks = %i\n", nd);
  printf(" crystal size  : %10.3f\n", 2*cal->caloGeomInfo().crystalHalfTrans());
  printf(" crystal length: %10.3f\n", 2*cal->caloGeomInfo().crystalHalfLength());

  for (int i=0; i<nd; i++) {
    disk = &cal->disk(i);
    printf(" ---- disk # %i\n",i);
    printf(" Rin  : %10.3f  Rout : %10.3f\n", disk->innerRadius(),disk->outerRadius());
    printf(" X : %12.3f Y : %12.3f Z : %12.3f\n",
	   disk->origin().x(),
	   disk->origin().y(),
	   disk->origin().z());
    printf(" Xsize : %10.3f Ysize : %10.3f Zsize : %10.3f\n", 
	   disk->size().x(),
	   disk->size().y(),
	   disk->size().z()
	   );
  }
}


//-----------------------------------------------------------------------------
void TAnaDump::printCaloCrystalHits(const char* ModuleLabel, 
				    const char* ProductName,
				    const char* ProcessName) {

  art::Selector  selector(art::ProductInstanceNameSelector(ProductName) &&
			  art::ProcessNameSelector(ProcessName)         && 
			  art::ModuleLabelSelector(ModuleLabel)            );

  art::Handle<mu2e::CaloCrystalHitCollection> caloCrystalHitsHandle;

  fEvent->get(selector,caloCrystalHitsHandle);

  const mu2e::CaloCrystalHitCollection* caloCrystalHits;

  caloCrystalHits = caloCrystalHitsHandle.operator->();

  int nhits = caloCrystalHits->size();

  const mu2e::CaloCrystalHit* hit;

  printf("----------------------------------------------------------------\n");
  printf("CrystalID      Time   Energy    EnergyTot  NRoids               \n");
  printf("----------------------------------------------------------------\n");

  for (int ic=0; ic<nhits; ic++) {
    hit  = &caloCrystalHits->at(ic);

    printf("%7i  %10.3f %10.3f %10.3f %5i\n",
	   hit->id(),
	   hit->time(),
	   hit->energyDep(),
	   hit->energyDepTotal(),
	   hit->numberOfROIdsUsed());
  }
}

//------------------------------------------------------------------
// void TAnaDump::printTrackClusterLink(const char* ModuleLabel, 
// 			     const char* ProductName,
// 			     const char* ProcessName) {

//   printf(">>>> ModuleLabel = %s\n",ModuleLabel);

//   //data about hits in the calorimeter crystals

//   art::Handle<mu2e::TrackClusterLink> trkCluLinkHandle;
//   const mu2e::TrackClusterLink* trkCluLink;

//   if (ProductName[0] != 0) {
//     art::Selector  selector(art::ProductInstanceNameSelector(ProductName) &&
// 			    art::ProcessNameSelector(ProcessName)         && 
// 			    art::ModuleLabelSelector(ModuleLabel)            );
//     fEvent->get(selector, trkCluLinkHandle);
//   }
//   else {
//     art::Selector  selector(art::ProcessNameSelector(ProcessName)         && 
// 			    art::ModuleLabelSelector(ModuleLabel)            );
//     fEvent->get(selector, trkCluLinkHandle);
//   }
  
//   trkCluLink = trkCluLinkHandle.operator ->();

//   int nhits = trkCluLink->size();

//   const mu2e::CaloCluster* clu;
//   //const KalRep *trk;

//   int banner_printed = 0;
//   const mu2e::TrkToCaloExtrapol* extrk;
  


//   for (int i=0; i<nhits; i++) {
//     clu = &(*(trkCluLink->at(i).second.get()) );
//     extrk = &(*(trkCluLink->at(i).first));
//     if (banner_printed == 0) {
//       printCaloCluster(clu, "banner");
//       //banner_printed = 1;
//     }
//     printCaloCluster(clu,"data");
    
//     if (banner_printed == 0) {
//       printTrkToCaloExtrapol(extrk,"banner");
//       //printKalRep(trk,"banner");
//       banner_printed = 1;
//     }
//     // printKalRep(trk,"data");
//     printTrkToCaloExtrapol(extrk,"data");
//   }
  
// }

////////////////////////////////////////////////////////////////////////////////

void TAnaDump::printTrkToCaloExtrapol(const mu2e::TrkToCaloExtrapol* trkToCalo,
				      const char* Opt) {
 TString opt = Opt;

  if ((opt == "") || (opt == "banner")) {
    printf("-------------------------------------------------------------------------------------------------------\n");
    printf("sectionId      Time     ExtPath     Ds       FitCon      t0          X           Y        Z          Mom  \n");
    printf("-------------------------------------------------------------------------------------------------------\n");
  }
  
  if ((opt == "") || (opt.Index("data") >= 0)) {

    double ds = trkToCalo->pathLengthExit()-trkToCalo->pathLengthEntrance();
  
    printf("%6i %10.3f %10.3f %8.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f \n",
	   trkToCalo->sectionId(),
	   trkToCalo->time(),
	   trkToCalo->pathLengthEntrance(),
	   ds,
	   trkToCalo->fitConsistency(),
	   trkToCalo->t0(),
	   trkToCalo->entrancePosition().x(),
	   trkToCalo->entrancePosition().y(),
	   trkToCalo->entrancePosition().z(),
	   trkToCalo->momentum().mag() );
  }
  
}

////////////////////////////////////////////////////////////////////////////////

void TAnaDump::printTrkToCaloExtrapolCollection(const char* ModuleLabel, 
						const char* ProductName,
						const char* ProcessName) {

  printf(">>>> ModuleLabel = %s\n",ModuleLabel);

  //data about hits in the calorimeter crystals

  art::Handle<mu2e::TrkToCaloExtrapolCollection> trkToCaloHandle;
  const mu2e::TrkToCaloExtrapolCollection* trkToCalo;

  if (ProductName[0] != 0) {
    art::Selector  selector(art::ProductInstanceNameSelector(ProductName) &&
			    art::ProcessNameSelector(ProcessName)         && 
			    art::ModuleLabelSelector(ModuleLabel)            );
    fEvent->get(selector, trkToCaloHandle);
  }
  else {
    art::Selector  selector(art::ProcessNameSelector(ProcessName)         && 
			    art::ModuleLabelSelector(ModuleLabel)            );
    fEvent->get(selector, trkToCaloHandle);
  }

  trkToCalo = trkToCaloHandle.operator ->();

  int nhits = trkToCalo->size();

  const mu2e::TrkToCaloExtrapol* hit;
  
  int banner_printed = 0;
  for (int i=0; i<nhits; i++) {
    hit = &trkToCalo->at(i);
    if (banner_printed == 0) {
      printTrkToCaloExtrapol(hit, "banner");
      banner_printed = 1;
    }
    printTrkToCaloExtrapol(hit,"data");
  }
  
}

////////////////////////////////////////////////////////////////////////////////


//-----------------------------------------------------------------------------
void TAnaDump::printStrawHit(const mu2e::StrawHit* Hit, const mu2e::StepPointMC* Step, const char* Opt, int IHit, int Flags) {
    TString opt = Opt;
    opt.ToLower();

    if ((opt == "") || (opt.Index("banner") >= 0)) {
      printf("-----------------------------------------------------------------------------------");
      printf("-------------------------------------------------\n");
      printf("   I   SHID  Flags      Station Sector Layer Straw     Time          dt       eDep ");
      printf("     PDG PDG(M)       GENID       ID         p   \n");
      printf("-----------------------------------------------------------------------------------");
      printf("-------------------------------------------------\n");
    }

    if (opt == "banner") return;

    mu2e::GeomHandle<mu2e::TTracker> ttHandle;
    const mu2e::TTracker* tracker = ttHandle.get();

    const mu2e::Straw* straw;

    straw = &tracker->getStraw(Hit->strawIndex());

// 12 - 11 -2013 giani added some MC info of the straws
    //Hep3Vector mom = Step->momentum();
    // double pt = std::sqrt(mom.mag()*mom.mag() - mom.z()*mom.z());

    const mu2e::SimParticle * sim (0);
    
    int      pdg_id(-1), mother_pdg_id(-1), gen_index(-1), sim_id(-1);
    double   mc_mom(-1.);

    mu2e::GenId gen_id;

    if (Step) {
      art::Ptr<mu2e::SimParticle> const& simptr = Step->simParticle(); 
      art::Ptr<mu2e::SimParticle> mother = simptr;

      while(mother->hasParent()) mother = mother->parent();

      sim = mother.operator ->();

      pdg_id        = simptr->pdgId();
      mother_pdg_id = sim->pdgId();

      if (simptr->fromGenerator()) gen_index = simptr->genParticle()->generatorId().id();
      else                         gen_index = -1;

      sim_id        = simptr->id().asInt();
      mc_mom        = Step->momentum().mag();
    }
    
    if ((opt == "") || (opt == "data")) {
      if (IHit  >= 0) printf("%5i " ,IHit);
      else            printf("    ");

      printf("%5i",Hit->strawIndex().asInt());

      if (Flags >= 0) printf(" %08x",Flags);
      else            printf("        ");
      printf("  %5i  %5i   %5i   %5i   %8.3f   %8.3f   %9.6f   %4i   %4i  %10i  %10i %8.3f\n",
	     straw->id().getDevice(),
	     straw->id().getSector(),
	     straw->id().getLayer(),
	     straw->id().getStraw(),
	     Hit->time(),
	     Hit->dt(),
	     Hit->energyDep(),
	     pdg_id,
	     mother_pdg_id,
	     gen_index,
	     sim_id,
	     mc_mom);
    }
  }


//-----------------------------------------------------------------------------
void TAnaDump::printStrawHitCollection(const char* ModuleLabel, 
				       const char* ProductName,
				       const char* ProcessName,
				       double TMin, double TMax) {

  const char* oname = "TAnaDump::printStrawHitCollection";

  art::Handle<mu2e::StrawHitCollection> shcHandle;
  const mu2e::StrawHitCollection*       shc;

  art::Handle<mu2e::StrawHitFlagCollection> shflagH;
  const mu2e::StrawHitFlagCollection*       shfcol;
  

//-----------------------------------------------------------------------------
// get straw hits
//-----------------------------------------------------------------------------
  if (ProductName[0] != 0) {
    art::Selector  selector(art::ProductInstanceNameSelector(ProductName) &&
			    art::ProcessNameSelector(ProcessName)         && 
			    art::ModuleLabelSelector(ModuleLabel)            );
    fEvent->get(selector, shcHandle);
  }
  else {
    art::Selector  selector(art::ProcessNameSelector(ProcessName)         && 
			    art::ModuleLabelSelector(ModuleLabel)            );
    fEvent->get(selector, shcHandle);
  }

  if (shcHandle.isValid()) shc = shcHandle.product();
  else {
    printf(">>> ERROR in %s: Straw Hit Collection by \"%s\" doesn't exist. Bail Out.\n",
	   oname,ModuleLabel);
    return;
  }
//-----------------------------------------------------------------------------
// get straw hit flags (half-hack)
//-----------------------------------------------------------------------------
  fEvent->getByLabel(fFlagBgrHitsModuleLabel.Data(),shflagH);

  if (shflagH.isValid()) shfcol = shflagH.product();
  else {
    printf(">>> ERROR in %s: Straw Hit Flag Collection by \"%s\" doesn't exist. Bail Out.\n",
	   oname,fFlagBgrHitsModuleLabel.Data());
    return;
  }

// 12 - 11 -2013 giani added some MC info of the straws
  art::Handle<mu2e::PtrStepPointMCVectorCollection> mcptrHandleStraw;
  fEvent->getByLabel(ModuleLabel,"StrawHitMCPtr",mcptrHandleStraw);
  mu2e::PtrStepPointMCVectorCollection const* hits_mcptrStraw = mcptrHandleStraw.product();
 
  //------------------------------------------------------------

  int nhits = shc->size();

  const mu2e::StrawHit* hit;

  int flags;

  int banner_printed = 0;
  for (int i=0; i<nhits; i++) {
    mu2e::PtrStepPointMCVector const& mcptr(hits_mcptrStraw->at(i ) );
    const mu2e::StepPointMC* Step = mcptr[0].operator ->();
    
    hit   = &shc->at(i);
					// assuming it doesn't move beyond 32 bits
    flags = *((int*) &shfcol->at(i));
    if (banner_printed == 0) {
      printStrawHit(hit, Step, "banner");
      banner_printed = 1;
    }
    if ((hit->time() >= TMin) && (hit->time() <= TMax)) {
      printStrawHit(hit, Step, "data", i, flags);
    }
  }
 
}



//-----------------------------------------------------------------------------
void TAnaDump::printSimParticle(const mu2e::SimParticle* P, const char* Opt) {

    TString opt = Opt;

    if ((opt == "") || (opt == "banner")) {
      printf("----------------------------------------------------------------------------------\n");
      printf("Index  Primary     ID Parent   PDG      X      Y       Z      T      Px      Py     Pz      E \n");
      printf("----------------------------------------------------------------------------------\n");
    }
 
    if ((opt == "") || (opt == "data")) {
      int  id        = P->id().asInt();
      
      int  parent_id (-1);

      if (P->parent()) {
	parent_id = P->parent()->id().asInt();
      }
      int  pdg_id    = P->pdgId();
      int  primary   = P->isPrimary();

      printf("%5i %10i %8i %5i %10i  %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f \n",
	     -1, primary, id, parent_id, pdg_id, 
	     P->startPosition().x(),
	     P->startPosition().y(),
	     P->startPosition().z(),
	     P->startGlobalTime(),
	     P->startMomentum().x(),
	     P->startMomentum().y(),
	     P->startMomentum().z(),
	     P->startMomentum().e());
    }
  }



//-----------------------------------------------------------------------------
void TAnaDump::printSimParticleCollection(const char* ModuleLabel, 
					  const char* ProductName, 
					  const char* ProcessName) {

  art::Handle<mu2e::SimParticleCollection> handle;
  const mu2e::SimParticleCollection*       coll(0);
  const mu2e::SimParticle*                 simp(0);

  if (ProductName[0] != 0) {
    art::Selector  selector(art::ProductInstanceNameSelector(ProductName) &&
			    art::ProcessNameSelector        (ProcessName) && 
			    art::ModuleLabelSelector        (ModuleLabel)    );
    fEvent->get(selector, handle);
  }
  else {
    art::Selector  selector(art::ProcessNameSelector(ProcessName)         && 
			    art::ModuleLabelSelector(ModuleLabel)            );
    fEvent->get(selector, handle);
  }

  if (handle.isValid()) coll = handle.product();
  else {
    printf(">>> ERROR in TAnaDump::printSimParticleCollection: failed to locate collection");
    printf(". BAIL OUT. \n");
    return;
  }

  int banner_printed(0);

  for ( mu2e::SimParticleCollection::const_iterator j=coll->begin(); j != coll->end(); ++j ){
    simp = &j->second;

    if (banner_printed == 0) {
      printSimParticle(simp,"banner");
      banner_printed = 1;
    }
    printSimParticle(simp,"data");
  }
}


//-----------------------------------------------------------------------------
void TAnaDump::printStepPointMC(const mu2e::StepPointMC* Step, const char* Opt) {
  const char* oname = "TAnaDump::printStepPointMC";
    TString opt = Opt;

    if ((opt == "") || (opt.Index("banner") >= 0)) {
      printf("------------------------------------------------------------------------------------------");
      printf("----------------------------");
      printf("-------------------------------------------------------------------------------------------------------------\n");
      printf("  Vol          PDG    ID GenIndex PPdg ParID     X          Y         Z            T      ");
      printf("  X0          Y0         Z0 ");
      printf("  Edep(Tot) Edep(NI) Edep(I)  Step   EndCode  Energy      EKin     Mom      Pt    doca  Creation      StopProc\n");
      printf("------------------------------------------------------------------------------------------");
      printf("----------------------------");
      printf("-------------------------------------------------------------------------------------------------------------\n");
    }

    mu2e::GeomHandle<mu2e::TTracker> ttHandle;
    const mu2e::TTracker* tracker = ttHandle.get();
  
    art::Ptr<mu2e::SimParticle> const& simptr = Step->simParticle();
    const mu2e::SimParticle* sim  = simptr.operator ->();
    if (sim == NULL) {
      printf(">>> ERROR: %s sim == NULL\n",oname);
    }

    art::Ptr<mu2e::SimParticle> const& parentptr = sim->parent();

    int parent_pdg_id (-1), parent_id(-1);

    const mu2e::SimParticle* par = parentptr.get();
    if (par != NULL) {
      parent_pdg_id = (int) par->pdgId();
      parent_id     = (int) par->id().asInt();
    }

    const mu2e::Straw* straw = &tracker->getStraw(mu2e::StrawIndex(Step->volumeId()));

    const Hep3Vector* v1 = &straw->getMidPoint();
    HepPoint p1(v1->x(),v1->y(),v1->z());

    const Hep3Vector* v2 = &Step->position();
    HepPoint    p2(v2->x(),v2->y(),v2->z());

    TrkLineTraj trstraw(p1,straw->getDirection()  ,0.,0.);
    TrkLineTraj trstep (p2,Step->momentum().unit(),0.,0.);

    // 2015-02-16 G. Pezzu and Murat change in the print out to be finished
    // 2015-02-25 P.Murat: fix sign - trajectory is the first !
    //  however, the sign of the disptance of closest approach is invariant
    // wrt the order
    TrkPoca poca(trstep, 0., trstraw, 0.);
    
    double doca = poca.doca();
    
    //    art::Ptr<mu2e::GenParticle> const& apgen = sim->genParticle();
    //    mu2e::GenParticle* gen = (mu2e::GenParticle*) &(*sim->genParticle());

    double mass = sim->startMomentum().m();

    double pabs = Step->momentum().mag();
    double energy = std::sqrt(pabs*pabs+mass*mass);
    double ekin  = energy-mass;
        
    Hep3Vector mom = Step->momentum();
    double pt = std::sqrt(pabs*pabs - mom.z()*mom.z());

    art::Handle<mu2e::PhysicalVolumeInfoCollection> volumes;
    fEvent->getRun().getByLabel("g4run", volumes);

//2014-26-11 gianipez added the timeoffsets to the steppoints time
    fgTimeOffsets->updateMap(*fEvent);
    
    double stepTime = fgTimeOffsets->timeWithOffsetsApplied(*Step);
    //    const mu2e::PhysicalVolumeInfo& pvinfo = volumes->at(sim->startVolumeIndex());
    //    const mu2e::PhysicalVolumeInfo& pvinfo = volumes->at(Step->volumeId()); - sometimes crashes..

    if ((opt == "") || (opt.Index("data") >= 0)) {
      printf("%5i %12i %5i %5i %5i %5i %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %8.2f %8.2f %8.2f %8.3f %4i %10.3f %8.3f %8.3f %8.3f %6.2f %-12s %-s\n",
	     (int) Step->volumeId(),
	     //	     pvinfo.name().data(), // smth is wrong with the name defined by volumeId()....
	     (int) sim->pdgId(),
	     (int) sim->id().asInt(),
	     (int) sim->generatorIndex(),
	     parent_pdg_id,
	     parent_id,
	     Step->position().x(),
	     Step->position().y(),
	     Step->position().z(),
	     stepTime,                             // Step->time(),
	     sim->startPosition().x(),
	     sim->startPosition().y(),
	     sim->startPosition().z(),
	     Step->totalEDep(),
	     Step->nonIonizingEDep(),
	     Step->ionizingEdep(),
	     Step->stepLength(),
	     Step->endProcessCode().id(),
	     energy,
	     ekin,
	     pabs,
	     pt,
	     doca,
	     sim->creationCode().name().data(),
	     Step->endProcessCode().name().data());
    }
}


//-----------------------------------------------------------------------------
void TAnaDump::printStepPointMCCollection(const char* ModuleLabel, 
					  const char* ProductName,
					  const char* ProcessName) {

  art::Handle<mu2e::StepPointMCCollection> handle;
  const mu2e::StepPointMCCollection*       coll(0);

  if (ProductName[0] != 0) {
    art::Selector  selector(art::ProductInstanceNameSelector(ProductName) &&
			    art::ProcessNameSelector        (ProcessName) && 
			    art::ModuleLabelSelector        (ModuleLabel)    );
    fEvent->get(selector, handle);
  }
  else {
    art::Selector  selector(art::ProcessNameSelector(ProcessName)         && 
			    art::ModuleLabelSelector(ModuleLabel)            );
    fEvent->get(selector, handle);
  }

  if (handle.isValid()) coll = handle.product();
  else {
    printf(">>> ERROR in TAnaDump::printStepPointMCCollection: failed to locate collection");
    printf(". BAIL OUT. \n");
    return;
  }

  int nsteps = coll->size();

  const mu2e::StepPointMC* step;


  int banner_printed = 0;
  for (int i=0; i<nsteps; i++) {
    step = &coll->at(i);
    if (banner_printed == 0) {
      printStepPointMC(step, "banner");
      banner_printed = 1;
    }
    printStepPointMC(step,"data");
  }
 
}

//-----------------------------------------------------------------------------
// this is 
//-----------------------------------------------------------------------------
void TAnaDump::printStepPointMCVectorCollection(const char* ModuleLabel, 
						const char* ProductName,
						const char* ProcessName) {
  int banner_printed(0);

  art::Handle<mu2e::PtrStepPointMCVectorCollection> mcptrHandle;
  if (ModuleLabel[0] != 0) {
    fEvent->getByLabel(ModuleLabel,"StrawHitMCPtr", mcptrHandle);
    if (mcptrHandle.isValid()) {
      fgListOfMCStrawHits = (mu2e::PtrStepPointMCVectorCollection*) mcptrHandle.product();
    }
    else {
      printf(">>> ERROR in TAnaDump::printStepPointMCVectorCollection: failed to locate collection");
      printf(". BAIL OUT. \n");
      return;
    }
  }
  else {
    printf(">>> ERROR in TAnaDump::printStepPointMCVectorCollection: ModuleLabel undefined");
    printf(". BAIL OUT. \n");
    return;
  }
  
  const mu2e::StepPointMC  *step;

  int nv = fgListOfMCStrawHits->size();
  printf("[TAnaDump::printStepPointMCVectorCollection] nv = %i\n",nv);

  for (int i=0; i<nv; i++) {
    mu2e::PtrStepPointMCVector  const& mcptr(fgListOfMCStrawHits->at(i));
    int nj = mcptr.size();
    printf("[TAnaDump::printStepPointMCVectorCollection] i = %i, nj=%i\n",i,nj);
    for (int j=0; j<nj; j++) {
      step = &(*mcptr.at(j));
      if (banner_printed == 0) {
	printStepPointMC(step, "banner");
	banner_printed = 1;
      }
      printStepPointMC(step,"data");
    }
  }
}
 
//-----------------------------------------------------------------------------
void TAnaDump::printStrawHitMCTruth(const mu2e::StrawHitMCTruth* Hit, const char* Opt) {
  TString opt = Opt;
  
  if ((opt == "") || (opt == "banner")) {
    printf("--------------------------------------------------------------------\n");
    printf(" Time Distance DistToMid         dt       eDep \n");
    printf("--------------------------------------------------------------------\n");
  }

  if ((opt == "") || (opt == "data")) {
    printf("%12.5f  %12.5f  %12.5f\n",
	   Hit->driftTime(),
	   Hit->driftDistance(),
	   Hit->distanceToMid());
  }
}


//-----------------------------------------------------------------------------
void TAnaDump::printStrawHitMCTruthCollection(const char* ModuleLabel, 
					      const char* ProductName,
					      const char* ProcessName) {

  art::Handle<mu2e::StrawHitMCTruthCollection> shcHandle;
  const mu2e::StrawHitMCTruthCollection*       shc;

  if (ProductName[0] != 0) {
    art::Selector  selector(art::ProductInstanceNameSelector(ProductName) &&
			    art::ProcessNameSelector(ProcessName)         && 
			    art::ModuleLabelSelector(ModuleLabel)            );
    fEvent->get(selector, shcHandle);
  }
  else {
    art::Selector  selector(art::ProcessNameSelector(ProcessName)         && 
			    art::ModuleLabelSelector(ModuleLabel)            );
    fEvent->get(selector, shcHandle);
  }

  shc = shcHandle.product();

  int nhits = shc->size();

  const mu2e::StrawHitMCTruth* hit;


  int banner_printed = 0;
  for (int i=0; i<nhits; i++) {
    hit = &shc->at(i);
    if (banner_printed == 0) {
      printStrawHitMCTruth(hit, "banner");
      banner_printed = 1;
    }
    printStrawHitMCTruth(hit,"data");
  }
 
}

//-----------------------------------------------------------------------------
void TAnaDump::printStrawHitPosition(const mu2e::StrawHitPosition* Hit, const char* Opt) {
  TString opt = Opt;
  
  if ((opt == "") || (opt == "banner")) {
    printf("---------------------------------------------------------------------------\n");
    printf("  STIndex     X        Y        Z        Wdist     Pres     RRes      Flag \n");
    printf("---------------------------------------------------------------------------\n");
  }

  int flag = *((int*) &Hit->flag());

  double wres = Hit->posRes(mu2e::StrawHitPosition::phi);
  if (wres > 999.) wres = 999.;

  double rres = Hit->posRes(mu2e::StrawHitPosition::rho);
  if (rres > 999.) rres = 999.;

  if ((opt == "") || (opt == "data")) {
    printf("   %6i %8.3f %8.3f %9.3f %8.3f  %7.2f  %7.2f  0x%08x\n",
	   Hit->stereoHitIndex(),
	   Hit->pos().x(),
	   Hit->pos().y(),
	   Hit->pos().z(),
	   Hit->wireDist(),
	   wres,
	   rres,
	   flag);
  }
}


//-----------------------------------------------------------------------------
void TAnaDump::printStrawHitPositionCollection(const char* ModuleLabel, 
					       const char* ProductName,
					       const char* ProcessName) {

  art::Handle<mu2e::StrawHitPositionCollection> spcHandle;
  const mu2e::StrawHitPositionCollection*       spc;

  if (ProductName[0] != 0) {
    art::Selector  selector(art::ProductInstanceNameSelector(ProductName) &&
			    art::ProcessNameSelector(ProcessName)         && 
			    art::ModuleLabelSelector(ModuleLabel)            );
    fEvent->get(selector, spcHandle);
  }
  else {
    art::Selector  selector(art::ProcessNameSelector(ProcessName)         && 
			    art::ModuleLabelSelector(ModuleLabel)            );
    fEvent->get(selector, spcHandle);
  }

  spc = spcHandle.product();

  int nhits = spc->size();

  const mu2e::StrawHitPosition* pos;


  int banner_printed = 0;
  for (int i=0; i<nhits; i++) {
    pos = &spc->at(i);
    if (banner_printed == 0) {
      printStrawHitPosition(pos, "banner");
      banner_printed = 1;
    }
    printStrawHitPosition(pos,"data");
  }
 
}


//-----------------------------------------------------------------------------
void TAnaDump::printTrackClusterMatch(const mu2e::TrackClusterMatch* Tcm, const char* Opt) {

  TString opt = Opt;
  
  if ((opt == "") || (opt == "banner")) {
    printf("------------------------------------------------------------------------------------------\n");
    printf("  Disk   Track    Cluster   chi2                                                         \n");
    printf("------------------------------------------------------------------------------------------\n");
  }

  if ((opt == "") || (opt == "data")) {

    const mu2e::CaloCluster*       cl  = Tcm->caloCluster();
    const mu2e::TrkToCaloExtrapol* tex = Tcm->textrapol  ();

    int disk     = cl->sectionId();
    double chi2  = Tcm->chi2();

    printf("%5i %16p  %16p  %10.3f\n",
	   disk,cl,tex,chi2);
  }
}



//-----------------------------------------------------------------------------
void TAnaDump::printTrackClusterMatchCollection(const char* ModuleLabel, 
						const char* ProductName,
						const char* ProcessName) {

  printf(">>>> ModuleLabel = %s\n",ModuleLabel);

  art::Handle<mu2e::TrackClusterMatchCollection> handle;
  const mu2e::TrackClusterMatchCollection*       coll;

  art::Selector  selector(art::ProductInstanceNameSelector(ProductName) &&
			  art::ProcessNameSelector(ProcessName)         && 
			  art::ModuleLabelSelector(ModuleLabel)            );

  fEvent->get(selector,handle);

  if (handle.isValid()) coll = handle.product();
  else {
    printf(">>> ERROR in TAnaDump::printTrackClusterMatchCollection: failed to locate collection");
    printf(". BAIL OUT. \n");
    return;
  }

  int nm = coll->size();

  const mu2e::TrackClusterMatch* obj;

  int banner_printed = 0;

  for (int i=0; i<nm; i++) {
    obj = &coll->at(i);
    if (banner_printed == 0) {
      printTrackClusterMatch(obj, "banner");
      banner_printed = 1;
    }
    printTrackClusterMatch(obj,"data");
  }
 
}


//-----------------------------------------------------------------------------
void TAnaDump::refitTrack(void* Trk, double NSig) {
  KalRep* trk = (KalRep*) Trk;

  const TrkHotList* hot_list = trk->hotList();

  for(TrkHotList::hot_iterator it=hot_list->begin(); it<hot_list->end(); it++) {

    // TrkStrawHit inherits from TrkHitOnTrk

    TrkHitOnTrk* hot = (TrkHitOnTrk*) &(*it);

    const mu2e::TrkStrawHit* hit = (const mu2e::TrkStrawHit*) &(*it);

    double res = hit->resid();

    if (fabs(res) > 0.1*NSig) {
      trk->deactivateHot(hot);
    }
  }

  trk->fit();

  printKalRep(trk, "hits");
}
