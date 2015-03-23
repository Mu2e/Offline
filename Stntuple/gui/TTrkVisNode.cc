///////////////////////////////////////////////////////////////////////////////
// May 04 2013 P.Murat
// 
// in 'XY' mode draw calorimeter clusters as circles with different colors 
// in 'Cal' mode draw every detail...
///////////////////////////////////////////////////////////////////////////////
#include "TVirtualX.h"
#include "TPad.h"
#include "TStyle.h"
#include "TVector3.h"
#include "TLine.h"
#include "TArc.h"
#include "TArrow.h"
#include "TBox.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"

#include "Stntuple/gui/TEvdTrack.hh"
#include "Stntuple/gui/TTrkVisNode.hh"
#include "Stntuple/gui/TEvdStraw.hh"
#include "Stntuple/gui/TEvdStrawHit.hh"
#include "Stntuple/gui/TEvdTrkStrawHit.hh"
#include "Stntuple/gui/TEvdStation.hh"
#include "Stntuple/gui/TEvdPanel.hh"
#include "Stntuple/gui/TEvdStrawTracker.hh"
#include "Stntuple/gui/TStnVisManager.hh"

#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"

#include "Stntuple/mod/TAnaDump.hh"

ClassImp(TTrkVisNode)

//_____________________________________________________________________________
TTrkVisNode::TTrkVisNode() : TVisNode("") {
}

//_____________________________________________________________________________
TTrkVisNode::TTrkVisNode(const char* name, const mu2e::TTracker* Tracker, TStnTrackBlock* TrackBlock): 
  TVisNode(name) {
  fTracker    = new TEvdStrawTracker(Tracker);
  fTrackBlock = TrackBlock;

  fArc        = new TArc;
  fEventTime  = 0;
  fTimeWindow = 1.e6;

  fListOfStrawHits = new TObjArray();
  fTimePeak        = NULL;
  fUseStereoHits   = 0;

  fListOfTracks    = new TObjArray();
}

//_____________________________________________________________________________
TTrkVisNode::~TTrkVisNode() {
  delete fArc;
  
  delete fListOfStrawHits;

  fListOfTracks->Delete();
  delete fListOfTracks;
}


//_____________________________________________________________________________
void TTrkVisNode::Paint(Option_t* Option) {
  //

  TStnVisManager* vm = TStnVisManager::Instance();

  const char* view = vm->GetCurrentView();

  if      (strstr(view,"trkxy") != 0) PaintXY  (Option);
  else if (strstr(view,"trkrz") != 0) PaintRZ  (Option);
  else if (strstr(view,"cal"  ) != 0) PaintCal (Option);
  else {
				// what is the default?
    printf(Form("[TTrkVisNode::Paint] >>> ERROR: Unknown option %s",Option));
  }
}

//-----------------------------------------------------------------------------
int TTrkVisNode::InitEvent() {

  const char* oname = "TTrkVisNode::InitEvent";

  mu2e::GeomHandle<mu2e::TTracker> ttHandle;
  const mu2e::TTracker* tracker = ttHandle.get();

  // Tracker calibration object.
  mu2e::ConditionsHandle<mu2e::TrackerCalibrations> trackCal("ignored");

  fListOfStrawHits->Delete();

  fListOfTracks->Delete();

  const mu2e::StrawHit              *hit;
  const mu2e::StrawHitPosition      *hit_pos;
  const mu2e::StrawHitFlag          *hit_id_word;
  const mu2e::PtrStepPointMCVector  *mcptr;

  TEvdStrawHit                      *evd_straw_hit; 
  const CLHEP::Hep3Vector           *mid, *w; 
  const mu2e::Straw                 *straw; 

  int                               n_straw_hits, display_hit, color, nl, ns; // , ipeak, ihit;
  bool                              isFromConversion, intime;
  size_t                            nmc;
  double                            sigw, vnorm, v, sigr; 
  CLHEP::Hep3Vector                 vx0, vx1, vx2;
//-----------------------------------------------------------------------------
// first, clear the cached hit information from the previous event
//-----------------------------------------------------------------------------
  int nst = fTracker->NStations();
  for (int ist=0; ist<nst; ist++) {
    TEvdStation* station = fTracker->Station(ist);
    int np = station->NPanels();
    for (int ip=0; ip<np; ip++) {
      TEvdPanel* panel = station->Panel(ip);
      nl = panel->NLayers();
      for (int il=0; il<nl; il++) {
	ns = panel->NStraws(il);
	for (int is=0; is<ns; is++) {
	  panel->Straw(il,is)->Clear();
	}
      }
    }
  }

  fListOfTracks->Delete();
//-----------------------------------------------------------------------------
// display hits corresponding to a given time peak, or all hits, 
// if the time peak is not found
//-----------------------------------------------------------------------------
  TEvdStraw* evd_straw;
  n_straw_hits = (*fStrawHitColl)->size();

  for (int ihit=0; ihit<n_straw_hits; ihit++ ) {

    hit         = &(*fStrawHitColl)    ->at(ihit);
    hit_pos     = &(*fStrawHitPosColl) ->at(ihit);
    hit_id_word = &(*fStrawHitFlagColl)->at(ihit);

    straw       = &tracker->getStraw(hit->strawIndex());
    display_hit = 1;

    if (! display_hit)                                      continue;
//-----------------------------------------------------------------------------
// deal with MC information - later
//-----------------------------------------------------------------------------
    mcptr = &(*fMcPtrColl)->at(ihit);
	
    // Get the straw information:

    mid   = &straw->getMidPoint();
    w     = &straw->getDirection();

    isFromConversion = false;

    nmc = mcptr->size();
    for (size_t j=0; j<nmc; ++j ){
      const mu2e::StepPointMC& step = *mcptr->at(j);
      art::Ptr<mu2e::SimParticle> const& simptr = step.simParticle();
      mu2e::SimParticleCollection::key_type trackId(step.trackId());
      const mu2e::SimParticle* sim  = simptr.operator ->();
      if (sim == NULL) {
	printf(">>> ERROR: %s sim == NULL\n",oname);
      }
      else {
	if ( sim->fromGenerator() ){
	  mu2e::GenParticle* gen = (mu2e::GenParticle*) &(*sim->genParticle());
	  //	    if ( gen->generatorId() == mu2e::GenId::conversionGun ){
	  if ( gen->generatorId() == mu2e::GenId::StoppedParticleReactionGun ){
	    isFromConversion = true;
	    break;
	  }
	}
      }
    }
//-----------------------------------------------------------------------------
// Position along wire displayed are StrawHitPosition's
//-----------------------------------------------------------------------------	
    v     = trackCal->TimeDiffToDistance( straw->index(), hit->dt() );
    vnorm = v/straw->getHalfLength();

    if (fUseStereoHits) {
//-----------------------------------------------------------------------------
// new default, hit position errors come from StrawHitPositionCollection
//-----------------------------------------------------------------------------
      sigw  = hit_pos->posRes(mu2e::StrawHitPosition::phi); 
      sigr  = hit_pos->posRes(mu2e::StrawHitPosition::rho); 
    }
    else {
//-----------------------------------------------------------------------------
// old default, draw semi-random errors
//-----------------------------------------------------------------------------
      sigw  = trackCal->TimeDivisionResolution( straw->index(), vnorm )/2.; // P.Murat
      sigr  = 5.; // in mm
    }
	
    intime = fabs(hit->time()-fEventTime) < fTimeWindow;
	
    if ( isFromConversion ) {
      if (intime) color = kRed;
      else        color = kBlue;
    }
    else          color = kBlack;
//-----------------------------------------------------------------------------
// add pointer to the hit to the straw 
//-----------------------------------------------------------------------------
    int mask = 0;
    if (intime          ) mask |= TEvdStrawHit::kInTimeBit;
    if (isFromConversion) mask |= TEvdStrawHit::kConversionBit;
    
    evd_straw_hit = new TEvdStrawHit(hit,
				     hit_pos->pos().x(),
				     hit_pos->pos().y(),
				     hit_pos->pos().z(),
				     w->x(),w->y(),
				     sigw,sigr,
				     mask,color);
    int ist, ip, il, is;

    ist = straw->id().getDevice();
    ip  = straw->id().getSector();
    il  = straw->id().getLayer();
    is  = straw->id().getStraw();
      
    evd_straw = fTracker->Station(ist)->Panel(ip)->Straw(il,is);
      



    evd_straw->AddHit(evd_straw_hit);

    fListOfStrawHits->Add(evd_straw_hit);
  }
//-----------------------------------------------------------------------------
// now initialize tracks
//-----------------------------------------------------------------------------
  TEvdTrack                *trk;
  const KalRep             *krep;  
  const mu2e::TrkStrawHit  *track_hit;

  int ntrk = (*fKalRepPtrCollection)->size();
  
  for (int i=0; i<ntrk; i++) {
    krep = (*fKalRepPtrCollection)->at(i).get();
    trk  = new TEvdTrack(i,krep);
//-----------------------------------------------------------------------------
// add hits
//-----------------------------------------------------------------------------
    const TrkHotList*        hot_list = krep->hotList();
    for (TrkHotList::hot_iterator it=hot_list->begin(); it<hot_list->end(); it++) {
      track_hit = (const mu2e::TrkStrawHit*) &(*it);
      TEvdTrkStrawHit* h = new TEvdTrkStrawHit(track_hit);
      trk->AddHit(h);
    }

    fListOfTracks->Add(trk);
  }
  
  return 0;
}

//-----------------------------------------------------------------------------
// draw reconstructed tracks and straw hits
//-----------------------------------------------------------------------------
void TTrkVisNode::PaintXY(Option_t* Option) {

  double        time;
  int           station, display_hit, ntrk;
  TEvdStrawHit  *hit;

  const mu2e::StrawHit   *straw_hit;
  const mu2e::Straw      *straw; 

  //  const char* view = TVisManager::Instance()->GetCurrentView();

  mu2e::GeomHandle<mu2e::TTracker> ttHandle;
  const mu2e::TTracker* tracker = ttHandle.get();

  TStnVisManager* vm = TStnVisManager::Instance();

  int ipeak = vm->TimePeak();

  if (ipeak >= 0) {
    if ((*fCalTimePeakColl) != NULL) {
      int ntp = (*fCalTimePeakColl)->size();
      if (ipeak < ntp) fTimePeak = &(*fCalTimePeakColl)->at(ipeak);
      else             fTimePeak = NULL;
    }
  }

  int nhits = fListOfStrawHits->GetEntries();
  for (int i=0; i<nhits; i++) {
    hit       = (TEvdStrawHit*) fListOfStrawHits->At(i);
    straw_hit = hit->StrawHit();
    straw     = &tracker->getStraw(straw_hit->strawIndex());
    station   = straw->id().getDevice();
    time      = straw_hit->time();

    if ((station >= vm->MinStation()) && (station <= vm->MaxStation())) { 
      display_hit = 1;
      if (fTimePeak != NULL) {
	if ((time < fTimePeak->TMin()) || (time > fTimePeak->TMax())) display_hit = 0;
      }
      if (display_hit) {
	hit->Paint(Option);
      }
    }
  }
//-----------------------------------------------------------------------------
// now - tracks
//-----------------------------------------------------------------------------
  TEvdTrack* evd_trk;
  TAnaDump::Instance()->printKalRep(0,"banner");

  ntrk = fListOfTracks->GetEntriesFast();
  for (int i=0; i<ntrk; i++ ) {
    evd_trk = (TEvdTrack*) fListOfTracks->At(i);
    evd_trk->Paint(Option);
  }

  gPad->Modified();
}

//-----------------------------------------------------------------------------
void TTrkVisNode::PaintRZ(Option_t* Option) {
  int             ntrk, nhits;
  TEvdTrack*      evd_trk;

  //  TStnVisManager* vm = TStnVisManager::Instance();

  fTracker->PaintRZ(Option);
//-----------------------------------------------------------------------------
// do not draw all straw hits - just redraw straws in different color instead
//-----------------------------------------------------------------------------
//   int nhits = fListOfStrawHits->GetEntries();
//   for (int i=0; i<nhits; i++) {
//     hit       = (TEvdStrawHit*) fListOfStrawHits->At(i);

//     if ((station >= vm->MinStation()) && (station <= vm->MaxStation())) continue;
//     if ((time    <  vm->TMin()      ) || (time     > vm->TMax()      )) continue; 

//     hit->Paint(Option);
//   }
//-----------------------------------------------------------------------------
// display tracks and track hits
//-----------------------------------------------------------------------------
  ntrk = fListOfTracks->GetEntriesFast();
  for (int i=0; i<ntrk; i++ ) {
    evd_trk = (TEvdTrack*) fListOfTracks->At(i);
    evd_trk->Paint(Option);

    nhits = evd_trk->NHits();
    for (int ih=0; ih<nhits; ih++) {
      TEvdTrkStrawHit* hit = evd_trk->Hit(ih);
      hit->PaintRZ(Option);
    }
  }

  gPad->Modified();
}

//_____________________________________________________________________________
void TTrkVisNode::PaintCal(Option_t* option) {
  // so far do nothing
 
}

//_____________________________________________________________________________
Int_t TTrkVisNode::DistancetoPrimitive(Int_t px, Int_t py) {
  return 9999;
}

//_____________________________________________________________________________
Int_t TTrkVisNode::DistancetoPrimitiveXY(Int_t px, Int_t py) {

  Int_t dist = 9999;

  static TVector3 global;
  //  static TVector3 local;

  //  Double_t    dx1, dx2, dy1, dy2, dx_min, dy_min, dr;

  global.SetXYZ(gPad->AbsPixeltoX(px),gPad->AbsPixeltoY(py),0);

  return dist;
}

//_____________________________________________________________________________
Int_t TTrkVisNode::DistancetoPrimitiveRZ(Int_t px, Int_t py) {
  return 9999;
}

