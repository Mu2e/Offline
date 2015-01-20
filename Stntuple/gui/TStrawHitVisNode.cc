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
#include "TMath.h"
#include "TBox.h"
#include "TObjArray.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"

#include "Stntuple/gui/TEvdStrawHit.hh"
#include "Stntuple/gui/TStrawHitVisNode.hh"
#include "Stntuple/gui/TStnVisManager.hh"

#include "TrackerGeom/inc/Straw.hh"
#include "TTrackerGeom/inc/TTracker.hh"

#include "Mu2eUtilities/inc/TrackTool.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"

#include "RecoDataProducts/inc/StrawHit.hh"

#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"

ClassImp(TStrawHitVisNode)

//_____________________________________________________________________________
TStrawHitVisNode::TStrawHitVisNode(const char* name): TVisNode(name) {
  fArc        = new TArc;
  fEventTime  = 0;
  fTimeWindow = 1.e6;

  fListOfStrawHits = new TObjArray();
  fTimePeak        = NULL;
  fUseStereoHits   = 0;
}

//_____________________________________________________________________________
TStrawHitVisNode::~TStrawHitVisNode() {
  delete fArc;
}

//_____________________________________________________________________________
int TStrawHitVisNode::InitEvent() {
  const char* oname = "TStrawHitVisNode::InitEvent";

  mu2e::GeomHandle<mu2e::TTracker> ttHandle;
  const mu2e::TTracker* tracker = ttHandle.get();

  // Tracker calibration object.
  mu2e::ConditionsHandle<mu2e::TrackerCalibrations> trackCal("ignored");

  fListOfStrawHits->Delete();

  const mu2e::StrawHit              *hit;
  const mu2e::StrawHitPosition      *hitpos;
  const mu2e::StrawHitFlag          *hit_id_word;
  const mu2e::PtrStepPointMCVector  *mcptr;

  TEvdStrawHit                      *evd_straw_hit; 
  const CLHEP::Hep3Vector           *mid, *w; 
  const mu2e::Straw                 *straw; 

  int                               n_straw_hits, display_hit, color; // , ipeak, ihit;
  bool                              isFromConversion, intime;
  size_t                            nmc;
  double                            sigv, vnorm, v, sigr; 
  CLHEP::Hep3Vector                 vx0, vx1, vx2; 
//-----------------------------------------------------------------------------
// display hits corresponding to a given time peak, or all if it is not found
//-----------------------------------------------------------------------------
  n_straw_hits = (*fStrawHitColl)->size();

  for (int ihit=0; ihit<n_straw_hits; ihit++ ) {

    hit         = &(*fStrawHitColl)    ->at(ihit);
    hitpos      = &(*fStrawHitPosColl) ->at(ihit);
    hit_id_word = &(*fStrawHitFlagColl)->at(ihit);

    straw = &tracker->getStraw(hit->strawIndex());

    //    int station = straw->id().getDevice();

    display_hit = 1;

//     if (fDisplayBackgroundHits == false) {
//       if (! hit_id_word->hasAllProperties(fGoodHitMask)) display_hit = 0;
//       if (hit_id_word->hasAnyProperty(fBadHitMask)     ) display_hit = 0;
//     }

    if (display_hit) {
      // StrawHitMCTruth const& truth(hitsTruth.at(ihit));

      evd_straw_hit = new TEvdStrawHit(hit);

      fListOfStrawHits->Add(evd_straw_hit);

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
	
      // Position along wire, from delta t.
	
      v     = trackCal->TimeDiffToDistance( straw->index(), hit->dt() );
      vnorm = v/straw->getHalfLength();

      if (fUseStereoHits) {
//-----------------------------------------------------------------------------
// new default, hit position errors come from StrawHitPositionCollection
//-----------------------------------------------------------------------------
	sigv  = hitpos->posRes(mu2e::StrawHitPosition::phi); 
	sigr  = hitpos->posRes(mu2e::StrawHitPosition::rho); 
      }
      else {
//-----------------------------------------------------------------------------
// old default, draw semi-random errors
//-----------------------------------------------------------------------------
	sigv  = trackCal->TimeDivisionResolution( straw->index(), vnorm )/2.; // P.Murat
	sigr  = 5.; // in mm
      }
	
      vx0 = (*mid) + v   *(*w);
      vx1 = vx0    + sigv*(*w);
      vx2 = vx0    - sigv*(*w);
	
      intime = fabs(hit->time()-fEventTime) < fTimeWindow;
	
      if ( isFromConversion ) {
	if (intime) color = kRed;
	else        color = kBlue;
      }
      else          color = kBlack;

      evd_straw_hit->SetColor   (color );
      evd_straw_hit->SetPos     (vx0.x(),vx0.y());
      evd_straw_hit->SetStrawDir(w->x() ,w->y());
      evd_straw_hit->SetSigW    (sigv  );
      evd_straw_hit->SetSigR    (sigr  );

      int mask = 0;
      if (intime          ) mask |= TEvdStrawHit::kInTimeBit;
      if (isFromConversion) mask |= TEvdStrawHit::kConversionBit;

      evd_straw_hit->SetMask    (mask  );  
    }
  }




  return 0;
}


//_____________________________________________________________________________
void TStrawHitVisNode::Paint(Option_t* option) {
  //
  const char oname[] = "TStrawHitVisNode::Paint";

				// parse option list

  const char* view = TVisManager::Instance()->GetCurrentView();

  if      (strstr(view,"trkxy" ) != 0) PaintXY (option);
  if      (strstr(view,"cal"   ) != 0) PaintCal(option);
  else {
    printf("[%s] >>> ERROR: unknown VIEW %s, DO NOTHING\n",oname,view);
  }

  gPad->Modified();
}


//_____________________________________________________________________________
void TStrawHitVisNode::PaintXY(Option_t* Option) {
  // draw calorimeter

  double        time;
  int           station, display_hit;
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
}

//_____________________________________________________________________________
void TStrawHitVisNode::PaintCal(Option_t* Option) {
  // draw calorimeter

}


//_____________________________________________________________________________
void TStrawHitVisNode::PaintRZ(Option_t* option) {
  // draw calorimeter
}

//_____________________________________________________________________________
Int_t TStrawHitVisNode::DistancetoPrimitive(Int_t px, Int_t py) {
  return 9999;
}

//_____________________________________________________________________________
Int_t TStrawHitVisNode::DistancetoPrimitiveXY(Int_t px, Int_t py) {

  Int_t dist = 9999;

  static TVector3 global;
  //  static TVector3 local;

  //  Double_t    dx1, dx2, dy1, dy2, dx_min, dy_min, dr;

  global.SetXYZ(gPad->AbsPixeltoX(px),gPad->AbsPixeltoY(py),0);

  return dist;
}

//_____________________________________________________________________________
Int_t TStrawHitVisNode::DistancetoPrimitiveRZ(Int_t px, Int_t py) {
  return 9999;
}

