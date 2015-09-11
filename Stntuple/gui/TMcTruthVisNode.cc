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

#include "Stntuple/gui/TMcTruthVisNode.hh"
#include "Stntuple/gui/TStnVisManager.hh"

#include "Mu2eUtilities/inc/SimParticlesWithHits.hh"
#include "Mu2eUtilities/inc/SortedStepPoints.hh"
#include "Mu2eUtilities/inc/TrackTool.hh"

#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"

#include <vector>

ClassImp(TMcTruthVisNode)

//_____________________________________________________________________________
TMcTruthVisNode::TMcTruthVisNode(const char* name): TVisNode(name) {
  fArc                  = new TArc;
  fEventTime            = 0;
  fTimeWindow           = 1.e6;
  fSimParticlesWithHits = 0;
  fGraph                = 0;
  fSteps                = 0;
  fGenpColl             = 0;
  fMinEnergyDep         = 0.;
}

//_____________________________________________________________________________
TMcTruthVisNode::~TMcTruthVisNode() {
  delete fArc;
}

//_____________________________________________________________________________
int TMcTruthVisNode::InitEvent() {
  return 0;
}


//_____________________________________________________________________________
void TMcTruthVisNode::Paint(Option_t* option) {
  //
				// parse option list
  const char* view = TVisManager::Instance()->GetCurrentView();

  if      (strstr(view,"trkxy" ) != 0) PaintXY (option);
  if      (strstr(view,"cal"   ) != 0) PaintCal(option);
  else {
				// what is the default?
    //    Warning("Paint",Form("Unknown option %s",option));
  }
}


//_____________________________________________________________________________
void TMcTruthVisNode::PaintXY(Option_t* option) {

  const char oname[] = "TMcTruthVisNode::PaintXY";
  // draw particles

  mu2e::SimParticleCollection::key_type key(1);
  const mu2e::GenParticle* gen;
 
  mu2e::SimParticleInfo const* info((*fSimParticlesWithHits)->findOrNull(key));

  const mu2e::StepPointMC  *firstStep; 
  static mu2e::TrackTool   *tt(0), *tg(0), *tmid(0);
  TArrow arrow;

  std::vector<double> xStep, yStep;

    if ( !info ) {
      printf("[%s] >>> ERROR: SimParticleInfo missing, SKIPPING EVENT\n",oname);
      firstStep = 0;
      //      return;
    }
    else {
      firstStep = &info->firstStepPointMCinTracker();
    }
    
    mu2e::SortedStepPoints sortedSteps(key,*(*fSteps));
    
    const mu2e::StepPointMC* midStep = & sortedSteps.middleByZ();
    
    if ((midStep == 0) || (firstStep ==0)) {
      printf("[%s] >>> ERROR : firstStep = %8p midstep = %8p, BAIL OUT\n",
	     oname,firstStep, midStep);
      //      goto END_OF_ROUTINE;
    }

    CLHEP::Hep3Vector pos1, mom1, pos2, mom2;
    if (firstStep ) {
      pos1 = firstStep->position();
      mom1 = firstStep->momentum();
    }
    else {
      pos1.set(1.,1.,1.);
      mom1.set(1.,1.,1.);
    }

    if (midStep) {
      pos2 = midStep->position();
      mom2 = midStep->momentum();
    }
    else {
      pos2.set(1.,1.,1.);
      mom2.set(1.,1.,1.);
    }
					// The generated particle - the first one! ... inherited...
    gen = &((*fGenpColl)->at(0));

    if (tt) {
      delete tt;
      delete tg;
      delete tmid;
    }

    tt   = new mu2e::TrackTool(gen->pdgId(),-1.,pos1           ,mom1           ,1.,CLHEP::Hep3Vector());
    tg   = new mu2e::TrackTool(gen->pdgId(),-1.,gen->position(),gen->momentum(),1.,CLHEP::Hep3Vector());
    tmid = new mu2e::TrackTool(gen->pdgId(),-1.,pos2           ,mom2           ,1.,CLHEP::Hep3Vector());

    int npt = (*fSteps)->size();
    for (int ipt=0; ipt<npt; ++ipt){
      mu2e::StepPointMC const& step =  (*fSteps)->at(ipt);
      if (step.totalEDep() > fMinEnergyDep) {
        xStep.push_back( step.position().x() );
        yStep.push_back( step.position().y() );
      }
    }

    fArc->SetFillStyle(0);
//-----------------------------------------------------------------------------
// draw StepPointMC's - as a graph - fix that!
//-----------------------------------------------------------------------------
    if (xStep.size() > 0) {
      if (fGraph) delete fGraph;
      fGraph = new TGraph( xStep.size(), &xStep[0], &yStep[0]);
      fGraph->SetMarkerStyle(kOpenCircle);
      fGraph->SetMarkerColor(kRed+1);
      fGraph->SetMarkerSize(0.6);
      fGraph->Draw("PSAME");
    }
//-----------------------------------------------------------------------------
// red marker and arrow: first point on the track.
//-----------------------------------------------------------------------------
    double xf1 = pos1.x();
    double yf1 = pos1.y();
    TGraph genPoint( 1, &xf1, &yf1 );
    genPoint.SetMarkerColor(kRed);
    genPoint.SetMarkerSize(1.0);
    genPoint.SetMarkerStyle(kFullCircle);
    genPoint.Draw("PSAME");
    
    double arrowLength(200.);
    double xf2 = xf1 + arrowLength*mom1.x()/mom1.perp();
    double yf2 = yf1 + arrowLength*mom1.y()/mom1.perp();
    arrow.SetLineColor(kRed);
    arrow.DrawArrow( xf1, yf1, xf2, yf2, 0.01, ">");
      
    double d0x  = tt->d0x();
    double d0y  = tt->d0y();
    double d0x2 =  tt->d0x() + arrowLength*tt->u0();
    double d0y2 =  tt->d0y() + arrowLength*tt->v0();
      
    arrow.SetLineColor(kBlue);
    arrow.DrawArrow(d0x, d0y, d0x2, d0y2, 0.01, ">");

	
  gPad->Modified();
}

//_____________________________________________________________________________
void TMcTruthVisNode::PaintCal(Option_t* option) {
  // draw calorimeter

  //  TEllipse e;

//   art::ServiceHandle<mu2e::GeometryService> geom;
//   mu2e::GeomHandle<mu2e::TTracker> tt;


  gPad->Modified();
}


//_____________________________________________________________________________
void TMcTruthVisNode::PaintRZ(Option_t* option) {
  // draw calorimeter
}

//_____________________________________________________________________________
Int_t TMcTruthVisNode::DistancetoPrimitive(Int_t px, Int_t py) {
  return 9999;
}

//_____________________________________________________________________________
Int_t TMcTruthVisNode::DistancetoPrimitiveXY(Int_t px, Int_t py) {

  Int_t dist = 9999;

  static TVector3 global;
  static TVector3 local;

  //  Double_t    dx1, dx2, dy1, dy2, dx_min, dy_min, dr;

  global.SetXYZ(gPad->AbsPixeltoX(px),gPad->AbsPixeltoY(py),0);

  return dist;
}

//_____________________________________________________________________________
Int_t TMcTruthVisNode::DistancetoPrimitiveRZ(Int_t px, Int_t py) {
  return 9999;
}

