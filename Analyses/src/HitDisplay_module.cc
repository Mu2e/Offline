//
// A sandbox for playing with tracks, including transformations to different representations.
// This is not production code but feel free to look at it.
//
// $Id: HitDisplay_module.cc,v 1.6 2011/05/20 19:18:44 wb Exp $
// $Author: wb $
// $Date: 2011/05/20 19:18:44 $
//
// Original author Rob Kutschke.
//

// C++ includes.
#include <iostream>
#include <string>

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes.
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "Mu2eUtilities/inc/SimParticlesWithHits.hh"
#include "Mu2eUtilities/inc/SortedStepPoints.hh"
#include "Mu2eUtilities/inc/TrackTool.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"
#include "Mu2eUtilities/inc/resolveDPIndices.hh"
#include "TTrackerGeom/inc/TTracker.hh"

#include "ToyDP/inc/DPIndexVectorCollection.hh"
#include "ToyDP/inc/GenParticleCollection.hh"
#include "ToyDP/inc/SimParticleCollection.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "ToyDP/inc/StrawHitCollection.hh"
#include "ToyDP/inc/StrawHitMCTruthCollection.hh"

// ROOT includes
#include "TApplication.h"
#include "TArc.h"
#include "TArrow.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TLine.h"
#include "TNtuple.h"

// Other includes
#include "CLHEP/Units/SystemOfUnits.h"

using namespace std;
using CLHEP::Hep3Vector;

namespace mu2e {

  class HitDisplay : public art::EDAnalyzer {
  public:
    explicit HitDisplay(fhicl::ParameterSet const& pset);
    virtual ~HitDisplay(){}

    void beginJob( );
    void analyze( art::Event const& e );

  private:

    // The module label of this instance of this module.
    std::string moduleLabel_;

    // Label of the modules that created the data products.
    std::string generatorModuleLabel_;
    std::string g4ModuleLabel_;
    std::string hitMakerModuleLabel_;

    // Name of the tracker StepPoint collection
    std::string trackerStepPoints_;

    // Cuts used inside SimParticleWithHits:
    //  - drop hits with too little energy deposited.
    //  - drop SimParticles with too few hits.
    double minEnergyDep_;
    size_t minHits_;

    bool doDisplay_;

    auto_ptr<TApplication> application_;
    TDirectory*   directory_;
    TCanvas*      canvas_;

    TH1F* _hnHits;
    TH1F* _hEnergyDep;
    TH1F* _hDeltaT;
    TH1F* _hTime;
    TH1F* _hx;
    TH1F* _hxnorm;

    TH1F* _hMissDist;

    TNtuple* _ntTrack;
    TNtuple* _ntHit;

  };

  HitDisplay::HitDisplay(fhicl::ParameterSet const& pset):
    moduleLabel_(pset.get<string>("@module_label")),
    generatorModuleLabel_(pset.get<std::string>("generatorModuleLabel")),
    g4ModuleLabel_(pset.get<std::string>("g4ModuleLabel")),
    hitMakerModuleLabel_(pset.get<std::string>("hitMakerModuleLabel")),
    trackerStepPoints_(pset.get<std::string>("trackerStepPoints")),
    minEnergyDep_(pset.get<double>("minEnergyDep")),
    minHits_(pset.get<unsigned>("minHits")),
    doDisplay_(pset.get<bool>("doDisplay",true)),
    application_(0),
    directory_(0),
    canvas_(0),
    _hnHits(0),
    _hEnergyDep(0),
    _hDeltaT(0),
    _hTime(0),
    _hx(0),
    _hxnorm(0),
    _hMissDist(0),
    _ntTrack(0),
    _ntHit(0){
  }

  void HitDisplay::beginJob(){

    // Get access to the TFile service.
    art::ServiceHandle<art::TFileService> tfs;

    // Create some 1D histograms.
    _hnHits     =  tfs->make<TH1F>( "hnHits",      "StrawHits per Event",              100,    0.,    200. );
    _hEnergyDep =  tfs->make<TH1F>( "hEnergyDep",  "Energy Deposition per Hit;(keV)",  100,    0.,     20. );
    _hDeltaT    =  tfs->make<TH1F>( "hDeltaT",     "Delta(time);(ns)",                 100,   -5.,      5. );
    _hTime      =  tfs->make<TH1F>( "hTime",       "Time;(ns)",                        100,    0.,   2000. );
    _hx         =  tfs->make<TH1F>( "hx",          "Displacement;(mm)",                100, -1000.,  1000. );
    _hxnorm     =  tfs->make<TH1F>( "hxnorm",      "Displacement;(HalfLength)",        100,    -1.,     1. );
    _hMissDist  =  tfs->make<TH1F>( "hMissDist",   "Distance to wire",                 100,   -50.,    50. );
    _ntTrack    =  tfs->make<TNtuple>( "ntTrack",  "TrackInfo", "d0gen:d0first" );
    _ntHit      =  tfs->make<TNtuple>( "ntHit",    "HitInfo",   "dca:z" );

    if ( !doDisplay_ ) return;

    // If needed, create the ROOT interactive environment. See note 1.
    if ( !gApplication ){
      int    tmp_argc(0);
      char** tmp_argv(0);
      application_ = auto_ptr<TApplication>(new TApplication( "noapplication", &tmp_argc, tmp_argv ));
    }

    // Create a canvas with a guaranteed unique name; the module label is unique within a job.
    TString name  = "canvas_"     + moduleLabel_;
    TString title = "Canvas for " + moduleLabel_;
    int window_size(800);
    canvas_ = tfs->make<TCanvas>(name,title,window_size,window_size);

    directory_ = gDirectory;

  }


  void
  HitDisplay::analyze(art::Event const& event) {

    // Tracker geometry.
    GeomHandle<TTracker> ttracker;

    // Geometry of tracker envelope.
    TubsParams envelope(ttracker->getTrackerEnvelopeParams());

    // Tracker calibration object.
    ConditionsHandle<TrackerCalibrations> trackCal("ignored");

    // Get information from the event.
    art::Handle<GenParticleCollection> gensHandle;
    event.getByLabel(generatorModuleLabel_,gensHandle);
    GenParticleCollection const& gens = *gensHandle;

    art::Handle<SimParticleCollection> simsHandle;
    event.getByLabel(g4ModuleLabel_,simsHandle);
    SimParticleCollection const& sims = *simsHandle;

    art::Handle<StrawHitCollection> hitsHandle;
    event.getByLabel(hitMakerModuleLabel_,hitsHandle);
    StrawHitCollection const& hits = *hitsHandle;

    /*
      art::Handle<StrawHitMCTruthCollection> hitsTruthHandle;
      event.getByLabel(hitMakerModuleLabel_,hitsTruthHandle);
      StrawHitMCTruthCollection const& hitsTruth = *hitsTruthHandle;
    */

    art::Handle<DPIndexVectorCollection> mcptrHandle;
    event.getByLabel(hitMakerModuleLabel_,"StrawHitMCPtr",mcptrHandle);
    DPIndexVectorCollection const& hits_mcptr = *mcptrHandle;

    art::Handle<StepPointMCCollection> stepsHandle;
    event.getByLabel(g4ModuleLabel_,trackerStepPoints_,stepsHandle);
    StepPointMCCollection const& steps = *stepsHandle;

    // Construct an object that ties together all of the simulated particle and hit info.
    SimParticlesWithHits simsInfo( event,
                                   g4ModuleLabel_,
                                   hitMakerModuleLabel_,
                                   trackerStepPoints_,
                                   minEnergyDep_,
                                   minHits_ );

    SimParticleCollection::key_type key(1);
    SimParticleInfo const* info(simsInfo.findOrNull(key));
    if ( !info ){
      cout << "Skipping event: " << event.id() << endl;
      return;
    }

    StepPointMC const& firstStep(info->firstStepPointMCinTracker());

    SortedStepPoints sortedSteps(key,steps);
    StepPointMC const& midStep(sortedSteps.middleByZ());

    // The generated particle.
    GenParticle const& gen = gens.at(0);
    TrackTool tt  ( gen.pdgId(), -1., firstStep.position(), firstStep.momentum(), 1., Hep3Vector() );
    TrackTool tg  ( gen.pdgId(), -1., gen.position(), gen.momentum(), 1., Hep3Vector() );
    TrackTool tmid( gen.pdgId(), -1., midStep.position(), midStep.momentum(), 1., Hep3Vector() );

    _hnHits->Fill( hits.size() );

    vector<double> xStep,yStep;
    for ( size_t ipt=0; ipt<steps.size(); ++ipt){
      StepPointMC const& step =  steps.at(ipt);
      if ( step.totalEDep() > minEnergyDep_ ) {
        xStep.push_back( step.position().x() );
        yStep.push_back( step.position().y() );
      }
    }

    static TLine*  line  = new TLine;
    static TArc*   arc   = new TArc;
    static TArrow* arrow = new TArrow;

    arc->SetFillStyle(0);

    if ( doDisplay_ ) {
      canvas_->cd(0);
      canvas_->Clear();

      // Draw the frame
      double plotLimits(850.);
      canvas_->DrawFrame(-plotLimits,-plotLimits,plotLimits,plotLimits);

      // Draw the inner and outer arc of the tracker.
      arc->SetLineColor(kBlack);
      arc->DrawArc(0.,0., envelope.outerRadius());
      arc->DrawArc(0.,0., envelope.innerRadius());
      arc->SetLineColor(kRed);
      arc->DrawArc( tt.xc(), tt.yc(), tt.rho());
      arc->SetLineColor(kMagenta);
      arc->DrawArc( tmid.xc(), tmid.yc(), tmid.rho());
      arc->SetLineColor(kRed);

    }

    float ntTrack[_ntTrack->GetNvar()];
    ntTrack[0] = tg.da();
    ntTrack[1] = tt.da();
    _ntTrack->Fill(ntTrack);

    // Loop over all straw hits.
    for ( size_t ihit=0; ihit<hits.size(); ++ihit ) {

      // Data and MC truth for this hit.
      StrawHit        const&   hit(hits.at(ihit));
      //StrawHitMCTruth const& truth(hitsTruth.at(ihit));
      DPIndexVector   const& mcptr(hits_mcptr.at(ihit));

      _hEnergyDep->Fill( hit.energyDep()/CLHEP::keV );

      // Skip hits with too little energy deposited in the straw.
      if ( hit.energyDep() < minEnergyDep_ ){
        continue;
      }

      // Get the straw information:
      const Straw&             straw = ttracker->getStraw( hit.strawIndex() );
      const CLHEP::Hep3Vector& mid   = straw.getMidPoint();
      const CLHEP::Hep3Vector& w     = straw.getDirection();

      // Is this hit from a conversion electron?
      bool isFromConversion(false);
      for ( size_t j=0; j<mcptr.size(); ++j ){
        StepPointMC const& step = *resolveDPIndex<StepPointMCCollection>( event, mcptr.at(j) );
        SimParticle const& sim  = sims[step.trackId()];
        if ( sim.fromGenerator() ){
          GenParticle const& gen = gens.at(sim.generatorIndex());
          if ( gen.generatorId() == GenId::conversionGun ){
            isFromConversion = true;
            break;
          }
        }
      }

      _hDeltaT->Fill( hit.dt() );
      _hTime->Fill( hit.time() );

      Hep3Vector pos( tt.positionAtZ( mid.z() ) );
      Hep3Vector mom( tt.momentumAtZ( mid.z() ) );
      TwoLinePCA pca( pos, mom.unit(), mid, w );
      double dca = pca.dca();

      Hep3Vector posMid( tmid.positionAtZ( mid.z() ) );
      Hep3Vector momMid( tmid.momentumAtZ( mid.z() ) );
      TwoLinePCA pcaMid( posMid, momMid.unit(), mid, w);

      /*
      double dcaMid = pcaMid.dca();
      cout << "Compare: "
           << mid.z() << " "
           << dca << " "
           << dcaMid << " | "
           << dca-dcaMid << " | "
           << mom.unit().dot(w) << " "
           << momMid.unit().dot(w) << " | "
           << pos << " "
           << posMid <<  " "
           << endl;
      */

      if ( isFromConversion ) {
        float ntHit[_ntHit->GetNvar()];
        _hMissDist->Fill( dca );
        ntHit[0] = dca;
        ntHit[1] = mid.z();
        _ntHit->Fill(ntHit);
      }

      // Position along wire, from delta t.
      double v     = trackCal->TimeDiffToDistance( straw.Index(), hit.dt() );
      double vnorm = v/straw.getHalfLength();
      double sigv = trackCal->TimeDivisionResolution( straw.Index(), vnorm );

      _hx->Fill(v);
      _hxnorm->Fill(vnorm);

      CLHEP::Hep3Vector x0 = mid + v*w;
      CLHEP::Hep3Vector x1 = x0 + sigv*w;
      CLHEP::Hep3Vector x2 = x0 - sigv*w;

      if ( doDisplay_  ){
        if ( isFromConversion ) line->DrawLine( x1.x(), x1.y(), x2.x(), x2.y() );
      }

    }


    if ( doDisplay_  ){

      // Draw the generated hits.
      TGraph graph( xStep.size(), &xStep[0], &yStep[0]);
      graph.SetMarkerStyle(kOpenTriangleUp);
      graph.Draw("PSAME");

      // Draw the first point on the track.
      double xf1 = firstStep.position().x();
      double yf1 = firstStep.position().y();
      TGraph genPoint( 1, &xf1, &yf1 );
      cout << "Size: " << genPoint.GetMarkerSize() << endl;
      genPoint.SetMarkerColor(kRed);
      genPoint.SetMarkerSize(1.5);
      genPoint.SetMarkerStyle(kFullCircle);
      genPoint.Draw("PSAME");

      CLHEP::Hep3Vector const& v(firstStep.momentum());
      double arrowLength(200.);
      double xf2 = xf1 + arrowLength*v.x()/v.perp();
      double yf2 = yf1 + arrowLength*v.y()/v.perp();
      arrow->SetLineColor(kRed);
      arrow->DrawArrow( xf1, yf1, xf2, yf2, 0.01, ">");

      double d0x  = tt.d0x();
      double d0y  = tt.d0y();
      double d0x2 =  tt.d0x() + arrowLength*tt.u0();
      double d0y2 =  tt.d0y() + arrowLength*tt.v0();

      double dot = tt.d0x()*tt.u0() + tt.d0y()*tt.v0();

      cout << "D0x etc: "
           << d0x << " "
           << d0y << " "
           << tt.u0() << " "
           << tt.v0() << " "
           << dot
           << endl;

      arrow->SetLineColor(kBlue);
      arrow->DrawArrow(tt.d0x(), tt.d0y(), d0x2, d0y2, 0.01, ">");


      // Plot origin.
      double xo = 0.;
      double yo = 0.;
      TGraph genPointo( 1, &xo, &yo );
      genPointo.SetMarkerColor(kBlue);
      genPointo.SetMarkerSize(2.5);
      genPointo.SetMarkerStyle(kPlus);
      genPointo.Draw("PSAME");



      canvas_->Modified();
      canvas_->Update();

      cerr << "Double click in the Canvas " << moduleLabel_ << " to continue:" ;
      gPad->WaitPrimitive();
      cerr << endl;
    }

    cout << "TubsParams: "
         << envelope.innerRadius() << " "
         << envelope.outerRadius() << " "
         << envelope.zHalfLength() << " "
         << endl;

  } // end of ::analyze.

}

using mu2e::HitDisplay;
DEFINE_ART_MODULE(HitDisplay);
