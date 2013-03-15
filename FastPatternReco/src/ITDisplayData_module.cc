
// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "ITrackerGeom/inc/Cell.hh"
//#include "DataProducts/inc/DPIndexVectorCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "ITrackerGeom/inc/ITracker.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "FastPatternReco/inc/ITHitPerTrackData.hh"
#include "FastPatternReco/inc/GenTrackData.hh"
#include "RecoDataProducts/inc/TrackerHitByID.hh"

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Provenance.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

// Root includes.
#include "TApplication.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TClonesArray.h"
#include "TEllipse.h"
#include "TLatex.h"

// Other external packages
#include "CLHEP/Vector/TwoVector.h"

// C++ includes.
#include <iostream>
#include <string>
#include <memory>
#include <map>
#include <utility>

using namespace std;

namespace mu2e {

  class Straw;

  class ITDisplayData : public art::EDAnalyzer {
  public:
    
    explicit ITDisplayData(fhicl::ParameterSet const& pset);
    virtual ~ITDisplayData() {delete _fakeCanvas; }

    virtual void beginJob();
    void endJob();

    // This is called for each event.
    void analyze(art::Event const& e);

  private:

    // Start: run time parameters

    // The module label of this module.
    std::string _moduleLabel;

    // Module label of the g4 module that made the hits.
    std::string _g4ModuleLabel;

    // Name of the tracker StepPoint collection
    std::string _trackerStepPoints;

    // Label of the module that made the hits.
    std::string _makerModuleLabel;

    // Label of the generator.
    std::string _generatorModuleLabel;

    // Label of the process for the remaping of the Tracker hit by Cell/Straw ID.
    std::string _mapTrackerHitByID;

    // Use efficincy value to simulate the proton rejection by ADC value
    float _protonRejecEff;
    float _pulseCut;
    bool  _useProtonRejec;
    bool  _perfectRejec;

    // drift vevocity
    double _driftVelocity;

    int   _hitClClSizeCut;

//    // Number of events to accumulate between prompts.
//    int _nAccumulate;

    // End: run time parameters

    // Pointers to histograms, ntuples, TGraphs.
    TH2F*         _hHitAtCenter;
    TH2F*         _hHitAtProperZ;
    TH2F*         _hHitStereoPl;
    TH2F*         _hHitStereoMi;
    TH2F*         _hCMapHitStereoPl;
    TH2F*         _hCMapHitStereoMi;
    TH2F*         _hRCMapHitStereoPl;
    TH2F*         _hRCMapHitStereoMi;
    TH1F*         _hHitZ;
    TH1I*         _hHitNloop;
    TH1F*         _hMCHitDTime;
    TH1F*         _hHitTime;
    TH1F*         _hHitClustTime;
    TH2F*         _hHitClustLayerPl;
    TH2F*         _hHitClustLayerMin;
    TCanvas*      _canvas;
    TCanvas*      _canvasPl;
    TCanvas*      _canvasClL;
    TCanvas*      _fakeCanvas;

    int   ntimeBin;
    float maxTimeHist; //ns
    float timeBinDim;  //ns


    // Some ugly but necessary ROOT related bookkeeping:

    // The job needs exactly one instance of TApplication.  See note 1.
    unique_ptr<TApplication> _application;

    // Save directory from beginJob so that we can go there in endJob. See note 3.
    TDirectory* _directory;

  };

  ITDisplayData::ITDisplayData(fhicl::ParameterSet const& pset) :

    // Run time parameters
    _moduleLabel(pset.get<string>("module_label")),
    _g4ModuleLabel(pset.get<string>("g4ModuleLabel")),
    _trackerStepPoints(pset.get<string>("trackerStepPoints","tracker")),
    _makerModuleLabel(pset.get<string>("makerModuleLabel")),
    _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
    _mapTrackerHitByID(pset.get<string>("mapTrackerHitByID")),
    _protonRejecEff(pset.get<float>("protonRejecEff",-1.0)),
    _pulseCut(pset.get<float>("pulseCut",-1.0)),
   /*_nAccumulate(pset.get<int>("nAccumulate",20)),*/
    _driftVelocity(pset.get<double>("driftVelocity",0.035)),   // mm/ns
    _hitClClSizeCut(pset.get<int>("hitClClSizeCut",15)),
    // ROOT objects that are the main focus of this example.
    _hHitAtCenter(0),
    _hHitAtProperZ(0),
    _hHitStereoPl(0),
    _hHitStereoMi(0),
    _hCMapHitStereoPl(0),
    _hCMapHitStereoMi(0),
    _hRCMapHitStereoPl(0),
    _hRCMapHitStereoMi(0),
    _hHitZ(0),
    _hHitNloop(0),
    _hMCHitDTime(0),
    _hHitTime(0),
    _hHitClustTime(0),
    _hHitClustLayerPl(0),
    _hHitClustLayerMin(0),
    _canvas(0),
    _canvasPl(0),
    _canvasClL(0),
    _fakeCanvas(0),

    ntimeBin(0),
    maxTimeHist(2500.0),
    timeBinDim(10.0),
    // Some ugly but necessary ROOT related bookkeeping.
    _application(nullptr),
    _directory(0){

          if (_protonRejecEff>0.0 && _pulseCut>0.0) {
                  _useProtonRejec=true;
                  if (_protonRejecEff>=1.0) {
                          _protonRejecEff=1.0;
                          _perfectRejec=true;
                  }
                  else {
                          _perfectRejec=false;
                  }
          }
          else {
                  _useProtonRejec=false;
                  _perfectRejec=false;
          }

 }

  void ITDisplayData::beginJob(){

    // Get access to the TFile service and save current directory for later use.
    art::ServiceHandle<art::TFileService> tfs;

    ntimeBin    = (int) maxTimeHist/timeBinDim;

    // Create a histogram.
    _hHitAtCenter  = tfs->make<TH2F>( "hHitAtCenter",  "Hits per Event at z=0", 1500, -75.0, 75.0, 1500, -75.0, 75.0 );
    _hHitAtProperZ = tfs->make<TH2F>( "hHitAtProperZ", "Hits per Event at true z", 1500, -75.0, 75.0, 1500, -75.0, 75.0 );
    _hHitStereoPl  = tfs->make<TH2F>( "hHitStereoPl",  "Stereo+ Hits per Event at z=0", 1500, -75.0, 75.0, 1500, -75.0, 75.0 );
    _hHitStereoMi  = tfs->make<TH2F>( "hHitStereoMi",  "Stereo- Hits per Event at z=0", 1500, -75.0, 75.0, 1500, -75.0, 75.0 );
    _hCMapHitStereoPl  = tfs->make<TH2F>( "hCMapHitStereoPl",  "Conformal Mapping of Stereo+ Hits per Event at z=0", 1000, -0.5, 0.5, 1000, -0.5, 0.5 );
    _hCMapHitStereoMi  = tfs->make<TH2F>( "hCMapHitStereoMi",  "Conformal Mapping of Stereo- Hits per Event at z=0", 1000, -0.5, 0.5, 1000, -0.5, 0.5 );
    _hRCMapHitStereoPl  = tfs->make<TH2F>( "hRCMapHitStereoPl",  "Relative Conformal Mapping of Stereo+ Hits per Event at z=0", 2000, -2, 2, 2000, -2, 2 );
    _hRCMapHitStereoMi  = tfs->make<TH2F>( "hRCMapHitStereoMi",  "Relative Conformal Mapping of Stereo- Hits per Event at z=0", 2000, -2, 2, 2000, -2, 2 );
    _hHitZ         = tfs->make<TH1F>( "hHitZ",         "Z of the Hits per Event", 1200, -300.0, 300.0 );
    _hHitNloop     = tfs->make<TH1I>( "hHitNloop",     "Track number of loops per Event", 5, 0, 5 );
    _hMCHitDTime   = tfs->make<TH1F>( "hMCHitDTime",   "Delta Time of the MC Step hits per Event", 200, 0.0, 50.0 );
    _hHitTime      = tfs->make<TH1F>( "hHitTime",      "Time of the Hits per Event", ntimeBin, 0.0, maxTimeHist );
    _hHitClustTime = tfs->make<TH1F>( "hHitClustTime", "Cluster of Time of the Hits per Event", ntimeBin, 0.0, maxTimeHist );
    _hHitClustLayerPl = tfs->make<TH2F>( "hHitClustLayerPl", "Cluster of consecutive Hits per Layer per Event", 50, 0.5, 50.5, 100, 0.5, 100.5 );
    _hHitClustLayerMin = tfs->make<TH2F>( "hHitClustLayerMin", "Cluster of consecutive Hits per Layer per Event", 50, 0.5, 50.5, 100, 0.5, 100.5 );

    _hHitAtCenter ->SetXTitle("cm");
    _hHitAtCenter ->SetYTitle("cm");
    _hHitAtProperZ->SetXTitle("cm");
    _hHitAtProperZ->SetYTitle("cm");
    _hHitStereoPl ->SetXTitle("cm");
    _hHitStereoPl ->SetYTitle("cm");
    _hHitStereoMi ->SetXTitle("cm");
    _hHitStereoMi ->SetYTitle("cm");
    _hCMapHitStereoPl ->SetXTitle("1/cm");
    _hCMapHitStereoPl ->SetYTitle("1/cm");
    _hCMapHitStereoMi ->SetXTitle("1/cm");
    _hCMapHitStereoMi ->SetYTitle("1/cm");
    _hRCMapHitStereoPl ->SetXTitle("1/cm");
    _hRCMapHitStereoPl ->SetYTitle("1/cm");
    _hRCMapHitStereoMi ->SetXTitle("1/cm");
    _hRCMapHitStereoMi ->SetYTitle("1/cm");
    _hHitZ        ->SetXTitle("cm");
    _hHitNloop    ->SetXTitle("N loops");
    _hMCHitDTime  ->SetXTitle("ns");
    _hHitTime     ->SetXTitle("ns");
    _hHitClustTime->SetXTitle("ns");
    _hHitClustLayerPl->SetXTitle("Layer");
    _hHitClustLayerPl->SetYTitle("Cl size");
    _hHitClustLayerMin->SetXTitle("Layer");
    _hHitClustLayerMin->SetYTitle("Cl size");

    // If needed, create the ROOT interactive environment. See note 1.
    if ( !gApplication ){
      int    tmp_argc(0);
      char** tmp_argv(0);
      _application = unique_ptr<TApplication>(new TApplication( "noapplication", &tmp_argc, tmp_argv ));
    }

    // Create a canvas with a unique name.  See note 2.
    TString name  = "canvas_"     + _moduleLabel;
    TString title = "Canvas for " + _moduleLabel;
    int window_size(860);
    _canvas = tfs->make<TCanvas>(name,title,window_size,window_size);
    _canvas->Divide(2,2);
    name  = "canvasPl_"     + _moduleLabel;
    title = "Canvas for Plots " + _moduleLabel;
    _canvasPl = tfs->make<TCanvas>(name,title,window_size,window_size);
    _canvasPl->Divide(2,2);
    name  = "canvasClL_"     + _moduleLabel;
    title = "Canvas for contiguous hits cluster " + _moduleLabel;
    _canvasClL = tfs->make<TCanvas>(name,title,window_size,window_size);
    _canvasClL->Divide(2,2);
    _fakeCanvas = new TCanvas("canvas_Fake","double click for next event",300,100);
//    _fakeCanvas = tfs->make<TCanvas>("canvas_Fake","double click for next event",400,200);

    // Draw the still empty histogram. It will be updated later.
    _canvas->cd(1);
    _hHitAtCenter->SetStats(kFALSE);
    _hHitAtCenter->Draw();
    _canvas->cd(2);
    _hHitAtProperZ->SetStats(kFALSE);
    _hHitAtProperZ->Draw();
    _canvas->cd(3);
    _hHitStereoMi->SetStats(kFALSE);
    _hHitStereoMi->Draw();
    _canvas->cd(4);
    _hHitStereoPl->SetStats(kFALSE);
    _hHitStereoPl->Draw();

    _canvasPl->cd(1);
    _hHitTime->SetStats(kTRUE);
    _hHitTime->Draw();
    _canvasPl->cd(2);
    _hHitClustTime->SetStats(kTRUE);
    _hHitClustTime->Draw();
    _canvasPl->cd(3);
    _hHitZ->SetStats(kFALSE);
    _hHitZ->Draw();
    _canvasPl->cd(4);
    _hHitNloop->SetStats(kFALSE);
    _hHitNloop->Draw();
//    _hMCHitDTime->SetStats(kFALSE);
//    _hMCHitDTime->Draw();

    // See note 3.
    _directory = gDirectory;

  }

  void ITDisplayData::analyze(art::Event const& event) {

/*
    // Ask the event to give us a "handle" to the requested hits.
    art::Handle<StepPointMCCollection> hitsHandle;
    event.getByLabel(_g4ModuleLabel,_trackerStepPoints,hitsHandle);
    StepPointMCCollection const& hits = *hitsHandle;

    // Fill histogram with number of hits per event.
    _hHitAtCenter->Fill(hits.size());

    // Periodically update the displayed histogram.

    if ( event.id().event()%_nAccumulate==0 ){
      _canvas->Modified();
      _canvas->Update();

      cerr << "Double click in the Canvas " << _moduleLabel << " to continue:" ;
      _canvas->WaitPrimitive();
      cerr << endl;

    }
*/

    //--------------------------------------------

    _hHitAtCenter->Reset();
    _hHitAtProperZ->Reset();
    _hHitStereoMi->Reset();
    _hHitStereoPl->Reset();
    _hCMapHitStereoPl->Reset();
    _hCMapHitStereoMi->Reset();
    _hRCMapHitStereoPl->Reset();
    _hRCMapHitStereoMi->Reset();
    _hHitTime->Reset();
    _hMCHitDTime->Reset();
    _hHitClustTime->Reset();
    _hHitClustLayerPl->Reset();
    _hHitClustLayerMin->Reset();
    _hHitZ->Reset();
    _hHitNloop->Reset();

    const Tracker& tracker = getTrackerOrThrow();
    const ITracker &itr = static_cast<const ITracker&>( tracker );
    CellGeometryHandle *itwp = itr.getCellGeometryHandle();

    itwp->SelectCell(0,0,0);
    float intTimeWind = 0.0;
    cout<<"-- cel 0,0,0 rad "<<itwp->GetCellRad()<<endl;
    intTimeWind += itwp->GetCellRad();
    itwp->SelectCell(itr.nSuperLayers()-1,0,0);
    cout<<"-- cel "<<itr.nSuperLayers()-1<<",0,0 rad "<<itwp->GetCellRad()<<endl;
    intTimeWind += itwp->GetCellRad();
    intTimeWind *= 0.5;
    intTimeWind = intTimeWind/_driftVelocity + 20.0;   //max drift time + max TOF of a signal electron
    cout<<"-- intTimeWind "<<intTimeWind<<endl;

    double rIn  = itr.r0()+itr.getWalls()->find(Wall::inner)->second->getTotalThickness();
    double rOut = itr.rOut()-itr.getWalls()->find(Wall::outer)->second->getTotalThickness();
    rIn/=CLHEP::cm;
    rOut/=CLHEP::cm;

    double tmpX, tmpY, cmapU, cmapV;

    art::Handle<StrawHitCollection> pdataHandle;
    event.getByLabel(_makerModuleLabel,pdataHandle);
    StrawHitCollection const* hits = pdataHandle.product();

    // Get the persistent data about the StrawHitsMCTruth.
    art::Handle<StrawHitMCTruthCollection> truthHandle;
    event.getByLabel(_makerModuleLabel,truthHandle);
    StrawHitMCTruthCollection const* hits_truth = truthHandle.product();

    // Get the persistent data about pointers to StepPointMCs
    art::Handle<PtrStepPointMCVectorCollection> mcptrHandle;
    event.getByLabel(_makerModuleLabel,"StrawHitMCPtr",mcptrHandle);
    PtrStepPointMCVectorCollection const* hits_mcptr = mcptrHandle.product();

    // Get the persistent data about the StepPointMCs. More correct implementation
    // should look for product ids in DPIndexVectorCollection, rather than
    // use producer name directly ("g4run").

    if (!(hits->size() == hits_truth->size() &&
          hits_mcptr->size() == hits->size() ) ) {
      throw cet::exception("RANGE")
        << "Strawhits: " << hits->size()
        << " MCTruthStrawHits: " << hits_truth->size()
        << " MCPtr: " << hits_mcptr->size();
    }

    // Get handles to the generated and simulated particles.
    art::Handle<GenParticleCollection> genParticles;
    event.getByLabel(_generatorModuleLabel, genParticles);

//    // Handle to information about G4 physical volumes.
//    art::Handle<PhysicalVolumeInfoCollection> volumes;
//    event.getRun().getByType(volumes);
//
//    // Some files might not have the SimParticle and volume information.
//    bool haveSimPart = ( simParticles.isValid() && volumes.isValid() );
//
//    // Other files might have empty collections.
//    if ( haveSimPart ){
//      haveSimPart = !(simParticles->empty() || volumes->empty());
//    }


    art::Handle<TrackerHitByID> hitByIDHandle;
    event.getByLabel(_mapTrackerHitByID,hitByIDHandle);
    TrackerHitByID const* hitByID = hitByIDHandle.product();
    TrackerHitByID::const_iterator hitByID_it;
    std::pair<TrackerHitByID::const_iterator, TrackerHitByID::const_iterator> rangeHitByID_it;
    std::set<size_t> hitLoked;


    size_t nStrawPerEvent = hits->size();

    TClonesArray *hitDraws = new TClonesArray("TEllipse");
    hitDraws->ExpandCreateFast(nStrawPerEvent);
    TClonesArray *hitDrawsAtZ = new TClonesArray("TEllipse");
    hitDrawsAtZ->ExpandCreateFast(nStrawPerEvent);

    //vector<int> idHitStereoMin;
    bool *isStereoMin = new bool [nStrawPerEvent];

    std::map<SimParticleCollection::key_type,std::pair<CLHEP::Hep3Vector,CLHEP::HepLorentzVector> > genElectrons;
    std::map<SimParticleCollection::key_type,std::pair<CLHEP::Hep3Vector,CLHEP::Hep3Vector> > firstHitOfGenEl;
    std::map<SimParticleCollection::key_type,double > firstHitTime;
    std::map<SimParticleCollection::key_type,std::pair<CLHEP::Hep3Vector,CLHEP::Hep3Vector> >::iterator firstHitOfGenEl_it;
    std::map<SimParticleCollection::key_type,double >::iterator firstHitTime_it;

    std::map<SimParticleCollection::key_type, GenTrackData<ITHitPerTrackData> > genTracks;
    std::map<SimParticleCollection::key_type, GenTrackData<ITHitPerTrackData> >::iterator genTracks_it;

    float aveZ=0.0;
    double mchittime=0.0;

    bool overlapped = false;
    bool isFirst = false;

    typedef std::vector<std::pair<size_t, std::pair<double,double> > > largeHitCls;
    largeHitCls potentialtrackSP_Pl;
    largeHitCls potentialtrackSP_Min;

    for (size_t i=0; i<nStrawPerEvent; ++i) {

      // Access data
      StrawHit        const&      hit(hits->at(i));
      StrawHitMCTruth const&    truth(hits_truth->at(i));
      PtrStepPointMCVector   const&    mcptr(hits_mcptr->at(i));

      //double hitEnergy = hit.energyDep();

      //Skip the straw if the energy of the hit is smaller than the minimum required
      //if (hitEnergy < _minimumEnergy) continue;

      if (_useProtonRejec) {
              if (hit.energyDep()>_pulseCut) {
                      continue;
                      //if (_perfectRejec) continue;
                      //else if (_randFlat.fire()<_protonRejecEff) continue;
              }
      }


      //Get hit straw
      StrawIndex si = hit.strawIndex();
      const Straw & str = tracker.getStraw(si);
      const Cell & cell = static_cast<const Cell&>( str );


      // cout << "Getting informations about cells" << endl;

      //int sid = cell.Id().getCell();
      //int lid = cell.Id().getLayer();
      //int did = cell.Id().getLayerId().getSuperLayer();
      //itwp->SelectCell(did,lid,sid);
      itwp->SelectCellDet(si.asUint());
      int sid = itwp->GetWire();
      int lid = itwp->GetCelRing();
      int did = itwp->GetSuperLayer();

      const CLHEP::Hep3Vector stMidPoint3 = itwp->GetCellCenter();//str.getMidPoint();

      //time of the hit
      double hitTime = hit.time();

      //direction of the straw
      const CLHEP::Hep3Vector stDirection3 = itwp->GetCellDirection();//str.getDirection();

      // cout << "Reading MCtruth info" << endl;

      // Get MC truth data
      //double driftTime = truth.driftTime();
      double driftDistance = truth.driftDistance();

      //Position along the wire using mctruth info
      double vMC = truth.distanceToMid();

      const CLHEP::Hep3Vector MCHitPoint = stMidPoint3 + (vMC/stDirection3.mag())*stDirection3;

      tmpX = stMidPoint3.getX()/CLHEP::cm;
      tmpY = stMidPoint3.getY()/CLHEP::cm;
      cmapV = 1.0/(tmpX*tmpX+tmpY*tmpY);
      cmapU = tmpX*cmapV;
      cmapV*= tmpY;
      ((TEllipse *) hitDraws->At(i))->SetX1(tmpX);
      ((TEllipse *) hitDraws->At(i))->SetY1(tmpY);
      ((TEllipse *) hitDraws->At(i))->SetR1(driftDistance/CLHEP::cm);
      ((TEllipse *) hitDraws->At(i))->SetR2( ( driftDistance/cos(cell.getWire()->getEpsilon()) )/CLHEP::cm );
//      ((TEllipse *) hitDraws->At(i))->SetTheta( TMath::ATan2(stMidPoint3.getY(),stMidPoint3.getX())/CLHEP::degree );
      ((TEllipse *) hitDraws->At(i))->SetTheta( stMidPoint3.phi()/CLHEP::degree );
      ((TEllipse *) hitDraws->At(i))->SetFillStyle(0);
      ((TEllipse *) hitDraws->At(i))->SetLineWidth(2);
      int color = 2;
      isStereoMin[i]=false;
      if (cell.getWire()->getEpsilon()<0.0) {
              color = 4;
              isStereoMin[i]=true;
              //idHitStereoMin.push_back(i);
              _hCMapHitStereoMi->Fill(cmapU,cmapV);
      }
      else {
              _hCMapHitStereoPl->Fill(cmapU,cmapV);
      }
      ((TEllipse *) hitDraws->At(i))->SetLineColor(color);
//      ((TEllipse *) hitDraws->At(i))->Draw("same");

      _hHitTime->Fill(hitTime);


      //common index for vectors
      //size_t trackIdx(0);

      aveZ=0.0;

      overlapped = false;
      isFirst = false;

      for (size_t j = 0; j < mcptr.size(); ++j) {

        StepPointMC const& mchit = *mcptr[j];

        // The simulated particle that made this hit.
        SimParticleCollection::key_type trackId(mchit.trackId());
        art::Ptr<SimParticle> const& simptr = mchit.simParticle();
        art::Ptr<SimParticle>::key_type simKey(simptr.key());
        //SimParticle const& sim = simParticles->at(trackId);
        SimParticle const& sim = *simptr;

        isFirst=( j==0 );
        if ( isFirst ){
                for ( size_t jj = 1; jj < mcptr.size(); ++jj) {
//                        if ( trackId != (*mchits)[mcptr[jj].index].trackId() ) {
                          if ( simKey!=mcptr[jj]->simParticle().key() && trackId != (*(mcptr[jj])).trackId() ) {
                                overlapped=true;
                                break;
                        }
                }
        }


        //cout<<"Hit time :"<<mchit.time()<<" proper "<<mchit.properTime()<<endl;
        if ( sim.pdgId()==11 && sim.isPrimary() ){
          if ( genElectrons.find(trackId)==genElectrons.end() ) {
            genElectrons.insert ( pair<SimParticleCollection::key_type,std::pair<CLHEP::Hep3Vector,CLHEP::HepLorentzVector> > ( trackId, std::pair<CLHEP::Hep3Vector,CLHEP::HepLorentzVector>( sim.startPosition(), sim.startMomentum() ) ) );
            firstHitOfGenEl.insert ( pair<SimParticleCollection::key_type,std::pair<CLHEP::Hep3Vector,CLHEP::Hep3Vector> > ( trackId, std::pair<CLHEP::Hep3Vector,CLHEP::Hep3Vector>( mchit.position(), mchit.momentum() ) ) );
            firstHitTime.insert ( pair<SimParticleCollection::key_type,double>(trackId,mchit.time()) );
          }
          else {
            firstHitTime_it = firstHitTime.find(trackId);
            firstHitOfGenEl_it = firstHitOfGenEl.find(trackId);
            if (mchit.time()<firstHitTime_it->second){
               firstHitTime.erase(firstHitTime_it);
               firstHitOfGenEl.erase(firstHitOfGenEl_it);
               firstHitOfGenEl.insert ( pair<SimParticleCollection::key_type,std::pair<CLHEP::Hep3Vector,CLHEP::Hep3Vector> > ( trackId, std::pair<CLHEP::Hep3Vector,CLHEP::Hep3Vector>( mchit.position(), mchit.momentum() ) ) );
               firstHitTime.insert ( pair<SimParticleCollection::key_type,double>(trackId,mchit.time()) );
            }
          }

          genTracks_it=genTracks.find(trackId);
          mchittime=mchit.time();
          if ( genTracks_it==genTracks.end() ) {
                  genTracks.insert( pair<SimParticleCollection::key_type,GenTrackData<ITHitPerTrackData> > ( trackId, GenTrackData<ITHitPerTrackData>(trackId, /*std::pair<CLHEP::Hep3Vector,CLHEP::HepLorentzVector>*/make_pair( sim.startPosition(), sim.startMomentum() ) ) ) );
                  (--(genTracks.end()))->second.addHitData( ITHitPerTrackData( i,mchittime,overlapped,isFirst,sid,lid,did,mchit.position(),mchit.momentum() ) );
          }
          else {
                  if (genTracks_it->second.isNotPresent(i) ) genTracks_it->second.addHitData( ITHitPerTrackData( i,mchittime,overlapped,isFirst,sid,lid,did,mchit.position(),mchit.momentum() ) );
          }
        }

        _hHitAtProperZ->Fill( mchit.position().getX()/CLHEP::cm, mchit.position().getY()/CLHEP::cm );
        aveZ+=mchit.position().getZ();
      }
      aveZ/= ((float) mcptr.size());
      float wpos[3];
      //itwp->WirePosAtZ(aveZ, wpos);
      itwp->WirePosAtZ(MCHitPoint.getZ(), wpos);

      aveZ/=CLHEP::cm;
      _hHitZ->Fill(aveZ);

      ((TEllipse *) hitDrawsAtZ->At(i))->SetX1(wpos[0]/CLHEP::cm);
      ((TEllipse *) hitDrawsAtZ->At(i))->SetY1(wpos[1]/CLHEP::cm);
      ((TEllipse *) hitDrawsAtZ->At(i))->SetR1(driftDistance/CLHEP::cm);
      ((TEllipse *) hitDrawsAtZ->At(i))->SetR2( ( driftDistance/cos(cell.getWire()->getEpsilon()) )/CLHEP::cm);
      ((TEllipse *) hitDrawsAtZ->At(i))->SetTheta( TMath::ATan2(wpos[1],wpos[0])/CLHEP::degree );
//      ((TEllipse *) hitDrawsAtZ->At(i))->SetTheta( stMidPoint3.phi()/CLHEP::degree );
      ((TEllipse *) hitDrawsAtZ->At(i))->SetFillStyle(0);
      ((TEllipse *) hitDrawsAtZ->At(i))->SetLineWidth(2);
//      color = 2;
//      if (cell.getWire()->getEpsilon()<0.0) color = 4;
      ((TEllipse *) hitDrawsAtZ->At(i))->SetLineColor(color);


      //-------------------------------------------------
      //Hit cluster per layer analysis

      if ( hitLoked.find(i)==hitLoked.end() ){
              unsigned long int tmpIndex;
              //tmpIndex = si.asUint();

              //cout<<"Hit n "<<i<<" Cell Id "<<tmpIndex<<" hits with same id:"<<endl;
              //rangeHitByID_it= hitByID->equal_range(tmpIndex);
              //art::Ptr<mu2e::StrawHit> c;
              //c.key()
              //for (hitByID_it=rangeHitByID_it.first; hitByID_it!=rangeHitByID_it.second; ++hitByID_it){
              //        cout<<"\t id "<<hitByID_it->first<<" hit n "<<hitByID_it->second.key()<<endl;
              //}
              int nCellPerLayer = itwp->GetITLayer()->nCells();
              int tmpSid=0, lost=0, hitClSizeInLayer=1;
              bool lostHit=true;
              //int lid = cell.Id().getLayer();
              //int did = cell.Id().getLayerId().getSuperLayer();
              for (int icel=1; icel<nCellPerLayer; icel++){
                      tmpSid = (sid+icel)%nCellPerLayer;
                      tmpIndex = itwp->computeDet(did,lid,tmpSid);
                      if ( hitByID->count(tmpIndex)>0 ){
                              lostHit=true;
                              rangeHitByID_it= hitByID->equal_range(tmpIndex);
                              for (hitByID_it=rangeHitByID_it.first; hitByID_it!=rangeHitByID_it.second; ++hitByID_it){
                                      if ( hitLoked.find(hitByID_it->second.key())==hitLoked.end() && fabs(hitByID_it->second->time()-hitTime)<=intTimeWind ){
                                              ++hitClSizeInLayer;
                                              hitLoked.insert( hitByID_it->second.key() );
                                              lostHit=false;
                                              break;
                                      }
                              }
                              if (lostHit) { ++lost; }
                      }
                      else { ++lost; }
                      if (lost>0) { break; }
              }
              lost=0;
              for (int icel=1; icel<nCellPerLayer; icel++){
                      tmpSid = sid-icel;
                      if (tmpSid<0) tmpSid+=nCellPerLayer;
                      tmpSid=tmpSid%nCellPerLayer;
                      tmpIndex = itwp->computeDet(did,lid,tmpSid);
                      if ( hitByID->count(tmpIndex)>0 ){
                              lostHit=true;
                              rangeHitByID_it= hitByID->equal_range(tmpIndex);
                              for (hitByID_it=rangeHitByID_it.first; hitByID_it!=rangeHitByID_it.second; ++hitByID_it){
                                      if ( hitLoked.find(hitByID_it->second.key())==hitLoked.end() && fabs(hitByID_it->second->time()-hitTime)<=intTimeWind ){
                                              ++hitClSizeInLayer;
                                              hitLoked.insert( hitByID_it->second.key() );
                                              lostHit=false;
                                              break;
                                      }
                              }
                              if (lostHit) { ++lost; }
                      }
                      else { ++lost; }
                      if (lost>0) { break; }
              }

              if (isStereoMin[i]){
                      _hHitClustLayerMin->Fill(did,hitClSizeInLayer);
                      if (hitClSizeInLayer>=_hitClClSizeCut) {
                              potentialtrackSP_Min.push_back(largeHitCls::value_type(i,std::pair<double,double>(tmpX,tmpY) ) );
                      }
              } else {
                      _hHitClustLayerPl->Fill(did,hitClSizeInLayer);
                      if (hitClSizeInLayer>=_hitClClSizeCut) {
                              potentialtrackSP_Pl.push_back(largeHitCls::value_type(i,std::pair<double,double>(tmpX,tmpY) ) );
                      }
              }
      }
      //-------------------------------------------------


    }

    //-------------------------------------------------
    cout<<"***** nHit "<<nStrawPerEvent<<endl;
    //Conformal Mapping relative to large Hit cluster
    double invTmpR2;
    //std::vector<double *> Uvec_Pl, Vvec_Pl, Uvec_Min, Vvec_Min;
    //std::vector<TGraph *> rCMapGraphs_Pl, rCMapGraphs_Min;
    std::vector<TH2F *> rCMap_Pl, rCMap_Min;
    int jCl = 0;
    for (largeHitCls::iterator lrgCl_it = potentialtrackSP_Min.begin(); lrgCl_it != potentialtrackSP_Min.end();  ++lrgCl_it) {
            //Uvec_Min.push_back(new double[nStrawPerEvent-1]);
            //Vvec_Min.push_back(new double[nStrawPerEvent-1]);
            rCMap_Min.push_back(new TH2F( Form("hRCMapHitStereoMin_%i",jCl), "Relative Conformal Mapping of Stereo- Hits per Event at z=0", 2000, -2, 2, 2000, -2, 2 ));
            for (size_t i=0; i<nStrawPerEvent; ++i) {
                    if(i==lrgCl_it->first || !isStereoMin[i]) continue;
                    // Access data
                    StrawHit        const&      hit(hits->at(i));

                    //Get hit straw
                    StrawIndex si = hit.strawIndex();
                    //const Straw & str = tracker.getStraw(si);
                    //const Cell & cell = static_cast<const Cell&>( str );

                    // cout << "Getting informations about cells" << endl;

                    //int sid = cell.Id().getCell();
                    //int lid = cell.Id().getLayer();
                    //int did = cell.Id().getLayerId().getSuperLayer();
                    //itwp->SelectCell(did,lid,sid);
                    itwp->SelectCellDet(si.asUint());
                    /*int sid = itwp->GetWire();
                    int lid = itwp->GetCelRing();
                    int did = itwp->GetSuperLayer();*/

                    const CLHEP::Hep3Vector stMidPoint3 = itwp->GetCellCenter();//str.getMidPoint();
                    tmpX = stMidPoint3.getX()/CLHEP::cm;
                    tmpY = stMidPoint3.getY()/CLHEP::cm;
                    cmapU = tmpX-lrgCl_it->second.first;
                    cmapV = tmpY-lrgCl_it->second.second;
                    invTmpR2  = 1.0/(cmapU*cmapU+cmapV*cmapV);
                    cmapU*= invTmpR2;
                    cmapV*= invTmpR2;

                    rCMap_Min.at(jCl)->Fill(cmapU,cmapV);
            }
            ++jCl;
    }
    cout<<"---- potentialtrackSP_Min size "<<jCl<<endl;
    jCl=0;
    for (largeHitCls::iterator lrgCl_it = potentialtrackSP_Pl.begin(); lrgCl_it != potentialtrackSP_Pl.end();  ++lrgCl_it) {
            //Uvec_Pl.push_back(new double[nStrawPerEvent-1]);
            //Vvec_Pl.push_back(new double[nStrawPerEvent-1]);
            rCMap_Pl.push_back(new TH2F( Form("hRCMapHitStereoPl_%i",jCl), "Relative Conformal Mapping of Stereo- Hits per Event at z=0", 2000, -2, 2, 2000, -2, 2 ));
            for (size_t i=0; i<nStrawPerEvent; ++i) {
                    if(i==lrgCl_it->first || isStereoMin[i]) continue;
                    // Access data
                    StrawHit        const&      hit(hits->at(i));

                    //Get hit straw
                    StrawIndex si = hit.strawIndex();
                    //const Straw & str = tracker.getStraw(si);
                    //const Cell & cell = static_cast<const Cell&>( str );

                    // cout << "Getting informations about cells" << endl;

                    //int sid = cell.Id().getCell();
                    //int lid = cell.Id().getLayer();
                    //int did = cell.Id().getLayerId().getSuperLayer();
                    //itwp->SelectCell(did,lid,sid);
                    itwp->SelectCellDet(si.asUint());
                    /*int sid = itwp->GetWire();
                    int lid = itwp->GetCelRing();
                    int did = itwp->GetSuperLayer();*/

                    const CLHEP::Hep3Vector stMidPoint3 = itwp->GetCellCenter();//str.getMidPoint();
                    tmpX = stMidPoint3.getX()/CLHEP::cm;
                    tmpY = stMidPoint3.getY()/CLHEP::cm;
                    cmapU = tmpX-lrgCl_it->second.first;
                    cmapV = tmpY-lrgCl_it->second.second;
                    invTmpR2  = 1.0/(cmapU*cmapU+cmapV*cmapV);
                    cmapU*= invTmpR2;
                    cmapV*= invTmpR2;

                    rCMap_Pl.at(jCl)->Fill(cmapU,cmapV);
            }
            ++jCl;
    }
    cout<<"---- potentialtrackSP_Pl size "<<jCl<<endl;

    //-------------------------------------------------


    TEllipse innerWall (0.0,0.0,rIn,rIn);
    innerWall.SetFillStyle(0);
    innerWall.SetLineWidth(1.5);

    TEllipse outerWall (0.0,0.0,rOut,rOut);
    outerWall.SetFillStyle(0);
    outerWall.SetLineWidth(1.5);

    std::map<SimParticleCollection::key_type,std::pair<CLHEP::Hep3Vector,CLHEP::HepLorentzVector> >::iterator genEl_it     = genElectrons.begin();
    /*std::map<SimParticleCollection::key_type,std::pair<CLHEP::Hep3Vector,CLHEP::Hep3Vector> >::iterator*/ firstHitOfGenEl_it = firstHitOfGenEl.begin();

    //TClonesArray *genElDraws = new TClonesArray("TEllipse");
    //int nGenEl = genElectrons.size();
    //genElDraws->ExpandCreateFast( nGenEl );
    while ( genEl_it != genElectrons.end() ){
      cout<<"Generated el at "<<genEl_it->second.first<<" of "<<genEl_it->second.second<<endl;
      cout<<"First hit in tracker of gen el at "<<firstHitOfGenEl_it->second.first<<" of "<<firstHitOfGenEl_it->second.second<<endl;
      ++genEl_it;
      ++firstHitOfGenEl_it;
    }


    //TClonesArray *hitDrawsAtZ = new TClonesArray("TEllipse");
    //hitDrawsAtZ->ExpandCreateFast(nStrawPerEvent);
    std::map<SimParticleCollection::key_type, TClonesArray* > genTracksCrircles;
    std::map<SimParticleCollection::key_type, TClonesArray* >::iterator genTracksCrircles_it;
    TEllipse *trckCircl;

    double ptMeV, rho;
    double B=1.0;
    CLHEP::Hep2Vector radDir;
    HepGeom::Point3D<double> CirCenter;

    for ( genTracks_it= genTracks.begin(); genTracks_it!= genTracks.end(); ++genTracks_it ){
//            genTracks_it->second.sort();
            _hHitNloop->Fill (genTracks_it->second.FindNTurns());
            genTracksCrircles.insert( std::pair<SimParticleCollection::key_type, TClonesArray* >( genTracks_it->first,  new TClonesArray("TEllipse") ) );
            ((TClonesArray *) genTracksCrircles[genTracks_it->first])->ExpandCreateFast(genTracks_it->second.getNumOfLoops());
            //cout<<genTracks_it->second;
            for ( unsigned int ilp=0; ilp<genTracks_it->second.getNumOfLoops(); ilp++ ){
                    ITHitPerTrackData hdil = genTracks_it->second.getithLoopHit(ilp);
                    ptMeV = TMath::Sqrt( pow(hdil.hitMomentum[0],2) + pow(hdil.hitMomentum[1],2) );
                    rho   = ptMeV/(B*0.3);
                    cout<<ilp<<" -th loop: p_t "<<ptMeV<<" rho mm "<<rho<<endl;
                    CirCenter.set(hdil.hitPoint.getX(),hdil.hitPoint.getY(),hdil.hitPoint.getZ());
                    radDir.setX(hdil.hitMomentum.getX());
                    radDir.setY(hdil.hitMomentum.getY());
                    radDir.rotate( ( (hdil.hitMomentum.getZ()>=0.0) ? 90.0 : -90.0 )*CLHEP::degree );
                    radDir=radDir.unit();
//                    radDir.rotateZ( ( (radDir[2]>=0.0) ? 90.0 : -90.0 )*CLHEP::degree );
//                    HepGeom::Transform3D CCenterTras = HepGeom::Translate3D(rho*radDir);
//                    CirCenter=CCenterTras*CirCenter;
                    CirCenter=CirCenter+rho*radDir;
                    trckCircl = ((TEllipse *) genTracksCrircles[genTracks_it->first]->At(ilp));
                    trckCircl->SetX1(CirCenter[0]/CLHEP::cm);
                    trckCircl->SetY1(CirCenter[1]/CLHEP::cm);
                    trckCircl->SetR1(rho/CLHEP::cm);
                    trckCircl->SetR2(rho/CLHEP::cm);
                    trckCircl->SetFillStyle(0);
                    trckCircl->SetLineWidth(1.5);
                    if ( genTracks_it->second.getTrkLrntzVec().rho()>=103.0 ) trckCircl->SetLineColor(kGreen+2);
                    else trckCircl->SetLineColor(kMagenta-2);
                    trckCircl->SetLineStyle(4);
                    cout<<" Hit Pos "<<hdil.hitPoint<<" Circ center "<<CirCenter<<endl;

            }
//            for (unsigned int ih=1; ih<genTracks_it->second.getNumOfHit(); ih++){
//                    _hMCHitDTime->Fill( genTracks_it->second.getHit(ih).mcHitTime - genTracks_it->second.getHit(ih-1).mcHitTime );
//            }
    }

    _canvas->cd(1);
    _hHitAtCenter->Draw();
    innerWall.Draw("same");
    outerWall.Draw("same");
    for (size_t i=0; i<nStrawPerEvent; ++i) ((TEllipse *) hitDraws->At(i))->Draw("same");
    _hHitAtCenter->Draw("same");

    _canvas->cd(2);
    _hHitAtProperZ->Draw();
    innerWall.Draw("same");
    outerWall.Draw("same");
    for (size_t i=0; i<nStrawPerEvent; ++i) ((TEllipse *) hitDrawsAtZ->At(i))->Draw("same");
    for ( genTracksCrircles_it=genTracksCrircles.begin(); genTracksCrircles_it!=genTracksCrircles.end(); genTracksCrircles_it++ ){
            for ( unsigned int ilp=0; ilp<((unsigned int)genTracksCrircles_it->second->GetEntries()); ilp++ )  ((TEllipse *) genTracksCrircles_it->second->At(ilp))-> Draw("same");
    }
    _hHitAtProperZ->Draw("same");

    _canvas->cd(3);
    _hHitStereoMi->Draw();
    innerWall.Draw("same");
    outerWall.Draw("same");
    for (size_t i=0; i<nStrawPerEvent; ++i) if ( isStereoMin[i] ) ((TEllipse *) hitDraws->At(i))->Draw("same");
    _hHitStereoMi->Draw("same");

    _canvas->cd(4);
    _hHitStereoPl->Draw();
    innerWall.Draw("same");
    outerWall.Draw("same");
    for (size_t i=0; i<nStrawPerEvent; ++i) if ( !isStereoMin[i] ) ((TEllipse *) hitDraws->At(i))->Draw("same");
    _hHitStereoPl->Draw("same");

    _canvasPl->cd(1);
    //_hHitTime->SetStats(kTRUE);
    _hHitTime->Draw();

    int nbingroup= (int) intTimeWind/timeBinDim;
    float hit=0.0;
    int binHalf = (int) nbingroup/2;
    int jbin;
    //int nBin=_hHitTime->GetNbinsX();
    for (int it=1; it<=ntimeBin; it++) {
      hit=0.0;
      for (int is=0; is<nbingroup; is++) {
        jbin = it -binHalf +is;
        if (jbin<1) continue;
        if (jbin>ntimeBin) break;
        hit+=_hHitTime->GetBinContent(jbin);
      }
      _hHitClustTime->SetBinContent(it,hit);
    }

    _canvasPl->cd(2);
    _hHitClustTime->Draw();

    _canvasPl->cd(3);
    _hHitZ->Draw();

    _canvasPl->cd(4);
    _hHitNloop->Draw();
//    _hMCHitDTime->Draw();

    _canvasClL->cd(1);
    _hHitClustLayerMin->Draw("colz");

    _canvasClL->cd(2);
    _hHitClustLayerPl->Draw("colz");

    _canvasClL->cd(3);
    _hCMapHitStereoMi->Draw("colz");

    _canvasClL->cd(4);
    _hCMapHitStereoPl->Draw("colz");

    TCanvas *tmpRCP_Min = new TCanvas();
    int nlrgCl_Min = rCMap_Min.size();
    int iPad = 0, maxPad = 1;
    if (nlrgCl_Min>1) {
      tmpRCP_Min->Divide(2,nlrgCl_Min/2);
      maxPad = 2*(nlrgCl_Min/2);
      iPad=1;
    }
    int iColor = 0;
    for (std::vector<TH2F*>::iterator hCmapMin_it = rCMap_Min.begin(); hCmapMin_it != rCMap_Min.end(); ++hCmapMin_it ) {
            if (iPad>maxPad) break;
            tmpRCP_Min->cd(iPad);
            (*hCmapMin_it)->SetMarkerStyle(8);
            (*hCmapMin_it)->SetMarkerSize(0.8);
            (*hCmapMin_it)->SetMarkerColor(++iColor);
            (*hCmapMin_it)->Draw();
            ++iPad;
    }

    TCanvas *tmpRCP_Pl = new TCanvas();
    int nlrgCl_Pl = rCMap_Pl.size();
    iPad = 0;
    maxPad = 1;
    if (nlrgCl_Pl>1) {
      tmpRCP_Pl->Divide(2,nlrgCl_Pl/2);
      maxPad = 2*(nlrgCl_Pl/2);
      iPad=1;
    }
    iColor = 0;
    for (std::vector<TH2F*>::iterator hCmapPl_it = rCMap_Pl.begin(); hCmapPl_it != rCMap_Pl.end(); ++hCmapPl_it ) {
            if (iPad>maxPad) break;
            tmpRCP_Pl->cd(iPad);
            (*hCmapPl_it)->SetMarkerStyle(8);
            (*hCmapPl_it)->SetMarkerSize(0.8);
            (*hCmapPl_it)->SetMarkerColor(++iColor);
            (*hCmapPl_it)->Draw();
            ++iPad;
    }

    _canvas->Modified();
    _canvas->Update();

    _canvasPl->Modified();
    _canvasPl->Update();

    _canvasClL->Modified();
    _canvasClL->Update();

    tmpRCP_Min->Modified();
    tmpRCP_Min->Update();

    tmpRCP_Pl->Modified();
    tmpRCP_Pl->Update();


    //-------------------------------------------------

//    TrackerHitByID::const_iterator hitByID_it;
//    const ITracker &itr = static_cast<const ITracker&>( tracker );
//    CellGeometryHandle *itwp = itr.getCellGeometryHandle();
//
//    itwp->SelectCell(0,0,0);
//    float intTimeWind = 0.0;
//    cout<<"-- cel 0,0,0 rad "<<itwp->GetCellRad()<<endl;
//    intTimeWind += itwp->GetCellRad();
//    itwp->SelectCell(itr.nSuperLayers()-1,0,0);
//    cout<<"-- cel "<<itr.nSuperLayers()-1<<",0,0 rad "<<itwp->GetCellRad()<<endl;
//    intTimeWind += itwp->GetCellRad();
//    for ( hitByID_it = hitByID->begin(); hitByID_it!= hitByID->end(); ++hitByID_it ){
//            cout<<"Cell/Straw id "<<hitByID_it->first<<" hit Data :"<<endl;
//            hitByID_it->second->print();
//    }

//    for (int isl=0; isl<itr.nSuperLayers(); isl++) {
//
//    }
//    const boost::shared_array<mu2e::SuperLayer> sprlr=itracker.getSuperLayersArray();
//    for ( int iS=0; iS<itracker.nSuperLayers(); iS++) {
//            for ( int iL=0; iL<sprlr[iS].nLayers(); iL++ ) {
//                    boost::shared_ptr<mu2e::ITLayer> ilr = sprlr[iS].getLayer(iL);
//                    for ( int iC=0; iC<ilr->nCells(); iC++ ) {
//                            mu2e::Cell *s = ilr->getCell(iC).get();
//
//                            int idLayer =  s->Id().getLayer();
//                            //itwp->SelectCell(iS,idLayer,iC);
//                            unsigned long int index;
//                            if (itracker.isDumbbell()){
//                                    index = itwp->computeDet(iS,idLayer,iC,true);
//                            } else {
//                                    index = itwp->computeDet(iS,idLayer,iC);
//                            }
//                    }
//            }
//    }


    //-------------------------------------------------


    //char t;
    //cerr << "Double click in the Canvas " << _moduleLabel << " to continue:" ;
    //_canvas->cd(0);
    //_canvas->WaitPrimitive();
    //gPad->WaitPrimitive();
    //cerr << "Press enter to continue:" ;
    //cin>>t;
    //cerr << endl;

    cerr << "Double click in the canvas_Fake to continue:" ;
    _fakeCanvas->cd();
    TLatex *printEvN = new TLatex(0.15,0.4,Form("Current Event: %d",event.id().event()));
    printEvN->SetTextFont(62);
    printEvN->SetTextSizePixels(180);
    printEvN->Draw();
    _fakeCanvas->Update();
    _fakeCanvas->WaitPrimitive();
    cerr << endl;
    delete printEvN;

    //for (size_t i=0; i<nStrawPerEvent; ++i) { delete ((TEllipse *) hitDraws->At(i));  delete ((TEllipse *) hitDrawsAtZ->At(i)); }
    hitDraws->Delete();
    hitDrawsAtZ->Delete();
    delete hitDraws;
    delete hitDrawsAtZ;
    delete isStereoMin;

    for (std::vector<TH2F*>::iterator hCmapMin_it = rCMap_Min.begin(); hCmapMin_it != rCMap_Min.end(); ++hCmapMin_it ) {
            (*hCmapMin_it)->Delete();
    }
    for (std::vector<TH2F*>::iterator hCmapPl_it = rCMap_Pl.begin(); hCmapPl_it != rCMap_Pl.end(); ++hCmapPl_it ) {
            (*hCmapPl_it)->Delete();
    }

    delete tmpRCP_Pl;
    delete tmpRCP_Min;

//    if ( event.id().event() ){
//      _canvas->Modified();
//      _canvas->Update();
//
//      cerr << "Double click in the Canvas " << _moduleLabel << " to continue:" ;
//      _canvas->WaitPrimitive();
//      cerr << endl;
//
//    }


    //-------------------------------------------



  } // end analyze

  void ITDisplayData::endJob(){

    // cd() to correct root directory. See note 3.
    TDirectory* save = gDirectory;
    _directory->cd();

    // Write canvas.  See note 4.
    _canvas->Write();

    // cd() back to where we were.  See note 3.
    save->cd();

  }

}  // end namespace mu2e

using mu2e::ITDisplayData;
DEFINE_ART_MODULE(ITDisplayData);
