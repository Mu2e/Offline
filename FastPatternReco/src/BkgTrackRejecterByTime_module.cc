//
// Fast Patter recognition bck rejection algorithm based on time peak analysis
//
// $Id: BkgTrackRejecterByTime_module.cc,v 1.11 2012/05/24 13:44:28 ignatov Exp $
// $Author: ignatov $
// $Date: 2012/05/24 13:44:28 $
//
// Original author G. Tassielli
//

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "ITrackerGeom/inc/Cell.hh"
//#include "DataProducts/inc/DPIndexVectorCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
//#include "MCDataProducts/inc/GenParticleCollection.hh"
//#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
//#include "MCDataProducts/inc/StepPointMCCollection.hh"
//#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "ITrackerGeom/inc/ITracker.hh"
#include "TTrackerGeom/inc/TTracker.hh"
//#include "FastPatternReco/inc/TTHitPerTrackData.hh"
//#include "FastPatternReco/inc/GenTrackData.hh"
#include "RecoDataProducts/inc/TrackerHitTimeCluster.hh"
#include "RecoDataProducts/inc/TrackerHitTimeClusterCollection.hh"
#include "SeedService/inc/SeedService.hh"

// From art and its toolchain.
#include "art/Framework/Core/EDProducer.h"
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
//#include "TCanvas.h"
//#include "TDirectory.h"
#include "TMath.h"
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TClonesArray.h"
#include "TObjArray.h"
#include "TLine.h"
#include "TEllipse.h"
#include "TMarker.h"
#include "TBox.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TF1.h"
#include "TSpectrum.h"
#include "TSpectrumFit.h"
#include "TPolyMarker.h"
#include "TLatex.h"

// From CLHEP
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Random/RandFlat.h"

// C++ includes.
#include <iostream>
#include <string>
#include <memory>
#include <map>
#include <utility>
#include <limits>
#include <ctime>

//using namespace std;
//using art::Event;

namespace mu2e {

typedef art::Ptr<StrawHit> StrawHitPtr;
typedef std::multimap<unsigned int, StrawHitPtr, std::less<unsigned int> > stbrel;

  //class Straw;

  class BkgTrackRejecterByTime : public art::EDProducer {
  public:

    explicit BkgTrackRejecterByTime(fhicl::ParameterSet const& pset);
    virtual ~BkgTrackRejecterByTime() {
            //_peakFinder->Delete();
//            _fg->Delete();
//            if (_peakFinder!=0x0)        delete _peakFinder;
//            if (_fg!=0x0)                delete _fg;
//            if (_fakeCanvas!=0x0)        delete _fakeCanvas;
//            if (_peaksCanvases!=0x0)     delete _peaksCanvases;
//            if (_peaksCanvHistos!=0x0)   delete _peaksCanvHistos;
            if (_hPkStDistTrs!=0x0)      delete _hPkStDistTrs;
//            if (_hPkStDistanceTrs!=0x0)  delete _hPkStDistanceTrs;
            if (_hPkStDist2W!=0x0)       delete _hPkStDist2W;
            if (_hPkStClustDist2W!=0x0)  delete _hPkStClustDist2W;
            if (_hPkSecDist2W!=0x0)      delete _hPkSecDist2W;
            if (_hPkSecClustDist2W!=0x0) delete _hPkSecClustDist2W;
    }

    virtual void beginJob();
//    void endJob();

    // This is called for each event.
    void produce(art::Event & e);

  private:

    // Start: run time parameters

    // The module label of this module.
    std::string _moduleLabel;

    // Label of the G4 module
    std::string _g4ModuleLabel;

    // Name of the tracker StepPoint collection
    std::string _trackerStepPoints;

    // Label of the module that made the hits.
    std::string _makerModuleLabel;

    // Use the peak sigma or the time window (drift time + TOF) for the hit selection near the peak.
    float _nsigmaForTimeSel;
    bool  _useSigmaForTimeSel;

    // Use efficincy value to simulate the proton rejection by ADC value
    float _protonRejecEff;
    float _pulseCut;
    bool  _useProtonRejec;
    bool  _perfectRejec;
    int   _ADCSttnPeriod;
    int   _factorASP;
    bool  _useADCSttnPer;

    // Label of the generator.
    std::string _generatorModuleLabel;

    // drift vevocity
    double _driftVelocity;

//    // Number of events to accumulate between prompts.
//    int _nAccumulate;

    // End: run time parameters

    // Pointers to histograms, ntuples, TGraphs.
    TH1F*         _hMCHitDTime;
    TH1F*         _hHitTime;
    TH1F*         _hHitClustTime;
    TH1F*         _hHitTime_UpStrm;
    TH1F*         _hHitClustTime_UpStrm;
//    TCanvas*      _canvasPl;
//    TObjArray*    _peaksCanvases;
//    TObjArray*    _peaksCanvHistos;
    TObjArray*    _hPkStDistTrs;
    TObjArray*    _hPkStDist2W;
    TObjArray*    _hPkStClustDist2W;
    TObjArray*    _hPkSecDist2W;
    TObjArray*    _hPkSecClustDist2W;

    TH1I*         _hClockCicles;
    TH1F*         _hExecTime;

//    TF1*          _fg;
//    TCanvas*      _cnvForPeakstudy;

//    TCanvas*      _fakeCanvas;

//    TSpectrum *   _peakFinder;

    int   ntimeBin;
    float maxTimeHist; //ns
    float timeBinDim;  //ns

    // Random number distributions.
    //std::auto_ptr<CLHEP::RandFlat>    _randFlat;
    CLHEP::RandFlat    _randFlat;

    void smooth(TH1F *tmpHhitTime, TH1F *tmpHhitClustTime, float integrWindow );
    int  peakFinder(TH1F *tmpHhitTime, TH1F *tmpHhitClustTime, float &sigmaFitted, double *timepeakPos, double *timepeakHei);
    void extractHitForTimePeak (std::auto_ptr<TrackerHitTimeClusterCollection> &thcc, TH1F *tmpHhitTime, stbrel &timeBin_Straw_rel, int nfound, float integrWindow, float sigmaFitted, double *timepeakPos, double *timepeakHei);
    int  rndup(float n);

    // Some ugly but necessary ROOT related bookkeeping:

    // The job needs exactly one instance of TApplication.  See note 1.
//    std::auto_ptr<TApplication> _application;

    // Save directory from beginJob so that we can go there in endJob. See note 3.
//    TDirectory* _directory;

  };

  BkgTrackRejecterByTime::BkgTrackRejecterByTime(fhicl::ParameterSet const& pset) :

    // Run time parameters
    _moduleLabel(pset.get<std::string>("module_label")),/*@module_label*/
    _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel")),
    _trackerStepPoints(pset.get<std::string>("trackerStepPoints","tracker")),
    _makerModuleLabel(pset.get<std::string>("makerModuleLabel")),
    _nsigmaForTimeSel(pset.get<float>("nsigmaForTimeSel")),
    _protonRejecEff(pset.get<float>("protonRejecEff",-1.0)),
    _pulseCut(pset.get<float>("pulseCut",-1.0)),
    _ADCSttnPeriod(pset.get<int>("ADCSttnPeriod",0)),
    _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
   /*_nAccumulate(pset.get<int>("nAccumulate",20)),*/
    _driftVelocity(pset.get<double>("driftVelocity",0.05)),   // mm/ns

    // ROOT objects that are the main focus of this example.
    _hMCHitDTime(0),
    _hHitTime(0),
    _hHitClustTime(0),
    _hHitTime_UpStrm(0),
    _hHitClustTime_UpStrm(0),
//    _canvasPl(0),
//    _peaksCanvases(0),
//    _peaksCanvHistos(0),
    _hPkStDistTrs(0),
    _hPkStDist2W(0),
    _hPkStClustDist2W(0),
    _hPkSecDist2W(0),
    _hPkSecClustDist2W(0),
    _hClockCicles(0),
    _hExecTime(0),

//    _fg(0),
//    _cnvForPeakstudy(0),

//    _fakeCanvas(0),

//    _peakFinder(0),
    ntimeBin(0),
    maxTimeHist(2500.0),
    timeBinDim(10.0),
    _randFlat( createEngine( art::ServiceHandle<SeedService>()->getSeed() ) )//,

    // Some ugly but necessary ROOT related bookkeeping.
//    _application(0),
//    _directory(0)
  {
          if (_nsigmaForTimeSel>0.0) _useSigmaForTimeSel=true;
          else _useSigmaForTimeSel=false;

          if (_protonRejecEff>0.0 && _pulseCut>0.0) {
                  _useProtonRejec=true;
                  if (_protonRejecEff>=1.0) {
                          _protonRejecEff=1.0;
                          _perfectRejec=true;
                  }
                  else {
                          _perfectRejec=false;
                  }
                  if (_ADCSttnPeriod>0 && _ADCSttnPeriod<18) {
                          _useADCSttnPer=true;
                          _factorASP=2*_ADCSttnPeriod;
                  }
                  else _useADCSttnPer=false;
          }
          else {
                  _useProtonRejec=false;
                  _perfectRejec=false;
          }

          // Tell the framework what we make.
          produces<TrackerHitTimeClusterCollection>();

  }

  void BkgTrackRejecterByTime::beginJob(){

    std::cout<<"Bkg rejection jos!"<<std::endl;

    // Get access to the TFile service and save current directory for later use.
    art::ServiceHandle<art::TFileService> tfs;

    ntimeBin    = (int) maxTimeHist/timeBinDim;

    // Create a histogram.
    _hMCHitDTime       = tfs->make<TH1F>( "hMCHitDTime",   "Delta Time of the MC Step hits per Event", 200, 0.0, 50.0  );
    _hHitTime          = tfs->make<TH1F>( "hHitTime",      "Time of the Hits per Event", ntimeBin, 0.0, maxTimeHist  );
    _hHitClustTime     = tfs->make<TH1F>( "hHitClustTime", "Cluster of Time of the Hits per Event", ntimeBin, 0.0, maxTimeHist  );
    _hHitTime_UpStrm          = tfs->make<TH1F>( "hHitTime_UpStrm",      "Time of the Hits per Event", ntimeBin, 0.0, maxTimeHist  );
    _hHitClustTime_UpStrm     = tfs->make<TH1F>( "hHitClustTime_UpStrm", "Cluster of Time of the Hits per Event", ntimeBin, 0.0, maxTimeHist  );

    _hMCHitDTime      ->SetXTitle("ns");
    _hHitTime         ->SetXTitle("ns");
    _hHitClustTime    ->SetXTitle("ns");
    _hHitTime_UpStrm         ->SetXTitle("ns");
    _hHitClustTime_UpStrm    ->SetXTitle("ns");

    _hClockCicles       = tfs->make<TH1I>( "hClockCicles",   "N clock cicles needed to analyze one Event by BkgTrackRejecterByTime", 1000, 200000, 400000  );
    _hExecTime          = tfs->make<TH1F>( "hExecTime",   "Execution time to analyze one Event by BkgTrackRejecterByTime", 1000, 0.0, 10.0  );

//    // If needed, create the ROOT interactive environment. See note 1.
//    if ( !gApplication ){
//      int    tmp_argc(0);
//      char** tmp_argv(0);
//      _application = std::auto_ptr<TApplication>(new TApplication( "noapplication", &tmp_argc, tmp_argv ));
//    }
//
//    gStyle->SetPalette(1);
//    gROOT->SetStyle("Plain");

//    _peakFinder = new TSpectrum(20);

//    _fg = new TF1("fg","gaus");
//    _cnvForPeakstudy = tfs->make<TCanvas>("cnvForPeakstudy","Peaks studies container");
//    TString name;
//    TString title;
//    int window_size(860);
    // Create a canvas with a unique name.  See note 2.

//    _peaksCanvases     = new TObjArray();
//    _peaksCanvHistos   = new TObjArray();
    _hPkStDistTrs      = new TObjArray();
    _hPkStDist2W       = new TObjArray();
    _hPkStClustDist2W  = new TObjArray();
    _hPkSecDist2W      = new TObjArray();
    _hPkSecClustDist2W = new TObjArray();

//    name  = "canvasPl_"     + _moduleLabel;
//    title = "Canvas for Plots " + _moduleLabel;
//    _canvasPl = tfs->make<TCanvas>(name,title,window_size,window_size);
//    _canvasPl->Divide(2,2);
//
//    _fakeCanvas = new TCanvas("canvas_Fake","double click for next event",300,100);
//
//    _canvasPl->cd(1);
//    _hHitTime->SetStats(kTRUE);
//    _hHitTime->Draw();
//    _canvasPl->cd(2);
//    _hHitClustTime->SetStats(kFALSE);
//    _hHitClustTime->Draw();

    // See note 3.
//    _directory = gDirectory;

//    if (_useProtonRejec) {
//            // Get the engine associated with this module instance.
//            art::ServiceHandle<art::RandomNumberGenerator> rng;
//            _randFlat = std::auto_ptr<CLHEP::RandFlat>  ( new CLHEP::RandFlat( rng->getEngine() ) );
//    }

  }

  void BkgTrackRejecterByTime::produce(art::Event & event ) {

/*
    // Ask the event to give us a "handle" to the requested hits.
    art::Handle<StepPointMCCollection> hitsHandle;
    event.getByLabel(_g4ModuleLabel,_trackerStepPoints,hitsHandle);
    StepPointMCCollection const& hits = *hitsHandle;

    // Fill histogram with number of hits per event.
    _hHitTransverse->Fill(hits.size());

    // Periodically update the displayed histogram.

    if ( event.id().event()%_nAccumulate==0 ){
      _canvas->Modified();
      _canvas->Update();

      cerr << "Double click in the Canvas " << _moduleLabel << " to continue:" ;
      _canvas->WaitPrimitive();
      cerr << std::endl;

    }
*/

    //--------------------------------------------

    _hMCHitDTime->Reset();
    _hHitTime->Reset();
    _hHitClustTime->Reset();
    _hHitTime_UpStrm->Reset();
    _hHitClustTime_UpStrm->Reset();

//    _peakFinder->Clear();
//    _peaksCanvases->Clear();

    std::auto_ptr<TrackerHitTimeClusterCollection> thcc(new TrackerHitTimeClusterCollection);

    float intTimeWind = 0.0;
    bool isBumbbell = false, isITracker = false;
    CellGeometryHandle *itwp=0x0;

    art::ServiceHandle<mu2e::GeometryService> geom;
    const Tracker& tracker = getTrackerOrThrow();
    if(geom->hasElement<mu2e::TTracker>()) {
            const TTracker &ttr = static_cast<const TTracker&>( tracker );
            //const std::vector<Device> ttrdev = ttr.getDevices();
            intTimeWind = ttr.strawRadius();
    } else if(geom->hasElement<mu2e::ITracker>()) {
            const ITracker &itr = static_cast<const ITracker&>( tracker );
            itwp = itr.getCellGeometryHandle();
            itwp->SelectCell(0,0,0);
            intTimeWind += itwp->GetCellRad();
            itwp->SelectCell(itr.nSuperLayers()-1,0,0);
            intTimeWind += itwp->GetCellRad();
            intTimeWind *= 0.5;
            isITracker = true;
            isBumbbell = itr.isDumbbell();
    }

    intTimeWind = intTimeWind/_driftVelocity + 20.0;   //max drift time + max TOF of a signal electron
    stbrel timeBin_Straw_rel;
    stbrel timeBin_Straw_UpStrm_rel;
    int tmpiTimeBin;


    art::Handle<StrawHitCollection> pdataHandle;
    event.getByLabel(_makerModuleLabel,pdataHandle);
    StrawHitCollection const* hits = pdataHandle.product();

    //    // Get the persistent data about the StrawHitsMCTruth.
    //    art::Handle<StrawHitMCTruthCollection> truthHandle;
    //    event.getByLabel(_makerModuleLabel,truthHandle);
    //    StrawHitMCTruthCollection const* hits_truth = truthHandle.product();

    // Get the persistent data about pointers to StepPointMCs
//    art::Handle<PtrStepPointMCVectorCollection> mcptrHandle;
//    event.getByLabel(_makerModuleLabel,"StrawHitMCPtr",mcptrHandle);
//    PtrStepPointMCVectorCollection const* hits_mcptr = mcptrHandle.product();

//    // Get the persistent data about pointers to StepPointMCs
//    art::Handle<DPIndexVectorCollection> mcptrHandle;
//    event.getByLabel(_makerModuleLabel,"StrawHitMCPtr",mcptrHandle);
//    DPIndexVectorCollection const* hits_mcptr = mcptrHandle.product();

    // Get the persistent data about the StepPointMCs. More correct implementation
    // should look for product ids in DPIndexVectorCollection, rather than
    // use producer name directly ("g4run").

//    art::Handle<StepPointMCCollection> mchitsHandle;
//    event.getByLabel(_g4ModuleLabel,_trackerStepPoints,mchitsHandle);
//    StepPointMCCollection const* mchits = mchitsHandle.product();

//    if (!(hits->size() == hits_truth->size() &&
//          hits_mcptr->size() == hits->size() ) ) {
//      throw cet::exception("RANGE")
//        << "Strawhits: " << hits->size()
//        << " MCTruthStrawHits: " << hits_truth->size()
//        << " MCPtr: " << hits_mcptr->size();
//    }
//
//    // Get handles to the generated and simulated particles.
//    art::Handle<GenParticleCollection> genParticles;
//    event.getByLabel(_generatorModuleLabel, genParticles);

//    art::Handle<SimParticleCollection> simParticles;
//    event.getByLabel(_g4ModuleLabel, simParticles);


    clock_t startClock = clock();
    clock_t clocksForProtonRej, starProtRejec, stopProtRejec;
    clocksForProtonRej=0;

    size_t nStrawPerEvent = hits->size();

//    double mchittime=0.0;

//    bool overlapped = false;
//    bool isFirst = false;

//    bool isThereProton=false;
//    bool isElFromTarget=false;
    bool isStationADCEquipped=false;

    for (size_t i=0; i<nStrawPerEvent; ++i) {
      // Access data
      StrawHit        const&      hit(hits->at(i));
//      StrawHitMCTruth const&    truth(hits_truth->at(i));
//      DPIndexVector   const&    mcptr(hits_mcptr->at(i));

//      PtrStepPointMCVector const&    mcptr(hits_mcptr->at(i));

      starProtRejec = clock();
      if (_useProtonRejec) {
              //StrawIndex si = hit.strawIndex();
              //const Straw & str = tracker.getStraw(si);
              if ( !isITracker && _useADCSttnPer ) {
                      //Get hit straw
                      StrawIndex si = hit.strawIndex();
                      const Straw & str = tracker.getStraw(si);
                      isStationADCEquipped=(((int) str.id().getDevice()/2)%_factorASP)<_ADCSttnPeriod;
              }
              if ( isITracker || (!_useADCSttnPer || isStationADCEquipped) ) {
                      //std::cout<<"ADC on Station "<<((int) str.id().getDevice()/2)<<std::endl;
                      //isThereProton=false;
                      //isElFromTarget=false;
                      //for (size_t j = 0; j < mcptr.size(); ++j) {

                      //        //        StepPointMC const& mchit = (*mchits)[mcptr[j].index];
                      //        StepPointMC const& mchit = *mcptr[j];

                      //        // The simulated particle that made this hit.
                      //        SimParticleCollection::key_type trackId(mchit.trackId());
                      //        SimParticle const& sim = simParticles->at(trackId);
                      //        if ( sim.pdgId()==11 && sim.fromGenerator() ) isElFromTarget=true;
                      //        if ( sim.pdgId()==2212 ){
                      //                isThereProton=true;
                      //                //break;
                      //        }
                      //}
                      //if (isElFromTarget && isThereProton) continue;
                      //if (isThereProton) {
                      if (hit.energyDep()>_pulseCut) {
                              if (_perfectRejec) continue;
                              else if (_randFlat.fire()<_protonRejecEff) continue;
                      }
              }
      }
      stopProtRejec = clock();
      clocksForProtonRej+=(stopProtRejec-starProtRejec);

      //double hitEnergy = hit.energyDep();

      //Skip the straw if the energy of the hit is smaller than the minimum required
      //if (hitEnergy < _minimumEnergy) continue;

      //Get hit straw
//      StrawIndex si = hit.strawIndex();
//      const Straw & str = tracker.getStraw(si);

      // std::cout << "Getting informations about cells" << std::endl;

//      StrawId sid = str.id();
//      int stn     = sid.getStraw();
//      int layern  = sid.getLayer();
//      int devicen = sid.getDevice();
//      int sectorn = sid.getSector();

      //time of the hit
      double hitTime = hit.time();

      // std::cout << "Reading MCtruth info" << std::endl;
      if (isITracker && isBumbbell) {
              StrawIndex si = hit.strawIndex();
              //const Straw & str = tracker.getStraw(si);
              //const Cell & cell = static_cast<const Cell&>( str );
              itwp->SelectCellDet(si.asUint());
              if ( itwp->isUpStream() ) {
                      tmpiTimeBin=_hHitTime_UpStrm->Fill(hitTime);
                      if (tmpiTimeBin>0) timeBin_Straw_UpStrm_rel.insert( stbrel::value_type((unsigned int) tmpiTimeBin, StrawHitPtr( pdataHandle, i) ) );
              } else {
                      tmpiTimeBin=_hHitTime->Fill(hitTime);
                      if (tmpiTimeBin>0) timeBin_Straw_rel.insert( stbrel::value_type((unsigned int) tmpiTimeBin, StrawHitPtr( pdataHandle, i) ) );
              }
      } else {
              tmpiTimeBin=_hHitTime->Fill(hitTime);
              if (tmpiTimeBin>0) timeBin_Straw_rel.insert( stbrel::value_type((unsigned int) tmpiTimeBin, StrawHitPtr( pdataHandle, i) ) );
      }

      //common index for vectors
      //size_t trackIdx(0);


    }


//    _canvasPl->cd(1);
//    //_hHitTime->SetStats(kTRUE);
//    _hHitTime->Draw();
////    _hSelHitTime->Draw("same");

    float sigmaFitted=-1.0;
    double *timepeakPos = new double[50];
    double *timepeakHei = new double[50];
    int nfound = 0;
    if (_hHitTime->GetEntries()>0 ) {
            std::cout<<"Searching into the Downstream part"<<std::endl;
            smooth( _hHitTime, _hHitClustTime, intTimeWind );

            //    _canvasPl->cd(2);

            int nfound = peakFinder(_hHitTime,_hHitClustTime, sigmaFitted, timepeakPos, timepeakHei);

            //    _fg->SetParameter(1,_hSelHitClustTime->GetMean());
            //    _fg->SetParameter(2,_hSelHitClustTime->GetRMS());
            //    _hSelHitClustTime->Fit("fg","0");
            //    _hSelHitClustTimeSigma->Fill(_fg->GetParameter(2));


            //    for ( stbrel::const_iterator stb_it=timeBin_Straw_rel.begin(); stb_it!=timeBin_Straw_rel.end(); ++stb_it ) {
            //            std::cout<<"i-th time bin "<<stb_it->first<<" hit id "<<stb_it->second<<std::endl;
            //    }

            if (nfound>0){
                    extractHitForTimePeak (thcc, _hHitTime, timeBin_Straw_rel, nfound, intTimeWind, sigmaFitted, timepeakPos, timepeakHei);
            }
            std::cout<<"Downstream part found n time peak "<<nfound<<std::endl;
    }
    if (isBumbbell && _hHitTime_UpStrm->GetEntries()>0) {
            std::cout<<"Searching into the Upstream part"<<std::endl;
            smooth( _hHitTime_UpStrm, _hHitClustTime_UpStrm, intTimeWind );
            sigmaFitted=-1.0;
            nfound = peakFinder(_hHitTime_UpStrm,_hHitClustTime_UpStrm, sigmaFitted, timepeakPos, timepeakHei);
            if (nfound>0){
                    extractHitForTimePeak (thcc, _hHitTime_UpStrm, timeBin_Straw_UpStrm_rel, nfound, intTimeWind, sigmaFitted, timepeakPos, timepeakHei);
            }
            std::cout<<"Upstream part found n time peak "<<nfound<<std::endl;
    }

    event.put(thcc);

    delete [] timepeakPos;
    delete [] timepeakHei;

    clock_t stopClock = clock();
    _hClockCicles->Fill((unsigned long)(stopClock-startClock-clocksForProtonRej));
    _hExecTime->Fill( (float)(stopClock-startClock-clocksForProtonRej)/((float) CLOCKS_PER_SEC ) );
    std::cout<<"-------- N clok to analyze 1 ev by BkgTrackRejecterByTime "<<stopClock-startClock-clocksForProtonRej<<" @ "<<CLOCKS_PER_SEC<<std::endl;


//    _canvasPl->Modified();
//    _canvasPl->Update();
//
//
//    cerr << "Double click in the canvas_Fake to continue:" ;
//    _fakeCanvas->cd();
//    TLatex *printEvN = new TLatex(0.15,0.4,Form("Current Event: %d",event.id().event()));
//    printEvN->SetTextFont(62);
//    printEvN->SetTextSizePixels(180);
//    printEvN->Draw();
//    _fakeCanvas->Update();
//    _fakeCanvas->WaitPrimitive();
//    cerr << std::endl;
//    delete printEvN;

//    _peaksCanvases->Delete();
//    _peaksCanvHistos->Delete();


    //-------------------------------------------



  } // end produce

//  void BkgTrackRejecterByTime::endJob(){
//
//    // cd() to correct root directory. See note 3.
//    TDirectory* save = gDirectory;
//    _directory->cd();
//
//    // Write canvas.  See note 4.
////    _canvas->Write();
//
//    // cd() back to where we were.  See note 3.
//    save->cd();
//
//  }

  void BkgTrackRejecterByTime::smooth(TH1F *tmpHhitTime, TH1F *tmpHhitClustTime, float integrWindow ) {
          int nbingroup= (int) integrWindow/timeBinDim;
          float nhit=0.0;
          //float Selhit=0.0;
          int binHalf = (int) nbingroup/2;
          int jbin;
          for (int it=1; it<=ntimeBin; it++) {
            nhit=0.0;
            //Selhit=0.0;
            for (int is=0; is<nbingroup; is++) {
              jbin = it -binHalf +is;
              if (jbin<1) continue;
              if (jbin>ntimeBin) break;
              nhit+=tmpHhitTime->GetBinContent(jbin);
            }
            tmpHhitClustTime->SetBinContent(it,nhit);
          }

  }

  int BkgTrackRejecterByTime::peakFinder(TH1F *tmpHhitTime, TH1F *tmpHhitClustTime, float &sigmaFitted, double *timepeakPos, double *timepeakHei) {
          //---------------------------------------- peak finder -----------------------------

          //    float thr = 10.0;
          //    int nfound = 0;
          //    if (tmpHhitTime->GetEntries()>0){
          //            thr/=tmpHhitClustTime->GetMaximum();
          //            if (thr>0.0 && thr<1.0) nfound = _peakFinder->Search( tmpHhitClustTime,2.61,"",thr );
          //    }
          //
          //    tmpHhitClustTime->Draw();
          //    _hSelHitClustTime->Draw("same");

          float thr = 5.0;
          int nfound = 0, oldNfound = 0;
          int ipkfinderIter=1;
          float oldChiRes=0.0, ChiRes;
          ChiRes=10000000000.0;
          float newSigma, sigma=2.61;
          double sigmaOut=-1.0, errSigmaOut=-1.0;
          float /*sigmaFitted=-1.0,*/ oldSigmaFitted=-1.0;
          float * BaseSource = new float[ntimeBin];
          float * source = new float[ntimeBin];
          float * dest = new float[ntimeBin];
          size_t ntbBytes=ntimeBin*sizeof(float);
          //    bool noFittingError=true;

          TH1F *tmpPeakFit = new TH1F("tmpPeakFit", "", ntimeBin, 0.0, maxTimeHist);
          //double *timepeakPos = new double[50];
          //double *timepeakHei = new double[50];

          //    std::vector<float> validPeakPos;
          //    float tmpPeakPos=0.0;

          if (tmpHhitTime->GetEntries()>0){
                  //memcpy(BaseSource,&(tmpHhitClustTime->GetArray()[1]),ntbBytes);
                  for (int ibin = 0; ibin < ntimeBin; ibin++) BaseSource[ibin]=tmpHhitClustTime->GetBinContent(ibin + 1);

                  thr/=tmpHhitClustTime->GetMaximum();
                  if (thr>0.0 && thr<1.0) {

                          while (ipkfinderIter<21) {
                                  std::cout<<"Number of iteration for peak finder "<<ipkfinderIter<<std::endl;
                                  oldChiRes=ChiRes;

                                  memcpy(source,BaseSource,ntbBytes);

                                  //std::cout<<"----------------- source arr -----------------"<<std::endl;
                                  //for (int isb=0; isb<ntimeBin; isb++) std::cout<<"i bin "<<isb<<" "<<source[isb]<<std::endl;
                                  //std::cout<<"----------------------------------------------"<<std::endl;

                                  newSigma = 0.1*((float)ipkfinderIter)*sigma;
                                  //std::cout<<"newSigma "<<newSigma<<std::endl;

                                  //_peakFinder->Clear();
                                  //nfound = _peakFinder->Search( tmpHhitClustTime,newSigma,"",thr );
                                  //nfound = _peakFinder->SearchHighRes(source, dest, ntimeBin, sigma, thr, kFALSE, 10000, kFALSE, 0);
                                  TSpectrum peaksSearcher;
                                  nfound =  peaksSearcher.Search( tmpHhitClustTime,newSigma,"",thr );

                                  if (nfound>0) {

                                          //                            //check if the peaks are inside the time window
                                          //                            validPeakPos.clear();
                                          //                            for (int iPF = 0; iPF < nfound; iPF++) {
                                          //                                    tmpPeakPos=peaksSearcher.GetPositionX()[iPF];
                                          //                                    if ( tmpPeakPos>0.0 && tmpPeakPos<maxTimeHist ) validPeakPos.push_back(tmpPeakPos);
                                          //                            }
                                          //                            nfound = validPeakPos.size();
                                          //                            //end checking

                                          //std::cout<<"N peaks found "<<nfound<<std::endl;

                                          bool *FixPos = new bool[nfound];
                                          bool *FixAmp = new bool[nfound];
                                          //filling in the initial estimates of the input parameters
                                          float *PosX = new float[nfound];
                                          float *PosX1 = new float[nfound];
                                          float *PosY = new float[nfound];
                                          //PosX = _peakFinder->GetPositionX();
                                          //PosX = peaksSearcher.GetPositionX();
                                          for (int iPF = 0; iPF < nfound; iPF++) {
                                                  PosX[iPF]   = peaksSearcher.GetPositionX()[iPF];
                                                  //                                    PosX[iPF]   = validPeakPos[iPF];
                                                  FixPos[iPF] = false;
                                                  FixAmp[iPF] = false;
                                                  PosX1[iPF]  = PosX[iPF]/timeBinDim;
                                                  PosY[iPF]   = BaseSource[(int)(PosX1[iPF])];
                                          }

                                          //TSpectrumFit *pfit=new TSpectrumFit(nfound);
                                          TSpectrumFit pfit(nfound);
                                          pfit.Clear();
                                          //TSpectrumFit pfit(nfound);
                                          //std::cout<<"Fitting peaks"<<std::endl;
                                          pfit.SetFitParameters(0.0, ntimeBin-1, 100, 0.1, pfit.kFitOptimChiCounts, pfit.kFitAlphaHalving, pfit.kFitPower2, pfit.kFitTaylorOrderFirst);
                                          pfit.SetPeakParameters(sigma, false, PosX1, FixPos, PosY, FixAmp);
                                          //pfit.FitStiefel(source);
                                          pfit.FitAwmi(source);
                                          //std::cout<<"Peaks fitted"<<std::endl;

                                          //std::cout<<"----------------- source arr -----------------"<<std::endl;
                                          //for (int isb=0; isb<ntimeBin; isb++) std::cout<<"i bin "<<isb<<" "<<source[isb]<<std::endl;
                                          //std::cout<<"----------------------------------------------"<<std::endl;
                                          ChiRes=abs(1.0-pfit.GetChi());
                                          if ( ( (std::numeric_limits<float>::has_infinity && ChiRes == std::numeric_limits<float>::infinity()) || ChiRes!=ChiRes ) ) {
                                                  //noFittingError=false;
                                                  sigmaFitted=newSigma;
                                                  for (int iPF = 0; iPF < nfound; iPF++) timepeakPos[iPF] = PosX[iPF];
                                                  delete [] FixPos;
                                                  delete [] FixAmp;
                                                  delete [] PosX;
                                                  delete [] PosX1;
                                                  delete [] PosY;
                                                  break;

                                          }
                                          pfit.GetSigma(sigmaOut,errSigmaOut);
                                          //std::cout<<"Fitted sigma "<<sigmaOut<<" err "<<errSigmaOut<<std::endl;
                                          sigmaFitted=sqrt(pow(sigmaOut,2)+pow(errSigmaOut,2));
                                          //std::cout<<"OldChiRes "<<oldChiRes<<" ChiRes "<<ChiRes<<std::endl;
                                          if ( ChiRes>oldChiRes ) {
                                                  nfound=oldNfound;
                                                  sigmaFitted=oldSigmaFitted;
                                                  delete [] FixPos;
                                                  delete [] FixAmp;
                                                  delete [] PosX;
                                                  delete [] PosX1;
                                                  delete [] PosY;
                                                  break;
                                          }
                                          oldNfound=nfound;
                                          oldSigmaFitted=sigmaFitted;
                                          for (int itb = 0; itb < ntimeBin; itb++) tmpPeakFit->SetBinContent(itb + 1,source[itb]);
                                          memcpy(timepeakPos,pfit.GetPositions(),nfound*sizeof(double));
                                          memcpy(timepeakHei,pfit.GetAmplitudes(),nfound*sizeof(double));
                                          for (int iPF = 0; iPF < nfound; iPF++) timepeakPos[iPF]*=timeBinDim;

                                          //delete pfit;
                                          //pfit->Delete();
                                          delete [] FixPos;
                                          delete [] FixAmp;
                                          delete [] PosX;
                                          delete [] PosX1;
                                          delete [] PosY;
                                          //peaksearcher->Delete();
                                          //delete peaksearcher;
                                          //std::cout<<"---- Ending of Number of iteration for peak finder "<<ipkfinderIter<<std::endl;
                                  }
                                  ipkfinderIter++;
                          }
                  }
          }

          std::cout<<"End of peak finder "<<ipkfinderIter<<std::endl;

          delete [] BaseSource;
          delete [] source;
          delete [] dest;

          //    tmpHhitClustTime->Draw();
          //    if ( noFittingError ) {
          //            TPolyMarker * pm = new TPolyMarker(nfound, timepeakPos, timepeakHei);
          //            tmpHhitClustTime->GetListOfFunctions()->Add(pm);
          //            pm->SetMarkerStyle(23);
          //            pm->SetMarkerColor(kBlue);
          //            pm->SetMarkerSize(1);
          //            tmpHhitClustTime->Draw();
          //
          //            tmpPeakFit->SetLineColor(kBlue);
          //            tmpPeakFit->SetLineStyle(2);
          //            tmpPeakFit->SetLineWidth(1);
          //            tmpPeakFit->Draw("SAME L");
          //    }
          //    _hSelHitClustTime->Draw("same");

          //---------------------------------------- end peak finder -----------------------------

          return nfound;
  }

  void BkgTrackRejecterByTime::extractHitForTimePeak (std::auto_ptr<TrackerHitTimeClusterCollection> &thcc, TH1F *tmpHhitTime, stbrel &timeBin_Straw_rel, int nfound, float integrWindow, float sigmaFitted, double *timepeakPos, double *timepeakHei) {

          //Float_t *timepeakPos=_peakFinder->GetPositionX();
          //Float_t *timepeakHei=_peakFinder->GetPositionY();
          unsigned int *timepeakPosId = new unsigned int[nfound];
          unsigned int frstTimeBinInP, lastTimeBinInP;
          //size_t ihit;
          //std::map <size_t, unsigned int>::iterator stHitTrkMrk_rel_it;

          int i1peak;
          int binHalf = (int) (integrWindow/timeBinDim)/2;

          StrawId sid;
          //int stn, layern, devicen, sectorn;
          //unsigned int absSect, devStId;

          TH2I *tmpStClustDist    = new TH2I("tmpStClustDist","tmp Smoothing of Device vs Straw multiplicity Distribution",36,0,36,1000,-200,800);
          TH2I *tmpSecClustDist   = new TH2I("tmpSecClustDist","tmp Smoothing of Device vs Sector multiplicity Distribution",36,0,36,20,-4,16);

          for (int ipeak=0; ipeak<nfound; ipeak++){

                  thcc->push_back(TrackerHitTimeCluster());

                  i1peak=ipeak;
                  ++i1peak;
//                    _peaksCanvases->AddAtAndExpand(new TCanvas(Form("canvasFor_%d-th_peak",i1peak),Form("Data of the peak at %f ns",timepeakPos[ipeak]),860,860),ipeak);
//                    _peaksCanvHistos->AddAtAndExpand(new TCanvas(Form("canvasForHistos_%d-th_peak",i1peak),Form("Histograms of the peak at %f ns",timepeakPos[ipeak]),860,860),ipeak);
                  _hPkStDistTrs->AddAtAndExpand(new TH2I(Form("h%d-Pk_StDistTrs",i1peak),Form("Sector-Straw Distribution in Transverse plane for the %d-th peak",i1peak),20,-4,16,50,0,50),ipeak);
//                    _hPkStDistanceTrs->AddAtAndExpand(new TH2I(Form("h%d-Pk_StDistanceTrs",i1peak),Form("Sector-Straw distance Distribution in Transverse plane for the %d-th peak",i1peak),20,-4,16,50,0,50),ipeak);
                  _hPkStDist2W->AddAtAndExpand(new TH2I(Form("h%d-Pk_StDist2W",i1peak),Form("Device vs Straw multiplicity Distribution for the %d-th peak",i1peak),36,0,36,1000,-200,800),ipeak);
                  _hPkStClustDist2W->AddAtAndExpand(new TH2I(Form("h%d-Pk_StClustDist2W",i1peak),Form("Smoothing of Device vs Straw multiplicity Distribution for the %d-th peak",i1peak),36,0,36,1000,-200,800),ipeak);
                  _hPkSecDist2W->AddAtAndExpand(new TH2I(Form("h%d-Pk_SecDist2W",i1peak),Form("Device vs Sector multiplicity Distribution for the %d-th peak",i1peak),36,0,36,20,-4,16),ipeak);
                  _hPkSecClustDist2W->AddAtAndExpand(new TH2I(Form("h%d-Pk_SecClustDist2W",i1peak),Form("Smoothing of Device vs Sector multiplicity Distribution for the %d-th peak",i1peak),36,0,36,20,-4,16),ipeak);

//                    TCanvas * iCanv=((TCanvas *) _peaksCanvases->At(ipeak));
//                    iCanv->Divide(2,2);

                  timepeakPosId[ipeak]=(unsigned int) (((unsigned int)timepeakPos[ipeak])/timeBinDim)+1;

                  std::cout<<"Peak pos "<<timepeakPos[ipeak]<<" = "<<timepeakPosId[ipeak]<<" bin center "<<tmpHhitTime->GetBinCenter(timepeakPosId[ipeak])<<std::endl;
                  std::cout<<"peak height "<<timepeakHei[ipeak]<<std::endl;
                  std::cout<<" n hit found at paek pos "<<timeBin_Straw_rel.count(timepeakPosId[ipeak])<<std::endl;
                  std::cout<<"First hit at peak pos "<<timeBin_Straw_rel.find(timepeakPosId[ipeak])->second.key()<<std::endl;

                  if (_useSigmaForTimeSel){
                          unsigned int width=rndup(_nsigmaForTimeSel*sigmaFitted);
                          frstTimeBinInP = timepeakPosId[ipeak]>width ? (timepeakPosId[ipeak] - width) : 0;
                          lastTimeBinInP = timepeakPosId[ipeak] + width;
                 } else {
                          frstTimeBinInP = (timepeakPosId[ipeak] > (((unsigned int) binHalf)+1) ) ? (timepeakPosId[ipeak] - binHalf -1) : 0;
                          lastTimeBinInP = timepeakPosId[ipeak] + binHalf +1;
                  }
                  std::cout<<"Peak range limits: "<<frstTimeBinInP<<" "<<lastTimeBinInP<<std::endl;

//                    timeBin_Straw_rel.begin();
//                    stbrel::const_iterator stb_it=timeBin_Straw_rel.find(frstTimeBinInP);
//                    int findFTB=frstTimeBinInP;
//                    while (stb_it==timeBin_Straw_rel.end()) {
//                            stb_it=timeBin_Straw_rel.begin();
//                            ++findFTB;
//                            stb_it=timeBin_Straw_rel.find(findFTB);
//                    }
                  stbrel::const_iterator stb_it=timeBin_Straw_rel.lower_bound(frstTimeBinInP);
                  thcc->back()/*at(ipeak)*/._minHitTime=tmpHhitTime->GetBinLowEdge(stb_it->first);//findFTB);
                  thcc->back()/*at(ipeak)*/._meanTime=timepeakPos[ipeak];

                  int hitFoundInPeak=0;
                  TH2I *ihPkStDistTrs      = ((TH2I *) _hPkStDistTrs->At(ipeak));
                  ihPkStDistTrs->SetStats(kFALSE);
                  ihPkStDistTrs->SetXTitle("absSect");
                  ihPkStDistTrs->SetYTitle("Straw");
                  ihPkStDistTrs->SetTitleOffset(1.2,"y");
                  TH2I *ihPkStDist2W       = ((TH2I *) _hPkStDist2W->At(ipeak));
                  ihPkStDist2W->SetStats(kFALSE);
                  ihPkStDist2W->SetXTitle("device");
                  ihPkStDist2W->SetYTitle("absStraw");
                  ihPkStDist2W->SetTitleOffset(1.35,"y");
                  TH2I *ihPkStClustDist2W  = ((TH2I *) _hPkStClustDist2W->At(ipeak));
                  ihPkStClustDist2W->SetStats(kFALSE);
                  ihPkStClustDist2W->SetXTitle("device");
                  ihPkStClustDist2W->SetYTitle("absStraw");
                  ihPkStClustDist2W->SetTitleOffset(1.35,"y");
                  TH2I *ihPkSecDist2W      = ((TH2I *) _hPkSecDist2W->At(ipeak));
                  ihPkSecDist2W->SetStats(kFALSE);
                  ihPkSecDist2W->SetXTitle("device");
                  ihPkSecDist2W->SetYTitle("absSect");
                  ihPkSecDist2W->SetTitleOffset(1.2,"y");
                  TH2I *ihPkSecClustDist2W = ((TH2I *) _hPkSecClustDist2W->At(ipeak));
                  ihPkSecClustDist2W->SetStats(kFALSE);
                  ihPkSecClustDist2W->SetXTitle("device");
                  ihPkSecClustDist2W->SetYTitle("absSect");
                  ihPkSecClustDist2W->SetTitleOffset(1.2,"y");

                  tmpStClustDist->Reset();
                  tmpSecClustDist->Reset();

//                    TCanvas * iCanvHist=((TCanvas *) _peaksCanvHistos->At(ipeak));
//                    iCanvHist->Divide(2,2);

//                    while (stb_it!=timeBin_Straw_rel.end()){
                  while (stb_it!=timeBin_Straw_rel.upper_bound(lastTimeBinInP)){
                          //if (stb_it->first>lastTimeBinInP) break;

                          //ihit=stb_it->second;
                          ++hitFoundInPeak;

                          thcc->back()/*at(ipeak)*/._maxHitTime=tmpHhitTime->GetBinLowEdge(stb_it->first)+timeBinDim;
                          thcc->back()/*at(ipeak)*/._selectedTrackerHits.push_back( stb_it->second );//art::Ptr<StrawHit> ( hits->at(ihit), ihit ) );

//                            StrawHit        const&      hit(hits->at(ihit));
//                            //Get hit straw
//                            StrawIndex si = hit.strawIndex();
//                            const Straw & str = tracker.getStraw(si);
//
//                            // std::cout << "Getting informations about cells" << std::endl;
//
//                            sid = str.id();
//                            stn     = sid.getStraw();
//                            layern  = sid.getLayer();
//                            devicen = sid.getDevice();
//                            sectorn = sid.getSector();
//
//                            absSect = (devicen%2)+2*sectorn;
//                            devStId = absSect*50+stn;
//
//                            //ihPkStDistTrs->Fill(devStId);
//                            ihPkStDistTrs->Fill(absSect,stn);
//                            if (absSect>=8 && absSect<=11)  ihPkStDistTrs->Fill(absSect-12,stn);
//                            if (absSect>=0 && absSect<=3)   ihPkStDistTrs->Fill(absSect+12,stn);
//
//                            ihPkSecDist2W->Fill(devicen,absSect);
//                            if (absSect>=8 && absSect<=11)  ihPkSecDist2W->Fill(devicen,absSect-12);
//                            if (absSect>=0 && absSect<=3)   ihPkSecDist2W->Fill(devicen,absSect+12);
//
//                            ihPkStDist2W->Fill(devicen,absSect*50+stn);
//                            if (absSect>=8 && absSect<=11)  ihPkStDist2W->Fill(devicen,(absSect-12)*50+stn);
//                            if (absSect>=0 && absSect<=3)   ihPkStDist2W->Fill(devicen,(absSect+12)*50+stn);
////                            ihPkStDist2W->Fill(devicen,devStId);
///*
//                            int nYbingroup=250;//5*50
//                            float Yhit=0.0;
//                            int YbinHalf = (int) nYbingroup/2;
//                            int jYbin;
//                            int nXBin=ihPkStDist2W->GetNbinsX();
//                            int nYBin=ihPkStDist2W->GetNbinsY();
//                            for (int itX=1; itX<=nXBin; itX++) {
//                                    for (int itY=1; itY<=nYBin; itY++) {
//                                            Yhit=0.0;
//                                            for (int is=0; is<nYbingroup; is++) {
//                                                    jYbin = itY -YbinHalf +is;
//                                                    if (jYbin<1) continue;
//                                                    if (jYbin>nYBin) break;
//                                                    Yhit+=ihPkStDist2W->GetBinContent(itX,jYbin);
//                                            }
//                                            tmpStClustDist->SetBinContent(itX,itY,Yhit);
//                                    }
//                            }
//                            int nXbingroup=9;
//                            float Xhit=0.0;
//                            int XbinHalf = (int) nXbingroup/2;
//                            int jXbin;
//                            for (int itY=1; itY<=nYBin; itY++) {
//                                    for (int itX=1; itX<=nXBin; itX++) {
//                                            Xhit=0.0;
//                                            for (int is=0; is<nXbingroup; is++) {
//                                                    jXbin = itX -XbinHalf +is;
//                                                    if (jXbin<1) continue;
//                                                    if (jXbin>nXBin) break;
//                                                    Xhit+=tmpStClustDist->GetBinContent(jXbin,itY);
//                                            }
//                                            ihPkStClustDist2W->SetBinContent(itX,itY,Xhit);
//                                    }
//                            }
//*/
//
//                            int nYbingroup=5;//5*50
//                            float Yhit=0.0;
//                            int YbinHalf = (int) nYbingroup/2;
//                            int jYbin;
//                            int nXBin=ihPkStDist2W->GetNbinsX();
//                            int nYBin=ihPkStDist2W->GetNbinsY();
//                            for (int itX=1; itX<=nXBin; itX++) {
//                                    for (int itY=1; itY<=nYBin; itY++) {
//                                            Yhit=0.0;
//                                            for (int is=0; is<nYbingroup; is++) {
//                                                    jYbin = itY -YbinHalf +is;
//                                                    if (jYbin<1) continue;
//                                                    if (jYbin>nYBin) break;
//                                                    Yhit+=ihPkSecDist2W->GetBinContent(itX,jYbin);
//                                            }
//                                            tmpSecClustDist->SetBinContent(itX,itY,Yhit);
//                                    }
//                            }
//                            int nXbingroup=9;
//                            float Xhit=0.0;
//                            int XbinHalf = (int) nXbingroup/2;
//                            int jXbin;
//                            for (int itY=1; itY<=nYBin; itY++) {
//                                    for (int itX=1; itX<=nXBin; itX++) {
//                                            Xhit=0.0;
//                                            for (int is=0; is<nXbingroup; is++) {
//                                                    jXbin = itX -XbinHalf +is;
//                                                    if (jXbin<1) continue;
//                                                    if (jXbin>nXBin) break;
//                                                    Xhit+=tmpSecClustDist->GetBinContent(jXbin,itY);
//                                            }
//                                            ihPkSecClustDist2W->SetBinContent(itX,itY,Xhit);
//                                    }
//                            }
//
//
//
////                            iCanv->cd(3);
////                            ihPkStDistTrs->Draw("col z");
////
////                            iCanv->cd(4);
////                            ihPkSecDist2W->Draw("col z");
////
////                            iCanvHist->cd(2);
////                            ihPkStDist2W->Draw("col z");
////
//////                            iCanvHist->cd(3);
//////                            tmpStClustDist->Draw("col z");
////
////                            iCanvHist->cd(3);
////                            ihPkSecClustDist2W->Draw("col z");
////
////                            iCanvHist->cd(4);
////                            ihPkStClustDist2W->Draw("col z");

                          ++stb_it;
                 }
                 std::cout<<"Total Hit found in "<<ipeak<<"-th peak: "<<hitFoundInPeak<<std::endl;

//                   iCanv->cd(1);
//                   _hHitTransverse->Draw("same");
//
//                   iCanv->cd(2);
//                   _hHitLongit->Draw("same");


//                    iCanv->Modified();
//                    iCanv->Update();

          }

          if (tmpStClustDist!=0x0) delete tmpStClustDist;
          if (tmpSecClustDist!=0x0) delete tmpSecClustDist;

          delete timepeakPosId;

          _hPkStDistTrs->Delete();
      //    _hPkStDistanceTrs->Delete();
          _hPkStDist2W->Delete();
          _hPkStClustDist2W->Delete();
          _hPkSecDist2W->Delete();
          _hPkSecClustDist2W->Delete();

  }

  int BkgTrackRejecterByTime::rndup(float n)//round up a float type and show one decimal place
  {
          float t;
          t=n-floor(n);
          if (t>=0.5)
          {
                  n*=10.00000;//where n is the multi-decimal float
                  n=ceil(n);
                  n/=10.00000;
          }
          else
          {
                  n*=10.00000;//where n is the multi-decimal float
                  n=floor(n);
                  n/=10.00000;
          }
          return (int)n;
  }

}  // end namespace mu2e

using mu2e::BkgTrackRejecterByTime;
DEFINE_ART_MODULE(BkgTrackRejecterByTime);
