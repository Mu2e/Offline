//
// Fast Patter recognition bck rejection algorithm based on time peak analysis
//
// $Id: BkgTrackRejecterByTime_module.cc,v 1.3 2011/10/11 17:31:06 tassiell Exp $
// $Author: tassiell $
// $Date: 2011/10/11 17:31:06 $
//
// Original author G. Tassielli
//


// C++ includes.
#include <iostream>
#include <string>
#include <memory>
#include <map>
#include <utility>
#include <limits>
#include <ctime>

// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/Handle.h"
#include "art/Persistency/Provenance/Provenance.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

// From CLHEP
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Random/RandFlat.h"

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

using namespace std;
using art::Event;

namespace mu2e {

typedef art::Ptr<StrawHit> StrawHitPtr;
typedef std::multimap<unsigned int, StrawHitPtr, less<unsigned int> > stbrel;

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
    bool  _useProtonRejec;
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
    //auto_ptr<CLHEP::RandFlat>    _randFlat;
    CLHEP::RandFlat    _randFlat;

    int rndup(float n);

    // Some ugly but necessary ROOT related bookkeeping:

    // The job needs exactly one instance of TApplication.  See note 1.
//    auto_ptr<TApplication> _application;

    // Save directory from beginJob so that we can go there in endJob. See note 3.
//    TDirectory* _directory;

  };

  BkgTrackRejecterByTime::BkgTrackRejecterByTime(fhicl::ParameterSet const& pset) :

    // Run time parameters
    _moduleLabel(pset.get<string>("module_label")),/*@module_label*/
    _g4ModuleLabel(pset.get<string>("g4ModuleLabel")),
    _trackerStepPoints(pset.get<string>("trackerStepPoints","tracker")),
    _makerModuleLabel(pset.get<string>("makerModuleLabel")),
    _nsigmaForTimeSel(pset.get<float>("nsigmaForTimeSel")),
    _protonRejecEff(pset.get<float>("protonRejecEff",-1.0)),
    _ADCSttnPeriod(pset.get<int>("ADCSttnPeriod",0)),
    _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
   /*_nAccumulate(pset.get<int>("nAccumulate",20)),*/
    _driftVelocity(pset.get<double>("driftVelocity",0.05)),   // mm/ns

    // ROOT objects that are the main focus of this example.
    _hMCHitDTime(0),
    _hHitTime(0),
    _hHitClustTime(0),
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
    _randFlat( createEngine( get_seed_value(pset)) )//,

    // Some ugly but necessary ROOT related bookkeeping.
//    _application(0),
//    _directory(0)
  {
          if (_nsigmaForTimeSel>0.0) _useSigmaForTimeSel=true;
          else _useSigmaForTimeSel=false;

          if (_protonRejecEff>0.0) {
                  _useProtonRejec=true;
                  if (_protonRejecEff>1.0) _protonRejecEff=1.0;
                  if (_ADCSttnPeriod>0 && _ADCSttnPeriod<18) {
                          _useADCSttnPer=true;
                          _factorASP=2*_ADCSttnPeriod;
                  }
                  else _useADCSttnPer=false;
          }
          else _useProtonRejec=false;

          // Tell the framework what we make.
          produces<TrackerHitTimeClusterCollection>();

  }

  void BkgTrackRejecterByTime::beginJob(){

    cout<<"Bkg rejection jos!"<<endl;

    // Get access to the TFile service and save current directory for later use.
    art::ServiceHandle<art::TFileService> tfs;

    ntimeBin    = (int) maxTimeHist/timeBinDim;

    // Create a histogram.
    _hMCHitDTime       = tfs->make<TH1F>( "hMCHitDTime",   "Delta Time of the MC Step hits per Event", 200, 0.0, 50.0  );
    _hHitTime          = tfs->make<TH1F>( "hHitTime",      "Time of the Hits per Event", ntimeBin, 0.0, maxTimeHist  );
    _hHitClustTime     = tfs->make<TH1F>( "hHitClustTime", "Cluster of Time of the Hits per Event", ntimeBin, 0.0, maxTimeHist  );

    _hMCHitDTime      ->SetXTitle("ns");
    _hHitTime         ->SetXTitle("ns");
    _hHitClustTime    ->SetXTitle("ns");

    _hClockCicles       = tfs->make<TH1I>( "hClockCicles",   "N clock cicles needed to analyze one Event by BkgTrackRejecterByTime", 1000, 200000, 400000  );
    _hExecTime          = tfs->make<TH1F>( "hExecTime",   "Execution time to analyze one Event by BkgTrackRejecterByTime", 1000, 0.0, 10.0  );

//    // If needed, create the ROOT interactive environment. See note 1.
//    if ( !gApplication ){
//      int    tmp_argc(0);
//      char** tmp_argv(0);
//      _application = auto_ptr<TApplication>(new TApplication( "noapplication", &tmp_argc, tmp_argv ));
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
//            _randFlat = auto_ptr<CLHEP::RandFlat>  ( new CLHEP::RandFlat( rng->getEngine() ) );
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
      cerr << endl;

    }
*/

    //--------------------------------------------

    _hMCHitDTime->Reset();
    _hHitTime->Reset();
    _hHitClustTime->Reset();

//    _peakFinder->Clear();
//    _peaksCanvases->Clear();

    auto_ptr<TrackerHitTimeClusterCollection> thcc(new TrackerHitTimeClusterCollection);


    const Tracker& tracker = getTrackerOrThrow();
    const TTracker &ttr = static_cast<const TTracker&>( tracker );
    const std::vector<Device> ttrdev = ttr.getDevices();

    float intTimeWind = ttr.strawRadius()/_driftVelocity + 20.0;   //max drift time + max TOF of a signal electron
    stbrel timeBin_Straw_rel;
    int tmpiTimeBin;


    art::Handle<StrawHitCollection> pdataHandle;
    event.getByLabel(_makerModuleLabel,pdataHandle);
    StrawHitCollection const* hits = pdataHandle.product();

    //    // Get the persistent data about the StrawHitsMCTruth.
    //    art::Handle<StrawHitMCTruthCollection> truthHandle;
    //    event.getByLabel(_makerModuleLabel,truthHandle);
    //    StrawHitMCTruthCollection const* hits_truth = truthHandle.product();

    // Get the persistent data about pointers to StepPointMCs
    art::Handle<PtrStepPointMCVectorCollection> mcptrHandle;
    event.getByLabel(_makerModuleLabel,"StrawHitMCPtr",mcptrHandle);
    PtrStepPointMCVectorCollection const* hits_mcptr = mcptrHandle.product();

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

    art::Handle<SimParticleCollection> simParticles;
    event.getByLabel(_g4ModuleLabel, simParticles);


    clock_t startClock = clock();
    clock_t clocksForProtonRej, starProtRejec, stopProtRejec;
    clocksForProtonRej=0;

    size_t nStrawPerEvent = hits->size();

//    double mchittime=0.0;

//    bool overlapped = false;
//    bool isFirst = false;

    bool isThereProton=false;
    bool isElFromTarget=false;
    bool isStationADCEquipped=false;

    for (size_t i=0; i<nStrawPerEvent; ++i) {
      // Access data
      StrawHit        const&      hit(hits->at(i));
//      StrawHitMCTruth const&    truth(hits_truth->at(i));
//      DPIndexVector   const&    mcptr(hits_mcptr->at(i));

      PtrStepPointMCVector const&    mcptr(hits_mcptr->at(i));

      starProtRejec = clock();
      if (_useProtonRejec) {
              //StrawIndex si = hit.strawIndex();
              //const Straw & str = tracker.getStraw(si);
              if ( _useADCSttnPer ) {
                      //Get hit straw
                      StrawIndex si = hit.strawIndex();
                      const Straw & str = tracker.getStraw(si);
                      isStationADCEquipped=(((int) str.id().getDevice()/2)%_factorASP)<_ADCSttnPeriod;
              }
              if ( !_useADCSttnPer || isStationADCEquipped ) {
                      //cout<<"ADC on Station "<<((int) str.id().getDevice()/2)<<endl;
                      isThereProton=false;
                      isElFromTarget=false;
                      for (size_t j = 0; j < mcptr.size(); ++j) {

                              //        StepPointMC const& mchit = (*mchits)[mcptr[j].index];
                              StepPointMC const& mchit = *mcptr[j];

                              // The simulated particle that made this hit.
                              SimParticleCollection::key_type trackId(mchit.trackId());
                              SimParticle const& sim = simParticles->at(trackId);
                              if ( sim.pdgId()==11 && sim.fromGenerator() ) isElFromTarget=true;
                              if ( sim.pdgId()==2212 ){
                                      isThereProton=true;
                                      //break;
                              }
                      }
                      if (isElFromTarget && isThereProton) continue;
                      if (isThereProton) {
                              if (_randFlat.fire()<_protonRejecEff) continue;
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

      // cout << "Getting informations about cells" << endl;

//      StrawId sid = str.id();
//      int stn     = sid.getStraw();
//      int layern  = sid.getLayer();
//      int devicen = sid.getDevice();
//      int sectorn = sid.getSector();

      //time of the hit
      double hitTime = hit.time();

      // cout << "Reading MCtruth info" << endl;

      tmpiTimeBin=_hHitTime->Fill(hitTime);
      if (tmpiTimeBin>0) timeBin_Straw_rel.insert( stbrel::value_type((unsigned int) tmpiTimeBin, StrawHitPtr( pdataHandle, i) ) );

      //common index for vectors
      //size_t trackIdx(0);


    }


//    _canvasPl->cd(1);
//    //_hHitTime->SetStats(kTRUE);
//    _hHitTime->Draw();
////    _hSelHitTime->Draw("same");

    int nbingroup= (int) intTimeWind/timeBinDim;
    float hit=0.0;
    float Selhit=0.0;
    int binHalf = (int) nbingroup/2;
    int jbin;
    for (int it=1; it<=ntimeBin; it++) {
      hit=0.0;
      Selhit=0.0;
      for (int is=0; is<nbingroup; is++) {
        jbin = it -binHalf +is;
        if (jbin<1) continue;
        if (jbin>ntimeBin) break;
        hit+=_hHitTime->GetBinContent(jbin);
      }
      _hHitClustTime->SetBinContent(it,hit);
    }

//    _canvasPl->cd(2);

    //---------------------------------------- peak finder -----------------------------

//    float thr = 10.0;
//    int nfound = 0;
//    if (_hHitTime->GetEntries()>0){
//            thr/=_hHitClustTime->GetMaximum();
//            if (thr>0.0 && thr<1.0) nfound = _peakFinder->Search( _hHitClustTime,2.61,"",thr );
//    }
//
//    _hHitClustTime->Draw();
//    _hSelHitClustTime->Draw("same");

    float thr = 5.0;
    int nfound = 0, oldNfound = 0;
    int ipkfinderIter=1;
    float oldChiRes=0.0, ChiRes;
    ChiRes=10000000000.0;
    float newSigma, sigma=2.61;
    double sigmaOut=-1.0, errSigmaOut=-1.0;
    float sigmaFitted=-1.0, oldSigmaFitted=-1.0;
    float * BaseSource = new float[ntimeBin];
    float * source = new float[ntimeBin];
    float * dest = new float[ntimeBin];
    size_t ntbBytes=ntimeBin*sizeof(float);
    bool noFittingError=true;

    TH1F *tmpPeakFit = new TH1F("tmpPeakFit", "", ntimeBin, 0.0, maxTimeHist);
    double *timepeakPos = new double[50];
    double *timepeakHei = new double[50];

//    std::vector<float> validPeakPos;
//    float tmpPeakPos=0.0;

    if (_hHitTime->GetEntries()>0){
            //memcpy(BaseSource,&(_hHitClustTime->GetArray()[1]),ntbBytes);
            for (int ibin = 0; ibin < ntimeBin; ibin++) BaseSource[ibin]=_hHitClustTime->GetBinContent(ibin + 1);

            thr/=_hHitClustTime->GetMaximum();
            if (thr>0.0 && thr<1.0) {

                    while (ipkfinderIter<21) {
                            cout<<"Number of iteration for peak finder "<<ipkfinderIter<<endl;
                            oldChiRes=ChiRes;

                            memcpy(source,BaseSource,ntbBytes);

                            //cout<<"----------------- source arr -----------------"<<endl;
                            //for (int isb=0; isb<ntimeBin; isb++) cout<<"i bin "<<isb<<" "<<source[isb]<<endl;
                            //cout<<"----------------------------------------------"<<endl;

                            newSigma = 0.1*((float)ipkfinderIter)*sigma;
                            //cout<<"newSigma "<<newSigma<<endl;

                            //_peakFinder->Clear();
                            //nfound = _peakFinder->Search( _hHitClustTime,newSigma,"",thr );
                            //nfound = _peakFinder->SearchHighRes(source, dest, ntimeBin, sigma, thr, kFALSE, 10000, kFALSE, 0);
                            TSpectrum peaksSearcher;
                            nfound =  peaksSearcher.Search( _hHitClustTime,newSigma,"",thr );

                            if (nfound>0) {

                                    //                            //check if the peaks are inside the time window
                                    //                            validPeakPos.clear();
                                    //                            for (int iPF = 0; iPF < nfound; iPF++) {
                                    //                                    tmpPeakPos=peaksSearcher.GetPositionX()[iPF];
                                    //                                    if ( tmpPeakPos>0.0 && tmpPeakPos<maxTimeHist ) validPeakPos.push_back(tmpPeakPos);
                                    //                            }
                                    //                            nfound = validPeakPos.size();
                                    //                            //end checking

                                    //cout<<"N peaks found "<<nfound<<endl;

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
                                    //cout<<"Fitting peaks"<<endl;
                                    pfit.SetFitParameters(0.0, ntimeBin-1, 100, 0.1, pfit.kFitOptimChiCounts, pfit.kFitAlphaHalving, pfit.kFitPower2, pfit.kFitTaylorOrderFirst);
                                    pfit.SetPeakParameters(sigma, false, PosX1, FixPos, PosY, FixAmp);
                                    //pfit.FitStiefel(source);
                                    pfit.FitAwmi(source);
                                    //cout<<"Peaks fitted"<<endl;

                                    //cout<<"----------------- source arr -----------------"<<endl;
                                    //for (int isb=0; isb<ntimeBin; isb++) cout<<"i bin "<<isb<<" "<<source[isb]<<endl;
                                    //cout<<"----------------------------------------------"<<endl;
                                    ChiRes=abs(1.0-pfit.GetChi());
                                    if ( ( (std::numeric_limits<float>::has_infinity && ChiRes == std::numeric_limits<float>::infinity()) || ChiRes!=ChiRes ) ) {
                                            noFittingError=false;
                                            sigmaFitted=newSigma;
                                            for (int iPF = 0; iPF < nfound; iPF++) timepeakPos[iPF] = PosX[iPF];
                                            delete FixPos;
                                            delete FixAmp;
                                            delete PosX;
                                            delete PosX1;
                                            delete PosY;
                                            break;

                                    }
                                    pfit.GetSigma(sigmaOut,errSigmaOut);
                                    //cout<<"Fitted sigma "<<sigmaOut<<" err "<<errSigmaOut<<endl;
                                    sigmaFitted=sqrt(pow(sigmaOut,2)+pow(errSigmaOut,2));
                                    //cout<<"OldChiRes "<<oldChiRes<<" ChiRes "<<ChiRes<<endl;
                                    if ( ChiRes>oldChiRes ) {
                                            nfound=oldNfound;
                                            sigmaFitted=oldSigmaFitted;
                                            delete FixPos;
                                            delete FixAmp;
                                            delete PosX;
                                            delete PosX1;
                                            delete PosY;
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
                                    delete FixPos;
                                    delete FixAmp;
                                    delete PosX;
                                    delete PosX1;
                                    delete PosY;
                                    //peaksearcher->Delete();
                                    //delete peaksearcher;
                                    //cout<<"---- Ending of Number of iteration for peak finder "<<ipkfinderIter<<endl;
                            }
                            ipkfinderIter++;
                    }
            }
    }

    cout<<"End of peak finder "<<ipkfinderIter<<endl;

    delete BaseSource;
    delete source;
    delete dest;

//    _hHitClustTime->Draw();
//    if ( noFittingError ) {
//            TPolyMarker * pm = new TPolyMarker(nfound, timepeakPos, timepeakHei);
//            _hHitClustTime->GetListOfFunctions()->Add(pm);
//            pm->SetMarkerStyle(23);
//            pm->SetMarkerColor(kBlue);
//            pm->SetMarkerSize(1);
//            _hHitClustTime->Draw();
//
//            tmpPeakFit->SetLineColor(kBlue);
//            tmpPeakFit->SetLineStyle(2);
//            tmpPeakFit->SetLineWidth(1);
//            tmpPeakFit->Draw("SAME L");
//    }
//    _hSelHitClustTime->Draw("same");

    //---------------------------------------- end peak finder -----------------------------

//    _fg->SetParameter(1,_hSelHitClustTime->GetMean());
//    _fg->SetParameter(2,_hSelHitClustTime->GetRMS());
//    _hSelHitClustTime->Fit("fg","0");
//    _hSelHitClustTimeSigma->Fill(_fg->GetParameter(2));


//    for ( stbrel::const_iterator stb_it=timeBin_Straw_rel.begin(); stb_it!=timeBin_Straw_rel.end(); ++stb_it ) {
//            cout<<"i-th time bin "<<stb_it->first<<" hit id "<<stb_it->second<<endl;
//    }

    TH2I *tmpSecClustDist=0x0;
    TH2I *tmpStClustDist=0x0;
    if (nfound>0){
            //Float_t *timepeakPos=_peakFinder->GetPositionX();
            //Float_t *timepeakHei=_peakFinder->GetPositionY();
            unsigned int *timepeakPosId = new unsigned int[nfound];
            unsigned int frstTimeBinInP, lastTimeBinInP;
            //size_t ihit;
            std::map <size_t, unsigned int>::iterator stHitTrkMrk_rel_it;

            int i1peak;

            StrawId sid;
            //int stn, layern, devicen, sectorn;
            //unsigned int absSect, devStId;

            tmpStClustDist    = new TH2I("tmpStClustDist","tmp Smooting of Device vs Straw multiplicity Distribution",36,0,36,1000,-200,800);
            tmpSecClustDist   = new TH2I("tmpSecClustDist","tmp Smooting of Device vs Sector multiplicity Distribution",36,0,36,20,-4,16);

            for (int ipeak=0; ipeak<nfound; ipeak++){

            	    thcc->push_back(TrackerHitTimeCluster());

                    i1peak=ipeak;
                    ++i1peak;
//                    _peaksCanvases->AddAtAndExpand(new TCanvas(Form("canvasFor_%d-th_peak",i1peak),Form("Data of the peak at %f ns",timepeakPos[ipeak]),860,860),ipeak);
//                    _peaksCanvHistos->AddAtAndExpand(new TCanvas(Form("canvasForHistos_%d-th_peak",i1peak),Form("Histograms of the peak at %f ns",timepeakPos[ipeak]),860,860),ipeak);
                    _hPkStDistTrs->AddAtAndExpand(new TH2I(Form("h%d-Pk_StDistTrs",i1peak),Form("Sector-Straw Distribution in Transverse plane for the %d-th peak",i1peak),20,-4,16,50,0,50),ipeak);
//                    _hPkStDistanceTrs->AddAtAndExpand(new TH2I(Form("h%d-Pk_StDistanceTrs",i1peak),Form("Sector-Straw distance Distribution in Transverse plane for the %d-th peak",i1peak),20,-4,16,50,0,50),ipeak);
                    _hPkStDist2W->AddAtAndExpand(new TH2I(Form("h%d-Pk_StDist2W",i1peak),Form("Device vs Straw multiplicity Distribution for the %d-th peak",i1peak),36,0,36,1000,-200,800),ipeak);
                    _hPkStClustDist2W->AddAtAndExpand(new TH2I(Form("h%d-Pk_StClustDist2W",i1peak),Form("Smooting of Device vs Straw multiplicity Distribution for the %d-th peak",i1peak),36,0,36,1000,-200,800),ipeak);
                    _hPkSecDist2W->AddAtAndExpand(new TH2I(Form("h%d-Pk_SecDist2W",i1peak),Form("Device vs Sector multiplicity Distribution for the %d-th peak",i1peak),36,0,36,20,-4,16),ipeak);
                    _hPkSecClustDist2W->AddAtAndExpand(new TH2I(Form("h%d-Pk_SecClustDist2W",i1peak),Form("Smooting of Device vs Sector multiplicity Distribution for the %d-th peak",i1peak),36,0,36,20,-4,16),ipeak);

//                    TCanvas * iCanv=((TCanvas *) _peaksCanvases->At(ipeak));
//                    iCanv->Divide(2,2);

                    timepeakPosId[ipeak]=(unsigned int) (((unsigned int)timepeakPos[ipeak])/timeBinDim)+1;

                    cout<<"Peak pos "<<timepeakPos[ipeak]<<" = "<<timepeakPosId[ipeak]<<" bin center "<<_hHitTime->GetBinCenter(timepeakPosId[ipeak])<<endl;
                    cout<<"peak height "<<timepeakHei[ipeak]<<endl;
                    cout<<" n hit found at paek pos "<<timeBin_Straw_rel.count(timepeakPosId[ipeak])<<endl;
                    cout<<"First hit at peak pos "<<timeBin_Straw_rel.find(timepeakPosId[ipeak])->second.key()<<endl;

                    if (_useSigmaForTimeSel){
                            int width=rndup(_nsigmaForTimeSel*sigmaFitted);
                            frstTimeBinInP = timepeakPosId[ipeak] - width;
                            lastTimeBinInP = timepeakPosId[ipeak] + width;
                   } else {
                            frstTimeBinInP = timepeakPosId[ipeak] - binHalf -1;
                            lastTimeBinInP = timepeakPosId[ipeak] + binHalf +1;
                    }
                    cout<<"Peak range limits: "<<frstTimeBinInP<<" "<<lastTimeBinInP<<endl;

//                    timeBin_Straw_rel.begin();
//                    stbrel::const_iterator stb_it=timeBin_Straw_rel.find(frstTimeBinInP);
//                    int findFTB=frstTimeBinInP;
//                    while (stb_it==timeBin_Straw_rel.end()) {
//                            stb_it=timeBin_Straw_rel.begin();
//                            ++findFTB;
//                            stb_it=timeBin_Straw_rel.find(findFTB);
//                    }
                    stbrel::const_iterator stb_it=timeBin_Straw_rel.lower_bound(frstTimeBinInP);
                    thcc->at(ipeak)._minHitTime=_hHitTime->GetBinLowEdge(stb_it->first);//findFTB);
                    thcc->at(ipeak)._meanTime=timepeakPos[ipeak];

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

                            thcc->at(ipeak)._maxHitTime=_hHitTime->GetBinLowEdge(stb_it->first)+timeBinDim;
                            thcc->at(ipeak)._selectedTrackerHits.push_back( stb_it->second );//art::Ptr<StawHit> ( hits->at(ihit), ihit ) );

//                            StrawHit        const&      hit(hits->at(ihit));
//                            //Get hit straw
//                            StrawIndex si = hit.strawIndex();
//                            const Straw & str = tracker.getStraw(si);
//
//                            // cout << "Getting informations about cells" << endl;
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
                   cout<<"Total Hit found in "<<ipeak<<"-th peak: "<<hitFoundInPeak<<endl;

//                   iCanv->cd(1);
//                   _hHitTransverse->Draw("same");
//
//                   iCanv->cd(2);
//                   _hHitLongit->Draw("same");


//                    iCanv->Modified();
//                    iCanv->Update();

            }
            delete timepeakPosId;
    }

    event.put(thcc);

    clock_t stopClock = clock();
    _hClockCicles->Fill((unsigned long)(stopClock-startClock-clocksForProtonRej));
    _hExecTime->Fill( (float)(stopClock-startClock-clocksForProtonRej)/((float) CLOCKS_PER_SEC ) );
    cout<<"-------- N clok to analyze 1 ev by BkgTrackRejecterByTime "<<stopClock-startClock-clocksForProtonRej<<" @ "<<CLOCKS_PER_SEC<<endl;


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
//    _fakeCanvas->WaitPrimitive();
//    cerr << endl;
//    delete printEvN;

//    _peaksCanvases->Delete();
//    _peaksCanvHistos->Delete();
    _hPkStDistTrs->Delete();
//    _hPkStDistanceTrs->Delete();
    _hPkStDist2W->Delete();
    _hPkStClustDist2W->Delete();
    _hPkSecDist2W->Delete();
    _hPkSecClustDist2W->Delete();

    if (tmpStClustDist!=0x0) delete tmpStClustDist;
    if (tmpSecClustDist!=0x0) delete tmpSecClustDist;


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
