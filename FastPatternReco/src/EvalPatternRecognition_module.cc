//
// performance evaluation of the Pattern Recognition modules
//
// $Id: EvalPatternRecognition_module.cc,v 1.5 2013/10/21 21:01:22 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/10/21 21:01:22 $
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

#include <boost/shared_ptr.hpp>


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


#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Units/SystemOfUnits.h"

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "ITrackerGeom/inc/ITracker.hh"
#include "MCDataProducts/inc/VisibleGenElTrack.hh"
#include "MCDataProducts/inc/VisibleGenElTrackCollection.hh"
#include "RecoDataProducts/inc/TrackerHitTimeCluster.hh"
#include "RecoDataProducts/inc/TrackerHitTimeClusterCollection.hh"
#include "RecoDataProducts/inc/SectorStationCluster.hh"
#include "RecoDataProducts/inc/TrackSeed.hh"
#include "RecoDataProducts/inc/TrackSeedCollection.hh"
// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrationsI.hh"

// Root includes.
#include "TApplication.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TF1.h"
#include "TTree.h"

#define invSqrt12 0.2886751345948129

namespace mu2e {

  typedef art::Ptr<StrawHit> StrawHitPtr;

  class EvalPatternRecognition : public art::EDAnalyzer {
  public:

    explicit EvalPatternRecognition(fhicl::ParameterSet const& pset);
    virtual ~EvalPatternRecognition() {
//            if (_fakeCanvas)        delete _fakeCanvas;
    }

    virtual void beginJob();
    void endJob();

    // This is called for each event.
    void analyze(art::Event const& e);

  private:

    // Start: run time parameters

//    // The module label of this module.
//    std::string _moduleLabel;
//
//    // Label of the G4 module
//    std::string _g4ModuleLabel;
//
//    // Name of the tracker StepPoint collection
//    std::string _trackerStepPoints;

    // Label of the module that made the hits.
    std::string _makerModuleLabel;

//    // Label of the generator.
//    std::string _generatorModuleLabel;

    // Label of the module that made the hits.
    std::string _extractElectronsData;

    // Label of the module that made the hits.
    std::string _timeRejecterModuleLabel;
    float _nsigmaForTimeSel;

    // Label of the module that made the hits.
    std::string _patternRecoModuleLabel;

    // End: run time parameters

//    TCanvas*      _fakeCanvas;
    // Pointers to histograms, ntuples, TGraphs.
    bool          _evalTimeSel;
    TH1F*         _hNtrackableEv;
    TH1F*         _hNhitTrackableEv;
    TH1F*         _hNBestTimePeakEv;
    TH1F*         _hNhitBestTimePeakEv;
    TH2F*         _hLostHitBstTmPkEv;
    TH1F*         _hNoiseHitBstTmPkEv;
    TH1F*         _hT0resBstTmPkEv;
    TH1F*         _hT0pullBstTmPkEv;
    TH1F*         _hT0res1BstTmPkEv;
    TH1F*         _hT0pull1BstTmPkEv;
    TH1F*         _hT0resAveBstTmPkEv;
    TH1F*         _hT0pullAveBstTmPkEv;

    bool          _evalPRSel;
    TH1F*         _hNBestSeedInTmPkEv;
    TH1F*         _hSGSeedHitInBstTmPkEv;
    TH2F*         _hLostHitBstSGSeedBstTmPkEv;
    TH2F*         _hLostHitBstSGSeedEv;
    TH1F*         _hNoiseHitBstSGSeedEv;
    TH2F*         _hNLoopResBstSGSeedEv;

    TTree*        _dataEvalPR;
    unsigned int  runID, eventID, evNHit, convElNHit, convElNLoop;
    float         sel_pt, sel_omega, sel_omega_start, sel_omega_end, sel_tanDip, sel_tanDip_start, sel_tanDip_end, convElFHitTime;
    float         sel_d0, sel_phi0, sel_z0;
    int           bestTPcElNHit, bestTPNHit, nPeaksFound;
    float         bestTPcMinTime, bestTPcMeanTime, bestTPcSigmaTime, bestTPcNominalTimeWdt;
    int           bestSeed_idx, bestSeedcElNHit, bestSeedNHit, nPotentTracks;
    float         bestSeed_omega, bestSeed_tanDip, bestSeed_d0, bestSeed_phi0, bestSeed_z0, bestSeed_T0;

    // Some ugly but necessary ROOT related bookkeeping:

    // The job needs exactly one instance of TApplication.  See note 1.
    //unique_ptr<TApplication> _application;

    // Save directory from beginJob so that we can go there in endJob. See note 3.
    TDirectory* _directory;

  };

  EvalPatternRecognition::EvalPatternRecognition(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),

    // Run time parameters
    _makerModuleLabel(pset.get<std::string>("makerModuleLabel")),
    _extractElectronsData(pset.get<std::string>("elextractModuleLabel")),
    _timeRejecterModuleLabel(pset.get<std::string>("tRejecterModuleLabel","off")),
    _nsigmaForTimeSel(pset.get<float>("nsigmaForTimeSel")),
    _patternRecoModuleLabel(pset.get<std::string>("patternRecoModuleLabel","off")),
//    _moduleLabel(pset.get<std::string>("module_label")),/*@module_label*/
//    _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel")),
//    _trackerStepPoints(pset.get<std::string>("trackerStepPoints","tracker")),
//    _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),

//    _fakeCanvas(0),
    _hNtrackableEv(0),
    _hNhitTrackableEv(0),
    _hNBestTimePeakEv(0),
    _hNhitBestTimePeakEv(0),
    _hLostHitBstTmPkEv(0),
    _hNoiseHitBstTmPkEv(0),
    _hT0resBstTmPkEv(0),
    _hT0pullBstTmPkEv(0),
    _hT0res1BstTmPkEv(0),
    _hT0pull1BstTmPkEv(0),
    _hT0resAveBstTmPkEv(0),
    _hT0pullAveBstTmPkEv(0),

    _hNBestSeedInTmPkEv(0),
    _hSGSeedHitInBstTmPkEv(0),
    _hLostHitBstSGSeedBstTmPkEv(0),
    _hLostHitBstSGSeedEv(0),
    _hNoiseHitBstSGSeedEv(0),
    _hNLoopResBstSGSeedEv(0),
    _dataEvalPR(0),

    // Some ugly but necessary ROOT related bookkeeping.
    //_application(nullptr),
    _directory(0){
          runID=eventID=evNHit=convElNHit=convElNLoop=0;
          sel_omega=sel_omega_start=sel_omega_end=sel_tanDip=sel_tanDip_start=sel_tanDip_end=convElFHitTime=0.0;
          sel_d0=sel_phi0=sel_z0=0.0;
          bestTPcElNHit=bestTPNHit=nPeaksFound=0;
          bestTPcMinTime=bestTPcMeanTime=bestTPcSigmaTime=bestTPcNominalTimeWdt=0.0;
          bestSeed_idx=-1;
          bestSeedcElNHit=bestSeedNHit=nPotentTracks=0;
          bestSeed_omega=bestSeed_tanDip=bestSeed_d0=bestSeed_phi0=bestSeed_z0=bestSeed_T0=0.0;
          if (_timeRejecterModuleLabel.compare("off")==0) { _evalTimeSel=false; }
          else { _evalTimeSel=true; }
          if (_patternRecoModuleLabel.compare("off")==0) { _evalPRSel=false; }
          else { _evalPRSel=true; }
          std::cout<<"_evalTimeSel "<<_evalTimeSel<<" _evalPRSel "<<_evalPRSel<<std::endl;
 }

  void EvalPatternRecognition::beginJob(){

          std::cout<<"Starting EvalPatternRecognition jos!"<<std::endl;

          // Get access to the TFile service and save current directory for later use.
          art::ServiceHandle<art::TFileService> tfs;

          // If needed, create the ROOT interactive environment. See note 1.
          //    if ( !gApplication ){
          //      int    tmp_argc(0);
          //      char** tmp_argv(0);
          //      _application = unique_ptr<TApplication>(new TApplication( "noapplication", &tmp_argc, tmp_argv ));
          //    }

          gStyle->SetPalette(1);
          gROOT->SetStyle("Plain");

          // Create a histogram.
          if (_evalTimeSel) {
                  _hNtrackableEv       = tfs->make<TH1F>( "hNtrackableEv",   "N of trackable signal electrons", 111, 49.75, 105.25  );
                  _hNhitTrackableEv    = tfs->make<TH1F>( "hNhitTrackableEv",   "N of hits for each signal electrons track", 111, 49.75, 105.25  );
                  _hNBestTimePeakEv    = tfs->make<TH1F>( "hNBestTimePeakEv",   "N of peak time that has the best agreement with trackable signal electrons", 111, 49.75, 105.25  );
                  _hNhitBestTimePeakEv = tfs->make<TH1F>( "hNhitBestTimePeakEv",   "N of hits of the signal electrons track that are found in the best time peak", 111, 49.75, 105.25  );
                  _hLostHitBstTmPkEv   = tfs->make<TH2F>( "hLostHitBstTmPkEv",   "N of lost hits of the signal electrons track that in the best time peak", 111, 49.75, 105.25, 200, 0, 200  );
                  _hNoiseHitBstTmPkEv  = tfs->make<TH1F>( "hNoiseHitBstTmPkEv",   "N of hits of noise present in the the best time peak over signal electrons track hits", 111, 49.75, 105.25  );
                  _hT0resBstTmPkEv     = tfs->make<TH1F>( "hT0resBstTmPkEv",   "T0 resolution by using the peak time of the best time peak over signal electrons track hits", 400, -40, 40  );
                  _hT0pullBstTmPkEv    = tfs->make<TH1F>( "hT0pullBstTmPkEv",   "T0 pull by using the peak time of the best time peak over signal electrons track hits", 120, -6, 6  );
                  _hT0res1BstTmPkEv     = tfs->make<TH1F>( "hT0res1BstTmPkEv",   "T0 resolution by using a correction to the peak time of the best time peak over signal electrons track hits", 400, -40, 40  );
                  _hT0pull1BstTmPkEv    = tfs->make<TH1F>( "hT0pull1BstTmPkEv",   "T0 pull by using a correction to the peak time of the best time peak over signal electrons track hits", 120, -6, 6  );
                  _hT0resAveBstTmPkEv     = tfs->make<TH1F>( "hT0resAveBstTmPkEv",   "T0 resolution by averaging the two measure from the best time peak over signal electrons track hits", 400, -40, 40  );
                  _hT0pullAveBstTmPkEv    = tfs->make<TH1F>( "hT0pullAveBstTmPkEv",   "T0 pull by averaging the two measure from the best time peak over signal electrons track hits", 120, -6, 6  );

                  _hNtrackableEv      ->SetXTitle("pt [MeV]");
                  _hNhitTrackableEv   ->SetXTitle("pt [MeV]");
                  _hNBestTimePeakEv   ->SetXTitle("pt [MeV]");
                  _hNhitBestTimePeakEv->SetXTitle("pt [MeV]");
                  _hLostHitBstTmPkEv  ->SetXTitle("pt [MeV]");
                  _hNoiseHitBstTmPkEv ->SetXTitle("pt [MeV]");
                  _hT0resBstTmPkEv    ->SetXTitle("ns");
                  //_hT0pullBstTmPkEv
                  _hT0res1BstTmPkEv   ->SetXTitle("ns");
                  //_hT0pull1BstTmPkEv
                  _hT0resAveBstTmPkEv   ->SetXTitle("ns");
                  //_hT0pullAveBstTmPkEv
          }

          if(_evalPRSel) {
                  _hNBestSeedInTmPkEv  = tfs->make<TH1F>( "hNBestSeedInTmPkEv",   "N of Best Selected Seed in the best peak time", 111, 49.75, 105.25  );
                  _hSGSeedHitInBstTmPkEv = tfs->make<TH1F>( "hSGSeedHitInBstTmPkEv",   "N of signal electron hit selected by best seed from the best time peak", 111, 49.75, 105.25  );
                  _hLostHitBstSGSeedBstTmPkEv = tfs->make<TH2F>( "hLostHitBstSGSeedBstTmPkEv",   "N of lost hits of the signal electrons track that in the best seed from the best time peak", 111, 49.75, 105.25, 200, 0, 200  );
                  _hLostHitBstSGSeedEv   = tfs->make<TH2F>( "hLostHitBstSGSeedEv",   "N of lost hits of the signal electrons track that in the best seed", 111, 49.75, 105.25, 200, 0, 200  );
                  _hNoiseHitBstSGSeedEv  = tfs->make<TH1F>( "hNoiseHitBstSGSeedEv",   "N of hits of noise present in the the best seed over signal electrons track hits", 111, 49.75, 105.25  );
                  _hNLoopResBstSGSeedEv  = tfs->make<TH2F>( "hNLoopResBstSGSeedEv",   "Res of N of track loop measured by the best seed", 111, 49.75, 105.25, 21, -10.5, 10.5  );

                  _hNBestSeedInTmPkEv  ->SetXTitle("pt [MeV]");
                  _hSGSeedHitInBstTmPkEv ->SetXTitle("pt [MeV]");
                  _hLostHitBstSGSeedBstTmPkEv ->SetXTitle("pt [MeV]");
                  _hLostHitBstSGSeedEv ->SetXTitle("pt [MeV]");
                  _hNoiseHitBstSGSeedEv ->SetXTitle("pt [MeV]");
                  _hNLoopResBstSGSeedEv ->SetXTitle("pt [MeV]");
          }

          _dataEvalPR = tfs->make<TTree>( "dataEvalPR", "data for Pattern Recognition performance evaluation" );
          _dataEvalPR->Branch("Run",&runID,"runID/i");
          _dataEvalPR->Branch("Event",&eventID,"eventID/i");
          _dataEvalPR->Branch("EvTotNHit",&evNHit,"evNHit/i");
          _dataEvalPR->Branch("ConvElNHit",&convElNHit,"convElNHit/i");
          _dataEvalPR->Branch("ConvElNLoop",&convElNLoop,"convElNLoop/i");
          _dataEvalPR->Branch("ConvEl_omega_Tracker",&sel_omega,"sel_omega/F");
          _dataEvalPR->Branch("ConvEl_omega_start",&sel_omega_start,"sel_omega_start/F");
          _dataEvalPR->Branch("ConvEl_omega_end",&sel_omega_end,"sel_omega_end/F");
          _dataEvalPR->Branch("ConvEl_tanDip_Tracker",&sel_tanDip,"sel_tanDip/F");
          _dataEvalPR->Branch("ConvEl_tanDip_start",&sel_tanDip_start,"sel_tanDip_start/F");
          _dataEvalPR->Branch("ConvEl_tanDip_end",&sel_tanDip_end,"sel_tanDip_end/F");
          _dataEvalPR->Branch("ConvEl_d0_Tracker",&sel_d0,"sel_d0/F");
          _dataEvalPR->Branch("ConvEl_phi0_Tracker",&sel_phi0,"sel_phi0/F");
          _dataEvalPR->Branch("ConvEl_z0_Tracker",&sel_z0,"sel_z0/F");
          _dataEvalPR->Branch("ConvEl_FrstHit_Time",&convElFHitTime,"convElFHitTime/F");
          if (_evalTimeSel) {
                  _dataEvalPR->Branch("BestTPcElNHit",&bestTPcElNHit,"bestTPcElNHit/I");
                  _dataEvalPR->Branch("BestTPNHit",&bestTPNHit,"bestTPNHit/I");
                  _dataEvalPR->Branch("TotNPeaksFound",&nPeaksFound,"nPeaksFound/I");
                  _dataEvalPR->Branch("BestTPcMinTime",&bestTPcMinTime,"bestTPcMinTime/F");
                  _dataEvalPR->Branch("BestTPcMeanTime",&bestTPcMeanTime,"bestTPcMeanTime/F");
                  _dataEvalPR->Branch("BestTPcSigmaTime",&bestTPcSigmaTime,"bestTPcSigmaTime/F");
                  _dataEvalPR->Branch("BestTPcNominalTimeWdt",&bestTPcNominalTimeWdt,"bestTPcNominalTimeWdt/F");
          }
          if (_evalPRSel) {
                  _dataEvalPR->Branch("BestSeed_Idxt",&bestSeed_idx,"bestSeed_idx/I");
                  _dataEvalPR->Branch("BestSeedcElNHit",&bestSeedcElNHit,"bestSeedcElNHit/I");
                  _dataEvalPR->Branch("BestSeedNHit",&bestSeedNHit,"bestSeedNHit/I");
                  _dataEvalPR->Branch("NPotentTracksFound",&nPotentTracks,"nPotentTracks/I");
                  _dataEvalPR->Branch("BestSeed_omega",&bestSeed_omega,"bestSeed_omega/F");
                  _dataEvalPR->Branch("BestSeed_tanDip",&bestSeed_tanDip,"bestSeed_tanDip/F");
                  _dataEvalPR->Branch("BestSeed_d0",&bestSeed_d0,"bestSeed_d0/F");
                  _dataEvalPR->Branch("BestSeed_phi0",&bestSeed_phi0,"bestSeed_phi0/F");
                  _dataEvalPR->Branch("BestSeed_z0",&bestSeed_z0,"bestSeed_z0/F");
                  _dataEvalPR->Branch("BestSeed_T0",&bestSeed_T0,"bestSeed_T0/F");
          }
          //_fakeCanvas = new TCanvas("canvas_Fake","double click for next event",300,100);

          // See note 3.
          _directory = gDirectory;

  }

  void EvalPatternRecognition::analyze(art::Event const& event ) {


    //const Tracker& tracker = getTrackerOrThrow();
    //const TTracker &ttr = static_cast<const TTracker&>( tracker );
    //const std::vector<Device> ttrdev = ttr.getDevices();

    art::Handle<StrawHitCollection> pdataHandle;
    event.getByLabel(_makerModuleLabel,pdataHandle);
    StrawHitCollection const* hits = pdataHandle.product();


    art::Handle<VisibleGenElTrackCollection> genEltrksHandle;
    event.getByLabel(_extractElectronsData,genEltrksHandle);
    VisibleGenElTrackCollection const* genEltrks = genEltrksHandle.product();

    TrackerHitTimeClusterCollection * tclusts=0x0;
    if (_evalTimeSel) {
            art::Handle<TrackerHitTimeClusterCollection> tclustHandle;
            event.getByLabel(_timeRejecterModuleLabel,tclustHandle);
            tclusts = const_cast<TrackerHitTimeClusterCollection *> (tclustHandle.product());
    }

    TrackSeedCollection * tracksSeed=0x0;
    if (_evalPRSel) {
            art::Handle<TrackSeedCollection> tracksSeedHandle;
            event.getByLabel(_patternRecoModuleLabel,tracksSeedHandle);
            tracksSeed = const_cast<TrackSeedCollection *> (tracksSeedHandle.product());
    }

    std::vector<mu2e::VisibleGenElTrack>::const_iterator genEltrk_it;

    double intTimeWind(0.0);
    ConditionsHandle<TrackerCalibrations> tcal("ignored");
    D2T d2t;
    static CLHEP::Hep3Vector zdir(0.0,0.0,1.0);
    StrawIndex fakeSInd(0);
    CellGeometryHandle *itwp=0x0;

    art::ServiceHandle<mu2e::GeometryService> geom;
    const Tracker& tracker = getTrackerOrThrow();
    if(geom->hasElement<mu2e::TTracker>()) {
            const TTracker &ttr = static_cast<const TTracker&>( tracker );
            tcal->DistanceToTime(fakeSInd,0.5*ttr.strawRadius(),zdir,d2t);
            intTimeWind = d2t._tdrift;
    } else if(geom->hasElement<mu2e::ITracker>()) {
            const ITracker &itr = static_cast<const ITracker&>( tracker );
            itwp = itr.getCellGeometryHandle();
            itwp->SelectCell(0,0,0);
            tcal->DistanceToTime(fakeSInd,0.5*itwp->GetCellInsideRad(),zdir,d2t);
            intTimeWind += d2t._tdrift;
            itwp->SelectCell(itr.nSuperLayers()-1,0,0);
            tcal->DistanceToTime(fakeSInd,0.5*itwp->GetCellInsideRad(),zdir,d2t);
            intTimeWind += d2t._tdrift;
            intTimeWind *= 0.5;
    }
    intTimeWind = intTimeWind+20;   //to invert the time integration window that was the average max drift time + max TOF of a signal electron
    intTimeWind*=0.5;
    //intTimeWind+= 20.0;

    std::cout<<"--------------------------- Analysis ---------------------------"<<std::endl;
    std::cout<<"Event: "<<event.id().event()<<std::endl;
    std::cout<<"----------------------------------------------------------------"<<std::endl;

    StrawId sid;
    bool goodEvent=false;
    //VisibleGenElTrack *signaLEltrk;
    std::vector<StrawHitPtr> signElGoodHit;

    double rho;
    //double pt, omega, phi0, z0, d0, tanDip;

    double pt_start, rho_start;
    double pt_end, rho_end;
    //double sel_rho;
    double B=1.0;
    CLHEP::Hep2Vector radDir, centerDir, PCA;
    HepGeom::Point3D<double> CirCenter;

    runID=event.run();
    eventID=event.event();
    evNHit=hits->size();
    bestSeed_idx=bestTPcElNHit=bestTPNHit=nPeaksFound=-1;
    bestSeedcElNHit=bestSeedNHit=nPotentTracks=-1;

    for ( genEltrk_it = genEltrks->begin(); genEltrk_it!= genEltrks->end(); ++genEltrk_it ){
            VisibleGenElTrack &iEltrk = const_cast<VisibleGenElTrack &>(*genEltrk_it);
            unsigned short &nloops = iEltrk.getNumOfLoops();
            std::cout<<"N loops "<<nloops<<std::endl;
            /*for ( unsigned int ilp=0; ilp<nloops; ilp++ ){
                    GenElHitData& hdil = iEltrk.getithLoopHit(ilp);
                    pt = sqrt( pow(hdil._hitMomentum[0],2) + pow(hdil._hitMomentum[1],2) );
                    rho   = pt/(B*0.3);
                    omega = 1.0/rho;
                    //iEltrk.
                    std::cout<<ilp<<" -th loop: p_t "<<pt<<" rho mm "<<rho<<std::endl;
                    CirCenter.set(hdil._hitPoint.getX(),hdil._hitPoint.getY(),hdil._hitPoint.getZ());
                    radDir.setX(hdil._hitMomentum.getX());
                    radDir.setY(hdil._hitMomentum.getY());
                    radDir.rotate( ( (hdil._hitMomentum.getZ()>=0.0) ? 90.0 : -90.0 )*CLHEP::degree ); //only valid for normal direction of CE
                    radDir=radDir.unit();
                    CirCenter=CirCenter+rho*radDir;
                    std::cout<<" Hit Pos "<<hdil._hitPoint<<" Circ center "<<CirCenter<<std::endl;
                    centerDir.set(CirCenter.x(),CirCenter.y());
                    centerDir = centerDir.unit();
                    d0 = CirCenter.rho() - rho;
                    PCA = d0*centerDir;
                    phi0 = PCA.phi();
                    centerDir.rotate(180.0*CLHEP::degree);
                    radDir.rotate(180.0*CLHEP::degree);
                    double phi = centerDir.angle(radDir); //wrong if it will be > pi FIXME
                    tanDip = hdil._hitMomentum.getZ()/pt;
                    z0 = hdil._hitPoint.getZ() - rho*phi*tanDip;

            }*/

            if ( iEltrk.isConversionEl() ) {
                    if ( iEltrk.getNumOfHit()>=6 ) {
                            GenElHitData& hdil = iEltrk.getFirstHit();//getithLoopHit(0);
                            convElFHitTime = hdil._mcHitTime;
                            if (convElFHitTime>2400.0) goodEvent=false;
                            else {
                                    goodEvent=true;
                                    sel_pt = sqrt( pow(hdil._hitMomentum[0],2) + pow(hdil._hitMomentum[1],2) );
                                    rho = sel_pt/(B*0.3);
                                    sel_omega = 1.0/rho;
                                    sel_tanDip = hdil._hitMomentum[2]/sel_pt;

                                    CirCenter.set(hdil._hitPoint.getX(),hdil._hitPoint.getY(),hdil._hitPoint.getZ());
                                    radDir.setX(hdil._hitMomentum.getX());
                                    radDir.setY(hdil._hitMomentum.getY());
                                    radDir.rotate( ( (hdil._hitMomentum.getZ()>=0.0) ? 90.0 : -90.0 )*CLHEP::degree ); //only valid for normal direction of CE
                                    radDir=radDir.unit();
                                    CirCenter=CirCenter+rho*radDir;
                                    std::cout<<" Hit Pos "<<hdil._hitPoint<<" Circ center "<<CirCenter<<std::endl;
                                    centerDir.set(CirCenter.x(),CirCenter.y());
                                    centerDir = centerDir.unit();
                                    sel_d0 = CirCenter.rho() - rho;
                                    PCA = centerDir;
                                    PCA*=sel_d0;
                                    sel_phi0 = PCA.phi();
                                    centerDir.rotate(180.0*CLHEP::degree);
                                    radDir.rotate(180.0*CLHEP::degree);
                                    double phi = centerDir.angle(radDir); //wrong if it will be > pi FIXME
                                    std::cout<<"phi first hit "<<phi<<std::endl;
                                    sel_z0 = hdil._hitPoint.getZ();// - rho*phi*sel_tanDip;

                                    _hNtrackableEv->Fill(sel_pt);
                                    for ( unsigned int iElHit=0; iElHit<iEltrk.getNumOfHit(); iElHit++) {
                                            GenElHitData& genElhit = iEltrk.getHit((int) iElHit);
                                            if (genElhit._isFirst) signElGoodHit.push_back(genElhit._iHit);
                                    }
                                    convElNHit = signElGoodHit.size();
                                    convElNLoop = nloops;
                                    pt_start = sqrt( pow(iEltrk.getTrkLrntzVec().px(),2) + pow(iEltrk.getTrkLrntzVec().py(),2) );
                                    rho_start = pt_start/(B*0.3); //the magnetic field at production point is not 1 FIXME
                                    sel_omega_start = 1.0/rho_start;
                                    sel_tanDip_start = iEltrk.getTrkLrntzVec().pz()/pt_start;
                                    GenElHitData& lasth = iEltrk.getHit((int)(iEltrk.getNumOfHit()-1));
                                    pt_end = sqrt( pow(lasth._hitMomentum[0],2) + pow(lasth._hitMomentum[1],2) );
                                    rho_end = pt_end/(B*0.3);
                                    sel_omega_end = 1.0/ rho_end;
                                    sel_tanDip_end = lasth._hitMomentum[2]/pt_end;
                                    _hNhitTrackableEv->Fill(sel_pt,signElGoodHit.size());
                                    //signaLEltrk=&iEltrk;
                            }

                    }
                    //std::cout<<"CE track data:"<<std::endl;
                    //std::cout<<iEltrk;
            }
            //std::cout<<"All track data:"<<std::endl;
            //std::cout<<iEltrk;

    }


    unsigned int elHitFound;
    std::multimap< unsigned int, size_t, std::less<unsigned int> > foundElHitInPeak;
    std::multimap< unsigned int, size_t, std::less<unsigned int> >::reverse_iterator BestFoundElInPeak_it;
    std::map<size_t, std::vector<StrawHitPtr> > signElGoodHitsInTmpks;
    //std::map<size_t, double > signElFHitTimeInTmpks;

    //---------------- for time algorithm ----------------
    if (_evalTimeSel) {
            size_t nTimeClusPerEvent = tclusts->size();
            nPeaksFound=nTimeClusPerEvent;

            if (goodEvent) {


                    for (size_t ipeak=0; ipeak<nTimeClusPerEvent; ipeak++) {
                            TrackerHitTimeCluster const&  tclust(tclusts->at(ipeak));
                            std::cout<<"\t Hit in peak "<<tclust._selectedTrackerHits.size()<<std::endl;

                            elHitFound=0;
                            std::vector<StrawHitPtr> tmpElHits;

                            for (std::vector<StrawHitPtr>::const_iterator iTCHit=tclust._selectedTrackerHits.begin(); iTCHit!=tclust._selectedTrackerHits.end(); ++iTCHit){
                                    // Access data

                                    //for ( int iElHit=0; iElHit<signaLEltrk->getNumOfHit(); iElHit++) {
                                    //        GenElHitData& genEl = signaLEltrk->getHit(iElHit);
                                    for ( std::vector<StrawHitPtr>::iterator iGoodSelHit_it=signElGoodHit.begin(); iGoodSelHit_it!=signElGoodHit.end(); ++iGoodSelHit_it ) {
                                            if ( iTCHit->key()==iGoodSelHit_it->key() ) {
                                                    elHitFound++;
                                                    tmpElHits.push_back(*iTCHit);
                                                    break;
                                            }
                                    }

                                    //                //                std::cout<<"\t\t iHit in peak at "<<*iTCHit<<std::endl;
                                    //                //StrawHit        const&      hit(hits->at(*iTCHit));
                                    //                StrawHit        const&      hit=*(*iTCHit);
                                    //                StrawIndex si = hit.strawIndex();
                                    //                const Straw & str = tracker.getStraw(si);
                                    //
                                    //                // std::cout << "Getting informations about cells" <<std::endl;
                                    //
                                    //                sid = str.id();
                            }

                            if ( elHitFound>0 ) {
                                    foundElHitInPeak.insert( std::pair<unsigned int, size_t> (elHitFound,ipeak) );
                                    signElGoodHitsInTmpks.insert( std::pair< size_t, std::vector<StrawHitPtr> > (ipeak,tmpElHits) );
                                    //signElFHitTimeInTmpks.insert( std::pair< size_t, double > (ipeak,) );
                            }
                    }

                    if ( !foundElHitInPeak.empty() ) {
                            _hNBestTimePeakEv->Fill(sel_pt);
                            BestFoundElInPeak_it=foundElHitInPeak.rbegin();
                            bestTPcElNHit=BestFoundElInPeak_it->first;
                            TrackerHitTimeCluster const&  tclust(tclusts->at(BestFoundElInPeak_it->second));
                            bestTPNHit=tclust._selectedTrackerHits.size();
                            _hNhitBestTimePeakEv->Fill(sel_pt,bestTPcElNHit);
                            _hLostHitBstTmPkEv->Fill(sel_pt, (signElGoodHit.size() - bestTPcElNHit) );
                            _hNoiseHitBstTmPkEv->Fill(sel_pt, (bestTPNHit - bestTPcElNHit) );
                            double t0peak, errt0peak, t0res;
                            tclust.expectedT0(t0peak,errt0peak);
                            t0res =  t0peak - convElFHitTime;
                            _hT0resBstTmPkEv->Fill(t0res);
                            _hT0pullBstTmPkEv->Fill(t0res/errt0peak);

                            tclust.expectedT0(t0peak,errt0peak,1);
                            t0res =  t0peak - convElFHitTime;
                            _hT0resAveBstTmPkEv->Fill(t0res);
                            _hT0pullAveBstTmPkEv->Fill(t0res/errt0peak); //propagated sigma between sigma and 0.6*sigma

                            tclust.expectedT0(t0peak,errt0peak,2);
                            t0res =  t0peak - convElFHitTime;
                            _hT0res1BstTmPkEv->Fill(t0res);
                            _hT0pull1BstTmPkEv->Fill(t0res/errt0peak); //the variance of half gaussian is approximately equal to 0.6*sigma

                            bestTPcMinTime=tclust._minHitTime;
                            bestTPcMeanTime=tclust._meanTime;
                            bestTPcSigmaTime=tclust._sigma;
                            bestTPcNominalTimeWdt=tclust._nominalWidth;
                            //if (fabs(t0res)>10) std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! worng t peak res "<<t0res<<std::endl;
                    }
            }
    }

    //---------------- for PR algorithm ----------------
    if (_evalPRSel) {
            std::cout<<"N of Ttracker seeds found: "<<tracksSeed->size()<<std::endl;
            nPotentTracks=tracksSeed->size();

            TrackSeedCollection::const_iterator tracksSeed_it;
            int iseed=0;

            //_hSGSeedHitInBstTmPkEv

            std::multimap< unsigned int, int, std::less<unsigned int> > foundElHitSeed;
            std::multimap< unsigned int, int, std::less<unsigned int> >::reverse_iterator BestFoundElInSeed_it;
            //std::map< int, unsigned int > foundElSeedHit;
            std::map<int, std::vector<StrawHitPtr> > signElGoodHitsInSeeds;
            std::map<int, std::vector<StrawHitPtr> >::iterator sgnElGdHtsInSeeds_it;

            //int sSeedTotalHit, iHitInSeed;//, elHitFoundInCl;

            if (goodEvent && !foundElHitInPeak.empty()) {

                    for ( tracksSeed_it=tracksSeed->begin(); tracksSeed_it!=tracksSeed->end(); ++tracksSeed_it ) {
                            std::cout<<"i-th seed :"<<iseed<<" : "<<std::endl<<*tracksSeed_it;

                            elHitFound=0;
                            std::vector<StrawHitPtr> tmpElHits;
                            //sSeedTotalHit=0;

                            if (tracksSeed_it->_relatedTimeCluster.key()==BestFoundElInPeak_it->second) {
                                    //iHitInSeed=0;
                                    //elHitFoundInCl=0;
                                    for ( std::vector<StrawHitPtr>::const_iterator iSeedHit_it=tracksSeed_it->_selectedTrackerHits.begin(); iSeedHit_it!=tracksSeed_it->_selectedTrackerHits.end(); ++iSeedHit_it ) {
                                            //iHitInSeed++;
                                            for ( std::vector<StrawHitPtr>::iterator iGeTpHit_it=signElGoodHitsInTmpks[BestFoundElInPeak_it->second].begin(); iGeTpHit_it!=signElGoodHitsInTmpks[BestFoundElInPeak_it->second].end(); ++iGeTpHit_it ) {
                                                    //for ( std::vector<StrawHitPtr>::iterator iGoodSelHit_it=signElGoodHit.begin(); iGoodSelHit_it!=signElGoodHit.end(); ++iGoodSelHit_it ) {
                                                    //std::cout<<"Check in Cl "<<iSeed<<" i-th hit "<<iHitInSeed<<" with key "<<iSeedHit_it->key()<<" vs "<<iGeTpHit_it->key()/*iGoodSelHit_it->key()*/<<std::endl;
                                                    if ( iSeedHit_it->key()==iGeTpHit_it->key() ) {
                                                            //if ( iSeedHit_it->key()==iGoodSelHit_it->key() ) {
                                                            elHitFound++;
                                                            //elHitFoundInCl++;
                                                            tmpElHits.push_back(*iSeedHit_it);
                                                            break;
                                                    }
                                            }
                                    }

                                    std::cout<<"Seed elHitFound "<<elHitFound<<std::endl;
                                    if ( elHitFound>0 ) {
                                            std::cout<<"------------ el trk found "<<std::endl;
                                            foundElHitSeed.insert( std::pair<unsigned int, int> (elHitFound,iseed) );
                                            sgnElGdHtsInSeeds_it=signElGoodHitsInSeeds.find(iseed);
                                            //foundElSeedHit.insert( std::pair<int,unsigned int> (iseed,sSeedTotalHit) );
                                            if ( sgnElGdHtsInSeeds_it==signElGoodHitsInSeeds.end() ) signElGoodHitsInSeeds.insert( std::pair< int, std::vector<StrawHitPtr> > (iseed,tmpElHits) );
                                            else if ( elHitFound>sgnElGdHtsInSeeds_it->second.size() ) {
                                                    std::cout<<"------------ removing previous el trk found "<<std::endl;
                                                    signElGoodHitsInSeeds.erase(sgnElGdHtsInSeeds_it);
                                                    signElGoodHitsInSeeds.insert( std::pair< int, std::vector<StrawHitPtr> > (iseed,tmpElHits) );
                                            }
                                    }
                            }

                            iseed++;
                    }
                    //bestSeed_omega, bestSeed_tanDip, bestSeed_d0, bestSeed_phi0, bestSeed_z0, bestSeed_T0;
                    if ( !foundElHitSeed.empty() ) {
                            std::cout<<"------------ histogramming"<<std::endl;
                            _hNBestSeedInTmPkEv->Fill(sel_pt);
                            BestFoundElInSeed_it=foundElHitSeed.rbegin();
                            bestSeedcElNHit=BestFoundElInSeed_it->first;
                            //bestSeedNHit=foundElSeedHit.find(BestFoundElInSeed_it->second)->second;
                            bestSeed_idx=BestFoundElInSeed_it->second;
                            mu2e::TrackSeed &bestSeed = tracksSeed->at(BestFoundElInSeed_it->second);
                            bestSeedNHit=bestSeed._selectedTrackerHits.size();
                            bestSeed_omega = bestSeed.omega();
                            bestSeed_tanDip = bestSeed.tanDip();
                            bestSeed_d0 = bestSeed.d0();
                            bestSeed_phi0 = bestSeed.phi0();
                            bestSeed_z0 = bestSeed.z0();
                            bestSeed_T0 = bestSeed.t0();

                            std::cout<<"--- convElNHit "<<convElNHit<<" bestTPcElNHit "<<bestTPcElNHit<<" bestSeedcElNHit "<<bestSeedcElNHit<<std::endl;
                            //std::cout<<"---  bestSeedNHit "<<foundElSeedHit.find(BestFoundElInSeed_it->second)->second<<" = "<<bestSeedNHit<<std::endl;
                            _hSGSeedHitInBstTmPkEv->Fill(sel_pt,bestSeedcElNHit);
                            _hLostHitBstSGSeedBstTmPkEv->Fill(sel_pt, (BestFoundElInPeak_it->first - bestSeedcElNHit) );
                            std::cout<<"---     !!!!!!!!!  "<<sel_pt<<" - "<<BestFoundElInPeak_it->first<<" - "<<bestSeedcElNHit<<std::endl;
                            _hLostHitBstSGSeedEv->Fill(sel_pt, (signElGoodHit.size() - bestSeedcElNHit) );
                            _hNoiseHitBstSGSeedEv->Fill(sel_pt, (bestSeedNHit - bestSeedcElNHit) );
                            _hNLoopResBstSGSeedEv->Fill(sel_pt, (convElNLoop - bestSeed.nLoops() ) );
                    }
            }
    }

    if (goodEvent) _dataEvalPR->Fill();

//    cerr << "Double click in the canvas_Fake to continue:" ;
//    _fakeCanvas->cd();
//    TLatex *printEvN = new TLatex(0.15,0.4,Form("Current Event: %d",event.id().event()));
//    printEvN->SetTextFont(62);
//    printEvN->SetTextSizePixels(180);
//    printEvN->Draw();
//    _fakeCanvas->Update();
//    _fakeCanvas->WaitPrimitive();
//    cerr << endl;
//    delete printEvN;

    std::cout<<"--------------------------- End of Analysis --------------------"<<std::endl;
    std::cout<<"Event: "<<event.id().event()<<std::endl;
    std::cout<<"----------------------------------------------------------------"<<std::endl;


  } // end analyze

  void EvalPatternRecognition::endJob(){

    // cd() to correct root directory. See note 3.
    TDirectory* save = gDirectory;
    _directory->cd();

    // Write canvas.  See note 4.
//    _canvas->Write();

    // cd() back to where we were.  See note 3.
    save->cd();

  }


}  // end namespace mu2e

using mu2e::EvalPatternRecognition;
DEFINE_ART_MODULE(EvalPatternRecognition);
