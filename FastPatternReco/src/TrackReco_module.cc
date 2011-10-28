//
// identification and track parameter extraction modules
//
// $Id: TrackReco_module.cc,v 1.1 2011/10/28 00:21:08 tassiell Exp $
// $Author: tassiell $
// $Date: 2011/10/28 00:21:08 $
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
#include <cmath>

#include <boost/shared_ptr.hpp>


// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/Handle.h"
#include "art/Persistency/Provenance/Provenance.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


//#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Units/SystemOfUnits.h"

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"
//#include "ITrackerGeom/inc/Cell.hh"
//#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
//#include "MCDataProducts/inc/GenParticleCollection.hh"
//#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
//#include "MCDataProducts/inc/SimParticleCollection.hh"
//#include "MCDataProducts/inc/StepPointMCCollection.hh"
//#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
//#include "ITrackerGeom/inc/ITracker.hh"
#include "TTrackerGeom/inc/TTracker.hh"
//#include "FastPatternReco/inc/TTHitPerTrackData.hh"
//#include "FastPatternReco/inc/GenTrackData.hh"
//#include "MCDataProducts/inc/GenId.hh"
//#include "MCDataProducts/inc/VisibleGenElTrack.hh"
//#include "MCDataProducts/inc/VisibleGenElTrackCollection.hh"
#include "RecoDataProducts/inc/TrackerHitTimeCluster.hh"
#include "RecoDataProducts/inc/TrackerHitTimeClusterCollection.hh"
#include "RecoDataProducts/inc/SectorStationCluster.hh"
#include "RecoDataProducts/inc/SctrSttnClusterGroup.hh"
#include "RecoDataProducts/inc/ZRotStrawHitMap.hh"
#include "RecoDataProducts/inc/ZRotStrawHitMapCollection.hh"
#include "RecoDataProducts/inc/SctrSttnClusterGroupCollection.hh"

// Root includes.
#include "TApplication.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
//#include "TClonesArray.h"
//#include "TObjArray.h"
//#include "TLine.h"
//#include "TArrow.h"
//#include "TEllipse.h"
//#include "TMarker.h"
//#include "TBox.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TF1.h"
#include "TSpectrum2.h"
#include "TLatex.h"
#include "TTree.h"

using namespace std;

namespace mu2e {

  double zo=-1500.0*CLHEP::mm;
  double targetRadMax=100.0*CLHEP::mm;

  typedef art::Ptr<StrawHit> StrawHitPtr;
  //typedef std::map<unsigned int, pair<double, double> > absSectAngleTable;  //key = absSectId, val=( cos(theta), sin(theta) )
  //absSectAngleTable _absSectAngleTablep;
  double _wireLengthStep=2.0*10.0*CLHEP::mm;

  struct HitData{
          const StrawHitPtr _shit;
          const CLHEP::Hep3Vector _direct;
          const CLHEP::Hep3Vector _midp;
          double _radius;
          double _halfLength;
          double _theta;   // ]-pi,pi]
          unsigned int _absZId;
          unsigned int _absSectId;
          vector< pair<double, double> > _wireXYTable;
          /*HitData():
                  _direct(0),
                  _midp(0),
                  _halfLength(0){
          }*/
          HitData(const StrawHitPtr &shit, const CLHEP::Hep3Vector & direct, const CLHEP::Hep3Vector & midp, double halfLength, unsigned int absZId, unsigned int absSectId):
                  _shit(shit),
                  _direct(direct),
                  _midp(midp),
                  _halfLength(halfLength),
                  _absZId(absZId),
                  _absSectId(absSectId){

                  _radius = sqrt( pow(_midp.getX(),2) + pow(_midp.getY(),2) );
                  _theta=_direct.getPhi()-CLHEP::halfpi;
                  if ( _theta<=-CLHEP::pi ) _theta+=CLHEP::twopi;
                  if ( _theta>  CLHEP::pi ) _theta-=CLHEP::twopi;
                  _theta+=CLHEP::pi;

                  //absSectAngleTable::iterator _absSectAngleTablep_it = _absSectAngleTablep.find(_absSectId);
                  //if ( _absSectAngleTablep_it==_absSectAngleTablep.end() ) {
                  //        //double theta=_direct.getPhi();
                  //        //_absSectAngleTablep_it = (_absSectAngleTablep.insert( absSectAngleTable::value_type( _absSectId, pair<double, double>(cos(theta),sin(theta) ) ) ) ).first;
                  //        _absSectAngleTablep_it = (_absSectAngleTablep.insert( absSectAngleTable::value_type( _absSectId, pair<double, double>(_direct.getX(),_direct.getY() ) ) ) ).first;
                  //}
                  _wireXYTable.clear();
                  //_wireXYTable.push_back( pair<double, double>(_midp.getX(), _midp.getY()) );
                  //double tmpL, tmpX, tmpY;
                  //for (int j=1; j<=(int)floor(_halfLength/_wireLengthStep); ++j ) {
                  //        tmpL = ((double) j)*_wireLengthStep;
                  //        tmpX = tmpL*_direct.getX()/*_absSectAngleTablep_it->second.first*/;
                  //        tmpY = tmpL*_direct.getY()/*_absSectAngleTablep_it->second.second*/;
                  //        _wireXYTable.push_back( pair<double, double>(_midp.getX()+tmpX, _midp.getY()+tmpY) );
                  //        _wireXYTable.push_back( pair<double, double>(_midp.getX()-tmpX, _midp.getY()-tmpY) );
                  //}
          }

          void RecalcWXYTab (double stepSize) {
                  //absSectAngleTable::iterator _absSectAngleTablep_it = _absSectAngleTablep.find(_absSectId);
                  _wireLengthStep=stepSize;
                  _wireXYTable.clear();
                  _wireXYTable.push_back( pair<double, double>(_midp.getX(), _midp.getY()) );
                  double tmpL, tmpX, tmpY;
                  for (int j=1; j<=(int)floor(_halfLength/_wireLengthStep); ++j ) {
                          tmpL = ((double) j)*_wireLengthStep;
                          tmpX = tmpL*_direct.getX()/*_absSectAngleTablep_it->second.first*/;
                          tmpY = tmpL*_direct.getY()/*_absSectAngleTablep_it->second.second*/;
                          _wireXYTable.push_back( pair<double, double>(_midp.getX()+tmpX, _midp.getY()+tmpY) );
                          _wireXYTable.push_back( pair<double, double>(_midp.getX()-tmpX, _midp.getY()-tmpY) );
                  }
          }
  };

  typedef std::map<int, HitData *> strawData;  //key = StrawIdx.asInt()

  class TrackReco : public art::EDAnalyzer {
  public:

    explicit TrackReco(fhicl::ParameterSet const& pset);
    virtual ~TrackReco() {
            if (_plotCanvas)        delete _plotCanvas;
            if (_plotCanvas_1)      delete _plotCanvas_1;
            if (_fakeCanvas)        delete _fakeCanvas;
            delete _iXHelStep;
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
//    std::string _makerModuleLabel;

//    // Label of the generator.
//    std::string _generatorModuleLabel;

    // Label of the module that made the hits.
//    std::string _extractElectronsData;

    // Label of the module that made the hits.
//    std::string _timeRejecterModuleLabel;

    // Label of the module that made the hits.
//    std::string _geomRejecterModuleLabel;

    // Label of the module that made the hits.
    std::string _remappingModuleLabel;

    bool _doDisplay;

    // End: run time parameters
    double removeAnglePeriod(double &angle);

    std::vector< std::pair<int, int> > _hitsCouplings; //hits couples by StrawIdx.asInt()
    void computeCombination( double minHStep=0.0, double maxHStep=0.0 );

    std::map< int, std::multimap<int, int> > _hitsTripleCouplings; //hits couples by StrawIdx.asInt()
    void computeTripleCombination( double minHStep=0.0, double maxHStep=0.0 );
    void measureRbyTripleCombination( double HStep, double HPhi0 );

    int butterflyFilter( float *in, float *out, int &nBinX, int &nBinY, float minCountCut );
    int butterflyFilterRot( float *in, float *out, int &nBinX, int &nBinY, float minCountCut );
    int butterflyFilterRot45( float *in, float *out, int &nBinX, int &nBinY, float minCountCut );

    double _minR;
    double _maxR;
    double _stepR;
    double _minLambda;
    double _maxLambda;
    double _stepLambda;
    double _minHelStep;
    double _maxHelStep;
    double _stepHelpStep;
    double _minMom;
    double _maxMom;
    double _rStraw;
    int    _nBinR;
    int    _nBinTanLambda;
    int    _nBinHelpStep;
    double _Bfield;

    double *_iXHelStep;
    strawData *strdat;

    TCanvas*      _plotCanvas;
    TCanvas*      _plotCanvas_1;
    TCanvas*      _fakeCanvas;

    TH2F*         _hRHelStep;
    TH1F*         _hR_1;
    TH1F*         _hR;
    TH1F*         _hPhi0;
    TH2F*         _hPhi0HelStep_1;
    TH2F*         _hPhi0HelStep;
    TH2F*         _hPhi0HelStepL;
    TH2F*         _BF_hPhi0HelStepL;
    TH2F*         _BFRot_hPhi0HelStepL;

    //    // Pointers to histograms, ntuples, TGraphs.
//    TH1F*         _hNtrackableEv;
//    TH1F*         _hNhitTrackableEv;
//    TH1F*         _hNBestTimePeakEv;
//    TH1F*         _hNhitBestTimePeakEv;
//    TH2F*         _hLostHitBstTmPkEv;
//    TH1F*         _hNoiseHitBstTmPkEv;
//
//    TH1F*         _hNBestSGClInTmPkEv;
//    TH1F*         _hSGClHitInBstTmPkEv;
//    TH2F*         _hLostHitBstSGClBstTmPkEv;
//    TH2F*         _hLostHitBstSGClEv;
//    TH1F*         _hNoiseHitBstSGClEv;
//    TH2F*         _hNLoopResBstSGClEv;
//
//    TTree*        _dataEvalBkgRejec;
    unsigned int  runID, eventID;// evNHit, convElNHit, convElNLoop;
//    float         sel_ptMeV, sel_ptMeV_start, sel_ptMeV_end, sel_plMeV, sel_plMeV_start, sel_plMeV_end, convElFHitTime;
//    int           bestTPcElNHit, bestTPNHit, nPeaksFound;
//    int           bestClGcElNHit, bestClGNHit, bestClGcNCls, nPotentTracks;

    // Some ugly but necessary ROOT related bookkeeping:

    // The job needs exactly one instance of TApplication.  See note 1.
    auto_ptr<TApplication> _application;

    // Save directory from beginJob so that we can go there in endJob. See note 3.
    TDirectory* _directory;

  };

  TrackReco::TrackReco(fhicl::ParameterSet const& pset) :

    // Run time parameters
    //_makerModuleLabel(pset.get<string>("makerModuleLabel")),
    //_extractElectronsData(pset.get<string>("elextractModuleLabel")),
    //_timeRejecterModuleLabel(pset.get<string>("tRejecterModuleLabel")),
    //_geomRejecterModuleLabel(pset.get<string>("gRejecterModuleLabel")),
    _remappingModuleLabel(pset.get<string>("reMapperSHModuleLabel")),
//    _moduleLabel(pset.get<string>("module_label")),/*@module_label*/
//    _g4ModuleLabel(pset.get<string>("g4ModuleLabel")),
//    _trackerStepPoints(pset.get<string>("trackerStepPoints","tracker")),
//    _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
    _doDisplay(pset.get<bool>("doDisplay",false)),

    _minR(280.0/2.0*CLHEP::mm),
    _maxR(800.0/2.0*CLHEP::mm),
    _stepR(2.5*1.0*CLHEP::mm),
    _minLambda(30.0*CLHEP::deg),
    _maxLambda(65.0*CLHEP::deg),
    _stepLambda(5.0*CLHEP::mrad),
    _stepHelpStep(2.0*5.0/sqrt(12)*CLHEP::mm),
    _minMom(50.0*CLHEP::MeV),
    _maxMom(150.0*CLHEP::MeV),
    _rStraw(2.50*CLHEP::mm),
    _nBinR(0),
    _nBinTanLambda(0),
    _nBinHelpStep(0),
    _Bfield(1.0),
    _iXHelStep(0x0),

    _plotCanvas(0),
    _plotCanvas_1(0),
    _fakeCanvas(0),
    _hRHelStep(0),
    _hR_1(0),
    _hR(0),
    _hPhi0(0),
    _hPhi0HelStep_1(0),
    _hPhi0HelStep(0),
    _hPhi0HelStepL(0),
    _BF_hPhi0HelStepL(0),
    _BFRot_hPhi0HelStepL(0),
//    _hNtrackableEv(0),
//    _hNhitTrackableEv(0),
//    _hNBestTimePeakEv(0),
//    _hNhitBestTimePeakEv(0),
//    _hLostHitBstTmPkEv(0),
//    _hNoiseHitBstTmPkEv(0),
//
//    _hNBestSGClInTmPkEv(0),
//    _hSGClHitInBstTmPkEv(0),
//    _hLostHitBstSGClBstTmPkEv(0),
//    _hLostHitBstSGClEv(0),
//    _hNoiseHitBstSGClEv(0),
//    _hNLoopResBstSGClEv(0),
//    _dataEvalBkgRejec(0),

    // Some ugly but necessary ROOT related bookkeeping.
    //_application(0),
    _directory(0){
          runID=eventID=/*evNHit=*/0;
          _minHelStep=CLHEP::twopi*_minMom/CLHEP::GeV*sin(_minLambda)/(0.3*_Bfield);
          _minHelStep*=CLHEP::m;
          _minHelStep=floor(_minHelStep);
          _maxHelStep=CLHEP::twopi*_maxMom/CLHEP::GeV*sin(_maxLambda)/(0.3*_Bfield);
          _maxHelStep*=CLHEP::m;
          _maxHelStep=ceil(_maxHelStep);
          //sel_ptMeV=sel_ptMeV_start=sel_ptMeV_end=sel_plMeV=sel_plMeV_start=sel_plMeV_end=convElFHitTime=0.0;
          //bestTPcElNHit=bestTPNHit=nPeaksFound=0;
          //bestClGcElNHit=bestClGNHit=bestClGcNCls=nPotentTracks=0;
 }

  void TrackReco::beginJob(){

          cout<<"Starting TrackReco jos!"<<endl;

    // Get access to the TFile service and save current directory for later use.
    art::ServiceHandle<art::TFileService> tfs;

    if(_doDisplay) {
            // If needed, create the ROOT interactive environment. See note 1.
            if ( !gApplication ){
                    int    tmp_argc(0);
                    char** tmp_argv(0);
                    _application = auto_ptr<TApplication>(new TApplication( "noapplication", &tmp_argc, tmp_argv ));
            }

            gStyle->SetPalette(1);
            gROOT->SetStyle("Plain");


            //_peaksCanvHistos   = new TObjArray();

            _plotCanvas = new TCanvas("plots","Hough Transform plots container",1290,860);
            _plotCanvas->Divide(2,2);
            _plotCanvas_1 = new TCanvas("plots1","Hough Transform plots container",1290,860);
            _plotCanvas_1->Divide(2,2);
           _fakeCanvas = new TCanvas("canvas_Fake","double click for next event",300,100);

    }

    // Create a histogram.
    double tmpNBin = (_maxR-_minR)/_stepR;
    _nBinR = (int)floor(tmpNBin+0.5);
    tmpNBin = (_maxR-_minR)/_stepR;
    _nBinTanLambda = (int)floor(tmpNBin+0.5);
    tmpNBin = (_maxHelStep-_minHelStep)/_stepHelpStep;
    _nBinHelpStep = floor(tmpNBin+0.5);

    cout<<"------------ Helix Step: "<<_minHelStep<<" "<<_maxHelStep<<" "<<_nBinHelpStep<<endl;
    _hRHelStep              = tfs->make<TH2F>( "hRHelStep",   "R vs Helix Step", _nBinHelpStep, _minHelStep, _maxHelStep, /*1200, -100,500*/_nBinR, _minR, _maxR );
    _iXHelStep = new double[_nBinHelpStep];
    for ( int ibin=0; ibin<_nBinHelpStep; ++ibin ) {
            _iXHelStep[ibin] = _hRHelStep->GetXaxis()->GetBinCenter(ibin+1);
    }
    _hPhi0                  = tfs->make<TH1F>( "hPhi0",   "Phi0", 400/*523*//*1047*/, 0.0, CLHEP::twopi/*-CLHEP::pi, CLHEP::pi*/ );
    _hR                     = tfs->make<TH1F>( "hR",   "R", _nBinR, _minR, _maxR );
    _hR_1                   = tfs->make<TH1F>( "hR_1",   "R", _nBinR, _minR, _maxR );
    _hPhi0HelStep_1         = tfs->make<TH2F>( "hPhi0HelStep_1",   "Phi0 vs Helix Step", _nBinHelpStep, _minHelStep, _maxHelStep, 200/*523*//*1047*/, 0.0, CLHEP::twopi/*-CLHEP::pi, CLHEP::pi*/ );


    _hPhi0HelStep           = tfs->make<TH2F>( "hPhi0HelStep",   "Phi0 vs Helix Step", _nBinHelpStep, _minHelStep, _maxHelStep, 100/*523*//*1047*/, -CLHEP::halfpi, CLHEP::halfpi );
    _hPhi0HelStepL          = tfs->make<TH2F>( "hPhi0HelStepL",  "Phi0 vs Helix Step (enlarged)", _nBinHelpStep, _minHelStep, _maxHelStep, 125/*523*//*1047*/, -CLHEP::halfpi, 1.5*CLHEP::halfpi );
    _BF_hPhi0HelStepL       = tfs->make<TH2F>( "BF_hPhi0HelStepL",  "Butterfly Filter of Phi0 vs Helix Step (enlarged)", _nBinHelpStep, _minHelStep, _maxHelStep, 125/*523*//*1047*/, -CLHEP::halfpi, 1.5*CLHEP::halfpi );
    _BFRot_hPhi0HelStepL    = tfs->make<TH2F>( "BFRot_hPhi0HelStepL",  "Butterfly Filter of Phi0 vs Helix Step (enlarged)", _nBinHelpStep, _minHelStep, _maxHelStep, 125/*523*//*1047*/, -CLHEP::halfpi, 1.5*CLHEP::halfpi );

    _hRHelStep->SetXTitle("Step [mm]");
    _hRHelStep->SetYTitle("R [mm]");
    _hPhi0->SetXTitle("#phi_{0} [rad]");
    _hPhi0HelStepL->SetXTitle("Step [mm]");
    _hPhi0HelStepL->SetYTitle("#phi_{0} [rad]");
    _BF_hPhi0HelStepL->SetXTitle("Step [mm]");
    _BF_hPhi0HelStepL->SetYTitle("#phi_{0} [rad]");
    _hR_1->SetXTitle("R [mm]");
    _hR_1->SetYTitle("Entries");



//    _hNtrackableEv       = tfs->make<TH1F>( "hNtrackableEv",   "N of trackable signal electrons", 111, 49.75, 105.25  );
//    _hNhitTrackableEv    = tfs->make<TH1F>( "hNhitTrackableEv",   "N of hits for each signal electrons track", 111, 49.75, 105.25  );
//    _hNBestTimePeakEv    = tfs->make<TH1F>( "hNBestTimePeakEv",   "N of peak time that has the best agreement with trackable signal electrons", 111, 49.75, 105.25  );
//    _hNhitBestTimePeakEv = tfs->make<TH1F>( "hNhitBestTimePeakEv",   "N of hits of the signal electrons track that are found in the best time peak", 111, 49.75, 105.25  );
//    _hLostHitBstTmPkEv   = tfs->make<TH2F>( "hLostHitBstTmPkEv",   "N of lost hits of the signal electrons track that in the best time peak", 111, 49.75, 105.25, 20, 0, 20  );
//    _hNoiseHitBstTmPkEv  = tfs->make<TH1F>( "hNoiseHitBstTmPkEv",   "N of hits of noise present in the the best time peak over signal electrons track hits", 111, 49.75, 105.25  );
//
//    _hNtrackableEv      ->SetXTitle("pt [MeV]");
//    _hNhitTrackableEv   ->SetXTitle("pt [MeV]");
//    _hNBestTimePeakEv   ->SetXTitle("pt [MeV]");
//    _hNhitBestTimePeakEv->SetXTitle("pt [MeV]");
//    _hLostHitBstTmPkEv  ->SetXTitle("pt [MeV]");
//    _hNoiseHitBstTmPkEv ->SetXTitle("pt [MeV]");
//
//    _hNBestSGClInTmPkEv  = tfs->make<TH1F>( "hNBestSGClInTmPkEv",   "N of Best Selected Clusters Group in the best peak time", 111, 49.75, 105.25  );
//    _hSGClHitInBstTmPkEv = tfs->make<TH1F>( "hSGClHitInBstTmPkEv",   "N of signal electron hit selected by best geom clustering algo from the best time peak", 111, 49.75, 105.25  );
//    _hLostHitBstSGClBstTmPkEv = tfs->make<TH2F>( "hLostHitBstSGClBstTmPkEv",   "N of lost hits of the signal electrons track that in the best geom cluster from the best time peak", 111, 49.75, 105.25, 20, 0, 20  );
//    _hLostHitBstSGClEv   = tfs->make<TH2F>( "hLostHitBstSGClEv",   "N of lost hits of the signal electrons track that in the best geom cluster", 111, 49.75, 105.25, 20, 0, 20  );
//    _hNoiseHitBstSGClEv  = tfs->make<TH1F>( "hNoiseHitBstSGClEv",   "N of hits of noise present in the the best geom cluster over signal electrons track hits", 111, 49.75, 105.25  );
//    _hNLoopResBstSGClEv  = tfs->make<TH2F>( "hNLoopResBstSGClEv",   "Res of N of track loop measured by the best geom cluster", 111, 49.75, 105.25, 21, -10.5, 10.5  );
//
//    _hNBestSGClInTmPkEv  ->SetXTitle("pt [MeV]");
//    _hSGClHitInBstTmPkEv ->SetXTitle("pt [MeV]");
//    _hLostHitBstSGClBstTmPkEv ->SetXTitle("pt [MeV]");
//    _hLostHitBstSGClEv ->SetXTitle("pt [MeV]");
//    _hNoiseHitBstSGClEv ->SetXTitle("pt [MeV]");
//    _hNLoopResBstSGClEv ->SetXTitle("pt [MeV]");
//
//    _dataEvalBkgRejec = tfs->make<TTree>( "dataEvalBkgRejec", "data for Bkg Rejection performance evaluation" );
//    _dataEvalBkgRejec->Branch("Run",&runID,"runID/i");
//    _dataEvalBkgRejec->Branch("Event",&eventID,"eventID/i");
//    _dataEvalBkgRejec->Branch("EvTotNHit",&evNHit,"evNHit/i");
//    _dataEvalBkgRejec->Branch("ConvElNHit",&convElNHit,"convElNHit/i");
//    _dataEvalBkgRejec->Branch("ConvElNLoop",&convElNLoop,"convElNLoop/i");
//    _dataEvalBkgRejec->Branch("ConvEl_ptMeV_Tracker",&sel_ptMeV,"sel_ptMeV/F");
//    _dataEvalBkgRejec->Branch("ConvEl_ptMeV_start",&sel_ptMeV_start,"sel_ptMeV_start/F");
//    _dataEvalBkgRejec->Branch("ConvEl_ptMeV_end",&sel_ptMeV_end,"sel_ptMeV_end/F");
//    _dataEvalBkgRejec->Branch("ConvEl_plMeV",&sel_plMeV,"sel_plMeV/F");
//    _dataEvalBkgRejec->Branch("ConvEl_plMeV_start",&sel_plMeV_start,"sel_plMeV_start/F");
//    _dataEvalBkgRejec->Branch("ConvEl_plMeV_end",&sel_plMeV_end,"sel_plMeV_end/F");
//    _dataEvalBkgRejec->Branch("ConvEl_FrstHit_Time",&convElFHitTime,"convElFHitTime/F");
//    _dataEvalBkgRejec->Branch("BestTPcElNHit",&bestTPcElNHit,"bestTPcElNHit/I");
//    _dataEvalBkgRejec->Branch("BestTPNHit",&bestTPNHit,"bestTPNHit/I");
//    _dataEvalBkgRejec->Branch("TotNPeaksFound",&nPeaksFound,"nPeaksFound/I");
//    _dataEvalBkgRejec->Branch("BestClGcElNHit",&bestClGcElNHit,"bestClGcElNHit/I");
//    _dataEvalBkgRejec->Branch("BestClGNHit",&bestClGNHit,"bestClGNHit/I");
//    _dataEvalBkgRejec->Branch("BestClGcNCls",&bestClGcNCls,"bestClGcNCls/I");
//    _dataEvalBkgRejec->Branch("NPotentTracksFound",&nPotentTracks,"nPotentTracks/I");


    //_fakeCanvas = new TCanvas("canvas_Fake","double click for next event",300,100);

    // See note 3.
    _directory = gDirectory;

  }

  void TrackReco::analyze(art::Event const& event ) {


    const Tracker& tracker = getTrackerOrThrow();
    const TTracker &ttr = static_cast<const TTracker&>( tracker );
    const std::vector<Device> ttrdev = ttr.getDevices();

//    art::Handle<StrawHitCollection> pdataHandle;
//    event.getByLabel(_makerModuleLabel,pdataHandle);
//    StrawHitCollection const* hits = pdataHandle.product();
//
//    art::Handle<VisibleGenElTrackCollection> genEltrksHandle;
//    event.getByLabel(_extractElectronsData,genEltrksHandle);
//    VisibleGenElTrackCollection const* genEltrks = genEltrksHandle.product();
//
//    art::Handle<TrackerHitTimeClusterCollection> tclustHandle;
//    event.getByLabel(_timeRejecterModuleLabel,tclustHandle);
//    TrackerHitTimeClusterCollection const* tclusts = tclustHandle.product();
//
//    art::Handle<SctrSttnClusterGroupCollection> gclusgtHandle;
//    event.getByLabel(_geomRejecterModuleLabel,gclusgtHandle);
//    SctrSttnClusterGroupCollection const* gclustgs = gclusgtHandle.product();
//
//    std::vector<mu2e::VisibleGenElTrack>::const_iterator genEltrk_it;

    art::Handle<ZRotStrawHitMapCollection> rmapdataHandle;
    event.getByLabel(_remappingModuleLabel,rmapdataHandle);
    ZRotStrawHitMapCollection const* mhits = rmapdataHandle.product();

    cout<<"--------------------------- Reco Track  ------------------------"<<endl;
    cout<<"Event: "<<event.id().event()<<endl;
    cout<<"----------------------------------------------------------------"<<endl;

    //StrawId sid;

    runID=event.run();
    eventID=event.event();
    //evNHit=hits->size();

    cout<<"N of group of Ttracker cluster of Hit found that could be tracks: "<<mhits->size()<<endl;
    //nPotentTracks=mhits->size();

    ZRotStrawHitMapCollection::const_iterator mhits_it;
    ZSectStrawHitMap::const_iterator zhitmap_tmp_it;
    AbsSectStrawHitMap::const_iterator scthitmap_tmp_it;

    ZStrawHitMap::const_iterator zInscthitmap_tmp1_it, zInscthitmap_tmp2_it;

    int maxStationDist=5, maxContSect=5;
    unsigned int negativeSectOver=2;
    int snegativeSectOver=-((int)negativeSectOver);
    int tmpSecDist;
    bool goodCoupling;

    unsigned int nCoupling, nCouplingInGroup, checkNCoupling;
    std::vector<unsigned int> nCouplingForGroups;

    unsigned int firstStation, tmpStation, lastStation, lastGroupStation;

    /*strawData **/strdat = new strawData;
    //StrawIndex firstStrawIdx, secondStrawIdx;

    for ( mhits_it=mhits->begin(); mhits_it!=mhits->end(); ++mhits_it ) {
            cout<<*mhits_it;
            //for ( ZSectStrawHitMap::const_iterator zhitmap_it = mhits_it->_zsctTrackerHits.begin(); zhitmap_it != mhits_it->_zsctTrackerHits.end(); ++zhitmap_it ) {
            //        for ( AbsSectStrawHitMap::const_iterator scthitmap_it =  zhitmap_it->second.begin(); scthitmap_it !=  zhitmap_it->second.end(); ++scthitmap_it ) {
            //                cout<<"Hit (z,sec): "<<zhitmap_it->first<<" "<<scthitmap_it->first<<endl;
            //        }
            //}
            //cout<<"---------------------------------"<<endl;
            _hitsCouplings.clear();
            _hitsTripleCouplings.clear();
            strdat->clear();

            _hRHelStep->Reset();
            _hPhi0->Reset();
            _hR->Reset();
            _hR_1->Reset();
            _hPhi0HelStep_1->Reset();

            _hPhi0HelStep->Reset();
            _hPhi0HelStepL->Reset();
            _BF_hPhi0HelStepL->Reset();
            _BFRot_hPhi0HelStepL->Reset();

            nCoupling=0;
            nCouplingForGroups.clear();
            for ( SectZStrawHitMap::const_iterator secthitmap_it = mhits_it->_sctZTrackerHits.begin(); secthitmap_it != mhits_it->_sctZTrackerHits.end(); ++secthitmap_it ) {
                    nCouplingInGroup=0;
                    for ( ZStrawHitMap::const_iterator zInscthitmap_it =  secthitmap_it->second.begin(); zInscthitmap_it !=  secthitmap_it->second.end(); ++zInscthitmap_it ) {
                            StrawIndex const&firstStrawIdx=((StrawHitPtr) zInscthitmap_it->second)->strawIndex();
                            const Straw & fstr = ttr.getStraw(firstStrawIdx);
                            if (strdat->count(firstStrawIdx.asInt())==0) strdat->insert(
                                            strawData::value_type(
                                                            firstStrawIdx.asInt(),
                                                            new HitData(zInscthitmap_it->second, fstr.getDirection(), fstr.getMidPoint(), fstr.getHalfLength(), zInscthitmap_it->first, secthitmap_it->first)
                                            )
                            );

                            std::multimap<int, int> tmpCoupling;
                            zInscthitmap_tmp1_it = zInscthitmap_it;
                            ++zInscthitmap_tmp1_it;
                            for ( ; zInscthitmap_tmp1_it !=  secthitmap_it->second.end(); ++zInscthitmap_tmp1_it ) {
                                    StrawIndex const&secondStrawIdx=((StrawHitPtr) zInscthitmap_tmp1_it->second)->strawIndex();
                                    const Straw & sstr = ttr.getStraw(secondStrawIdx);
                                    if ( sstr.getMidPoint().getZ() < fstr.getMidPoint().getZ() ) continue;
                                    if (strdat->count(secondStrawIdx.asInt())==0) strdat->insert(
                                                    strawData::value_type(
                                                                    secondStrawIdx.asInt(),
                                                                    new HitData(zInscthitmap_tmp1_it->second, sstr.getDirection(), sstr.getMidPoint(), sstr.getHalfLength(), zInscthitmap_tmp1_it->first, secthitmap_it->first)
                                                    )
                                    );
                                    _hitsCouplings.push_back(std::pair<int,int>(firstStrawIdx.asInt(),secondStrawIdx.asInt()));

                                    zInscthitmap_tmp2_it = zInscthitmap_tmp1_it;
                                    ++zInscthitmap_tmp2_it;

                                    for ( ; zInscthitmap_tmp2_it !=  secthitmap_it->second.end(); ++zInscthitmap_tmp2_it ) {
                                            StrawIndex const&thirdStrawIdx=((StrawHitPtr) zInscthitmap_tmp2_it->second)->strawIndex();
                                            const Straw & tstr = ttr.getStraw(thirdStrawIdx);
                                            if ( tstr.getMidPoint().getZ() < sstr.getMidPoint().getZ() ) continue;
                                            if (strdat->count(thirdStrawIdx.asInt())==0) strdat->insert(
                                                            strawData::value_type(
                                                                            thirdStrawIdx.asInt(),
                                                                            new HitData(zInscthitmap_tmp2_it->second, tstr.getDirection(), tstr.getMidPoint(), tstr.getHalfLength(), zInscthitmap_tmp2_it->first, secthitmap_it->first)
                                                            )
                                            );

                                            tmpCoupling.insert( std::pair<int,int>( secondStrawIdx.asInt(),thirdStrawIdx.asInt() ) );
                                            ++nCouplingInGroup;
                                            ++nCoupling;
                                    }

                            }
                            _hitsTripleCouplings.insert( std::pair< int,std::multimap<int, int> >( firstStrawIdx.asInt(), tmpCoupling) );
                    }
                    if ( nCouplingInGroup>0 ) {
                            cout<<"N Coupling in previos group "<<nCouplingInGroup<<endl;
                            nCouplingForGroups.push_back(nCouplingInGroup);
                            nCouplingForGroups.push_back(nCouplingInGroup);
                    }
            }

            cout<<"Tot N Coupling "<<nCoupling<<endl;

            computeCombination( mhits_it->_min_HStep, mhits_it->_max_HStep);

            computeTripleCombination( mhits_it->_min_HStep, mhits_it->_max_HStep);

            /*int binMaxX, binMaxY, maxVal;
	    _hPhi0HelStep->GetMaximumBin(binMaxX,binMaxY,maxVal);
	    maxVal=(int)_hPhi0HelStep->GetBinContent(binMaxX,binMaxY);
	    cout<<"Maximun in the Histo:"<<maxVal<<" at bin "<<binMaxX-1<<" "<<binMaxY-1<<endl;
	    measureRbyTripleCombination( _hPhi0HelStep->GetXaxis()->GetBinCenter(binMaxX), _hPhi0HelStep->GetYaxis()->GetBinCenter(binMaxY) );*/

	    Float_t *hPhi0HelStep_arr = _hPhi0HelStep->GetArray();
            Float_t *hPhi0HelStepL_arr = _hPhi0HelStepL->GetArray();
            int nBinsX  = _hPhi0HelStep->GetNbinsX();
            int nBinsY  = _hPhi0HelStep->GetNbinsY();
            int nBinsYL = _hPhi0HelStepL->GetNbinsY();
            memcpy(hPhi0HelStepL_arr+(1*(nBinsX+2)),hPhi0HelStep_arr+(1*(nBinsX+2)),nBinsY*(nBinsX+2)*sizeof(Float_t));
            memcpy(hPhi0HelStepL_arr+((nBinsY+1)*(nBinsX+2)),hPhi0HelStep_arr+1*(nBinsX+2),(nBinsYL-nBinsY)*(nBinsX+2)*sizeof(Float_t));
            _hPhi0HelStepL->SetEntries(_hPhi0HelStep->GetEntries()+(Int_t)_hPhi0HelStep->Integral(1,nBinsX,1,(nBinsYL-nBinsY)));

            Float_t *BF_hPhi0HelStepL_arr = _BF_hPhi0HelStepL->GetArray();
            Float_t *BFRot_hPhi0HelStepL_arr = _BFRot_hPhi0HelStepL->GetArray();
            _BF_hPhi0HelStepL->SetEntries( butterflyFilterRot45( hPhi0HelStepL_arr, BF_hPhi0HelStepL_arr, nBinsX, nBinsYL, 30.0 ) );
            _BFRot_hPhi0HelStepL->SetEntries( butterflyFilterRot( hPhi0HelStepL_arr, BFRot_hPhi0HelStepL_arr, nBinsX, nBinsYL, 30.0 ) );

            int binMaxX, binMaxY, maxVal;
            _BF_hPhi0HelStepL->GetMaximumBin(binMaxX,binMaxY,maxVal);
            maxVal=(int)_BF_hPhi0HelStepL->GetBinContent(binMaxX,binMaxY);
            cout<<"Maximun in the Histo:"<<maxVal<<" at bin "<<binMaxX-1<<" "<<binMaxY-1<<endl;
            measureRbyTripleCombination( _BF_hPhi0HelStepL->GetXaxis()->GetBinCenter(binMaxX), _BF_hPhi0HelStepL->GetYaxis()->GetBinCenter(binMaxY) );

            //TSpectrum2 peaksSearcher(10);
            //int nfound =  peaksSearcher.Search( _hPhi0HelStepL, 2.0, "", 0.1 );
            //TSpectrum2 peaksSearcherBF(10);
            //int nfoundBF =  peaksSearcherBF.Search( _BF_hPhi0HelStepL, 2.0, "", 0.1 );

            //cout<<"Peaks found in Voting Array "<<nfound<<" in BFiltered "<<nfoundBF<<endl;

            //TH1D *_BF_hPhi0HelStepL_Y = _BF_hPhi0HelStepL->ProjectionY();



//            _hRHelStep->Reset();

//            nCoupling=0;
//            nCouplingInGroup=0;
//            lastGroupStation=0;
//            nCouplingForGroups.clear();
//            for ( ZSectStrawHitMap::const_iterator zhitmap_it = mhits_it->_zsctTrackerHits.begin(); zhitmap_it != mhits_it->_zsctTrackerHits.end(); ++zhitmap_it ) {
//                    firstStation = (unsigned int)zhitmap_it->first/4;
//                    //cout<<"******** firstStation "<<firstStation<<" lastGroupStation "<<lastGroupStation<<endl;
//                    if ( nCoupling>0 && ((int)(firstStation-lastGroupStation))>0 ) {
//                            cout<<"N Coupling in previos group "<<nCouplingInGroup<<endl;
//                            nCouplingForGroups.push_back(nCouplingInGroup);
//                            nCouplingInGroup=0;
//                    }
//                    //goodCoupling=false;
//                    for ( AbsSectStrawHitMap::const_iterator scthitmap_it =  zhitmap_it->second.begin(); scthitmap_it !=  zhitmap_it->second.end(); ++scthitmap_it ) {
//                            //cout<<"Starting hit at z "<<zhitmap_it->first<<" sec "<<scthitmap_it->first<<endl;
//                            lastStation=firstStation;
//                            //if (lastStation>lastGroupStation) lastGroupStation=lastStation;
//                            zhitmap_tmp_it  = zhitmap_it;
//                            scthitmap_tmp_it = scthitmap_it;
//                            ++scthitmap_tmp_it;
//                            goodCoupling=false;
//                            //cout<<"Coupling Z: "<<zhitmap_it->first<<" - "<<zhitmap_tmp_it->first<<":"<<endl;
//                            for ( ; scthitmap_tmp_it != zhitmap_tmp_it->second.end(); ++scthitmap_tmp_it) {
//                                    if (scthitmap_tmp_it->first != scthitmap_it->first) break;
//                                    //tmpSecDist = scthitmap_tmp_it->first - scthitmap_it->first;
//                                    //if ( tmpSecDist<snegativeSectOver || tmpSecDist>=maxContSect ) break;
//                                    //cout<<"\t Sect: "<<scthitmap_it->first<<" - "<<scthitmap_tmp_it->first<<endl;
//                                    goodCoupling=true;
//                                    nCoupling++;
//                                    nCouplingInGroup++;
//
//                                    StrawIndex const&firstStrawIdx=((StrawHitPtr) scthitmap_it->second)->strawIndex();
//                                    StrawIndex const&secondStrawIdx=((StrawHitPtr) scthitmap_tmp_it->second)->strawIndex();
//                                    _hitsCouplings.push_back(std::pair<int,int>(firstStrawIdx.asInt(),secondStrawIdx.asInt()));
//                                    const Straw & fstr = ttr.getStraw(firstStrawIdx);
//                                    if (strdat->count(firstStrawIdx.asInt())==0) strdat->insert(
//                                                    strawData::value_type(
//                                                                    firstStrawIdx.asInt(),
//                                                                    new HitData(scthitmap_it->second, fstr.getDirection(), fstr.getMidPoint(), fstr.getHalfLength(), zhitmap_it->first, scthitmap_it->first)
//                                                    )
//                                    );
//                                    const Straw & sstr = ttr.getStraw(secondStrawIdx);
//                                    if (strdat->count(secondStrawIdx.asInt())==0) strdat->insert(
//                                                    strawData::value_type(
//                                                                    secondStrawIdx.asInt(),
//                                                                    new HitData(scthitmap_it->second, sstr.getDirection(), sstr.getMidPoint(), sstr.getHalfLength(), zhitmap_tmp_it->first, scthitmap_tmp_it->first)
//                                                    )
//                                    );
//
//                            }
//                            ++zhitmap_tmp_it;
//
//                            for ( ; zhitmap_tmp_it != mhits_it->_zsctTrackerHits.end(); ++zhitmap_tmp_it ) {
//                                    tmpStation = (unsigned int)zhitmap_tmp_it->first/4;
//                                    if ( ((int)(tmpStation - lastStation))>1 || ((int)(tmpStation - firstStation))>maxStationDist ) break;
//                                    goodCoupling=false;
//                                    //cout<<"Coupling Z: "<<zhitmap_it->first<<" - "<<zhitmap_tmp_it->first<<":"<<endl;
//                                    //cout<<"("<<"tmpStation "<<tmpStation<<" lastStation "<<lastStation<<" firstStation "<<firstStation<<")"<<endl;
//
//                                    if (scthitmap_it->first<negativeSectOver) {
//                                            scthitmap_tmp_it = zhitmap_tmp_it->second.lower_bound( 0 );
//                                            for ( ; scthitmap_tmp_it != zhitmap_tmp_it->second.end(); ++scthitmap_tmp_it) {
//                                                    tmpSecDist = scthitmap_tmp_it->first - scthitmap_it->first;
//                                                    if ( tmpSecDist<snegativeSectOver || tmpSecDist>=maxContSect ) break;
//                                                    //cout<<"\t Sect: "<<scthitmap_it->first<<" - "<<scthitmap_tmp_it->first<<endl;
//                                                    goodCoupling=true;
//                                                    nCoupling++;
//                                                    nCouplingInGroup++;
//
//                                                    StrawIndex const&firstStrawIdx=((StrawHitPtr) scthitmap_it->second)->strawIndex();
//                                                    StrawIndex const&secondStrawIdx=((StrawHitPtr) scthitmap_tmp_it->second)->strawIndex();
//                                                    _hitsCouplings.push_back(std::pair<int,int>(firstStrawIdx.asInt(),secondStrawIdx.asInt()));
//                                                    const Straw & fstr = ttr.getStraw(firstStrawIdx);
//                                                    if (strdat->count(firstStrawIdx.asInt())==0) strdat->insert(
//                                                                    strawData::value_type(
//                                                                                    firstStrawIdx.asInt(),
//                                                                                    new HitData(scthitmap_it->second, fstr.getDirection(), fstr.getMidPoint(), fstr.getHalfLength(), zhitmap_it->first, scthitmap_it->first)
//                                                                    )
//                                                    );
//                                                    const Straw & sstr = ttr.getStraw(secondStrawIdx);
//                                                    if (strdat->count(secondStrawIdx.asInt())==0) strdat->insert(
//                                                                    strawData::value_type(
//                                                                                    secondStrawIdx.asInt(),
//                                                                                    new HitData(scthitmap_it->second, sstr.getDirection(), sstr.getMidPoint(), sstr.getHalfLength(), zhitmap_tmp_it->first, scthitmap_tmp_it->first)
//                                                                    )
//                                                    );
//
//                                           }
//                                            scthitmap_tmp_it = zhitmap_tmp_it->second.lower_bound( 12-negativeSectOver );
//                                            for ( ; scthitmap_tmp_it != zhitmap_tmp_it->second.end(); ++scthitmap_tmp_it) {
//                                                    tmpSecDist = scthitmap_it->first /*+1*/-(scthitmap_tmp_it->first-12);
//                                                    if ( tmpSecDist<snegativeSectOver || tmpSecDist>=maxContSect ) break;
//                                                    //cout<<"\t Sect: "<<scthitmap_it->first<<" - "<<scthitmap_tmp_it->first<<endl;
//                                                    goodCoupling=true;
//                                                    nCoupling++;
//                                                    nCouplingInGroup++;
//
//                                                    StrawIndex const&firstStrawIdx=((StrawHitPtr) scthitmap_it->second)->strawIndex();
//                                                    StrawIndex const&secondStrawIdx=((StrawHitPtr) scthitmap_tmp_it->second)->strawIndex();
//                                                    _hitsCouplings.push_back(std::pair<int,int>(firstStrawIdx.asInt(),secondStrawIdx.asInt()));
//                                                    const Straw & fstr = ttr.getStraw(firstStrawIdx);
//                                                    if (strdat->count(firstStrawIdx.asInt())==0) strdat->insert(
//                                                                    strawData::value_type(
//                                                                                    firstStrawIdx.asInt(),
//                                                                                    new HitData(scthitmap_it->second, fstr.getDirection(), fstr.getMidPoint(), fstr.getHalfLength(), zhitmap_it->first, scthitmap_it->first)
//                                                                    )
//                                                    );
//                                                    const Straw & sstr = ttr.getStraw(secondStrawIdx);
//                                                    if (strdat->count(secondStrawIdx.asInt())==0) strdat->insert(
//                                                                    strawData::value_type(
//                                                                                    secondStrawIdx.asInt(),
//                                                                                    new HitData(scthitmap_it->second, sstr.getDirection(), sstr.getMidPoint(), sstr.getHalfLength(), zhitmap_tmp_it->first, scthitmap_tmp_it->first)
//                                                                    )
//                                                    );
//
//                                           }
//                                    }
//                                    else {
//                                            scthitmap_tmp_it = zhitmap_tmp_it->second.lower_bound( scthitmap_it->first-negativeSectOver );
//                                            for ( ; scthitmap_tmp_it != zhitmap_tmp_it->second.end(); ++scthitmap_tmp_it) {
//                                                    tmpSecDist = scthitmap_tmp_it->first - scthitmap_it->first;
//                                                    if ( tmpSecDist<snegativeSectOver || tmpSecDist>=maxContSect ) break;
//                                                    //cout<<"\t Sect: "<<scthitmap_it->first<<" - "<<scthitmap_tmp_it->first<<endl;
//                                                    goodCoupling=true;
//                                                    nCoupling++;
//                                                    nCouplingInGroup++;
//
//                                                    StrawIndex const&firstStrawIdx=((StrawHitPtr) scthitmap_it->second)->strawIndex();
//                                                    StrawIndex const&secondStrawIdx=((StrawHitPtr) scthitmap_tmp_it->second)->strawIndex();
//                                                    _hitsCouplings.push_back(std::pair<int,int>(firstStrawIdx.asInt(),secondStrawIdx.asInt()));
//                                                    const Straw & fstr = ttr.getStraw(firstStrawIdx);
//                                                    if (strdat->count(firstStrawIdx.asInt())==0) strdat->insert(
//                                                                    strawData::value_type(
//                                                                                    firstStrawIdx.asInt(),
//                                                                                    new HitData(scthitmap_it->second, fstr.getDirection(), fstr.getMidPoint(), fstr.getHalfLength(), zhitmap_it->first, scthitmap_it->first)
//                                                                    )
//                                                    );
//                                                    const Straw & sstr = ttr.getStraw(secondStrawIdx);
//                                                    if (strdat->count(secondStrawIdx.asInt())==0) strdat->insert(
//                                                                    strawData::value_type(
//                                                                                    secondStrawIdx.asInt(),
//                                                                                    new HitData(scthitmap_it->second, sstr.getDirection(), sstr.getMidPoint(), sstr.getHalfLength(), zhitmap_tmp_it->first, scthitmap_tmp_it->first)
//                                                                    )
//                                                    );
//
//                                            }
//                                    }
//                                    if ( scthitmap_it->first>(12-maxStationDist) ) {
//                                            scthitmap_tmp_it = zhitmap_tmp_it->second.lower_bound( 0 );
//                                            for ( ; scthitmap_tmp_it != zhitmap_tmp_it->second.end() && (scthitmap_tmp_it->first - (scthitmap_it->first-12))<maxContSect ; ++scthitmap_tmp_it) {
//                                                    //cout<<"\t Sect: "<<scthitmap_it->first<<" - "<<scthitmap_tmp_it->first<<endl;
//                                                    goodCoupling=true;
//                                                    nCoupling++;
//                                                    nCouplingInGroup++;
//
//                                                    StrawIndex const&firstStrawIdx=((StrawHitPtr) scthitmap_it->second)->strawIndex();
//                                                    StrawIndex const&secondStrawIdx=((StrawHitPtr) scthitmap_tmp_it->second)->strawIndex();
//                                                    _hitsCouplings.push_back(std::pair<int,int>(firstStrawIdx.asInt(),secondStrawIdx.asInt()));
//                                                    const Straw & fstr = ttr.getStraw(firstStrawIdx);
//                                                    if (strdat->count(firstStrawIdx.asInt())==0) strdat->insert(
//                                                                    strawData::value_type(
//                                                                                    firstStrawIdx.asInt(),
//                                                                                    new HitData(scthitmap_it->second, fstr.getDirection(), fstr.getMidPoint(), fstr.getHalfLength(), zhitmap_it->first, scthitmap_it->first)
//                                                                    )
//                                                    );
//                                                    const Straw & sstr = ttr.getStraw(secondStrawIdx);
//                                                    if (strdat->count(secondStrawIdx.asInt())==0) strdat->insert(
//                                                                    strawData::value_type(
//                                                                                    secondStrawIdx.asInt(),
//                                                                                    new HitData(scthitmap_it->second, sstr.getDirection(), sstr.getMidPoint(), sstr.getHalfLength(), zhitmap_tmp_it->first, scthitmap_tmp_it->first)
//                                                                    )
//                                                    );
//
//                                            }
//                                   }
//                                   if (goodCoupling) {
//                                           lastStation=tmpStation;
//                                           if (lastStation>lastGroupStation) lastGroupStation=lastStation;
//                                           //cout<<"--------- lastStation "<<lastStation<<" lastGroupStation "<<lastGroupStation<<endl;
//                                   }
//                            }
//                    }
//                    //if (goodCoupling) lastStation=firstStation;
//            }
//            if ( nCouplingInGroup>0 ) {
//                    cout<<"N Coupling in previos group "<<nCouplingInGroup<<endl;
//                    nCouplingForGroups.push_back(nCouplingInGroup);
//            }
//
//            cout<<endl<<"N of total hit coupling for hits gorup: "<<nCoupling<<" = to n poit "<<(1+sqrt(1+8*nCoupling))/2<<endl<<endl;
//            checkNCoupling=0;
//            for (std::vector<unsigned int>::iterator nCouplingForGroups_it=nCouplingForGroups.begin(); nCouplingForGroups_it!=nCouplingForGroups.end(); ++nCouplingForGroups_it){
//                    checkNCoupling+=* nCouplingForGroups_it;
//            }
//            if (nCouplingForGroups.size()>1 && checkNCoupling==nCoupling) cout<<"Good separtion in groups!!!!"<<endl;
//
//            //for ( strawData::iterator  strdat_it=strdat->begin(); strdat_it!=strdat->end(); ++strdat_it ) {
//            //        cout<<strdat_it->first <<" absZ "<<strdat_it->second->_absZId<<" absSect "<<strdat_it->second->_absSectId<<" dir: "<<strdat_it->second->_direct<<" mid "<<strdat_it->second->_midp<<endl;
//            //        cout<<"wire points (x,y):"<<endl;
//            //        for ( vector< pair<double, double> >::iterator wxy_it = strdat_it->second->_wireXYTable.begin(); wxy_it != strdat_it->second->_wireXYTable.end(); ++wxy_it ) {
//            //                cout<<"\t"<<wxy_it->first<<" , "<<wxy_it->second<<endl;
//            //        }
//            //}
//
//            computeCombination( mhits_it->_min_HStep, mhits_it->_max_HStep);
//            int binMaxX, binMaxY, maxVal;
//            _hRHelStep->GetMaximumBin(binMaxX,binMaxY,maxVal);
//            maxVal=(int)_hRHelStep->GetBinContent(binMaxX,binMaxY);
//            cout<<"Maximun in the Histo:"<<maxVal<<" at bin "<<binMaxX-1<<" "<<binMaxY-1<<endl;

            if (_doDisplay) {

                    if (mhits_it->_min_HStep>0.0 && mhits_it->_max_HStep>0.0) {
                            _hRHelStep->GetXaxis()->SetRangeUser( mhits_it->_min_HStep, mhits_it->_max_HStep);
                            _hPhi0HelStep_1->GetXaxis()->SetRangeUser( mhits_it->_min_HStep, mhits_it->_max_HStep);

                            _hPhi0HelStep->GetXaxis()->SetRangeUser( mhits_it->_min_HStep, mhits_it->_max_HStep);
                            _hPhi0HelStepL->GetXaxis()->SetRangeUser( mhits_it->_min_HStep, mhits_it->_max_HStep);
                            _BF_hPhi0HelStepL->GetXaxis()->SetRangeUser( mhits_it->_min_HStep, mhits_it->_max_HStep);
                            _BFRot_hPhi0HelStepL->GetXaxis()->SetRangeUser( mhits_it->_min_HStep, mhits_it->_max_HStep);
                    }

                    //_plotCanvas->cd();
                    //_hRHelStep->Draw("col z");
                    _plotCanvas->cd(1);
                    //_hRHelStep->Draw("col z");
                    //_hPhi0HelStep->Draw("col z");
                    _hPhi0HelStepL->Draw("col z");
                    _plotCanvas->cd(2);
                    //_hPhi0HelStepL->Draw("col z");
                    _BF_hPhi0HelStepL->Draw("col z");
                    _plotCanvas->cd(3);
                    _hR_1->Draw();
                    //_hPhi0->Draw();
                    //_BF_hPhi0HelStepL_Y->Draw();
                    //_BF_hPhi0HelStepL->Draw("col z");
                    _plotCanvas->cd(4);
                    //_BF_hPhi0HelStepL->Draw("col z");
                    //_BFRot_hPhi0HelStepL->Draw("col z");
                    _plotCanvas->Update();

                    _plotCanvas_1->cd(1);
                    _hRHelStep->Draw("col z");
                    _plotCanvas_1->cd(2);
                    _hPhi0HelStep_1->Draw("col z");
                    _plotCanvas_1->cd(3);
                    _hR->Draw();
                    _plotCanvas_1->cd(4);
                    _hPhi0->Draw();
                    _plotCanvas_1->Update();

                    cerr << "Double click in the canvas_Fake to continue:" ;
                    _fakeCanvas->Clear();
                    _fakeCanvas->cd();
                    TLatex *printEvN = new TLatex(0.15,0.4,Form("Current Event: %d",event.id().event()));
                    printEvN->SetTextFont(62);
                    printEvN->SetTextSizePixels(180);
                    printEvN->Draw();
                    _fakeCanvas->Update();
                    _fakeCanvas->WaitPrimitive();
                    cerr << endl;
                    delete printEvN;
            }


    }

    strdat->clear();
    delete strdat;


    cout<<"--------------------------- End of Reco Track ------------------"<<endl;
    cout<<"Event: "<<event.id().event()<<endl;
    cout<<"----------------------------------------------------------------"<<endl;


  } // end analyze

  void TrackReco::endJob(){

    // cd() to correct root directory. See note 3.
    TDirectory* save = gDirectory;
    _directory->cd();

    // Write canvas.  See note 4.
//    _canvas->Write();

    // cd() back to where we were.  See note 3.
    save->cd();

  }

  double TrackReco::removeAnglePeriod(double &angle){
          // remove the 2pi period and return the angle in the range ]-pi,pi]
          double intpart, fractpart;
          double tmp=angle/CLHEP::twopi;
          fractpart=std::modf (tmp , &intpart);
          fractpart*=CLHEP::twopi;
          if (fractpart<0.0) fractpart+=CLHEP::twopi;
          //if (fractpart<-CLHEP::pi) fractpart+=CLHEP::twopi;
          //if (fractpart> CLHEP::pi) fractpart-=CLHEP::twopi;
          //fractpart+=CLHEP::pi;

          //cout<<"Rescaling angle in: "<<angle/CLHEP::degree<<" out: "<<fractpart/CLHEP::degree<<endl;
          return fractpart;
  }

  void TrackReco::computeCombination( double minHStep, double maxHStep ){
          cout<<"In TrackReco::computeCombination"<<endl;
          cout<<"-------- minHStep "<<minHStep<<" maxHStep "<<maxHStep<<endl;
          int firstXbin=0, lastXBin=_nBinHelpStep;
          if ( minHStep>_minHelStep) {
                  for (int ibin=0; ibin<_nBinHelpStep; ++ibin){
                          if (minHStep<=_iXHelStep[ibin]) {
                                  firstXbin=ibin;
                                  break;
                          }
                  }
          }
          if ( maxHStep<_maxHelStep ) {
                  lastXBin=firstXbin;
                  for (int ibin=firstXbin; ibin<_nBinHelpStep; ++ibin){
                          if (maxHStep<=_iXHelStep[ibin]) {
                                  break;
                          }
                          ++lastXBin;
                  }
          }
          cout<<"-------- firstXbin "<<firstXbin<<" lastXBin "<<lastXBin<<endl;
          //firstXbin=0;


          double tmpR1, tmpDen1, tmpNum1, tmpR2, tmpDen2, tmpNum2;
          double tmpXo, tmpYo, tmpXC, tmpYC, tmpCdist, tmpCspread;
          double DeltaZ1, DeltaZ2;
          double Phi01, Phi0lTheta1, Phi1, Phi02, Phi0lTheta2, Phi2;
          double cosPhi01, sinPhi01, cosPhi02, sinPhi02;
          double tmpVal;
          double tmpNBin, stepXY, maxXY = _maxR-targetRadMax;
          tmpNBin = 2.0*maxXY/_stepR;
          int nXYStep = (int)floor(tmpNBin+0.5);
          stepXY=2.0*maxXY/((double)nXYStep);
          double tmpTwopiTanLambda1, tmpTwopiTanLambda2, twopiTanLambdaMin = CLHEP::twopi*tan(_minLambda), twoPitanLambdaMax = CLHEP::twopi*tan(_maxLambda);

          for ( vector< std::pair<int, int> >::iterator hitsCouplings_it = _hitsCouplings.begin(); hitsCouplings_it != _hitsCouplings.end(); ++hitsCouplings_it ) {
                  //cout<<"I'm analyzing the couple: "<<hitsCouplings_it->first<<" "<<hitsCouplings_it->second<<endl;

                  HitData *first  = strdat->find(hitsCouplings_it->first)->second;
                  HitData *second = strdat->find(hitsCouplings_it->second)->second;
                  DeltaZ1 = first->_midp.getZ() - zo;
                  DeltaZ2 = second->_midp.getZ() - zo;
                  //cout<<"AbsSect "<<first->_absSectId<<" theta "<<first->_theta<<endl;
                  for ( int ibin=firstXbin; ibin<lastXBin; ++ibin ) {
                          tmpVal      = CLHEP::twopi/_iXHelStep[ibin];

                          Phi1        = tmpVal*DeltaZ1;
                          Phi1        = removeAnglePeriod(Phi1);
                          //if (Phi1==0.0 || /*Phi1==-CLHEP::pi ||*/Phi1==CLHEP::twopi || Phi1==CLHEP::pi) {Phi0lTheta1 = Phi1+CLHEP::halfpi; cout<<"No der"<<endl;}
                          if (Phi1==0.0 || /*Phi1==-CLHEP::pi ||*/Phi1==CLHEP::twopi ) Phi0lTheta1 = CLHEP::halfpi;
                          else if (Phi1==CLHEP::pi ) Phi0lTheta1 = -CLHEP::halfpi;
                          else {
                                  //Phi0lTheta1  = std::atan(2.0*std::tan(Phi1+CLHEP::halfpi));
                                  //if (Phi1>=0.0 && Phi1<CLHEP::halfpi) Phi0lTheta1+=CLHEP::pi;
                                  //if (Phi1>CLHEP::halfpi && Phi1<CLHEP::pi) Phi0lTheta1-=CLHEP::pi;
                                  Phi0lTheta1  = std::atan(2.0*std::tan(Phi1+CLHEP::halfpi));
                                  if (Phi1<=CLHEP::pi) Phi0lTheta1+=CLHEP::pi;
                          }
                          Phi0lTheta1=removeAnglePeriod(Phi0lTheta1);
                          tmpVal = -Phi1-Phi0lTheta1;
                          tmpVal = removeAnglePeriod(tmpVal);
                          if (tmpVal>=CLHEP::pi) {
                                  if (Phi1>CLHEP::pi/*Phi1>=0.0&&Phi1<CLHEP::pi*/) continue;
                                  Phi0lTheta1=0.0;
                          }
                          Phi01 = Phi0lTheta1+first->_theta;
                          Phi01 = removeAnglePeriod(Phi01);
                          _hPhi0HelStep_1->Fill(_iXHelStep[ibin],Phi01);
                          cosPhi01 = cos(Phi01);
                          sinPhi01 = sin(Phi01);

                          Phi2        = tmpVal*DeltaZ2;
                          Phi2        = removeAnglePeriod(Phi2);
                          //if (Phi2==0.0 || /*Phi2==-CLHEP::pi ||*/Phi2==CLHEP::twopi || Phi2==CLHEP::pi) {Phi0lTheta2 = Phi2+CLHEP::halfpi; cout<<"No der"<<endl;}
                          if (Phi2==0.0 || /*Phi1==-CLHEP::pi ||*/Phi2==CLHEP::twopi ) Phi0lTheta2 = CLHEP::halfpi;
                          else if (Phi2==CLHEP::pi ) Phi0lTheta2 = -CLHEP::halfpi;
                          else {
                                  //Phi0lTheta2  = std::atan(2.0*std::tan(Phi2+CLHEP::halfpi));
                                  //if (Phi2>=0.0 && Phi2<CLHEP::halfpi) Phi0lTheta2+=CLHEP::pi;
                                  //if (Phi2>CLHEP::halfpi && Phi2<CLHEP::pi) Phi0lTheta2-=CLHEP::pi;
                                  Phi0lTheta2  = std::atan(2.0*std::tan(Phi2+CLHEP::halfpi));
                                  if (Phi2<=CLHEP::pi) Phi0lTheta2+=CLHEP::pi;
                          }
                          Phi0lTheta2=removeAnglePeriod(Phi0lTheta2);
                          tmpVal = -Phi2-Phi0lTheta2;
                          tmpVal = removeAnglePeriod(tmpVal);
                          if (tmpVal>=CLHEP::pi) {
                                  if (Phi2>CLHEP::pi/*Phi2>=0.0*/) continue;
                                  Phi0lTheta2=0.0;
                          }
                          Phi02 = Phi0lTheta2+second->_theta;
                          Phi02 = removeAnglePeriod(Phi02);
                          _hPhi0HelStep_1->Fill(_iXHelStep[ibin],Phi02);
                          cosPhi02 = cos(Phi02);
                          sinPhi02 = sin(Phi02);

                          tmpDen1 = 1.0/( cos(Phi0lTheta1)*(cos(Phi1)-1.0) - sin(Phi0lTheta1)*sin(Phi1) );
                          tmpDen2 = 1.0/( cos(Phi0lTheta2)*(cos(Phi2)-1.0) - sin(Phi0lTheta2)*sin(Phi2) );
                          tmpXo=-maxXY+0.5*stepXY;
/*                          for ( int ix=0; ix<nXYStep; ++ix ) {
                                  tmpYo=-maxXY+0.5*stepXY;
                                  for ( int iy=0; iy<nXYStep; ++iy ) {
                                          tmpVal = tmpXo*first->_direct.getX() + tmpYo*first->_direct.getX();
                                          tmpNum1 = first->_radius - tmpVal;
                                          tmpR1   = tmpNum1*tmpDen1;
                                          if (tmpR1<=0.0) continue;
                                          tmpTwopiTanLambda1 =  _iXHelStep[ibin]/tmpR1;
                                          if ( tmpTwopiTanLambda1>twoPitanLambdaMax || (tmpTwopiTanLambda1>-twopiTanLambdaMin && tmpTwopiTanLambda1<twopiTanLambdaMin) || tmpTwopiTanLambda1<-twoPitanLambdaMax ) continue;
                                          tmpXC = tmpXo - tmpR1*cosPhi01;
                                          tmpYC = tmpYo - tmpR1*sinPhi01;
                                          tmpCdist = sqrt( pow(tmpXC,2) + pow(tmpYC,2) );
                                          tmpCspread = tmpCdist - tmpR1;
                                          if ( tmpCspread>-targetRadMax && tmpCspread<targetRadMax ) {
                                                  _hRHelStep->Fill(_iXHelStep[ibin],tmpR1);
                                                  _hPhi0->Fill(Phi01);
                                                  _hR->Fill(tmpR1);
                                                  _hPhi0HelStep_1->Fill(_iXHelStep[ibin],Phi01);

                                                  tmpNum2 = second->_radius - tmpVal;
                                                  tmpR2   = tmpNum2*tmpDen2;
                                                  if (tmpR2<=0.0) continue;
                                                  tmpTwopiTanLambda2 =  _iXHelStep[ibin]/tmpR2;
                                                  if ( tmpTwopiTanLambda2>twoPitanLambdaMax || (tmpTwopiTanLambda2>-twopiTanLambdaMin && tmpTwopiTanLambda2<twopiTanLambdaMin) || tmpTwopiTanLambda2<-twoPitanLambdaMax ) continue;
                                                  tmpXC = tmpXo - tmpR2*cosPhi02;
                                                  tmpYC = tmpYo - tmpR2*sinPhi02;
                                                  tmpCdist = sqrt( pow(tmpXC,2) + pow(tmpYC,2) );
                                                  tmpCspread = tmpCdist - tmpR2;
                                                  if ( tmpCspread>-targetRadMax && tmpCspread<targetRadMax ) {
                                                          _hRHelStep->Fill(_iXHelStep[ibin],tmpR2);
                                                          _hPhi0->Fill(Phi02);
                                                          _hR->Fill(tmpR2);
                                                          _hPhi0HelStep_1->Fill(_iXHelStep[ibin],Phi02);

                                                  }
                                          }
                                          tmpYo+=stepXY;
                                  }
                                  tmpXo+=stepXY;
                          }*/
                  }
          }
  }

//  void TrackReco::computeCombination( double minHStep, double maxHStep ){
//          cout<<"-------- minHStep "<<minHStep<<" maxHStep "<<maxHStep<<endl;
//          int firstXbin=0, lastXBin=_nBinHelpStep;
//          if ( minHStep>_minHelStep) {
//                  for (int ibin=0; ibin<_nBinHelpStep; ++ibin){
//                          if (minHStep<=_iXHelStep[ibin]) {
//                                  firstXbin=ibin;
//                                  break;
//                          }
//                  }
//          }
//          if ( maxHStep<_maxHelStep ) {
//                  lastXBin=firstXbin;
//                  for (int ibin=firstXbin; ibin<_nBinHelpStep; ++ibin){
//                          if (maxHStep<=_iXHelStep[ibin]) {
//                                  break;
//                          }
//                          ++lastXBin;
//                  }
//          }
//          cout<<"-------- firstXbin "<<firstXbin<<" lastXBin "<<lastXBin<<endl;
//          //firstXbin=0;
//
//
//          double DeltaZ, chord, tmpR;
//          for ( vector< std::pair<int, int> >::iterator hitsCouplings_it = _hitsCouplings.begin(); hitsCouplings_it != _hitsCouplings.end(); ++hitsCouplings_it ) {
//                  //cout<<"I'm analyzing the couple: "<<hitsCouplings_it->first<<" "<<hitsCouplings_it->second<<endl;
//
//                  HitData *first  = strdat->find(hitsCouplings_it->first)->second;
//                  HitData *second = strdat->find(hitsCouplings_it->second)->second;
//                  DeltaZ = second->_midp.getZ() - first->_midp.getZ();
//                  for ( vector< pair<double, double> >::iterator fWXY_it = first->_wireXYTable.begin(); fWXY_it != first->_wireXYTable.end(); ++fWXY_it ) {
//                          for ( vector< pair<double, double> >::iterator sWXY_it = second->_wireXYTable.begin(); sWXY_it != second->_wireXYTable.end(); ++sWXY_it ) {
//                                  chord = sqrt( pow( sWXY_it->first - fWXY_it->first ,2) + pow( sWXY_it->second - fWXY_it->second ,2) );
//                                  for ( int ibin=firstXbin; ibin<lastXBin; ++ibin ) {
//                                          tmpR=2.0*sin(CLHEP::pi*DeltaZ/_iXHelStep[ibin]);
//                                          tmpR=chord/tmpR;
//                                          _hRHelStep->Fill(_iXHelStep[ibin],tmpR);                                  }
//                          }
//                  }
//          }
//  }

  void TrackReco::computeTripleCombination( double minHStep, double maxHStep ){
                  cout<<"In TrackReco::computeTripleCombination"<<endl;
          cout<<"-------- minHStep "<<minHStep<<" maxHStep "<<maxHStep<<endl;
          int firstXbin=0, lastXBin=_nBinHelpStep;
          if ( minHStep>_minHelStep) {
                  for (int ibin=0; ibin<_nBinHelpStep; ++ibin){
                          if (minHStep<=_iXHelStep[ibin]) {
                                  firstXbin=ibin;
                                  break;
                          }
                  }
          }
          if ( maxHStep<_maxHelStep ) {
                  lastXBin=firstXbin;
                  for (int ibin=firstXbin; ibin<_nBinHelpStep; ++ibin){
                          if (maxHStep<=_iXHelStep[ibin]) {
                                  break;
                          }
                          ++lastXBin;
                  }
          }
          cout<<"-------- firstXbin "<<firstXbin<<" lastXBin "<<lastXBin<<endl;
          //firstXbin=0;

          double DeltaZ21, DeltaZ31, DeltaZ1;
          double DeltaPhi21, cosDPhi21_1, sinDPhi21, DeltaPhi31, cosDPhi31_1, sinDPhi31;
          double Phi1, tanPhi1;
          double X2X1, X3X1, XsRatio;
          double tmpVal;
          double turnsCut = 5.0*CLHEP::twopi;
          for ( std::map< int, std::multimap<int, int> >::iterator hitsTripleCouplings_it = _hitsTripleCouplings.begin();
                  hitsTripleCouplings_it != _hitsTripleCouplings.end(); ++hitsTripleCouplings_it ) 
          {
                  HitData *first  = strdat->find(hitsTripleCouplings_it->first)->second;
                  DeltaZ1 = first->_midp.getZ() - zo;
                  for ( std::multimap<int, int>::iterator secondHitsPairs_it = hitsTripleCouplings_it->second.begin();
                          secondHitsPairs_it != hitsTripleCouplings_it->second.end(); ++secondHitsPairs_it )
                  {
                          HitData *second = strdat->find(secondHitsPairs_it->first)->second;
                          HitData *third  = strdat->find(secondHitsPairs_it->second)->second;
                          DeltaZ21 = second->_midp.getZ() - first->_midp.getZ();
                          DeltaZ31 = third->_midp.getZ() - first->_midp.getZ();
                          X2X1 = second->_radius - first->_radius;
                          X3X1 = third->_radius - first->_radius;
                          XsRatio = X3X1/X2X1;
                          for ( int ibin=firstXbin; ibin<lastXBin; ++ibin ) {
                                  tmpVal      = CLHEP::twopi/_iXHelStep[ibin];
                                  DeltaPhi21  = tmpVal*DeltaZ21;
                                  DeltaPhi31  = tmpVal*DeltaZ31;
                                  Phi1        = tmpVal*DeltaZ1;
                                  if (Phi1>turnsCut || DeltaPhi21>turnsCut|| DeltaPhi31>turnsCut) continue;
                                  cosDPhi21_1 = cos(DeltaPhi21) - 1.0;
                                  cosDPhi31_1 = cos(DeltaPhi31) - 1.0;
                                  sinDPhi21   = sin(DeltaPhi21);
                                  sinDPhi31   = sin(DeltaPhi31);
                                  tanPhi1     = tan(Phi1);

                                  tmpVal = XsRatio*( first->_direct.getY()*cosDPhi21_1 -first->_direct.getX()*sinDPhi21 );
                                  tmpVal -= ( first->_direct.getY()*cosDPhi31_1 -first->_direct.getX()*sinDPhi31 );
                                  tmpVal /= ( XsRatio*( first->_direct.getX()*cosDPhi21_1 +first->_direct.getY()*sinDPhi21 )
                                                  -( first->_direct.getX()*cosDPhi31_1 +first->_direct.getY()*sinDPhi31 ) );

                                  _hPhi0HelStep->Fill( _iXHelStep[ibin], /*atan2( tmpVal-tanPhi1, 1.0+tmpVal*tanPhi1 )*/ atan( (tmpVal-tanPhi1)/(1.0+tmpVal*tanPhi1) ) );
                          }
                  }
          }
  }

  void TrackReco::measureRbyTripleCombination( double HStep, double HPhi0 ){
          double DeltaZ21, DeltaZ31, DeltaZ1;
          double DeltaPhi21, cosDPhi21_1, sinDPhi21, DeltaPhi31, cosDPhi31_1, sinDPhi31;
          double Phi1, tanPhi1;
          double X2X1, X3X1, XsRatio;
          double tmpR21, tmpR31;
          double tmpVal = CLHEP::twopi/HStep;
          double cosPhi01, sinPhi01;
          double tmpAngle, HPhi0compl;
          HPhi0compl=HPhi0+CLHEP::pi;

          for ( std::map< int, std::multimap<int, int> >::iterator hitsTripleCouplings_it = _hitsTripleCouplings.begin();
                  hitsTripleCouplings_it != _hitsTripleCouplings.end(); ++hitsTripleCouplings_it )
          {
                  HitData *first  = strdat->find(hitsTripleCouplings_it->first)->second;
                  DeltaZ1 = first->_midp.getZ() - zo;
                  for ( std::multimap<int, int>::iterator secondHitsPairs_it = hitsTripleCouplings_it->second.begin();
                          secondHitsPairs_it != hitsTripleCouplings_it->second.end(); ++secondHitsPairs_it )
                  {
                          HitData *second = strdat->find(secondHitsPairs_it->first)->second;
                          HitData *third  = strdat->find(secondHitsPairs_it->second)->second;
                          DeltaZ21 = second->_midp.getZ() - first->_midp.getZ();
                          DeltaZ31 = third->_midp.getZ() - first->_midp.getZ();
                          X2X1 = second->_radius - first->_radius;
                          X3X1 = third->_radius - first->_radius;
                          DeltaPhi21  = tmpVal*DeltaZ21;
                          DeltaPhi31  = tmpVal*DeltaZ31;
                          Phi1        = tmpVal*DeltaZ1;
                          cosDPhi21_1 = cos(DeltaPhi21) - 1.0;
                          cosDPhi31_1 = cos(DeltaPhi31) - 1.0;
                          sinDPhi21   = sin(DeltaPhi21);
                          sinDPhi31   = sin(DeltaPhi31);

                          tmpAngle=HPhi0+Phi1;
                          cosPhi01=cos(tmpAngle);
                          sinPhi01=sin(tmpAngle);
                          tmpR21 = X2X1/( first->_direct.getY()*( cosDPhi21_1*cosPhi01 - sinDPhi21*sinPhi01 ) -first->_direct.getX()*( cosDPhi21_1*sinPhi01 + sinDPhi21*cosPhi01 ) );
                          tmpR31 = X3X1/( first->_direct.getY()*( cosDPhi31_1*cosPhi01 - sinDPhi31*sinPhi01 ) -first->_direct.getX()*( cosDPhi31_1*sinPhi01 + sinDPhi31*cosPhi01 ) );
                          _hR_1->Fill(tmpR21);
                          _hR_1->Fill(tmpR31);
                          tmpAngle=HPhi0compl+Phi1;
                          cosPhi01=cos(tmpAngle);
                          sinPhi01=sin(tmpAngle);
                          tmpR21 = X2X1/( first->_direct.getY()*( cosDPhi21_1*cosPhi01 - sinDPhi21*sinPhi01 ) -first->_direct.getX()*( cosDPhi21_1*sinPhi01 + sinDPhi21*cosPhi01 ) );
                          tmpR31 = X3X1/( first->_direct.getY()*( cosDPhi31_1*cosPhi01 - sinDPhi31*sinPhi01 ) -first->_direct.getX()*( cosDPhi31_1*sinPhi01 + sinDPhi31*cosPhi01 ) );
                          _hR_1->Fill(tmpR21);
                          _hR_1->Fill(tmpR31);
                  }
          }
  }

  int TrackReco::butterflyFilter( float *in, float *out, int &nBinX, int &nBinY, float minCountCut ){
          // the 3x3 mask for this application is:
          //    0 -3  0
          //    2  2  2
          //    0 -3  0

          float maskY = -3.0;
          // the mask on X is {2.0, 2.0, 2.0} but for optimization is better to directly sum the 3 bins and after multiply by 2

          int nOutEntries=0;
          int tmpRow, tmpRowUp, tmpRowDown, tmpCol;
          int nBinXeff = nBinX+2;
          int nBinYeff = nBinY+2;

          float *buff = new float [nBinXeff*nBinYeff];
          for (int ib=0; ib<nBinXeff*nBinYeff; ib++) buff[ib]=0.00000;

          for (int j=1; j<(nBinY-1); j++){
                  tmpRowDown = j*nBinXeff;
                  tmpRow     = tmpRowDown + nBinXeff; //(j+1)*nBinXeff
                  tmpRowUp   = tmpRow + nBinXeff;     //(j+2)*nBinXeff
                  for (int i=0; i<nBinX; i++){
                          tmpCol=i+1;
                          //buff[tmpRow+tmpCol] = maskY[0]*in[tmpRowDown+tmpCol] /*+ maskY[1]*in[tmpRow+tmpCol]*/ + maskY[2]*in[tmpRowUp+tmpCol];  //the center element will be computed during the sum on rows
                          buff[tmpRow+tmpCol] = maskY*(in[tmpRowDown+tmpCol] + in[tmpRowUp+tmpCol]);  //the center element will be computed during the sum on rows
                  }
          }

          int tmpColR, tmpColL;
          for (int i=1; i<(nBinX-1); i++){
                  tmpColL = i;
                  tmpCol  = tmpColL+1;
                  tmpColR = tmpCol+1;
                  for (int j=1; j<(nBinY-1); j++){
                          tmpRow     = (j+1)*nBinXeff;
                          out[tmpRow+tmpCol] = 2.0*( in[tmpRow+tmpColL] + in[tmpRow+tmpCol] + in[tmpRow+tmpColR] ) + buff[tmpRow+tmpCol];
                          if (out[tmpRow+tmpCol]<minCountCut) out[tmpRow+tmpCol]=0.0;
                          else nOutEntries+=(int)out[tmpRow+tmpCol];
                  }
          }

          delete buff;
          return nOutEntries;
  }

  int TrackReco::butterflyFilterRot( float *in, float *out, int &nBinX, int &nBinY, float minCountCut ){
          // the 3x3 mask for this application is:
          //    0  2  0
          //   -3  2 -3
          //    0  2  0

          float maskX = -3.0;
          // the mask on Y is {2.0, 2.0, 2.0} but for optimization is better to directly sum the 3 bins and after multiply by 2

          int nOutEntries=0;
          int tmpRow, tmpRowUp, tmpRowDown, tmpCol;
          int nBinXeff = nBinX+2;
          int nBinYeff = nBinY+2;

          float *buff = new float [nBinXeff*nBinYeff];
          for (int ib=0; ib<nBinXeff*nBinYeff; ib++) buff[ib]=0.00000;

          for (int j=1; j<(nBinY-1); j++){
                  tmpRowDown = j*nBinXeff;
                  tmpRow     = tmpRowDown + nBinXeff; //(j+1)*nBinXeff
                  tmpRowUp   = tmpRow + nBinXeff;     //(j+2)*nBinXeff
                  for (int i=0; i<nBinX; i++){
                          tmpCol=i+1;
                          buff[tmpRow+tmpCol] = 2.0*(in[tmpRowDown+tmpCol] +in[tmpRow+tmpCol] + in[tmpRowUp+tmpCol]);  //the center element will be computed during the sum on rows
                  }
          }

          int tmpColR, tmpColL;
          for (int i=1; i<(nBinX-1); i++){
                  tmpColL = i;
                  tmpCol  = tmpColL+1;
                  tmpColR = tmpCol+1;
                  for (int j=1; j<(nBinY-1); j++){
                          tmpRow     = (j+1)*nBinXeff;
                          out[tmpRow+tmpCol] = maskX*(in[tmpRow+tmpColL] + in[tmpRow+tmpColR] ) + buff[tmpRow+tmpCol];
                          if (out[tmpRow+tmpCol]<minCountCut) out[tmpRow+tmpCol]=0.0;
                          else nOutEntries+=(int)out[tmpRow+tmpCol];
                  }
          }

          delete buff;
          return nOutEntries;
  }

  int TrackReco::butterflyFilterRot45( float *in, float *out, int &nBinX, int &nBinY, float minCountCut ){
          // the 3x3 mask for this application is:
          //   -3  0  2
          //    0  2  0
          //    2  0 -3
          ////    2  0 -3
          ////    0  2  0
          ////   -3  0  2

          float maskD = -3.0;
          // the mask on Y is {2.0, 2.0, 2.0} but for optimization is better to directly sum the 3 bins and after multiply by 2

          int nOutEntries=0;
          int tmpRow, tmpRowUp, tmpRowDown, tmpCol;
          int nBinXeff = nBinX+2;
          int nBinYeff = nBinY+2;
          int tmpColR, tmpColL;

          float *buff = new float [nBinXeff*nBinYeff];
          for (int ib=0; ib<nBinXeff*nBinYeff; ib++) buff[ib]=0.00000;

          for (int j=1; j<(nBinY-1); j++){
                  tmpRowDown = j*nBinXeff;
                  tmpRow     = tmpRowDown + nBinXeff; //(j+1)*nBinXeff
                  tmpRowUp   = tmpRow + nBinXeff;     //(j+2)*nBinXeff
                  for (int i=0; i<nBinX; i++){
                          tmpColL = i;
                          tmpCol  = tmpColL+1;
                          tmpColR = tmpCol+1;
                          buff[tmpRow+tmpCol] = 2.0*(in[tmpRowDown+tmpColL] +in[tmpRow+tmpCol] + in[tmpRowUp+tmpColR]);  //the center element will be computed during the sum on rows
                          //buff[tmpRow+tmpCol] = 2.0*(in[tmpRowDown+tmpColR] +in[tmpRow+tmpCol] + in[tmpRowUp+tmpColL]);  //the center element will be computed during the sum on rows
                  }
          }

          for (int i=1; i<(nBinX-1); i++){
                  tmpColL = i;
                  tmpCol  = tmpColL+1;
                  tmpColR = tmpCol+1;
                  for (int j=1; j<(nBinY-1); j++){
                          tmpRowDown = j*nBinXeff;
                          tmpRow     = tmpRowDown + nBinXeff; //(j+1)*nBinXeff
                          tmpRowUp   = tmpRow + nBinXeff;     //(j+2)*nBinXeff
                          out[tmpRow+tmpCol] = maskD*(in[tmpRowUp+tmpColL] + in[tmpRowDown+tmpColR] ) + buff[tmpRow+tmpCol];
                          //out[tmpRow+tmpCol] = maskD*(in[tmpRowUp+tmpColR] + in[tmpRowDown+tmpColL] ) + buff[tmpRow+tmpCol];
                          if (out[tmpRow+tmpCol]<minCountCut) out[tmpRow+tmpCol]=0.0;
                          else nOutEntries+=(int)out[tmpRow+tmpCol];
                  }
          }

          delete buff;
          return nOutEntries;
  }



}  // end namespace mu2e

using mu2e::TrackReco;
DEFINE_ART_MODULE(TrackReco);
