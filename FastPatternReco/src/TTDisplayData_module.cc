//
// this is a old version, visualization and embedded implementation of the bck rejection algorithm
//
// $Id: TTDisplayData_module.cc,v 1.1 2011/06/23 21:56:11 tassiell Exp $
// $Author: tassiell $
// $Date: 2011/06/23 21:56:11 $
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
#include "art/Framework/Core/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/Handle.h"
#include "art/Persistency/Provenance/Provenance.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


#include "CLHEP/Vector/TwoVector.h"

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
#include "FastPatternReco/inc/TTHitPerTrackData.hh"
#include "FastPatternReco/inc/GenTrackData.hh"

// Root includes.
#include "TApplication.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TClonesArray.h"
#include "TObjArray.h"
#include "TLine.h"
#include "TArrow.h"
#include "TEllipse.h"
#include "TMarker.h"
#include "TBox.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TF1.h"
#include "TSpectrum.h"
#include "TLatex.h"

using namespace std;

typedef std::multimap<unsigned int, size_t, less<unsigned int> > stbrel;                                        //"Straw-histogram bin" relation (navigation)

typedef std::multimap<unsigned short, unsigned short, less<unsigned short> > ptcsctrrel;                        //"Station Pitch-Station" relation
typedef std::map<unsigned short, ptcsctrrel> ptcbnrel;                                                          //container for each row "Station Pitch-Station" relation
typedef std::set<unsigned short, less<unsigned short> > goodsttnrel;                                                               //"good Station" that has at least one vote in both scan direction per row
typedef std::map<unsigned short, goodsttnrel> goodsctrstatnrel;                                                 //container for each row "good Station"
typedef std::multimap<unsigned short, std::pair<size_t, unsigned short>, less<unsigned short> > ptcclsbrel;     //"Station Pitch-row Cluster" relation ("row Cluster" = cluster along each row in the Station-Sector map)
typedef std::multimap<unsigned short, ptcclsbrel, less<unsigned short> > ptcmaprowclsrel;                       //container of the "Station Pitch-row Cluster" relation for each sector (row should be see like an alias for sector)
typedef std::pair<size_t, size_t> rwclclcpl;                                                                    //"row Cluster-row Cluster" coupling along a row
typedef std::map<rwclclcpl, std::vector<unsigned short> > avptcclscpl;                                          //list of all available station pitch for each "row Cluster-row Cluster" pair
typedef std::multimap<unsigned short, avptcclscpl , less<unsigned short> > rwavptcclscpl;                       //container of the data for "row Cluster-row Cluster" pair for each row

namespace mu2e {

  class Straw;

  struct rowClust{
          friend struct Clust;

          rowClust():_firstStationID(0),
                     _lastStationID(0),
                     _nHit(0),
                     _mean(0.00000),
                     _MS(0.00000),
                     _sigma(0.00000),
                     tmpData(0.00000),
                     tmpMean(0.00000),
                     tmpMS(0.00000)
          {}

          rowClust(unsigned short iSttID, unsigned short iHitMult):
                     _nHit(0),
                     _mean(0.00000),
                     _MS(0.00000),
                     _sigma(0.00000),
                     tmpData(0.00000),
                     tmpMean(0.00000),
                     tmpMS(0.00000)

          {
                  _firstStationID=iSttID;
                  addHit(iSttID, iHitMult);
          }
          unsigned short _firstStationID;
          unsigned short _lastStationID;
          unsigned int _nHit;
          float _mean;
          float _MS;
          float _sigma;

          void addHit(unsigned short iSttID, unsigned short iHitMult) {
                  _lastStationID=iSttID;
                  _nHit+=iHitMult;
                  tmpData=(float)iSttID + 0.50000;
                  tmpMean+= ((float)iHitMult)*tmpData;
                  tmpMS+= ((float)iHitMult)*pow(tmpData,2);
                  _mean=tmpMean/((float) _nHit);
                  _MS=tmpMS/((float) _nHit);
                  if (_lastStationID!=_firstStationID) _sigma = sqrt(_MS-pow(_mean,2));
          }

  protected:
          float tmpData;
          float tmpMean;
          float tmpMS;
  };

  typedef boost::shared_ptr<rowClust> rwClustPtr;
  typedef std::vector< rwClustPtr > rwclstvec;
  typedef std::map<unsigned short, rwclstvec> rwclinrwrel;                                                         //"row Clusters" for each row contained
  typedef std::multimap<unsigned short, rwClustPtr > rwclincl;                                                   //"row Clusters" for each row contained in Cluster


  struct Clust{
          Clust():_nHit(0),
                  _mean_Sttn(0.00000),
                  _MS_Sttn(0.00000),
                  _sigma_Sttn(0.00000),
                  _mean_Sctr(0.00000),
                  _MS_Sctr(0.00000),
                  _sigma_Sctr(0.00000),
                  _m(0.00000),
                  _q(0.00000),
                  _errm(0.00000),
                  _errq(0.00000),
                  _firstSectorID(0),
                  _lastSectorID(0),
                  _minStationID(0),
                  _maxStationID(0),
                  tmpMeanX(0.00000),
                  tmpMSX(0.00000),
                  tmpDataY(0.00000),
                  tmpMeanY(0.00000),
                  tmpMSY(0.00000),
                  tmpMeanXY(0.00000)
          {}
          Clust( unsigned short iRow, rwClustPtr &firstRwclust ):
                  _nHit(0),
                  _mean_Sttn(0.00000),
                  _MS_Sttn(0.00000),
                  _sigma_Sttn(0.00000),
                  _m(0.00000),
                  _q(0.00000),
                  _errm(0.00000),
                  _errq(0.00000),
                  _mean_Sctr(0.00000),
                  _MS_Sctr(0.00000),
                  _sigma_Sctr(0.00000),
                  tmpMeanX(0.00000),
                  tmpMSX(0.00000),
                  tmpDataY(0.00000),
                  tmpMeanY(0.00000),
                  tmpMSY(0.00000),
                  tmpMeanXY(0.00000)
          {
                  _rClusts.clear();
                  _rClusts.insert( rwclincl::value_type( iRow, firstRwclust ) );
                  _nHit+=firstRwclust->_nHit;
                  _firstSectorID=iRow;
                  _lastSectorID=_firstSectorID;
                  _minStationID=firstRwclust->_firstStationID;
                  _maxStationID=firstRwclust->_lastStationID;

                  tmpMeanX+=firstRwclust->tmpMean;
                  tmpMSX+=firstRwclust->tmpMS;
                  _mean_Sttn=tmpMeanX/((float) _nHit);
                  _MS_Sttn=tmpMSX/((float) _nHit);

                  tmpDataY=(float)iRow + 0.50000;
                  tmpMeanY+= ((float)firstRwclust->_nHit)*tmpDataY;
                  tmpMSY+= ((float)firstRwclust->_nHit)*pow(tmpDataY,2);
                  tmpMeanXY+=firstRwclust->tmpMean*tmpDataY;
                  _mean_Sctr=tmpDataY;//tmpMeanY/((float) _nHit);
                  _MS_Sctr=tmpMSY/((float) _nHit);
                  //if (_lastSectorID!=_firstSectorID) _sigma_Sctr = sqrt(_MS_Sctr-pow(_mean_Sctr,2));
                  if (_maxStationID!=_minStationID) {
                          _sigma_Sttn = sqrt(_MS_Sttn-pow(_mean_Sttn,2));
//                          _q=tmpDataY;
                  }
//                  else {
//                          _m=std::numeric_limits<float>::infinity();
//                          _q=(float)_maxStationID + 0.50000;
//                  }

                  _q=tmpDataY;


          }
          void addRwClust( unsigned short iRow, rwClustPtr &iRwclust ) {
                  _rClusts.insert( rwclincl::value_type( iRow, iRwclust ) );
                  _nHit+=iRwclust->_nHit;
                  if ( iRow<_firstSectorID ) _firstSectorID=iRow;
                  if ( iRow>_lastSectorID ) _lastSectorID=iRow;
                  if ( iRwclust->_firstStationID<_minStationID ) _minStationID=iRwclust->_firstStationID;
                  if ( iRwclust->_lastStationID>_maxStationID ) _maxStationID=iRwclust->_lastStationID;

                  tmpMeanX+=iRwclust->tmpMean;
                  tmpMSX+=iRwclust->tmpMS;
                  _mean_Sttn=tmpMeanX/((float) _nHit);
                  _MS_Sttn=tmpMSX/((float) _nHit);
                  if (_maxStationID!=_minStationID) _sigma_Sttn = sqrt(_MS_Sttn-pow(_mean_Sttn,2));
                  else {
                          _m=std::numeric_limits<float>::infinity();
                          _q=(float)_maxStationID + 0.50000;
                  }

                  tmpDataY=(float)iRow + 0.50000;
                  tmpMeanY+= ((float)iRwclust->_nHit)*tmpDataY;
                  tmpMSY+= ((float)iRwclust->_nHit)*pow(tmpDataY,2);
                  tmpMeanXY+=iRwclust->tmpMean*tmpDataY;
                  _mean_Sctr=tmpMeanY/((float) _nHit);
                  _MS_Sctr=tmpMSY/((float) _nHit);
                  if (_lastSectorID!=_firstSectorID) {
                          _sigma_Sctr = sqrt(_MS_Sctr-pow(_mean_Sctr,2));

                          if ( _maxStationID!=_minStationID ) {
                                  tmpNsCovX=_nHit*tmpMSX-pow(tmpMeanX,2);
                                  tmpNsCovY=_nHit*tmpMSY-pow(tmpMeanY,2);
                                  _m=(_nHit*tmpMeanXY-tmpMeanX*tmpMeanY)/tmpNsCovX;
                                  _q=(tmpMSX*tmpMeanY-tmpMeanX*tmpMeanXY)/tmpNsCovX;
                                  if (_nHit>2) {
                                          tmpSigmaY=sqrt( tmpNsCovY/(_nHit*(_nHit-2)) );
                                          _errm=tmpSigmaY/sqrt( tmpNsCovX/_nHit );
                                          _errq=tmpSigmaY*sqrt( 1.00000/_nHit+pow(tmpMeanX,2)/(_nHit*tmpNsCovX) );
                                  }
                          }
                  }
          }
          rwclincl _rClusts;
          unsigned int _nHit;
          float _mean_Sttn;
          float _MS_Sttn;
          float _sigma_Sttn;
          float _mean_Sctr;
          float _MS_Sctr;
          float _sigma_Sctr;
          float _m;
          float _q;
          float _errm;
          float _errq;
          unsigned short _firstSectorID;
          unsigned short _lastSectorID;
          unsigned short _minStationID;
          unsigned short _maxStationID;

  private:
          float tmpMeanX;
          float tmpMSX;
          float tmpNsCovX;
          float tmpDataY;
          float tmpMeanY;
          float tmpMSY;
          float tmpNsCovY;
          float tmpSigmaY;
          float tmpMeanXY;
  };




  class TTDisplayData : public art::EDAnalyzer {
  public:

    explicit TTDisplayData(fhicl::ParameterSet const& pset);
    virtual ~TTDisplayData() {
            _peakFinder->Delete();
            _fg->Delete();
            //if (_peakFinder)        delete _peakFinder;
            //if (_fg)                delete _fg;
            if (_fakeCanvas)        delete _fakeCanvas;
            if (_peaksCanvases)     delete _peaksCanvases;
            if (_peaksCanvHistos)   delete _peaksCanvHistos;
            if (_hPkStDistTrs)      delete _hPkStDistTrs;
            if (_hPkStDistanceTrs)  delete _hPkStDistanceTrs;
            if (_hPkStDist2W)       delete _hPkStDist2W;
            if (_hPkStClustDist2W)  delete _hPkStClustDist2W;
            if (_hPkSecDist2W)      delete _hPkSecDist2W;
            if (_hPkSecClustDist2W) delete _hPkSecClustDist2W;
            if (_hPkSecVsStationDist2W)      delete _hPkSecVsStationDist2W;
            if (_hPkSecVsStationDist2WGood)      delete _hPkSecVsStationDist2WGood;
            if (_hPkAccArrSecStation)      delete _hPkAccArrSecStation;
    }

    virtual void beginJob();
    void endJob();

    // This is called for each event.
    void analyze(art::Event const& e);

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

    // Label of the generator.
    std::string _generatorModuleLabel;

    // drift vevocity
    double _driftVelocity;

//    // Number of events to accumulate between prompts.
//    int _nAccumulate;

    // End: run time parameters

    // Pointers to histograms, ntuples, TGraphs.
    TH2F*         _hHitTransverse;
    TH2F*         _hHitLongit;
    TH2F*         _hHitStereoPl;
    TH2F*         _hHitStereoMi;
    TH1F*         _hHitZ;
    TH1I*         _hHitNloop;
    TH1F*         _hMCHitDTime;
    TH1F*         _hHitTime;
    TH1F*         _hHitClustTime;
    TH1F*         _hSelHitTime;
    TH1F*         _hSelHitClustTime;
    TH1F*         _hSelHitClustTimeSigma;
    TCanvas*      _canvasTLview;
    TCanvas*      _canvasStLview;
    TCanvas*      _canvasPl;
    TObjArray*    _peaksCanvases;
    TObjArray*    _peaksCanvHistos;
    TObjArray*    _hPkStDistTrs;
    TObjArray*    _hPkStDistanceTrs;
    TObjArray*    _hPkStDist2W;
    TObjArray*    _hPkStClustDist2W;
    TObjArray*    _hPkSecDist2W;
    TObjArray*    _hPkSecClustDist2W;
    TObjArray*    _hPkSecVsStationDist2W;
    TObjArray*    _hPkSecVsStationDist2WGood;
    TObjArray*    _hPkAccArrSecStation;

    TF1*          _fg;
    TCanvas*      _cnvForPeakstudy;

    TCanvas*      _fakeCanvas;

    TSpectrum *   _peakFinder;

    int   ntimeBin;
    float maxTimeHist; //ns
    float timeBinDim;  //ns

    void fillVotingArray(TH2I *startDist, TH2I *votingArray, int minPitch=3, int maxPitch=11);
    rwclinrwrel rwClst_forRw_rel;
    std::vector<Clust> clustersList;
    void findCluster(unsigned short tmpRowId, rwClustPtr &startingRowClust);

    // Some ugly but necessary ROOT related bookkeeping:

    // The job needs exactly one instance of TApplication.  See note 1.
    auto_ptr<TApplication> _application;

    // Save directory from beginJob so that we can go there in endJob. See note 3.
    TDirectory* _directory;

  };

  TTDisplayData::TTDisplayData(fhicl::ParameterSet const& pset) :

    // Run time parameters
    _moduleLabel(pset.get<string>("module_label")),/*@module_label*/
    _g4ModuleLabel(pset.get<string>("g4ModuleLabel")),
    _trackerStepPoints(pset.get<string>("trackerStepPoints","tracker")),
    _makerModuleLabel(pset.get<string>("makerModuleLabel")),
    _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
   /*_nAccumulate(pset.get<int>("nAccumulate",20)),*/
    _driftVelocity(pset.get<double>("driftVelocity",0.05)),   // mm/ns

    // ROOT objects that are the main focus of this example.
    _hHitTransverse(0),
    _hHitLongit(0),
    _hHitStereoPl(0),
    _hHitStereoMi(0),
    _hHitZ(0),
    _hHitNloop(0),
    _hMCHitDTime(0),
    _hHitTime(0),
    _hHitClustTime(0),
    _hSelHitTime(0),
    _hSelHitClustTime(0),
    _hSelHitClustTimeSigma(0),
    _canvasTLview(0),
    _canvasStLview(0),
    _canvasPl(0),
    _peaksCanvases(0),
    _peaksCanvHistos(0),
    _hPkStDistTrs(0),
    _hPkStDistanceTrs(0),
    _hPkStDist2W(0),
    _hPkStClustDist2W(0),
    _hPkSecDist2W(0),
    _hPkSecClustDist2W(0),
    _hPkSecVsStationDist2W(0),
    _hPkSecVsStationDist2WGood(0),

    _fg(0),
    _cnvForPeakstudy(0),

    _fakeCanvas(0),

    _peakFinder(0),
    ntimeBin(0),
    maxTimeHist(2500.0),
    timeBinDim(10.0),

    // Some ugly but necessary ROOT related bookkeeping.
    _application(0),
    _directory(0){

 }

  void TTDisplayData::beginJob(){

	  cout<<"Starting DisplayData jos!"<<endl;

    // Get access to the TFile service and save current directory for later use.
    art::ServiceHandle<art::TFileService> tfs;

    ntimeBin    = (int) maxTimeHist/timeBinDim;

    // Create a histogram.
    _hHitTransverse    = tfs->make<TH2F>( "hHitTransverse",  "Hits per Event in transverse view", 1500, -75.0, 75.0, 1500, -75.0, 75.0  );
    _hHitLongit        = tfs->make<TH2F>( "hHitLongit", "Hits per Event in Z view ", 3200, -160.0, 160.0, 1500, -75.0, 75.0  );
    _hHitStereoPl      = tfs->make<TH2F>( "hHitStereoPl",  "Stereo+ Hits per Event at z=0", 1500, -75.0, 75.0, 1500, -75.0, 75.0  );
    _hHitStereoMi      = tfs->make<TH2F>( "hHitStereoMi",  "Stereo- Hits per Event at z=0", 1500, -75.0, 75.0, 1500, -75.0, 75.0  );
    _hHitZ             = tfs->make<TH1F>( "hHitZ",         "Z of the Hits per Event", 1280, -160.0, 160.0  );
    _hHitNloop         = tfs->make<TH1I>( "hHitNloop",     "Track number of loops per Event", 5, 0, 5  );
    _hMCHitDTime       = tfs->make<TH1F>( "hMCHitDTime",   "Delta Time of the MC Step hits per Event", 200, 0.0, 50.0  );
    _hHitTime          = tfs->make<TH1F>( "hHitTime",      "Time of the Hits per Event", ntimeBin, 0.0, maxTimeHist  );
    _hHitClustTime     = tfs->make<TH1F>( "hHitClustTime", "Cluster of Time of the Hits per Event", ntimeBin, 0.0, maxTimeHist  );
    _hSelHitTime       = tfs->make<TH1F>( "hSelHitTime",      "Time of the Signal e^{-} Hits per Event", ntimeBin, 0.0, maxTimeHist  );
    _hSelHitClustTime  = tfs->make<TH1F>( "hSelHitClustTime", "Cluster of Time of the Signal e^{-} Hits per Event", ntimeBin, 0.0, maxTimeHist  );
    _hSelHitClustTimeSigma  = tfs->make<TH1F>( "hSelHitClustTimeSigma", "#sigma of the Cluster of Time of the Signal e^{-} Hits per Event", 200, 0.0, 100.0  );

    _hHitTransverse   ->SetXTitle("cm");
    _hHitTransverse   ->SetYTitle("cm");
    _hHitLongit       ->SetXTitle("cm");
    _hHitLongit       ->SetYTitle("cm");
    _hHitStereoPl     ->SetXTitle("cm");
    _hHitStereoPl     ->SetYTitle("cm");
    _hHitStereoMi     ->SetXTitle("cm");
    _hHitStereoMi     ->SetYTitle("cm");
    _hHitZ            ->SetXTitle("cm");
    _hHitNloop        ->SetXTitle("N loops");
    _hMCHitDTime      ->SetXTitle("ns");
    _hHitTime         ->SetXTitle("ns");
    _hHitClustTime    ->SetXTitle("ns");
    _hSelHitTime      ->SetXTitle("ns");
    _hSelHitClustTime ->SetXTitle("ns");
    _hSelHitClustTimeSigma ->SetXTitle("ns");

    _hSelHitTime      ->SetFillColor(kRed);
    _hSelHitTime      ->SetFillStyle(3001);
    _hSelHitClustTime ->SetFillColor(kRed);
    _hSelHitClustTime ->SetFillStyle(3001);

    // If needed, create the ROOT interactive environment. See note 1.
    if ( !gApplication ){
      int    tmp_argc(0);
      char** tmp_argv(0);
      _application = auto_ptr<TApplication>(new TApplication( "noapplication", &tmp_argc, tmp_argv ));
    }

    gStyle->SetPalette(1);
    gROOT->SetStyle("Plain");

    _peakFinder = new TSpectrum(20);

    _fg = new TF1("fg","gaus");
    _cnvForPeakstudy = tfs->make<TCanvas>("cnvForPeakstudy","Peaks studies container");

    // Create a canvas with a unique name.  See note 2.
    TString name  = "canvasTLview_"     + _moduleLabel;
    TString title = "Canvas Transverse and Longitudinal view for " + _moduleLabel;
    int window_size(860);
    _canvasTLview = tfs->make<TCanvas>(name,title,2*window_size,window_size);
    _canvasTLview->Divide(2,1);

    name  = "canvasLview_"     + _moduleLabel;
    title = "Canvas Longitudinal view for " + _moduleLabel;
    _canvasStLview = tfs->make<TCanvas>(name,title,window_size,window_size);

//    name  = "canvas_"     + _moduleLabel;
//    title = "Canvas for " + _moduleLabel;
//    _canvas = tfs->make<TCanvas>(name,title,window_size,window_size);
//    _canvas->Divide(2,2);
    _peaksCanvases     = new TObjArray();
    _peaksCanvHistos   = new TObjArray();
    _hPkStDistTrs      = new TObjArray();
    _hPkStDistanceTrs  = new TObjArray();
    _hPkStDist2W       = new TObjArray();
    _hPkStClustDist2W  = new TObjArray();
    _hPkSecDist2W      = new TObjArray();
    _hPkSecClustDist2W = new TObjArray();
    _hPkSecVsStationDist2W      = new TObjArray();
    _hPkSecVsStationDist2WGood      = new TObjArray();
    _hPkAccArrSecStation      = new TObjArray();

    name  = "canvasPl_"     + _moduleLabel;
    title = "Canvas for Plots " + _moduleLabel;
    _canvasPl = tfs->make<TCanvas>(name,title,window_size,window_size);
    _canvasPl->Divide(2,2);

    _fakeCanvas = new TCanvas("canvas_Fake","double click for next event",300,100);
//    _fakeCanvas = tfs->make<TCanvas>("canvas_Fake","double click for next event",400,200);

    // Draw the still empty histogram. It will be updated later.
    _canvasTLview->cd(1);
    _hHitTransverse->SetStats(kFALSE);
    _hHitTransverse->Draw();
//    _canvasStLview->cd();
    _canvasTLview->cd(2);
    _hHitLongit->SetStats(kFALSE);
    _hHitLongit->Draw();
//    _canvas->cd(1);
    _hHitTransverse->SetStats(kFALSE);
//    _hHitTransverse->Draw();
//    _canvas->cd(2);
    _hHitLongit->SetStats(kFALSE);
//    _hHitLongit->Draw();
//    _canvas->cd(3);
    _hHitStereoMi->SetStats(kFALSE);
//    _hHitStereoMi->Draw();
//    _canvas->cd(4);
    _hHitStereoPl->SetStats(kFALSE);
//    _hHitStereoPl->Draw();

    _canvasPl->cd(1);
    _hHitTime->SetStats(kTRUE);
    _hHitTime->Draw();
    _canvasPl->cd(2);
    _hHitClustTime->SetStats(kFALSE);
    _hHitClustTime->Draw();
    _canvasPl->cd(3);
    _hHitZ->SetStats(kFALSE);
    _hHitZ->Draw();
    _canvasPl->cd(4);
    _hHitNloop->SetStats(kFALSE);
    _hHitNloop->Draw();
//    _hMCHitDTime->SetStats(kFALSE);
//    _hMCHitDTime->Draw();

    _hSelHitTime->SetStats(kFALSE);
    _hSelHitClustTime->SetStats(kFALSE);
    _hSelHitClustTimeSigma->SetStats(kTRUE);

    // See note 3.
    _directory = gDirectory;

  }

  void TTDisplayData::analyze(art::Event const& event ) {

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

    _hHitTransverse->Reset();
    _hHitLongit->Reset();
    _hHitStereoMi->Reset();
    _hHitStereoPl->Reset();
    _hMCHitDTime->Reset();
    _hHitTime->Reset();
    _hHitClustTime->Reset();
    _hSelHitTime->Reset();
    _hSelHitClustTime->Reset();
    _hHitZ->Reset();
    _hHitNloop->Reset();

//    _hHitClustTime->SetMaximum(200.0);

    _peakFinder->Clear();
    //_peaksCanvases->Clear();
    typedef boost::shared_ptr<TArrow> clssegdr;
    std::vector< std::vector< clssegdr > > clstSegments;
    TF1 clustSegline("clustSegline","pol1", 0.0, 20.0);

    const Tracker& tracker = getTrackerOrThrow();
    const TTracker &ttr = static_cast<const TTracker&>( tracker );
    const std::vector<Device> ttrdev = ttr.getDevices();

    float intTimeWind = ttr.strawRadius()/_driftVelocity + 20.0;   //max drift time + max TOF of a signal electron
    stbrel timeBin_Straw_rel;
    int tmpiTimeBin;

    int nDevice = ttr.nDevices();
    double deviceHalfLength = ttr.getDeviceEnvelopeParams().zHalfLength();
    double rIn  = ttr.getTrackerEnvelopeParams().innerRadius(); //380.00;
    double rOut = ttr.getSupportParams().innerRadius(); //700.00;
    rIn/=CLHEP::cm;
    rOut/=CLHEP::cm;

    stbrel ssmap_Bin_Straw_rel;  //relation between straw and the AbsSector (Rotation) vs Station mao
    int tmpiGeomBin;

    art::Handle<StrawHitCollection> pdataHandle;
    event.getByLabel(_makerModuleLabel,pdataHandle);
    StrawHitCollection const* hits = pdataHandle.product();

    // Get the persistent data about the StrawHitsMCTruth.
    art::Handle<StrawHitMCTruthCollection> truthHandle;
    event.getByLabel(_makerModuleLabel,truthHandle);
    StrawHitMCTruthCollection const* hits_truth = truthHandle.product();

//    // Get the persistent data about pointers to StepPointMCs
//    art::Handle<DPIndexVectorCollection> mcptrHandle;
//    event.getByLabel(_makerModuleLabel,"StrawHitMCPtr",mcptrHandle);
//    DPIndexVectorCollection const* hits_mcptr = mcptrHandle.product();
//
//    // Get the persistent data about the StepPointMCs. More correct implementation
//    // should look for product ids in DPIndexVectorCollection, rather than
//    // use producer name directly ("g4run").
//
//    art::Handle<StepPointMCCollection> mchitsHandle;
//    event.getByLabel(_g4ModuleLabel,_trackerStepPoints,mchitsHandle);
//    StepPointMCCollection const* mchits = mchitsHandle.product();

    // Get the persistent data about pointers to StepPointMCs
    art::Handle<PtrStepPointMCVectorCollection> mcptrHandle;
    event.getByLabel(_makerModuleLabel,"StrawHitMCPtr",mcptrHandle);
    PtrStepPointMCVectorCollection const* hits_mcptr = mcptrHandle.product();

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

    art::Handle<SimParticleCollection> simParticles;
    event.getByLabel(_g4ModuleLabel, simParticles);

//    // Handle to information about G4 physical volumes.
//    art::Handle<PhysicalVolumeInfoCollection> volumes;
//    event.getRun().getByLabel(_g4ModuleLabelvolumes);
//
//    // Some files might not have the SimParticle and volume information.
//    bool haveSimPart = ( simParticles.isValid() && volumes.isValid() );
//
//    // Other files might have empty collections.
//    if ( haveSimPart ){
//      haveSimPart = !(simParticles->empty() || volumes->empty());
//    }

    size_t nStrawPerEvent = hits->size();

    TClonesArray *hitDrawsTrans = new TClonesArray("TLine");
    hitDrawsTrans->ExpandCreateFast(nStrawPerEvent);
    TClonesArray *hitStrawDrawsTrans = new TClonesArray("TLine");
    hitStrawDrawsTrans->ExpandCreateFast(4*nStrawPerEvent);
    TClonesArray *hitDrawsTransPoint = new TClonesArray("TMarker");  // only for signal electron
    TClonesArray *hitDrawsLongPoint  = new TClonesArray("TMarker");  // only for signal electron
    int hDP_prevDim;
    TClonesArray *hitDrawsZ = new TClonesArray("TEllipse");
    hitDrawsZ->ExpandCreateFast(nStrawPerEvent);

    TClonesArray *hitDrawsStZprj = new TClonesArray("TEllipse");
    hitDrawsStZprj->ExpandCreateFast(nStrawPerEvent);

    //vector<int> idHitStereoMin;
    bool *isStereoMin = new bool [nStrawPerEvent];

    std::map<SimParticleCollection::key_type,std::pair<CLHEP::Hep3Vector,CLHEP::HepLorentzVector> > genElectrons;  //Generation data of electrons that come from target
    std::map<SimParticleCollection::key_type,std::pair<CLHEP::Hep3Vector,CLHEP::Hep3Vector> > firstHitOfGenEl;     //Position and momentum data of the first (in TOF) interaction hit (Geant4 hit) inside active volume
    std::map<SimParticleCollection::key_type,double > firstHitTime;                                                //Time (TOF) of the first interaction hit (Geant4 hit) inside active volume
    std::map<SimParticleCollection::key_type,std::pair<CLHEP::Hep3Vector,CLHEP::Hep3Vector> >::iterator firstHitOfGenEl_it;  //_it=iterator
    std::map<SimParticleCollection::key_type,double >::iterator firstHitTime_it;

    std::map<SimParticleCollection::key_type, GenTrackData<TTHitPerTrackData> > genTracks;                         //Data (postion etc.) of the interaction hit (Geant4 hit) of the electrons tracks (from target) that should give a signal
    std::map<SimParticleCollection::key_type, GenTrackData<TTHitPerTrackData> >::iterator genTracks_it;

    float aveZ=0.0;
    double mchittime=0.0;

    bool overlapped = false;
    bool isFirst = false;

    int il;
    int StrawColor = kCyan-6;
    int noiseColor = kYellow+3;

    for (size_t i=0; i<nStrawPerEvent; ++i) {

      // Access data
      StrawHit        const&      hit(hits->at(i));
      StrawHitMCTruth const&    truth(hits_truth->at(i));
//      DPIndexVector   const&    mcptr(hits_mcptr->at(i));
      PtrStepPointMCVector const&    mcptr(hits_mcptr->at(i));

      double hitEnergy = hit.energyDep();

      //Skip the straw if the energy of the hit is smaller than the minimum required
      //if (hitEnergy < _minimumEnergy) continue;

      //Get hit straw
      StrawIndex si = hit.strawIndex();
      const Straw & str = tracker.getStraw(si);

      // cout << "Getting informations about cells" << endl;

      StrawId sid = str.id();
      int stn     = sid.getStraw();
      int layern  = sid.getLayer();
      int devicen = sid.getDevice();
      int sectorn = sid.getSector();

      const CLHEP::Hep3Vector stMidPoint3 = str.getMidPoint();

      //time of the hit
      double hitTime = hit.time();

      //direction of the straw
      const CLHEP::Hep3Vector stDirection3 = str.getDirection();

      // cout << "Reading MCtruth info" << endl;

      // Get MC truth data
      double driftTime = truth.driftTime();
      double driftDistance = truth.driftDistance();

      //Position along the wire using mctruth info
      double vMC = truth.distanceToMid();

      const CLHEP::Hep3Vector MCHitPoint = stMidPoint3 + (vMC/stDirection3.mag())*stDirection3;

      CLHEP::Hep2Vector StYNormalDir;
      StYNormalDir.set(stDirection3.getX(),stDirection3.getY());
      StYNormalDir.rotate(90.0*CLHEP::degree);
      StYNormalDir=StYNormalDir.unit();
//      cout<<"Normal dir "<<StYNormalDir<<endl;
      CLHEP::Hep3Vector tmpUprPoint  = driftDistance*StYNormalDir+MCHitPoint;
      CLHEP::Hep3Vector tmpDownPoint = -driftDistance*StYNormalDir+MCHitPoint;

//      cout<<"Point "<<MCHitPoint<<endl;
//      cout<<"Up Point  "<<tmpUprPoint<<endl;
//      cout<<"Dow Point "<<tmpDownPoint<<endl;

      ((TLine *) hitDrawsTrans->At(i))->SetX1(tmpUprPoint.getX()/CLHEP::cm);
      ((TLine *) hitDrawsTrans->At(i))->SetY1(tmpUprPoint.getY()/CLHEP::cm);
      ((TLine *) hitDrawsTrans->At(i))->SetX2(tmpDownPoint.getX()/CLHEP::cm);
      ((TLine *) hitDrawsTrans->At(i))->SetY2(tmpDownPoint.getY()/CLHEP::cm);
      ((TLine *) hitDrawsTrans->At(i))->SetLineWidth(2.0);
      ((TLine *) hitDrawsTrans->At(i))->SetLineColor(noiseColor);

      CLHEP::Hep2Vector StTPDir;
      StTPDir.set(stDirection3.getX(),stDirection3.getY());
      StTPDir=StTPDir.unit();
      CLHEP::Hep3Vector tmpStDUpPoint1  = stMidPoint3+StTPDir*str.getHalfLength()+str.getRadius()*StYNormalDir;
      CLHEP::Hep3Vector tmpStDUpPoint2  = stMidPoint3-StTPDir*str.getHalfLength()+str.getRadius()*StYNormalDir;
      CLHEP::Hep3Vector tmpStDDownPoint1 = stMidPoint3+StTPDir*str.getHalfLength()-str.getRadius()*StYNormalDir;
      CLHEP::Hep3Vector tmpStDDownPoint2 = stMidPoint3-StTPDir*str.getHalfLength()-str.getRadius()*StYNormalDir;
      il=4*i;
      ((TLine *) hitStrawDrawsTrans->At(il))->SetX1(tmpStDUpPoint1.getX()/CLHEP::cm);
      ((TLine *) hitStrawDrawsTrans->At(il))->SetY1(tmpStDUpPoint1.getY()/CLHEP::cm);
      ((TLine *) hitStrawDrawsTrans->At(il))->SetX2(tmpStDUpPoint2.getX()/CLHEP::cm);
      ((TLine *) hitStrawDrawsTrans->At(il))->SetY2(tmpStDUpPoint2.getY()/CLHEP::cm);
      ((TLine *) hitStrawDrawsTrans->At(il))->SetLineWidth(0.1);
      ((TLine *) hitStrawDrawsTrans->At(il))->SetLineStyle(3);
      ((TLine *) hitStrawDrawsTrans->At(il))->SetLineColor(StrawColor);
      ((TLine *) hitStrawDrawsTrans->At(il+1))->SetX1(tmpStDDownPoint1.getX()/CLHEP::cm);
      ((TLine *) hitStrawDrawsTrans->At(il+1))->SetY1(tmpStDDownPoint1.getY()/CLHEP::cm);
      ((TLine *) hitStrawDrawsTrans->At(il+1))->SetX2(tmpStDDownPoint2.getX()/CLHEP::cm);
      ((TLine *) hitStrawDrawsTrans->At(il+1))->SetY2(tmpStDDownPoint2.getY()/CLHEP::cm);
      ((TLine *) hitStrawDrawsTrans->At(il+1))->SetLineWidth(0.1);
      ((TLine *) hitStrawDrawsTrans->At(il+1))->SetLineStyle(3);
      ((TLine *) hitStrawDrawsTrans->At(il+1))->SetLineColor(StrawColor);
      ((TLine *) hitStrawDrawsTrans->At(il+2))->SetX1(tmpStDUpPoint1.getX()/CLHEP::cm);
      ((TLine *) hitStrawDrawsTrans->At(il+2))->SetY1(tmpStDUpPoint1.getY()/CLHEP::cm);
      ((TLine *) hitStrawDrawsTrans->At(il+2))->SetX2(tmpStDDownPoint1.getX()/CLHEP::cm);
      ((TLine *) hitStrawDrawsTrans->At(il+2))->SetY2(tmpStDDownPoint1.getY()/CLHEP::cm);
      ((TLine *) hitStrawDrawsTrans->At(il+2))->SetLineWidth(0.1);
      ((TLine *) hitStrawDrawsTrans->At(il+2))->SetLineStyle(3);
      ((TLine *) hitStrawDrawsTrans->At(il+2))->SetLineColor(StrawColor);
      ((TLine *) hitStrawDrawsTrans->At(il+3))->SetX1(tmpStDUpPoint2.getX()/CLHEP::cm);
      ((TLine *) hitStrawDrawsTrans->At(il+3))->SetY1(tmpStDUpPoint2.getY()/CLHEP::cm);
      ((TLine *) hitStrawDrawsTrans->At(il+3))->SetX2(tmpStDDownPoint2.getX()/CLHEP::cm);
      ((TLine *) hitStrawDrawsTrans->At(il+3))->SetY2(tmpStDDownPoint2.getY()/CLHEP::cm);
      ((TLine *) hitStrawDrawsTrans->At(il+3))->SetLineWidth(0.1);
      ((TLine *) hitStrawDrawsTrans->At(il+3))->SetLineStyle(3);
      ((TLine *) hitStrawDrawsTrans->At(il+3))->SetLineColor(StrawColor);
      //((TLine *) hitStrawDrawsTrans->At(i))->SetFillStyle(0);



      tmpiTimeBin=_hHitTime->Fill(hitTime);
      if (tmpiTimeBin>0) timeBin_Straw_rel.insert( stbrel::value_type((unsigned int) tmpiTimeBin, i) );

      //common index for vectors
      //size_t trackIdx(0);

      aveZ=0.0;

      overlapped = false;   // more that one track are crossing the cell/straw in the max drift time
      isFirst = false;      // is the hit of a track that produces the signal (the first in signal arrival time)?

      for (size_t j = 0; j < mcptr.size(); ++j) {

//        StepPointMC const& mchit = (*mchits)[mcptr[j].index];
        StepPointMC const& mchit = *mcptr[j];

        // The simulated particle that made this hit.
        SimParticleCollection::key_type trackId(mchit.trackId());
        SimParticle const& sim = simParticles->at(trackId);

        isFirst=( j==0 );
        if ( isFirst ){
                for ( size_t jj = 1; jj < mcptr.size(); ++jj) {
//                        if ( trackId != (*mchits)[mcptr[jj].index].trackId() ) {
                          if ( trackId != (*(mcptr[jj])).trackId() ) {
                                overlapped=true;
                                break;
                        }
                }
        }

        //cout<<"Hit time :"<<mchit.time()<<" proper "<<mchit.properTime()<<endl;

        if ( sim.pdgId()==11 && sim.fromGenerator() /*&& sim.startMomentum().rho()>=103.0*/ ){
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
                  genTracks.insert( pair<SimParticleCollection::key_type,GenTrackData<TTHitPerTrackData> > ( trackId, GenTrackData<TTHitPerTrackData> (trackId, /*std::pair<CLHEP::Hep3Vector,CLHEP::HepLorentzVector>*/make_pair( sim.startPosition(), sim.startMomentum() ) ) ) );
                  (--(genTracks.end()))->second.addHitData( TTHitPerTrackData( i,mchittime,overlapped,isFirst,stn,layern,sectorn,devicen,mchit.position(),mchit.momentum() ) );
          }
          else {
                  if (genTracks_it->second.isNotPresent(i) ) genTracks_it->second.addHitData( TTHitPerTrackData( i,mchittime,overlapped,isFirst,stn,layern,sectorn,devicen,mchit.position(),mchit.momentum() ) );
          }
        }

//        _hHitLongit->Fill( mchit.position().getX()/CLHEP::cm, mchit.position().getY()/CLHEP::cm );

        aveZ+=mchit.position().getZ();
      }
      aveZ/= ((float) mcptr.size());
      float wpos[3];
//      itwp->WirePosAtZ(aveZ, wpos);

      aveZ/=CLHEP::cm;
      _hHitZ->Fill(aveZ);

      ((TEllipse *) hitDrawsZ->At(i))->SetX1(MCHitPoint.getZ()/CLHEP::cm);
      ((TEllipse *) hitDrawsZ->At(i))->SetY1(MCHitPoint.getY()/CLHEP::cm);
      ((TEllipse *) hitDrawsZ->At(i))->SetR1(driftDistance/CLHEP::cm);
      ((TEllipse *) hitDrawsZ->At(i))->SetR2( ( driftDistance*cos( TMath::ATan2(MCHitPoint.getY(),MCHitPoint.getX()) ) )/CLHEP::cm);
      ((TEllipse *) hitDrawsZ->At(i))->SetFillStyle(0);
      ((TEllipse *) hitDrawsZ->At(i))->SetLineWidth(2.0);
      ((TEllipse *) hitDrawsZ->At(i))->SetLineColor(noiseColor);

      ((TEllipse *) hitDrawsStZprj->At(i))->SetX1(stMidPoint3.getZ()/CLHEP::cm);
      ((TEllipse *) hitDrawsStZprj->At(i))->SetY1(stMidPoint3.getRho()/CLHEP::cm);
      ((TEllipse *) hitDrawsStZprj->At(i))->SetR1(str.getRadius()/CLHEP::cm);
      ((TEllipse *) hitDrawsStZprj->At(i))->SetR2(str.getRadius()/CLHEP::cm);
      ((TEllipse *) hitDrawsStZprj->At(i))->SetFillStyle(1);
      ((TEllipse *) hitDrawsStZprj->At(i))->SetLineWidth(2.0);
      ((TEllipse *) hitDrawsStZprj->At(i))->SetLineColor(noiseColor);

    }

    TEllipse innerWall (0.0,0.0,rIn,rIn);
    innerWall.SetFillStyle(0);
    innerWall.SetLineWidth(1.5);

    TEllipse outerWall (0.0,0.0,rOut,rOut);
    outerWall.SetFillStyle(0);
    outerWall.SetLineWidth(1.5);

    TClonesArray *deviceDrawsLong = new TClonesArray("TBox");
    deviceDrawsLong->ExpandCreateFast(nDevice);
    for ( int idev=0; idev<nDevice; idev++){
            ((TBox *) deviceDrawsLong->At(idev))->SetX1( (ttrdev.at(idev).origin().getZ()-deviceHalfLength)/CLHEP::cm );
            ((TBox *) deviceDrawsLong->At(idev))->SetX2( (ttrdev.at(idev).origin().getZ()+deviceHalfLength)/CLHEP::cm );
            ((TBox *) deviceDrawsLong->At(idev))->SetY1( -rOut );
            ((TBox *) deviceDrawsLong->At(idev))->SetY2( rOut );
            ((TBox *) deviceDrawsLong->At(idev))->SetLineWidth(0.1);
            ((TBox *) deviceDrawsLong->At(idev))->SetFillStyle(0);
            ((TBox *) deviceDrawsLong->At(idev))->SetLineStyle(3);
            ((TBox *) deviceDrawsLong->At(idev))->SetLineColor(kGray/*StrawColor*/);
   }

    std::map<SimParticleCollection::key_type,std::pair<CLHEP::Hep3Vector,CLHEP::HepLorentzVector> >::iterator genEl_it     = genElectrons.begin();
    /*std::map<SimParticleCollection::key_type,std::pair<CLHEP::Hep3Vector,CLHEP::Hep3Vector> >::iterator*/ firstHitOfGenEl_it = firstHitOfGenEl.begin();

    //TClonesArray *genElDraws = new TClonesArray("TEllipse");
    int nGenEl = genElectrons.size();
    //genElDraws->ExpandCreateFast( nGenEl );
    while ( genEl_it != genElectrons.end() ){
      cout<<"Generated el at "<<genEl_it->second.first<<" of "<<genEl_it->second.second<<endl;
      cout<<"First hit in tracker of gen el at "<<firstHitOfGenEl_it->second.first<<" of "<<firstHitOfGenEl_it->second.second<<endl;
      ++genEl_it;
      ++firstHitOfGenEl_it;
    }


    //TClonesArray *hitDrawsZ = new TClonesArray("TEllipse");
    //hitDrawsZ->ExpandCreateFast(nStrawPerEvent);
    std::map<SimParticleCollection::key_type, TClonesArray* > genTracksCrircles;
    std::map<SimParticleCollection::key_type, TClonesArray* >::iterator genTracksCrircles_it;
    TEllipse *trckCircl;

    double ptMeV, rho;
    double B=1.0;
    CLHEP::Hep2Vector radDir;
    HepGeom::Point3D<double> CirCenter;

    int NSignalEle=0;
    std::map <size_t, unsigned int> stHitTrkMrk_rel;
    for ( genTracks_it= genTracks.begin(); genTracks_it!= genTracks.end(); ++genTracks_it ){
//            genTracks_it->second.sort();
            _hHitNloop->Fill (genTracks_it->second.FindNTurns());
            genTracksCrircles.insert( std::pair<SimParticleCollection::key_type, TClonesArray* >( genTracks_it->first,  new TClonesArray("TEllipse") ) );
            ((TClonesArray *) genTracksCrircles[genTracks_it->first])->ExpandCreateFast(genTracks_it->second.getNumOfLoops());
            //cout<<genTracks_it->second;
            for ( unsigned int ilp=0; ilp<genTracks_it->second.getNumOfLoops(); ilp++ ){
                    TTHitPerTrackData hdil = genTracks_it->second.getithLoopHit(ilp);
                    ptMeV = sqrt( pow(hdil.hitMomentum[0],2) + pow(hdil.hitMomentum[1],2) );
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

            if ( genTracks_it->second.getTrkLrntzVec().rho()>=103.0 ) {
                    NSignalEle++;
                    hDP_prevDim = hitDrawsTransPoint->GetEntries();
                    hitDrawsTransPoint->ExpandCreateFast( hDP_prevDim+genTracks_it->second.getNumOfHit() );
                    hitDrawsLongPoint->ExpandCreateFast( hDP_prevDim+genTracks_it->second.getNumOfHit() );
                    for (unsigned int ih=0; ih<genTracks_it->second.getNumOfHit(); ih++){
                            ((TMarker *) hitDrawsTransPoint->At(hDP_prevDim+ih))->SetX(genTracks_it->second.getHit(ih).hitPoint.getX()/CLHEP::cm);
                            ((TMarker *) hitDrawsTransPoint->At(hDP_prevDim+ih))->SetY(genTracks_it->second.getHit(ih).hitPoint.getY()/CLHEP::cm);
                            ((TMarker *) hitDrawsTransPoint->At(hDP_prevDim+ih))->SetMarkerStyle(5+NSignalEle);
                            ((TMarker *) hitDrawsTransPoint->At(hDP_prevDim+ih))->SetMarkerSize(2.0);
                            ((TMarker *) hitDrawsTransPoint->At(hDP_prevDim+ih))->SetMarkerColor(1+NSignalEle);

                            ((TMarker *) hitDrawsLongPoint->At(hDP_prevDim+ih))->SetX(genTracks_it->second.getHit(ih).hitPoint.getZ()/CLHEP::cm);
                            ((TMarker *) hitDrawsLongPoint->At(hDP_prevDim+ih))->SetY(genTracks_it->second.getHit(ih).hitPoint.getY()/CLHEP::cm);
                            ((TMarker *) hitDrawsLongPoint->At(hDP_prevDim+ih))->SetMarkerStyle(5+NSignalEle);
                            ((TMarker *) hitDrawsLongPoint->At(hDP_prevDim+ih))->SetMarkerSize(2.0);
                            ((TMarker *) hitDrawsLongPoint->At(hDP_prevDim+ih))->SetMarkerColor(1+NSignalEle);

                            stHitTrkMrk_rel.insert( pair<size_t,unsigned int>(genTracks_it->second.getHit(ih).iHit,hDP_prevDim+ih));

                            if (genTracks_it->second.getHit(ih).isFirst) {
                                    ((TLine *) hitDrawsTrans->At(genTracks_it->second.getHit(ih).iHit))->SetLineColor(3+NSignalEle);
                                    ((TLine *) hitDrawsTrans->At(genTracks_it->second.getHit(ih).iHit))->SetLineWidth(2.5);
                                    ((TEllipse *) hitDrawsZ->At(genTracks_it->second.getHit(ih).iHit))->SetLineColor(3+NSignalEle);
                                    ((TEllipse *) hitDrawsZ->At(genTracks_it->second.getHit(ih).iHit))->SetLineWidth(2.5);

                            }
                            _hSelHitTime->Fill( hits->at(genTracks_it->second.getHit(ih).iHit).time() );

                            ((TEllipse *) hitDrawsStZprj->At(genTracks_it->second.getHit(ih).iHit))->SetLineColor(2+NSignalEle+ (genTracks_it->second.getHit(ih).devicen%2)+2*genTracks_it->second.getHit(ih).sectorn );

                    }
            }

    }

    _canvasTLview->cd(1);
    _hHitTransverse->Draw();
    innerWall.Draw("same");
    outerWall.Draw("same");
    for ( genTracksCrircles_it=genTracksCrircles.begin(); genTracksCrircles_it!=genTracksCrircles.end(); genTracksCrircles_it++ ){
            for ( unsigned int ilp=0; ilp<genTracksCrircles_it->second->GetEntries(); ilp++ )  ((TEllipse *) genTracksCrircles_it->second->At(ilp))-> Draw("same");
    }
    for (size_t i=0; i<hitDrawsTransPoint->GetEntries(); ++i) ((TMarker *) hitDrawsTransPoint->At(i))->Draw("same");
    for (size_t i=0; i<4*nStrawPerEvent; ++i) ((TLine *) hitStrawDrawsTrans->At(i))->Draw("same");
    for (size_t i=0; i<nStrawPerEvent; ++i) ((TLine *) hitDrawsTrans->At(i))->Draw("same");
    _hHitTransverse->Draw("same");

    _canvasTLview->cd(2);
    _hHitLongit->Draw();
    for (int id=0; id<nDevice; id++) ((TBox *) deviceDrawsLong->At(id))->Draw("same");
    for (size_t i=0; i<hitDrawsLongPoint->GetEntries(); ++i) ((TMarker *) hitDrawsLongPoint->At(i))->Draw("same");
    for (size_t i=0; i<nStrawPerEvent; ++i) ((TEllipse *) hitDrawsZ->At(i))->Draw("same");
    _hHitLongit->Draw("same");

    _canvasStLview->cd();
    _hHitLongit->Draw();
    for (int id=0; id<nDevice; id++) ((TBox *) deviceDrawsLong->At(id))->Draw("same");
    //for (size_t i=0; i<hitDrawsLongPoint->GetEntries(); ++i) ((TMarker *) hitDrawsLongPoint->At(i))->Draw("same");
    for (size_t i=0; i<nStrawPerEvent; ++i) ((TEllipse *) hitDrawsStZprj->At(i))->Draw("same");
    _hHitLongit->Draw("same");


    _canvasPl->cd(1);
    //_hHitTime->SetStats(kTRUE);
    _hHitTime->Draw();
    _hSelHitTime->Draw("same");

    int nbingroup= (int) intTimeWind/timeBinDim;
    float hit=0.0;
    float Selhit=0.0;
    int binHalf = (int) nbingroup/2;
    int jbin;
    //int nBin=_hHitTime->GetNbinsX();
    for (int it=1; it<=ntimeBin; it++) {
      hit=0.0;
      Selhit=0.0;
      for (int is=0; is<nbingroup; is++) {
        jbin = it -binHalf +is;
        if (jbin<1) continue;
        if (jbin>ntimeBin) break;
        hit+=_hHitTime->GetBinContent(jbin);
        Selhit+=_hSelHitTime->GetBinContent(jbin);
      }
      _hHitClustTime->SetBinContent(it,hit);
      _hSelHitClustTime->SetBinContent(it,Selhit);
    }

    _canvasPl->cd(2);

    float thr = 10.0;
    int nfound = 0;
    if (_hHitTime->GetEntries()>0){
            thr/=_hHitClustTime->GetMaximum();
            if (thr>0.0 && thr<1.0) nfound = _peakFinder->Search( _hHitClustTime,2.61,"",thr );
    }

    _hHitClustTime->Draw();
    _hSelHitClustTime->Draw("same");

    _fg->SetParameter(1,_hSelHitClustTime->GetMean());
    _fg->SetParameter(2,_hSelHitClustTime->GetRMS());
    _hSelHitClustTime->Fit("fg","0");
    _hSelHitClustTimeSigma->Fill(_fg->GetParameter(2));

//    _hSelHitClustTimeSigma->GetFunction("fg")->Delete();


    for ( stbrel::const_iterator stb_it=timeBin_Straw_rel.begin(); stb_it!=timeBin_Straw_rel.end(); ++stb_it ) {
            cout<<"i-th time bin "<<stb_it->first<<" hit id "<<stb_it->second<<endl;
    }

    TH2I *tmpSecClustDist=0x0;
    TH2I *tmpStClustDist=0x0;
    if (nfound>0){
            Float_t *timepeakPos=_peakFinder->GetPositionX();
            Float_t *timepeakHei=_peakFinder->GetPositionY();
            unsigned int *timepeakPosId = new unsigned int[nfound];
            unsigned int frstTimeBinInP, lastTimeBinInP;
            size_t ihit;
            std::map <size_t, unsigned int>::iterator stHitTrkMrk_rel_it;

            int i1peak;

            StrawId sid;
            int stn, layern, devicen, sectorn;
            unsigned int absSect, devStId;

            tmpStClustDist    = new TH2I("tmpStClustDist","tmp Smoothing of Device vs Straw multiplicity Distribution",36,0,36,1000,-200,800);
            tmpSecClustDist   = new TH2I("tmpSecClustDist","tmp Smoothing of Device vs Sector multiplicity Distribution",36,0,36,20,-4,16);

            for (int ipeak=0; ipeak<nfound; ipeak++){

                    i1peak=ipeak;
                    ++i1peak;
                    _peaksCanvases->AddAtAndExpand(new TCanvas(Form("canvasFor_%d-th_peak",i1peak),Form("Data of the peak at %f ns",timepeakPos[ipeak]),860,860),ipeak);
                    _peaksCanvHistos->AddAtAndExpand(new TCanvas(Form("canvasForHistos_%d-th_peak",i1peak),Form("Histograms of the peak at %f ns",timepeakPos[ipeak]),860,860),ipeak);
                    _hPkStDistTrs->AddAtAndExpand(new TH2I(Form("h%d-Pk_StDistTrs",i1peak),Form("Sector-Straw Distribution in Transverse plane for the %d-th peak",i1peak),20,-4,16,50,0,50),ipeak);
                    _hPkStDistanceTrs->AddAtAndExpand(new TH2I(Form("h%d-Pk_StDistanceTrs",i1peak),Form("Sector-Straw distance Distribution in Transverse plane for the %d-th peak",i1peak),20,-4,16,50,0,50),ipeak);
                    _hPkStDist2W->AddAtAndExpand(new TH2I(Form("h%d-Pk_StDist2W",i1peak),Form("Device vs Straw multiplicity Distribution for the %d-th peak",i1peak),36,0,36,1000,-200,800),ipeak);
                    _hPkStClustDist2W->AddAtAndExpand(new TH2I(Form("h%d-Pk_StClustDist2W",i1peak),Form("Smoothing of Device vs Straw multiplicity Distribution for the %d-th peak",i1peak),36,0,36,1000,-200,800),ipeak);
                    _hPkSecDist2W->AddAtAndExpand(new TH2I(Form("h%d-Pk_SecDist2W",i1peak),Form("Device vs Sector multiplicity Distribution for the %d-th peak",i1peak),36,0,36,20,-4,16),ipeak);
                    _hPkSecClustDist2W->AddAtAndExpand(new TH2I(Form("h%d-Pk_SecClustDist2W",i1peak),Form("Smoothing of Device vs Sector multiplicity Distribution for the %d-th peak",i1peak),36,0,36,20,-4,16),ipeak);
                    _hPkSecVsStationDist2W->AddAtAndExpand(new TH2I(Form("h%d-Pk_SecDist2W",i1peak),Form("Station vs Sector multiplicity Distribution for the %d-th peak",i1peak),18,0,18,16,0,16),ipeak);
                    _hPkSecVsStationDist2WGood->AddAtAndExpand(new TH2I(Form("h%d-Pk_SecDist2WG",i1peak),Form("Selected Station vs Sector multiplicity Distribution for the %d-th peak",i1peak),18,0,18,16,0,16),ipeak);
                    _hPkAccArrSecStation->AddAtAndExpand(new TH2I(Form("h%d-Pk_AccArrSecStation",i1peak),Form("Accumulator Array for Station distance vs Sector for the %d-th peak",i1peak),9,3,12,32,-16,16),ipeak);

                    TCanvas * iCanv=((TCanvas *) _peaksCanvases->At(ipeak));
                    iCanv->Divide(2,2);
                    iCanv->cd(1);
                    _hHitTransverse->Draw();
                    innerWall.Draw("same");
                    outerWall.Draw("same");

                    iCanv->cd(2);
                    _hHitLongit->Draw();
                    for (int id=0; id<nDevice; id++) ((TBox *) deviceDrawsLong->At(id))->Draw("same");

                    cout<<"----------"<<endl;
                    timepeakPosId[ipeak]=(unsigned int) (((unsigned int)timepeakPos[ipeak])/timeBinDim)+1;

                    cout<<"Peak pos "<<timepeakPos[ipeak]<<" = "<<timepeakPosId[ipeak]<<" bin center "<<_hHitTime->GetBinCenter(timepeakPosId[ipeak])<<endl;
                    cout<<"peak height "<<timepeakHei[ipeak]<<endl;
                    cout<<" n hit found at paek pos "<<timeBin_Straw_rel.count(timepeakPosId[ipeak])<<endl;
                    cout<<"First hit at peak pos "<<timeBin_Straw_rel.find(timepeakPosId[ipeak])->second<<endl;
                    //timeBin_Straw_rel

                    frstTimeBinInP = timepeakPosId[ipeak] - binHalf -1;
                    lastTimeBinInP = timepeakPosId[ipeak] + binHalf +1;
                    cout<<"Peak range limits: "<<frstTimeBinInP<<" "<<lastTimeBinInP<<endl;
                    timeBin_Straw_rel.begin();
                    stbrel::const_iterator stb_it=timeBin_Straw_rel.find(frstTimeBinInP);
                    int findFTB=frstTimeBinInP;
                    while (stb_it==timeBin_Straw_rel.end()) {
                            stb_it=timeBin_Straw_rel.begin();
                            ++findFTB;
                            stb_it=timeBin_Straw_rel.find(findFTB);
                    }

                    int hitFoundInPeak=0;
                    TH2I *ihPkStDistTrs      = ((TH2I *) _hPkStDistTrs->At(ipeak));
                    ihPkStDistTrs->SetStats(kFALSE);
                    ihPkStDistTrs->SetXTitle("absSect");
                    ihPkStDistTrs->SetYTitle("Straw");
                    ihPkStDistTrs->SetTitleOffset(1.2,"y");
//                    TH2I *ihPkStDistanceTrs  = ((TH2I *) _hPkStDistanceTrs->At(ipeak));
//                    ihPkStDistanceTrs->SetStats(kFALSE);
//                    ihPkStDistanceTrs->SetXTitle("absSect");
//                    ihPkStDistanceTrs->SetYTitle("Straw");
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

                    TH2I *ihPkSecVsStationDist2W      = ((TH2I *) _hPkSecVsStationDist2W->At(ipeak));
                    ihPkSecVsStationDist2W->SetStats(kFALSE);
                    ihPkSecVsStationDist2W->SetXTitle("Station");
                    ihPkSecVsStationDist2W->SetYTitle("absSect");
                    ihPkSecVsStationDist2W->SetTitleOffset(1.2,"y");

                    TH2I *ihPkSecVsStationDist2WG      = ((TH2I *) _hPkSecVsStationDist2WGood->At(ipeak));
                    ihPkSecVsStationDist2WG->SetStats(kFALSE);
                    ihPkSecVsStationDist2WG->SetXTitle("Station");
                    ihPkSecVsStationDist2WG->SetYTitle("absSect");
                    ihPkSecVsStationDist2WG->SetTitleOffset(1.2,"y");

                    TH2I *ihPkAccArrSecStation      = ((TH2I *) _hPkAccArrSecStation->At(ipeak));
                    ihPkAccArrSecStation->SetStats(kFALSE);
                    ihPkAccArrSecStation->SetXTitle("Pitch [Station]");
                    ihPkAccArrSecStation->SetYTitle("absSect");
                    ihPkAccArrSecStation->SetTitleOffset(1.2,"y");

                    tmpStClustDist->Reset();
                    tmpSecClustDist->Reset();

                    TCanvas * iCanvHist=((TCanvas *) _peaksCanvHistos->At(ipeak));
                    iCanvHist->Divide(2,2);

                    ssmap_Bin_Straw_rel.clear();
                    while (stb_it!=timeBin_Straw_rel.end()){
                            if (stb_it->first>lastTimeBinInP) break;
                            ihit=stb_it->second;
                            ++hitFoundInPeak;

                            iCanv->cd(1);
                            stHitTrkMrk_rel_it=stHitTrkMrk_rel.find(ihit);
                            if (stHitTrkMrk_rel_it!=stHitTrkMrk_rel.end()) ((TMarker *) hitDrawsTransPoint->At(stHitTrkMrk_rel_it->second))->Draw("same");
                            for (size_t isl=0; isl<4; ++isl) ((TLine *) hitStrawDrawsTrans->At(4*ihit+isl))->Draw("same");
                            ((TLine *) hitDrawsTrans->At(ihit))->Draw("same");

                            iCanv->cd(2);
                            if (stHitTrkMrk_rel_it!=stHitTrkMrk_rel.end()) ((TMarker *) hitDrawsLongPoint->At(stHitTrkMrk_rel_it->second))->Draw("same");
                            ((TEllipse *) hitDrawsZ->At(ihit))->Draw("same");

                            StrawHit        const&      hit(hits->at(ihit));
                            //Get hit straw
                            StrawIndex si = hit.strawIndex();
                            const Straw & str = tracker.getStraw(si);

                            // cout << "Getting informations about cells" << endl;

                            sid = str.id();
                            stn     = sid.getStraw();
                            layern  = sid.getLayer();
                            devicen = sid.getDevice();
                            sectorn = sid.getSector();

                            absSect = (devicen%2)+2*sectorn;
                            devStId = absSect*50+stn;

                            //ihPkStDistTrs->Fill(devStId);
                            ihPkStDistTrs->Fill(absSect,stn);
                            if (absSect>=8 && absSect<=11)  ihPkStDistTrs->Fill(absSect-12,stn);
                            if (absSect>=0 && absSect<=3)   ihPkStDistTrs->Fill(absSect+12,stn);

                            ihPkSecDist2W->Fill(devicen,absSect);
                            if (absSect>=8 && absSect<=11)  ihPkSecDist2W->Fill(devicen,absSect-12);
                            if (absSect>=0 && absSect<=3)   ihPkSecDist2W->Fill(devicen,absSect+12);

                            ihPkStDist2W->Fill(devicen,absSect*50+stn);
                            if (absSect>=8 && absSect<=11)  ihPkStDist2W->Fill(devicen,(absSect-12)*50+stn);
                            if (absSect>=0 && absSect<=3)   ihPkStDist2W->Fill(devicen,(absSect+12)*50+stn);
//                            ihPkStDist2W->Fill(devicen,devStId);
/*
                            int nYbingroup=250;//5*50
                            float Yhit=0.0;
                            int YbinHalf = (int) nYbingroup/2;
                            int jYbin;
                            int nXBin=ihPkStDist2W->GetNbinsX();
                            int nYBin=ihPkStDist2W->GetNbinsY();
                            for (int itX=1; itX<=nXBin; itX++) {
                                    for (int itY=1; itY<=nYBin; itY++) {
                                            Yhit=0.0;
                                            for (int is=0; is<nYbingroup; is++) {
                                                    jYbin = itY -YbinHalf +is;
                                                    if (jYbin<1) continue;
                                                    if (jYbin>nYBin) break;
                                                    Yhit+=ihPkStDist2W->GetBinContent(itX,jYbin);
                                            }
                                            tmpStClustDist->SetBinContent(itX,itY,Yhit);
                                    }
                            }
                            int nXbingroup=9;
                            float Xhit=0.0;
                            int XbinHalf = (int) nXbingroup/2;
                            int jXbin;
                            for (int itY=1; itY<=nYBin; itY++) {
                                    for (int itX=1; itX<=nXBin; itX++) {
                                            Xhit=0.0;
                                            for (int is=0; is<nXbingroup; is++) {
                                                    jXbin = itX -XbinHalf +is;
                                                    if (jXbin<1) continue;
                                                    if (jXbin>nXBin) break;
                                                    Xhit+=tmpStClustDist->GetBinContent(jXbin,itY);
                                            }
                                            ihPkStClustDist2W->SetBinContent(itX,itY,Xhit);
                                    }
                            }
*/

                            tmpiGeomBin = ihPkSecVsStationDist2W->Fill(devicen/2,absSect);
                            ssmap_Bin_Straw_rel.insert( stbrel::value_type((unsigned int) tmpiGeomBin, ihit) );
                            //if (absSect>=8 && absSect<=11)  ihPkSecVsStationDist2W->Fill(devicen/2,absSect-12);
                            if (absSect>=0 && absSect<=3)   ihPkSecVsStationDist2W->Fill(devicen/2,absSect+12);

                            int nYbingroup, YbinHalf, jYbin, nXBin, nYBin;
                            float Yhit;
                            int nXbingroup, XbinHalf, jXbin;
                            float Xhit;

////                            nYbingroup=250;//5*50
////                            Yhit=0.0;
////                            YbinHalf = (int) nYbingroup/2;
////                            jYbin;
////                            nXBin=ihPkStDist2W->GetNbinsX();
////                            nYBin=ihPkStDist2W->GetNbinsY();
////                            for (int itX=1; itX<=nXBin; itX++) {
////                                    for (int itY=1; itY<=nYBin; itY++) {
////                                            Yhit=0.0;
////                                            for (int is=0; is<nYbingroup; is++) {
////                                                    jYbin = itY -YbinHalf +is;
////                                                    if (jYbin<1) continue;
////                                                    if (jYbin>nYBin) break;
////                                                    Yhit+=ihPkStDist2W->GetBinContent(itX,jYbin);
////                                            }
////                                            tmpStClustDist->SetBinContent(itX,itY,Yhit);
////                                    }
////                            }
////                            nXbingroup=7;//9
////                            Xhit=0.0;
////                            XbinHalf = (int) nXbingroup/2;
////                            jXbin;
////                            for (int itY=1; itY<=nYBin; itY++) {
////                                    for (int itX=1; itX<=nXBin; itX++) {
////                                            Xhit=0.0;
////                                            for (int is=0; is<nXbingroup; is++) {
////                                                    jXbin = itX -XbinHalf +is;
////                                                    if (jXbin<1) continue;
////                                                    if (jXbin>nXBin) break;
////                                                    Xhit+=tmpStClustDist->GetBinContent(jXbin,itY);
////                                            }
////                                            ihPkStClustDist2W->SetBinContent(itX,itY,Xhit);
////                                    }
////                            }
//
//                            nYbingroup=5;//5*50
//                            Yhit=0.0;
//                            YbinHalf = (int) nYbingroup/2;
//                            jYbin;
//                            nXBin=ihPkStDist2W->GetNbinsX();
//                            nYBin=ihPkStDist2W->GetNbinsY();
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
//                            nXbingroup=7;//9
//                            Xhit=0.0;
//                            XbinHalf = (int) nXbingroup/2;
//                            jXbin;
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

                            ++stb_it;
                   }


                   cout<<"------- start research ------"<<endl;
                   fillVotingArray(ihPkSecVsStationDist2W, ihPkAccArrSecStation);





                   iCanv->cd(3);
                   ihPkStDistTrs->Draw("col z");

                   iCanv->cd(4);
                   ihPkSecDist2W->Draw("col z");

//                   iCanvHist->cd(1);
//                   iCanvHist->cd(1)->SetGrid();
//                   ihPkSecDist2W->Draw("col z");

                   iCanvHist->cd(1);//2
                   iCanvHist->cd(1)->SetGrid();//2
                   ihPkSecVsStationDist2W->Draw("col z");

                   clstSegments.push_back( std::vector< clssegdr >() );
                   int icl=0;
                   for ( std::vector<Clust>::iterator clstlst_it=clustersList.begin(); clstlst_it!=clustersList.end(); ++clstlst_it ) {
                           if ( clstlst_it->_m!=std::numeric_limits<float>::infinity() ) {
                                   clustSegline.SetParameter(0,clstlst_it->_q);
                                   clustSegline.SetParameter(1,clstlst_it->_m);
                                   clstSegments.back().push_back( clssegdr(new TArrow( (float)clstlst_it->_minStationID, clustSegline.Eval((float)clstlst_it->_minStationID),
                                                   (float)(clstlst_it->_maxStationID+1), clustSegline.Eval((float)(clstlst_it->_maxStationID+1)), 0.02 ) ) );
                           } else {
                                   clstSegments.back().push_back( clssegdr(new TArrow( clstlst_it->_q, (float)clstlst_it->_firstSectorID, clstlst_it->_q, (float) (clstlst_it->_lastSectorID+1), 0.02 ) ) );
                           }
                           //clstSegments.back().back()->SetLineStyle(2);
                           clstSegments.back().back()->SetLineColor(kBlack);
                           for ( rwclincl::iterator rwclincl_it=clstlst_it->_rClusts.begin(); rwclincl_it!=clstlst_it->_rClusts.end(); ++rwclincl_it) {
                                   for ( int ibin=rwclincl_it->second->_firstStationID; ibin<=rwclincl_it->second->_lastStationID; ibin++) {
                                           ihPkSecVsStationDist2WG->SetBinContent(ibin+1,rwclincl_it->first+1,
                                                           ihPkSecVsStationDist2W->GetBinContent(ibin+1,rwclincl_it->first+1));
                                   }
                           }
                           icl++;
                   }
                   iCanvHist->cd(2);
                   iCanvHist->cd(2)->SetGrid();
                   ihPkSecVsStationDist2WG->Draw("col z");
                   for ( std::vector< clssegdr >::iterator clsS_it=clstSegments.back().begin(); clsS_it!=clstSegments.back().end(); ++clsS_it) (*clsS_it)->Draw();


//                            iCanvHist->cd(2);
//                            ihPkStDist2W->Draw("col z");

//                            iCanvHist->cd(3);
//                            tmpStClustDist->Draw("col z");

//                   iCanvHist->cd(3);
//                   ihPkSecClustDist2W->Draw("col z");




//                            iCanvHist->cd(4);
//                            ihPkStClustDist2W->Draw("col z");
                   iCanvHist->cd(3);//4
                   ihPkAccArrSecStation->Draw("col z");


                   cout<<"Total Hit found in "<<ipeak<<"-th peak: "<<hitFoundInPeak<<endl;

                   iCanv->cd(1);
                   _hHitTransverse->Draw("same");

                   iCanv->cd(2);
                   _hHitLongit->Draw("same");


//                    for ( genTracksCrircles_it=genTracksCrircles.begin(); genTracksCrircles_it!=genTracksCrircles.end(); genTracksCrircles_it++ ){
//                            for ( unsigned int ilp=0; ilp<genTracksCrircles_it->second->GetEntries(); ilp++ )  ((TEllipse *) genTracksCrircles_it->second->At(ilp))-> Draw("same");
//                    }

                    //    iCanv->cd(3);
                    ////    _hHitStereoMi->Draw();
                    ////    innerWall.Draw("same");
                    ////    outerWall.Draw("same");
                    ////    for (size_t i=0; i<nStrawPerEvent; ++i) if ( isStereoMin[i] ) ((TEllipse *) hitDrawsTrans->At(i))->Draw("same");
                    ////    _hHitStereoMi->Draw("same");
                    //
                    //    iCanv->cd(4);
                    ////    _hHitStereoPl->Draw();
                    ////    innerWall.Draw("same");
                    ////    outerWall.Draw("same");
                    ////    for (size_t i=0; i<nStrawPerEvent; ++i) if ( !isStereoMin[i] ) ((TEllipse *) hitDrawsTrans->At(i))->Draw("same");
                    ////    _hHitStereoPl->Draw("same");

                    iCanv->Modified();
                    iCanv->Update();

            }
            delete timepeakPosId;
    }


    _canvasPl->cd(3);
    _hHitZ->Draw();

    _canvasPl->cd(4);
    _hHitNloop->Draw();
//    _hMCHitDTime->Draw();

    _canvasTLview->Modified();
    _canvasTLview->Update();

    _canvasStLview->Modified();
    _canvasStLview->Update();

    _canvasPl->Modified();
    _canvasPl->Update();

    _cnvForPeakstudy->cd();
    _hSelHitClustTimeSigma->Draw();
    _cnvForPeakstudy->Modified();
    _cnvForPeakstudy->Update();


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
    _fakeCanvas->WaitPrimitive();
    cerr << endl;
    delete printEvN;

    _peaksCanvases->Delete();
    _peaksCanvHistos->Delete();
    _hPkStDistTrs->Delete();
    _hPkStDistanceTrs->Delete();
    _hPkStDist2W->Delete();
    _hPkStClustDist2W->Delete();
    _hPkSecDist2W->Delete();
    _hPkSecClustDist2W->Delete();
    _hPkSecVsStationDist2W->Delete();
    _hPkSecVsStationDist2WGood->Delete();

    //for (size_t i=0; i<nStrawPerEvent; ++i) { delete ((TEllipse *) hitDrawsTrans->At(i));  delete ((TEllipse *) hitDrawsZ->At(i)); }
    hitDrawsTrans->Delete();
    hitStrawDrawsTrans->Delete();
    hitDrawsZ->Delete();
    hitDrawsStZprj->Delete();
    hitDrawsTransPoint->Delete();
    hitDrawsLongPoint->Delete();
    deviceDrawsLong->Delete();
    delete hitDrawsTrans;
    delete hitStrawDrawsTrans;
    delete hitDrawsZ;
    delete hitDrawsStZprj;
    delete hitDrawsTransPoint;
    delete hitDrawsLongPoint;
    delete isStereoMin;
    delete deviceDrawsLong;


    if (tmpStClustDist) delete tmpStClustDist;
    if (tmpSecClustDist) delete tmpSecClustDist;

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

  void TTDisplayData::endJob(){

    // cd() to correct root directory. See note 3.
    TDirectory* save = gDirectory;
    _directory->cd();

    // Write canvas.  See note 4.
//    _canvas->Write();

    // cd() back to where we were.  See note 3.
    save->cd();

  }

  void TTDisplayData::fillVotingArray(TH2I *startDist, TH2I *votingArray, int minPitch, int maxPitch){

          int nXbin = startDist->GetNbinsX();
          int nYbin = startDist->GetNbinsY();
          int lastXbin = nXbin-1;
          int *contArr = startDist->GetArray();
          int dataArr[nYbin][nXbin];
          int dataArrRev[nYbin][nXbin];
          int iY, iX, jXRange, jOutOfRange;
          int minAvaiPitch, tmpMinPitch, lastFullBin, tmpPitch;

          int nVotXbin=votingArray->GetNbinsX();
          int iVotX;
          TH1I tmpVoting ("tmpVoting","temporany container",nVotXbin,minPitch,maxPitch);
          TH1I tmpVotingRev ("tmpVotingRev","temporany container",nVotXbin,minPitch,maxPitch);
          int *votArr = tmpVoting.GetArray();
          int *votArrRev = tmpVotingRev.GetArray();

          ptcsctrrel tmp_ptc_sctr_rel;
          ptcbnrel sttn_sectPitch_rel;

          goodsttnrel tmp_good_sctr_val;
          goodsctrstatnrel good_sttn_sctr_rel;

          rwclstvec rowClusts;
          rwclstvec::iterator tmp_rwcl_it;

//          rwclinrwrel rwClst_forRw_rel;
          rwClst_forRw_rel.clear();

          stbrel rowClst_Column_rel;
          size_t nRowClust=-1;

          ptcclsbrel tmp_ptc_rowClst__rel;
          ptcmaprowclsrel ptcMap_rowCls_rel;
          rwavptcclscpl rwMap_avPtc_rwClrwCl;
          //rwavptcclscpl::iterator ptcrwClrwClmap_it;
          avptcclscpl tmp_avPtc_rwClrwCl;
          avptcclscpl::iterator ptcrwClrwCl_it;
          size_t couplingRwCl;
          rwclclcpl tmp_rwClrwClpair;


//          for (int i=0; i<nYbin; i++) { for(int j=0; j<nXbin; j++) cout<<contArr[(i+1)*(nXbin+2)+1+j]<<" "; cout<<endl;}

          for ( iY=0; iY<nYbin; iY++) {
                  memcpy( dataArr[iY], contArr+( ( (nXbin+2)*(iY+1) ) +1 ), nXbin*4);
                  std::reverse_copy(dataArr[iY],dataArr[iY]+nXbin,dataArrRev[iY]);

                  tmpVoting.Reset();
                  tmpVotingRev.Reset();

                  minAvaiPitch=minPitch;
                  lastFullBin=-1;
                  tmp_ptc_sctr_rel.clear();
                  tmp_ptc_rowClst__rel.clear();
                  tmp_avPtc_rwClrwCl.clear();

                  for ( iX=0; iX<nXbin; iX++) {
                          if ( dataArr[iY][iX] ) {
//                                  if ( iX==0 || (iX-lastFullBin)>1 ) {
//                                          rowClusts.push_back(rowClust( (unsigned short) iX, (unsigned short) dataArr[iY][iX]));
//                                          nRowClust++;
//                                          rowClst_Column_rel.insert( stbrel::value_type((unsigned int) iY, nRowClust) );
//                                  }
//                                  else rowClusts.back().addHit( (unsigned short) iX, (unsigned short) dataArr[iY][iX]);
                                  lastFullBin=iX;
                                  for ( jXRange=iX+minAvaiPitch; jXRange<(iX+maxPitch); jXRange++){
                                          if ( jXRange>lastXbin ) break;
//                                          if ( jXRange>=lastXbin ) {
//                                                  tmpMinPitch=jXRange-iX;
//                                                  jOutOfRange= (tmpMinPitch>minAvaiPitch) ? tmpMinPitch : minAvaiPitch;
//                                                  for ( ; jOutOfRange<=maxPitch; jOutOfRange++) votingArray->Fill(jOutOfRange,iY);
//                                                  break;
//                                          }
                                          if ( dataArr[iY][jXRange] ) {
                                                  tmpPitch=jXRange-iX;
                                                  votingArray->Fill(tmpPitch,iY);
                                                  tmpVoting.Fill(tmpPitch);
                                                  tmp_ptc_sctr_rel.insert ( ptcsctrrel::value_type(
                                                                  (unsigned short) tmpPitch,
                                                                  (unsigned short) iX )
                                                  );

//                                                  tmp_ptc_rowClst__rel.insert(
//                                                                  ptcclsbrel::value_type(
//                                                                                  (unsigned short) tmpPitch, std::pair<size_t,unsigned short>(nRowClust,(unsigned short) jXRange)
//                                                                                  )
//                                                  );

                                          }
                                  }
                                  //minAvaiPitch=
                          }
                          else {
                                  tmpMinPitch=iX-lastFullBin;
                                  if ( tmpMinPitch>minAvaiPitch ) minAvaiPitch=tmpMinPitch;
                          }
                          if (minAvaiPitch>=maxPitch) break;
                  }
                  //if (tmp_ptc_rowClst__rel.size()) ptcMap_rowCls_rel.insert( ptcmaprowclsrel::value_type((unsigned short) iY, tmp_ptc_rowClst__rel) );

                  //in reverse order
                  minAvaiPitch=minPitch;
                  lastFullBin=-1;
                  for ( iX=0; iX<nXbin; iX++) {
                          if ( dataArrRev[iY][iX] ) {
                                  lastFullBin=iX;
                                  for ( jXRange=iX+minAvaiPitch; jXRange<(iX+maxPitch); jXRange++){
                                          if ( jXRange>lastXbin ) break;
//                                          if ( jXRange>=lastXbin ) {
//                                                  tmpMinPitch=jXRange-iX;
//                                                  jOutOfRange= (tmpMinPitch>minAvaiPitch) ? tmpMinPitch : minAvaiPitch;
//                                                  for ( ; jOutOfRange<=maxPitch; jOutOfRange++) votingArray->Fill(jOutOfRange,iY-nYbin);
//                                                  break;
//                                          }
                                          if ( dataArrRev[iY][jXRange] ) {
                                                  tmpPitch=jXRange-iX;
                                                  votingArray->Fill(tmpPitch,iY-nYbin);
                                                  tmpVotingRev.Fill(tmpPitch);
                                                  tmp_ptc_sctr_rel.insert ( ptcsctrrel::value_type(
                                                                  (unsigned short) tmpPitch,
                                                                  (unsigned short) lastXbin-iX )
                                                  );
                                         }
                                  }
                                  //minAvaiPitch=
                          }
                          else {
                                  tmpMinPitch=iX-lastFullBin;
                                  if ( tmpMinPitch>minAvaiPitch ) minAvaiPitch=tmpMinPitch;
                          }
                          if (minAvaiPitch>=maxPitch) break;
                  }

                  for ( iVotX=1; iVotX<=nVotXbin; iVotX++) {
                          tmpMinPitch=iVotX-1+minPitch;  //is the pitch under check
                          //if ( votArr[iVotX] && votArrRev[iVotX] ) cout<<"found"<<endl;
                          if ( (votArr[iVotX] && !votArrRev[iVotX]) || (!votArr[iVotX] && votArrRev[iVotX]) ) {
                                  ptcsctrrel::iterator ptcst_it=tmp_ptc_sctr_rel.find(tmpMinPitch);
                                  while ( ptcst_it != tmp_ptc_sctr_rel.end() ){
                                          tmp_ptc_sctr_rel.erase(ptcst_it);
                                          ptcst_it=tmp_ptc_sctr_rel.find(tmpMinPitch);
                                  }
                          }
//                          if ( votArr[iVotX] && !votArrRev[iVotX] ) {
//
//
//                                  cout<<"\t !!!! I'm removing vote at row "<<iY<<" and pitch "<<tmpMinPitch<<endl;
//
//                                  ptcclsbrel::iterator prcc_it=tmp_ptc_rowClst__rel.find(tmpMinPitch);
//                                  while ( prcc_it != tmp_ptc_rowClst__rel.end() ){
//                                          tmp_ptc_rowClst__rel.erase(prcc_it);
//                                          prcc_it=tmp_ptc_rowClst__rel.find(tmpMinPitch);
//                                  }
//                          }
                  }
                  if (tmp_ptc_sctr_rel.size()) {
                          sttn_sectPitch_rel.insert( ptcbnrel::value_type((unsigned short) iY, tmp_ptc_sctr_rel) );
                          tmp_good_sctr_val.clear();
                          for ( ptcsctrrel::iterator ptcst_it=tmp_ptc_sctr_rel.begin(); ptcst_it!=tmp_ptc_sctr_rel.end(); ++ptcst_it ) {
                                  if (tmp_good_sctr_val.find(ptcst_it->second)==tmp_good_sctr_val.end()) tmp_good_sctr_val.insert( goodsttnrel::value_type(ptcst_it->second) );
                          }

                          unsigned short lastElement;
                          unsigned short iEl;
                          goodsttnrel::iterator ptcst_it=tmp_good_sctr_val.begin();

                          lastElement=*ptcst_it;
                          rowClusts.push_back( rwClustPtr( new rowClust(lastElement, (unsigned short) dataArr[iY][lastElement]) ) );
                          rwClst_forRw_rel.insert( rwclinrwrel::value_type( (unsigned short) iY, rwclstvec(1,rowClusts.back()) ) );
                          nRowClust++;
                          rowClst_Column_rel.insert( stbrel::value_type((unsigned int) iY, nRowClust) );
                          ptcst_it++;
                          for ( iEl=1; iEl<(unsigned short)tmp_good_sctr_val.size(); iEl++ ) {
                                  //cout<<"Element "<<(*ptcst_it)<<" previous "<<lastElement<<endl;
                                  if ( ((*ptcst_it)-lastElement)>1 ) {
                                  //        cout<<"New Cluster added"<<endl;
                                          rowClusts.push_back( rwClustPtr ( new rowClust( (*ptcst_it), (unsigned short) dataArr[iY][(*ptcst_it)]) ) );
                                          rwClst_forRw_rel[(unsigned short) iY].push_back(rowClusts.back());
                                          nRowClust++;
                                          rowClst_Column_rel.insert( stbrel::value_type((unsigned int) iY, nRowClust) );
                                  }
                                  else {
                                  //        cout<<"hit added to Cluster"<<endl;
                                          rowClusts.back()->addHit( (*ptcst_it), (unsigned short) dataArr[iY][(*ptcst_it)]);
                                  }
                                  lastElement=*ptcst_it;
                                  ptcst_it++;
                          }

                          good_sttn_sctr_rel.insert( goodsctrstatnrel::value_type((unsigned short) iY, tmp_good_sctr_val) );
                  }
//                  if (tmp_ptc_rowClst__rel.size()) {
//
//                          ptcMap_rowCls_rel.insert( ptcmaprowclsrel::value_type((unsigned short) iY, tmp_ptc_rowClst__rel) );
//
//                          for ( ptcclsbrel::iterator prcc_it=tmp_ptc_rowClst__rel.begin(); prcc_it!=tmp_ptc_rowClst__rel.end(); prcc_it++ ){
//                                  couplingRwCl = prcc_it->second.first;
//                                  tmp_rwClrwClpair.first=couplingRwCl;
//                                  couplingRwCl++;
//                                  while ( couplingRwCl<rowClusts.size()-1 ){
//                                          tmpRwClust=&(rowClusts.at(couplingRwCl));
//                                          if ( prcc_it->second.second>=tmpRwClust->_firstStationID && prcc_it->second.second<=tmpRwClust->_lastStationID ) break;
//                                          couplingRwCl++;
//                                  }
//                                  tmp_rwClrwClpair.second=couplingRwCl;
//
//                                  ptcrwClrwCl_it = tmp_avPtc_rwClrwCl.find(tmp_rwClrwClpair);
//                                  if ( ptcrwClrwCl_it==tmp_avPtc_rwClrwCl.end() )
//                                          tmp_avPtc_rwClrwCl.insert (
//                                                  avptcclscpl::value_type(
//                                                                  tmp_rwClrwClpair, std::vector<unsigned short> ( (size_t) 1, (unsigned short) prcc_it->first )
//                                                                  )
//                                          );
//                                  else {
//                                          ptcrwClrwCl_it->second.push_back( (unsigned short) prcc_it->first );
//                                  }
//
//                          }
//
//                          //tmp_avPtc_rwClrwCl.insert( );
//                          rwMap_avPtc_rwClrwCl.insert( rwavptcclscpl::value_type( (unsigned short) iY, tmp_avPtc_rwClrwCl ) );
//                  }

          }




          //std::vector<Clust> clustersList;
          clustersList.clear();
          rwClustPtr startingRowClust;

          for ( rwclinrwrel::iterator rcr_it=rwClst_forRw_rel.begin(); rcr_it!=rwClst_forRw_rel.end(); ++rcr_it ) {
//                  cout<<"i-th row "<<rcr_it->first<<endl;
//                  for ( rwclstvec::iterator cr_it=rcr_it->second.begin();  cr_it!=rcr_it->second.end();  ++cr_it  ) {
//                          rwClustPtr &irClus = *cr_it;
//                          cout<<"\t cluster data: "<<irClus->_firstStationID<<" "<<irClus->_lastStationID<<" "<<irClus->_nHit
//                                          <<" "<<irClus->_mean<<" "<<irClus->_sigma<<endl;
//
//                  }
                  ///-----------------------------------

                  if (rcr_it->second.empty()) continue;
                  rwclstvec::iterator cr_it=rcr_it->second.begin();
                  while (cr_it!=rcr_it->second.end()){
                  //if ( clustersList.empty() ) {
                          startingRowClust=*cr_it;
                          cout<<"New cluster start at row "<<rcr_it->first<<" with cluster of row bin "<<startingRowClust->_firstStationID<<" "<<startingRowClust->_lastStationID<<endl;
                          clustersList.push_back( Clust(  rcr_it->first, startingRowClust ) );
                          rcr_it->second.erase(cr_it);
                          findCluster(rcr_it->first, startingRowClust);
                          if (!rcr_it->second.empty()) cr_it=rcr_it->second.begin();
                  }

//                  //if ( .... )
//



//                  stbrel::iterator rcr_it=rowClst_Column_rel.find( tmpRowId );
//                  while( rcr_it!=rowClst_Column_rel.end() ){
//                          tmpRwClust=&(rowClusts.at(rcr_it->second));
//                          if ( (tmpRwClust->_firstStationID-startingRowClust->_lastStationID)<=1 &&
//                               (startingRowClust->_firstStationID-tmpRwClust->_lastStationID)<=1 ){
//
//                          }
//                          rcr_it=rowClst_Column_rel.find( tmpRowId );
//                  }
//                  rowClst_Column_rel
//
//                  tmpRwClust=&(rowClusts.at());

          }




          for ( stbrel::iterator rcc_it=rowClst_Column_rel.begin(); rcc_it!=rowClst_Column_rel.end(); ++rcc_it ) {
                  cout<<"i-th row "<<rcc_it->first<<" j-th clusters "<<rcc_it->second<<endl;
                  rwClustPtr &irClus = rowClusts.at(rcc_it->second);
                  cout<<"\t cluster data: "<<irClus->_firstStationID<<" "<<irClus->_lastStationID<<" "<<irClus->_nHit
                                  <<" "<<irClus->_mean<<" "<<irClus->_sigma<<endl;
          }
          for ( ptcmaprowclsrel::const_iterator pmrcc_it=ptcMap_rowCls_rel.begin(); pmrcc_it!=ptcMap_rowCls_rel.end(); ++pmrcc_it ) {
                  cout<<"i-th row "<<pmrcc_it->first<<endl;
                  for ( ptcclsbrel::const_iterator prcc_it=pmrcc_it->second.begin(); prcc_it!=pmrcc_it->second.end(); ++prcc_it ) {
                          cout<<"\t pitch: "<<prcc_it->first<<" rowClust "<<prcc_it->second.first<<" from bin "<<prcc_it->second.second<<endl;
                  }

                  cout<<"Row Clust pairs"<<endl;
                  rwavptcclscpl::iterator ptcrwClrwClmap_it = rwMap_avPtc_rwClrwCl.find(pmrcc_it->first);
                  for ( ptcrwClrwCl_it=ptcrwClrwClmap_it->second.begin(); ptcrwClrwCl_it!=ptcrwClrwClmap_it->second.end(); ++ptcrwClrwCl_it ){
                          cout<<"\t row Clusters pair "<<ptcrwClrwCl_it->first.first<<" - "<<ptcrwClrwCl_it->first.second<<endl;
                          cout<<"\t\t available pitches: ";
                          for ( std::vector<unsigned short>::iterator ptcList=ptcrwClrwCl_it->second.begin(); ptcList!=ptcrwClrwCl_it->second.end(); ++ptcList ){
                                  cout<<*ptcList<<" ";
                          }
                          cout<<endl;

                  }

          }


          for (int i=0; i<nYbin; i++) { for(int j=0; j<nXbin; j++) cout<<dataArr[i][j]<<" "; cout<<endl;}
//          cout<<"reverse"<<endl;
//          for (int i=0; i<nYbin; i++) { for(int j=0; j<nXbin; j++) cout<<dataArrRev[i][j]<<" "; cout<<endl;}
          for ( ptcbnrel::iterator sttnptc_it=sttn_sectPitch_rel.begin(); sttnptc_it!=sttn_sectPitch_rel.end(); ++sttnptc_it ){
                  cout<<"Row "<<sttnptc_it->first<<" : ";
                  for ( ptcsctrrel::iterator ptcst_it=sttnptc_it->second.begin(); ptcst_it!=sttnptc_it->second.end(); ++ptcst_it ) {
                          cout<<ptcst_it->second<<" ";
                  }
                  cout<<endl;
          }
          cout<<"--- without station repetitions ---"<<endl;
          for ( goodsctrstatnrel::iterator sttnptc_it=good_sttn_sctr_rel.begin(); sttnptc_it!=good_sttn_sctr_rel.end(); ++sttnptc_it ){
                  cout<<"Row "<<sttnptc_it->first<<" : ";
                  for ( goodsttnrel::iterator ptcst_it=sttnptc_it->second.begin(); ptcst_it!=sttnptc_it->second.end(); ++ptcst_it ) {
                          cout<<*ptcst_it<<" ";
                  }
                  cout<<endl;
          }





//          for ( rwavptcclscpl::iterator ptcrwClrwClmap_it=rwMap_avPtc_rwClrwCl.begin(); ptcrwClrwClmap_it!=rwMap_avPtc_rwClrwCl.end(); ++ptcrwClrwClmap_it ) {
//                  for ( ptcrwClrwCl_it=ptcrwClrwClmap_it->second.begin(); ptcrwClrwCl_it!=ptcrwClrwClmap_it->second.end(); ++ptcrwClrwCl_it ) {
//                          if ( clustersList.empty() ) {
//                                  assgnRwClus.push_back(ptcrwClrwCl_it->first.first);
//                                  startingRowClust=&(rowClusts.at(ptcrwClrwCl_it->first.first));
//                                  clustersList.push_back( Clust(  ptcrwClrwClmap_it->
//                                                  first, startingRowClust ) );
//                          }
//                          //if ( .... )
//
//                          tmpRowId=ptcrwClrwCl_it->first.first;
//                          tmpRowId++;
//                          stbrel::iterator rcc_it=rowClst_Column_rel.find( tmpRowId );
//                          while( rcc_it!=rowClst_Column_rel.end() ){
//                                  tmpRwClust=&(rowClusts.at(rcc_it->second));
//                                  if ( (tmpRwClust->_firstStationID-startingRowClust->_lastStationID)<=1 &&
//                                       (startingRowClust->_firstStationID-tmpRwClust->_lastStationID)<=1 ){
//
//                                  }
//                                  rcc_it=rowClst_Column_rel.find( tmpRowId );
//                          }
//                          rowClst_Column_rel
//
//                          tmpRwClust=&(rowClusts.at());
//
//                  }
//          }


          int iCl=0;
          for ( std::vector<Clust>::iterator clstlst_it=clustersList.begin(); clstlst_it!=clustersList.end(); ++clstlst_it ) {
                 cout<<"Cluster "<< iCl<<endl;
                 for ( rwclincl::iterator rwclincl_it=clstlst_it->_rClusts.begin(); rwclincl_it!=clstlst_it->_rClusts.end(); ++rwclincl_it) {
                         cout<<"\t row "<<rwclincl_it->first<<" rwclst bin "<<rwclincl_it->second->_firstStationID<<" "<<rwclincl_it->second->_lastStationID<<endl;

                 }
                 iCl++;
          }

          std::vector<Clust>::iterator clstlst_it=clustersList.begin();
          std::vector<Clust>::iterator tmp_clstlst_it;
          short tmpMaxSect_1, tmpMinSect_1, tmpMaxSect_2, tmpMinSect_2;
          bool nexClust;
          while ( clstlst_it!=clustersList.end() ) {
                  nexClust=true;
                  tmp_clstlst_it=clstlst_it;
                  tmp_clstlst_it++;
                  tmpMaxSect_1=clstlst_it->_lastSectorID;
                  tmpMinSect_1=clstlst_it->_firstSectorID;
                  if ( tmpMaxSect_1>=12 ) {
                          tmpMaxSect_1-=12;
                          tmpMinSect_1-=12;
                  }
                  for ( ; tmp_clstlst_it!=clustersList.end(); ++tmp_clstlst_it) {
                          tmpMaxSect_2=tmp_clstlst_it->_lastSectorID;
                          tmpMinSect_2=tmp_clstlst_it->_firstSectorID;
                          if ( tmpMaxSect_2>=12 ) {
                                  tmpMaxSect_2-=12;
                                  tmpMinSect_2-=12;
                          }
                          if ( ( tmpMaxSect_1>=tmpMaxSect_2 && tmpMinSect_1<=tmpMinSect_2 ) &&
                               ( clstlst_it->_minStationID<=tmp_clstlst_it->_minStationID && clstlst_it->_maxStationID>=tmp_clstlst_it->_maxStationID ) ) {
                                  clustersList.erase(tmp_clstlst_it);
                                  break;
                          }
                          if ( ( tmpMaxSect_2>=tmpMaxSect_1 && tmpMinSect_2<=tmpMinSect_1 ) &&
                               ( tmp_clstlst_it->_minStationID<=clstlst_it->_minStationID && tmp_clstlst_it->_maxStationID>=clstlst_it->_maxStationID ) ) {
                                  clustersList.erase(clstlst_it);
                                  nexClust=false;
                                  break;
                          }
                  }
                  if ( nexClust ) clstlst_it++;
                  //clstlst_it=clustersList.begin();
          }

          cout<<"------------ double cluster removed ------------"<<endl;
          iCl=0;
          for ( std::vector<Clust>::iterator clstlst_it=clustersList.begin(); clstlst_it!=clustersList.end(); ++clstlst_it ) {
                 cout<<"Cluster "<< iCl <<" Station mean "<<clstlst_it->_mean_Sttn<<" sigma "<<clstlst_it->_sigma_Sttn<<
                                 " Sector mean "<<clstlst_it->_mean_Sctr<<" sigma "<<clstlst_it->_sigma_Sctr<<endl;
                 cout<<"Cluster m "<<clstlst_it->_m<<" q "<<clstlst_it->_q<<" errm "<<clstlst_it->_errm<<" errq "<<clstlst_it->_errq<<endl;
                 for ( rwclincl::iterator rwclincl_it=clstlst_it->_rClusts.begin(); rwclincl_it!=clstlst_it->_rClusts.end(); ++rwclincl_it) {
                         cout<<"\t row "<<rwclincl_it->first<<" rwclst bin "<<rwclincl_it->second->_firstStationID<<" "<<rwclincl_it->second->_lastStationID<<
                                         " mean "<<rwclincl_it->second->_mean<<" sigma "<<rwclincl_it->second->_sigma<<endl;
                 }
                 iCl++;
          }


  }

  void TTDisplayData::findCluster(unsigned short tmpRowId, rwClustPtr &startingRowClust){

          rwClustPtr tmpRwClust;
          rwclinrwrel::iterator tmp_rcr_it;

          tmpRowId++;
          tmp_rcr_it=rwClst_forRw_rel.find(tmpRowId);
          if (tmp_rcr_it!=rwClst_forRw_rel.end()) {
                  //cout<<"Up) Lokking for row cluster to add in row "<<tmpRowId<<endl;
                  if (!tmp_rcr_it->second.empty()) {
                          for ( rwclstvec::iterator tmp_cr_it=tmp_rcr_it->second.begin(); tmp_cr_it!=tmp_rcr_it->second.end(); ++tmp_cr_it) {
                                  tmpRwClust=*tmp_cr_it;
                                  if ( (tmpRwClust->_firstStationID-startingRowClust->_lastStationID)<=1 &&
                                                  (startingRowClust->_firstStationID-tmpRwClust->_lastStationID)<=1 ) {
                                          clustersList.back().addRwClust(tmpRowId, tmpRwClust);
                                          tmp_rcr_it->second.erase(tmp_cr_it);
                                          findCluster(tmpRowId,tmpRwClust);
                                          break;
                                  }
                          }
                  }
          }
          tmpRowId--; tmpRowId--;
          tmp_rcr_it=rwClst_forRw_rel.find(tmpRowId);
          if (tmp_rcr_it!=rwClst_forRw_rel.end()) {
                  //cout<<"Down) Lokking for row cluster to add in row "<<tmpRowId<<endl;
                  if (!tmp_rcr_it->second.empty()) {
                          for ( rwclstvec::iterator tmp_cr_it=tmp_rcr_it->second.begin(); tmp_cr_it!=tmp_rcr_it->second.end(); ++tmp_cr_it) {
                                  tmpRwClust=*tmp_cr_it;
                                  if ( (tmpRwClust->_firstStationID-startingRowClust->_lastStationID)<=1 &&
                                                  (startingRowClust->_firstStationID-tmpRwClust->_lastStationID)<=1 ) {
                                          clustersList.back().addRwClust(tmpRowId, tmpRwClust);
                                          tmp_rcr_it->second.erase(tmp_cr_it);
                                          findCluster(tmpRowId,tmpRwClust);
                                          break;
                                  }
                          }
                  }
          }
  }


}  // end namespace mu2e

using mu2e::TTDisplayData;
DEFINE_ART_MODULE(TTDisplayData);
