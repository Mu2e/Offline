//
// Fast Patter recognition bck rejection algorithm based on geometry considerations
//
// $Id: BkgTrackRejecterByGeomTplg_module.cc,v 1.4 2011/07/14 16:38:54 tassiell Exp $
// $Author: tassiell $
// $Date: 2011/07/14 16:38:54 $
//
// Original author G. Tassielli
//

// C++ includes.
#include <iostream>
//#include <string>
#include <memory>
#include <map>
#include <utility>
#include <cstring>
#include <algorithm>
#include <limits>
#include <ctime>

#include <boost/shared_ptr.hpp>

// Framework includes.
#include "art/Framework/Core/EDProducer.h"
//#include "art/Framework/Core/EDAnalyzer.h"
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
//#include "MCDataProducts/inc/GenParticleCollection.hh"
//#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
//#include "MCDataProducts/inc/SimParticleCollection.hh"
//#include "MCDataProducts/inc/StepPointMCCollection.hh"
//#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "ITrackerGeom/inc/ITracker.hh"
#include "TTrackerGeom/inc/TTracker.hh"
//#include "FastPatternReco/inc/TTHitPerTrackData.hh"
//#include "FastPatternReco/inc/GenTrackData.hh"
#include "RecoDataProducts/inc/TrackerHitTimeCluster.hh"
#include "RecoDataProducts/inc/TrackerHitTimeClusterCollection.hh"
#include "RecoDataProducts/inc/SectorStationCluster.hh"
#include "RecoDataProducts/inc/SctrSttnClusterGroup.hh"
//#include "RecoDataProducts/inc/SectorStationClusterCollection.hh"
#include "RecoDataProducts/inc/SctrSttnClusterGroupCollection.hh"

// Root includes.
#include "TApplication.h"
#include "TCanvas.h"
//#include "TDirectory.h"
#include "TMath.h"
#include "TH1I.h"
#include "TH1F.h"
#include "TH2I.h"
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
using art::Event;

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

typedef art::Ptr<StrawHit> StrawHitPtr;
typedef std::multimap<unsigned int, StrawHitPtr, less<unsigned int> > stMaprel;
typedef art::Ptr<TrackerHitTimeCluster> TrackerHitTimeClusterPtr;

//  class Straw;

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
                  else _sigma = 0.288675135;   // 1/sqrt(12)
          }

          bool operator==( rowClust const& comp) const {
                  return ( _firstStationID==comp._firstStationID &&
                                  _lastStationID==comp._lastStationID &&
                                  _nHit==comp._nHit );
          }

  protected:
          float tmpData;
          float tmpMean;
          float tmpMS;
  };

  typedef boost::shared_ptr<rowClust> rwClustPtr;
  typedef std::vector< rwClustPtr > rwclstvec;
  typedef std::map<unsigned short, rwclstvec> rwclinrwrel;                                                         //"row Clusters" for each row contained
  typedef std::multimap<unsigned short, rwClustPtr, less<unsigned short> > rwclincl;                                                   //"row Clusters" for each row contained in Cluster


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
                  _mean_Sctr(0.00000),
                  _MS_Sctr(0.00000),
                  _sigma_Sctr(0.00000),
                  _m(0.00000),
                  _q(0.00000),
                  _errm(0.00000),
                  _errq(0.00000),
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
                  else {
                          _sigma_Sttn=0.288675135;   // 1/sqrt(12)
                  }
//                  else {
//                          _m=std::numeric_limits<float>::infinity();
//                          _q=(float)_maxStationID + 0.50000;
//                  }

                  _q=tmpDataY;


          }
          void addRwClust( unsigned short iRow, rwClustPtr &iRwclust ) {
                  //cout<<"\t start adding I row "<<iRow<<" to cluster"<<endl;
                  //cout<<"\t previous dim of _rClusts "<<_rClusts.size()<<endl;
                  _rClusts.insert( rwclincl::value_type( iRow, iRwclust ) );
                  //cout<<"\t --------------------- 1 ---------------------"<<endl;
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
                          _sigma_Sttn=0.288675135;   // 1/sqrt(12)
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
                                          //tmpSigmaY=sqrt( tmpNsCovY/(_nHit*(_nHit-2)) );
                                          tmpSigmaY=sqrt( ((tmpNsCovY-pow(_m,2)*tmpNsCovX)/(_nHit*_nHit-2)) );
                                          _errm=tmpSigmaY/sqrt( tmpNsCovX/_nHit );
                                          //_errq=tmpSigmaY*sqrt( 1.00000/_nHit+pow(tmpMeanX,2)/(_nHit*tmpNsCovX) );
                                          _errq=tmpSigmaY*tmpMSX/tmpNsCovX;
                                  }
                          }
                  }
                  //cout<<"\t I row "<<iRow<<" added to cluster"<<endl;
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


  class BkgTrackRejecterByGeomTplg : public art::EDProducer {
//  class BkgTrackRejecterByGeomTplg : public art::EDAnalyzer {

  public:

    explicit BkgTrackRejecterByGeomTplg(fhicl::ParameterSet const& pset);
    virtual ~BkgTrackRejecterByGeomTplg() {
//            _peakFinder->Delete();
//            _fg->Delete();
//            if (_peakFinder!=0x0)        delete _peakFinder;
//            if (_fg!=0x0)                delete _fg;
            if (_fakeCanvas!=0x0)        delete _fakeCanvas;
            if (_peaksCanvHistos!=0x0)   delete _peaksCanvHistos;
            if (_hPkStDistTrs!=0x0)      delete _hPkStDistTrs;
//            if (_hPkStDistanceTrs!=0x0)  delete _hPkStDistanceTrs;
            if (_hPkStDist2W!=0x0)       delete _hPkStDist2W;
            if (_hPkStClustDist2W!=0x0)  delete _hPkStClustDist2W;
            if (_hPkSecDist2W!=0x0)      delete _hPkSecDist2W;
            if (_hPkSecClustDist2W!=0x0) delete _hPkSecClustDist2W;
            if (_hPkSecVsStationDist2W!=0x0)      delete _hPkSecVsStationDist2W;
            if (_hPkSecVsStationDist2WGood!=0x0)      delete _hPkSecVsStationDist2WGood;
            if (_hPkAccArrSecStation!=0x0)      delete _hPkAccArrSecStation;
    }

    virtual void beginJob();
    //void endJob();

    // This is called for each event.
    void produce(art::Event & e);
    //void analyze(art::Event const& e);

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

    // Label of the module that made the hits.
    std::string _timeRejecterModuleLabel;

    // Use n times sigma for matching the clusters in the Sector plane
    float _nSigmaForClMatch;

    // Use n times sigma for matching the clusters slope with 0
    float _nSigmaForClSlope;

    // Label of the generator.
    std::string _generatorModuleLabel;

    bool _doDisplay;

//    // Number of events to accumulate between prompts.
//    int _nAccumulate;

    // End: run time parameters

    // Pointers to histograms, ntuples, TGraphs.

    TObjArray*    _peaksCanvHistos;
    TObjArray*    _hPkStDistTrs;
    TObjArray*    _hPkStDist2W;
    TObjArray*    _hPkStClustDist2W;
    TObjArray*    _hPkSecDist2W;
    TObjArray*    _hPkSecClustDist2W;
    TObjArray*    _hPkSecVsStationDist2W;
    TObjArray*    _hPkSecVsStationDist2WGood;
    TObjArray*    _hPkAccArrSecStation;

    TH1I*         _hClockCicles;
    TH1F*         _hExecTime;

//    TF1*          _fg;
//    TCanvas*      _cnvForPeakstudy;

    TCanvas*      _fakeCanvas;

//    TSpectrum *   _peakFinder;

    int   ntimeBin;
    float maxTimeHist; //ns
    float timeBinDim;  //ns

    inline unsigned int iRot(int device, int sector){
            return (unsigned int) (device%2)+2*sector;
    }

    void SctrSttnMapAnalyze(TH2I *startDist, TH2I *votingArray, int minPitch=3, int maxPitch=11);
    rwclinrwrel rwClst_forRw_rel;
    std::vector<Clust> clustersList;
    void findCluster(unsigned short tmpRowId, rwClustPtr &startingRowClust);

    // Some ugly but necessary ROOT related bookkeeping:

    // The job needs exactly one instance of TApplication.  See note 1.
    auto_ptr<TApplication> _application;

    // Save directory from beginJob so that we can go there in endJob. See note 3.
//    TDirectory* _directory;

  };

  BkgTrackRejecterByGeomTplg::BkgTrackRejecterByGeomTplg(fhicl::ParameterSet const& pset) :

    // Run time parameters
    _moduleLabel(pset.get<string>("module_label")),/*@module_label*/
    _g4ModuleLabel(pset.get<string>("g4ModuleLabel")),
    _trackerStepPoints(pset.get<string>("trackerStepPoints","tracker")),
    _makerModuleLabel(pset.get<string>("makerModuleLabel")),
    _timeRejecterModuleLabel(pset.get<string>("tRejecterModuleLabel")),
    _nSigmaForClMatch(pset.get<float>("nSigmaForClMatch")),
    _nSigmaForClSlope(pset.get<float>("nSigmaForClSlope")),
    _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
   /*_nAccumulate(pset.get<int>("nAccumulate",20)),*/
    _doDisplay(pset.get<bool>("doDisplay",false)),

    // ROOT objects that are the main focus of this example.
    _peaksCanvHistos(0),
    _hPkStDistTrs(0),
    _hPkStDist2W(0),
    _hPkStClustDist2W(0),
    _hPkSecDist2W(0),
    _hPkSecClustDist2W(0),
    _hPkSecVsStationDist2W(0),
    _hPkSecVsStationDist2WGood(0),
    _hPkAccArrSecStation(0),
    _hClockCicles(0),
    _hExecTime(0),

//    _fg(0),
//    _cnvForPeakstudy(0),

    _fakeCanvas(0),

//    _peakFinder(0),
    ntimeBin(0),
    maxTimeHist(2500.0),
    timeBinDim(10.0),

    // Some ugly but necessary ROOT related bookkeeping.
    _application(0)
//    _directory(0)
  {
          cout<<"Constructed"<<endl;
          if ( _nSigmaForClMatch<0.0 ) _nSigmaForClMatch*=-1.00000;
          if ( _nSigmaForClSlope<0.0 ) _nSigmaForClSlope*=-1.00000;
          // Tell the framework what we make.
          produces<SctrSttnClusterGroupCollection>();

  }

  void BkgTrackRejecterByGeomTplg::beginJob(){

    cout<<"Bkg rejection by Geom job!"<<endl;

    // Get access to the TFile service and save current directory for later use.
    art::ServiceHandle<art::TFileService> tfs;
    _hClockCicles       = tfs->make<TH1I>( "hClockCicles",   "N clock cicles needed to analyze one Event by BkgTrackRejecterByGeomTplg", 2000, 5000, 15000  );
    _hExecTime          = tfs->make<TH1F>( "hExecTime",   "Execution time to analyze one Event by BkgTrackRejecterByTime", 1000, 0.0, 1.0  );


    ntimeBin    = (int) maxTimeHist/timeBinDim;

    // Create a histogram.

    _hPkStDistTrs      = new TObjArray();
    _hPkStDist2W       = new TObjArray();
    _hPkStClustDist2W  = new TObjArray();
    _hPkSecDist2W      = new TObjArray();
    _hPkSecClustDist2W = new TObjArray();
    _hPkSecVsStationDist2W      = new TObjArray();
    _hPkSecVsStationDist2WGood      = new TObjArray();
    _hPkAccArrSecStation      = new TObjArray();

    if(_doDisplay) {
            // If needed, create the ROOT interactive environment. See note 1.
            if ( !gApplication ){
                    int    tmp_argc(0);
                    char** tmp_argv(0);
                    _application = auto_ptr<TApplication>(new TApplication( "noapplication", &tmp_argc, tmp_argv ));
            }

            gStyle->SetPalette(1);
            gROOT->SetStyle("Plain");


            _peaksCanvHistos   = new TObjArray();

            _fakeCanvas = new TCanvas("canvas_Fake","double click for next event",300,100);

    }

    //_peakFinder = new TSpectrum(20);

    // See note 3.
//    _directory = gDirectory;

  }

  void BkgTrackRejecterByGeomTplg::produce(art::Event & event ) {
//  void BkgTrackRejecterByGeomTplg::analyze(art::Event const& event ) {

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


    //_peakFinder->Clear();
    //_peaksCanvases->Clear();
    typedef boost::shared_ptr<TArrow> clssegdr;
    std::vector< std::vector< clssegdr > > clstSegments;
    TF1 clustSegline("clustSegline","pol1", 0.0, 20.0);

    //auto_ptr<SectorStationClusterCollection> sscc(new SectorStationClusterCollection);
    auto_ptr<SctrSttnClusterGroupCollection> sscc(new SctrSttnClusterGroupCollection);


    const Tracker& tracker = getTrackerOrThrow();
    const TTracker &ttr = static_cast<const TTracker&>( tracker );
    const std::vector<Device> ttrdev = ttr.getDevices();

    stMaprel ssmap_Bin_Straw_rel;           //relation between straw and the AbsSector (Rotation) vs Station map
    stMaprel::iterator ssmpbnstrw_rel_it;   //iterator for relation between straw and the AbsSector (Rotation) vs Station map
    int tmpiGeomBin, tmpiGeomBinRow;

    art::Handle<StrawHitCollection> pdataHandle;
    event.getByLabel(_makerModuleLabel,pdataHandle);
    //StrawHitCollection const* hits = pdataHandle.product();

    art::Handle<TrackerHitTimeClusterCollection> tclustHandle;
    event.getByLabel(_timeRejecterModuleLabel,tclustHandle);
    TrackerHitTimeClusterCollection const* tclusts = tclustHandle.product();


    clock_t startClock = clock();


    size_t nTimeClusPerEvent = tclusts->size();
    cout<<"----------------------------------------------------------------------------------"<<endl;
    cout<<"event "<<event.id().event()<<" N peak found "<<nTimeClusPerEvent<<endl;
    cout<<"----------------------------------------------------------------------------------"<<endl;
    int i1peak;

    StrawId sid;
    int stn, layern, devicen, sectorn, stationn;
    unsigned int absSect, devStId;
    unsigned int ihit;
    bool storeNewClstrsGrpFortPeak;

    int nMapXBin;
    bool storeClusterHit, skip, recheckprevious, addingBadCluster;
    unsigned int icl;
    int iclgr, groupStoredInEv;
    typedef std::multimap< int, unsigned int, less<int>  > clugrouprel;
    clugrouprel storedClutID;
    clugrouprel::iterator lastsCID_it;

    float clSctrDist=0.0, errClSctrDist=0.0;
    float meanPitch=0.0, errPitch=0.0, invSigma2Pitch=0.0;
    float tmpPitch=0.0, tmpPitchSigma2=0.0, tmpInvNloop=0.0;
    bool clMcompatibleWith0, hugeCluster;

    SctrSttnClusterGroupCollection::iterator sscc_it;

    TH2I *tmpStClustDist    = new TH2I("tmpStClustDist","tmp Smoothing of Device vs Straw multiplicity Distribution",36,0,36,1000,-200,800);
    TH2I *tmpSecClustDist   = new TH2I("tmpSecClustDist","tmp Smoothing of Device vs Sector multiplicity Distribution",36,0,36,20,-4,16);
    TH2I *ihPkStDistTrs, *ihPkStDist2W, *ihPkStClustDist2W, *ihPkSecDist2W, *ihPkSecClustDist2W, *ihPkSecVsStationDist2WG, *ihPkSecVsStationDist2W, *ihPkAccArrSecStation;
    TCanvas * iCanvHist=0x0;

    for (size_t ipeak=0; ipeak<nTimeClusPerEvent; ipeak++) {
    	TrackerHitTimeCluster const&  tclust(tclusts->at(ipeak));
    	cout<<"\t Hit in peak "<<tclust._selectedTrackerHits.size()<<endl;

    	i1peak=ipeak;
        ++i1peak;

        _hPkStDistTrs->AddAtAndExpand(new TH2I(Form("h%d-Pk_StDistTrs",i1peak),Form("Sector-Straw Distribution in Transverse plane for the %d-th peak",i1peak),16,0,16,50,0,50),ipeak);
//                    _hPkStDistanceTrs->AddAtAndExpand(new TH2I(Form("h%d-Pk_StDistanceTrs",i1peak),Form("Sector-Straw distance Distribution in Transverse plane for the %d-th peak",i1peak),20,-4,16,50,0,50),ipeak);
        _hPkStDist2W->AddAtAndExpand(new TH2I(Form("h%d-Pk_StDist2W",i1peak),Form("Device vs Straw multiplicity Distribution for the %d-th peak",i1peak),36,0,36,1000,-200,800),ipeak);
        _hPkStClustDist2W->AddAtAndExpand(new TH2I(Form("h%d-Pk_StClustDist2W",i1peak),Form("Smoothing of Device vs Straw multiplicity Distribution for the %d-th peak",i1peak),36,0,36,1000,-200,800),ipeak);
        _hPkSecDist2W->AddAtAndExpand(new TH2I(Form("h%d-Pk_SecDist2W",i1peak),Form("Device vs Sector multiplicity Distribution for the %d-th peak",i1peak),36,0,36,16,0,16),ipeak);
        _hPkSecClustDist2W->AddAtAndExpand(new TH2I(Form("h%d-Pk_SecClustDist2W",i1peak),Form("Smoothing of Device vs Sector multiplicity Distribution for the %d-th peak",i1peak),36,0,36,16,0,16),ipeak);
        _hPkSecVsStationDist2W->AddAtAndExpand(new TH2I(Form("h%d-Pk_SecVsStationDist2W",i1peak),Form("Station vs Sector multiplicity Distribution for the %d-th peak",i1peak),18,0,18,16,0,16),ipeak);
        _hPkSecVsStationDist2WGood->AddAtAndExpand(new TH2I(Form("h%d-Pk_SecDist2WG",i1peak),Form("Selected Station vs Sector multiplicity Distribution for the %d-th peak",i1peak),18,0,18,28,0,28),ipeak);
        _hPkAccArrSecStation->AddAtAndExpand(new TH2I(Form("h%d-Pk_AccArrSecStation",i1peak),Form("Accumulator Array for Station distance vs Sector for the %d-th peak",i1peak),9,3,12,32,-16,16),ipeak);

        ihPkStDistTrs      = ((TH2I *) _hPkStDistTrs->At(ipeak));
        ihPkStDist2W       = ((TH2I *) _hPkStDist2W->At(ipeak));
        ihPkStClustDist2W  = ((TH2I *) _hPkStClustDist2W->At(ipeak));
        ihPkSecDist2W      = ((TH2I *) _hPkSecDist2W->At(ipeak));
        ihPkSecClustDist2W = ((TH2I *) _hPkSecClustDist2W->At(ipeak));

        ihPkSecVsStationDist2W  = ((TH2I *) _hPkSecVsStationDist2W->At(ipeak));
        ihPkSecVsStationDist2WG = ((TH2I *) _hPkSecVsStationDist2WGood->At(ipeak));
        ihPkAccArrSecStation    = ((TH2I *) _hPkAccArrSecStation->At(ipeak));

        nMapXBin = ihPkSecVsStationDist2W->GetNbinsX();

        tmpStClustDist->Reset();
        tmpSecClustDist->Reset();


        ssmap_Bin_Straw_rel.clear();
        ihit=0;
        for (std::vector<StrawHitPtr>::const_iterator iTCHit=tclust._selectedTrackerHits.begin(); iTCHit!=tclust._selectedTrackerHits.end(); ++iTCHit){
                // Access data
//                cout<<"\t\t iHit in peak at "<<*iTCHit<<endl;
                //StrawHit        const&      hit(hits->at(*iTCHit));
                StrawHit        const&      hit=*(*iTCHit);
                StrawIndex si = hit.strawIndex();
                const Straw & str = tracker.getStraw(si);

                // cout << "Getting informations about cells" << endl;

                sid = str.id();
                stn     = sid.getStraw();
                layern  = sid.getLayer();
                devicen = sid.getDevice();
                sectorn = sid.getSector();

                absSect = iRot(devicen,sectorn);
                devStId = absSect*50+stn;

                //ihPkStDistTrs->Fill(devStId);
                ihPkStDistTrs->Fill(absSect,stn);
                //if (absSect>=8 && absSect<=11)  ihPkStDistTrs->Fill(absSect-12,stn);
                if (absSect>=0 && absSect<=3)   ihPkStDistTrs->Fill(absSect+12,stn);

                ihPkSecDist2W->Fill(devicen,absSect);
                //if (absSect>=8 && absSect<=11)  ihPkSecDist2W->Fill(devicen,absSect-12);
                if (absSect>=0 && absSect<=3)   ihPkSecDist2W->Fill(devicen,absSect+12);

                ihPkStDist2W->Fill(devicen,absSect*50+stn);
                //if (absSect>=8 && absSect<=11)  ihPkStDist2W->Fill(devicen,(absSect-12)*50+stn);
                if (absSect>=0 && absSect<=3)   ihPkStDist2W->Fill(devicen,(absSect+12)*50+stn);
//                ihPkStDist2W->Fill(devicen,devStId);

                stationn=devicen/2;
                tmpiGeomBin = ihPkSecVsStationDist2W->Fill(stationn,absSect);
                ssmap_Bin_Straw_rel.insert( stMaprel::value_type((unsigned int) tmpiGeomBin, *iTCHit/*StrawHitPtr( pdataHandle, ihit)*/ ) );
                //if (absSect>=8 && absSect<=11)  ihPkSecVsStationDist2W->Fill(stationn,absSect-12);
                if (absSect>=0 && absSect<=3)   {
                        tmpiGeomBin = ihPkSecVsStationDist2W->Fill(stationn,absSect+12);
                        ssmap_Bin_Straw_rel.insert( stMaprel::value_type((unsigned int) tmpiGeomBin, *iTCHit/*ihit*/) );
                }

                ihit++;

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

//                int nYbingroup=5;//5*50
//                float Yhit=0.0;
//                int YbinHalf = (int) nYbingroup/2;
//                int jYbin;
//                int nXBin=ihPkStDist2W->GetNbinsX();
//                int nYBin=ihPkStDist2W->GetNbinsY();
//                for (int itX=1; itX<=nXBin; itX++) {
//                        for (int itY=1; itY<=nYBin; itY++) {
//                                Yhit=0.0;
//                                for (int is=0; is<nYbingroup; is++) {
//                                        jYbin = itY -YbinHalf +is;
//                                        if (jYbin<1) continue;
//                                        if (jYbin>nYBin) break;
//                                        Yhit+=ihPkSecDist2W->GetBinContent(itX,jYbin);
//                                }
//                                tmpSecClustDist->SetBinContent(itX,itY,Yhit);
//                        }
//                }
//                int nXbingroup=9;
//                float Xhit=0.0;
//                int XbinHalf = (int) nXbingroup/2;
//                int jXbin;
//                for (int itY=1; itY<=nYBin; itY++) {
//                        for (int itX=1; itX<=nXBin; itX++) {
//                                Xhit=0.0;
//                                for (int is=0; is<nXbingroup; is++) {
//                                        jXbin = itX -XbinHalf +is;
//                                        if (jXbin<1) continue;
//                                        if (jXbin>nXBin) break;
//                                        Xhit+=tmpSecClustDist->GetBinContent(jXbin,itY);
//                                }
//                                ihPkSecClustDist2W->SetBinContent(itX,itY,Xhit);
//                        }
//                }

        }

        cout<<"------- start research ------"<<endl;
        SctrSttnMapAnalyze(ihPkSecVsStationDist2W, ihPkAccArrSecStation);
        storeNewClstrsGrpFortPeak=false;
        icl=0;
        iclgr=-1;
        groupStoredInEv=sscc->size();
        storedClutID.clear();
        std::vector<Clust>::iterator clstlst_it=clustersList.begin();
        //for ( std::vector<Clust>::iterator clstlst_it=clustersList.begin(); clstlst_it!=clustersList.end(); ++clstlst_it ) {
        while ( icl<clustersList.size() ) {
                skip=false;
                recheckprevious=false;
                for ( clugrouprel::iterator sCID_it=storedClutID.begin(); sCID_it!=storedClutID.end(); ++sCID_it ) if (sCID_it->second==icl) {skip=true; break;}
                if ( skip ) {
                        icl++;
                        ++clstlst_it;
                        continue;
                }
                storeClusterHit=false;
                clMcompatibleWith0 = abs(clstlst_it->_m)<_nSigmaForClSlope*clstlst_it->_errm;
                hugeCluster = ( (clstlst_it->_maxStationID-clstlst_it->_minStationID)>12 && (clstlst_it->_lastSectorID-clstlst_it->_firstSectorID)>2 );
                if ( clstlst_it->_m>0.00000 || (clMcompatibleWith0/*clstlst_it->_m==0.00000*/ && clstlst_it->_lastSectorID>clstlst_it->_firstSectorID ) || hugeCluster ) {
                        /*if ( clstlst_it->_m==std::numeric_limits<float>::infinity() && clstlst_it->_lastSectorID==clstlst_it->_firstSectorID ) storeNewClstrsGrpFortPeak=false;
                        else*/ storeNewClstrsGrpFortPeak=true;
                        addingBadCluster=false;
                        //for ( sscc_it=sscc->begin(); sscc_it!=sscc->end(); ++sscc_it ) {

                        //cout<<" ---- storedClutID size "<<storedClutID.size()<<endl;
                        for ( clugrouprel::iterator sCID_it=storedClutID.begin(); sCID_it!=storedClutID.end(); ++sCID_it ) {
                                clSctrDist=clstlst_it->_mean_Sctr - clustersList[sCID_it->second]._mean_Sctr;
                                //                - ((clstlst_it->_lastSectorID>11) ? 12 : 0);
                                //clSctrDist-=(clustersList[sCID_it->second]._mean_Sctr - ((clustersList[sCID_it->second]._lastSectorID>11) ? 12 : 0) );
                                errClSctrDist=sqrt( pow(clstlst_it->_sigma_Sctr,2) + pow(clustersList[sCID_it->second]._sigma_Sctr,2) );
                                errClSctrDist*=_nSigmaForClMatch;
                                if ( clSctrDist>=-errClSctrDist && clSctrDist<=errClSctrDist ) {
                                        storeClusterHit=true;
                                        storeNewClstrsGrpFortPeak=false;
                                        lastsCID_it=storedClutID.insert( clugrouprel::value_type(sCID_it->first,icl) );
                                        break;
                                }
                                else if ( (clstlst_it->_lastSectorID>11) && ( (clSctrDist-12.00000)>=-errClSctrDist && (clSctrDist-12.00000)<=errClSctrDist )){
                                        storeClusterHit=true;
                                        storeNewClstrsGrpFortPeak=false;
                                        lastsCID_it=storedClutID.insert( clugrouprel::value_type(sCID_it->first,icl) );
                                        break;
                                }
                                else if ( (clustersList[sCID_it->second]._lastSectorID>11) && ( (clSctrDist+12.00000)>=-errClSctrDist && (clSctrDist+12.00000)<=errClSctrDist )){
                                        storeClusterHit=true;
                                        storeNewClstrsGrpFortPeak=false;
                                        lastsCID_it=storedClutID.insert( clugrouprel::value_type(sCID_it->first,icl) );
                                        break;
                                }

                        }
                        if (storeNewClstrsGrpFortPeak) {
                                storeClusterHit=true;
                                ++iclgr;
                                lastsCID_it=storedClutID.insert( clugrouprel::value_type(groupStoredInEv+iclgr,icl) );
                                //sscc->push_back(SectorStationCluster());
                                sscc->push_back(SctrSttnClusterGroup());
                                //cout<<"Before "<<endl;
                                sscc->back()._relatedTimeCluster = TrackerHitTimeClusterPtr( tclustHandle, ipeak );
                                if ( clMcompatibleWith0/*clstlst_it->_m==0.00000*/ || clstlst_it->_m==std::numeric_limits<float>::infinity() || hugeCluster ) addingBadCluster=true;
                                //cout<<"After "<<endl;
                                //storeNewClstrsGrpFortPeak=false;
                                recheckprevious=true;
                        }
                }
                else if ( clstlst_it->_nHit>1 ) {
                        for ( clugrouprel::iterator sCID_it=storedClutID.begin(); sCID_it!=storedClutID.end(); ++sCID_it ){
                                clSctrDist=clstlst_it->_mean_Sctr - clustersList[sCID_it->second]._mean_Sctr;
                                //                - ((clstlst_it->_lastSectorID>11) ? 12 : 0);
                                //clSctrDist-=(clustersList[sCID_it->second]._mean_Sctr - ((clustersList[sCID_it->second]._lastSectorID>11) ? 12 : 0) );
                                errClSctrDist=sqrt( pow(clstlst_it->_sigma_Sctr,2) + pow(clustersList[sCID_it->second]._sigma_Sctr,2) );
                                if ( clSctrDist>=-errClSctrDist && clSctrDist<=errClSctrDist ) {
                                        storeClusterHit=true;
                                        addingBadCluster=true;
                                        lastsCID_it=storedClutID.insert( clugrouprel::value_type(sCID_it->first,icl) );
                                        break;
                                }
                                else if ( (clstlst_it->_lastSectorID>11) && ( (clSctrDist-12.00000)>=-errClSctrDist && (clSctrDist-12.00000)<=errClSctrDist )){
                                        storeClusterHit=true;
                                        addingBadCluster=true;
                                        lastsCID_it=storedClutID.insert( clugrouprel::value_type(sCID_it->first,icl) );
                                        break;
                                }
                                else if ( (clustersList[sCID_it->second]._lastSectorID>11) && ( (clSctrDist+12.00000)>=-errClSctrDist && (clSctrDist+12.00000)<=errClSctrDist )){
                                        storeClusterHit=true;
                                        addingBadCluster=true;
                                        lastsCID_it=storedClutID.insert( clugrouprel::value_type(sCID_it->first,icl) );
                                        break;
                                }
                        }
                }
                //cout<<"----------------------- 0 -----------------------"<<endl;
                if (storeClusterHit){
                        //cout<<"Cluster "<<icl<<" with hit "<<clstlst_it->_nHit<<endl;
                        SctrSttnClusterGroup &clugrp = sscc->at(lastsCID_it->first);
                        clugrp._selectedClusters.push_back(
                                        SectorStationCluster(clstlst_it->_mean_Sttn, clstlst_it->_sigma_Sttn, clstlst_it->_mean_Sctr, clstlst_it->_sigma_Sctr,
                                                        clstlst_it->_m, clstlst_it->_q, clstlst_it->_errm, clstlst_it->_errq,
                                                        clstlst_it->_firstSectorID, clstlst_it->_lastSectorID, clstlst_it->_minStationID, clstlst_it->_maxStationID)
                                        );
                        if ( addingBadCluster ) {
                                if ( clugrp._coupling==SctrSttnClusterGroup::good ) clugrp._coupling=SctrSttnClusterGroup::mixed;
                                else if ( clugrp._coupling==SctrSttnClusterGroup::mixed ) clugrp._coupling=SctrSttnClusterGroup::bad;
                        }
                        else if ( clugrp._coupling==SctrSttnClusterGroup::mixed ) clugrp._coupling=SctrSttnClusterGroup::good;
                        for ( rwclincl::iterator rwclincl_it=clstlst_it->_rClusts.begin(); rwclincl_it!=clstlst_it->_rClusts.end(); ++rwclincl_it) {
                                //cout<<"i row "<<rwclincl_it->first<<" station "<<rwclincl_it->second->_firstStationID<<" - "<<rwclincl_it->second->_lastStationID<<endl;
                                ssmap_Bin_Straw_rel.begin();
                                tmpiGeomBinRow = ( ( (nMapXBin+2)*(rwclincl_it->first+1-((rwclincl_it->first>15)?12:0)) ) +1 );
                                tmpiGeomBin= tmpiGeomBinRow +  rwclincl_it->second->_firstStationID;
                                ssmpbnstrw_rel_it=ssmap_Bin_Straw_rel.find(tmpiGeomBin);
                                while ( ssmpbnstrw_rel_it->first<= (tmpiGeomBinRow +  rwclincl_it->second->_lastStationID) ) {
                                        //cout<<"I'm storing bin i-th "<<ssmpbnstrw_rel_it->first<<endl;
                                        clugrp._selectedClusters.back()._selectedTrackerHits.push_back( ssmpbnstrw_rel_it->second );
                                        ++ssmpbnstrw_rel_it;
                                        if ( ssmpbnstrw_rel_it==ssmap_Bin_Straw_rel.end() ) break;
                                }
                        }

                        //cout<<"A) in Time Peak Cluster Group: "<<lastsCID_it->first<<" - "<<clugrp<<endl;
                }
                if (recheckprevious){
                        icl=0;
                        clstlst_it=clustersList.begin();
                }else {
                  icl++;
                  ++clstlst_it;
                }
         }

        //cout<<"----------------------- 1 -----------------------"<<endl;

        if (!storedClutID.empty()) {
                std::vector<unsigned int> clinGrp;
                bool addCl;
                pair<clugrouprel::iterator,clugrouprel::iterator> sCgID_it;
                for ( int jclgr=0; jclgr<=iclgr; ++jclgr ) {
                        clinGrp.clear();
                        meanPitch=0.00000;
                        errPitch=0.00000;

                        sCgID_it=storedClutID.equal_range(groupStoredInEv+jclgr);
                        SctrSttnClusterGroup &clugrp = sscc->at(sCgID_it.first->first);
                        for ( clugrouprel::iterator sCID_it=sCgID_it.first; sCID_it!=sCgID_it.second; ++sCID_it ) {
                                addCl=true;
                                for ( std::vector<unsigned int>::iterator clinGrp_it=clinGrp.begin(); clinGrp_it!=clinGrp.end(); ++clinGrp_it ){
                                        if ( clustersList[sCID_it->second]._mean_Sttn<clustersList[*clinGrp_it]._mean_Sttn ) {
                                                clinGrp.insert(clinGrp_it,1,sCID_it->second);
                                                addCl=false;
                                                break;
                                        }
                                }//store the index of the clusters of a group ordered for cluster position in z plane
                                if (addCl) clinGrp.push_back(sCID_it->second);

                                cout<<"!!!!!!!!!!!!!!! Cluster saved at: min station "<<clustersList[sCID_it->second]._minStationID<<" max station "<<clustersList[sCID_it->second]._maxStationID<<" first sect "
                                                <<clustersList[sCID_it->second]._firstSectorID<<" last sect "<<clustersList[sCID_it->second]._lastSectorID<<endl;
                        }
                        if ( clinGrp.size()>1 ) {
                                for ( unsigned int i_cl=0; i_cl<(clinGrp.size()-1); i_cl++ ) {
                                        for ( unsigned int j_cl=i_cl+1; j_cl<clinGrp.size(); j_cl++ ) {
                                                tmpInvNloop=1.00000/((float) (j_cl-i_cl));
                                                tmpPitch=abs(clustersList[clinGrp[i_cl]]._mean_Sttn-clustersList[clinGrp[j_cl]]._mean_Sttn)*tmpInvNloop;
                                                tmpPitchSigma2= (pow(clustersList[clinGrp[i_cl]]._sigma_Sttn,2) + pow(clustersList[clinGrp[j_cl]]._sigma_Sttn,2)) * pow(tmpInvNloop,2);

                                                invSigma2Pitch=1.00000/tmpPitchSigma2;
                                                meanPitch+=invSigma2Pitch*tmpPitch;
                                                errPitch+=invSigma2Pitch;
                                        }
                                }
                                errPitch=1.00000/errPitch;
                                meanPitch*=errPitch;
                                errPitch=sqrt(errPitch);
                                clugrp._meanPitch=meanPitch;
                                clugrp._sigmaPitch=errPitch;
                                 }
                        else {
                                clugrp._meanPitch=clustersList[clinGrp[0]]._mean_Sttn;
                                clugrp._sigmaPitch=clustersList[clinGrp[0]]._sigma_Sttn;
                        }
                       //}
                        //cout<<"B) in Time Peak Cluster Group: "<<jclgr<<" - "<<clugrp<<endl;

                }
        }

        //cout<<"----------------------- 2 -----------------------"<<endl;

        if (_doDisplay) {
                _peaksCanvHistos->AddAtAndExpand(new TCanvas(Form("canvasForHistos_%d-th_peak",i1peak),Form("Histograms of the peak at %f ns",tclust._meanTime),1290,860),ipeak);
                ihPkStDistTrs->SetStats(kFALSE);
                ihPkStDistTrs->SetXTitle("absSect");
                ihPkStDistTrs->SetYTitle("Straw");
                ihPkStDistTrs->SetTitleOffset(1.2,"y");
                ihPkStDist2W->SetStats(kFALSE);
                ihPkStDist2W->SetXTitle("device");
                ihPkStDist2W->SetYTitle("absStraw");
                ihPkStDist2W->SetTitleOffset(1.35,"y");
                ihPkStClustDist2W->SetStats(kFALSE);
                ihPkStClustDist2W->SetXTitle("device");
                ihPkStClustDist2W->SetYTitle("absStraw");
                ihPkStClustDist2W->SetTitleOffset(1.35,"y");
                ihPkSecDist2W->SetStats(kFALSE);
                ihPkSecDist2W->SetXTitle("device");
                ihPkSecDist2W->SetYTitle("absSect");
                ihPkSecDist2W->SetTitleOffset(1.2,"y");
                ihPkSecClustDist2W->SetStats(kFALSE);
                ihPkSecClustDist2W->SetXTitle("device");
                ihPkSecClustDist2W->SetYTitle("absSect");
                ihPkSecClustDist2W->SetTitleOffset(1.2,"y");

                ihPkSecVsStationDist2W->SetStats(kFALSE);
                ihPkSecVsStationDist2W->SetXTitle("Station");
                ihPkSecVsStationDist2W->SetYTitle("absSect");
                ihPkSecVsStationDist2W->SetTitleOffset(1.2,"y");
                ihPkSecVsStationDist2WG->SetStats(kFALSE);
                ihPkSecVsStationDist2WG->SetXTitle("Station");
                ihPkSecVsStationDist2WG->SetYTitle("absSect");
                ihPkSecVsStationDist2WG->SetTitleOffset(1.2,"y");
                ihPkAccArrSecStation->SetStats(kFALSE);
                ihPkAccArrSecStation->SetXTitle("Station");
                ihPkAccArrSecStation->SetYTitle("absSect");
                ihPkAccArrSecStation->SetTitleOffset(1.2,"y");

                iCanvHist=((TCanvas *) _peaksCanvHistos->At(ipeak));
                iCanvHist->Divide(3,2);

                iCanvHist->cd(1);
                ihPkStDistTrs->Draw("col z");

//                iCanvHist->cd(2);
//                ihPkSecDist2W->Draw("col z");

//                iCanvHist->cd(3);
//                ihPkStDist2W->Draw("col z");

                iCanvHist->cd(2);//3
                iCanvHist->cd(2)->SetGrid();
                ihPkSecVsStationDist2W->Draw("col z");

                clstSegments.push_back( std::vector< clssegdr >() );
                icl=0;
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
                                                        ihPkSecVsStationDist2W->GetBinContent( ibin+1, rwclincl_it->first+1-((rwclincl_it->first>15) ? 12 : 0) ) );
                                }
                        }
                        icl++;
                }
                iCanvHist->cd(3);
                iCanvHist->cd(3)->SetGrid();
                ihPkSecVsStationDist2WG->Draw("col z");
                for ( std::vector< clssegdr >::iterator clsS_it=clstSegments.back().begin(); clsS_it!=clstSegments.back().end(); ++clsS_it) (*clsS_it)->Draw();

                //                iCanvHist->cd(3);
                //                tmpStClustDist->Draw("col z");

//                iCanvHist->cd(5);
//                ihPkSecClustDist2W->Draw("col z");

//                iCanvHist->cd(6);
//                ihPkStClustDist2W->Draw("col z");

                iCanvHist->cd(5);//6
                ihPkAccArrSecStation->Draw("col z");



                iCanvHist->Modified();
                iCanvHist->Update();
        }

    }

    event.put(sscc);

    //for (unsigned long i=0; i<1000000; i++) cout<<"lost time "<<endl;
    clock_t stopClock = clock();
    _hClockCicles->Fill((unsigned long)(stopClock-startClock));
    _hExecTime->Fill( (float)(stopClock-startClock)/((float) CLOCKS_PER_SEC ) );
    cout<<"-------- N clok to analyze 1 ev by BkgTrackRejecterByGeomTplg "<<stopClock-startClock<<" @ "<<CLOCKS_PER_SEC<<endl;

    if (_doDisplay) {
            cerr << "Double click in the canvas_Fake to continue:" ;
            _fakeCanvas->cd();
            TLatex *printEvN = new TLatex(0.15,0.4,Form("Current Event: %d",event.id().event()));
            printEvN->SetTextFont(62);
            printEvN->SetTextSizePixels(180);
            printEvN->Draw();
            _fakeCanvas->WaitPrimitive();
            cerr << endl;
            delete printEvN;
            _peaksCanvHistos->Delete();
    }

    _hPkStDistTrs->Delete();
//    _hPkStDistanceTrs->Delete();
    _hPkStDist2W->Delete();
    _hPkStClustDist2W->Delete();
    _hPkSecDist2W->Delete();
    _hPkSecClustDist2W->Delete();
    _hPkSecVsStationDist2W->Delete();
    _hPkSecVsStationDist2WGood->Delete();
    _hPkAccArrSecStation->Delete();

    if (tmpStClustDist!=0x0) delete tmpStClustDist;
    if (tmpSecClustDist!=0x0) delete tmpSecClustDist;

    //-------------------------------------------

  } // end produce

//  void BkgTrackRejecterByGeomTplg::endJob(){
////
////    // cd() to correct root directory. See note 3.
////    TDirectory* save = gDirectory;
////    _directory->cd();
////
////    // Write canvas.  See note 4.
//////    _canvas->Write();
////
////    // cd() back to where we were.  See note 3.
////    save->cd();
////
//  }

  void BkgTrackRejecterByGeomTplg::SctrSttnMapAnalyze(TH2I *startDist, TH2I *votingArray, int minPitch, int maxPitch) {

          int nXbin = startDist->GetNbinsX();
          int nYbin = startDist->GetNbinsY();
          int lastXbin = nXbin-1;
          int *contArr = startDist->GetArray();
          int dataArr[nYbin][nXbin];
          int dataArrRev[nYbin][nXbin];
          int iY, iX, jXRange;//, jOutOfRange;
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
//          size_t couplingRwCl;
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
          std::vector<Clust>::iterator tmp_clstlst_it, end_clstlst_it;
          short tmpMaxSect_1, tmpMinSect_1, tmpMaxSect_2, tmpMinSect_2;
          bool nexClust=true, noCl2erased=true;
          end_clstlst_it=clustersList.end();
          rwclincl nonMatching_Cl1, nonMatching_Cl2, matching_Cl1_2;
          bool alreadyChecked=false, Cl1_2_NoMatching=true;

          while ( clstlst_it<end_clstlst_it/*clstlst_it!=end_clstlst_it*/ ) {
                  if (clustersList.size()<2) break;
                  nexClust=true;
                  tmp_clstlst_it=clstlst_it;
                  tmp_clstlst_it++;
                  tmpMaxSect_1=clstlst_it->_lastSectorID;
                  tmpMinSect_1=clstlst_it->_firstSectorID;
                  if ( tmpMaxSect_1>=12 ) {
                          tmpMaxSect_1-=12;
                          tmpMinSect_1-=12;
                  }
                  //for ( ; tmp_clstlst_it!=clustersList.end(); ++tmp_clstlst_it) {
                  //end_clstlst_it=clustersList.end();
                  while ( tmp_clstlst_it<end_clstlst_it/*!=clustersList.end()*/ ) {
                          if (clustersList.size()<2) break;
                          noCl2erased=true;
                          tmpMaxSect_2=tmp_clstlst_it->_lastSectorID;
                          tmpMinSect_2=tmp_clstlst_it->_firstSectorID;
                          if ( tmpMaxSect_2>=12 ) {
                                  tmpMaxSect_2-=12;
                                  tmpMinSect_2-=12;
                          }
                          cout<<"Clust 1 "<<clstlst_it->_firstSectorID<<" - "<<clstlst_it->_lastSectorID<<" renorm "<<tmpMinSect_1<<" - "<<tmpMaxSect_1<<" - "<<clstlst_it->_minStationID<<" - "<<clstlst_it->_maxStationID<<endl;
                          cout<<"Clust 2 "<<tmp_clstlst_it->_firstSectorID<<" - "<<tmp_clstlst_it->_lastSectorID<<" renorm "<<tmpMinSect_2<<" - "<<tmpMaxSect_2<<" - "<<tmp_clstlst_it->_minStationID<<" - "<<tmp_clstlst_it->_maxStationID<<endl;
//                          if ( ( tmpMaxSect_1>=tmpMaxSect_2 && tmpMinSect_1<=tmpMinSect_2 ) &&
//                               ( clstlst_it->_minStationID<=tmp_clstlst_it->_minStationID && clstlst_it->_maxStationID>=tmp_clstlst_it->_maxStationID ) ) {
//                                  clustersList.erase(tmp_clstlst_it);
//                                  cout<<"Cluster 2 removed"<<endl;
//                                  noCl2erased=false;
//                                  end_clstlst_it=clustersList.end();
//                                  //break;
//                          }
//                          //if (erased) end_clstlst_it=clustersList.end();
//                          else  {
//                                  if ( ( tmpMaxSect_2>=tmpMaxSect_1 && tmpMinSect_2<=tmpMinSect_1 ) &&
//                                                  ( tmp_clstlst_it->_minStationID<=clstlst_it->_minStationID && tmp_clstlst_it->_maxStationID>=clstlst_it->_maxStationID ) ) {
//                                          clustersList.erase(clstlst_it);
//                                          cout<<"Cluster 1 removed"<<endl;
//                                          nexClust=false;
//                                          end_clstlst_it=clustersList.end();
//                                          break;
//                                  }
//                                  else {
                                          //identify and remove cluster that partly match
                                          cout<<"--- Cl dist "<<abs(clstlst_it->_mean_Sttn-tmp_clstlst_it->_mean_Sttn)<<" err "<<3.0*sqrt(pow(clstlst_it->_sigma_Sttn,2)+pow(tmp_clstlst_it->_sigma_Sttn,2))<<endl;
                                          if ( abs(clstlst_it->_mean_Sttn-tmp_clstlst_it->_mean_Sttn)<
                                               3.0*sqrt(pow(clstlst_it->_sigma_Sttn,2)+pow(tmp_clstlst_it->_sigma_Sttn,2)) ) {
                                                  nonMatching_Cl1.clear();
                                                  nonMatching_Cl2.clear();
                                                  matching_Cl1_2.clear();
                                                  for ( rwclincl::iterator rwclincl_it=clstlst_it->_rClusts.begin(); rwclincl_it!=clstlst_it->_rClusts.end(); ++rwclincl_it) {

                                                          pair<rwclincl::iterator,rwclincl::iterator> potentMatch = tmp_clstlst_it->_rClusts.equal_range(rwclincl_it->first+12);
                                                          if (potentMatch.first==tmp_clstlst_it->_rClusts.end()) nonMatching_Cl1.insert( rwclincl::value_type( rwclincl_it->first, rwclincl_it->second ) );
                                                          else {
                                                                  Cl1_2_NoMatching=true;
                                                                  for ( rwclincl::iterator tmp_rwclincl_it=potentMatch.first; tmp_rwclincl_it!=potentMatch.second; ++tmp_rwclincl_it ) {
                                                                          if ( (*(rwclincl_it->second))==(*(tmp_rwclincl_it->second))/*tmp_rwclincl_it->second._internal_equiv( tmp_match_it->second)*/ ) {
                                                                                  matching_Cl1_2.insert( rwclincl::value_type( rwclincl_it->first, rwclincl_it->second ) );
                                                                                  Cl1_2_NoMatching=false;
                                                                                  break;
                                                                          }
                                                                          //else nonMatching_Cl2.insert( rwclincl::value_type( tmp_rwclincl_it->first, tmp_rwclincl_it->second ) );
                                                                  }
                                                                  if (Cl1_2_NoMatching) nonMatching_Cl1.insert( rwclincl::value_type( rwclincl_it->first, rwclincl_it->second ) );
                                                          }
                                                  }
                                                  for ( rwclincl::iterator tmp_rwclincl_it=tmp_clstlst_it->_rClusts.begin(); tmp_rwclincl_it!=tmp_clstlst_it->_rClusts.end(); ++tmp_rwclincl_it) {
                                                          alreadyChecked=false;
                                                          for ( rwclincl::iterator tmp_match_it=matching_Cl1_2.begin(); tmp_match_it!=matching_Cl1_2.end(); ++tmp_match_it ) {
                                                                  if ( (*(tmp_rwclincl_it->second))==(*(tmp_match_it->second)) ) {
                                                                          alreadyChecked=true;
                                                                          //Cl1_2_NoMatching=false;
                                                                          break;
                                                                  }
                                                          }
                                                          if (alreadyChecked) continue;
                                                          else nonMatching_Cl2.insert( rwclincl::value_type( tmp_rwclincl_it->first, tmp_rwclincl_it->second ) );

                                                  }

//                                                  for ( rwclincl::iterator rwclincl_it=clstlst_it->_rClusts.begin(); rwclincl_it!=clstlst_it->_rClusts.end(); ++rwclincl_it) {
//                                                          for ( rwclincl::iterator tmp_rwclincl_it=tmp_clstlst_it->_rClusts.begin(); tmp_rwclincl_it!=tmp_clstlst_it->_rClusts.end(); ++tmp_rwclincl_it) {
//                                                                  Cl1_2_NoMatching=true;
//                                                                  alreadyChecked=false;
//                                                                  cout<<"--- rw clust 1 "<<rwclincl_it->second->_firstStationID<<" - "<<rwclincl_it->second->_lastStationID<<" - "<<rwclincl_it->second->_nHit<<endl;
//                                                                  cout<<"--- rw clust 2 "<<tmp_rwclincl_it->second->_firstStationID<<" - "<<tmp_rwclincl_it->second->_lastStationID<<" - "<<tmp_rwclincl_it->second->_nHit<<endl;
//                                                                  cout<<"--- matching? "<<((*(rwclincl_it->second))==(*(tmp_rwclincl_it->second)))<<endl;
//                                                                  for ( rwclincl::iterator tmp_match_it=matching_Cl1_2.begin(); tmp_match_it!=matching_Cl1_2.end(); ++tmp_match_it ) {
//                                                                          if ( (*(tmp_rwclincl_it->second))==(*(tmp_match_it->second))/*tmp_rwclincl_it->second._internal_equiv( tmp_match_it->second)*/ ) {
//                                                                                  alreadyChecked=true;
//                                                                                  Cl1_2_NoMatching=false;
//                                                                                  break;
//                                                                          }
//                                                                  }
//                                                                  if (alreadyChecked) continue;
//                                                                  for ( rwclincl::iterator tmp_match_it=nonMatching_Cl2.begin(); tmp_match_it!=nonMatching_Cl2.end(); ++tmp_match_it ) {
//                                                                          if ( (*(tmp_rwclincl_it->second))==(*(tmp_match_it->second))/*tmp_rwclincl_it->second._internal_equiv( tmp_match_it->second)*/ ) {
//                                                                                  alreadyChecked=true;
//                                                                                  Cl1_2_NoMatching=false;
//                                                                                  break;
//                                                                          }
//                                                                  }
//                                                                  if (alreadyChecked) continue;
//
//                                                                  if ( (*(rwclincl_it->second))==(*(tmp_rwclincl_it->second))/*rwclincl_it->second._internal_equiv( tmp_rwclincl_it->second)*/ ) {
//                                                                          matching_Cl1_2.insert( rwclincl::value_type( rwclincl_it->first, rwclincl_it->second ) );
//                                                                          Cl1_2_NoMatching=false;
//                                                                          break;
//                                                                  }
//                                                                  nonMatching_Cl2.insert( rwclincl::value_type( tmp_rwclincl_it->first, tmp_rwclincl_it->second ) );
//                                                          }
//                                                          if (Cl1_2_NoMatching) nonMatching_Cl1.insert( rwclincl::value_type( rwclincl_it->first, rwclincl_it->second ) );
//                                                  }

                                                  /*cout<<"--------------------------------------------------------"<<endl;
                                                  for ( rwclincl::iterator tmp_match_it=matching_Cl1_2.begin(); tmp_match_it!=matching_Cl1_2.end(); ++tmp_match_it ) {
                                                          cout<<"--- rw clust in Cl1 that matches Cl2: "<<"rw "<<tmp_match_it->first<<" @ "<<tmp_match_it->second->_firstStationID<<" - "<<tmp_match_it->second->_lastStationID<<" n Hit "<<tmp_match_it->second->_nHit<<endl;
                                                  }
                                                  cout<<"--------------------------------------------------------"<<endl;
                                                  for ( rwclincl::iterator tmp_match_it=nonMatching_Cl1.begin(); tmp_match_it!=nonMatching_Cl1.end(); ++tmp_match_it ) {
                                                          cout<<"--- rw clust in Cl1 that doens't match Cl2: "<<"rw "<<tmp_match_it->first<<" @ "<<tmp_match_it->second->_firstStationID<<" - "<<tmp_match_it->second->_lastStationID<<" n Hit "<<tmp_match_it->second->_nHit<<endl;
                                                  }
                                                  cout<<"--------------------------------------------------------"<<endl;
                                                  for ( rwclincl::iterator tmp_match_it=nonMatching_Cl2.begin(); tmp_match_it!=nonMatching_Cl2.end(); ++tmp_match_it ) {
                                                          cout<<"--- rw clust in Cl2 that doens't match Cl1: "<<"rw "<<tmp_match_it->first<<" @ "<<tmp_match_it->second->_firstStationID<<" - "<<tmp_match_it->second->_lastStationID<<" n Hit "<<tmp_match_it->second->_nHit<<endl;
                                                  }
                                                  cout<<"--------------------------------------------------------"<<endl;
                                                 */
                                                  if ( !matching_Cl1_2.empty() ) {
                                                          if ( !nonMatching_Cl1.empty() && nonMatching_Cl2.empty() ) {
                                                                  clustersList.erase(tmp_clstlst_it);
                                                                  cout<<"Cluster 2 removed"<<endl;
                                                                  noCl2erased=false;
                                                                  end_clstlst_it=clustersList.end();
                                                                  //break;
                                                          }
                                                          else if ( nonMatching_Cl1.empty() && !nonMatching_Cl2.empty() ) {
                                                                  clustersList.erase(clstlst_it);
                                                                  cout<<"Cluster 1 removed"<<endl;
                                                                  nexClust=false;
                                                                  end_clstlst_it=clustersList.end();
                                                                  break;
                                                          }
                                                          else {
                                                                  if ( tmpMinSect_2<0 ) {
                                                                          cout<<"Cluster 1 removed and Cluster 2 modified"<<endl;
                                                                          //cout<<"no Matching rw cluster in Cl1 "<<nonMatching_Cl1.size()<<endl;
                                                                          for ( rwclincl::iterator tmp_match_it=nonMatching_Cl1.begin(); tmp_match_it!=nonMatching_Cl1.end(); ++tmp_match_it ) {
                                                                                  //cout<<"--- rw clust in Cl1 that doens't match Cl2 "<<tmp_match_it->second->_firstStationID<<" - "<<tmp_match_it->second->_lastStationID<<" - "<<tmp_match_it->second->_nHit<<endl;
                                                                                  /*if ( (tmp_match_it->first + 12) <16 )*/ tmp_clstlst_it->addRwClust( tmp_match_it->first+12, tmp_match_it->second );
                                                                                  //cout<<"--- row cl added to Cl 2"<<endl;
                                                                          }
                                                                          //cout<<"--- All no matching row cls of Cl 1 added to Cl 2"<<endl;
                                                                          clustersList.erase(clstlst_it);
                                                                          nexClust=false;
                                                                          end_clstlst_it=clustersList.end();
                                                                          break;
                                                                  }
                                                                  else {
                                                                          cout<<"Cluster 2 removed and Cluster 1 modified"<<endl;
                                                                          //cout<<"no Matching rw cluster in Cl2 "<<nonMatching_Cl1.size()<<endl;
                                                                          for ( rwclincl::iterator tmp_match_it=nonMatching_Cl2.begin(); tmp_match_it!=nonMatching_Cl2.end(); ++tmp_match_it ) {
                                                                                  if ( (tmp_match_it->first - 12) >0 ) clstlst_it->addRwClust( tmp_match_it->first-12, tmp_match_it->second );
                                                                                  //cout<<"--- row cl added to Cl 1"<<endl;
                                                                          }
                                                                          //cout<<"--- All no matching row cls of Cl 2 added to Cl 1"<<endl;
                                                                          clustersList.erase(tmp_clstlst_it);
                                                                          noCl2erased=false;
                                                                          end_clstlst_it=clustersList.end();
                                                                          //break;
                                                                  }
                                                          }
                                                  }
                                                  //cout<<"New Clusters group stored n: "<<clustersList.size()<<endl;
                                                  //if (clustersList.size()<2) break;
                                          }
//                                  }//end of fixing the partly match clusters issue
//
//                                  //++tmp_clstlst_it;
//                          }
                          if (noCl2erased) ++tmp_clstlst_it;
                  }
                  if ( nexClust ) {
                          clstlst_it++;
                          //end_clstlst_it=clustersList.end();
                  }
                  //clstlst_it=clustersList.begin();
          }// periodic (in sector) clusters removed

          cout<<"------------ double cluster removed ------------"<<endl;
          cout<<"Cluster List size "<<clustersList.size()<<endl;
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

  void BkgTrackRejecterByGeomTplg::findCluster(unsigned short tmpRowId, rwClustPtr &startingRowClust){

          rwClustPtr tmpRwClust;
          rwclinrwrel::iterator tmp_rcr_it;

          tmpRowId++;
          tmp_rcr_it=rwClst_forRw_rel.find(tmpRowId);
          bool notFound;
          if (tmp_rcr_it!=rwClst_forRw_rel.end()) {
                  //cout<<"Up) Lokking for row cluster to add in row "<<tmpRowId<<endl;
                  if (!tmp_rcr_it->second.empty()) {
                          //for ( rwclstvec::iterator tmp_cr_it=tmp_rcr_it->second.begin(); tmp_cr_it!=tmp_rcr_it->second.end(); ++tmp_cr_it) {
                          rwclstvec::iterator tmp_cr_it=tmp_rcr_it->second.begin();
                          rwclstvec::iterator end_cr_it=tmp_rcr_it->second.end();
                          while ( tmp_cr_it!=end_cr_it ) {
                                  notFound=true;
                                  tmpRwClust=*tmp_cr_it;
                                  if ( (tmpRwClust->_firstStationID-startingRowClust->_lastStationID)<=1 &&
                                                  (startingRowClust->_firstStationID-tmpRwClust->_lastStationID)<=1 ) {
                                          clustersList.back().addRwClust(tmpRowId, tmpRwClust);
                                          tmp_rcr_it->second.erase(tmp_cr_it);
                                          findCluster(tmpRowId,tmpRwClust);
                                          //break;
                                          notFound=false;
                                  }
                                  if (notFound) ++tmp_cr_it;
                                  else if (tmp_rcr_it->second.empty()) break;
                                  end_cr_it=tmp_rcr_it->second.end();
                                  if (distance(tmp_cr_it,end_cr_it)<0) break;
                          }
                  }
          }
          tmpRowId--; tmpRowId--;
          tmp_rcr_it=rwClst_forRw_rel.find(tmpRowId);
          if (tmp_rcr_it!=rwClst_forRw_rel.end()) {
                  //cout<<"Down) Lokking for row cluster to add in row "<<tmpRowId<<endl;
                  if (!tmp_rcr_it->second.empty()) {
                          //for ( rwclstvec::iterator tmp_cr_it=tmp_rcr_it->second.begin(); tmp_cr_it!=tmp_rcr_it->second.end(); ++tmp_cr_it) {
                          rwclstvec::iterator tmp_cr_it=tmp_rcr_it->second.begin();
                          rwclstvec::iterator end_cr_it=tmp_rcr_it->second.end();
                          while ( tmp_cr_it!=end_cr_it ) {
                                  notFound=true;
                                  tmpRwClust=*tmp_cr_it;
                                  if ( (tmpRwClust->_firstStationID-startingRowClust->_lastStationID)<=1 &&
                                                  (startingRowClust->_firstStationID-tmpRwClust->_lastStationID)<=1 ) {
                                          clustersList.back().addRwClust(tmpRowId, tmpRwClust);
                                          tmp_rcr_it->second.erase(tmp_cr_it);
                                          findCluster(tmpRowId,tmpRwClust);
                                          //break;
                                          notFound=false;
                                  }
                                  if (notFound) ++tmp_cr_it;
                                  else if (tmp_rcr_it->second.empty()) break;
                                  end_cr_it=tmp_rcr_it->second.end();
                                  if (distance(tmp_cr_it,end_cr_it)<0) break;
                          }
                 }
          }
  }


}  // end namespace mu2e

using mu2e::BkgTrackRejecterByGeomTplg;
DEFINE_ART_MODULE(BkgTrackRejecterByGeomTplg);
