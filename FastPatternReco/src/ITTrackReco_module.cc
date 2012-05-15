//
// Fast Patter recognition for the ITracker
//
// $Id: ITTrackReco_module.cc,v 1.1 2012/05/15 07:52:41 tassiell Exp $
// $Author: tassiell $
// $Date: 2012/05/15 07:52:41 $
//
// Original author G. Tassielli
//

// C++ includes.
#include <iostream>
//#include <string>
#include <new>
#include <memory>
#include <set>
#include <map>
#include <utility>
#include <cstring>
#include <algorithm>
#include <limits>
#include <ctime>
#include <cmath>

#include <boost/shared_ptr.hpp>

// Framework includes.
//#include "art/Framework/Core/EDProducer.h"
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

// CLHEP includes.
#include "CLHEP/Vector/TwoVector.h"

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "ITrackerGeom/inc/ITracker.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "ITrackerGeom/inc/Cell.hh"
//#include "DataProducts/inc/DPIndexVectorCollection.hh"
//#include "MCDataProducts/inc/GenParticleCollection.hh"
//#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
//#include "MCDataProducts/inc/SimParticleCollection.hh"
//#include "MCDataProducts/inc/StepPointMCCollection.hh"
//#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
//#include "FastPatternReco/inc/GenTrackData.hh"
#include "RecoDataProducts/inc/TrackerHitTimeCluster.hh"
#include "RecoDataProducts/inc/TrackerHitTimeClusterCollection.hh"
#include "RecoDataProducts/inc/TrackerHitByID.hh"
#include "FastPatternReco/inc/FastPatRecoUtilsAndDataDef.hh"

// Root includes.
#include "TApplication.h"
#include "TCanvas.h"
//#include "TDirectory.h"
#include "TMath.h"
#include "TH1I.h"
#include "TH1F.h"
#include "TH2I.h"
#include "TH2F.h"
#include "TGraphErrors.h"
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

#define invSqrt3 0.577350269
#define invSqrt12 0.288675135

using namespace std;

namespace mu2e {

  struct confMapDraw{

          confMapDraw():
                  _cmap(0x0),
                  _cmapCHT(0x0)
          {}

          confMapDraw( TGraphErrors *cmap, TH2F *cmapCHT):
                  _cmap(cmap),
                  _cmapCHT(cmapCHT)
          {}

          ~confMapDraw() {
                  cout<<" pointers "<<_cmap<<" "<<_cmapCHT<<endl;
                  if (_cmap!=0x0) delete _cmap;
                  if (_cmapCHT!=0x0) delete _cmapCHT;
          }

          //TH2F *_cmap;
          TGraphErrors *_cmap;
          TH2F *_cmapCHT;

  };


  class ITTrackReco : public art::EDAnalyzer {

  public:

    explicit ITTrackReco(fhicl::ParameterSet const& pset);
    virtual ~ITTrackReco() {
//            _peakFinder->Delete();
//            _fg->Delete();
//            if (_peakFinder!=0x0)        delete _peakFinder;
//            if (_fg!=0x0)                delete _fg;
            if (_fakeCanvas!=0x0)        delete _fakeCanvas;
//            if (_peaksCanvHistos!=0x0)   delete _peaksCanvHistos;
//            if (_hPkStDistTrs!=0x0)      delete _hPkStDistTrs;
////            if (_hPkStDistanceTrs!=0x0)  delete _hPkStDistanceTrs;
//            if (_hPkStDist2W!=0x0)       delete _hPkStDist2W;
//            if (_hPkStClustDist2W!=0x0)  delete _hPkStClustDist2W;
//            if (_hPkSecDist2W!=0x0)      delete _hPkSecDist2W;
//            if (_hPkSecClustDist2W!=0x0) delete _hPkSecClustDist2W;
//            if (_hPkSecVsStationDist2W!=0x0)      delete _hPkSecVsStationDist2W;
//            if (_hPkSecVsStationDist2WGood!=0x0)      delete _hPkSecVsStationDist2WGood;
//            if (_hPkAccArrSecStation!=0x0)      delete _hPkAccArrSecStation;
    }

    virtual void beginJob();
    void endJob();

    // This is called for each event.
    //void produce(art::Event & e);
    void analyze(art::Event const& e);

  private:

    // Start: run time parameters

    // The module label of this module.
    std::string _moduleLabel;

    // Label of the module that made the hits.
    std::string _makerModuleLabel;

    // Label of the module that made the hits.
    std::string _timeRejecterModuleLabel;

    // Label of the process for the remaping of the Tracker hit by Cell/Straw ID.
    std::string _mapTrackerHitByID;

    // Threshold on the multiplicity of Closest Hit Clusters to identify a potential bending point
    int _minClusMultBend;

    // Max acceptable value for the Chi2/NDOF for a circle
    float _circChi2cut;

    // Use n times sigma for matching the clusters in the Sector plane
    //float _nSigmaForClMatch;

    // Use n times sigma for matching the clusters slope with 0
    //float _nSigmaForClSlope;

    // Diagnostics level.
    int _diagLevel;

    bool _doDisplay;

//    // Number of events to accumulate between prompts.
//    int _nAccumulate;

    // End: run time parameters

    // Pointers to histograms, ntuples, TGraphs.

    //TObjArray*    _peaksCanvHistos;
    //TObjArray*    _hPkStDistTrs;
    //TObjArray*    _hPkStDist2W;
    //TObjArray*    _hPkStClustDist2W;
    //TObjArray*    _hPkSecDist2W;
    //TObjArray*    _hPkSecClustDist2W;
    //TObjArray*    _hPkSecVsStationDist2W;
    //TObjArray*    _hPkSecVsStationDist2WGood;
    //TObjArray*    _hPkAccArrSecStation;

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

    //void SctrSttnMapAnalyze(TH2I *startDist, TH2I *votingArray, int minPitch=3, int maxPitch=11);
    //rwclinrwrel rwClst_forRw_rel;
    //std::vector<Clust> clustersList;
    bool  findCluster(int &did, int &lid, int &sid, size_t &iglbHit, std::set<size_t> &hitLoked, closClstcol &clusCont, TrackerHitByID const* hitByID, TrackerHitTimeCluster const&  tclust, CellGeometryHandle *itwp);
    bool  findClosClust( unsigned int tmpX, unsigned int tmpY, std::map<unsigned int, bool> &notLooked, CHTVotArr &cmapCHTVotArr, ClosClustCHTVot &votClus );
    bool  findCircles(closClinRadLay &radCls, CellGeometryHandle *itwp, std::vector<confMapDraw*> &rCMap, circlesCont & potcircles, std::string type="", unsigned int thrCHTaddVot=5);
    void  mergeCircles( circlesCont &potCircles );
    void  printPotCircles( circlesCont &potCircles );

    // Some ugly but necessary ROOT related bookkeeping:

    // The job needs exactly one instance of TApplication.  See note 1.
    auto_ptr<TApplication> _application;

    // Save directory from beginJob so that we can go there in endJob. See note 3.
//    TDirectory* _directory;

  };

  ITTrackReco::ITTrackReco(fhicl::ParameterSet const& pset) :

    // Run time parameters
    _moduleLabel(pset.get<string>("module_label")),/*@module_label*/
    _makerModuleLabel(pset.get<string>("makerModuleLabel")),
    _timeRejecterModuleLabel(pset.get<string>("tRejecterModuleLabel")),
    _mapTrackerHitByID(pset.get<string>("mapTrackerHitByID")),
    _minClusMultBend(pset.get<int>("minClusMultBend")),
    _circChi2cut(pset.get<float>("circChi2cut")),
    //_nSigmaForClMatch(pset.get<float>("nSigmaForClMatch")),
    //_nSigmaForClSlope(pset.get<float>("nSigmaForClSlope")),
    /*_nAccumulate(pset.get<int>("nAccumulate",20)),*/
    _diagLevel(pset.get<int>("diagLevel",0)),
    _doDisplay(pset.get<bool>("doDisplay",false)),

    // ROOT objects that are the main focus of this example.
    //_peaksCanvHistos(0),
    //_hPkStDistTrs(0),
    //_hPkStDist2W(0),
    //_hPkStClustDist2W(0),
    //_hPkSecDist2W(0),
    //_hPkSecClustDist2W(0),
    //_hPkSecVsStationDist2W(0),
    //_hPkSecVsStationDist2WGood(0),
    //_hPkAccArrSecStation(0),
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
          //if ( _nSigmaForClMatch<0.0 ) _nSigmaForClMatch*=-1.00000;
          //if ( _nSigmaForClSlope<0.0 ) _nSigmaForClSlope*=-1.00000;
          // Tell the framework what we make.
          //produces<SctrSttnClusterGroupCollection>();

  }

  void ITTrackReco::beginJob(){

    cout<<"Bkg rejection by Geom job!"<<endl;

    // Get access to the TFile service and save current directory for later use.
    art::ServiceHandle<art::TFileService> tfs;
    _hClockCicles       = tfs->make<TH1I>( "hClockCicles",   "N clock cicles needed to analyze one Event by ITTrackReco", 2000, 5000, 15000  );
    _hExecTime          = tfs->make<TH1F>( "hExecTime",   "Execution time to analyze one Event by BkgTrackRejecterByTime", 1000, 0.0, 1.0  );


    ntimeBin    = (int) maxTimeHist/timeBinDim;

    // Create a histogram.

    //_hPkStDistTrs      = new TObjArray();
    //_hPkStDist2W       = new TObjArray();
    //_hPkStClustDist2W  = new TObjArray();
    //_hPkSecDist2W      = new TObjArray();
    //_hPkSecClustDist2W = new TObjArray();
    //_hPkSecVsStationDist2W      = new TObjArray();
    //_hPkSecVsStationDist2WGood      = new TObjArray();
    //_hPkAccArrSecStation      = new TObjArray();

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

            _fakeCanvas = new TCanvas("canvas_Fake","double click for next event",500,100);

    }

    //_peakFinder = new TSpectrum(20);

    // See note 3.
//    _directory = gDirectory;

  }

//  void ITTrackReco::produce(art::Event & event ) {
  void ITTrackReco::analyze(art::Event const& event ) {

/*

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

    //auto_ptr<SctrSttnClusterGroupCollection> sscc(new SctrSttnClusterGroupCollection);


    const Tracker& tracker = getTrackerOrThrow();
    const ITracker &itr = static_cast<const ITracker&>( tracker );
    CellGeometryHandle *itwp = itr.getCellGeometryHandle();

    double rIn  = itr.r0()+itr.getWalls()->find(Wall::inner)->second->getTotalThickness();
    double rOut = itr.rOut()-itr.getWalls()->find(Wall::outer)->second->getTotalThickness();
    rIn/=CLHEP::cm;
    rOut/=CLHEP::cm;

    stMaprel ssmap_Bin_Straw_rel;           //relation between straw and the AbsSector (Rotation) vs Station map
    stMaprel::iterator ssmpbnstrw_rel_it;   //iterator for relation between straw and the AbsSector (Rotation) vs Station map
    int tmpiGeomBin, tmpiGeomBinRow;

    art::Handle<StrawHitCollection> pdataHandle;
    event.getByLabel(_makerModuleLabel,pdataHandle);
    StrawHitCollection const* hits = pdataHandle.product();

    art::Handle<TrackerHitTimeClusterCollection> tclustHandle;
    event.getByLabel(_timeRejecterModuleLabel,tclustHandle);
    TrackerHitTimeClusterCollection const* tclusts = tclustHandle.product();

    art::Handle<TrackerHitByID> hitByIDHandle;
    event.getByLabel(_mapTrackerHitByID,hitByIDHandle);
    TrackerHitByID const* hitByID = hitByIDHandle.product();
    TrackerHitByID::const_iterator hitByID_it;
    std::pair<TrackerHitByID::const_iterator, TrackerHitByID::const_iterator> rangeHitByID_it;
    //std::set<size_t> hitLoked;


    clock_t startClock = clock();

    size_t nStrawPerEvent = hits->size();
    size_t nTimeClusPerEvent = tclusts->size();
    cout<<"----------------------------------------------------------------------------------"<<endl;
    cout<<"event "<<event.id().event()<<" tot N hit "<<nStrawPerEvent<<" N peak found "<<nTimeClusPerEvent<<endl;
    cout<<"----------------------------------------------------------------------------------"<<endl;
    int i1peak;

    int sid, lid, did, radId;
    unsigned int absSect/*, devStId*/;
    unsigned int ihit;

    //double refX, refY, /*tmpX, tmpY,*/ cmapU, cmapV, invTmpR2;
    std::vector<confMapDraw*> rCMap_Pl, rCMap_Min;// rCMap_Pl_Ups, rCMap_Min_Ups;
    size_t iglbHit=0;
    bool *isStereoMin = new (nothrow) bool [nStrawPerEvent];
    for (int jhit=0; jhit<nStrawPerEvent; ++jhit) {
            isStereoMin[jhit] = false;
    }

    //TH2I *tmpStClustDist    = new TH2I("tmpStClustDist","tmp Smoothing of Device vs Straw multiplicity Distribution",36,0,36,1000,-200,800);
    //TH2I *tmpSecClustDist   = new TH2I("tmpSecClustDist","tmp Smoothing of Device vs Sector multiplicity Distribution",36,0,36,20,-4,16);
    //TH2I *ihPkStDistTrs, *ihPkStDist2W, *ihPkStClustDist2W, *ihPkSecDist2W, *ihPkSecClustDist2W, *ihPkSecVsStationDist2WG, *ihPkSecVsStationDist2W, *ihPkAccArrSecStation;
    //TCanvas * iCanvHist=0x0;

    //int nCellPerLayer = 0;
    int tmpSid=0, lost=0, hitClSizeInLayer=1;
    int nClus=0, nTotHitInCls=0, nTotHitInCls_Pl=0, nTotHitInCls_MIn=0, nTotHit_Pl=0, nTotHit_MIn=0;
    bool lostHit=true;
    //closClstcol
    closClinRadLay radCls_Pl, radCls_Min;// radCls_Pl_Ups, radCls_Min_Ups;
    closClinRadLay::iterator radCls_it, sndRadCls_it;

    for (size_t ipeak=0; ipeak<nTimeClusPerEvent; ipeak++) {
    	TrackerHitTimeCluster const&  tclust(tclusts->at(ipeak));
    	cout<<"\t Hit in peak "<<tclust._selectedTrackerHits.size()<<endl;

    	i1peak=ipeak;
        ++i1peak;

        std::set<size_t> hitLoked;
        ihit=0;
        nClus=nTotHitInCls=0;
        nTotHit_Pl=nTotHit_MIn=0;
        for (std::vector<StrawHitPtr>::const_iterator iTCHit=tclust._selectedTrackerHits.begin(); iTCHit!=tclust._selectedTrackerHits.end(); ++iTCHit){
                // Access data
//                cout<<"\t\t iglbHit in peak at "<<*iTCHit<<endl;
                //StrawHit        const&      hit(hits->at(*iTCHit));
                iglbHit = iTCHit->key();
                StrawHit        const&      hit=*(*iTCHit);
                StrawIndex si = hit.strawIndex();
                const Straw & str = tracker.getStraw(si);
                const Cell & cell = static_cast<const Cell&>( str );


                // cout << "Getting informations about cells" << endl;

                sid = cell.Id().getCell();
                lid = cell.Id().getLayer();
                did = cell.Id().getLayerId().getSuperLayer();
                itwp->SelectCell(did,lid,sid);
                radId = itwp->GetCellAbsRadID();

                //const CLHEP::Hep3Vector stMidPoint3 = itwp->GetCellCenter();//str.getMidPoint();

                //tmpX = stMidPoint3.getX();
                //tmpY = stMidPoint3.getY();
                //cmapV = 1.0/(tmpX*tmpX+tmpY*tmpY);
                //cmapU = tmpX*cmapV;
                //cmapV*= tmpY;
                //-------------------------------------------------
                //Hit cluster per layer analysis
                if (isStereoMin[iglbHit] || cell.getWire()->getEpsilon()<0.0) {
                        isStereoMin[iglbHit]=true;
                        ++nTotHit_MIn;
                        //idHitStereoMin.push_back(iglbHit);
                        //_hCMapHitStereoMi->Fill(cmapU,cmapV);

                        radCls_it=radCls_Min.find(radId);
                        if (radCls_it==radCls_Min.end()) {
                                radCls_it=( radCls_Min.insert( closClinRadLay::value_type( radId, closClstcol() ) ) ).first;
                        }
                        findCluster(did, lid, sid, iglbHit, hitLoked, radCls_it->second, hitByID, tclust, itwp);
                }
                else {
                        //isStereoMin[iglbHit]=false;
                        ++nTotHit_Pl;

                        radCls_it=radCls_Pl.find(radId);
                        if (radCls_it==radCls_Pl.end()) {
                                radCls_it=( radCls_Pl.insert( closClinRadLay::value_type( radId, closClstcol() ) ) ).first;
                        }
                        findCluster(did, lid, sid, iglbHit, hitLoked, radCls_it->second, hitByID, tclust, itwp);
                }
                //-------------------------------------------------

                ihit++;
        }

        if(_diagLevel > 1){
                nClus=0;
                nTotHitInCls=0;
                nTotHitInCls_MIn=nTotHitInCls_Pl=0;
                cout<<"Clusters for Min streo:"<<endl;

                for ( closClinRadLay::iterator radCls_it = radCls_Min.begin(); radCls_it != radCls_Min.end(); ++radCls_it ) {
                        cout<<"***** cluster for Radial layer "<<radCls_it->first<<" :"<<endl;
                        for ( closClstcol::iterator closCls_it = radCls_it->second.begin(); closCls_it != radCls_it->second.end(); ++closCls_it ) {
                                cout<<"-- cluster n "<<nClus<<endl;
                                ++nClus;
                                nTotHitInCls_MIn+=(*closCls_it)->_nHit;
                                nTotHitInCls+=(*closCls_it)->_nHit;
                                cout<<"\t cluster size "<<(*closCls_it)->_nHit<<" minCell "<<(*closCls_it)->_minCellID<<" maxCell "<<(*closCls_it)->_maxCellID<<" centralCell "<<(*closCls_it)->_centerCellID<<endl;
                                cout<<"\t hits in cluster: ";
                                for (hitsInClsID::const_iterator clsHit_it = (*closCls_it)->_hitIdx.begin(); clsHit_it != (*closCls_it)->_hitIdx.end(); ++clsHit_it){
                                        cout<<"("<<clsHit_it->first<<", "<<clsHit_it->second<<") ";
                                }
                                cout<<endl;
                        }
                }
                cout<<"Clusters for Plus streo:"<<endl;

                for ( closClinRadLay::iterator radCls_it = radCls_Pl.begin(); radCls_it != radCls_Pl.end(); ++radCls_it ) {
                        cout<<"***** cluster for Radial layer "<<radCls_it->first<<" :"<<endl;
                        for ( closClstcol::iterator closCls_it = radCls_it->second.begin(); closCls_it != radCls_it->second.end(); ++closCls_it ) {
                                cout<<"-- cluster n "<<nClus<<endl;
                                ++nClus;
                                nTotHitInCls_Pl+=(*closCls_it)->_nHit;
                                nTotHitInCls+=(*closCls_it)->_nHit;
                                cout<<"\t cluster size "<<(*closCls_it)->_nHit<<" minCell "<<(*closCls_it)->_minCellID<<" maxCell "<<(*closCls_it)->_maxCellID<<" centralCell "<<(*closCls_it)->_centerCellID<<endl;
                                cout<<"\t hits in cluster: ";
                                for (hitsInClsID::const_iterator clsHit_it = (*closCls_it)->_hitIdx.begin(); clsHit_it != (*closCls_it)->_hitIdx.end(); ++clsHit_it){
                                        cout<<"("<<clsHit_it->first<<", "<<clsHit_it->second<<") ";
                                }
                                cout<<endl;
                        }
                }

                cout<<"In Time peak :"<<endl;
                cout<<" n. Hit for Stereo Min "<<nTotHit_MIn<<" n. Hit for Stereo Pl "<<nTotHit_Pl<<endl;
                cout<<"NClus "<<nClus<<" n. Tot hit grouped in clusters "<<nTotHitInCls<<" n. Hit in Clusters for Stereo Min "<<nTotHitInCls_MIn<<" n. Hit in Clusters for Stereo Pl "<<nTotHitInCls_Pl<<endl;
        }

        cout<<"------- start research ------"<<endl;
        isHitIDUsed associedHits;

        cout<<"   ----  stereo - ----   "<<endl;
        circlesCont potCirclesMin;
        bool foundCircles_Min = findCircles(radCls_Min,itwp,rCMap_Min,potCirclesMin,"Min");
        if ( foundCircles_Min ) {
                printPotCircles(potCirclesMin);
                mergeCircles(potCirclesMin);
                printPotCircles(potCirclesMin);
        }
        cout<<"   ----  stereo + ----   "<<endl;
        circlesCont potCirclesPl;
        bool foundCircles_Pl = findCircles(radCls_Pl,itwp,rCMap_Pl,potCirclesPl,"Pl");
        if ( foundCircles_Pl ) {
                printPotCircles(potCirclesPl);
                mergeCircles(potCirclesPl);
                printPotCircles(potCirclesPl);

        }
        //findCircles(radCls_Min_Up,itwp,rCMap_Min_Up,"Min_Up");
        //findCircles(radCls_Pl_Up,itwp,rCMap_Pl_Up,"Pl_Up");

        if ( foundCircles_Min && foundCircles_Pl) {

        }

        if (_doDisplay) {

                TCanvas *tmpRCP_Min = new TCanvas();
                int nlrgCl_Min = rCMap_Min.size();
                int iPad = 0, maxPad = 1;
                if (nlrgCl_Min>=1) {
                        tmpRCP_Min->Divide(2,nlrgCl_Min);
                        maxPad = 2*nlrgCl_Min;
                        iPad=1;
                }
                int iColor = 0;
                for (std::vector<confMapDraw*>::iterator hCmapMin_it = rCMap_Min.begin(); hCmapMin_it != rCMap_Min.end(); ++hCmapMin_it ) {
                        if (iPad>maxPad) break;
                        tmpRCP_Min->cd(iPad);
                        (*hCmapMin_it)->_cmap->SetMarkerStyle(8);
                        (*hCmapMin_it)->_cmap->SetMarkerSize(0.8);
                        (*hCmapMin_it)->_cmap->SetMarkerColor(++iColor);
                        (*hCmapMin_it)->_cmap->Draw("AP");
                        ++iPad;
                        tmpRCP_Min->cd(iPad);
                        (*hCmapMin_it)->_cmapCHT->Draw("colz");
                        ++iPad;
                }

                TCanvas *tmpRCP_Pl = new TCanvas();
                int nlrgCl_Pl = rCMap_Pl.size();
                iPad = 0;
                maxPad = 1;
                if (nlrgCl_Pl>=1) {
                        tmpRCP_Pl->Divide(2,nlrgCl_Pl);
                        maxPad = 2*nlrgCl_Pl;
                        iPad=1;
                }
                iColor = 0;
                for (std::vector<confMapDraw*>::iterator hCmapPl_it = rCMap_Pl.begin(); hCmapPl_it != rCMap_Pl.end(); ++hCmapPl_it ) {
                        if (iPad>maxPad) break;
                        tmpRCP_Pl->cd(iPad);
                        (*hCmapPl_it)->_cmap->SetMarkerStyle(8);
                        (*hCmapPl_it)->_cmap->SetMarkerSize(0.8);
                        (*hCmapPl_it)->_cmap->SetMarkerColor(++iColor);
                        (*hCmapPl_it)->_cmap->Draw("AP");
                        ++iPad;
                        tmpRCP_Pl->cd(iPad);
                        (*hCmapPl_it)->_cmapCHT->Draw("colz");
                        ++iPad;
                }

                tmpRCP_Min->Modified();
                tmpRCP_Min->Update();

                tmpRCP_Pl->Modified();
                tmpRCP_Pl->Update();



                TEllipse innerWall (0.0,0.0,rIn,rIn);
                innerWall.SetFillStyle(0);
                innerWall.SetLineWidth(1.5);

                TEllipse outerWall (0.0,0.0,rOut,rOut);
                outerWall.SetFillStyle(0);
                outerWall.SetLineWidth(1.5);

                TCanvas *tmpCircleFound = new TCanvas("PotCirclesDraw","",860,430);
                tmpCircleFound->Divide(2,1);
                TH2F hHitStereoPl( "hHitStereoPl",  "Stereo+ Hits per Event at z=0", 1500, -75.0, 75.0, 1500, -75.0, 75.0 );
                TH2F hHitStereoMin( "hHitStereoMin",  "Stereo- Hits per Event at z=0", 1500, -75.0, 75.0, 1500, -75.0, 75.0 );
                hHitStereoPl.SetXTitle("cm");
                hHitStereoPl.SetYTitle("cm");
                hHitStereoMin.SetXTitle("cm");
                hHitStereoMin.SetYTitle("cm");
                tmpCircleFound->cd(1);
                hHitStereoMin.SetStats(kFALSE);
                hHitStereoMin.Draw();
                innerWall.Draw("same");
                outerWall.Draw("same");
                int iCol = 0;
                std::vector<TEllipse *> drawPotCirclesMin;
                std::vector<TEllipse *> drawPotKrmCirclesMin;
                std::vector<TClonesArray *> drawPotCircPntsMin;
                for (circlesCont::iterator potCirclesMin_it = potCirclesMin.begin(); potCirclesMin_it != potCirclesMin.end(); ++potCirclesMin_it ) {
                        drawPotCirclesMin.push_back( new TEllipse( potCirclesMin_it->_center.x()/CLHEP::cm, potCirclesMin_it->_center.y()/CLHEP::cm, potCirclesMin_it->_radius/CLHEP::cm, potCirclesMin_it->_radius/CLHEP::cm) );
                        drawPotCirclesMin.back()->SetFillStyle(0);
                        drawPotCirclesMin.back()->SetLineWidth(2);
                        drawPotCirclesMin.back()->SetLineColor(++iCol);
                        drawPotCirclesMin.back()->SetLineStyle(4);
                        //drawPotCirclesMin.back()->Draw("same");
                        CLHEP::Hep2Vector radDir;
                        radDir.setX( cos(potCirclesMin_it->_krmCircFit.phi) );
                        radDir.setY( sin(potCirclesMin_it->_krmCircFit.phi) );
                        radDir.rotate( ( (potCirclesMin_it->_krmCircFit.rho>=0.0) ? -90.0 : 90.0 )*CLHEP::degree );
                        radDir=radDir.unit();
                        HepGeom::Point3D<double> CirCenter (0.0,0.0,0.0);
                        float krmRad = 1.0/fabs(potCirclesMin_it->_krmCircFit.rho);
                        CirCenter = CirCenter +  (krmRad-potCirclesMin_it->_krmCircFit.dca)*radDir;
                        drawPotKrmCirclesMin.push_back( new TEllipse( CirCenter.x()/CLHEP::cm, CirCenter.y()/CLHEP::cm, krmRad/CLHEP::cm, krmRad/CLHEP::cm) );
                        //radDir *= (krmRad-potCirclesMin_it->_krmCircFit.dca);
                        //drawPotKrmCirclesMin.push_back( new TEllipse( radDir.x()/CLHEP::cm, radDir.y()/CLHEP::cm, krmRad/CLHEP::cm, krmRad/CLHEP::cm) );
                        drawPotKrmCirclesMin.back()->SetFillStyle(0);
                        drawPotKrmCirclesMin.back()->SetLineWidth(2);
                        drawPotKrmCirclesMin.back()->SetLineColor(iCol);
                        drawPotKrmCirclesMin.back()->SetLineStyle(4);
                        drawPotKrmCirclesMin.back()->Draw("same");

                        drawPotCircPntsMin.push_back( new TClonesArray("TEllipse") );
                        drawPotCircPntsMin.back()->ExpandCreateFast(potCirclesMin_it->_listHitptrs.size());
                        int iHit = 0;
                        for ( std::vector< ptrHitInClosCust >::iterator potCirclesPntsMin_it = potCirclesMin_it->_listHitptrs.begin(); potCirclesPntsMin_it != potCirclesMin_it->_listHitptrs.end(); ++potCirclesPntsMin_it ) {
                                itwp->SelectCell( potCirclesPntsMin_it->getRadLayID(), potCirclesPntsMin_it->getinLayerCellID() );
                                ((TEllipse *) drawPotCircPntsMin.back()->At(iHit))->SetX1(itwp->GetCellCenter().x()/CLHEP::cm);
                                ((TEllipse *) drawPotCircPntsMin.back()->At(iHit))->SetY1(itwp->GetCellCenter().y()/CLHEP::cm);
                                ((TEllipse *) drawPotCircPntsMin.back()->At(iHit))->SetR1(itwp->GetCellRad()/CLHEP::cm);
                                ((TEllipse *) drawPotCircPntsMin.back()->At(iHit))->SetR2( ( itwp->GetCellRad()/cos(itwp->GetWireEpsilon()) )/CLHEP::cm );
                                ((TEllipse *) drawPotCircPntsMin.back()->At(iHit))->SetTheta( itwp->GetCellCenter().phi()/CLHEP::degree );
                                ((TEllipse *) drawPotCircPntsMin.back()->At(iHit))->SetFillStyle(0);
                                ((TEllipse *) drawPotCircPntsMin.back()->At(iHit))->SetLineWidth(2);
                                ((TEllipse *) drawPotCircPntsMin.back()->At(iHit))->SetLineColor(iCol);
                                ((TEllipse *) drawPotCircPntsMin.back()->At(iHit))->Draw("same");
                                ++iHit;
                        }
                }

                tmpCircleFound->cd(2);
                hHitStereoPl.SetStats(kFALSE);
                hHitStereoPl.Draw();
                innerWall.Draw("same");
                outerWall.Draw("same");
                iCol = 0;
                std::vector<TEllipse *> drawPotCirclesPl;
                std::vector<TEllipse *> drawPotKrmCirclesPl;
                std::vector<TClonesArray *> drawPotCircPntsPl;
                for (circlesCont::iterator potCirclesPl_it = potCirclesPl.begin(); potCirclesPl_it != potCirclesPl.end(); ++potCirclesPl_it ) {
                        drawPotCirclesPl.push_back( new TEllipse( potCirclesPl_it->_center.x()/CLHEP::cm, potCirclesPl_it->_center.y()/CLHEP::cm, potCirclesPl_it->_radius/CLHEP::cm, potCirclesPl_it->_radius/CLHEP::cm) );
                        drawPotCirclesPl.back()->SetFillStyle(0);
                        drawPotCirclesPl.back()->SetLineWidth(2);
                        drawPotCirclesPl.back()->SetLineColor(++iCol);
                        drawPotCirclesPl.back()->SetLineStyle(4);
                        //drawPotCirclesPl.back()->Draw("same");
                        CLHEP::Hep2Vector radDir;
                        radDir.setX( cos(potCirclesPl_it->_krmCircFit.phi) );
                        radDir.setY( sin(potCirclesPl_it->_krmCircFit.phi) );
                        radDir.rotate( ( (potCirclesPl_it->_krmCircFit.rho>=0.0) ? -90.0 : 90.0 )*CLHEP::degree );
                        radDir=radDir.unit();
                        HepGeom::Point3D<double> CirCenter (0.0,0.0,0.0);
                        float krmRad = 1.0/fabs(potCirclesPl_it->_krmCircFit.rho);
                        CirCenter = CirCenter + (krmRad-potCirclesPl_it->_krmCircFit.dca)*radDir;
                        drawPotKrmCirclesPl.push_back( new TEllipse( CirCenter.x()/CLHEP::cm, CirCenter.y()/CLHEP::cm, krmRad/CLHEP::cm, krmRad/CLHEP::cm) );
                        //radDir *= (krmRad-potCirclesPl_it->_krmCircFit.dca);
                        //drawPotKrmCirclesMin.push_back( new TEllipse( radDir.x()/CLHEP::cm, radDir.y()/CLHEP::cm, krmRad/CLHEP::cm, krmRad/CLHEP::cm) );
                        drawPotKrmCirclesPl.back()->SetFillStyle(0);
                        drawPotKrmCirclesPl.back()->SetLineWidth(2);
                        drawPotKrmCirclesPl.back()->SetLineColor(iCol);
                        drawPotKrmCirclesPl.back()->SetLineStyle(4);
                        drawPotKrmCirclesPl.back()->Draw("same");

                        drawPotCircPntsPl.push_back( new TClonesArray("TEllipse") );
                        drawPotCircPntsPl.back()->ExpandCreateFast(potCirclesPl_it->_listHitptrs.size());
                        int iHit = 0;
                        for ( std::vector< ptrHitInClosCust >::iterator potCirclesPntsPl_it = potCirclesPl_it->_listHitptrs.begin(); potCirclesPntsPl_it != potCirclesPl_it->_listHitptrs.end(); ++potCirclesPntsPl_it ) {
                                itwp->SelectCell( potCirclesPntsPl_it->getRadLayID(), potCirclesPntsPl_it->getinLayerCellID() );
                                ((TEllipse *) drawPotCircPntsPl.back()->At(iHit))->SetX1(itwp->GetCellCenter().x()/CLHEP::cm);
                                ((TEllipse *) drawPotCircPntsPl.back()->At(iHit))->SetY1(itwp->GetCellCenter().y()/CLHEP::cm);
                                ((TEllipse *) drawPotCircPntsPl.back()->At(iHit))->SetR1(itwp->GetCellRad()/CLHEP::cm);
                                ((TEllipse *) drawPotCircPntsPl.back()->At(iHit))->SetR2( ( itwp->GetCellRad()/cos(itwp->GetWireEpsilon()) )/CLHEP::cm );
                                ((TEllipse *) drawPotCircPntsPl.back()->At(iHit))->SetTheta( itwp->GetCellCenter().phi()/CLHEP::degree );
                                ((TEllipse *) drawPotCircPntsPl.back()->At(iHit))->SetFillStyle(0);
                                ((TEllipse *) drawPotCircPntsPl.back()->At(iHit))->SetLineWidth(2);
                                ((TEllipse *) drawPotCircPntsPl.back()->At(iHit))->SetLineColor(iCol);
                                ((TEllipse *) drawPotCircPntsPl.back()->At(iHit))->Draw("same");
                                ++iHit;
                        }
               }
                tmpCircleFound->Modified();
                tmpCircleFound->Update();


                cerr << "Double click in the canvas_Fake to continue:" ;
                _fakeCanvas->cd();
                TLatex *printEvN = new TLatex(0.15,0.4,Form("Current Event: %d - Tpeak %d",event.id().event(),ipeak));
                printEvN->SetTextFont(62);
                printEvN->SetTextSizePixels(160);
                printEvN->Draw();
                _fakeCanvas->Modified();
                _fakeCanvas->Update();
                _fakeCanvas->WaitPrimitive();
                cerr << endl;
                delete printEvN;

                for (std::vector<confMapDraw*>::iterator hCmapMin_it = rCMap_Min.begin(); hCmapMin_it != rCMap_Min.end(); ++hCmapMin_it ) {
                        delete (*hCmapMin_it);
                }
                delete tmpRCP_Min;
                for (std::vector<confMapDraw*>::iterator hCmapPl_it = rCMap_Pl.begin(); hCmapPl_it != rCMap_Pl.end(); ++hCmapPl_it ) {
                        delete (*hCmapPl_it);
                }
                delete tmpRCP_Pl;

                rCMap_Min.clear();
                rCMap_Pl.clear();

                delete tmpCircleFound;
                for (std::vector<TEllipse *>::iterator drawPotCircles_it = drawPotCirclesMin.begin(); drawPotCircles_it != drawPotCirclesMin.end(); ++drawPotCircles_it ){
                        delete *drawPotCircles_it;
                }
                for (std::vector<TEllipse *>::iterator drawPotCircles_it = drawPotCirclesPl.begin(); drawPotCircles_it != drawPotCirclesPl.end(); ++drawPotCircles_it ){
                        delete *drawPotCircles_it;
                }

                for (std::vector<TClonesArray *>::iterator drawPotCircPntsMin_it = drawPotCircPntsMin.begin(); drawPotCircPntsMin_it != drawPotCircPntsMin.end(); ++drawPotCircPntsMin_it ) {
                        (*drawPotCircPntsMin_it)->Delete();
                }
                for (std::vector<TClonesArray *>::iterator drawPotCircPntsPl_it = drawPotCircPntsPl.begin(); drawPotCircPntsPl_it != drawPotCircPntsPl.end(); ++drawPotCircPntsPl_it ) {
                        (*drawPotCircPntsPl_it)->Delete();
                }

        }

    }

    if (isStereoMin!=0x0) delete [] isStereoMin;


    //event.put(sscc);

    //for (unsigned long i=0; i<1000000; i++) cout<<"lost time "<<endl;
    clock_t stopClock = clock();
    _hClockCicles->Fill((unsigned long)(stopClock-startClock));
    _hExecTime->Fill( (float)(stopClock-startClock)/((float) CLOCKS_PER_SEC ) );
    cout<<"-------- N clok to analyze 1 ev by ITTrackReco "<<stopClock-startClock<<" @ "<<CLOCKS_PER_SEC<<endl;

    /*if (_doDisplay) {
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
            //_peaksCanvHistos->Delete();
    }*/

//    _hPkStDistTrs->Delete();
////    _hPkStDistanceTrs->Delete();
//    _hPkStDist2W->Delete();
//    _hPkStClustDist2W->Delete();
//    _hPkSecDist2W->Delete();
//    _hPkSecClustDist2W->Delete();
//    _hPkSecVsStationDist2W->Delete();
//    _hPkSecVsStationDist2WGood->Delete();
//    _hPkAccArrSecStation->Delete();

//    if (tmpStClustDist!=0x0) delete tmpStClustDist;
//    if (tmpSecClustDist!=0x0) delete tmpSecClustDist;

    //-------------------------------------------

  } // end produce

  void ITTrackReco::endJob(){
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
  }

  bool ITTrackReco::findCluster(int &did, int &lid, int &sid, size_t &iglbHit, std::set<size_t> &hitLoked, closClstcol &clusCont ,TrackerHitByID const* hitByID, TrackerHitTimeCluster const&  tclust, CellGeometryHandle *itwp){
          if ( hitLoked.find(iglbHit)==hitLoked.end() ){
                  unsigned long int tmpIndex;
                  //tmpIndex = si.asUint();
                  TrackerHitByID::const_iterator hitByID_it;
                  std::pair<TrackerHitByID::const_iterator, TrackerHitByID::const_iterator> rangeHitByID_it;

                  nCellPerLayer = itwp->GetITLayer()->nCells();
                  int tmpSid=0, lost=0, hitClSizeInLayer=1;

                  closHitClust *tmpClust = new closHitClust(sid, iglbHit);
                  hitLoked.insert( iglbHit );

                  bool lostHit=true;
                  for (int icel=1; icel<nCellPerLayer; icel++){
                          tmpSid = (sid+icel)%nCellPerLayer;
                          tmpIndex = itwp->computeDet(did,lid,tmpSid);
                          if ( hitByID->count(tmpIndex)>0 ){
                                  lostHit=true;
                                  rangeHitByID_it= hitByID->equal_range(tmpIndex);
                                  for (hitByID_it=rangeHitByID_it.first; hitByID_it!=rangeHitByID_it.second; ++hitByID_it){
                                          if ( hitLoked.find(hitByID_it->second.key())==hitLoked.end() && hitByID_it->second->time()>tclust._minHitTime && hitByID_it->second->time()<tclust._maxHitTime ){
                                                  ++hitClSizeInLayer;
                                                  hitLoked.insert( hitByID_it->second.key() );
                                                  tmpClust->addHit(tmpSid,icel,hitByID_it->second.key());
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
                                          if ( hitLoked.find(hitByID_it->second.key())==hitLoked.end() && hitByID_it->second->time()>tclust._minHitTime && hitByID_it->second->time()<tclust._maxHitTime ){
                                                  ++hitClSizeInLayer;
                                                  hitLoked.insert( hitByID_it->second.key() );
                                                  tmpClust->addHit(tmpSid,-icel,hitByID_it->second.key());
                                                  lostHit=false;
                                                  break;
                                          }
                                  }
                                  if (lostHit) { ++lost; }
                          }
                          else { ++lost; }
                          if (lost>0) { break; }
                  }
                  if(_diagLevel > 1){
                          cout<<"***** tmpClust ***** starting hitId "<<iglbHit<<" at SL "<<did<<" lay "<<lid<<" cell "<<sid<<endl;
                          cout<<"before "<<tmpClust->_nHit<<" min "<<tmpClust->_minCellID<<" max "<<tmpClust->_maxCellID<<" center "<<tmpClust->_centerCellID<<endl;
                          closClstcol::iterator tmpClust_it = ( clusCont.insert( boost::shared_ptr<mu2e::closHitClust>(tmpClust) ) ).first;
                          cout<<"after "<<(*tmpClust_it)->_nHit<<" min "<<(*tmpClust_it)->_minCellID<<" max "<<(*tmpClust_it)->_maxCellID<<" center "<<(*tmpClust_it)->_centerCellID<<endl;
                          cout<<"are they equal? "<<((*tmpClust)==(*(*tmpClust_it)))<<endl;
                          cout<<"are they not equal? "<<((*tmpClust)!=(*(*tmpClust_it)))<<endl;
                  }
                  return true;
          }
          else return false;
  }

  bool ITTrackReco::findClosClust( unsigned int tmpX, unsigned int tmpY, std::map<unsigned int, bool> &notLooked, CHTVotArr &cmapCHTVotArr, ClosClustCHTVot &votClus ) {
          unsigned int closTmpX, closTmpY, closTmpuID;
          bool addHit = false;
          for (int ix=-1; ix<2; ++ix) {
                  closTmpX=tmpX+ix;
                  for (int iy=-1; iy<2; ++iy) {
                          if (ix==0 && iy==0) { continue; }
                          closTmpY=tmpY+iy;
                          if ( cmapCHTVotArr.xyTouID(closTmpX,closTmpY,closTmpuID) ) {
                                  std::set<unsigned int>::const_iterator closOverThrs_it = cmapCHTVotArr._overTHRs.find(closTmpuID);
                                  if ( closOverThrs_it!=cmapCHTVotArr._overTHRs.end() && notLooked[*closOverThrs_it] ) {
                                          notLooked[*closOverThrs_it]=false;
                                          votClus.addhit(closTmpX, closTmpY, *closOverThrs_it, cmapCHTVotArr.get(closTmpX,closTmpY));
                                          addHit = true;
                                          findClosClust( closTmpX, closTmpY, notLooked, cmapCHTVotArr, votClus );
                                  }
                          }
                  }
          }
          return addHit;
  }

  bool ITTrackReco::findCircles(closClinRadLay &radCls, CellGeometryHandle *itwp, std::vector<confMapDraw*> &rCMap, circlesCont & potcircles, std::string type, unsigned int thrCHTaddVot) {
          int nPotBending=0;
          bool foundCircles = false;

          float refX, refY, refRes2, closeCntCellX, closeCntCellY, cmapU, cmapV, TmpR2, invTmpR2, minGoodR=150.0, maxGoodR=450.0 ;
          int closeCntCellID, closeCntCellShift, halfwayCellID;
          float cmapU_clsCntCell, cmapV_clsCntCell, CHTtanTheta0_DX, CHTtanTheta0_DY, CHTtanTheta0, CHTTheta0, CHTr0, CHTcosArcTanT0;
          unsigned int votArrBin, halfwayVotArrBin;
          CHTVotArr_HitPtrrel votArr_HitPtrrel;
          ptrHitInClosCust refHitptr, clsCntCellHitptr;
          closeCntCellShift = (_minClusMultBend+1)/4+1;
          halfwayCellID = closeCntCellShift/2;
          if (halfwayCellID==0) halfwayCellID=1;
          if (halfwayCellID==closeCntCellShift) halfwayCellID=-1;

          float tmpCirChi2;

          //loop over all one view hits clusters, starting by the external radius, and make a CHT respect to the center of the starting cluster
          for ( closClinRadLay::iterator radCls_it = radCls.begin(); radCls_it != radCls.end(); ++radCls_it ) {
                  for ( closClstcol::iterator closCls_it = radCls_it->second.begin(); closCls_it != radCls_it->second.end(); ++closCls_it ) {
                          if ((*closCls_it)->_nHit>_minClusMultBend){
                                  ++nPotBending;
                                  std::vector<confMapPoint > confMap;
                                  std::vector<std::pair<double, double> > CHTconfMap;  // Hough Transform value (theta, r) for the point of the Conformal Mapping
                                  CHTVotArr cmapCHTVotArr;
                                  cmapCHTVotArr.setTHR(thrCHTaddVot);

                                  itwp->SelectCell(radCls_it->first, (*closCls_it)->_centerCellID);
                                  refX = itwp->GetCellCenter().getX();
                                  refY = itwp->GetCellCenter().getY();
                                  refRes2 = itwp->GetCellRad()*itwp->GetCellRad()/3.0;
                                  closeCntCellID = (*closCls_it)->_centerCellID+closeCntCellShift;//-1;
                                  nCellPerLayer = itwp->GetITLayer()->nCells();
                                  if (closeCntCellID<0) closeCntCellID+=nCellPerLayer;
                                  closeCntCellID=closeCntCellID%nCellPerLayer;
                                  itwp->SelectCell(radCls_it->first, closeCntCellID);
                                  closeCntCellX = itwp->GetCellCenter().getX();
                                  closeCntCellY = itwp->GetCellCenter().getY();
                                  cmapU_clsCntCell = closeCntCellX - refX;
                                  cmapV_clsCntCell = closeCntCellY - refY;
                                  //cout<<"------ ref (ID,x,y) "<<(*closCls_it)->_centerCellID<<", "<<refX<<", "<<refY<<" close "<<closeCntCellID<<", "<<closeCntCellX<<", "<<closeCntCellY<<endl;

                                  invTmpR2  = 1.0/(cmapU_clsCntCell*cmapU_clsCntCell+cmapV_clsCntCell*cmapV_clsCntCell);
                                  cmapU_clsCntCell*= invTmpR2;
                                  cmapV_clsCntCell*= invTmpR2;

                                  halfwayCellID = (*closCls_it)->_centerCellID + halfwayCellID;
                                  if (halfwayCellID<0) halfwayCellID+=nCellPerLayer;
                                  halfwayCellID=halfwayCellID%nCellPerLayer;

                                  for (hitsInClsID::const_iterator clHitsIDs_it=(*closCls_it)->_hitIdx.begin(); clHitsIDs_it!=(*closCls_it)->_hitIdx.end(); ++clHitsIDs_it) {
                                          if ((*closCls_it)->_centerCellID==clHitsIDs_it->first) {
                                                  refHitptr.setHitInClosCust(clHitsIDs_it,closCls_it,radCls_it);
                                                  continue;
                                          }
                                          itwp->SelectCell(itwp->GetSuperLayer(),itwp->GetCelRing(),clHitsIDs_it->first);
                                          cmapU = itwp->GetCellCenter().getX() - refX;
                                          cmapV = itwp->GetCellCenter().getY() - refY;
                                          TmpR2 = cmapU*cmapU+cmapV*cmapV;
                                          invTmpR2  = 1.0/TmpR2;
                                          cmapU*= invTmpR2;
                                          cmapV*= invTmpR2;
                                          confMap.push_back( confMapPoint(cmapU,cmapV,invTmpR2,/*15.0*invSqrt12*/itwp->GetCellRad()*invSqrt3) );
                                          if (closeCntCellID==clHitsIDs_it->first) {
                                                  clsCntCellHitptr.setHitInClosCust(clHitsIDs_it,closCls_it,radCls_it);
                                                  continue;
                                          }
                                          CHTtanTheta0_DX = cmapU-cmapU_clsCntCell;
                                          CHTtanTheta0_DY = cmapV_clsCntCell-cmapV;
                                          CHTtanTheta0 = CHTtanTheta0_DX / CHTtanTheta0_DY;
                                          CHTTheta0 = atan(CHTtanTheta0);//atan2(CHTtanTheta0_DX,CHTtanTheta0_DY);
                                          //CHTcosArcTanT0 = 1.0/sqrt(1.0+CHTtanTheta0*CHTtanTheta0);
                                          //CHTr0 = cmapU_clsCntCell*CHTcosArcTanT0 + cmapV_clsCntCell*CHTtanTheta0*CHTcosArcTanT0;
                                          //CHTr0 = cmapU_clsCntCell*cos(CHTTheta0) + cmapV_clsCntCell*sin(CHTTheta0);
                                          CHTr0 = cmapU*cos(CHTTheta0) + cmapV*sin(CHTTheta0);
                                          CHTconfMap.push_back(std::pair<double, double>( CHTTheta0/*atan2(CHTtanTheta0_DY,CHTtanTheta0_DX)*/, CHTr0 ) );
                                          //votArrBin = cmapCHTVotArr.Fill( CHTconfMap.back().first, CHTconfMap.back().second );
                                          //patch to avoid that the point of the ref Cluster will not be selected!!! FIXME
                                          //The points that should be in the same bin of halfwayCellID will be added to all CHT vot clusters that will be found!!! FIXME
                                          cmapCHTVotArr.xyTouID( CHTconfMap.back().first, CHTconfMap.back().second, votArrBin );
                                          if (halfwayCellID==clHitsIDs_it->first) {
                                                  halfwayVotArrBin = votArrBin;
                                          }
                                          votArr_HitPtrrel.insert(
                                                          CHTVotArr_HitPtrrel::value_type(
                                                                          votArrBin,
                                                                          ptrHitInClosCust(clHitsIDs_it,closCls_it,radCls_it)
                                                                          )
                                          );
                                          //cout<<"Hit of the ref cluster goes into hit of CHT vot "<<votArrBin<<endl;
                                  }
                                  closClinRadLay::iterator sndRadCls_it = radCls_it;
                                  ++sndRadCls_it;
                                  for ( ; sndRadCls_it != radCls.end(); ++sndRadCls_it ) {
                                          itwp->SelectCell(sndRadCls_it->first, 0);
                                          for ( closClstcol::iterator sndClosCls_it = sndRadCls_it->second.begin(); sndClosCls_it != sndRadCls_it->second.end(); ++sndClosCls_it ) {
                                                  for (hitsInClsID::const_iterator clHitsIDs_it=(*sndClosCls_it)->_hitIdx.begin(); clHitsIDs_it!=(*sndClosCls_it)->_hitIdx.end(); ++clHitsIDs_it) {
                                                          itwp->SelectCell(itwp->GetSuperLayer(),itwp->GetCelRing(),clHitsIDs_it->first);
                                                          cmapU = itwp->GetCellCenter().getX() - refX;
                                                          cmapV = itwp->GetCellCenter().getY() - refY;
                                                          TmpR2 = cmapU*cmapU+cmapV*cmapV;
                                                          invTmpR2  = 1.0/TmpR2;
                                                          cmapU*= invTmpR2;
                                                          cmapV*= invTmpR2;
                                                          confMap.push_back( confMapPoint(cmapU,cmapV,invTmpR2,/*15.0*invSqrt12*/itwp->GetCellRad()*invSqrt3) );

                                                          CHTtanTheta0_DX = cmapU-cmapU_clsCntCell;
                                                          CHTtanTheta0_DY = cmapV_clsCntCell-cmapV;
                                                          CHTtanTheta0 = CHTtanTheta0_DX / CHTtanTheta0_DY;
                                                          CHTTheta0 = atan(CHTtanTheta0);//atan2(CHTtanTheta0_DX,CHTtanTheta0_DY);
                                                          //CHTcosArcTanT0 = 1.0/sqrt(1.0+CHTtanTheta0*CHTtanTheta0);
                                                          //CHTr0 = cmapU_clsCntCell*CHTcosArcTanT0 + cmapV_clsCntCell*CHTtanTheta0*CHTcosArcTanT0;
                                                          //CHTr0 = cmapU_clsCntCell*cos(CHTTheta0) + cmapV_clsCntCell*sin(CHTTheta0);
                                                          CHTr0 = cmapU*cos(CHTTheta0) + cmapV*sin(CHTTheta0);
                                                          CHTconfMap.push_back(std::pair<double, double>( CHTTheta0/*atan2(CHTtanTheta0_DY,CHTtanTheta0_DX)*/, CHTr0 ) );
                                                          votArrBin = cmapCHTVotArr.Fill( CHTconfMap.back().first, CHTconfMap.back().second );
                                                          votArr_HitPtrrel.insert(
                                                                          CHTVotArr_HitPtrrel::value_type(
                                                                                          votArrBin,
                                                                                          ptrHitInClosCust(clHitsIDs_it,sndClosCls_it,sndRadCls_it)
                                                                                          )
                                                          );
                                                 }
                                          }
                                  }

                                  //cout<<"CHT vot hit over thr "<<cmapCHTVotArr._overTHRs.size()<<endl;
                                  //group closest over threshold votes
                                  std::map<unsigned int, bool> notLookedClos;
                                  std::vector<ClosClustCHTVot> closClustCHTVots;
                                  for ( std::set<unsigned int>::iterator overThrs_it = cmapCHTVotArr._overTHRs.begin(); overThrs_it != cmapCHTVotArr._overTHRs.end(); ++overThrs_it ){
                                          notLookedClos.insert(std::pair<unsigned int, bool>(*overThrs_it,true) );
                                  }
                                  unsigned int tmpCHTvot_X=0, tmpCHTvot_Y=0;
                                  for ( std::set<unsigned int>::iterator overThrs_it = cmapCHTVotArr._overTHRs.begin(); overThrs_it != cmapCHTVotArr._overTHRs.end(); ++overThrs_it ){
                                          if ( notLookedClos[*overThrs_it] ) {
                                                  notLookedClos[*overThrs_it] = false;
                                                  cmapCHTVotArr.uIDToxy(*overThrs_it, tmpCHTvot_X, tmpCHTvot_Y);
                                                  ClosClustCHTVot tmpVotClus;
                                                  tmpVotClus.addhit(tmpCHTvot_X, tmpCHTvot_Y, *overThrs_it, cmapCHTVotArr.get(tmpCHTvot_X,tmpCHTvot_Y));
                                                  findClosClust( tmpCHTvot_X, tmpCHTvot_Y, notLookedClos, cmapCHTVotArr, tmpVotClus );
                                                  //std::set<ClosClustCHTVot>::iterator closClustCHTVots_it = closClustCHTVots.insert(tmpVotClus).first;

                                                  //patch to avoid that the point of the ref Cluster will not be selected!!! FIXME
                                                  //The points that should be in the same bin of halfwayCellID will be added to all CHT vot clusters that will be found!!! FIXME
                                                  cmapCHTVotArr.uIDToxy(halfwayVotArrBin, tmpCHTvot_X, tmpCHTvot_Y);
                                                  tmpVotClus.addhit(tmpCHTvot_X, tmpCHTvot_Y, halfwayVotArrBin, votArr_HitPtrrel.count(halfwayVotArrBin));

                                                  closClustCHTVots.push_back(tmpVotClus);
                                          }
                                  }

                                  //evaluate the circles by the clusters and store them if found
                                  float tmpTheta, tmpCHT_SinTheta, tmpCHT_CosTheta2, tmpinvCHT_r, tmpCHT_sigmaTheta2, tmpCHT_sigmar2;
                                  float tmpCMap_B, tmpCMap_A, tmpCMap_B2, tmpCMap_A2, tmpCMap_sigmaB2, tmpCMap_sigmaA2;
                                  float tmpRad;
                                  float cellAveRes;
                                  for ( std::vector<ClosClustCHTVot>::iterator closClustCHTVots_it=closClustCHTVots.begin(); closClustCHTVots_it!=closClustCHTVots.end(); ++closClustCHTVots_it ) {

                                          SimpleCircle2D tmpCircle;
                                          tmpTheta            = closClustCHTVots_it->_meanX*cmapCHTVotArr.getBinWdtX() + cmapCHTVotArr.getMinX();
                                          tmpCHT_SinTheta     = sin( tmpTheta );
                                          tmpCHT_CosTheta2    = 1.0 - tmpCHT_SinTheta*tmpCHT_SinTheta;
                                          tmpinvCHT_r         = 1.0/( closClustCHTVots_it->_meanY*cmapCHTVotArr.getBinWdtY() + cmapCHTVotArr.getMinY() );
                                          //cout<<"---- Theta "<<tmpTheta<<" r "<<1.0/tmpinvCHT_r<<endl;
                                          tmpCHT_sigmaTheta2  = closClustCHTVots_it->_sigmaX*cmapCHTVotArr.getBinWdtX();
                                          tmpCHT_sigmaTheta2 *= tmpCHT_sigmaTheta2;
                                          tmpCHT_sigmar2      = closClustCHTVots_it->_sigmaY*cmapCHTVotArr.getBinWdtY();
                                          tmpCHT_sigmar2     *= tmpCHT_sigmar2;
                                          tmpCMap_B           = 0.5 * tmpCHT_SinTheta * tmpinvCHT_r;
                                          tmpCMap_B2          = tmpCMap_B * tmpCMap_B;
                                          tmpCMap_sigmaB2     = tmpinvCHT_r * tmpinvCHT_r * ( 0.25*tmpCHT_CosTheta2*tmpCHT_sigmaTheta2 + tmpCMap_B2*tmpCHT_sigmar2 );
                                          tmpCMap_sigmaA2     = 1.0 / tan(tmpTheta);
                                          tmpCMap_A           = tmpCMap_B * tmpCMap_sigmaA2;
                                          tmpCMap_A2          = tmpCMap_A * tmpCMap_A;
                                          tmpCMap_sigmaA2    *= tmpCMap_sigmaA2;
                                          tmpCMap_sigmaA2    *= ( tmpCMap_sigmaB2  +  tmpCMap_A2*tmpCHT_sigmaTheta2 / ( tmpCHT_CosTheta2*tmpCHT_CosTheta2 ) );

                                          tmpRad = sqrt( tmpCMap_A2 + tmpCMap_B2 );
                                          if (tmpRad>=minGoodR && tmpRad<=maxGoodR ) {
                                                  tmpCircle.SetCenter( tmpCMap_A + refX, tmpCMap_B + refY, sqrt( tmpCMap_sigmaA2 + refRes2 ), sqrt( tmpCMap_sigmaB2 + refRes2 ) );
                                                  tmpCircle.SetRadius( tmpRad, sqrt( tmpCMap_A2*tmpCMap_sigmaA2 + tmpCMap_B2*tmpCMap_sigmaB2 )/tmpRad );
                                                  //for ( std::set<unsigned int>::iterator votInClus_it = closClustCHTVots_it->_listCHTVotIDs.begin();  votInClus_it != closClustCHTVots_it->_listCHTVotIDs.end(); ++votInClus_it ) {
                                                  //        cmapCHTVotArr.uIDToxy(*votInClus_it, tmpCHTvot_X, tmpCHTvot_Y);
                                                  //        cout<<" ("<<tmpCHTvot_X<<", "<<tmpCHTvot_Y<<") with val "<<cmapCHTVotArr.get(tmpCHTvot_X, tmpCHTvot_Y);
                                                  //}
                                                  tmpCircle._listHitptrs.push_back(refHitptr);//insert(refHitptr);
                                                  tmpCircle._listHitptrs.push_back(clsCntCellHitptr);//insert(clsCntCellHitptr);

                                                  /*itwp->SelectCell(refHitptr.getRadLayID(), refHitptr.getinLayerCellID());
                                                  cellAveRes = itwp->GetCellRad() * invSqrt3;
                                                  tmpCircle._krmCircFit.addHit( itwp->GetCellCenter().getX(), itwp->GetCellCenter().getY(), cellAveRes, cellAveRes );
                                                  itwp->SelectCell(clsCntCellHitptr.getRadLayID(), clsCntCellHitptr.getinLayerCellID());
                                                  cellAveRes = itwp->GetCellRad() * invSqrt3;
                                                  tmpCircle._krmCircFit.addHit( itwp->GetCellCenter().getX(), itwp->GetCellCenter().getY(), cellAveRes, cellAveRes );*/
                                                  cellAveRes = sqrt(refRes2);
                                                  tmpCircle._krmCircFit.addHit( refX, refY, cellAveRes, cellAveRes );
                                                  tmpCircle._krmCircFit.addHit( closeCntCellX, closeCntCellY, cellAveRes, cellAveRes );

                                                  //cout<<"For Potential circles"<<endl;
                                                  for (std::set< unsigned int >::const_iterator listCHTVotIDs_it=closClustCHTVots_it->_listCHTVotIDs.begin();
                                                                  listCHTVotIDs_it!=closClustCHTVots_it->_listCHTVotIDs.end(); ++listCHTVotIDs_it) {
                                                          //cout<<"bins in clusters "<<*listCHTVotIDs_it<<endl;
                                                          pair<CHTVotArr_HitPtrrel::iterator,CHTVotArr_HitPtrrel::iterator>
                                                                  votArr_HitPtrrel_its = votArr_HitPtrrel.equal_range(*listCHTVotIDs_it);
                                                          for ( CHTVotArr_HitPtrrel::iterator tmpCHTVotArr_HitPtrrel_it=votArr_HitPtrrel_its.first;
                                                                          tmpCHTVotArr_HitPtrrel_it!=votArr_HitPtrrel_its.second;
                                                                          ++tmpCHTVotArr_HitPtrrel_it ) {
                                                                  tmpCircle._listHitptrs.push_back(tmpCHTVotArr_HitPtrrel_it->second);//insert(tmpCHTVotArr_HitPtrrel_it->second);
                                                                  itwp->SelectCell(tmpCHTVotArr_HitPtrrel_it->second.getRadLayID(), tmpCHTVotArr_HitPtrrel_it->second.getinLayerCellID());
                                                                  cellAveRes = itwp->GetCellRad() * invSqrt3;
                                                                  tmpCircle._krmCircFit.addHit( itwp->GetCellCenter().getX(), itwp->GetCellCenter().getY(), cellAveRes, cellAveRes );
                                                         }
                                                  }
                                                  foundCircles = true;
                                                  if ( tmpCircle._krmCircFit.computeBestCirc(/*closeCntCellX, closeCntCellY*/) ) {
                                                          tmpCirChi2 = tmpCircle._krmCircFit.chicir/((float) (tmpCircle._krmCircFit.np()-3));
                                                          if ( tmpCirChi2>0.001 && tmpCirChi2<_circChi2cut ) {
                                                                  potcircles.push_back( tmpCircle );
                                                          }
                                                  }
                                          }
                                  }

                                  if (_diagLevel>1) {
                                          int icmpel=0;
                                          cout<<"for pot bend "<<nPotBending<<" :"<<endl;
                                          for (std::vector<confMapPoint >::iterator confMap_it=confMap.begin(); confMap_it!=confMap.end(); ++confMap_it) {
                                                  cout<<"\t u "<<confMap_it->_u<<" v "<<confMap_it->_v<<" r "<<sqrt(confMap_it->_rSq)<<endl;
                                                  ++icmpel;
                                          }

                                          cout<<"Point over threshold in CHT voting array of the conformal mapping: "<<endl;
                                          cout<<"\t points: ";
                                          for ( std::set<unsigned int>::iterator overThrs_it = cmapCHTVotArr._overTHRs.begin(); overThrs_it != cmapCHTVotArr._overTHRs.end(); ++overThrs_it ){
                                                  cmapCHTVotArr.uIDToxy(*overThrs_it, tmpCHTvot_X, tmpCHTvot_Y);
                                                  cout<<" ("<<tmpCHTvot_X<<", "<<tmpCHTvot_Y<<") with val "<<cmapCHTVotArr.get(tmpCHTvot_X, tmpCHTvot_Y);
                                          }
                                          cout<<endl;

                                          int clusOfVot = 0;
                                          //for ( std::set<ClosClustCHTVot>::iterator closClustCHTVots_it=closClustCHTVots.begin(); closClustCHTVots_it!=closClustCHTVots.end(); ++closClustCHTVots_it ) {
                                          for ( std::vector<ClosClustCHTVot>::iterator closClustCHTVots_it=closClustCHTVots.begin(); closClustCHTVots_it!=closClustCHTVots.end(); ++closClustCHTVots_it ) {
                                                  ++clusOfVot;
                                                  cout<<"clusOfVot "<<clusOfVot<<", meanX "<<closClustCHTVots_it->_meanX<<", meanY "<<closClustCHTVots_it->_meanY
                                                                  <<", sigmaX "<<closClustCHTVots_it->_sigmaX<<", sigmaY "<<closClustCHTVots_it->_sigmaY<<endl;
                                                  cout<<"\t points: ";
                                                  for ( std::set<unsigned int>::iterator votInClus_it = closClustCHTVots_it->_listCHTVotIDs.begin();  votInClus_it != closClustCHTVots_it->_listCHTVotIDs.end(); ++votInClus_it ) {
                                                          cmapCHTVotArr.uIDToxy(*votInClus_it, tmpCHTvot_X, tmpCHTvot_Y);
                                                          cout<<" ("<<tmpCHTvot_X<<", "<<tmpCHTvot_Y<<") with val "<<cmapCHTVotArr.get(tmpCHTvot_X, tmpCHTvot_Y);
                                                  }
                                                  cout<<endl;
                                          }

                                  }
                                  if (_doDisplay) {
                                          //rCMap.push_back(new TH2F( Form("hRCMapHitStereo%s_%i",type.c_str(),nPotBending), "Relative Conformal Mapping of Stereo- Hits per Event at z=0", 2000, -2, 2, 2000, -2, 2 ));
                                          float *Uarr=new float[confMap.size()];
                                          float *Varr=new float[confMap.size()];
                                          float *errUarr=new float[confMap.size()];
                                          float *errVarr=new float[confMap.size()];
                                          int icmpel=0;
                                          for (std::vector<confMapPoint >::iterator confMap_it=confMap.begin(); confMap_it!=confMap.end(); ++confMap_it) {
                                                  //rCMap.back()->Fill(confMap_it->_u,confMap_it->_v);
                                                  Uarr[icmpel] = confMap_it->_u;
                                                  Varr[icmpel] = confMap_it->_v;
                                                  errUarr[icmpel] = 0.0;//confMap_it->_errU;
                                                  errVarr[icmpel] = 0.0;//confMap_it->_errV;
                                                  ++icmpel;
                                          }
                                          rCMap.push_back(new confMapDraw());
                                          rCMap.back()->_cmap = new TGraphErrors( icmpel, Uarr, Varr, errUarr, errVarr );
                                          delete Uarr;
                                          delete Varr;
                                          delete errUarr;
                                          delete errVarr;

                                          rCMap.back()->_cmapCHT = new TH2F( Form("hRCHT_CMapHitStereo%s_%i",type.c_str(),nPotBending), Form("Relative CHT of Conformal Mapping of Stereo%s Hits per Event at z=0",type.c_str()), 2513, -TMath::Pi(), TMath::Pi(), 500, -0.25, 0.25 );
                                          for (std::vector<std::pair<double, double> >::iterator CHTconfMap_it = CHTconfMap.begin(); CHTconfMap_it != CHTconfMap.end(); ++CHTconfMap_it){
                                                  rCMap.back()->_cmapCHT->Fill(CHTconfMap_it->first,CHTconfMap_it->second);
                                          }
                                  }
                          }
                  }
          }

          return foundCircles;
  }

  void ITTrackReco::mergeCircles( circlesCont &potCircles ) {
          typedef std::pair<circlesCont::iterator,circlesCont::iterator> circlePair;
          std::vector<circlePair > matchingCirc;
          std::set<circlesCont::iterator, greater<circlesCont::iterator> > circToErase;

          //Check if two or more circles match into one
          for ( circlesCont::iterator potCircles_it = potCircles.begin(); potCircles_it != potCircles.end(); ++potCircles_it ) {
                  circlesCont::iterator potCircles2_it = potCircles_it;
                  ++potCircles2_it;
                  for (; potCircles2_it != potCircles.end(); ++potCircles2_it ) {
                          //if (potCircles_it->isCompatibleByCmapWith(*potCircles2_it,2.0)) {
                          //if (potCircles_it->isCompatibleByKrmWith(*potCircles2_it,10.0)) {
                                 matchingCirc.push_back( circlePair(potCircles_it,potCircles2_it) );
                          //}
                  }
          }
          //Merge the matching circles into one
          cout<<"Circle couples found "<<matchingCirc.size()<<endl;
          if ( matchingCirc.size()>0 ) {
                  bool removeCirc = false;
                  int nAddedPoints = 0;
                  for ( std::vector<circlePair >::iterator matchingCirc_it = matchingCirc.begin(); matchingCirc_it != matchingCirc.end(); ++matchingCirc_it ) {
                          nAddedPoints = 0;
                          //if ( matchingCirc_it->first->_listHitptrs.size()>=matchingCirc_it->second->_listHitptrs.size() ) {
                                  nAddedPoints = matchingCirc_it->first->mergeCirc( *matchingCirc_it->second, removeCirc);
                                  if (removeCirc) {
                                          circToErase.insert(matchingCirc_it->second);
                                  }
                          /*} else {
                                  nAddedPoints = matchingCirc_it->second->mergeCirc ( *matchingCirc_it->first, removeCirc);
                                  if (removeCirc) {
                                          circToErase.insert(matchingCirc_it->first);
                                  }
                          }*/


                          /*if (matchingCirc_it->first->_listHitptrs.begin()->getinEventHitID()==matchingCirc_it->second->_listHitptrs.begin()->getinEventHitID()) {
                                  matchingCirc_it->first->summCirc(*matchingCirc_it->second,2);
                                  circToErase.insert(matchingCirc_it->second);
                                  cout<<"Merging couple with equal ref"<<endl;
                          } else {
                                  bool refPntNotShared = true;
                                  for (circlesCont::iterator potCircles_it = potCircles.begin(); potCircles_it != potCircles.end(); ++potCircles_it ) {
                                          if ( potCircles_it == matchingCirc_it->first || potCircles_it == matchingCirc_it->second ) {continue;}
                                          //cout<<"ref hit id i-th circ "<<potCircles_it->_listHitptrs.begin()->getinEventHitID()<<" first in coue "<<matchingCirc_it->first->_listHitptrs.begin()->getinEventHitID();
                                          //cout<<" second in couple "<<matchingCirc_it->second->_listHitptrs.begin()->getinEventHitID()<<endl;
                                          if ( potCircles_it->_listHitptrs.begin()->getinEventHitID() == matchingCirc_it->first->_listHitptrs.begin()->getinEventHitID() ) {
                                                  //cout<<"First has duplicate ref point"<<endl;
                                                  matchingCirc_it->second->summCirc(*matchingCirc_it->first,2);
                                                  circToErase.insert(matchingCirc_it->first);
                                                  cout<<"Merging couple with different ref but ref of first was scared"<<endl;
                                                  refPntNotShared = false;
                                                  break;
                                          }
                                          if ( potCircles_it->_listHitptrs.begin()->getinEventHitID() == matchingCirc_it->second->_listHitptrs.begin()->getinEventHitID() ) {
                                                  //cout<<"Second has duplicate ref point"<<endl;
                                                  matchingCirc_it->first->summCirc(*matchingCirc_it->second,2);
                                                  circToErase.insert(matchingCirc_it->second);
                                                  cout<<"Merging couple with different ref but ref of second was scared"<<endl;
                                                  refPntNotShared = false;
                                                  break;
                                          }

                                  }
                                  if (refPntNotShared) {
                                          //is it right? FIXME
                                          *matchingCirc_it->first += *matchingCirc_it->second;
                                          circToErase.insert(matchingCirc_it->second);
                                          cout<<"Merging couple with different ref"<<endl;
                                  }
                          }*/
                  }

                  //clear Circles that have been added
                  for (std::set<circlesCont::iterator >::iterator circToErase_it = circToErase.begin(); circToErase_it != circToErase.end(); ++circToErase_it ) {
                          potCircles.erase(*circToErase_it);
                  }
          }
  }

  void ITTrackReco::printPotCircles( circlesCont &potCircles ) {
          if (_diagLevel>1) {
                  cout<<"Found "<< potCircles.size() <<" potential circles: "<<endl;
                  for (circlesCont::iterator potCircles_it = potCircles.begin(); potCircles_it != potCircles.end(); ++potCircles_it ) {
                          cout<<"\t Center "<<potCircles_it->_center<<" Rasius "<<potCircles_it->_radius
                                          <<" sigmas (x,y,r) "<<potCircles_it->_sigmaCenter[0]<<", "<<potCircles_it->_sigmaCenter[1]<<", "<<potCircles_it->_sigmaRad<<endl;
                          cout<<"\t by Karimaki fit: rho "<<potCircles_it->_krmCircFit.rho<<" +- "<<sqrt(potCircles_it->_krmCircFit.covrfd[0])
                                          <<" phi "<< potCircles_it->_krmCircFit.phi<<" +- "<< sqrt(potCircles_it->_krmCircFit.covrfd[2])
                                          <<" dca "<< potCircles_it->_krmCircFit.dca<<" +- "<< sqrt(potCircles_it->_krmCircFit.covrfd[5])
                                          <<" chi2/ndof "<< potCircles_it->_krmCircFit.chicir/((float) (potCircles_it->_krmCircFit.np()-3))<<endl;
                          /*circlesCont::iterator potCircles2_it = potCircles_it;
                          ++potCircles2_it;
                          for (; potCircles2_it != potCircles.end(); ++potCircles2_it ) {
                                  if (potCircles_it->isCompatibleByCmapWith(*potCircles2_it)) {
                                          cout<<"found match at 1 sigma"<<endl;
                                          if (potCircles_it->_listHitptrs.begin()->getinEventHitID() == potCircles2_it->_listHitptrs.begin()->getinEventHitID() ) {
                                                  cout<<"Same ref point"<<endl;
                                          }
                                  }
                                  if (potCircles_it->isCompatibleByCmapWith(*potCircles2_it,2.0)) {
                                          cout<<"found match at 2 sigma"<<endl;
                                          if (potCircles_it->_listHitptrs.begin()->getinEventHitID() == potCircles2_it->_listHitptrs.begin()->getinEventHitID() ) {
                                                  cout<<"Same ref point"<<endl;
                                          }
                                  }
                                  if (potCircles_it->isCompatibleByCmapWith(*potCircles2_it,3.0)) {
                                          cout<<"found match at 3 sigma"<<endl;
                                          if (potCircles_it->_listHitptrs.begin()->getinEventHitID() == potCircles2_it->_listHitptrs.begin()->getinEventHitID() ) {
                                                  cout<<"Same ref point"<<endl;
                                          }
                                  }
                          }*/
                          if (_diagLevel>3) {
                          cout<<"\t Points in circle : "<<potCircles_it->_listHitptrs.size()<<endl;
                                  for ( /*std::set*/std::vector< ptrHitInClosCust >::iterator hitsInCircles_it=potCircles_it->_listHitptrs.begin();
                                                  hitsInCircles_it!=potCircles_it->_listHitptrs.end(); ++hitsInCircles_it ) {
                                          cout<<"\t\t"<<*hitsInCircles_it<<endl;
                                  }
                          }
                  }
          }
  }

}  // end namespace mu2e

using mu2e::ITTrackReco;
DEFINE_ART_MODULE(ITTrackReco);
