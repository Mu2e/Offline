//
// Fast Patter recognition for the ITracker
//
// $Id: ITTrackReco_module.cc,v 1.20 2013/10/21 21:01:22 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/10/21 21:01:22 $
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
#include "RecoDataProducts/inc/HelixVal.hh"
#include "RecoDataProducts/inc/TrackSeed.hh"
#include "RecoDataProducts/inc/TrackSeedCollection.hh"
#include "FastPatternReco/inc/FastPatRecoUtilsAndDataDef.hh"
#include "FastPatternReco/inc/TrkHelixFitIT.hh"
#include "FastPatternReco/inc/TrackSeedUtils.hh"

//temp El Visual
#include "MCDataProducts/inc/VisibleGenElTrack.hh"
#include "MCDataProducts/inc/VisibleGenElTrackCollection.hh"


//For track fit
// BaBar
#include "BTrk/BaBar/BaBar.hh"
#include "KalmanTests/inc/TrkDef.hh"
//#include "KalmanTests/inc/TrkStrawHit.hh"
//#include "KalmanTests/inc/KalFit.hh"
//#include "KalmanTests/inc/KalFitMC.hh"
//#include "TrkPatRec/inc/TrkHitFilter.hh"
//#include "TrkPatRec/inc/HelixFit.hh"
//#include "BTrk/TrkBase/TrkPoca.hh"
//#include "TrkPatRec/src/TrkPatRec_module.cc"
#include "TrkPatRec/inc/TrkPatRec.hh"



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
#define sqrt2 1.414213562

using namespace std;

namespace mu2e {

//typedef std::map<std::vector<XYZP>::iterator, std::pair<size_t, size_t> > RbstFitPntsLpsPntsRel;
typedef std::map<size_t, std::pair<size_t, size_t> > RbstFitPntsLpsPntsRel;


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


class ITTrackReco : public art::EDProducer {

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
        void produce(art::Event & e);
        //void analyze(art::Event const& e);

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

        float  _minGoodR, _maxGoodR;
        float  _minZ, _maxZ;
        double _rIn;
        double _rOut;
        bool   _hitsInUpstream;

        // Use n times sigma for matching the clusters in the Sector plane
        //float _nSigmaForClMatch;

        // Use n times sigma for matching the clusters slope with 0
        //float _nSigmaForClSlope;

        // Diagnostics level.
        int _diagLevel;

        bool _doDisplay;

        //temp El Visual
        // Label of the module that made the hits.
        std::string _extractElectronsData;


        // robust helix fitter
        TrkHelixFitIT _hfit;


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

        double _maxTdrift;
        double _tpeakerr;

        //void SctrSttnMapAnalyze(TH2I *startDist, TH2I *votingArray, int minPitch=3, int maxPitch=11);
        //rwclinrwrel rwClst_forRw_rel;
        //std::vector<Clust> clustersList;
        bool *notAssociatedHits;
        size_t nAssociatedHit;
        hitCrossingList hitCrossingChecked;
        hitCrossingList hitCrossingCheckedT;
        hitCrossingList segmantHitCrossingChecked;
        hitCrossingList segmantHitCrossingCheckedT;

        bool  findCluster(int &did, int &lid, int &sid, bool isUpStream, size_t &iglbHit, std::set<size_t> &hitLoked, closClstcol &clusCont, TrackerHitByID const* hitByID, TrackerHitTimeCluster const&  tclust, CellGeometryHandle *itwp);
        bool  findClosClust( unsigned int tmpX, unsigned int tmpY, std::map<unsigned int, bool> &notLooked, CHTVotArr &cmapCHTVotArr, ClosClustCHTVot &votClus );
        bool  findCircles(closClinRadLay &radCls, CellGeometryHandle *itwp, circlesCont & potcircles, std::string type);
        bool  findCircles(closClinRadLay &radCls, CellGeometryHandle *itwp, std::vector<confMapDraw*> &rCMap, circlesCont & potcircles, std::string type="", int minClusMultBend =6, unsigned int thrCHTaddVot=5, bool skipAssociated=false);
        void  mergeCircles( circlesCont &potCircles );
        void  markUsedHit( circlesCont &potCircles );
        void  printPotCircles( circlesCont &potCircles );

        bool  findCircles3D(closClinRadLay &radCls_1v, closClinRadLay &radCls_2v, CellGeometryHandle *itwp, std::vector<points3D> &potLoops/*circlesCont & potcircles*/, int minClusMultBend=10, bool elForwDirection=true);
        inline bool findWireCross(int crossCellId1, int crossCellId2, closClinRadLay::iterator &radCls_1v_it, closClstcol::iterator &closCls_it, hitsInClsID::const_iterator &clHitsIDs_it,  closClinRadLay::iterator &radCls_2v_it, closClstcol::iterator &lowRadClosCls_it, hitsInClsID::const_iterator &lowRadClHitsIDs_it, int matchZdir, float prevZVal, CellGeometryHandle *itwp, points3D &potLoopPoints);
        bool  iteratMaxMinWireCross(closClinRadLay &radCls_1v, closClinRadLay &radCls_2v, closClinRadLay::iterator &radCls_1v_it, closClstcol::iterator &closCls_it, closClinRadLay::iterator &radCls_2v_it, closClstcol::iterator &lowRadClosCls_it, int matchZdir, float prevZVal, CellGeometryHandle *itwp, points3D &potLoopPoints);
        bool  iteratMinMaxWireCross(closClinRadLay &radCls_1v, closClinRadLay &radCls_2v, closClinRadLay::iterator &radCls_1v_it, closClstcol::iterator &closCls_it, closClinRadLay::iterator &radCls_2v_it, closClstcol::iterator &lowRadClosCls_it, int matchZdir, float prevZVal, CellGeometryHandle *itwp, points3D &potLoopPoints);
        int   findZclusters( std::vector<XYZP> &xyzp, float zdim=20.0, float zPeakSigma=4.0, float zPeakThr=0.01, TH1F *zHist=0x0 );
        int   findZclusters( HelixTraj &trkhelixe, std::vector<XYZP> &xyzp, std::vector< std::vector<size_t> > &pntIdsLoops/*, std::vector<TrkHelix> &helixeFitLoops, std::vector<HelixTraj> &trkhelixeLoops*/);

        void  searchTracks(closClinRadLay &radCls_Pl, closClinRadLay &radCls_Min, unique_ptr<TrackSeedCollection> &seeds,
                        CellGeometryHandle *itwp, std::vector<TrkTimePeak> &tpeaks, size_t &ipeak, art::Handle<TrackerHitTimeClusterCollection> &tclustHandle, art::Handle<StrawHitCollection> &pdataHandle
                        , VisibleGenElTrackCollection const* genEltrks
                        );

        // Some ugly but necessary ROOT related bookkeeping:

        // The job needs exactly one instance of TApplication.  See note 1.
        unique_ptr<TApplication> _application;

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
                    _minGoodR(150.0),
                    _maxGoodR(450.0),
                    _minZ(0.0),
                    _maxZ(0.0),
                    _rIn(0.0),
                    _rOut(0.0),
                    _hitsInUpstream(false),
                    //_nSigmaForClMatch(pset.get<float>("nSigmaForClMatch")),
                    //_nSigmaForClSlope(pset.get<float>("nSigmaForClSlope")),
                    /*_nAccumulate(pset.get<int>("nAccumulate",20)),*/
                    _diagLevel(pset.get<int>("diagLevel",0)),
                    _doDisplay(pset.get<bool>("doDisplay",false)),
                    _extractElectronsData(pset.get<string>("elextractModuleLabel")),//temp El Visual

                    _hfit(pset),

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

                    notAssociatedHits(0),

                    // Some ugly but necessary ROOT related bookkeeping.
                    _application(nullptr)
//    _directory(0)
{
        cout<<"Constructed"<<endl;
        //if ( _nSigmaForClMatch<0.0 ) _nSigmaForClMatch*=-1.00000;
        //if ( _nSigmaForClSlope<0.0 ) _nSigmaForClSlope*=-1.00000;

        // Tell the framework what we make.
        produces<TrackSeedCollection>();

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
                        _application = unique_ptr<TApplication>(new TApplication( "noapplication", &tmp_argc, tmp_argv ));
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

void ITTrackReco::produce(art::Event & event ) {
//void ITTrackReco::analyze(art::Event const& event ) {

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

        unique_ptr<TrackSeedCollection> seeds(new TrackSeedCollection);


        const Tracker& tracker = getTrackerOrThrow();
        const ITracker &itr = static_cast<const ITracker&>( tracker );
        CellGeometryHandle *itwp = itr.getCellGeometryHandle();
        _maxZ =  itr.maxEndCapDim()-itr.getWalls()->find(Wall::endcap)->second->getTotalThickness();
        _minZ = -_maxZ;

        _maxTdrift = 200.0;


        //_hfit._itwp = itwp;

        _rIn  = itr.r0()+itr.getWalls()->find(Wall::inner)->second->getTotalThickness();
        _rOut = itr.rOut()-itr.getWalls()->find(Wall::outer)->second->getTotalThickness();

        stMaprel ssmap_Bin_Straw_rel;           //relation between straw and the AbsSector (Rotation) vs Station map
        stMaprel::iterator ssmpbnstrw_rel_it;   //iterator for relation between straw and the AbsSector (Rotation) vs Station map

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


        //temp El Visual
        art::Handle<VisibleGenElTrackCollection> genEltrksHandle;
        event.getByLabel(_extractElectronsData,genEltrksHandle);
        VisibleGenElTrackCollection const* genEltrks = genEltrksHandle.product();



        clock_t startClock = clock();

        size_t nStrawPerEvent = hits->size();
        size_t nTimeClusPerEvent = tclusts->size();
        cout<<"----------------------------------------------------------------------------------"<<endl;
        cout<<"event "<<event.id().event()<<" tot N hit "<<nStrawPerEvent<<" N peak found "<<nTimeClusPerEvent<<endl;
        cout<<"----------------------------------------------------------------------------------"<<endl;
        int i1peak;

        int sid, lid, did, radId;
        unsigned int ihit;

        //double refX, refY, /*tmpX, tmpY,*/ cmapU, cmapV, invTmpR2;
        size_t iglbHit=0;
        bool *isStereoMin = new (nothrow) bool [nStrawPerEvent];
        for (size_t jhit=0; jhit<nStrawPerEvent; ++jhit) {
                isStereoMin[jhit] = false;
        }

        //TH2I *tmpStClustDist    = new TH2I("tmpStClustDist","tmp Smoothing of Device vs Straw multiplicity Distribution",36,0,36,1000,-200,800);
        //TH2I *tmpSecClustDist   = new TH2I("tmpSecClustDist","tmp Smoothing of Device vs Sector multiplicity Distribution",36,0,36,20,-4,16);
        //TH2I *ihPkStDistTrs, *ihPkStDist2W, *ihPkStClustDist2W, *ihPkSecDist2W, *ihPkSecClustDist2W, *ihPkSecVsStationDist2WG, *ihPkSecVsStationDist2W, *ihPkAccArrSecStation;
        //TCanvas * iCanvHist=0x0;

        //int nCellPerLayer = 0;
        int nClus=0, nTotHitInCls=0, nTotHitInCls_Pl=0, nTotHitInCls_MIn=0, nTotHit_Pl=0, nTotHit_MIn=0;
        //closClstcol
        closClinRadLay radCls_Pl, radCls_Min, radCls_Pl_Ups, radCls_Min_Ups;
        closClinRadLay::iterator radCls_it, sndRadCls_it;

        notAssociatedHits = new bool[nStrawPerEvent];
        for (size_t jhit=0; jhit<nStrawPerEvent; ++jhit) { notAssociatedHits[jhit]=true; }

        //find for ciclres using Dave's trkpaterc code
        std::vector<TrkTimePeak> tpeaks;

        for (size_t ipeak=0; ipeak<nTimeClusPerEvent; ipeak++) {

                nAssociatedHit=0;
                hitCrossingChecked.clear();
                hitCrossingCheckedT.clear();

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
                        itwp->SelectCellDet(si.asUint());
                        //itwp->SelectCell(did,lid,sid);
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

                                if (itwp->isDownStream()) {
                                        radCls_it=radCls_Min.find(radId);
                                        if (radCls_it==radCls_Min.end()) {
                                                radCls_it=( radCls_Min.insert( closClinRadLay::value_type( radId, closClstcol() ) ) ).first;
                                        }
                                        findCluster(did, lid, sid, false, iglbHit, hitLoked, radCls_it->second, hitByID, tclust, itwp);
                                } else if (itwp->isUpStream()) {
                                        radCls_it=radCls_Min_Ups.find(radId);
                                        if (radCls_it==radCls_Min_Ups.end()) {
                                                radCls_it=( radCls_Min_Ups.insert( closClinRadLay::value_type( radId, closClstcol() ) ) ).first;
                                        }
                                        findCluster(did, lid, sid, true, iglbHit, hitLoked, radCls_it->second, hitByID, tclust, itwp);
                                }
                        }
                        else {
                                //isStereoMin[iglbHit]=false;
                                ++nTotHit_Pl;

                                if (itwp->isDownStream()) {
                                        radCls_it=radCls_Pl.find(radId);
                                        if (radCls_it==radCls_Pl.end()) {
                                                radCls_it=( radCls_Pl.insert( closClinRadLay::value_type( radId, closClstcol() ) ) ).first;
                                        }
                                        findCluster(did, lid, sid, false, iglbHit, hitLoked, radCls_it->second, hitByID, tclust, itwp);
                                } else if (itwp->isUpStream()) {
                                        radCls_it=radCls_Pl_Ups.find(radId);
                                        if (radCls_it==radCls_Pl_Ups.end()) {
                                                radCls_it=( radCls_Pl_Ups.insert( closClinRadLay::value_type( radId, closClstcol() ) ) ).first;
                                        }
                                        findCluster(did, lid, sid, false, iglbHit, hitLoked, radCls_it->second, hitByID, tclust, itwp);
                                }
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

                        int nTotHitInCls_MIn_Ups=0, nTotHitInCls_Pl_Ups=0;
                        if (itr.isDumbbell()) {
                                cout<<"Clusters for Min Upstream streo:"<<endl;

                                for ( closClinRadLay::iterator radCls_it = radCls_Min_Ups.begin(); radCls_it != radCls_Min_Ups.end(); ++radCls_it ) {
                                        cout<<"***** cluster for Radial layer "<<radCls_it->first<<" :"<<endl;
                                        for ( closClstcol::iterator closCls_it = radCls_it->second.begin(); closCls_it != radCls_it->second.end(); ++closCls_it ) {
                                                cout<<"-- cluster n "<<nClus<<endl;
                                                ++nClus;
                                                nTotHitInCls_MIn_Ups+=(*closCls_it)->_nHit;
                                                nTotHitInCls+=(*closCls_it)->_nHit;
                                                cout<<"\t cluster size "<<(*closCls_it)->_nHit<<" minCell "<<(*closCls_it)->_minCellID<<" maxCell "<<(*closCls_it)->_maxCellID<<" centralCell "<<(*closCls_it)->_centerCellID<<endl;
                                                cout<<"\t hits in cluster: ";
                                                for (hitsInClsID::const_iterator clsHit_it = (*closCls_it)->_hitIdx.begin(); clsHit_it != (*closCls_it)->_hitIdx.end(); ++clsHit_it){
                                                        cout<<"("<<clsHit_it->first<<", "<<clsHit_it->second<<") ";
                                                }
                                                cout<<endl;
                                        }
                                }
                                cout<<"Clusters for Plus Upstream streo:"<<endl;

                                for ( closClinRadLay::iterator radCls_it = radCls_Pl_Ups.begin(); radCls_it != radCls_Pl_Ups.end(); ++radCls_it ) {
                                        cout<<"***** cluster for Radial layer "<<radCls_it->first<<" :"<<endl;
                                        for ( closClstcol::iterator closCls_it = radCls_it->second.begin(); closCls_it != radCls_it->second.end(); ++closCls_it ) {
                                                cout<<"-- cluster n "<<nClus<<endl;
                                                ++nClus;
                                                nTotHitInCls_Pl_Ups+=(*closCls_it)->_nHit;
                                                nTotHitInCls+=(*closCls_it)->_nHit;
                                                cout<<"\t cluster size "<<(*closCls_it)->_nHit<<" minCell "<<(*closCls_it)->_minCellID<<" maxCell "<<(*closCls_it)->_maxCellID<<" centralCell "<<(*closCls_it)->_centerCellID<<endl;
                                                cout<<"\t hits in cluster: ";
                                                for (hitsInClsID::const_iterator clsHit_it = (*closCls_it)->_hitIdx.begin(); clsHit_it != (*closCls_it)->_hitIdx.end(); ++clsHit_it){
                                                        cout<<"("<<clsHit_it->first<<", "<<clsHit_it->second<<") ";
                                                }
                                                cout<<endl;
                                        }
                                }
                        }

                        cout<<"In Time peak :"<<endl;
                        cout<<" n. Hit for Stereo Min "<<nTotHit_MIn<<" n. Hit for Stereo Pl "<<nTotHit_Pl<<endl;
                        cout<<"NClus "<<nClus<<" n. Tot hit grouped in clusters "<<nTotHitInCls<<" n. Hit in Clusters for Stereo Min "<<nTotHitInCls_MIn<<" n. Hit in Clusters for Stereo Pl "<<nTotHitInCls_Pl;
                        if (itr.isDumbbell()) {cout<<" n. Hit in Clusters for Stereo Min Ups "<<nTotHitInCls_MIn_Ups<<" n. Hit in Clusters for Stereo Pl Ups "<<nTotHitInCls_Pl_Ups;}
                        cout<<endl;
                }

                _tpeakerr = (tclust._maxHitTime - tclust._minHitTime)*invSqrt12;

                cout<<"------- start Tracks search for event "<<event.id().event()<<" ------"<<endl;
                //isHitIDUsed notAssociatedHits;
                //notAssociatedHits.clear();
                /*rCMap_Pl.clear();
                rCMap_Min.clear();
                rCMap_Pl2.clear();
                rCMap_Min2.clear();*/

                _hitsInUpstream=false;
                searchTracks(radCls_Pl, radCls_Min, seeds, itwp, tpeaks, ipeak, tclustHandle, pdataHandle, genEltrks);
                if (itr.isDumbbell()) {
                        _hitsInUpstream=true;
                        searchTracks(radCls_Pl_Ups, radCls_Min_Ups, seeds, itwp, tpeaks, ipeak, tclustHandle, pdataHandle, genEltrks);
                }


        }

        if (isStereoMin!=0x0) delete [] isStereoMin;
        if (notAssociatedHits!=0x0) delete [] notAssociatedHits;


        event.put(std::move(seeds));

        //for (unsigned long i=0; i<1000000; i++) cout<<"lost time "<<endl;
        clock_t stopClock = clock();
        _hClockCicles->Fill((unsigned long)(stopClock-startClock));
        _hExecTime->Fill( (float)(stopClock-startClock)/((float) CLOCKS_PER_SEC ) );
        cout<<"-------- N clok to analyze 1 ev by ITTrackReco "<<stopClock-startClock<<" @ "<<CLOCKS_PER_SEC<<endl;

        if (_doDisplay) {
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
        }

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

bool ITTrackReco::findCluster(int &did, int &lid, int &sid, bool isUpStream, size_t &iglbHit, std::set<size_t> &hitLoked, closClstcol &clusCont ,TrackerHitByID const* hitByID, TrackerHitTimeCluster const&  tclust, CellGeometryHandle *itwp){
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
                        tmpIndex = itwp->computeDet(did,lid,tmpSid,isUpStream);
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
                        tmpIndex = itwp->computeDet(did,lid,tmpSid,isUpStream);
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
                        cout<<"***** tmpClust ***** starting hitId "<<iglbHit<<" at SL "<<did<<" lay "<<lid<<" cell "<<sid<<" isUpstream "<<isUpStream<<endl;
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

bool ITTrackReco::findCircles3D(closClinRadLay &radCls_1v, closClinRadLay &radCls_2v, CellGeometryHandle *itwp, std::vector<points3D> &potLoops/*circlesCont & potcircles*/, int minClusMultBend, bool elForwDirection) {
        int nPotBending=0;
        bool foundCircles = false;

        //float refX, refY, refRes2, closeCntCellX, closeCntCellY, cmapU, cmapV, TmpR2, invTmpR2;
        //int closeCntCellID, closeCntCellShift, halfwayCellID;
        //float cmapU_clsCntCell, cmapV_clsCntCell, CHTtanTheta0_DX, CHTtanTheta0_DY, CHTtanTheta0, CHTTheta0, CHTr0, CHTcosArcTanT0;
        //unsigned int votArrBin, halfwayVotArrBin;
        //CHTVotArr_HitPtrrel votArr_HitPtrrel;
        //ptrHitInClosCust refHitptr, clsCntCellHitptr;
        //closeCntCellShift = (_minClusMultBend+1)/4+1;
        //halfwayCellID = closeCntCellShift/2;
        //if (halfwayCellID==0) halfwayCellID=1;
        //if (halfwayCellID==closeCntCellShift) halfwayCellID=-1;

        //float cellAveRes, cellAveZRes;
        //float tmpCirChi2, tmpRad;
        //bool firstInsert = false;
        bool foundWrCross=false;
        //float zCorssPos, wrsDistAtz, minZCorssPos, minWrsDistAtz, cutWrsDistAtz;
        int crossCellId1, crossCellId2;
        int tmpAbsRad;
        int matchZdir=1;

        //std::map<int, std::vector<std::> > wrsCrossingMap
        //notAssociatedHits

        //loop over all one view hits clusters, starting by the external radius, and make a CHT respect to the center of the starting cluster
        for ( closClinRadLay::iterator radCls_1v_it = radCls_1v.begin(); radCls_1v_it != radCls_1v.end(); ++radCls_1v_it ) {
                for ( closClstcol::iterator closCls_it = radCls_1v_it->second.begin(); closCls_it != radCls_1v_it->second.end(); ++closCls_it ) {
                        if ((*closCls_it)->_nHit>((unsigned int)minClusMultBend)){
                                ++nPotBending;
                                potLoops.push_back(points3D());
                                points3D &potLoopPoints = potLoops.back();

                                itwp->SelectCell(radCls_1v_it->first, (*closCls_it)->_maxCellID,_hitsInUpstream);
                                //cellAveRes = itwp->GetCellRad() * invSqrt3;
                                //cellAveZRes = cellAveRes/cos( itwp->GetWireEpsilon() );
                                //cutWrsDistAtz = 2.5*itwp->GetCellRad();
                                tmpAbsRad = (radCls_1v_it->first-1);
                                closClinRadLay::iterator radCls_2v_it = radCls_2v.find( tmpAbsRad );
                                if (radCls_2v_it!=radCls_2v.end()) {
                                        for ( closClstcol::iterator lowRadClosCls_it = radCls_2v_it->second.begin(); lowRadClosCls_it != radCls_2v_it->second.end(); ++lowRadClosCls_it) {
                                                foundWrCross=false;
                                                /*cout<<"----------- ("<< radCls_1v_it->first<<","<< (*closCls_it)->_maxCellID<<") and ("<< tmpAbsRad<<","<< (*lowRadClosCls_it)->_minCellID<<") is? "<<itwp->canIntersectInZ(zCorssPos, wrsDistAtz, tmpAbsRad, (*lowRadClosCls_it)->_minCellID )<<endl;
                                                  if ( itwp->canIntersectInZ(zCorssPos, wrsDistAtz, tmpAbsRad, (*lowRadClosCls_it)->_minCellID ) ) {
                                                          cout<<"tmp z cross "<<zCorssPos<<" wires dist "<<wrsDistAtz<<endl;
                                                          if (wrsDistAtz<cutWrsDistAtz && wrsDistAtz<minWrsDistAtz) {
                                                                  minWrsDistAtz= wrsDistAtz;
                                                                  minZCorssPos=zCorssPos;
                                                                  crossCellId=(*lowRadClosCls_it)->_minCellID;
                                                                  foundWrCross=true;
                                                          }
                                                  }
                                                  cout<<"----------- ("<< radCls_1v_it->first<<","<< (*closCls_it)->_maxCellID<<") and ("<< tmpAbsRad<<","<< (*lowRadClosCls_it)->_maxCellID<<") is? "<<itwp->canIntersectInZ(zCorssPos, wrsDistAtz, tmpAbsRad, (*lowRadClosCls_it)->_maxCellID )<<endl;
                                                  if ( itwp->canIntersectInZ(zCorssPos, wrsDistAtz, tmpAbsRad, (*lowRadClosCls_it)->_maxCellID ) ) {
                                                          cout<<" tmp z cross "<<zCorssPos<<" wires dist "<<wrsDistAtz<<endl;
                                                          if (wrsDistAtz<cutWrsDistAtz && wrsDistAtz<minWrsDistAtz) {
                                                                  minWrsDistAtz= wrsDistAtz;
                                                                  minZCorssPos=zCorssPos;
                                                                  crossCellId=(*lowRadClosCls_it)->_maxCellID;
                                                                  foundWrCross=true;
                                                          }
                                                  }*/
                                                if (elForwDirection) {
                                                        crossCellId1 = (*closCls_it)->_maxCellID;
                                                        crossCellId2 = (*lowRadClosCls_it)->_minCellID;
                                                        hitsInClsID::const_iterator clHitsIDs_it=((*closCls_it)->_hitIdx.end()-1);
                                                        hitsInClsID::const_iterator lowRadClHitsIDs_it = (*lowRadClosCls_it)->_hitIdx.begin();
                                                        if (/*notAssociatedHits[clHitsIDs_it->second] && notAssociatedHits[lowRadClHitsIDs_it->second] &&*/
                                                                        hitCrossingChecked.find(std::pair<size_t,size_t>(clHitsIDs_it->second,lowRadClHitsIDs_it->second))==hitCrossingChecked.end() &&
                                                                        hitCrossingCheckedT.find(std::pair<size_t,size_t>(lowRadClHitsIDs_it->second,clHitsIDs_it->second))==hitCrossingCheckedT.end() )
                                                        {
                                                                foundWrCross = findWireCross( crossCellId1, crossCellId2, radCls_1v_it, closCls_it, clHitsIDs_it, radCls_2v_it, lowRadClosCls_it, lowRadClHitsIDs_it, 0, 0.0, itwp, potLoopPoints);
                                                                if (foundWrCross) {
                                                                        matchZdir = (_hitsInUpstream) ? -1 : 1;
                                                                        iteratMaxMinWireCross(radCls_1v, radCls_2v, radCls_1v_it, closCls_it, radCls_2v_it, lowRadClosCls_it, /*0*/matchZdir, potLoopPoints.back()._pos.z(), itwp, potLoopPoints);
                                                                }
                                                                if ((*closCls_it)->_nHit>3) {
                                                                        matchZdir = (_hitsInUpstream) ? 1 : -1;
                                                                        crossCellId1 = (*closCls_it)->_minCellID;
                                                                        crossCellId2 = (*lowRadClosCls_it)->_maxCellID;
                                                                        hitsInClsID::const_iterator sClHitsIDs_it=(*closCls_it)->_hitIdx.begin();
                                                                        hitsInClsID::const_iterator sLowRadClHitsIDs_it = ((*lowRadClosCls_it)->_hitIdx.end()-1);
                                                                        if (/*notAssociatedHits[sClHitsIDs_it->second] && notAssociatedHits[sLowRadClHitsIDs_it->second] &&*/
                                                                                        hitCrossingChecked.find(std::pair<size_t,size_t>(sClHitsIDs_it->second,sLowRadClHitsIDs_it->second))==hitCrossingChecked.end() &&
                                                                                        hitCrossingCheckedT.find(std::pair<size_t,size_t>(sLowRadClHitsIDs_it->second,sClHitsIDs_it->second))==hitCrossingCheckedT.end() )
                                                                        {
                                                                                foundWrCross = findWireCross( crossCellId1, crossCellId2, radCls_1v_it, closCls_it, sClHitsIDs_it, radCls_2v_it, lowRadClosCls_it, sLowRadClHitsIDs_it, 0, 0.0, itwp, potLoopPoints);
                                                                                if (foundWrCross) {
                                                                                        iteratMinMaxWireCross(radCls_1v, radCls_2v, radCls_1v_it, closCls_it, radCls_2v_it, lowRadClosCls_it, /*0*/matchZdir, potLoopPoints.back()._pos.z(), itwp, potLoopPoints);
                                                                                }
                                                                        }
                                                                }
                                                        }
                                                } else {
                                                        crossCellId1 = (*closCls_it)->_minCellID;
                                                        crossCellId2=(*lowRadClosCls_it)->_maxCellID;
                                                        hitsInClsID::const_iterator clHitsIDs_it=(*closCls_it)->_hitIdx.begin();
                                                        hitsInClsID::const_iterator lowRadClHitsIDs_it = ((*lowRadClosCls_it)->_hitIdx.end()-1);
                                                        if (/*notAssociatedHits[clHitsIDs_it->second] && notAssociatedHits[lowRadClHitsIDs_it->second] &&*/
                                                                        hitCrossingChecked.find(std::pair<size_t,size_t>(clHitsIDs_it->second,lowRadClHitsIDs_it->second))==hitCrossingChecked.end() &&
                                                                        hitCrossingCheckedT.find(std::pair<size_t,size_t>(lowRadClHitsIDs_it->second,clHitsIDs_it->second))==hitCrossingCheckedT.end() )
                                                        {
                                                                foundWrCross = findWireCross( crossCellId1, crossCellId2, radCls_1v_it, closCls_it, clHitsIDs_it, radCls_2v_it, lowRadClosCls_it, lowRadClHitsIDs_it, 0, 0.0, itwp, potLoopPoints);
                                                                if (foundWrCross) {
                                                                        matchZdir = (_hitsInUpstream) ? -1 : 1;
                                                                        iteratMinMaxWireCross(radCls_1v, radCls_2v, radCls_1v_it, closCls_it, radCls_2v_it, lowRadClosCls_it, /*0*/matchZdir, potLoopPoints.back()._pos.z(), itwp, potLoopPoints);
                                                                }
                                                                if ((*closCls_it)->_nHit>3) {
                                                                        matchZdir = (_hitsInUpstream) ? 1 : -1;
                                                                        crossCellId1 = (*closCls_it)->_maxCellID;
                                                                        crossCellId2=(*lowRadClosCls_it)->_minCellID;
                                                                        hitsInClsID::const_iterator sClHitsIDs_it=((*closCls_it)->_hitIdx.end()-1);
                                                                        hitsInClsID::const_iterator sLowRadClHitsIDs_it = (*lowRadClosCls_it)->_hitIdx.begin();
                                                                        if (/*notAssociatedHits[sClHitsIDs_it->second] && notAssociatedHits[sLowRadClHitsIDs_it->second] &&*/
                                                                                        hitCrossingChecked.find(std::pair<size_t,size_t>(sClHitsIDs_it->second,sLowRadClHitsIDs_it->second))==hitCrossingChecked.end() &&
                                                                                        hitCrossingCheckedT.find(std::pair<size_t,size_t>(sLowRadClHitsIDs_it->second,sClHitsIDs_it->second))==hitCrossingCheckedT.end() )
                                                                        {
                                                                                foundWrCross = findWireCross( crossCellId1, crossCellId2, radCls_1v_it, closCls_it, sClHitsIDs_it, radCls_2v_it, lowRadClosCls_it, sLowRadClHitsIDs_it, 0, 0.0, itwp, potLoopPoints);
                                                                                if (foundWrCross) {
                                                                                        iteratMaxMinWireCross(radCls_1v, radCls_2v, radCls_1v_it, closCls_it, radCls_2v_it, lowRadClosCls_it, /*0*/matchZdir, potLoopPoints.back()._pos.z(), itwp, potLoopPoints);
                                                                                }
                                                                        }
                                                                }
                                                        }

                                                }

                                                //if (foundWrCross) {
                                                //        cout<<"Found cross between (AbsRad,CellId) ("<<radCls_1v_it->first<<","<< (*closCls_it)->_maxCellID<<") and ("<<tmpAbsRad<<","<<crossCellId<<") "<<endl;
                                                //        cout<<"\t at z "<<minZCorssPos<<" distance between wires "<<minWrsDistAtz<<endl;
                                                //}
                                        }
                                }

                                if (potLoops.back().size()<10) {
                                        for (points3D::iterator points_it = potLoops.back().begin(); points_it != potLoops.back().end(); ++points_it ) {
                                                std::pair<size_t,size_t> clearCross(points_it->getInEventHitID(),points_it->getCrossInEventHitID());
                                                std::pair<size_t,size_t> clearCrossT(points_it->getCrossInEventHitID(),points_it->getInEventHitID());
                                                hitCrossingChecked.erase(clearCross);
                                                hitCrossingCheckedT.erase(clearCrossT);
                                       }
                                        potLoops.pop_back();
                                }

                                //                                  (*closCls_it)->_minCellID;
                                //
                                //                                  SimpleCircle2D tmpCircle;
                                //
                                //                                  itwp->SelectCell(radCls_1v_it->first, (*closCls_it)->_centerCellID, _hitsInUpstream);
                                //                                  refX = itwp->GetCellCenter().getX();
                                //                                  refY = itwp->GetCellCenter().getY();
                                //                                  refRes2 = itwp->GetCellRad()*itwp->GetCellRad()/3.0;
                                //                                  cellAveRes = sqrt(refRes2);
                                //                                  tmpCircle._krmCircFit.addHit(refX,refY,cellAveRes,cellAveRes);
                                //
                                //                                  closeCntCellID = (*closCls_it)->_centerCellID+closeCntCellShift;//-1;
                                //                                  nCellPerLayer = itwp->GetITLayer()->nCells();
                                //                                  if (closeCntCellID<0) closeCntCellID+=nCellPerLayer;
                                //                                  closeCntCellID=closeCntCellID%nCellPerLayer;
                                //                                  itwp->SelectCell(radCls_1v_it->first, closeCntCellID, _hitsInUpstream);
                                //                                  closeCntCellX = itwp->GetCellCenter().getX();
                                //                                  closeCntCellY = itwp->GetCellCenter().getY();
                                //
                                //                                  tmpCircle._krmCircFit.addHit(closeCntCellX,closeCntCellY,cellAveRes,cellAveRes);
                                //                                  //cout<<"------ ref (ID,x,y) "<<(*closCls_it)->_centerCellID<<", "<<refX<<", "<<refY<<" close "<<closeCntCellID<<", "<<closeCntCellX<<", "<<closeCntCellY<<endl;
                                //
                                //
                                //                                  halfwayCellID = (*closCls_it)->_centerCellID + halfwayCellID;
                                //                                  if (halfwayCellID<0) halfwayCellID+=nCellPerLayer;
                                //                                  halfwayCellID=halfwayCellID%nCellPerLayer;
                                //
                                //                                  firstInsert = false;
                                //
                                //                                  for (hitsInClsID::const_iterator clHitsIDs_it=(*closCls_it)->_hitIdx.begin(); clHitsIDs_it!=(*closCls_it)->_hitIdx.end(); ++clHitsIDs_it) {
                                //                                          if ((*closCls_it)->_centerCellID==clHitsIDs_it->first) {
                                //                                                  tmpCircle._listHitptrs.insert( tmpCircle._listHitptrs.begin(), ptrHitInClosCust(clHitsIDs_it,closCls_it,radCls_1v_it) );
                                //                                                  firstInsert = true;
                                //                                                  continue;
                                //                                          }
                                //                                          if (closeCntCellID==clHitsIDs_it->first) {
                                //                                                  if (firstInsert) {
                                //                                                          tmpCircle._listHitptrs.insert( (tmpCircle._listHitptrs.begin()+1), ptrHitInClosCust(clHitsIDs_it,closCls_it,radCls_1v_it) );
                                //                                                  }
                                //                                                  else {
                                //                                                          tmpCircle._listHitptrs.insert( tmpCircle._listHitptrs.begin(), ptrHitInClosCust(clHitsIDs_it,closCls_it,radCls_1v_it) );
                                //                                                  }
                                //                                                 continue;
                                //                                          }
                                //                                          if (halfwayCellID==clHitsIDs_it->first) {
                                //                                          }
                                //
                                //                                          itwp->SelectCell(itwp->GetSuperLayer(),itwp->GetCelRing(),clHitsIDs_it->first, _hitsInUpstream);
                                //                                          tmpCircle._krmCircFit.addHit(itwp->GetCellCenter().getX(),itwp->GetCellCenter().getY(),cellAveRes,cellAveRes);
                                //                                          tmpCircle._listHitptrs.push_back( ptrHitInClosCust(clHitsIDs_it,closCls_it,radCls_1v_it) );
                                //
                                //                                          //confMap.push_back( confMapPoint(cmapU,cmapV,invTmpR2,/*15.0*invSqrt12*/itwp->GetCellRad()*invSqrt3) );
                                //                                          /*votArr_HitPtrrel.insert(
                                //                                                          CHTVotArr_HitPtrrel::value_type(
                                //                                                                          votArrBin,
                                //                                                                          ptrHitInClosCust(clHitsIDs_it,closCls_it,radCls_1v_it)
                                //                                                                          )
                                //                                          );*/
                                //                                          //cout<<"Hit of the ref cluster goes into hit of CHT vot "<<votArrBin<<endl;
                                //                                  }
                                //                                  if ( !tmpCircle._krmCircFit.computeBestCirc() ) {
                                //                                          cout<<"Circle not good"<<endl;
                                //                                          break;////////////////FIXME
                                //                                  }
                                //
                                //                                  closClinRadLay::iterator sndRadCls_1v_it = radCls_1v_it;
                                //                                  ++sndRadCls_1v_it;
                                //                                  for ( ; sndRadCls_1v_it != radCls_1v.end(); ++sndRadCls_1v_it ) {
                                //                                          itwp->SelectCell(sndRadCls_1v_it->first, 0, _hitsInUpstream);
                                //                                          cellAveRes = itwp->GetCellRad()*invSqrt3;
                                //
                                //                                          for ( closClstcol::iterator sndClosCls_it = sndRadCls_1v_it->second.begin(); sndClosCls_it != sndRadCls_1v_it->second.end(); ++sndClosCls_it ) {
                                //                                                  for (hitsInClsID::const_iterator clHitsIDs_it=(*sndClosCls_it)->_hitIdx.begin(); clHitsIDs_it!=(*sndClosCls_it)->_hitIdx.end(); ++clHitsIDs_it) {
                                //                                                          itwp->SelectCell(itwp->GetSuperLayer(),itwp->GetCelRing(),clHitsIDs_it->first, _hitsInUpstream);
                                //                                                          if ( tmpCircle._krmCircFit.checkBfrAddHit(itwp->GetCellCenter().getX(),itwp->GetCellCenter().getY(),cellAveRes,cellAveRes) ) {
                                //                                                                  tmpCircle._listHitptrs.push_back( ptrHitInClosCust(clHitsIDs_it,sndClosCls_it,sndRadCls_1v_it) );
                                //                                                                  cout<<"Point added"<<endl;
                                //                                                          }
                                //                                                 }
                                //                                          }
                                //                                  }
                                //
                                //                                  if ( tmpCircle._krmCircFit.ierror==0 ) {
                                //                                          tmpCirChi2 = tmpCircle._krmCircFit.chicir/((float) (tmpCircle._krmCircFit.np()-3));
                                //                                          if ( tmpCirChi2>0.001 && tmpCirChi2<_circChi2cut ) {
                                //                                                  cout<<"Circ Chi2 good"<<endl;
                                //                                                  tmpRad = 1.0/tmpCircle._krmCircFit.rho;
                                //                                                  //if (tmpRad>=_minGoodR && tmpRad<=_maxGoodR ) {
                                //                                                          cout<<"Circle added"<<endl;
                                //                                                          potcircles.push_back( tmpCircle );
                                //                                                          foundCircles = true;
                                //                                                  //}
                                //                                          }
                                //                                  }


                                /*
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
                                 */
                        }
                }
        }

        return foundCircles;
}

inline bool ITTrackReco::findWireCross(int crossCellId1, int crossCellId2, closClinRadLay::iterator &radCls_1v_it, closClstcol::iterator &closCls_it, hitsInClsID::const_iterator &clHitsIDs_it,  closClinRadLay::iterator &radCls_2v_it, closClstcol::iterator &lowRadClosCls_it, hitsInClsID::const_iterator &lowRadClHitsIDs_it, int matchZdir, float prevZVal, CellGeometryHandle *itwp, points3D &potLoopPoints) {
        itwp->SelectCell(radCls_1v_it->first, crossCellId1, _hitsInUpstream);
        float zCorssPos=0.0, wrsDistAtz=0.0, cutWrsDistAtz;
        float cellAveRes, cellAveZRes;
        cellAveRes = itwp->GetCellRad() * invSqrt3;
        cellAveZRes = sqrt2 * cellAveRes/cos( itwp->GetWireEpsilon() );
        cutWrsDistAtz = 2.5*itwp->GetCellRad();
        bool foundWrCross=false;
        bool condition = false;
        float deltaZ, deltaZRes = sqrt2*cellAveZRes;

        //cout<<"Looking cross between (AbsRad,CellId) ("<<radCls_1v_it->first<<","<< crossCellId1<<") and ("<<radCls_2v_it->first<<","<<crossCellId2<<") "<<endl;

        //if ( itwp->canIntersectInZ(zCorssPos, wrsDistAtz, radCls_2v_it->first, crossCellId2 ) ) {
        itwp->SelectComp_Cell(radCls_2v_it->first, crossCellId2, _hitsInUpstream);
        if ( itwp->canIntersectInZ(zCorssPos, wrsDistAtz ) ) {
                condition = wrsDistAtz<cutWrsDistAtz;
                if (matchZdir!=0) { //matchZdir = 0 disable z matching
                        deltaZ = zCorssPos-prevZVal;
                        if (matchZdir>0) {
                                condition =  condition && ( deltaZ>-7.0*deltaZRes && deltaZ<15.0*deltaZRes );
                                //condition =  condition && ( deltaZ>-3.0*deltaZRes && deltaZ<10.0*deltaZRes );
                        } else {
                                condition =  condition && ( deltaZ>-15.0*deltaZRes && deltaZ<7.0*deltaZRes );
                                //condition =  condition && ( deltaZ>-10.0*deltaZRes && deltaZ<3.0*deltaZRes );
                        }
                }

                if (condition) {
                        potLoopPoints.push_back( corssingPoints() );
                        itwp->WirePosAtZ( zCorssPos, potLoopPoints.back()._pos );
                        potLoopPoints.back()._sigmax = cellAveRes;
                        potLoopPoints.back()._sigmay = cellAveRes;
                        potLoopPoints.back()._sigmaz = cellAveZRes;

                        potLoopPoints.back()._hitRef.setHitInClosCust(clHitsIDs_it,closCls_it,radCls_1v_it);

                        potLoopPoints.back()._crossHitRef.setHitInClosCust(lowRadClHitsIDs_it,lowRadClosCls_it,radCls_2v_it);
                        hitCrossingChecked.insert( std::pair<size_t,size_t>( clHitsIDs_it->second, lowRadClHitsIDs_it->second ) );
                        hitCrossingCheckedT.insert( std::pair<size_t,size_t>( lowRadClHitsIDs_it->second, clHitsIDs_it->second ) );
                        //cout<<"\t found at z "<<zCorssPos<<" distance between wires "<<wrsDistAtz<<endl;
                        foundWrCross=true;
                }
        }

        return foundWrCross;
}

bool ITTrackReco::iteratMaxMinWireCross(closClinRadLay &radCls_1v, closClinRadLay &radCls_2v, closClinRadLay::iterator &radCls_1v_it, closClstcol::iterator &closCls_it,  closClinRadLay::iterator &radCls_2v_it, closClstcol::iterator &lowRadClosCls_it, int matchZdir, float prevZVal, CellGeometryHandle *itwp, points3D &potLoopPoints) {
        closClinRadLay::iterator newRadCls_1v_it = radCls_2v_it;
        int crossCellId1 = (*lowRadClosCls_it)->_maxCellID;
        hitsInClsID::const_iterator newClHitsIDs_it=((*lowRadClosCls_it)->_hitIdx.end()-1);
        closClstcol::iterator newClosCls_it = lowRadClosCls_it;
        int tmpAbsRad = (newRadCls_1v_it->first-1);
        closClinRadLay::iterator newRadCls_2v_it = radCls_1v.find( tmpAbsRad );
        //float oldZ = potLoopPoints.back()._pos.z();

        //cout<<"iteratMaxMinWireCross"<<endl;
        bool newFoundWrCross = false;
        if (newRadCls_2v_it!=radCls_1v.end()) {
                for ( closClstcol::iterator newLowRadClosCls_it = newRadCls_2v_it->second.begin(); newLowRadClosCls_it != newRadCls_2v_it->second.end(); ++newLowRadClosCls_it ) {
                        int crossCellId2=(*newLowRadClosCls_it)->_minCellID;
                        hitsInClsID::const_iterator newLowRadClHitsIDs_it = (*newLowRadClosCls_it)->_hitIdx.begin();
                        if ( /*(notAssociatedHits[newClHitsIDs_it->second] || notAssociatedHits[newLowRadClHitsIDs_it->second]) &&*/
                                        hitCrossingChecked.find(std::pair<size_t,size_t>(newClHitsIDs_it->second,newLowRadClHitsIDs_it->second))==hitCrossingChecked.end() &&
                                        hitCrossingCheckedT.find(std::pair<size_t,size_t>(newLowRadClHitsIDs_it->second,newClHitsIDs_it->second))==hitCrossingCheckedT.end() ) {
                                //cout<<"checking"<<endl;
                                newFoundWrCross = findWireCross( crossCellId1, crossCellId2, newRadCls_1v_it, newClosCls_it, newClHitsIDs_it, newRadCls_2v_it, newLowRadClosCls_it, newLowRadClHitsIDs_it, matchZdir, prevZVal, itwp, potLoopPoints);
                                if (newFoundWrCross) {
                                        float newPrevZ = potLoopPoints.back()._pos.z();
                                        iteratMaxMinWireCross( radCls_2v, radCls_1v, newRadCls_1v_it, newClosCls_it, newRadCls_2v_it, newLowRadClosCls_it, matchZdir, newPrevZ, itwp, potLoopPoints );
                                }
                        }
                }
        }
        return newFoundWrCross;

}

bool ITTrackReco::iteratMinMaxWireCross(closClinRadLay &radCls_1v, closClinRadLay &radCls_2v, closClinRadLay::iterator &radCls_1v_it, closClstcol::iterator &closCls_it,  closClinRadLay::iterator &radCls_2v_it, closClstcol::iterator &lowRadClosCls_it, int matchZdir, float prevZVal, CellGeometryHandle *itwp, points3D &potLoopPoints) {
        closClinRadLay::iterator newRadCls_1v_it = radCls_2v_it;
        int crossCellId1 = (*lowRadClosCls_it)->_minCellID;
        hitsInClsID::const_iterator newClHitsIDs_it=(*lowRadClosCls_it)->_hitIdx.begin();
        closClstcol::iterator newClosCls_it = lowRadClosCls_it;
        int tmpAbsRad = (newRadCls_1v_it->first-1);
        closClinRadLay::iterator newRadCls_2v_it = radCls_1v.find( tmpAbsRad );

        //cout<<"iteratMinMaxWireCross"<<endl;
        bool newFoundWrCross = false;
        if (newRadCls_2v_it!=radCls_1v.end()) {
                for ( closClstcol::iterator newLowRadClosCls_it = newRadCls_2v_it->second.begin(); newLowRadClosCls_it != newRadCls_2v_it->second.end(); ++newLowRadClosCls_it ) {
                        int crossCellId2=(*newLowRadClosCls_it)->_maxCellID;
                        hitsInClsID::const_iterator newLowRadClHitsIDs_it = ((*newLowRadClosCls_it)->_hitIdx.end()-1);
                        if ( /*(notAssociatedHits[newClHitsIDs_it->second] || notAssociatedHits[newLowRadClHitsIDs_it->second]) &&*/
                                        hitCrossingChecked.find(std::pair<size_t,size_t>(newClHitsIDs_it->second,newLowRadClHitsIDs_it->second))==hitCrossingChecked.end() &&
                                        hitCrossingCheckedT.find(std::pair<size_t,size_t>(newLowRadClHitsIDs_it->second,newClHitsIDs_it->second))==hitCrossingCheckedT.end() ) {
                                //cout<<"checking"<<endl;
                                newFoundWrCross = findWireCross( crossCellId1, crossCellId2, newRadCls_1v_it, newClosCls_it, newClHitsIDs_it, newRadCls_2v_it, newLowRadClosCls_it, newLowRadClHitsIDs_it, matchZdir, prevZVal, itwp, potLoopPoints);
                                if (newFoundWrCross) {
                                        float newPrevZ = potLoopPoints.back()._pos.z();
                                        iteratMinMaxWireCross( radCls_2v, radCls_1v, newRadCls_1v_it, newClosCls_it, newRadCls_2v_it, newLowRadClosCls_it, matchZdir, newPrevZ, itwp, potLoopPoints );
                                }
                        }
                }
        }
        return newFoundWrCross;
}

int   ITTrackReco::findZclusters( std::vector<XYZP> &xyzp, float zdim, float zPeakSigma, float zPeakThr, TH1F *zHist ) {

        int nBins;
        TH1F *tmpZhist;
        if (zHist!=0x0) {
                tmpZhist = zHist;
                nBins = tmpZhist->GetNbinsX();
        } else {
                nBins = TMath::Nint( fabs(_maxZ-_minZ)/zdim );
                tmpZhist = new TH1F("tmpZhist","tmpZhist", nBins, _minZ,_maxZ);
        }
        //cout<<"ZMin "<<_minZ<<" Zmax "<<_maxZ<<" nBins "<<nBins<<endl;
        std::multimap <int, size_t> hitBinsrel;
        size_t kHit=0;
        for (std::vector<XYZP>::iterator xyzp_it=xyzp.begin(); xyzp_it!=xyzp.end(); ++xyzp_it) {
                hitBinsrel.insert ( std::pair<int,size_t> (tmpZhist->Fill(xyzp_it->_pos.z()), kHit) );
                ++kHit;
        }
        tmpZhist->Smooth(4);
        TSpectrum peaksSearcher;
        int nfound =  peaksSearcher.Search( tmpZhist, zPeakSigma, "", zPeakThr );
        cout<<"n of z Peaks found "<<nfound<<endl;
        if (nfound>0) {
                float *PosX = new float[nfound];
                float *PosX1 = new float[nfound];
                float *PosY = new float[nfound];
                //PosX = _peakFinder->GetPositionX();
                //PosX = peaksSearcher.GetPositionX();
                for (int iPF = 0; iPF < nfound; iPF++) {
                        PosX[iPF]   = peaksSearcher.GetPositionX()[iPF];
                        //                                    PosX[iPF]   = validPeakPos[iPF];
                        ///FixPos[iPF] = false;
                        //FixAmp[iPF] = false;
                        PosX1[iPF]  = (PosX[iPF]-_minZ)/zdim;
                        PosY[iPF]   = hitBinsrel.count( (int)(PosX1[iPF]) );

                        cout<<"peak pos at z "<<PosX[iPF]<<" = bin "<<(int)(PosX1[iPF])<<" with nhit "<<PosY[iPF]<<endl;
                }
                delete [] PosX;
                delete [] PosX1;
                delete [] PosY;
        }
        return nfound;
}

int   ITTrackReco::findZclusters( HelixTraj &trkhelixe, std::vector<XYZP> &xyzp, std::vector< std::vector<size_t> > &pntIdsLoops/*, std::vector<TrkHelix> &helixeFitLoops, std::vector<HelixTraj> &trkhelixeLoops*/) {

        size_t lHit=0;
        listZHitIdrels pntsOdered;
        for (std::vector<XYZP>::iterator xyzp_it = xyzp.begin(); xyzp_it != xyzp.end(); ++xyzp_it) {
                pntsOdered.push_back( make_pair(xyzp_it->_pos.z(),lHit) );
                ++lHit;
        }
        std::sort( pntsOdered.begin(), pntsOdered.end() );

        std::vector< std::pair<double, double> > loopsZRange;
        double firstZ = pntsOdered.front().first;
        double lastZ  = pntsOdered.back().first;
        HepPoint fPos=trkhelixe.position( trkhelixe.zFlight(firstZ) );
        double startRange=0.0;
        bool   prevNegSign=false;
        if ( (std::hypot(fPos.x(),fPos.y()) - _rIn) >= 0 ) {
                startRange = firstZ;
                prevNegSign = false;
        } else {
                startRange = -1000000.0;
                prevNegSign = true;
        }
        double tmpZ = firstZ+50.0;
        while ( tmpZ<=lastZ ) {
                HepPoint tmpPos=trkhelixe.position( trkhelixe.zFlight(tmpZ) );
                if ( (std::hypot(tmpPos.x(),tmpPos.y()) - _rIn) > 0 ) {
                        if ( prevNegSign ) {
                                startRange = tmpZ - 50.0;
                        }
                        prevNegSign = false;
                } else {
                        if ( !prevNegSign ) {
                                loopsZRange.push_back( make_pair(startRange,tmpZ) );
                        }
                        prevNegSign = true;
                }
                tmpZ+=50.0;
        }
        if (!prevNegSign) {
                loopsZRange.push_back( make_pair(startRange,lastZ) );
        }

        listZHitIdrels::iterator pntsOdered_it = pntsOdered.begin();

        for ( std::vector< std::pair<double, double> >::iterator loopsZRange_it = loopsZRange.begin(); loopsZRange_it != loopsZRange.end(); ++loopsZRange_it) {
                std::vector<size_t> iHelLoopPntIds;
                //TrkHelix iHelLoop;
                //HelixTraj iHelTraj;
                for ( ; pntsOdered_it != pntsOdered.end(); ++pntsOdered_it) {
                        if (pntsOdered_it->first>loopsZRange_it->second) {
                                break;
                        }
                        if (pntsOdered_it->first>=loopsZRange_it->first) {
                                iHelLoopPntIds.push_back(pntsOdered_it->second);
                        }
                }
                pntIdsLoops.push_back(iHelLoopPntIds);
                //helixeFitLoops.push_back(iHelLoop);
                //trkhelixeLoops.push_back(iHelTraj);
        }

        return loopsZRange.size();
}


bool ITTrackReco::findCircles(closClinRadLay &radCls, CellGeometryHandle *itwp, circlesCont & potcircles, std::string type) {
        int nPotBending=0;
        bool foundCircles = false;

        float refX, refY, refRes2, closeCntCellX, closeCntCellY/*, cmapU, cmapV, TmpR2, invTmpR2*/;
        int closeCntCellID, closeCntCellShift, halfwayCellID;
        //float cmapU_clsCntCell, cmapV_clsCntCell, CHTtanTheta0_DX, CHTtanTheta0_DY, CHTtanTheta0, CHTTheta0, CHTr0, CHTcosArcTanT0;
        //unsigned int votArrBin, halfwayVotArrBin;
        CHTVotArr_HitPtrrel votArr_HitPtrrel;
        ptrHitInClosCust refHitptr, clsCntCellHitptr;
        closeCntCellShift = (_minClusMultBend+1)/4+1;
        halfwayCellID = closeCntCellShift/2;
        if (halfwayCellID==0) halfwayCellID=1;
        if (halfwayCellID==closeCntCellShift) halfwayCellID=-1;

        float cellAveRes;
        float tmpCirChi2/*, tmpRad*/;
        bool firstInsert = false;

        //loop over all one view hits clusters, starting by the external radius, and make a CHT respect to the center of the starting cluster
        for ( closClinRadLay::iterator radCls_it = radCls.begin(); radCls_it != radCls.end(); ++radCls_it ) {
                for ( closClstcol::iterator closCls_it = radCls_it->second.begin(); closCls_it != radCls_it->second.end(); ++closCls_it ) {
                        if ((*closCls_it)->_nHit>((unsigned int)_minClusMultBend)){
                                ++nPotBending;

                                SimpleCircle2D tmpCircle;

                                itwp->SelectCell(radCls_it->first, (*closCls_it)->_centerCellID, _hitsInUpstream);
                                refX = itwp->GetCellCenter().getX();
                                refY = itwp->GetCellCenter().getY();
                                refRes2 = itwp->GetCellRad()*itwp->GetCellRad()/3.0;
                                cellAveRes = sqrt(refRes2);
                                tmpCircle._krmCircFit.addHit(refX,refY,cellAveRes,cellAveRes);

                                closeCntCellID = (*closCls_it)->_centerCellID+closeCntCellShift;//-1;
                                nCellPerLayer = itwp->GetITLayer()->nCells();
                                if (closeCntCellID<0) closeCntCellID+=nCellPerLayer;
                                closeCntCellID=closeCntCellID%nCellPerLayer;
                                itwp->SelectCell(radCls_it->first, closeCntCellID, _hitsInUpstream);
                                closeCntCellX = itwp->GetCellCenter().getX();
                                closeCntCellY = itwp->GetCellCenter().getY();

                                tmpCircle._krmCircFit.addHit(closeCntCellX,closeCntCellY,cellAveRes,cellAveRes);
                                //cout<<"------ ref (ID,x,y) "<<(*closCls_it)->_centerCellID<<", "<<refX<<", "<<refY<<" close "<<closeCntCellID<<", "<<closeCntCellX<<", "<<closeCntCellY<<endl;


                                halfwayCellID = (*closCls_it)->_centerCellID + halfwayCellID;
                                if (halfwayCellID<0) halfwayCellID+=nCellPerLayer;
                                halfwayCellID=halfwayCellID%nCellPerLayer;

                                firstInsert = false;

                                for (hitsInClsID::const_iterator clHitsIDs_it=(*closCls_it)->_hitIdx.begin(); clHitsIDs_it!=(*closCls_it)->_hitIdx.end(); ++clHitsIDs_it) {
                                        if ((*closCls_it)->_centerCellID==clHitsIDs_it->first) {
                                                tmpCircle._listHitptrs.insert( tmpCircle._listHitptrs.begin(), ptrHitInClosCust(clHitsIDs_it,closCls_it,radCls_it) );
                                                firstInsert = true;
                                                continue;
                                        }
                                        if (closeCntCellID==clHitsIDs_it->first) {
                                                if (firstInsert) {
                                                        tmpCircle._listHitptrs.insert( (tmpCircle._listHitptrs.begin()+1), ptrHitInClosCust(clHitsIDs_it,closCls_it,radCls_it) );
                                                }
                                                else {
                                                        tmpCircle._listHitptrs.insert( tmpCircle._listHitptrs.begin(), ptrHitInClosCust(clHitsIDs_it,closCls_it,radCls_it) );
                                                }
                                                continue;
                                        }
                                        if (halfwayCellID==clHitsIDs_it->first) {
                                        }

                                        itwp->SelectCell(itwp->GetSuperLayer(),itwp->GetCelRing(),clHitsIDs_it->first, _hitsInUpstream);
                                        tmpCircle._krmCircFit.addHit(itwp->GetCellCenter().getX(),itwp->GetCellCenter().getY(),cellAveRes,cellAveRes);
                                        tmpCircle._listHitptrs.push_back( ptrHitInClosCust(clHitsIDs_it,closCls_it,radCls_it) );

                                        //confMap.push_back( confMapPoint(cmapU,cmapV,invTmpR2,/*15.0*invSqrt12*/itwp->GetCellRad()*invSqrt3) );
                                        /*votArr_HitPtrrel.insert(
                                                          CHTVotArr_HitPtrrel::value_type(
                                                                          votArrBin,
                                                                          ptrHitInClosCust(clHitsIDs_it,closCls_it,radCls_it)
                                                                          )
                                          );*/
                                        //cout<<"Hit of the ref cluster goes into hit of CHT vot "<<votArrBin<<endl;
                                }
                                if ( !tmpCircle._krmCircFit.computeBestCirc() ) {
                                        cout<<"Circle not good"<<endl;
                                        break;////////////////FIXME
                                }

                                closClinRadLay::iterator sndRadCls_it = radCls_it;
                                ++sndRadCls_it;
                                for ( ; sndRadCls_it != radCls.end(); ++sndRadCls_it ) {
                                        itwp->SelectCell(sndRadCls_it->first, 0, _hitsInUpstream);
                                        cellAveRes = itwp->GetCellRad()*invSqrt3;

                                        for ( closClstcol::iterator sndClosCls_it = sndRadCls_it->second.begin(); sndClosCls_it != sndRadCls_it->second.end(); ++sndClosCls_it ) {
                                                for (hitsInClsID::const_iterator clHitsIDs_it=(*sndClosCls_it)->_hitIdx.begin(); clHitsIDs_it!=(*sndClosCls_it)->_hitIdx.end(); ++clHitsIDs_it) {
                                                        itwp->SelectCell(itwp->GetSuperLayer(),itwp->GetCelRing(),clHitsIDs_it->first, _hitsInUpstream);
                                                        if ( tmpCircle._krmCircFit.checkBfrAddHit(itwp->GetCellCenter().getX(),itwp->GetCellCenter().getY(),cellAveRes,cellAveRes) ) {
                                                                tmpCircle._listHitptrs.push_back( ptrHitInClosCust(clHitsIDs_it,sndClosCls_it,sndRadCls_it) );
                                                                cout<<"Point added"<<endl;
                                                        }
                                                }
                                        }
                                }

                                if ( tmpCircle._krmCircFit.ierror==0 ) {
                                        tmpCirChi2 = tmpCircle._krmCircFit.chicir/((float) (tmpCircle._krmCircFit.np()-3));
                                        if ( tmpCirChi2>0.001 && tmpCirChi2<_circChi2cut ) {
                                                cout<<"Circ Chi2 good"<<endl;
                                                //tmpRad = 1.0/tmpCircle._krmCircFit.rho;
                                                //if (tmpRad>=_minGoodR && tmpRad<=_maxGoodR ) {
                                                cout<<"Circle added"<<endl;
                                                potcircles.push_back( tmpCircle );
                                                foundCircles = true;
                                                //}
                                        }
                                }


                                /*
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
                                 */
                        }
                }
        }

        return foundCircles;
}

bool ITTrackReco::findCircles(closClinRadLay &radCls, CellGeometryHandle *itwp, std::vector<confMapDraw*> &rCMap, circlesCont & potcircles, std::string type, int minClusMultBend, unsigned int thrCHTaddVot, bool skipAssociated) {
        int nPotBending=0;
        bool foundCircles = false;

        float refX, refY, refRes2, closeCntCellX, closeCntCellY, cmapU, cmapV, TmpR2, invTmpR2;
        int closeCntCellID, closeCntCellShift, halfwayCellID;
        float cmapU_clsCntCell, cmapV_clsCntCell, CHTtanTheta0_DX, CHTtanTheta0_DY, CHTtanTheta0, CHTTheta0, CHTr0/*, CHTcosArcTanT0*/;
        unsigned int votArrBin, halfwayVotArrBin;
        CHTVotArr_HitPtrrel votArr_HitPtrrel;
        ptrHitInClosCust refHitptr, clsCntCellHitptr;
        closeCntCellShift = (minClusMultBend+1)/4+1;
        halfwayCellID = closeCntCellShift/2;
        if (halfwayCellID==0) halfwayCellID=1;
        if (halfwayCellID==closeCntCellShift) halfwayCellID=-1;

        float tmpCirChi2;
        bool startPointNotUsed;

        //loop over all one view hits clusters, starting by the external radius, and make a CHT respect to the center of the starting cluster
        for ( closClinRadLay::iterator radCls_it = radCls.begin(); radCls_it != radCls.end(); ++radCls_it ) {
                for ( closClstcol::iterator closCls_it = radCls_it->second.begin(); closCls_it != radCls_it->second.end(); ++closCls_it ) {
                        startPointNotUsed = true;
                        if ((*closCls_it)->_nHit>((unsigned int)minClusMultBend) ){
                                if (skipAssociated) {
                                        for (hitsInClsID::const_iterator clHitsIDs_it=(*closCls_it)->_hitIdx.begin(); clHitsIDs_it!=(*closCls_it)->_hitIdx.end(); ++clHitsIDs_it) {
                                                if ((*closCls_it)->_centerCellID==clHitsIDs_it->first) {
                                                        if ( !notAssociatedHits[clHitsIDs_it->second] ) {
                                                                startPointNotUsed = false;
                                                                break;
                                                        }
                                                }
                                        }
                                }
                                if (startPointNotUsed) {
                                        //continue;
                                        //}

                                        ++nPotBending;
                                        std::vector<confMapPoint > confMap;
                                        std::vector<std::pair<double, double> > CHTconfMap;  // Hough Transform value (theta, r) for the point of the Conformal Mapping
                                        CHTVotArr cmapCHTVotArr;
                                        cmapCHTVotArr.setTHR(thrCHTaddVot);

                                        itwp->SelectCell(radCls_it->first, (*closCls_it)->_centerCellID, _hitsInUpstream);
                                        refX = itwp->GetCellCenter().getX();
                                        refY = itwp->GetCellCenter().getY();
                                        refRes2 = itwp->GetCellRad()*itwp->GetCellRad()/3.0;
                                        closeCntCellID = (*closCls_it)->_centerCellID+closeCntCellShift;//-1;
                                        nCellPerLayer = itwp->GetITLayer()->nCells();
                                        if (closeCntCellID<0) closeCntCellID+=nCellPerLayer;
                                        closeCntCellID=closeCntCellID%nCellPerLayer;
                                        itwp->SelectCell(radCls_it->first, closeCntCellID, _hitsInUpstream);
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
                                                if ( skipAssociated && !notAssociatedHits[clHitsIDs_it->second] ) {
                                                        continue;
                                                }

                                                if ((*closCls_it)->_centerCellID==clHitsIDs_it->first) {
                                                        refHitptr.setHitInClosCust(clHitsIDs_it,closCls_it,radCls_it);
                                                        continue;
                                                }
                                                itwp->SelectCell(itwp->GetSuperLayer(),itwp->GetCelRing(),clHitsIDs_it->first, _hitsInUpstream);
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
                                                itwp->SelectCell(sndRadCls_it->first, 0, _hitsInUpstream);
                                                for ( closClstcol::iterator sndClosCls_it = sndRadCls_it->second.begin(); sndClosCls_it != sndRadCls_it->second.end(); ++sndClosCls_it ) {
                                                        for (hitsInClsID::const_iterator clHitsIDs_it=(*sndClosCls_it)->_hitIdx.begin(); clHitsIDs_it!=(*sndClosCls_it)->_hitIdx.end(); ++clHitsIDs_it) {
                                                                if ( skipAssociated && !notAssociatedHits[clHitsIDs_it->second] ) {
                                                                        continue;
                                                                }

                                                                itwp->SelectCell(itwp->GetSuperLayer(),itwp->GetCelRing(),clHitsIDs_it->first, _hitsInUpstream);
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

                                        cout<<"CHT vot hit over thr "<<cmapCHTVotArr._overTHRs.size()<<endl;
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

                                        cout<<"************************ 1 ************************"<<endl;

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
                                                if (tmpRad>=_minGoodR && tmpRad<=_maxGoodR ) {
                                                        tmpCircle.SetCenter( tmpCMap_A + refX, tmpCMap_B + refY, sqrt( tmpCMap_sigmaA2 + refRes2 ), sqrt( tmpCMap_sigmaB2 + refRes2 ) );
                                                        tmpCircle.SetRadius( tmpRad, sqrt( tmpCMap_A2*tmpCMap_sigmaA2 + tmpCMap_B2*tmpCMap_sigmaB2 )/tmpRad );
                                                        //for ( std::set<unsigned int>::iterator votInClus_it = closClustCHTVots_it->_listCHTVotIDs.begin();  votInClus_it != closClustCHTVots_it->_listCHTVotIDs.end(); ++votInClus_it ) {
                                                        //        cmapCHTVotArr.uIDToxy(*votInClus_it, tmpCHTvot_X, tmpCHTvot_Y);
                                                        //        cout<<" ("<<tmpCHTvot_X<<", "<<tmpCHTvot_Y<<") with val "<<cmapCHTVotArr.get(tmpCHTvot_X, tmpCHTvot_Y);
                                                        //}
                                                        cellAveRes = sqrt(refRes2);
                                                        if ( skipAssociated ) {
                                                                if ( notAssociatedHits[refHitptr.getinEventHitID()] ) {
                                                                        tmpCircle._listHitptrs.push_back(refHitptr);//insert(refHitptr);
                                                                        tmpCircle._krmCircFit.addHit( refX, refY, cellAveRes, cellAveRes );
                                                                }
                                                                if ( notAssociatedHits[clsCntCellHitptr.getinEventHitID()] ) {
                                                                        tmpCircle._listHitptrs.push_back(clsCntCellHitptr);//insert(clsCntCellHitptr);
                                                                        tmpCircle._krmCircFit.addHit( closeCntCellX, closeCntCellY, cellAveRes, cellAveRes );
                                                                }
                                                        } else {
                                                                tmpCircle._listHitptrs.push_back(refHitptr);//insert(refHitptr);
                                                                tmpCircle._krmCircFit.addHit( refX, refY, cellAveRes, cellAveRes );
                                                                tmpCircle._listHitptrs.push_back(clsCntCellHitptr);//insert(clsCntCellHitptr);
                                                                tmpCircle._krmCircFit.addHit( closeCntCellX, closeCntCellY, cellAveRes, cellAveRes );
                                                        }

                                                        //cout<<"For Potential circles"<<endl;
                                                        for (std::set< unsigned int >::const_iterator listCHTVotIDs_it=closClustCHTVots_it->_listCHTVotIDs.begin();
                                                                        listCHTVotIDs_it!=closClustCHTVots_it->_listCHTVotIDs.end(); ++listCHTVotIDs_it) {
                                                                //cout<<"bins in clusters "<<*listCHTVotIDs_it<<endl;
                                                                pair<CHTVotArr_HitPtrrel::iterator,CHTVotArr_HitPtrrel::iterator>
                                                                votArr_HitPtrrel_its = votArr_HitPtrrel.equal_range(*listCHTVotIDs_it);
                                                                for ( CHTVotArr_HitPtrrel::iterator tmpCHTVotArr_HitPtrrel_it=votArr_HitPtrrel_its.first;
                                                                                tmpCHTVotArr_HitPtrrel_it!=votArr_HitPtrrel_its.second;
                                                                                ++tmpCHTVotArr_HitPtrrel_it ) {
                                                                        if ( skipAssociated ) {
                                                                                if ( notAssociatedHits[tmpCHTVotArr_HitPtrrel_it->second.getinEventHitID()] ) {
                                                                                        tmpCircle._listHitptrs.push_back(tmpCHTVotArr_HitPtrrel_it->second);//insert(tmpCHTVotArr_HitPtrrel_it->second);
                                                                                        itwp->SelectCell(tmpCHTVotArr_HitPtrrel_it->second.getRadLayID(), tmpCHTVotArr_HitPtrrel_it->second.getinLayerCellID(), _hitsInUpstream);
                                                                                        cellAveRes = itwp->GetCellRad() * invSqrt3;
                                                                                        tmpCircle._krmCircFit.addHit( itwp->GetCellCenter().getX(), itwp->GetCellCenter().getY(), cellAveRes, cellAveRes );
                                                                                }
                                                                        } else {
                                                                                tmpCircle._listHitptrs.push_back(tmpCHTVotArr_HitPtrrel_it->second);//insert(tmpCHTVotArr_HitPtrrel_it->second);
                                                                                itwp->SelectCell(tmpCHTVotArr_HitPtrrel_it->second.getRadLayID(), tmpCHTVotArr_HitPtrrel_it->second.getinLayerCellID(), _hitsInUpstream);
                                                                                cellAveRes = itwp->GetCellRad() * invSqrt3;
                                                                                tmpCircle._krmCircFit.addHit( itwp->GetCellCenter().getX(), itwp->GetCellCenter().getY(), cellAveRes, cellAveRes );
                                                                        }
                                                                }
                                                        }

                                                        cout<<"************************ 2 ************************"<<endl;
                                                        if ( tmpCircle._krmCircFit.np()>3 && tmpCircle._krmCircFit.computeBestCirc(/*closeCntCellX, closeCntCellY*/) ) {
                                                                tmpCirChi2 = tmpCircle._krmCircFit.chicir/((float) (tmpCircle._krmCircFit.np()-3));
                                                                if ( tmpCirChi2>0.001 && tmpCirChi2<_circChi2cut ) {
                                                                        foundCircles = true;
                                                                        size_t irej = 2;
                                                                        std::vector<circPoint>::iterator points_end = tmpCircle._krmCircFit.points.end();
                                                                        //refine the circle by rejecting points
                                                                        for (std::vector<circPoint>::iterator points_it =(tmpCircle._krmCircFit.points.begin()+2); points_it <points_end;) {
                                                                                //std::vector<circPoint>::iterator test_points_it = points_it;
                                                                                irej = points_it-tmpCircle._krmCircFit.points.begin();
                                                                                if ( tmpCircle._krmCircFit.rejectHits(points_it) ) {
                                                                                        cout<<"Rejected hit "<<irej<<endl;
                                                                                        tmpCircle._listHitptrs.erase(tmpCircle._listHitptrs.begin()+irej);
                                                                                        points_it = (tmpCircle._krmCircFit.points.begin()+irej);
                                                                                        //--points_it;
                                                                                        //--irej;
                                                                                        points_end = tmpCircle._krmCircFit.points.end();
                                                                                        continue;
                                                                                }
                                                                                ++points_it;
                                                                                //++irej;
                                                                                if (tmpCircle._krmCircFit.points.size()<=3) {
                                                                                        foundCircles = false;
                                                                                        break;
                                                                                }
                                                                        }
                                                                        if (foundCircles) {
                                                                                tmpRad = 1.0/fabs(tmpCircle._krmCircFit.rho);
                                                                                if (tmpRad<_minGoodR || tmpRad>_maxGoodR) {
                                                                                        foundCircles = false;
                                                                                }
                                                                        }

                                                                        //try to recover points for the circle
                                                                        if (foundCircles) {
                                                                                if (skipAssociated) {
                                                                                        for (std::vector< ptrHitInClosCust >::iterator listHitptrs_it = tmpCircle._listHitptrs.begin(); listHitptrs_it != tmpCircle._listHitptrs.end(); ++listHitptrs_it) {
                                                                                                notAssociatedHits[listHitptrs_it->getinEventHitID()]=false;
                                                                                                ++nAssociatedHit;
                                                                                        }
                                                                                }

                                                                                closClinRadLay::iterator sndRadCls_it = radCls_it;
                                                                                ++sndRadCls_it;
                                                                                for ( ; sndRadCls_it != radCls.end(); ++sndRadCls_it ) {
                                                                                        itwp->SelectCell(sndRadCls_it->first, 0, _hitsInUpstream);
                                                                                        cellAveRes = itwp->GetCellRad()*invSqrt3;

                                                                                        for ( closClstcol::iterator sndClosCls_it = sndRadCls_it->second.begin(); sndClosCls_it != sndRadCls_it->second.end(); ++sndClosCls_it ) {
                                                                                                for (hitsInClsID::const_iterator clHitsIDs_it=(*sndClosCls_it)->_hitIdx.begin(); clHitsIDs_it!=(*sndClosCls_it)->_hitIdx.end(); ++clHitsIDs_it) {
                                                                                                        if ( skipAssociated ) {
                                                                                                                if ( notAssociatedHits[clHitsIDs_it->second] ) {
                                                                                                                        itwp->SelectCell(itwp->GetSuperLayer(),itwp->GetCelRing(),clHitsIDs_it->first, _hitsInUpstream);
                                                                                                                        if ( tmpCircle._krmCircFit.checkBfrAddHit(itwp->GetCellCenter().getX(),itwp->GetCellCenter().getY(),cellAveRes,cellAveRes,_circChi2cut, 1 ) ) {
                                                                                                                                tmpCircle._listHitptrs.push_back( ptrHitInClosCust(clHitsIDs_it,sndClosCls_it,sndRadCls_it) );
                                                                                                                                notAssociatedHits[clHitsIDs_it->second]=false;
                                                                                                                                ++nAssociatedHit;
                                                                                                                                cout<<"Point added"<<endl;
                                                                                                                        }
                                                                                                                }
                                                                                                        } else {
                                                                                                                itwp->SelectCell(itwp->GetSuperLayer(),itwp->GetCelRing(),clHitsIDs_it->first, _hitsInUpstream);
                                                                                                                if ( tmpCircle._krmCircFit.checkBfrAddHit(itwp->GetCellCenter().getX(),itwp->GetCellCenter().getY(),cellAveRes,cellAveRes,_circChi2cut, 1 ) ) {
                                                                                                                        tmpCircle._listHitptrs.push_back( ptrHitInClosCust(clHitsIDs_it,sndClosCls_it,sndRadCls_it) );
                                                                                                                        cout<<"Point added"<<endl;
                                                                                                                }
                                                                                                        }
                                                                                                }
                                                                                        }
                                                                                }
                                                                                potcircles.push_back( tmpCircle );
                                                                        }
                                                                }
                                                        }
                                                }
                                        }
                                        cout<<"************************ 3 ************************"<<endl;

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
                //int nAddedPoints = 0;
                for ( std::vector<circlePair >::iterator matchingCirc_it = matchingCirc.begin(); matchingCirc_it != matchingCirc.end(); ++matchingCirc_it ) {
                        //nAddedPoints = 0;
                        //if ( matchingCirc_it->first->_listHitptrs.size()>=matchingCirc_it->second->_listHitptrs.size() ) {
                        //nAddedPoints = matchingCirc_it->first->mergeCirc( *matchingCirc_it->second, removeCirc);
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

void ITTrackReco::markUsedHit( circlesCont &potCircles ) {
        for ( circlesCont::iterator potCircles_it = potCircles.begin(); potCircles_it != potCircles.end(); ++potCircles_it ) {
                for (std::vector< ptrHitInClosCust >::iterator listHitptrs_it = potCircles_it->_listHitptrs.begin(); listHitptrs_it != potCircles_it->_listHitptrs.end(); ++listHitptrs_it) {
                        notAssociatedHits[listHitptrs_it->getinEventHitID()]=false;
                        ++nAssociatedHit;
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


void ITTrackReco::searchTracks(closClinRadLay &radCls_Pl, closClinRadLay &radCls_Min, unique_ptr<TrackSeedCollection> &seeds,
                CellGeometryHandle *itwp, std::vector<TrkTimePeak> &tpeaks, size_t &ipeak, art::Handle<TrackerHitTimeClusterCollection> &tclustHandle, art::Handle<StrawHitCollection> &pdataHandle
                , VisibleGenElTrackCollection const* genEltrks
                ) {

        TrackerHitTimeCluster const&  tclust(tclustHandle.product()->at(ipeak));
        StrawHitCollection const* hits = pdataHandle.product();

        std::vector<confMapDraw*> rCMap_Pl, rCMap_Min, rCMap_Pl2, rCMap_Min2;

        circlesCont potCirclesMin;
        circlesCont potCirclesPl;

        std::vector<points3D> potLoops;
        int minClusMultBend = _minClusMultBend;
        size_t prevFound=0;
        int nPotLoops=0;
        while (minClusMultBend>=1) {
                //hitCrossingChecked.clear();
                //hitCrossingCheckedT.clear();
                findCircles3D(radCls_Min, radCls_Pl, itwp, potLoops,minClusMultBend);
                //hitCrossingChecked.clear();
                //hitCrossingCheckedT.clear();
                findCircles3D(radCls_Pl, radCls_Min, itwp, potLoops,minClusMultBend);
                /*potLoops.push_back(points3D());
                findCircles3D(radCls_Min, radCls_Pl, itwp, potLoops.back(),minClusMultBend);
                if (potLoops.back().size()>3) {
                        potLoops.push_back(points3D());
                } else {
                        potLoops.back().clear();
                }

                findCircles3D(radCls_Pl, radCls_Min, itwp, potLoops.back(),minClusMultBend);

                if (potLoops.back().size()<4) {
                        potLoops.pop_back();
                }*/
                /*
                for (std::vector<points3D>::iterator potLoops_it = (potLoops.begin()+prevFound); potLoops_it < potLoops.end(); ++potLoops_it) {
                        for (points3D::iterator loopPoints_it = potLoops_it->begin(); loopPoints_it != potLoops_it->end(); ++loopPoints_it) {
                                //cout<<"----------- "<<endl;
                                //cout<<loopPoints_it->getInEventHitID()<<endl;
                                //cout<<loopPoints_it->getCrossInEventHitID()<<endl;
                                notAssociatedHits[loopPoints_it->getInEventHitID()] = false;
                                notAssociatedHits[loopPoints_it->getCrossInEventHitID()] = false;
                        }
                }*/

                cout<<"Using cluster size cut for starting search of "<<minClusMultBend<<" found:"<<endl;
                cout<<"\t n potential loop "<<potLoops.size()<<endl;
                for (std::vector<points3D>::iterator potLoops_it = (potLoops.begin()+prevFound); potLoops_it < potLoops.end(); ++potLoops_it) {

                  cout<<"\t\t nCrossing Point for "<<++nPotLoops<<"-th = "<<potLoops_it->size()<<endl;
                  for (points3D::iterator loopPoints_it = potLoops_it->begin();
                       loopPoints_it != potLoops_it->end(); ++loopPoints_it) {
                    std::cout<<"3dpoint in Pot Loop "<<nPotLoops<<" hid "<<loopPoints_it->getInEventHitID()<<" xyz "
                             <<loopPoints_it->_pos.x()<<" "<<loopPoints_it->_pos.y()<<" "<<loopPoints_it->_pos.z()<<std::endl;
                  }
                }

                prevFound = potLoops.size();
                minClusMultBend/=2;
        }


        points3D selectedPnts;

        // track fitting objects for this peak
        std::vector<HelixFitResult> helixesfit;
        std::vector<HelixTraj> trkhelixes;
        size_t iTrackFound=0;
        std::map<size_t, std::vector<HelixFitResult> > helixesfitLoops;
        std::map<size_t, std::vector<HelixTraj> > trkhelixesLoops;
        std::map<size_t,  std::vector< std::vector<size_t> > > helixezPntIdsLoops;
        RbstFitPntsLpsPntsRel xyzpLoopRel;
        std::vector<XYZP> xyzp;
        std::vector<XYZP> xyzpSel;
        std::vector<std::pair<size_t,size_t> > sel_AllxyzPntsRel;
        std::vector<XYZP> xyzpRej;
        std::vector<std::pair<size_t,size_t> > rej_AllxyzPntsRel;

        //TH1F tmpZhist ("tmpZhist","tmpZhist", TMath::Nint( fabs(_maxZ-_minZ)/40.0 ), _minZ,_maxZ);

        if (potLoops.size()>0) {
                tpeaks.push_back( TrkTimePeak(tclust._meanTime,tclust._maxHitTime) );
                size_t iPotLoop = 0, iPntInLoop=0;
                for (std::vector<points3D>::iterator potLoops_it = potLoops.begin(); potLoops_it != potLoops.end(); ++potLoops_it) {
                        iPntInLoop=0;
                        for (points3D::iterator loopPoints_it = potLoops_it->begin(); loopPoints_it != potLoops_it->end(); ++loopPoints_it) {
                                itwp->SelectCell(loopPoints_it->getRadLayID(),loopPoints_it->getInLayerCellID(), _hitsInUpstream);
                                xyzp.push_back(XYZP( loopPoints_it->getInEventHitID(),
				loopPoints_it->_pos, itwp->GetCellDirection(), loopPoints_it->_sigmaz,loopPoints_it->_sigmax)); //sigmax=sigmay=sigmaR
                                //std::vector<XYZP>::iterator lastPoint = xyzp.end();
                                //--lastPoint;
                                //xyzpLoopRel.insert( RbstFitPntsLpsPntsRel::value_type ( lastPoint,
                                xyzpLoopRel.insert( RbstFitPntsLpsPntsRel::value_type ( xyzp.size()-1,
                                                std::pair<size_t, size_t > (iPotLoop,iPntInLoop) ) );
                                tpeaks.back()._trkptrs.push_back( hitIndex(loopPoints_it->getInEventHitID()) );
                                ++iPntInLoop;
                        }
                        ++iPotLoop;
                }

                //TrkKalFit seedfit, kalfit;
                // create track definitions for the different fits from this initial information
                TrkDef helixdef(hits,/*tpeaks[ipeak]*/tpeaks.back()._trkptrs);
                // track fitting objects for this peak
                HelixFitResult helixfit(helixdef);
                // robust helix fit
		XYZPVector xyzpv(xyzp);
                bool isfitted=_hfit.findHelix(xyzpv,helixfit);
                if(isfitted){
                        cout<<"Using Robust Fit, Circle Found with center "<<helixfit._center<<" and Rad "<<helixfit._radius<<endl;
                        for (RbstFitPntsLpsPntsRel::iterator xyzpLoopRel_it = xyzpLoopRel.begin(); xyzpLoopRel_it != xyzpLoopRel.end(); ++xyzpLoopRel_it ) {
                                XYZP &tmpXYZP = xyzp.at(xyzpLoopRel_it->first);
                                if ( tmpXYZP.use() ) {
                                        double tmpDist = fabs( sqrt( pow(tmpXYZP._pos.x()-helixfit._center.x(),2) +
                                                        pow(tmpXYZP._pos.y()-helixfit._center.y(),2) ) -
                                                        helixfit._radius );
                                        if ( tmpDist<3.0*tmpXYZP._rerr ) {
                                                selectedPnts.push_back( potLoops.at(xyzpLoopRel_it->second.first).at(xyzpLoopRel_it->second.second) );
                                                xyzpSel.push_back(tmpXYZP);
                                                sel_AllxyzPntsRel.push_back( std::pair<size_t,size_t> (xyzpSel.size()-1, xyzpLoopRel_it->first) );
                                        } else {
                                                tmpXYZP.setUse(false);
                                                xyzpRej.push_back(tmpXYZP);
                                                rej_AllxyzPntsRel.push_back( std::pair<size_t,size_t> (xyzpRej.size()-1, xyzpLoopRel_it->first) );
                                        }
                                }
                                xyzpRej.push_back(tmpXYZP);
                                rej_AllxyzPntsRel.push_back( std::pair<size_t,size_t> (xyzpRej.size()-1, xyzpLoopRel_it->first) );
                        }

                        if (helixfit._radius>=_minGoodR && helixfit._radius<=_maxGoodR ) {
                                helixesfit.push_back(helixfit);
                        } else {
                                cout<<"Track rejected for radius"<<endl;
                                isfitted = false;
                        }
                        //try to use the rejected points to find a new circle!!! FIXME
                }

                if(isfitted) {
                        HepVector hpar;
                        HepVector hparerr;
                        _hfit.helixParams(helixfit,hpar,hparerr);
                        HepSymMatrix hcov = vT_times_v(hparerr);
                        HelixTraj seed(hpar,hcov);
                        std::cout<<"fitted "<<isfitted<<" ";
                        seed.printAll(std::cout);
                        trkhelixes.push_back(seed);

                        findZclusters( trkhelixes.back(), xyzpSel, helixezPntIdsLoops[iTrackFound]/*, helixesfitLoops[iTrackFound], trkhelixesLoops[iTrackFound]*/ );
                        bool isLoopFitted;
                        for ( size_t iTrkLoop=0; iTrkLoop<helixezPntIdsLoops[iTrackFound].size(); ++iTrkLoop ) {
                                std::vector<size_t> &loopPntIds = helixezPntIdsLoops[iTrackFound].at(iTrkLoop);
                                std::vector<XYZP> xyzpLoop;
                                for (std::vector<size_t>::iterator loopPntIds_it = loopPntIds.begin(); loopPntIds_it != loopPntIds.end(); ++loopPntIds_it) {
                                        xyzpLoop.push_back( xyzpSel.at(*loopPntIds_it) );
                                }
                                HelixFitResult helixfitLoop(helixdef);
				XYZPVector xyzplv(xyzpLoop);
                                isLoopFitted=_hfit.findHelix(xyzplv,helixfitLoop);
                                if (isLoopFitted) {
                                        helixesfitLoops[iTrackFound].push_back(helixfitLoop);
                                        HepVector hparLoop;
                                        HepVector hparerrLoop;
                                        _hfit.helixParams(helixfitLoop,hparLoop,hparerrLoop);
                                        HepSymMatrix hcovLoop = vT_times_v(hparerrLoop);
                                        HelixTraj seedLoop(hparLoop,hcovLoop);
                                        std::cout<<"fitted "<<iTrkLoop<<"-th loop:"<<endl;
                                        seedLoop.printAll(std::cout);
                                        trkhelixesLoops[iTrackFound].push_back(seedLoop);
                                }

                        }

                        if (_diagLevel>1) {
                                cout<<"For "<<iTrackFound<<"-th track n Loops "<<helixezPntIdsLoops[iTrackFound].size()<<endl;
                                for ( size_t iTrkLoop=0; iTrkLoop<helixezPntIdsLoops[iTrackFound].size(); ++iTrkLoop ) {
                                        cout<<"\t for "<<iTrkLoop<<"-th loop z points:"<<endl;
                                        std::vector<size_t> &loopPntIds = helixezPntIdsLoops[iTrackFound].at(iTrkLoop);
                                        for (std::vector<size_t>::iterator loopPntIds_it = loopPntIds.begin(); loopPntIds_it != loopPntIds.end(); ++loopPntIds_it) {
                                                cout<<"\t\t "<<xyzpSel.at(*loopPntIds_it)._pos<<endl;
                                        }
                                }
                        }

                        ++iTrackFound;
                }
        }



        int jTrackFound = 0;
        for(std::vector<HelixTraj>::iterator ittrk=trkhelixes.begin(); ittrk!=trkhelixes.end(); ittrk++){
                cout<<"For i track "<<jTrackFound<<" there are "<<trkhelixesLoops[jTrackFound].size()<<" Loops to draw"<<endl;


                TrackSeed tmpOutSeed;
                tmpOutSeed._relatedTimeCluster=TrackerHitTimeClusterPtr(tclustHandle,ipeak);
                HelixTraj2HelixVal (*ittrk, tmpOutSeed._fullTrkSeed);

                std::set<size_t > hitInserted;
                double minSeedTime = 1.e9, tmpSeedHitTime=0.0;
                for ( size_t jSeedPoint = 0; jSeedPoint<sel_AllxyzPntsRel.size(); ++jSeedPoint ) {
                        std::pair<size_t, size_t> &iPLoop_jPnt = xyzpLoopRel[sel_AllxyzPntsRel.at(jSeedPoint).second];
                        corssingPoints &crossPoint = potLoops[iPLoop_jPnt.first][iPLoop_jPnt.second];
                        //insert the two crossing points
                        if ( hitInserted.find(crossPoint.getInEventHitID())==hitInserted.end() ) {
                                tmpOutSeed._fullTrkSeed._selectedTrackerHitsIdx.push_back( HitIndex( crossPoint.getInEventHitID() ) );
                                tmpOutSeed._selectedTrackerHits.push_back( StrawHitPtr ( pdataHandle, crossPoint.getInEventHitID() ) );
                                hitInserted.insert(crossPoint.getInEventHitID());
                                tmpSeedHitTime = tmpOutSeed._selectedTrackerHits.back()->time();// correct for TOF and signal propagation!!! FIXME
                                if (tmpSeedHitTime<minSeedTime) { minSeedTime=tmpSeedHitTime; }
                        }
                        if ( hitInserted.find(crossPoint.getCrossInEventHitID())==hitInserted.end() ) {
                                tmpOutSeed._fullTrkSeed._selectedTrackerHitsIdx.push_back( HitIndex( crossPoint.getCrossInEventHitID() ) );
                                tmpOutSeed._selectedTrackerHits.push_back( StrawHitPtr ( pdataHandle, crossPoint.getCrossInEventHitID() ) );
                                hitInserted.insert(crossPoint.getCrossInEventHitID());
                                tmpSeedHitTime = tmpOutSeed._selectedTrackerHits.back()->time();// correct for TOF and signal propagation!!! FIXME
                                if (tmpSeedHitTime<minSeedTime) { minSeedTime=tmpSeedHitTime; }
                        }

                        //inserted all the points of the cluster that contains the upper of the two crossing points
                        for ( hitsInClsID::const_iterator hitsInCls_it = crossPoint.getHitClust()._hitIdx.begin();
                                        hitsInCls_it != crossPoint.getHitClust()._hitIdx.end(); ++hitsInCls_it ) {
                                if ( hitInserted.find(hitsInCls_it->second)==hitInserted.end() ) {
                                        tmpOutSeed._fullTrkSeed._selectedTrackerHitsIdx.push_back( HitIndex( hitsInCls_it->second ) );
                                        tmpOutSeed._selectedTrackerHits.push_back( StrawHitPtr ( pdataHandle, hitsInCls_it->second ) );
                                        hitInserted.insert(hitsInCls_it->second);
                                        tmpSeedHitTime = tmpOutSeed._selectedTrackerHits.back()->time();// correct for TOF and signal propagation!!! FIXME
                                        if (tmpSeedHitTime<minSeedTime) { minSeedTime=tmpSeedHitTime; }
                                }
                        }
                }
                tmpOutSeed._t0    = minSeedTime;
                tmpOutSeed._errt0 = _tpeakerr;

                int kTrckLoop=0;
                for (std::vector<HelixTraj>::iterator trkhelixesLoop_it = trkhelixesLoops[jTrackFound].begin(); trkhelixesLoop_it != trkhelixesLoops[jTrackFound].end(); ++trkhelixesLoop_it) {
                        cout<<"Drawing single Loop for track "<<jTrackFound<<endl;
                        HelixVal tmpSeedLoop;
                        HelixTraj2HelixVal (*trkhelixesLoop_it, tmpSeedLoop);
                        std::vector<size_t> &hitsInLoopId = helixezPntIdsLoops[jTrackFound].at(kTrckLoop);
                        std::set<size_t > loopHitInserted;

                        for (std::vector<size_t>::iterator hitsInLoopId_it = hitsInLoopId.begin(); hitsInLoopId_it != hitsInLoopId.end(); ++hitsInLoopId_it) {
                                std::pair<size_t, size_t> &iPLoop_jPnt = xyzpLoopRel[sel_AllxyzPntsRel.at(*hitsInLoopId_it).second];
                                corssingPoints &crossPoint = potLoops[iPLoop_jPnt.first][iPLoop_jPnt.second];
                                //insert the two crossing points
                                if ( loopHitInserted.find(crossPoint.getInEventHitID())==loopHitInserted.end() ) {
                                        tmpSeedLoop._selectedTrackerHitsIdx.push_back( HitIndex( crossPoint.getInEventHitID() ) );
                                        tmpOutSeed._selectedTrackerHits.push_back( StrawHitPtr ( pdataHandle, crossPoint.getInEventHitID() ) );
                                        loopHitInserted.insert(crossPoint.getInEventHitID());
                                        tmpSeedHitTime = tmpOutSeed._selectedTrackerHits.back()->time();// correct for TOF and signal propagation!!! FIXME
                                        if (tmpSeedHitTime<minSeedTime) { minSeedTime=tmpSeedHitTime; }
                                }
                                if ( loopHitInserted.find(crossPoint.getCrossInEventHitID())==loopHitInserted.end() ) {
                                        tmpSeedLoop._selectedTrackerHitsIdx.push_back( HitIndex( crossPoint.getCrossInEventHitID() ) );
                                        tmpOutSeed._selectedTrackerHits.push_back( StrawHitPtr ( pdataHandle, crossPoint.getCrossInEventHitID() ) );
                                        loopHitInserted.insert(crossPoint.getCrossInEventHitID());
                                        tmpSeedHitTime = tmpOutSeed._selectedTrackerHits.back()->time();// correct for TOF and signal propagation!!! FIXME
                                        if (tmpSeedHitTime<minSeedTime) { minSeedTime=tmpSeedHitTime; }
                                }

                                //inserted all the points of the cluster that contains the upper of the two crossing points
                                for ( hitsInClsID::const_iterator hitsInCls_it = crossPoint.getHitClust()._hitIdx.begin();
                                                hitsInCls_it != crossPoint.getHitClust()._hitIdx.end(); ++hitsInCls_it ) {
                                        if ( loopHitInserted.find(hitsInCls_it->second)==loopHitInserted.end() ) {
                                                tmpSeedLoop._selectedTrackerHitsIdx.push_back( HitIndex( hitsInCls_it->second ) );
                                                tmpOutSeed._selectedTrackerHits.push_back( StrawHitPtr ( pdataHandle, hitsInCls_it->second ) );
                                                loopHitInserted.insert(hitsInCls_it->second);
                                                tmpSeedHitTime = tmpOutSeed._selectedTrackerHits.back()->time();// correct for TOF and signal propagation!!! FIXME
                                                if (tmpSeedHitTime<minSeedTime) { minSeedTime=tmpSeedHitTime; }
                                        }
                                }
                        }

                        tmpOutSeed._loopSeeds.push_back(tmpSeedLoop);
                        ++kTrckLoop;
                }

                seeds->push_back(tmpOutSeed);
                ++jTrackFound;
        }

        /*
        for (std::vector<XYZP>::iterator xyzp_it = xyzp.begin(); xyzp_it != xyzp.end(); ++xyzp_it ) {
                cout<<" xyzp "<<xyzp_it->_pos<<" serr "<<xyzp_it->_rerr<<" use "<<xyzp_it->use()<<endl;
        }
        cout<<"-----------------------------------------"<<endl;
        for (RbstFitPntsLpsPntsRel::iterator xyzpLoopRel_it = xyzpLoopRel.begin(); xyzpLoopRel_it != xyzpLoopRel.end(); ++xyzpLoopRel_it ) {
                //cout<<" xyzp "<<xyzpLoopRel_it->first->_pos<<" serr "<<xyzpLoopRel_it->first->_rerr<<" use "<<xyzpLoopRel_it->first->use()<<endl;
                cout<<" xyzp "<<xyzp.at(xyzpLoopRel_it->first)._pos<<" serr "<<xyzp.at(xyzpLoopRel_it->first)._rerr<<" use "<<xyzp.at(xyzpLoopRel_it->first).use()<<endl;
        }
        */




        //        cout<<"Using cluster size cut for starting search of "<<_minClusMultBend<<" found:"<<endl;
        //        cout<<"\t n potential loop "<<potLoops.size()<<endl;
        //        int nPotLoops=0;
        //        for (std::vector<points3D>::iterator potLoops_it = potLoops.begin(); potLoops_it != potLoops.end(); ++potLoops_it) {
        //                cout<<"\t\t nCrossing Point for "<<++nPotLoops<<"-th = "<<potLoops_it->size()<<endl;
        //        }

        //        cout<<"   ----  stereo - ----   "<<endl;
        //        findCircles3D(radCls_Min, radCls_Pl, itwp, potCirclesMin);
        //
        //        bool foundCircles_Min = false ;/*findCircles(radCls_Min,itwp,rCMap_Min,potCirclesMin,"Min", _minClusMultBend, 10, true);
        //        //bool foundCircles_Min = findCircles(radCls_Min,itwp,potCirclesMin,"Min");
        //        if ( foundCircles_Min ) {
        //                printPotCircles(potCirclesMin);
        //                //mergeCircles(potCirclesMin);
        //                //markUsedHit(potCirclesMin);
        //                //cout<<"**** After merging ****"<<endl;
        //                //printPotCircles(potCirclesMin);
        //        }
        //        cout<<"   ----  second search stereo - ----   "<<endl;
        //        foundCircles_Min = findCircles(radCls_Min,itwp,rCMap_Min2,potCirclesMin,"Min", 6, 5, true);
        //        //bool foundCircles_Min = findCircles(radCls_Min,itwp,potCirclesMin,"Min");
        //        if ( foundCircles_Min ) {
        //                printPotCircles(potCirclesMin);
        //                //mergeCircles(potCirclesMin);
        //                //printPotCircles(potCirclesMin);
        //        }
        //        */
        //        cout<<"   ----  stereo + ----   "<<endl;
        //        findCircles3D(radCls_Pl, radCls_Min, itwp, potCirclesPl);
        //        bool foundCircles_Pl = false; /*findCircles(radCls_Pl,itwp,rCMap_Pl,potCirclesPl,"Pl", _minClusMultBend, 10,true);
        //        //bool foundCircles_Pl = findCircles(radCls_Pl,itwp,potCirclesPl,"Pl");
        //        if ( foundCircles_Pl ) {
        //                printPotCircles(potCirclesPl);
        //                //mergeCircles(potCirclesPl);
        //                //markUsedHit(potCirclesPl);
        //                //cout<<"**** After merging ****"<<endl;
        //                //printPotCircles(potCirclesPl);
        //
        //        }
        //        cout<<"   ----  second search stereo + ----   "<<endl;
        //        foundCircles_Pl = findCircles(radCls_Pl,itwp,rCMap_Pl2,potCirclesPl,"Pl", 6, 5, true);
        //        //bool foundCircles_Pl = findCircles(radCls_Pl,itwp,potCirclesPl,"Pl");
        //        if ( foundCircles_Pl ) {
        //                printPotCircles(potCirclesPl);
        //                //mergeCircles(potCirclesPl);
        //                //printPotCircles(potCirclesPl);
        //
        //        }*/
        //        //findCircles(radCls_Min_Up,itwp,rCMap_Min_Up,"Min_Up");
        //        //findCircles(radCls_Pl_Up,itwp,rCMap_Pl_Up,"Pl_Up");
        //
        //        if ( foundCircles_Min && foundCircles_Pl) {
        //
        //        }

        if (_doDisplay) {

                double rIn  = _rIn/CLHEP::cm;
                double rOut = _rOut/CLHEP::cm;

                //cout<<"************************ 1.1 ************************"<<endl;

                TEllipse innerWall (0.0,0.0,rIn,rIn);
                innerWall.SetFillStyle(0);
                innerWall.SetLineWidth(1.5);

                TEllipse outerWall (0.0,0.0,rOut,rOut);
                outerWall.SetFillStyle(0);
                outerWall.SetLineWidth(1.5);

                TCanvas *tmpRCP_Min=0,*tmpRCP_Pl=0,*tmpCircleFound=0;
                TH2F hHitStereoPl( "hHitStereoPl",  "Stereo+ Hits per Event at z=0", 1500, -75.0, 75.0, 1500, -75.0, 75.0 );
                TH2F hHitStereoMin( "hHitStereoMin",  "Stereo- Hits per Event at z=0", 1500, -75.0, 75.0, 1500, -75.0, 75.0 );
                std::vector<TEllipse *> drawPotCirclesMin;
                std::vector<TEllipse *> drawPotKrmCirclesMin;
                std::vector<TClonesArray *> drawPotCircPntsMin;

                std::vector<TEllipse *> drawPotCirclesPl;
                std::vector<TEllipse *> drawPotKrmCirclesPl;
                std::vector<TClonesArray *> drawPotCircPntsPl;

                bool drawStereoProj=false;
                if(drawStereoProj){//stereo projections
                        tmpRCP_Min = new TCanvas();
                        /*
        int nlrgCl_Min = rCMap_Min.size();
        if (nlrgCl_Min>0) {
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
        }
                         */
                        //cout<<"************************ 1.2 ************************"<<endl;
                        tmpRCP_Pl = new TCanvas();
                        /*
        int nlrgCl_Pl = rCMap_Pl.size();
        if (nlrgCl_Pl>0) {
        int iPad = 0, maxPad = 1;
        if (nlrgCl_Pl>=1) {
                tmpRCP_Pl->Divide(2,nlrgCl_Pl);
                maxPad = 2*nlrgCl_Pl;
                iPad=1;
        }
        int iColor = 0;
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
        }
                         */

                        tmpRCP_Min->Modified();
                        tmpRCP_Min->Update();

                        tmpRCP_Pl->Modified();
                        tmpRCP_Pl->Update();


                        //cout<<"************************ 1.3 ************************"<<endl;

                        tmpCircleFound = new TCanvas("PotCirclesDraw","",860,430);
                        tmpCircleFound->Divide(2,1);
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
                        for (circlesCont::iterator potCirclesMin_it = potCirclesMin.begin(); potCirclesMin_it != potCirclesMin.end(); ++potCirclesMin_it ) {
                                drawPotCirclesMin.push_back( new TEllipse( potCirclesMin_it->_center.x()/CLHEP::cm, potCirclesMin_it->_center.y()/CLHEP::cm, potCirclesMin_it->_radius/CLHEP::cm, potCirclesMin_it->_radius/CLHEP::cm) );
                                drawPotCirclesMin.back()->SetFillStyle(0);
                                drawPotCirclesMin.back()->SetLineWidth(2);
                                drawPotCirclesMin.back()->SetLineColor(++iCol);
                                drawPotCirclesMin.back()->SetLineStyle(4);
                                //drawPotCirclesMin.back()->Draw("same");
                                CLHEP::Hep2Vector radDir;
                                float rSign = (potCirclesMin_it->_krmCircFit.rho>=0.0) ? 1.0 : -1.0;
                                radDir.setX( cos(potCirclesMin_it->_krmCircFit.phi) );
                                radDir.setY( sin(potCirclesMin_it->_krmCircFit.phi) );
                                radDir.rotate( /*( (potCirclesMin_it->_krmCircFit.rho>=0.0) ? -90.0 : 90.0 )*/-90.0*rSign*CLHEP::degree );
                                radDir=radDir.unit();
                                HepGeom::Point3D<double> CirCenter (0.0,0.0,0.0);
                                float krmRad = 1.0/fabs(potCirclesMin_it->_krmCircFit.rho);
                                CirCenter = CirCenter +  (krmRad + rSign*potCirclesMin_it->_krmCircFit.dca)*radDir;
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
                                        itwp->SelectCell( potCirclesPntsMin_it->getRadLayID(), potCirclesPntsMin_it->getinLayerCellID(), _hitsInUpstream );
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

                        //cout<<"************************ 1.4 ************************"<<endl;

                        tmpCircleFound->cd(2);
                        hHitStereoPl.SetStats(kFALSE);
                        hHitStereoPl.Draw();
                        innerWall.Draw("same");
                        outerWall.Draw("same");
                        iCol = 0;
                        for (circlesCont::iterator potCirclesPl_it = potCirclesPl.begin(); potCirclesPl_it != potCirclesPl.end(); ++potCirclesPl_it ) {
                                drawPotCirclesPl.push_back( new TEllipse( potCirclesPl_it->_center.x()/CLHEP::cm, potCirclesPl_it->_center.y()/CLHEP::cm, potCirclesPl_it->_radius/CLHEP::cm, potCirclesPl_it->_radius/CLHEP::cm) );
                                drawPotCirclesPl.back()->SetFillStyle(0);
                                drawPotCirclesPl.back()->SetLineWidth(2);
                                drawPotCirclesPl.back()->SetLineColor(++iCol);
                                drawPotCirclesPl.back()->SetLineStyle(4);
                                //drawPotCirclesPl.back()->Draw("same");
                                CLHEP::Hep2Vector radDir;
                                float rSign = (potCirclesPl_it->_krmCircFit.rho>=0.0) ? 1.0 : -1.0;
                                radDir.setX( cos(potCirclesPl_it->_krmCircFit.phi) );
                                radDir.setY( sin(potCirclesPl_it->_krmCircFit.phi) );
                                radDir.rotate( /*( (potCirclesPl_it->_krmCircFit.rho>=0.0) ? -90.0 : 90.0 )*/-90.0*rSign*CLHEP::degree );
                                radDir=radDir.unit();
                                HepGeom::Point3D<double> CirCenter (0.0,0.0,0.0);
                                float krmRad = 1.0/fabs(potCirclesPl_it->_krmCircFit.rho);
                                CirCenter = CirCenter + (krmRad + rSign*potCirclesPl_it->_krmCircFit.dca)*radDir;
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
                                        itwp->SelectCell( potCirclesPntsPl_it->getRadLayID(), potCirclesPntsPl_it->getinLayerCellID(), _hitsInUpstream );
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
                }

                TCanvas *tmpCircleFound3D = new TCanvas("PotCirclesDraw3D","",860,860);
                tmpCircleFound3D->Divide(2,2);
                TH2F hHit3DT( "hHit3DT",  "Hits per Event at z crossing", 1500, -75.0, 75.0, 1500, -75.0, 75.0 );
                TH2F hHit3DL( "hHit3DL",  "Hits per Event at z", TMath::Nint((_maxZ-_minZ)/10.0), _minZ/CLHEP::cm, _maxZ/CLHEP::cm, 750, 0, 75.0 );
                hHit3DT.SetXTitle("X [cm]");
                hHit3DT.SetYTitle("Y [cm]");
                hHit3DL.SetXTitle("Z [cm]");
                hHit3DL.SetYTitle("R [cm]");
                tmpCircleFound3D->cd(1);
                hHit3DT.SetStats(kFALSE);
                hHit3DT.Draw();
                innerWall.Draw("same");
                outerWall.Draw("same");
                tmpCircleFound3D->cd(2);
                TLine innerWallL;
                TLine outerWallL;
                innerWallL.SetX1(hHit3DL.GetXaxis()->GetXmin());
                innerWallL.SetY1(rIn);
                innerWallL.SetX2(hHit3DL.GetXaxis()->GetXmax());
                innerWallL.SetY2(rIn);
                outerWallL.SetX1(hHit3DL.GetXaxis()->GetXmin());
                outerWallL.SetY1(rOut);
                outerWallL.SetX2(hHit3DL.GetXaxis()->GetXmax());
                outerWallL.SetY2(rOut);
                innerWallL.SetLineWidth(1.5);
                outerWallL.SetLineWidth(1.5);
                hHit3DL.SetStats(kFALSE);
                hHit3DL.Draw();
                innerWallL.Draw("same");
                outerWallL.Draw("same");
                tmpCircleFound3D->cd(3);
                hHit3DT.Draw();
                innerWall.Draw("same");
                outerWall.Draw("same");
                tmpCircleFound3D->cd(4);
                hHit3DL.Draw();
                innerWallL.Draw("same");
                outerWallL.Draw("same");
                tmpCircleFound3D->cd(3);

                int iCol = 0;
                std::vector<TClonesArray *> drawPotCircPnts3D;
                std::vector<TGraphErrors *> drawPotCircPnts3D_L;
                std::vector<float *> zarr;
                std::vector<float *> rarr;
                std::vector<float *> errzarr;
                std::vector<float *> errrarr;
                for (std::vector<points3D>::iterator potLoops_it = potLoops.begin(); potLoops_it != potLoops.end(); ++potLoops_it) {
                        ++iCol;
                        zarr.push_back(new float[potLoops_it->size()]);
                        rarr.push_back(new float[potLoops_it->size()]);
                        errzarr.push_back(new float[potLoops_it->size()]);
                        errrarr.push_back(new float[potLoops_it->size()]);
                        drawPotCircPnts3D.push_back( new TClonesArray("TEllipse") );
                        drawPotCircPnts3D.back()->ExpandCreateFast(potLoops_it->size());
                        int iHit = 0;
                        tmpCircleFound3D->cd(1);
                        for (points3D::iterator loopPoints_it = potLoops_it->begin(); loopPoints_it != potLoops_it->end(); ++loopPoints_it) {
                                //cout<<"hit "<<loopPoints_it->_pos.x()<<endl;
                                ((TEllipse *) drawPotCircPnts3D.back()->At(iHit))->SetX1(loopPoints_it->_pos.x()/CLHEP::cm);
                                ((TEllipse *) drawPotCircPnts3D.back()->At(iHit))->SetY1(loopPoints_it->_pos.y()/CLHEP::cm);
                                ((TEllipse *) drawPotCircPnts3D.back()->At(iHit))->SetR1(loopPoints_it->_sigmax/invSqrt3/CLHEP::cm);
                                ((TEllipse *) drawPotCircPnts3D.back()->At(iHit))->SetR2(loopPoints_it->_sigmay/invSqrt3/CLHEP::cm );
                                ((TEllipse *) drawPotCircPnts3D.back()->At(iHit))->SetTheta( loopPoints_it->_pos.phi()/CLHEP::degree );
                                ((TEllipse *) drawPotCircPnts3D.back()->At(iHit))->SetFillStyle(0);
                                ((TEllipse *) drawPotCircPnts3D.back()->At(iHit))->SetLineWidth(2);
                                ((TEllipse *) drawPotCircPnts3D.back()->At(iHit))->SetLineColor(iCol);
                                ((TEllipse *) drawPotCircPnts3D.back()->At(iHit))->Draw("same");
                                zarr.back()[iHit] = loopPoints_it->_pos.z()/CLHEP::cm;
                                errzarr.back()[iHit] = loopPoints_it->_sigmaz/CLHEP::cm;
                                rarr.back()[iHit] = loopPoints_it->_pos.perp()/CLHEP::cm;
                                errrarr.back()[iHit] = loopPoints_it->_sigmax/CLHEP::cm;
                                ++iHit;
                        }
                        tmpCircleFound3D->cd(2);
                        drawPotCircPnts3D_L.push_back( new TGraphErrors(potLoops_it->size(),zarr.back(),rarr.back(),errzarr.back(),errrarr.back()) );
                        drawPotCircPnts3D_L.back()->SetMarkerColor(iCol);
                        drawPotCircPnts3D_L.back()->SetMarkerSize(0.8);
                        drawPotCircPnts3D_L.back()->SetMarkerStyle(8);
                        drawPotCircPnts3D_L.back()->Draw("P");
                }
                for (size_t jLoop = 0; jLoop<zarr.size();++jLoop) {
                        delete [] zarr.at(jLoop);
                        delete [] rarr.at(jLoop);
                        delete [] errzarr.at(jLoop);
                        delete [] errrarr.at(jLoop);
                }


                tmpCircleFound3D->cd(1);
                int iRobCirCol = 0;
                std::vector<TEllipse *> drawPotCirclesRobustFit;
                for (std::vector<HelixFitResult>::iterator helixesfit_it = helixesfit.begin(); helixesfit_it != helixesfit.end(); ++helixesfit_it ) {
                        drawPotCirclesRobustFit.push_back( new TEllipse( helixesfit_it->_center.x()/CLHEP::cm, helixesfit_it->_center.y()/CLHEP::cm, helixesfit_it->_radius/CLHEP::cm, helixesfit_it->_radius/CLHEP::cm) );
                        drawPotCirclesRobustFit.back()->SetFillStyle(0);
                        drawPotCirclesRobustFit.back()->SetLineWidth(2);
                        drawPotCirclesRobustFit.back()->SetLineColor(++iRobCirCol);
                        drawPotCirclesRobustFit.back()->SetLineStyle(4);
                        drawPotCirclesRobustFit.back()->Draw("same");
                        tmpCircleFound3D->cd(3);
                        drawPotCirclesRobustFit.back()->Draw("same");

                }


                if(helixesfit.size()>0)
                        for (RbstFitPntsLpsPntsRel::iterator xyzpLoopRel_it = xyzpLoopRel.begin(); xyzpLoopRel_it != xyzpLoopRel.end(); ++xyzpLoopRel_it ) {
                                XYZP &tmpXYZP = xyzp.at(xyzpLoopRel_it->first);
                                if ( tmpXYZP.use() ) {
                                        double tmpDist = fabs( sqrt( pow(tmpXYZP._pos.x()-helixesfit.back()._center.x(),2) +
                                                        pow(tmpXYZP._pos.y()-helixesfit.back()._center.y(),2) ) -
                                                        helixesfit.back()._radius );
                                        if ( tmpDist<3.0*tmpXYZP._rerr ) {
                                                ((TEllipse *) drawPotCircPnts3D.at(xyzpLoopRel_it->second.first)->At(xyzpLoopRel_it->second.second))->Draw("same");
                                        }
                                }
                        }
                tmpCircleFound3D->cd(4);
                float * selZarr = new float [selectedPnts.size()];
                float * selRarr = new float [selectedPnts.size()];
                float * selErrZarr = new float [selectedPnts.size()];
                float * selErrRarr = new float [selectedPnts.size()];
                int iSelPnt=0;
                for (points3D::iterator selectedPnts_it = selectedPnts.begin(); selectedPnts_it != selectedPnts.end(); ++selectedPnts_it) {
                        selZarr[iSelPnt] = selectedPnts_it->_pos.z()/CLHEP::cm;
                        selRarr[iSelPnt] = selectedPnts_it->_pos.perp()/CLHEP::cm;
                        selErrZarr[iSelPnt] = selectedPnts_it->_sigmaz/CLHEP::cm;
                        selErrRarr[iSelPnt] = selectedPnts_it->_sigmax/CLHEP::cm;
                        ++iSelPnt;
                }
                TGraphErrors drawSelPotCircPnts3D_L( selectedPnts.size(),selZarr,selRarr,selErrZarr,selErrRarr );
                if(selectedPnts.size()){
                        //drawSelPotCircPnts3D_L.SetMarkerColor(iCol);
                        drawSelPotCircPnts3D_L.SetMarkerSize(0.8);
                        drawSelPotCircPnts3D_L.SetMarkerStyle(8);
                        drawSelPotCircPnts3D_L.Draw("P");
                }
                delete [] selZarr;
                delete [] selRarr;
                delete [] selErrZarr;
                delete [] selErrRarr;

                std::vector<TGraph *> drawTracksT;
                std::vector<TGraph *> drawTracksL;

                iRobCirCol = 0;
                iTrackFound = 0;
                for(std::vector<HelixTraj>::iterator ittrk=trkhelixes.begin();ittrk!=trkhelixes.end();ittrk++){
                        iRobCirCol++;
                        TGraph *xy=new TGraph(400);
                        TGraph *zr=new TGraph(400);
                        double flt0=(*ittrk).zFlight(_minZ);
                        double flt1=(*ittrk).zFlight(_maxZ);
                        // std::cout<<flt0<<" "<<flt1<<std::endl;
                        for(int i=0;i<1000;i++){
                                HepPoint pos=(*ittrk).position(flt0+(flt1-flt0)*i/1000);
                                // std::cout<<i<<" "<<pos<<endl;
                                xy->SetPoint(i,pos.x()/CLHEP::cm,pos.y()/CLHEP::cm);
                                zr->SetPoint(i,pos.z()/CLHEP::cm,std::hypot(pos.x()/CLHEP::cm,pos.y()/CLHEP::cm));
                        }
                        xy->SetMarkerColor(iRobCirCol);
                        zr->SetMarkerColor(iRobCirCol);
                        xy->SetMarkerSize(0.1);
                        xy->SetMarkerStyle(6);
                        zr->SetMarkerSize(0.1);
                        zr->SetMarkerStyle(6);
                        tmpCircleFound3D->cd(2);
                        zr->Draw("P");
                        tmpCircleFound3D->cd(4);
                        zr->Draw("P");

                        drawTracksT.push_back(xy);
                        drawTracksL.push_back(zr);

                        int iMarckerLoop = 7;
                        cout<<"For i track "<<iTrackFound<<" there are "<<trkhelixesLoops[iTrackFound].size()<<" Loops to draw"<<endl;
                        for (std::vector<HelixTraj>::iterator trkhelixesLoop_it = trkhelixesLoops[iTrackFound].begin(); trkhelixesLoop_it != trkhelixesLoops[iTrackFound].end(); ++trkhelixesLoop_it) {
                                cout<<"Drawing single Loop for track "<<iTrackFound<<endl;
                                TGraph *xyLoop=new TGraph(400);
                                TGraph *zrLoop=new TGraph(400);
                                double flt0Loop=(*trkhelixesLoop_it).zFlight(_minZ);
                                double flt1Loop=(*trkhelixesLoop_it).zFlight(_maxZ);
                                // std::cout<<flt0<<" "<<flt1<<std::endl;
                                for(int i=0;i<1000;i++){
                                        HepPoint pos=(*trkhelixesLoop_it).position(flt0Loop+(flt1Loop-flt0Loop)*i/1000);
                                        // std::cout<<i<<" "<<pos<<endl;
                                        xyLoop->SetPoint(i,pos.x()/CLHEP::cm,pos.y()/CLHEP::cm);
                                        zrLoop->SetPoint(i,pos.z()/CLHEP::cm,std::hypot(pos.x()/CLHEP::cm,pos.y()/CLHEP::cm));
                                }
                                xyLoop->SetMarkerColor(iRobCirCol);
                                zrLoop->SetMarkerColor(iRobCirCol);
                                xyLoop->SetMarkerSize(0.1);
                                xyLoop->SetMarkerStyle(iMarckerLoop);
                                zrLoop->SetMarkerSize(0.1);
                                zrLoop->SetMarkerStyle(iMarckerLoop);
                                tmpCircleFound3D->cd(2);
                                zrLoop->Draw("P");
                                tmpCircleFound3D->cd(4);
                                zrLoop->Draw("P");

                                drawTracksT.push_back(xyLoop);
                                drawTracksL.push_back(zrLoop);
                                ++iMarckerLoop;
                        }
                        ++iTrackFound;
                }



                //temp El Visual
                std::vector<float *> simXarr;
                std::vector<float *> simYarr;
                std::vector<float *> simZarr;
                std::vector<float *> simRarr;
                int ielTrack=0, trckCol=0;
                cout<<"Generate el tracks "<<genEltrks->size()<<endl;
                for ( std::vector<mu2e::VisibleGenElTrack>::const_iterator genEltrk_it = genEltrks->begin(); genEltrk_it!= genEltrks->end(); ++genEltrk_it ){
                        VisibleGenElTrack &iEltrk = const_cast<VisibleGenElTrack &>(*genEltrk_it);

                        simXarr.push_back(new float [iEltrk.getNumOfHit()] );
                        simYarr.push_back(new float [iEltrk.getNumOfHit()] );
                        simZarr.push_back(new float [iEltrk.getNumOfHit()] );
                        simRarr.push_back(new float [iEltrk.getNumOfHit()] );

                        /*
                        unsigned short &nloops = iEltrk.getNumOfLoops();
                        cout<<"N loops "<<nloops<<endl;
                        convElNLoop=nloops;
                        for ( unsigned int ilp=0; ilp<nloops; ilp++ ){
                                GenElHitData& hdil = iEltrk.getithLoopHit(ilp);
                                ptMeV = sqrt( pow(hdil._hitMomentum[0],2) + pow(hdil._hitMomentum[1],2) );
                                rho   = ptMeV/(B*0.3);
                                cout<<ilp<<" -th loop: p_t "<<ptMeV<<" rho mm "<<rho<<endl;
                                CirCenter.set(hdil._hitPoint.getX(),hdil._hitPoint.getY(),hdil._hitPoint.getZ());
                                radDir.setX(hdil._hitMomentum.getX());
                                radDir.setY(hdil._hitMomentum.getY());
                                radDir.rotate( ( (hdil._hitMomentum.getZ()>=0.0) ? 90.0 : -90.0 )*CLHEP::degree );
                                radDir=radDir.unit();
                                CirCenter=CirCenter+rho*radDir;
                                cout<<" Hit Pos "<<hdil._hitPoint<<" Circ center "<<CirCenter<<endl;

                        }*/
                        if ( iEltrk.isConversionEl() ) {
                                trckCol = kGreen+2;
                                for ( unsigned int iElHit=0; iElHit<iEltrk.getNumOfHit(); iElHit++) {
                                        GenElHitData& genElhit = iEltrk.getHit((int)iElHit);
                                        simXarr.at(ielTrack)[iElHit] = genElhit._hitPoint.x()/CLHEP::cm;
                                        simYarr.at(ielTrack)[iElHit] = genElhit._hitPoint.y()/CLHEP::cm;
                                        simZarr.at(ielTrack)[iElHit] = genElhit._hitPoint.z()/CLHEP::cm;
                                        simRarr.at(ielTrack)[iElHit] = genElhit._hitPoint.perp()/CLHEP::cm;

                                        //cout<<"El Hit data "<<genElhit._hitPoint<<endl;
                                        if (genElhit._isFirst) {
                                        }
                                }

                        }
                        else {
                                /*
                                convElHitOvrlppd=0;
                                convElHitOvrlppdByP=0;
                                for ( unsigned int iElHit=0; iElHit<iEltrk.getNumOfHit(); iElHit++) {
                                        GenElHitData& genElhit = iEltrk.getHit((int)iElHit);
                                        if (genElhit._isOverlapped) {
                                                convElHitOvrlppd++;
                                                if (genElhit._isOvrlpByProton) convElHitOvrlppdByP++;
                                        }
                                }
                                _hNhitOverlapNoConvElEv_nh->Fill(iEltrk.getNumOfHit(),convElHitOvrlppd);
                                _hNhitOverlapByPNoConvElEv_nh->Fill(iEltrk.getNumOfHit(),convElHitOvrlppdByP);
                                 */
                        }

                        drawTracksT.push_back(new TGraph( iEltrk.getNumOfHit(), simXarr.back(), simYarr.back() ) );
                        drawTracksL.push_back(new TGraph( iEltrk.getNumOfHit(), simZarr.back(), simRarr.back() ) );

                        tmpCircleFound3D->cd(1);
                        drawTracksT.back()->SetMarkerColor(trckCol);
                        drawTracksT.back()->SetMarkerSize(0.5);
                        drawTracksT.back()->SetMarkerStyle(6);
                        drawTracksT.back()->Draw("P");
                        tmpCircleFound3D->cd(2);
                        drawTracksL.back()->SetMarkerColor(trckCol);
                        drawTracksL.back()->SetMarkerSize(0.5);
                        drawTracksL.back()->SetMarkerStyle(6);
                        drawTracksL.back()->Draw("P");


                        ++ielTrack;

                }


                //TCanvas *zClustCanv = new TCanvas();
                //tmpZhist.Draw();


                cerr << "Double click in the canvas_Fake to continue:" ;
                _fakeCanvas->cd();
                //TLatex *printEvN = new TLatex(0.15,0.4,Form("Current Event: %d - Tpeak %d",event.id().event(),ipeak));
                TLatex *printEvN = new TLatex(0.15,0.4,"Double click here to continue:");
                //printEvN->SetTextFont(50);
                printEvN->SetTextSizePixels(160);
                printEvN->Draw();
                _fakeCanvas->Modified();
                _fakeCanvas->Update();
                _fakeCanvas->WaitPrimitive();
                cerr << endl;
                delete printEvN;

                delete tmpRCP_Min;
                delete tmpRCP_Pl;


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

                delete tmpCircleFound3D;
                ///*int*/ iLoop = 0;
                for (std::vector<TClonesArray *>::iterator drawPotCircPnts3D_it = drawPotCircPnts3D.begin(); drawPotCircPnts3D_it != drawPotCircPnts3D.end(); ++drawPotCircPnts3D_it ) {
                        (*drawPotCircPnts3D_it)->Delete();
                        //if ( drawPotCircPnts3D_L.at(iLoop)!=0x0 ) drawPotCircPnts3D_L.at(iLoop)->Delete();
                        //if ( zarr.at(iLoop)!=0x0 ) delete [] zarr.at(iLoop);
                        //if ( rarr.at(iLoop)!=0x0 ) delete [] rarr.at(iLoop);
                        //if ( errzarr.at(iLoop)!=0x0 ) delete [] errzarr.at(iLoop);
                        //if ( errrarr.at(iLoop)!=0x0 ) delete [] errrarr.at(iLoop);
                }

                for (std::vector<TEllipse *>::iterator drawPotCirclesRobustFit_it = drawPotCirclesRobustFit.begin(); drawPotCirclesRobustFit_it != drawPotCirclesRobustFit.end(); ++drawPotCirclesRobustFit_it ) {
                        delete *drawPotCirclesRobustFit_it;
                }

                for_each(drawTracksT.begin(), drawTracksT.end(), std::bind2nd(std::mem_fun(&TGraph::Delete),""));
                for_each(drawTracksL.begin(), drawTracksL.end(), std::bind2nd(std::mem_fun(&TGraph::Delete),""));
                for_each(drawPotCircPnts3D_L.begin(), drawPotCircPnts3D_L.end(), std::bind2nd(std::mem_fun(&TGraphErrors::Delete),""));

                //delete zClustCanv;

        }

        for (std::vector<confMapDraw*>::iterator hCmapMin_it = rCMap_Min.begin(); hCmapMin_it != rCMap_Min.end(); ++hCmapMin_it ) {
                delete (*hCmapMin_it);
        }
        for (std::vector<confMapDraw*>::iterator hCmapPl_it = rCMap_Pl.begin(); hCmapPl_it != rCMap_Pl.end(); ++hCmapPl_it ) {
                delete (*hCmapPl_it);
        }
        rCMap_Min.clear();
        rCMap_Pl.clear();

        for (std::vector<confMapDraw*>::iterator hCmapMin_it = rCMap_Min2.begin(); hCmapMin_it != rCMap_Min2.end(); ++hCmapMin_it ) {
                delete (*hCmapMin_it);
        }
        for (std::vector<confMapDraw*>::iterator hCmapPl_it = rCMap_Pl2.begin(); hCmapPl_it != rCMap_Pl2.end(); ++hCmapPl_it ) {
                delete (*hCmapPl_it);
        }
        rCMap_Min2.clear();
        rCMap_Pl2.clear();

}

}  // end namespace mu2e

using mu2e::ITTrackReco;
DEFINE_ART_MODULE(ITTrackReco);
