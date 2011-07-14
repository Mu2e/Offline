//
// performance evaluation of the Bkg rejection modules
//
// $Id: EvalBkgTrackRejection_module.cc,v 1.2 2011/07/14 16:38:54 tassiell Exp $
// $Author: tassiell $
// $Date: 2011/07/14 16:38:54 $
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
#include "RecoDataProducts/inc/VisibleGenElTrack.hh"
#include "RecoDataProducts/inc/VisibleGenElTrackCollection.hh"
#include "RecoDataProducts/inc/TrackerHitTimeCluster.hh"
#include "RecoDataProducts/inc/TrackerHitTimeClusterCollection.hh"
#include "RecoDataProducts/inc/SectorStationCluster.hh"
#include "RecoDataProducts/inc/SctrSttnClusterGroup.hh"
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
//#include "TSpectrum.h"
//#include "TLatex.h"
#include "TTree.h"

using namespace std;

namespace mu2e {

  typedef art::Ptr<StrawHit> StrawHitPtr;

  class EvalBkgTrackRejection : public art::EDAnalyzer {
  public:

    explicit EvalBkgTrackRejection(fhicl::ParameterSet const& pset);
    virtual ~EvalBkgTrackRejection() {
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

    // Label of the module that made the hits.
    std::string _geomRejecterModuleLabel;

    // End: run time parameters

//    TCanvas*      _fakeCanvas;
    // Pointers to histograms, ntuples, TGraphs.
    TH1F*         _hNtrackableEv;
    TH1F*         _hNhitTrackableEv;
    TH1F*         _hNBestTimePeakEv;
    TH1F*         _hNhitBestTimePeakEv;
    TH2F*         _hLostHitBstTmPkEv;
    TH1F*         _hNoiseHitBstTmPkEv;

    TH1F*         _hNBestSGClInTmPkEv;
    TH1F*         _hSGClHitInBstTmPkEv;
    TH2F*         _hLostHitBstSGClBstTmPkEv;
    TH2F*         _hLostHitBstSGClEv;
    TH1F*         _hNoiseHitBstSGClEv;
    TH2F*         _hNLoopResBstSGClEv;

    TTree*        _dataEvalBkgRejec;
    unsigned int  runID, eventID, evNHit, convElNHit, convElNLoop;
    float         sel_ptMeV, sel_ptMeV_start, sel_ptMeV_end, sel_plMeV, sel_plMeV_start, sel_plMeV_end, convElFHitTime;
    int           bestTPcElNHit, bestTPNHit, nPeaksFound;
    int           bestClGcElNHit, bestClGNHit, bestClGcNCls, nPotentTracks;

    // Some ugly but necessary ROOT related bookkeeping:

    // The job needs exactly one instance of TApplication.  See note 1.
    //auto_ptr<TApplication> _application;

    // Save directory from beginJob so that we can go there in endJob. See note 3.
    TDirectory* _directory;

  };

  EvalBkgTrackRejection::EvalBkgTrackRejection(fhicl::ParameterSet const& pset) :

    // Run time parameters
    _extractElectronsData(pset.get<string>("elextractModuleLabel")),
    _timeRejecterModuleLabel(pset.get<string>("tRejecterModuleLabel")),
    _geomRejecterModuleLabel(pset.get<string>("gRejecterModuleLabel")),
//    _moduleLabel(pset.get<string>("module_label")),/*@module_label*/
//    _g4ModuleLabel(pset.get<string>("g4ModuleLabel")),
//    _trackerStepPoints(pset.get<string>("trackerStepPoints","tracker")),
    _makerModuleLabel(pset.get<string>("makerModuleLabel")),
//    _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),

//    _fakeCanvas(0),
    _hNtrackableEv(0),
    _hNhitTrackableEv(0),
    _hNBestTimePeakEv(0),
    _hNhitBestTimePeakEv(0),
    _hLostHitBstTmPkEv(0),
    _hNoiseHitBstTmPkEv(0),

    _hNBestSGClInTmPkEv(0),
    _hSGClHitInBstTmPkEv(0),
    _hLostHitBstSGClBstTmPkEv(0),
    _hLostHitBstSGClEv(0),
    _hNoiseHitBstSGClEv(0),
    _hNLoopResBstSGClEv(0),
    _dataEvalBkgRejec(0),

    // Some ugly but necessary ROOT related bookkeeping.
    //_application(0),
    _directory(0){
          runID=eventID=evNHit=convElNHit=convElNLoop=0;
          sel_ptMeV=sel_ptMeV_start=sel_ptMeV_end=sel_plMeV=sel_plMeV_start=sel_plMeV_end=convElFHitTime=0.0;
          bestTPcElNHit=bestTPNHit=nPeaksFound=0;
          bestClGcElNHit=bestClGNHit=bestClGcNCls=nPotentTracks=0;
 }

  void EvalBkgTrackRejection::beginJob(){

	  cout<<"Starting EvalBkgTrackRejection jos!"<<endl;

    // Get access to the TFile service and save current directory for later use.
    art::ServiceHandle<art::TFileService> tfs;

    // If needed, create the ROOT interactive environment. See note 1.
//    if ( !gApplication ){
//      int    tmp_argc(0);
//      char** tmp_argv(0);
//      _application = auto_ptr<TApplication>(new TApplication( "noapplication", &tmp_argc, tmp_argv ));
//    }

    gStyle->SetPalette(1);
    gROOT->SetStyle("Plain");

    // Create a histogram.
    _hNtrackableEv       = tfs->make<TH1F>( "hNtrackableEv",   "N of trackable signal electrons", 111, 49.75, 105.25  );
    _hNhitTrackableEv    = tfs->make<TH1F>( "hNhitTrackableEv",   "N of hits for each signal electrons track", 111, 49.75, 105.25  );
    _hNBestTimePeakEv    = tfs->make<TH1F>( "hNBestTimePeakEv",   "N of peak time that has the best agreement with trackable signal electrons", 111, 49.75, 105.25  );
    _hNhitBestTimePeakEv = tfs->make<TH1F>( "hNhitBestTimePeakEv",   "N of hits of the signal electrons track that are found in the best time peak", 111, 49.75, 105.25  );
    _hLostHitBstTmPkEv   = tfs->make<TH2F>( "hLostHitBstTmPkEv",   "N of lost hits of the signal electrons track that in the best time peak", 111, 49.75, 105.25, 20, 0, 20  );
    _hNoiseHitBstTmPkEv  = tfs->make<TH1F>( "hNoiseHitBstTmPkEv",   "N of hits of noise present in the the best time peak over signal electrons track hits", 111, 49.75, 105.25  );

    _hNtrackableEv      ->SetXTitle("pt [MeV]");
    _hNhitTrackableEv   ->SetXTitle("pt [MeV]");
    _hNBestTimePeakEv   ->SetXTitle("pt [MeV]");
    _hNhitBestTimePeakEv->SetXTitle("pt [MeV]");
    _hLostHitBstTmPkEv  ->SetXTitle("pt [MeV]");
    _hNoiseHitBstTmPkEv ->SetXTitle("pt [MeV]");

    _hNBestSGClInTmPkEv  = tfs->make<TH1F>( "hNBestSGClInTmPkEv",   "N of Best Selected Clusters Group in the best peak time", 111, 49.75, 105.25  );
    _hSGClHitInBstTmPkEv = tfs->make<TH1F>( "hSGClHitInBstTmPkEv",   "N of signal electron hit selected by best geom clustering algo from the best time peak", 111, 49.75, 105.25  );
    _hLostHitBstSGClBstTmPkEv = tfs->make<TH2F>( "hLostHitBstSGClBstTmPkEv",   "N of lost hits of the signal electrons track that in the best geom cluster from the best time peak", 111, 49.75, 105.25, 20, 0, 20  );
    _hLostHitBstSGClEv   = tfs->make<TH2F>( "hLostHitBstSGClEv",   "N of lost hits of the signal electrons track that in the best geom cluster", 111, 49.75, 105.25, 20, 0, 20  );
    _hNoiseHitBstSGClEv  = tfs->make<TH1F>( "hNoiseHitBstSGClEv",   "N of hits of noise present in the the best geom cluster over signal electrons track hits", 111, 49.75, 105.25  );
    _hNLoopResBstSGClEv  = tfs->make<TH2F>( "hNLoopResBstSGClEv",   "Res of N of track loop measured by the best geom cluster", 111, 49.75, 105.25, 21, -10.5, 10.5  );

    _hNBestSGClInTmPkEv  ->SetXTitle("pt [MeV]");
    _hSGClHitInBstTmPkEv ->SetXTitle("pt [MeV]");
    _hLostHitBstSGClBstTmPkEv ->SetXTitle("pt [MeV]");
    _hLostHitBstSGClEv ->SetXTitle("pt [MeV]");
    _hNoiseHitBstSGClEv ->SetXTitle("pt [MeV]");
    _hNLoopResBstSGClEv ->SetXTitle("pt [MeV]");

    _dataEvalBkgRejec = tfs->make<TTree>( "dataEvalBkgRejec", "data for Bkg Rejection performance evaluation" );
    _dataEvalBkgRejec->Branch("Run",&runID,"runID/i");
    _dataEvalBkgRejec->Branch("Event",&eventID,"eventID/i");
    _dataEvalBkgRejec->Branch("EvTotNHit",&evNHit,"evNHit/i");
    _dataEvalBkgRejec->Branch("ConvElNHit",&convElNHit,"convElNHit/i");
    _dataEvalBkgRejec->Branch("ConvElNLoop",&convElNLoop,"convElNLoop/i");
    _dataEvalBkgRejec->Branch("ConvEl_ptMeV_Tracker",&sel_ptMeV,"sel_ptMeV/F");
    _dataEvalBkgRejec->Branch("ConvEl_ptMeV_start",&sel_ptMeV_start,"sel_ptMeV_start/F");
    _dataEvalBkgRejec->Branch("ConvEl_ptMeV_end",&sel_ptMeV_end,"sel_ptMeV_end/F");
    _dataEvalBkgRejec->Branch("ConvEl_plMeV",&sel_plMeV,"sel_plMeV/F");
    _dataEvalBkgRejec->Branch("ConvEl_plMeV_start",&sel_plMeV_start,"sel_plMeV_start/F");
    _dataEvalBkgRejec->Branch("ConvEl_plMeV_end",&sel_plMeV_end,"sel_plMeV_end/F");
    _dataEvalBkgRejec->Branch("ConvEl_FrstHit_Time",&convElFHitTime,"convElFHitTime/F");
    _dataEvalBkgRejec->Branch("BestTPcElNHit",&bestTPcElNHit,"bestTPcElNHit/I");
    _dataEvalBkgRejec->Branch("BestTPNHit",&bestTPNHit,"bestTPNHit/I");
    _dataEvalBkgRejec->Branch("TotNPeaksFound",&nPeaksFound,"nPeaksFound/I");
    _dataEvalBkgRejec->Branch("BestClGcElNHit",&bestClGcElNHit,"bestClGcElNHit/I");
    _dataEvalBkgRejec->Branch("BestClGNHit",&bestClGNHit,"bestClGNHit/I");
    _dataEvalBkgRejec->Branch("BestClGcNCls",&bestClGcNCls,"bestClGcNCls/I");
    _dataEvalBkgRejec->Branch("NPotentTracksFound",&nPotentTracks,"nPotentTracks/I");


    //_fakeCanvas = new TCanvas("canvas_Fake","double click for next event",300,100);

    // See note 3.
    _directory = gDirectory;

  }

  void EvalBkgTrackRejection::analyze(art::Event const& event ) {


    const Tracker& tracker = getTrackerOrThrow();
    const TTracker &ttr = static_cast<const TTracker&>( tracker );
    const std::vector<Device> ttrdev = ttr.getDevices();

    art::Handle<StrawHitCollection> pdataHandle;
    event.getByLabel(_makerModuleLabel,pdataHandle);
    StrawHitCollection const* hits = pdataHandle.product();


    art::Handle<VisibleGenElTrackCollection> genEltrksHandle;
    event.getByLabel(_extractElectronsData,genEltrksHandle);
    VisibleGenElTrackCollection const* genEltrks = genEltrksHandle.product();

    art::Handle<TrackerHitTimeClusterCollection> tclustHandle;
    event.getByLabel(_timeRejecterModuleLabel,tclustHandle);
    TrackerHitTimeClusterCollection const* tclusts = tclustHandle.product();

    art::Handle<SctrSttnClusterGroupCollection> gclusgtHandle;
    event.getByLabel(_geomRejecterModuleLabel,gclusgtHandle);
    SctrSttnClusterGroupCollection const* gclustgs = gclusgtHandle.product();

    std::vector<mu2e::VisibleGenElTrack>::const_iterator genEltrk_it;

    cout<<"--------------------------- Analysis ---------------------------"<<endl;
    cout<<"Event: "<<event.id().event()<<endl;
    cout<<"----------------------------------------------------------------"<<endl;

    StrawId sid;
    bool goodEvent=false;
    VisibleGenElTrack *signaLEltrk;
    std::vector<StrawHitPtr> signElGoodHit;

    double ptMeV, rho;
    double sel_rho;
    double B=1.0;
    CLHEP::Hep2Vector radDir;
    HepGeom::Point3D<double> CirCenter;

    runID=event.run();
    eventID=event.event();
    evNHit=hits->size();
    bestTPcElNHit=bestTPNHit=nPeaksFound=-1;
    bestClGcElNHit=bestClGNHit=bestClGcNCls=nPotentTracks=-1;

    for ( genEltrk_it = genEltrks->begin(); genEltrk_it!= genEltrks->end(); ++genEltrk_it ){
            VisibleGenElTrack &iEltrk = const_cast<VisibleGenElTrack &>(*genEltrk_it);
            unsigned short &nloops = iEltrk.getNumOfLoops();
            cout<<"N loops "<<nloops<<endl;
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

            }

            if ( iEltrk.isConversionEl() ) {
                    if ( iEltrk.getNumOfHit()>=6 ) {
                            GenElHitData& hdil = iEltrk.getithLoopHit(0);
                            convElFHitTime = hdil._mcHitTime;
                            if (convElFHitTime>2400.0) goodEvent=false;
                            else {
                                    goodEvent=true;
                                    sel_ptMeV = sqrt( pow(hdil._hitMomentum[0],2) + pow(hdil._hitMomentum[1],2) );
                                    sel_plMeV = hdil._hitMomentum[2];
                                    _hNtrackableEv->Fill(sel_ptMeV);
                                    for ( int iElHit=0; iElHit<iEltrk.getNumOfHit(); iElHit++) {
                                            GenElHitData& genElhit = iEltrk.getHit(iElHit);
                                            if (genElhit._isFirst) signElGoodHit.push_back(genElhit._iHit);
                                    }
                                    convElNHit = signElGoodHit.size();
                                    convElNLoop = nloops;
                                    sel_ptMeV_start = sqrt( pow(iEltrk.getTrkLrntzVec().px(),2) + pow(iEltrk.getTrkLrntzVec().py(),2) );
                                    sel_plMeV_start = iEltrk.getTrkLrntzVec().pz();
                                    GenElHitData& lasth = iEltrk.getHit((int)(iEltrk.getNumOfHit()-1));
                                    sel_ptMeV_end = sqrt( pow(lasth._hitMomentum[0],2) + pow(lasth._hitMomentum[1],2) );
                                    sel_plMeV_end = lasth._hitMomentum[2];
                                    _hNhitTrackableEv->Fill(sel_ptMeV,signElGoodHit.size());
                                    signaLEltrk=&iEltrk;
                            }
                    }
            }
            //cout<<"All track data:"<<endl;
            //cout<<iEltrk;

    }


    //---------------- for time algorithm ----------------
    size_t nTimeClusPerEvent = tclusts->size();
    nPeaksFound=nTimeClusPerEvent;

    unsigned int elHitFound;
    std::multimap< unsigned int, size_t, less<unsigned int> > foundElHitInPeak;
    std::multimap< unsigned int, size_t, less<unsigned int> >::reverse_iterator BestFoundElInPeak_it;
    std::map<size_t, std::vector<StrawHitPtr> > signElGoodHitsInTmpks;

    if (goodEvent) {


            for (size_t ipeak=0; ipeak<nTimeClusPerEvent; ipeak++) {
                    TrackerHitTimeCluster const&  tclust(tclusts->at(ipeak));
                    cout<<"\t Hit in peak "<<tclust._selectedTrackerHits.size()<<endl;

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

                            //                //                cout<<"\t\t iHit in peak at "<<*iTCHit<<endl;
                            //                //StrawHit        const&      hit(hits->at(*iTCHit));
                            //                StrawHit        const&      hit=*(*iTCHit);
                            //                StrawIndex si = hit.strawIndex();
                            //                const Straw & str = tracker.getStraw(si);
                            //
                            //                // cout << "Getting informations about cells" << endl;
                            //
                            //                sid = str.id();
                    }

                    if ( elHitFound>0 ) {
                            foundElHitInPeak.insert( pair<unsigned int, size_t> (elHitFound,ipeak) );
                            signElGoodHitsInTmpks.insert( pair< size_t, std::vector<StrawHitPtr> > (ipeak,tmpElHits) );
                    }
            }

            if ( !foundElHitInPeak.empty() ) {
                    _hNBestTimePeakEv->Fill(sel_ptMeV);
                    BestFoundElInPeak_it=foundElHitInPeak.rbegin();
                    bestTPcElNHit=BestFoundElInPeak_it->first;
                    bestTPNHit=tclusts->at(BestFoundElInPeak_it->second)._selectedTrackerHits.size();
                    _hNhitBestTimePeakEv->Fill(sel_ptMeV,bestTPcElNHit);
                    _hLostHitBstTmPkEv->Fill(sel_ptMeV, (signElGoodHit.size() - bestTPcElNHit) );
                    _hNoiseHitBstTmPkEv->Fill(sel_ptMeV, (bestTPNHit - bestTPcElNHit) );
            }
    }
    //---------------- for geom algorithm ----------------

    cout<<"N of group of Ttracker cluster of Hit found that could be tracks: "<<gclustgs->size()<<endl;
    nPotentTracks=gclustgs->size();

    SctrSttnClusterGroupCollection::const_iterator gclustgs_it;
    std::vector<SectorStationCluster>::const_iterator iclustg_it;
    int igroup=0, iclust;

    //_hSGClHitInBstTmPkEv

    std::multimap< unsigned int, int, less<unsigned int> > foundElHitSClGrp;
    std::multimap< unsigned int, int, less<unsigned int> >::reverse_iterator BestFoundElInSClGrp_it;
    std::map< int, unsigned int > foundElSClGrpHit;
    std::map<int, std::vector<StrawHitPtr> > signElGoodHitsInSClGrps;
    std::map<int, std::vector<StrawHitPtr> >::iterator sgnElGdHtsInSClGrps_it;

    int sClGrpTotalHit, iHitInCl, iCl;//, elHitFoundInCl;

    if (goodEvent && !foundElHitInPeak.empty()) {

            for ( gclustgs_it=gclustgs->begin(); gclustgs_it!=gclustgs->end(); ++gclustgs_it ) {
                    cout<<"i-th group :"<<igroup<<" : "<<endl<<*gclustgs_it;
                    iclust=0;
                    for ( iclustg_it=gclustgs_it->_selectedClusters.begin(); iclustg_it!=gclustgs_it->_selectedClusters.end(); ++iclustg_it ){
                            cout<<"i-th clust :"<<iclust<<" : "<<endl<<*iclustg_it;
                            iclust++;
                    }

                    elHitFound=0;
                    std::vector<StrawHitPtr> tmpElHits;
                    sClGrpTotalHit=0;

                    if (gclustgs_it->_relatedTimeCluster.key()==BestFoundElInPeak_it->second) {
                            iCl=0;
                            for ( iclustg_it=gclustgs_it->_selectedClusters.begin(); iclustg_it!=gclustgs_it->_selectedClusters.end(); ++iclustg_it ){
                                    sClGrpTotalHit+=iclustg_it->_selectedTrackerHits.size();
                                    cout<<"--- N Hits in Cls "<<sClGrpTotalHit<<" partial "<<iclustg_it->_selectedTrackerHits.size()<<endl;
                                    //cout<<"hit in Cls "
                                    iHitInCl=0;
                                    //elHitFoundInCl=0;
                                    for ( std::vector<StrawHitPtr>::const_iterator iSClGrpHit_it=iclustg_it->_selectedTrackerHits.begin(); iSClGrpHit_it!=iclustg_it->_selectedTrackerHits.end(); ++iSClGrpHit_it ) {
                                            iHitInCl++;
                                            for ( std::vector<StrawHitPtr>::iterator iGeTpHit_it=signElGoodHitsInTmpks[BestFoundElInPeak_it->second].begin(); iGeTpHit_it!=signElGoodHitsInTmpks[BestFoundElInPeak_it->second].end(); ++iGeTpHit_it ) {
                                            //for ( std::vector<StrawHitPtr>::iterator iGoodSelHit_it=signElGoodHit.begin(); iGoodSelHit_it!=signElGoodHit.end(); ++iGoodSelHit_it ) {
                                                    //cout<<"Check in Cl "<<iCl<<" i-th hit "<<iHitInCl<<" with key "<<iSClGrpHit_it->key()<<" vs "<<iGeTpHit_it->key()/*iGoodSelHit_it->key()*/<<endl;
                                                    if ( iSClGrpHit_it->key()==iGeTpHit_it->key() ) {
                                                    //if ( iSClGrpHit_it->key()==iGoodSelHit_it->key() ) {
                                                            elHitFound++;
                                                            //elHitFoundInCl++;
                                                            tmpElHits.push_back(*iSClGrpHit_it);
                                                            break;
                                                    }
                                            }
                                    }
                                    //cout<<"hit Found in Cl "<<elHitFoundInCl<<endl;
                                    iCl++;
                            }
                            cout<<"Clgroup elHitFound "<<elHitFound<<endl;
                            if ( elHitFound>0 ) {
                                    cout<<"------------ el trk found "<<endl;
                                    foundElHitSClGrp.insert( pair<unsigned int, int> (elHitFound,igroup) );
                                    sgnElGdHtsInSClGrps_it=signElGoodHitsInSClGrps.find(igroup);
                                    foundElSClGrpHit.insert( pair<int,unsigned int> (igroup,sClGrpTotalHit) );
                                    if ( sgnElGdHtsInSClGrps_it==signElGoodHitsInSClGrps.end() ) signElGoodHitsInSClGrps.insert( pair< int, std::vector<StrawHitPtr> > (igroup,tmpElHits) );
                                    else if ( elHitFound>sgnElGdHtsInSClGrps_it->second.size() ) {
                                            cout<<"------------ removing previous el trk found "<<endl;
                                            signElGoodHitsInSClGrps.erase(sgnElGdHtsInSClGrps_it);
                                            signElGoodHitsInSClGrps.insert( pair< int, std::vector<StrawHitPtr> > (igroup,tmpElHits) );
                                    }
                            }
                    }

                    igroup++;
            }

            if ( !foundElHitSClGrp.empty() ) {
                    cout<<"------------ histogramming"<<endl;
                    _hNBestSGClInTmPkEv->Fill(sel_ptMeV);
                    BestFoundElInSClGrp_it=foundElHitSClGrp.rbegin();
                    bestClGcElNHit=BestFoundElInSClGrp_it->first;
                    bestClGNHit=foundElSClGrpHit.find(BestFoundElInSClGrp_it->second)->second;
                    cout<<"--- convElNHit "<<convElNHit<<" bestTPcElNHit "<<bestTPcElNHit<<" bestClGcElNHit "<<bestClGcElNHit<<endl;
                    //cout<<"---  bestClGNHit "<<foundElSClGrpHit.find(BestFoundElInSClGrp_it->second)->second<<" = "<<bestClGNHit<<endl;
                    bestClGcNCls=gclustgs->at(BestFoundElInSClGrp_it->second)._selectedClusters.size();
                    _hSGClHitInBstTmPkEv->Fill(sel_ptMeV,bestClGcElNHit);
                    _hLostHitBstSGClBstTmPkEv->Fill(sel_ptMeV, (BestFoundElInPeak_it->first - bestClGcElNHit) );
                    cout<<"---     !!!!!!!!!  "<<sel_ptMeV<<" - "<<BestFoundElInPeak_it->first<<" - "<<bestClGcElNHit<<endl;
                    _hLostHitBstSGClEv->Fill(sel_ptMeV, (signElGoodHit.size() - bestClGcElNHit) );
                    _hNoiseHitBstSGClEv->Fill(sel_ptMeV, (bestClGNHit - bestClGcElNHit) );
                    _hNLoopResBstSGClEv->Fill(sel_ptMeV, (convElNLoop - bestClGcNCls ) );
            }
    }

    if (goodEvent) _dataEvalBkgRejec->Fill();

//    cerr << "Double click in the canvas_Fake to continue:" ;
//    _fakeCanvas->cd();
//    TLatex *printEvN = new TLatex(0.15,0.4,Form("Current Event: %d",event.id().event()));
//    printEvN->SetTextFont(62);
//    printEvN->SetTextSizePixels(180);
//    printEvN->Draw();
//    _fakeCanvas->WaitPrimitive();
//    cerr << endl;
//    delete printEvN;

    cout<<"--------------------------- End of Analysis --------------------"<<endl;
    cout<<"Event: "<<event.id().event()<<endl;
    cout<<"----------------------------------------------------------------"<<endl;


  } // end analyze

  void EvalBkgTrackRejection::endJob(){

    // cd() to correct root directory. See note 3.
    TDirectory* save = gDirectory;
    _directory->cd();

    // Write canvas.  See note 4.
//    _canvas->Write();

    // cd() back to where we were.  See note 3.
    save->cd();

  }


}  // end namespace mu2e

using mu2e::EvalBkgTrackRejection;
DEFINE_ART_MODULE(EvalBkgTrackRejection);
