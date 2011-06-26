//
// performance evaluation of the Bkg rejection modules
//
// $Id: EvalBkgTrackRejection_module.cc,v 1.1 2011/06/26 00:03:27 tassiell Exp $
// $Author: tassiell $
// $Date: 2011/06/26 00:03:27 $
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
//#include "TrackerGeom/inc/Tracker.hh"
//#include "TrackerGeom/inc/Straw.hh"
//#include "ITrackerGeom/inc/Cell.hh"
//#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
//#include "MCDataProducts/inc/GenParticleCollection.hh"
//#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
//#include "MCDataProducts/inc/SimParticleCollection.hh"
//#include "MCDataProducts/inc/StepPointMCCollection.hh"
//#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
//#include "RecoDataProducts/inc/StrawHitCollection.hh"
//#include "ITrackerGeom/inc/ITracker.hh"
//#include "TTrackerGeom/inc/TTracker.hh"
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
//
//    // Label of the module that made the hits.
//    std::string _makerModuleLabel;
//
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

    // Some ugly but necessary ROOT related bookkeeping:

    // The job needs exactly one instance of TApplication.  See note 1.
    auto_ptr<TApplication> _application;

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
//    _makerModuleLabel(pset.get<string>("makerModuleLabel")),
//    _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),

//    _fakeCanvas(0),

    // Some ugly but necessary ROOT related bookkeeping.
    _application(0),
    _directory(0){

 }

  void EvalBkgTrackRejection::beginJob(){

	  cout<<"Starting EvalBkgTrackRejection jos!"<<endl;

    // Get access to the TFile service and save current directory for later use.
    art::ServiceHandle<art::TFileService> tfs;

    // If needed, create the ROOT interactive environment. See note 1.
    if ( !gApplication ){
      int    tmp_argc(0);
      char** tmp_argv(0);
      _application = auto_ptr<TApplication>(new TApplication( "noapplication", &tmp_argc, tmp_argv ));
    }

    gStyle->SetPalette(1);
    gROOT->SetStyle("Plain");

    //_fakeCanvas = new TCanvas("canvas_Fake","double click for next event",300,100);

    // See note 3.
    _directory = gDirectory;

  }

  void EvalBkgTrackRejection::analyze(art::Event const& event ) {


//    const Tracker& tracker = getTrackerOrThrow();
//    const TTracker &ttr = static_cast<const TTracker&>( tracker );
//    const std::vector<Device> ttrdev = ttr.getDevices();

//    art::Handle<StrawHitCollection> pdataHandle;
//    event.getByLabel(_makerModuleLabel,pdataHandle);
//    StrawHitCollection const* hits = pdataHandle.product();


    art::Handle<VisibleGenElTrackCollection> genEltrksHandle;
    event.getByLabel(_extractElectronsData,genEltrksHandle);
    VisibleGenElTrackCollection const* genEltrks = genEltrksHandle.product();

    art::Handle<TrackerHitTimeClusterCollection> tclustHandle;
    event.getByLabel(_timeRejecterModuleLabel,tclustHandle);
    TrackerHitTimeClusterCollection const* tclusts = tclustHandle.product();

    art::Handle<SctrSttnClusterGroupCollection> gclusgtHandle;
    event.getByLabel(_geomRejecterModuleLabel,gclusgtHandle);
    SctrSttnClusterGroupCollection const* gclustgs = gclusgtHandle.product();



    size_t nTimeClusPerEvent = tclusts->size();

    std::vector<mu2e::VisibleGenElTrack>::const_iterator genEltrk_it;

    double ptMeV, rho;
    double B=1.0;
    CLHEP::Hep2Vector radDir;
    HepGeom::Point3D<double> CirCenter;

//    for ( genEltrk_it = genEltrks->begin(); genEltrk_it!= genEltrks->end(); ++genEltrk_it ){
//            VisibleGenElTrack &iEltrk = const_cast<VisibleGenElTrack &>(*genEltrk_it);
//            unsigned short &nloops = iEltrk.getNumOfLoops();
//            cout<<"N loops "<<nloops<<endl;
//            for ( unsigned int ilp=0; ilp<nloops; ilp++ ){
//                    GenElHitData& hdil = iEltrk.getithLoopHit(ilp);
//                    ptMeV = sqrt( pow(hdil._hitMomentum[0],2) + pow(hdil._hitMomentum[1],2) );
//                    rho   = ptMeV/(B*0.3);
//                    cout<<ilp<<" -th loop: p_t "<<ptMeV<<" rho mm "<<rho<<endl;
//                    CirCenter.set(hdil._hitPoint.getX(),hdil._hitPoint.getY(),hdil._hitPoint.getZ());
//                    radDir.setX(hdil._hitMomentum.getX());
//                    radDir.setY(hdil._hitMomentum.getY());
//                    radDir.rotate( ( (hdil._hitMomentum.getZ()>=0.0) ? 90.0 : -90.0 )*CLHEP::degree );
//                    radDir=radDir.unit();
//                    CirCenter=CirCenter+rho*radDir;
//                    cout<<" Hit Pos "<<hdil._hitPoint<<" Circ center "<<CirCenter<<endl;
//
//            }
//            //cout<<"All track data:"<<endl;
//            //cout<<iEltrk;
//
//    }


    cout<<"--------------------------- Analysis ---------------------------"<<endl;
    cout<<"Event: "<<event.id().event()<<endl;
    cout<<"N of group of Ttracker cluster of Hit found that could be tracks: "<<gclustgs->size()<<endl;
    SctrSttnClusterGroupCollection::const_iterator gclustgs_it;
    std::vector<SectorStationCluster>::const_iterator iclustg_it;
    int igroup=0, iclust;
    for ( gclustgs_it=gclustgs->begin(); gclustgs_it!=gclustgs->end(); ++gclustgs_it ) {
            cout<<"i-th group :"<<igroup<<" : "<<endl<<*gclustgs_it;
            iclust=0;
            for ( iclustg_it=gclustgs_it->_selectedClusters.begin(); iclustg_it!=gclustgs_it->_selectedClusters.end(); ++iclustg_it ){
                    cout<<"i-th clust :"<<iclust<<" : "<<endl<<*iclustg_it;
                    iclust++;
            }
            igroup++;
    }




//    cerr << "Double click in the canvas_Fake to continue:" ;
//    _fakeCanvas->cd();
//    TLatex *printEvN = new TLatex(0.15,0.4,Form("Current Event: %d",event.id().event()));
//    printEvN->SetTextFont(62);
//    printEvN->SetTextSizePixels(180);
//    printEvN->Draw();
//    _fakeCanvas->WaitPrimitive();
//    cerr << endl;
//    delete printEvN;


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
