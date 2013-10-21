//
// example of a module to read Data of the Electrons tracks that came from the targets
//
// $Id: ReadExtractElectronsData_module.cc,v 1.9 2013/10/21 21:01:22 kutschke Exp $
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
#include "MCDataProducts/inc/VisibleGenElTrack.hh"
#include "MCDataProducts/inc/VisibleGenElTrackCollection.hh"

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
#include "TProfile.h"
#include "TTree.h"

using namespace std;

namespace mu2e {

  typedef art::Ptr<StrawHit> StrawHitPtr;

  class ReadExtractElectronsData : public art::EDAnalyzer {
  public:

    explicit ReadExtractElectronsData(fhicl::ParameterSet const& pset);
    virtual ~ReadExtractElectronsData() {
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

    TProfile *_hNhitOverlapConvElEv_nh;
    TProfile *_hNhitOverlapByPConvElEv_nh;
    TProfile *_hNhitOverlapNoConvElEv_nh;
    TProfile *_hNhitOverlapByPNoConvElEv_nh;

    TTree    *_dataConvEl;
    unsigned int  runID, eventID, convElGoodNHit, convElTotNHit, convElNLoop, el_AtTracker_LID;
    float         el_Start_px, el_Start_py, el_Start_pz, el_AtTracker_px, el_AtTracker_py, el_AtTracker_pz, el_AtTracker_rad, el_AtTracker_z;
    float         el_ExitTracker_px, el_ExitTracker_py, el_ExitTracker_pz, el_ExitTracker_rad, el_ExitTracker_z;
    float         convElFHitTime/*, convElStartTime*/, convElLstHitTime;
    double        hitDelaE[1000];
    float         signalTimeCoin[1000], lastSignalTimeCoin[1000];

    // End: run time parameters

//    TCanvas*      _fakeCanvas;

    // Some ugly but necessary ROOT related bookkeeping:

    // The job needs exactly one instance of TApplication.  See note 1.
    unique_ptr<TApplication> _application;

    // Save directory from beginJob so that we can go there in endJob. See note 3.
    TDirectory* _directory;

  };

  ReadExtractElectronsData::ReadExtractElectronsData(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),

    // Run time parameters
    _extractElectronsData(pset.get<string>("elextractModuleLabel")),
//    _moduleLabel(pset.get<string>("module_label")),/*@module_label*/
//    _g4ModuleLabel(pset.get<string>("g4ModuleLabel")),
//    _trackerStepPoints(pset.get<string>("trackerStepPoints","tracker")),
//    _makerModuleLabel(pset.get<string>("makerModuleLabel")),
//    _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
    _hNhitOverlapConvElEv_nh(0),
    _hNhitOverlapByPConvElEv_nh(0),
    _hNhitOverlapNoConvElEv_nh(0),
    _hNhitOverlapByPNoConvElEv_nh(0),
    _dataConvEl(0),

//    _fakeCanvas(0),

    // Some ugly but necessary ROOT related bookkeeping.
    _application(nullptr),
    _directory(0){
          runID=eventID=convElGoodNHit=convElTotNHit=convElNLoop=0;
          el_Start_px=el_Start_py=el_Start_pz=el_AtTracker_px=el_AtTracker_py=el_AtTracker_pz=0.0;
          convElFHitTime=/*convElStartTime=*/0.0;
 }

  void ReadExtractElectronsData::beginJob(){

	  cout<<"Starting DisplayData jos!"<<endl;

    // Get access to the TFile service and save current directory for later use.
    art::ServiceHandle<art::TFileService> tfs;

    // If needed, create the ROOT interactive environment. See note 1.
    if ( !gApplication ){
      int    tmp_argc(0);
      char** tmp_argv(0);
      _application = unique_ptr<TApplication>(new TApplication( "noapplication", &tmp_argc, tmp_argv ));
    }

    gStyle->SetPalette(1);
    gROOT->SetStyle("Plain");

    _hNhitOverlapConvElEv_nh       = tfs->make<TProfile>( "hNhitOverlapConvElEv_nh", "Profile of the number of Hits overlapped by another particle for trackable conv. electrons", 120, 6, 126 );
    _hNhitOverlapByPConvElEv_nh    = tfs->make<TProfile>( "hNhitOverlapByPConvElEv_nh", "Profile of the number of Hits overlapped by proton for trackable conv. electrons", 120, 6, 126 );
    _hNhitOverlapNoConvElEv_nh     = tfs->make<TProfile>( "hNhitOverlapNoConvElEv_nh", "Profile of the number of Hits overlapped by another particle for not conv. electrons", 120, 6, 126 );
    _hNhitOverlapByPNoConvElEv_nh  = tfs->make<TProfile>( "hNhitOverlapByPNoConvElEv_nh", "Profile of the number of Hits overlapped by proton for not conv. electrons", 120, 6, 126 );

    _dataConvEl                    = tfs->make<TTree>( "dataConvEl", "data for Conversion Electron" );
    _dataConvEl->Branch("Run",&runID,"runID/i");
    _dataConvEl->Branch("Event",&eventID,"eventID/i");
    _dataConvEl->Branch("ConvElGoodNHit",&convElGoodNHit,"convElGoodNHit/i");
    _dataConvEl->Branch("ConvElTotNHit",&convElTotNHit,"convElTotNHit/i");
    _dataConvEl->Branch("ConvElNLoop",&convElNLoop,"convElNLoop/i");
    _dataConvEl->Branch("ConvEl_Start_px",&el_Start_px,"el_Start_px/F");
    _dataConvEl->Branch("ConvEl_Start_py",&el_Start_py,"el_Start_py/F");
    _dataConvEl->Branch("ConvEl_Start_pz",&el_Start_pz,"el_Start_pz/F");
    _dataConvEl->Branch("ConvEl_Tracker_px",&el_AtTracker_px,"el_AtTracker_px/F");
    _dataConvEl->Branch("ConvEl_Tracker_py",&el_AtTracker_py,"el_AtTracker_py/F");
    _dataConvEl->Branch("ConvEl_Tracker_pz",&el_AtTracker_pz,"el_AtTracker_pz/F");
    _dataConvEl->Branch("ConvEl_Tracker_layerID",&el_AtTracker_LID,"el_AtTracker_LID/i");
    _dataConvEl->Branch("ConvEl_Tracker_rad",&el_AtTracker_rad,"el_AtTracker_rad/F");
    _dataConvEl->Branch("ConvEl_Tracker_z",&el_AtTracker_z,"el_AtTracker_z/F");
    _dataConvEl->Branch("ConvEl_FrstHit_Time",&convElFHitTime,"convElFHitTime/F");
    _dataConvEl->Branch("ConvEl_ExitTracker_px",&el_ExitTracker_px,"el_ExitTracker_px/F");
    _dataConvEl->Branch("ConvEl_ExitTracker_py",&el_ExitTracker_py,"el_ExitTracker_py/F");
    _dataConvEl->Branch("ConvEl_ExitTracker_pz",&el_ExitTracker_pz,"el_ExitTracker_pz/F");
    _dataConvEl->Branch("ConvEl_ExitTracker_rad",&el_ExitTracker_rad,"el_ExitTracker_rad/F");
    _dataConvEl->Branch("ConvEl_ExitTracker_z",&el_ExitTracker_z,"el_ExitTracker_z/F");
    _dataConvEl->Branch("ConvEl_LstHit_Time",&convElLstHitTime,"convElLstHitTime/F");
    _dataConvEl->Branch("ConvEl_Eloss_btwn_Hits",hitDelaE,"hitDelaE[convElTotNHit]/D");
    //_dataConvEl->Branch("ConvEl_Start_Time",&convElStartTime,"convElStartTime/F");
    _dataConvEl->Branch("ConvEl_Frst_Time_Coinc",signalTimeCoin,"signalTimeCoin[convElTotNHit]/F");
    _dataConvEl->Branch("ConvEl_Last_Time_Coinc",lastSignalTimeCoin,"lastSignalTimeCoin[convElTotNHit]/F");

    //_fakeCanvas = new TCanvas("canvas_Fake","double click for next event",300,100);

    // See note 3.
    _directory = gDirectory;

  }

  void ReadExtractElectronsData::analyze(art::Event const& event ) {


//    const Tracker& tracker = getTrackerOrThrow();
//    const TTracker &ttr = static_cast<const TTracker&>( tracker );
//    const std::vector<Device> ttrdev = ttr.getDevices();

//    art::Handle<StrawHitCollection> pdataHandle;
//    event.getByLabel(_makerModuleLabel,pdataHandle);
//    StrawHitCollection const* hits = pdataHandle.product();


    art::Handle<VisibleGenElTrackCollection> genEltrksHandle;
    event.getByLabel(_extractElectronsData,genEltrksHandle);
    VisibleGenElTrackCollection const* genEltrks = genEltrksHandle.product();

    std::vector<mu2e::VisibleGenElTrack>::const_iterator genEltrk_it;

    runID=event.run();
    eventID=event.event();

    double ptMeV, rho;
    double B=1.0;
    CLHEP::Hep2Vector radDir;
    HepGeom::Point3D<double> CirCenter;
    int convElHitOvrlppd=0;
    int convElHitOvrlppdByP=0;

    for ( genEltrk_it = genEltrks->begin(); genEltrk_it!= genEltrks->end(); ++genEltrk_it ){
            VisibleGenElTrack &iEltrk = const_cast<VisibleGenElTrack &>(*genEltrk_it);
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

            }
            //cout<<"All track data:"<<endl;
            //cout<<iEltrk;
            if ( iEltrk.isConversionEl() ) {
                    convElHitOvrlppd=0;
                    convElHitOvrlppdByP=0;
                    convElTotNHit=iEltrk.getNumOfHit();
                    convElGoodNHit=0;
                    el_Start_px=iEltrk.getTrkLrntzVec().px();
                    el_Start_py=iEltrk.getTrkLrntzVec().py();
                    el_Start_pz=iEltrk.getTrkLrntzVec().pz();
                    //convElStartTime=;
                    GenElHitData& hdil = iEltrk.getFirstHit();//getHit(0);
                    el_AtTracker_px=hdil._hitMomentum[0];
                    el_AtTracker_py=hdil._hitMomentum[1];
                    el_AtTracker_pz=hdil._hitMomentum[2];
                    convElFHitTime=hdil._mcHitTime;
                    el_AtTracker_rad = sqrt(hdil._hitPoint.getX()*hdil._hitPoint.getX()+hdil._hitPoint.getY()*hdil._hitPoint.getY());
                    el_AtTracker_z = hdil._hitPoint.getZ();
                    el_AtTracker_LID = 0;

                    GenElHitData& hdil_l = iEltrk.getLastHit();//getHit(0);
                    el_ExitTracker_px=hdil_l._hitMomentum[0];
                    el_ExitTracker_py=hdil_l._hitMomentum[1];
                    el_ExitTracker_pz=hdil_l._hitMomentum[2];
                    convElLstHitTime=hdil_l._mcHitTime;
                    el_ExitTracker_rad = sqrt(hdil._hitPoint.getX()*hdil_l._hitPoint.getX()+hdil_l._hitPoint.getY()*hdil_l._hitPoint.getY());
                    el_ExitTracker_z = hdil_l._hitPoint.getZ();

                    double prevE(0.0), currE(0.0);
                    double pmas2 = iEltrk.getTrkLrntzVec().invariantMass2();
                    for ( unsigned int iElHit=0; iElHit<iEltrk.getNumOfHit(); iElHit++) {
                            //GenElHitData& genElhit = iEltrk.getHit((int)iElHit);
                            GenElHitData& genElhit = iEltrk.getHitTimeOrder((int)iElHit);
                            if (genElhit._isFirst) convElGoodNHit++;
                            if (genElhit._isOverlapped) {
                                    convElHitOvrlppd++;
                                    if (genElhit._isOvrlpByProton) convElHitOvrlppdByP++;
                            }
                            if (iElHit<1000) {
                                    signalTimeCoin[iElHit]=genElhit._iHit->time()-convElFHitTime;
                                    lastSignalTimeCoin[iElHit]=signalTimeCoin[iElHit]+genElhit._iHit->dt();
                            }
                            currE = sqrt( genElhit._hitMomentum.mag2() + pmas2 );
                            if (iElHit!=0) {
                                    hitDelaE[iElHit]=currE-prevE;
                            } else {
                                    hitDelaE[0]=0.0;
                            }
                            prevE = currE;
                    }
                    _hNhitOverlapConvElEv_nh->Fill(iEltrk.getNumOfHit(),convElHitOvrlppd);
                    _hNhitOverlapByPConvElEv_nh->Fill(iEltrk.getNumOfHit(),convElHitOvrlppdByP);
                    _dataConvEl->Fill();
            }
            else {
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
            }


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

  void ReadExtractElectronsData::endJob(){

    // cd() to correct root directory. See note 3.
    TDirectory* save = gDirectory;
    _directory->cd();

    // Write canvas.  See note 4.
//    _canvas->Write();

    // cd() back to where we were.  See note 3.
    save->cd();

  }


}  // end namespace mu2e

using mu2e::ReadExtractElectronsData;
DEFINE_ART_MODULE(ReadExtractElectronsData);
