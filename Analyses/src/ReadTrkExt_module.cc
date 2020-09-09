//
// Read the tracks added to the event by KalFitTest_module.
//
//
// Original author MyeongJae Lee
//

// Mu2e includes.

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"

// ROOT incldues
#include "TH1F.h"

// Need this for the BaBar headers.
using namespace CLHEP;

// BaBar includes
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/TrkBase/TrkParticle.hh"
// mu2e tracking
#include "RecoDataProducts/inc/TrkFitDirection.hh"
#include "TrkDiag/inc/KalDiag.hh"
// C++ includes.
#include <iostream>
#include <string>

// This is fragile and needs to be last until CLHEP is
// properly qualified and included in the BaBar classes.
#include "RecoDataProducts/inc/KalRepCollection.hh"


#include "TrkExt/inc/TrkExtInstanceName.hh"
#include "RecoDataProducts/inc/TrkExtTraj.hh"
#include "RecoDataProducts/inc/TrkExtTrajCollection.hh"
#include "TrkExt/inc/TrkExtDiag.hh"

using namespace std;

namespace mu2e {

  class ReadTrkExt : public art::EDAnalyzer {
  public:

    explicit ReadTrkExt(fhicl::ParameterSet const& pset);
    virtual ~ReadTrkExt() { }

    void beginJob();
    void beginRun(art::Run const& run) override;
    void beginSubRun(art::SubRun const & lblock ) override; 
    void analyze(const art::Event& e);


  private:

    // Module label of the module that performed the fits.
    std::vector<std::string> _fitterModuleLabelArray;
    std::vector<int> _fitparticleArray;
    std::vector<int> _fitdirectionArray;
    std::string _trkextModuleLabel;
    
    // diagnostic of Kalman fit
    bool _recordKalFit;
    KalDiag _kdiag;
    TrkExtDiag _trkext;

    // Control level of printout.
    int _verbosity;
    int _maxPrint;
    bool _processEmpty; // whether or not to include MC info for empty events

    // Histograms
    TTree* _trkdiag;
    TTree* _extdiag;

    // internal variables
    TrkExtInstanceName _trkPatRecInstanceName;
    Int_t _trkid,_eventid, _updown, _hepid;

  };

  ReadTrkExt::ReadTrkExt(fhicl::ParameterSet const& pset):
    art::EDAnalyzer(pset),
    _fitterModuleLabelArray(pset.get<std::vector<std::string> >("fitterModuleLabelArray")),
    _fitparticleArray(pset.get<std::vector<int> >("fitparticleArray")),
    _fitdirectionArray(pset.get<std::vector<int> >("fitdirectionArray")),
    _trkextModuleLabel(pset.get<string>("trkextModuleLabel")),
    _recordKalFit(pset.get<bool>("recordKalFit", false)),
    _kdiag(pset.get<fhicl::ParameterSet>("KalDiag",fhicl::ParameterSet())),
    _trkext(pset.get<fhicl::ParameterSet>("TrkExt", fhicl::ParameterSet())),
    _verbosity(pset.get<int>("verbosity",0)),
    _maxPrint(pset.get<int>("maxPrint",0)),
    _processEmpty(pset.get<bool>("processEmpty",false)),
    _trkdiag(0),
    _extdiag(0)

    {
    // construct the data product instance name
    if (   _fitterModuleLabelArray.size() <=0 
        || _fitparticleArray.size() <= 0  
        || _fitdirectionArray.size() <=0        
        || _fitterModuleLabelArray.size() != _fitparticleArray.size()
        || _fitterModuleLabelArray.size() != _fitdirectionArray.size()
        || _fitdirectionArray.size() != _fitparticleArray.size() ) { 
      throw cet::exception("CONFIGURATION") 
        << "ReadTrkExt : Error in fitterModuleNameArray, fitparticleArray, or fitdirectionArray. They are not given properly or never given.";
    }

    for (unsigned int i = 0 ; i <_fitterModuleLabelArray.size() ; ++i) {
      _trkPatRecInstanceName.construct(
          (TrkParticle::type)(_fitparticleArray[i]),
          _fitdirectionArray[i],      
          _fitterModuleLabelArray[i]);
      if (_verbosity>=1) cout << "ReadTrkExt : module=" << _trkPatRecInstanceName.fitterName(i) << ", instance=" << _trkPatRecInstanceName.name(i) << " created" << endl;
    }
  }

  void ReadTrkExt::beginJob( ){
    art::ServiceHandle<art::TFileService> tfs;
    if (_recordKalFit) _trkdiag    = _kdiag.createTrkDiag();
    _extdiag = _trkext.createTrkExtDiag();
// add local branches
    _extdiag->Branch("eventid",&_eventid,"eventid/I");
    _extdiag->Branch("trkid",&_trkid,"trkid/I");
    _extdiag->Branch("updown",&_updown,"updown/I");
    _extdiag->Branch("hepid",&_hepid,"hepid/I");
    _eventid = 0;
  }

  void ReadTrkExt::beginRun(art::Run const& run) {
    _trkext.setRunInfo();
  }

  void ReadTrkExt::beginSubRun(art::SubRun const & lblock ) {
    _trkext.setSubRunInfo();
  }

  void ReadTrkExt::analyze(const art::Event& event) {

    if (_recordKalFit) _kdiag.findMCData(event);
    art::Handle<KalRepCollection> trksHandle;
    art::Handle<TrkExtTrajCollection> trkextHandle;
    _eventid = event.id().event();

    for (unsigned int instanceIter = 0 ; instanceIter < _trkPatRecInstanceName.size() ; ++instanceIter) {
      TrkExtInstanceNameEntry & instance = _trkPatRecInstanceName.get(instanceIter);
      if (instance.updown) _updown = 1;
      else _updown = 0;
      _hepid = instance.hepid;

      event.getByLabel(instance.fitterName, instance.name, trksHandle);
      event.getByLabel(_trkextModuleLabel, instance.name, trkextHandle);

      bool trkvalidity = trksHandle.isValid();
      bool extvalidity = trkextHandle.isValid();
      if (!trkvalidity && !extvalidity) { continue; }
      else if (trkvalidity && extvalidity) {;}
      else {
        throw cet::exception("CONFIGURATION") 
          << "ReadTrkExt Error : strange recon status : trk " << (trkvalidity?"valid":"not valid") << ", ext " <<(extvalidity?"valid":"not valid") << endl << "Maybe Track fit did not run?" << endl;
        return;
      }
  
      KalRepCollection const& trks = *trksHandle;
      TrkExtTrajCollection const & trkexts = *trkextHandle;

      if (trks.size() != trkexts.size()) {
        throw cet::exception("DATA")
          << "ReadTrkExt Error : number of entries from trkfit and trkext are not same : N(trkfit)=" << trks.size() << ", N(trkext)=" << trkexts.size() << endl;
        return;
      }
  
      if ( _verbosity > 0 && _eventid <= _maxPrint ){
        cout << "ReadKalmanFits  for event: " << event.id() << "  Number of fitted tracks: " << trks.size() << endl;
      }
  
      _trkid = -1;
      for ( size_t i=0; i< trks.size(); ++i ){
        _trkid = i;
        if (_recordKalFit) {
          KalRep const* krep = trks.get(i);
          if ( !krep ) {
            throw cet::exception("DATA")
              <<"ReadTrkExt Error : krep object not exist" << endl;
          }
          else {
            _kdiag.kalDiag(krep);
          }
          TrkExtTraj const *trkext = &(trkexts[i]);
          if ( !trkext) {
            throw cet::exception("DATA")
              <<"ReadTrkExt Error : trkexttraj object not exist" << endl;
          }
          else {
            _trkext.trkExtDiag(event, *krep, *trkext);
            
          }
        }

  
      }
      // if there are no tracks, enter dummies
      if(trks.size() == 0 && _processEmpty){
        if (_recordKalFit) {
          _kdiag.kalDiag(0);
          _trkext.trkExtDiag();
        }
      }
    }
  }


}  // end namespace mu2e

// Part of the magic that makes this class a module.
// create an instance of the module.  It also registers
using mu2e::ReadTrkExt;
DEFINE_ART_MODULE(ReadTrkExt);
