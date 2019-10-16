#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Provenance.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Utilities/InputTag.h"

#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/CalorimeterCalibrations.hh"

// #include "MCDataProducts/inc/SimParticleCollection.hh"
// #include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"

#include "CalorimeterGeom/inc/Calorimeter.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"

#include "RecoDataProducts/inc/CaloTrigSeedCollection.hh"
#include "RecoDataProducts/inc/TriggerInfo.hh"

#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"

// Root includes
#include "TDirectory.h"
#include "TTree.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

// C++ includes
#include <vector>

using namespace std;

namespace mu2e {

  
  class FilterEcalMVATrigger : public art::EDFilter {
  public:
    explicit FilterEcalMVATrigger(fhicl::ParameterSet const& pset);
    virtual ~FilterEcalMVATrigger() { }

    bool filter( art::Event& event);

    //virtual bool beginRun(art::Run &run) override;

    virtual void beginJob() override;
    virtual void endJob() override;

  private:

    int _diagLevel;
    std::string    _trigPath;

    std::string                _MVAMethodLabel;
    std::string _caloTrigSeedModuleLabel;
    std::string                _weightsfile;

    float _TOFF; // time offset to align fast clustering with tracker

    float _ENEMIN;
    float _DTMAX;

  
    float                      _MVAhighcut0;
    float                      _MVArpivot0;
    float                      _MVApivotcut0;
    float                      _MVAlowcut0;
    float                      _MVAhighcut1;
    float                      _MVArpivot1;
    float                      _MVApivotcut1;
    float                      _MVAlowcut1;
    int                        _downscale500_factor;
    int                        _step;
    int                        _nProcessed;

    static const int nECALDISKs=2;
    float _MVAcutA[nECALDISKs];
    float _MVAcutB[nECALDISKs];
    float _MVArpivot[nECALDISKs];
    float _MVAlowcut[nECALDISKs];


    TMVA::Reader *reader;

    std::string _MVAmethod;
    float _MVA;
    float xpeak,ypeak;
    float _rpeak,_tpeak;
    float _Epeak,_E10,_E11,_E20;
    float _fevt,_fdiskpeak; // needed by TMVA reader
    int disk;
    float MVAcut;
  };

  FilterEcalMVATrigger::FilterEcalMVATrigger(fhicl::ParameterSet const& pset):
    art::EDFilter{pset},
    _diagLevel                  (pset.get<int>("diagLevel",0)),
    _trigPath                   (pset.get<std::string>("triggerPath")),
    _MVAMethodLabel             (pset.get<std::string>("MVAMethod","BDT")), 
    _caloTrigSeedModuleLabel    (pset.get<std::string>("caloTrigSeedModuleLabel")),  
    _weightsfile                (pset.get<string>("weightsfile")),
    _TOFF                       (pset.get<float>("TimeOFFSET",22.5)),
    _MVAhighcut0                (pset.get<float>("MVAhighcut0",0.5)),
    _MVArpivot0                 (pset.get<float>("MVArpivot0",445.)),
    _MVApivotcut0               (pset.get<float>("MVApivotcut0",0.2)),
    _MVAlowcut0                 (pset.get<float>("MVAlowcut0",0.2)),
    _MVAhighcut1                (pset.get<float>("MVAhighcut1",0.5)),
    _MVArpivot1                 (pset.get<float>("MVArpivot1",445.)),
    _MVApivotcut1               (pset.get<float>("MVApivotcut1",0.2)),
    _MVAlowcut1                 (pset.get<float>("MVAlowcut1",0.2)),
    _downscale500_factor        (pset.get<float>("downscale500factor",0)),
    _step                       (pset.get<float>("step",10)),
    _nProcessed(0)  
  {
    _MVAmethod= _MVAMethodLabel + " method"; 
 
    _MVAcutB[0]=(_MVAhighcut0-_MVApivotcut0)/(395.-_MVArpivot0); // high cut is at r=395 mm
    _MVAcutA[0]=_MVApivotcut0-_MVAcutB[0]*_MVArpivot0;
    _MVArpivot[0]= _MVArpivot0;
    _MVAlowcut[0]= _MVAlowcut0;
    _MVAcutB[1]=(_MVAhighcut1-_MVApivotcut1)/(395.-_MVArpivot1); // high cut is at r=395 mm
    _MVAcutA[1]=_MVApivotcut1-_MVAcutB[1]*_MVArpivot1;
    _MVArpivot[1]= _MVArpivot1;
    _MVAlowcut[1]= _MVAlowcut1;
    
    produces<TriggerInfo>();
  }


  void FilterEcalMVATrigger::beginJob(){
    // art::ServiceHandle<art::TFileService> tfs;
    // art::TFileDirectory tfdir = tfs->mkdir("diag");
    
    ConfigFileLookupPolicy configFile;
    // this loads the TMVA library
    TMVA::Tools::Instance();
    reader = new TMVA::Reader( "!Color:Silent" );  
    reader->AddVariable("rpeak",&_rpeak);
    reader->AddVariable("tpeak",&_tpeak);
    reader->AddVariable("Epeak",&_Epeak);
    reader->AddVariable("E10",&_E10);
    reader->AddVariable("E11",&_E11);
    reader->AddVariable("E20",&_E20);
    reader->AddSpectator("fdiskpeak",&_fdiskpeak);
    reader->AddSpectator("fevt",&_fevt);  
    reader->BookMVA(_MVAmethod ,configFile(_weightsfile));
  }

  bool FilterEcalMVATrigger::filter(art::Event& event) {
    unique_ptr<TriggerInfo> triginfo(new TriggerInfo);
    bool   retval(false);

    if (_step==0) return false;

    _fevt=(float)event.id().event();

    ++_nProcessed;
    if (_nProcessed%10==0 && _diagLevel > 0) std::cout<<"Processing event from FilterEcalMVATrigger =  "<<_nProcessed <<std::endl;
   
    //Handle to the calorimeter
    art::ServiceHandle<GeometryService> geom;
    if( ! geom->hasElement<Calorimeter>() ) return false;
    Calorimeter const & cal = *(GeomHandle<Calorimeter>());
    // ------------------------- Unpack Calo Trig Seeds -------------------------------
    art::Handle<CaloTrigSeedCollection> caloTrigSeedsHandle;
    event.getByLabel(_caloTrigSeedModuleLabel, caloTrigSeedsHandle);
    CaloTrigSeedCollection const& caloTrigSeeds(*caloTrigSeedsHandle);
    if (_step==4) return false;

    for (CaloTrigSeedCollection::const_iterator seedIt = caloTrigSeeds.begin(); seedIt != caloTrigSeeds.end(); ++seedIt){
      disk= cal.crystal((int)seedIt->crystalid()).diskId();
      _fdiskpeak   = (float) disk;
      _Epeak   = seedIt->epeak();
      _tpeak   = seedIt->tpeak()+_TOFF;
      _rpeak   = seedIt->rpeak();
      _E10   = seedIt->ring1max();
      _E11   = seedIt->ring1max2();
      _E20   = seedIt->ring2max();
      //
      if (_step>5){
	_MVA= reader->EvaluateMVA(_MVAmethod);
	if (_diagLevel > 0){
	  std::cout<<"EVENT " << event.id().event() << " DISK " << disk << " Epeak=" << _Epeak << " tpeak=" << _tpeak << " E10=" << _E10 << " E11=" << _E11 << " E20=" << _E20 << " MVA=" << _MVA << " rpeak=" << _rpeak ;
	    
	  if (_rpeak>_MVArpivot[disk]){
	    std::cout << " cut at " << _MVAlowcut[disk] << std::endl;
	  }
	  else{
	    MVAcut=_MVAcutA[disk]+_MVAcutB[disk]*_rpeak;
	    std::cout << " cut at " << MVAcut <<  std::endl;
	  }
	}
	if (_rpeak>_MVArpivot[disk]){
	  if (_MVA>_MVAlowcut[disk]) {
	    retval = true;
	    triginfo->_triggerBits.merge(TriggerFlag::caloTrigSeed);
	    triginfo->_triggerPath = _trigPath;
	    size_t index = std::distance(caloTrigSeeds.begin(),seedIt);
	    triginfo->_caloTrigSeed = art::Ptr<CaloTrigSeed>(caloTrigSeedsHandle,index);
	    break;
	  }
	}
	else{
	  MVAcut=_MVAcutA[disk]+_MVAcutB[disk]*_rpeak;
	  if (_MVA>MVAcut) {
	    retval = true;
	    triginfo->_triggerBits.merge(TriggerFlag::caloTrigSeed);
	    triginfo->_triggerPath = _trigPath;
	    size_t index = std::distance(caloTrigSeeds.begin(),seedIt);
	    triginfo->_caloTrigSeed = art::Ptr<CaloTrigSeed>(caloTrigSeedsHandle,index);
	    break;
	  }
	}
      }
    }
    event.put(std::move(triginfo));
    return retval;
  }
  void FilterEcalMVATrigger::endJob(){
    cout << "FilterEcalMVATrigger filter end job:" << _nProcessed << " events processed" << endl;
  }
  
}

using mu2e::FilterEcalMVATrigger;
DEFINE_ART_MODULE(FilterEcalMVATrigger);
