#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Optional/TFileService.h"
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

#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"

#include "CalorimeterGeom/inc/Calorimeter.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"

#include "RecoDataProducts/inc/CaloClusterCollection.hh"

#include "RecoDataProducts/inc/CaloDigi.hh"
#include "RecoDataProducts/inc/CaloDigiCollection.hh"

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

  struct EcalTrigDigi{
    int cryId;
    float x;
    float y;
    float phi;
    float tpeak;
    float amp;
  };

  struct strawHit{
    float x;
    float y;
    float z;
    float t;
    float dadz;
    float dfdz;
    float dt;
  };
  
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

    std::string                _MVAMethodLabel;
    std::string _caloDigiModuleLabel;
    std::string                _weightsfile;

    SimParticleTimeOffset _toff;  // time offset smearing

    float _ENEMIN;
    float _DTMAX;

    static const int nECALDISKs=2;
    std::multimap<float,EcalTrigDigi> _DigiMapE[nECALDISKs];
    std::multimap<float,EcalTrigDigi> _Ring1MapE;
    std::multimap<float,EcalTrigDigi> _Ring2MapE;
  

    const int                  _caloCrystals = 1356;
    float                      _ADC2MeV;
    float                      _PEAK2E;
    float                      _DIGITHRESHOLD;
    float                      _T0MIN;
    float                      _PEAKmin;
    float                      _DTmax;
    float                      _RINGRMIN0;
    float                      _RINGRMIN1;
    float                      _RINGRMIN2;
    float                      _wcry;
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

    float _MVAcutA[nECALDISKs];
    float _MVAcutB[nECALDISKs];
    float _MVArpivot[nECALDISKs];
    float _MVAlowcut[nECALDISKs];

    // digitizations
    std::vector < int > _DigiIds;
    std::vector < double > _DigiT0s;
    std::vector < std::vector < int > > _DigiWaves;

    map <int,double> _wave;
    vector <int> _cryIDs;
    vector <double> _TPEAKs;
    vector <double> _EPEAKs;

    TMVA::Reader *reader;

    std::string _MVAmethod;
    float _MVA;
    float _rpeak,_tpeak;
    float _Epeak,_E10,_E11,_E20;
    float _fevt,_fdiskpeak; // needed by TMVA reader
  };

  FilterEcalMVATrigger::FilterEcalMVATrigger(fhicl::ParameterSet const& pset):
    _diagLevel(pset.get<int>("diagLevel",0)),

    _MVAMethodLabel(pset.get<std::string>("MVAMethod","BDT")), 
    _caloDigiModuleLabel      (pset.get<std::string>("caloDigiModuleLabel","CaloDigiFromShower")), 
    _weightsfile               (pset.get<string>("weightsfile")),
    _toff(pset.get<fhicl::ParameterSet>("TimeOffsets", fhicl::ParameterSet())),
    _ADC2MeV                   (pset.get<float>("ADC2MeV",0.0076)),
    _PEAK2E                    (pset.get<float>("PEAK2E",6.04)),
    _DIGITHRESHOLD             (pset.get<float>("DigiThreshold",1.)),
    _T0MIN                     (pset.get<float>("T0min",500.)),
    _PEAKmin                   (pset.get<float>("PEAKmin",20.)),
    _DTmax                     (pset.get<float>("DTmax",6.)),
    _RINGRMIN0                 (pset.get<float>("RINGRmin0",0.)),
    _RINGRMIN1                 (pset.get<float>("RINGRmin1",0.)),
    _RINGRMIN2                 (pset.get<float>("RINGRmin2",0.)),
    _wcry                      (pset.get<float>("wcry",34.3)),
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
  }


  void FilterEcalMVATrigger::beginJob(){
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory tfdir = tfs->mkdir("diag");
    
    ConfigFileLookupPolicy configFile;
    // this loads the TMVA library
    TMVA::Tools::Instance();
    reader = new TMVA::Reader( "!Color:!Silent" );  
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

    if (_step==0) return false;

    _fevt=(float)event.id().event();

    ++_nProcessed;
    if (_nProcessed%10==0 && _diagLevel > 0) std::cout<<"Processing event from FilterEcalMVATrigger =  "<<_nProcessed <<std::endl;

    //Handle to the calorimeter
    art::ServiceHandle<GeometryService> geom;
    if( ! geom->hasElement<Calorimeter>() ) return false;
    Calorimeter const & cal = *(GeomHandle<Calorimeter>());

    // Microbunch time
    ConditionsHandle<AcceleratorParams> accPar("ignored");
    _toff.updateMap(event);
    
    //--------------------------  Do Calo Digis  --------------------------------    
    art::Handle<CaloDigiCollection> caloDigisHandle;
    event.getByLabel(_caloDigiModuleLabel, caloDigisHandle);
    //auto caloDigiFlag = event.getValidHandle<mu2e::CaloDigiCollection>(_caloDigiModuleLabel);
    //if (*caloDigiFlag.product() == 0){
    if (! caloDigisHandle.isValid() ){
      std::cout << "FILTERECALMVATRIGGER: CaloDigiHandle NOT VALID" << std::endl;
      return false; return false;
    }
    //    const CaloDigiCollection& caloDigis = *caloDigiFlag.product();
    const CaloDigiCollection& caloDigis(*caloDigisHandle);
    _DigiIds.clear();
    _DigiT0s.clear();
    _DigiWaves.clear();
    for (const auto& caloDigi : caloDigis){
      _DigiIds.push_back(caloDigi.roId());
      _DigiT0s.push_back(caloDigi.t0());
      _DigiWaves.push_back(caloDigi.waveform());
    }


    //******************************************************************************
    // Process Calo Digi
    int cryID;
    int index;
    float tnow,tmin;
    _cryIDs.clear();
    _TPEAKs.clear();
    _EPEAKs.clear();
    index=0;
    for (std::vector<vector<int> >::iterator it=_DigiWaves.begin(); it!=_DigiWaves.end(); ++it){
      if (_DigiIds.at(index)%2==1){
	index++;
	continue; // consider only even sensors
      }
      cryID=_DigiIds.at(index)/2;
      tmin=(float) _DigiT0s.at(index);
      for (unsigned int i=2;i<it->size()-2;++i){
	if ( it->at(i)*_ADC2MeV*_PEAK2E > _DIGITHRESHOLD){
	  tnow=tmin+(i+0.5)*5.; // 5 ns tick 
	  if ( tnow > _T0MIN - _DTmax) {
	    if ( it->at(i) > it->at(i-1) &&
		 it->at(i) > it->at(i-2) &&
		 it->at(i) > it->at(i+1) &&
		 it->at(i) > it->at(i+2)){
	      _cryIDs.push_back(cryID); 
	      _TPEAKs.push_back((double)tnow);
	      _EPEAKs.push_back((double)it->at(i)*_ADC2MeV*_PEAK2E);
	    }
	  }
	}
      }
      index++;
    }
  
    if (_step==1) return false;
    // Iterate on reconstructed digi
    float tpeak;
    float amp;
    int disk;
    CLHEP::Hep3Vector  crystalPos(0);
    EcalTrigDigi _EcalTrigDigi;
    for (int idisk=0;idisk<nECALDISKs;idisk++){
      _DigiMapE[idisk].clear();
    }
    index=0;
    for (std::vector<int>::iterator it=_cryIDs.begin(); it!=_cryIDs.end(); ++it){
      tpeak= _TPEAKs.at(index);
      if (tpeak<_T0MIN){
	index++;
	continue;
      }
      cryID=*it;
      amp= (float) _EPEAKs.at(index);
      //
      _EcalTrigDigi.cryId= cryID;
      crystalPos = cal.crystal(cryID).localPosition();
      _EcalTrigDigi.x=crystalPos.x();
      _EcalTrigDigi.y=crystalPos.y();
      _EcalTrigDigi.phi = atan2(_EcalTrigDigi.y,_EcalTrigDigi.x);
      _EcalTrigDigi.tpeak = tpeak; 
      _EcalTrigDigi.amp = amp; 
      disk = cal.crystal(cryID).diskId();
      _DigiMapE[disk].insert(std::make_pair(amp,_EcalTrigDigi));
      index++;
    }
    if (_step==2) return false;
    float xpeak,ypeak;
    float x,y,r;
    float dist;
    float MVAcut;
    // Loop on disks
    for (int idisk=0;idisk<nECALDISKs;idisk++){
      _fdiskpeak=(float) idisk;
      // Loop on digi map for this disk (amp decreasing ordered)
      for (std::multimap<float,EcalTrigDigi>::reverse_iterator rit=_DigiMapE[idisk].rbegin(); rit!=_DigiMapE[idisk].rend(); ++rit){
	if ( rit->second.amp < _PEAKmin) break; // peak amplitude threshold
	xpeak=rit->second.x;
	ypeak=rit->second.y;
	_rpeak=sqrt(xpeak*xpeak+ypeak*ypeak);
	if (_rpeak<_RINGRMIN0) continue; // skip crystals too close to beam
	_tpeak=rit->second.tpeak;
	_Epeak=rit->first;

	// Search crystals in close rings
	_Ring1MapE.clear();
	_Ring2MapE.clear();
	for (std::multimap<float,EcalTrigDigi>::reverse_iterator rit2=_DigiMapE[idisk].rbegin(); rit2!=_DigiMapE[idisk].rend(); ++rit2){
	  //
	  x=rit2->second.x;
	  y=rit2->second.y;
	  //
	  dist=sqrt((x-xpeak)*(x-xpeak)+(y-ypeak)*(y-ypeak));
	  if (dist>2.5*_wcry) continue; // too far
	  if (dist<0.5*_wcry) continue; // same as peak
	  //
	  if (abs(rit2->second.tpeak-_tpeak)>_DTmax) continue; // too early or too late
	  r=sqrt(x*x+y*y);
	  if (dist<1.5*_wcry){ // first ring
	    // cut on minimum r
	    if (r>_RINGRMIN1) _Ring1MapE.insert(std::make_pair(rit2->second.amp,rit2->second));
	  }
	  else{
	    if (r>_RINGRMIN2) _Ring2MapE.insert(std::make_pair(rit2->second.amp,rit2->second));
	  }
	}
	if (_Ring1MapE.size()>0){
	  std::multimap<float,EcalTrigDigi>::reverse_iterator ring1it=_Ring1MapE.rbegin();
	  _E10=ring1it->first;
	  if (_Ring1MapE.size()>1){
	    std::advance(ring1it,1);
	    _E11=ring1it->first;
	  }
	  else _E11=0.;
	}
	else{
	  _E10=0.;
	  _E11=0.;
	}
	if (_Ring2MapE.size()>0){
	  _E20=_Ring2MapE.rbegin()->first; 
	}
	else _E20=0.;
	//
	if (_step>3){
	  _MVA= reader->EvaluateMVA(_MVAmethod);
	  if (_diagLevel > 0){
	    std::cout<<"EVENT " << event.id().event() << " DISK " << idisk << " Epeak=" << _Epeak << " tpeak=" << _tpeak << " E10=" << _E10 << " E11=" << _E11 << " E20=" << _E20 << " MVA=" << _MVA << " rpeak=" << _rpeak ;
	  
	    if (_rpeak>_MVArpivot[idisk]){
	      std::cout << " cut at " << _MVAlowcut[idisk] << std::endl;
	    }
	    else{
	      MVAcut=_MVAcutA[idisk]+_MVAcutB[idisk]*_rpeak;
	      std::cout << " cut at " << MVAcut <<  std::endl;
	    }
	  }
	  if (_rpeak>_MVArpivot[idisk]){
	    if (_MVA>_MVAlowcut[idisk]) return true;
	  }
	  else{
	    MVAcut=_MVAcutA[idisk]+_MVAcutB[idisk]*_rpeak;
	    if (_MVA>MVAcut) return true;
	  }
	}
      }
    }

    return false;
    
  }
  void FilterEcalMVATrigger::endJob(){
    cout << "FilterEcalMVATrigger filter end job:" << _nProcessed << " events processed" << endl;
  }
  
}

using mu2e::FilterEcalMVATrigger;
DEFINE_ART_MODULE(FilterEcalMVATrigger);
