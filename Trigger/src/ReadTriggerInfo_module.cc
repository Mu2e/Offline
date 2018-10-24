//
// An EDAnalyzer module that reads the Trigger Info 
//
// $Id:  $
// $Author:  $
// $Date:  $
//
// Original author G. Pezzullo
//

// #include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
// #include "GlobalConstantsService/inc/ParticleDataTable.hh"
// #include "GlobalConstantsService/inc/unknownPDGIdName.hh"
// #include "ConditionsService/inc/AcceleratorParams.hh"
// #include "ConditionsService/inc/ConditionsHandle.hh"

// #include "CalorimeterGeom/inc/Calorimeter.hh"
// #include "CalorimeterGeom/inc/DiskCalorimeter.hh"

// #include "GeometryService/inc/GeomHandle.hh"
// #include "GeometryService/inc/GeometryService.hh"
// #include "GeometryService/inc/VirtualDetector.hh"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Provenance.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
// #include "canvas/Utilities/InputTag.h"

#include "RecoDataProducts/inc/TriggerInfo.hh"
// #include "TDirectory.h"
// #include "TNtuple.h"
// #include "TTree.h"
// #include "TH2F.h"
#include "TH1F.h"

#include <cmath>
// #include <iostream>
#include <string>
// #include <map>
#include <vector>



namespace mu2e {


  class ReadTriggerInfo : public art::EDAnalyzer {

  public:

    explicit ReadTriggerInfo(fhicl::ParameterSet const& pset);
    virtual ~ReadTriggerInfo() { }

    virtual void beginJob();
    virtual void endJob();
    virtual void endSubRun(const art::SubRun& sr);

    // This is called for each event.
    virtual void analyze(const art::Event& e);
    
    void   findTrigIndex(std::vector<std::string> Vec, std::string ModuleLabel, int &Index);
  private:

    int      _diagLevel;
    size_t   _nMaxTrig;
    int      _nProcess;
    int      _numEvents;
    
    std::vector<std::string>  _trigLabels;
    std::vector<int>          _trigCounts;
    
    std::vector<std::string>  _trigLabelsCaloOnly;
    std::vector<int>          _trigCountsCaloOnly;
    
    std::vector<std::string>  _trigLabelsTrack;
    std::vector<int>          _trigCountsTrack;
    
    std::vector<std::string>  _trigLabelsEvtPS;
    std::vector<int>          _trigCountsEvtPS;
    
    TH1F *_hTrigInfo[4];

  };


  ReadTriggerInfo::ReadTriggerInfo(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset), 
    _diagLevel     (pset.get<int>("diagLevel", 0)),
    _nMaxTrig      (pset.get<size_t>("nFilters", 20)),
    _nProcess(0), 
    _numEvents(0)
  {
    _trigLabels.resize(_nMaxTrig);
    _trigCounts.resize(_nMaxTrig);

    _trigLabelsCaloOnly.resize(_nMaxTrig);
    _trigCountsCaloOnly.resize(_nMaxTrig);

    _trigLabelsTrack.resize(_nMaxTrig);
    _trigCountsTrack.resize(_nMaxTrig);

    _trigLabelsEvtPS.resize(_nMaxTrig);
    _trigCountsEvtPS.resize(_nMaxTrig);
  }
  
  void ReadTriggerInfo::beginJob(){

    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory trigInfoDir = tfs->mkdir("trigInfo");
    art::TFileDirectory caloTrigDir = tfs->mkdir("caloOnly");
    art::TFileDirectory trkTrigDir  = tfs->mkdir("track");

    _hTrigInfo[0] = trigInfoDir.make<TH1F>("hTrigInfo_global", "Global Trigger rejection", (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));       
    _hTrigInfo[1] = trigInfoDir.make<TH1F>("hTrigInfo_track" , "Calo-only Triggers rejection", (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));       
    _hTrigInfo[2] = trigInfoDir.make<TH1F>("hTrigInfo_calo"  , "Track Triggers rejection", (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));       
    _hTrigInfo[3] = trigInfoDir.make<TH1F>("hTrigInfo_evtPS" , "Event prescaler Trigger bits distribution", (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));       
    

  }

  void ReadTriggerInfo::endJob(){
    //set the labels of the histrograms
    for (size_t i=0; i<_trigCounts.size(); ++i ){
      _hTrigInfo[0]->GetXaxis()->SetBinLabel(i+1, _trigLabels[i].c_str());
      if (_trigCounts[i] > 0) _hTrigInfo[0]->SetBinContent(i+1, _nProcess/_trigCounts[i]);

      _hTrigInfo[1]->GetXaxis()->SetBinLabel(i+1, _trigLabelsTrack[i].c_str());
      if (_trigCountsTrack[i] > 0) _hTrigInfo[1]->SetBinContent(i+1, _nProcess/_trigCountsTrack[i]);

      _hTrigInfo[2]->GetXaxis()->SetBinLabel(i+1, _trigLabelsCaloOnly[i].c_str());
      if (_trigCountsCaloOnly[i] > 0) _hTrigInfo[2]->SetBinContent(i+1, _nProcess/_trigCountsCaloOnly[i]);

      _hTrigInfo[3]->GetXaxis()->SetBinLabel(i+1, _trigLabelsEvtPS[i].c_str());
      if (_trigCountsEvtPS[i] > 0) _hTrigInfo[3]->SetBinContent(i+1, _trigCountsEvtPS[i]);
    }
  }

  //================================================================
  void ReadTriggerInfo::endSubRun(const art::SubRun& sr){}

  void ReadTriggerInfo::analyze(const art::Event& event) {

    ++_nProcess;
    
    std::vector<art::Handle<TriggerInfo> > hTrigInfoVec;
    event.getManyByType(hTrigInfoVec);
    
    art::Handle<TriggerInfo>       hTrigInfo;
    TriggerFlag                    prescalerFlag   = TriggerFlag::prescaleRandom;
    TriggerFlag                    trackFlag       = TriggerFlag::track;
    TriggerFlag                    caloFlag        = TriggerFlag::caloCluster;
    TriggerFlag                    caloOrTrackFlag = trackFlag; caloOrTrackFlag.merge(caloFlag);
    
    for (size_t i=0; i<hTrigInfoVec.size(); ++i){
      hTrigInfo = hTrigInfoVec.at(i);
      if (!hTrigInfo.isValid())         continue;
      const TriggerInfo* trigInfo  = hTrigInfo.product();
      const TriggerFlag  flag      = trigInfo->triggerBits();

      std::string    moduleLabel         = hTrigInfo.provenance()->moduleLabel();
      int            index(0);         
      
      //fill the Global Trigger bits info
      findTrigIndex(_trigLabels, moduleLabel, index);
      _trigLabels[index] = moduleLabel;
      if ( flag.hasAnyProperty(caloOrTrackFlag)){ 
	_trigCounts[index] = _trigCounts[index] + 1;
      }
      
      //fill the Calo-Only Trigger bits info
      findTrigIndex(_trigLabelsCaloOnly, moduleLabel, index);
      if ( flag.hasAnyProperty(caloFlag)){ 
	_trigLabelsCaloOnly[index] = moduleLabel;
	_trigCountsCaloOnly[index] = _trigCountsCaloOnly[index] + 1;
      }
      
      //fill the Track Trigger bits info
      findTrigIndex(_trigLabelsTrack, moduleLabel, index);
      if ( flag.hasAnyProperty(trackFlag)){ 
	_trigLabelsTrack[index] = moduleLabel;
	_trigCountsTrack[index] = _trigCountsTrack[index] + 1;
      }
      
      //fill the Event-Prescaler Trigger bits info
      findTrigIndex(_trigLabelsEvtPS, moduleLabel, index);
      if ( flag.hasAnyProperty(prescalerFlag)){ 
	_trigLabelsEvtPS[index] = moduleLabel;
	_trigCountsEvtPS[index] = _trigCountsEvtPS[index] + 1;
      }
      
      
    }//end loop over the TriggerInfo Handles
  }
  
  void ReadTriggerInfo::findTrigIndex(std::vector<std::string> Vec, std::string ModuleLabel, int &Index){
    //reset the index value
    Index = 0;
    for (size_t i=0; i<Vec.size(); ++i){
      if (Vec[i] == ModuleLabel) { 
	Index = i;
	break;
      }else if (Vec[i] != ""){
	Index = i+1;
      }
    }
  }
  
}  

DEFINE_ART_MODULE(mu2e::ReadTriggerInfo);


