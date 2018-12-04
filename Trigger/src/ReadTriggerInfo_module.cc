//
// An EDAnalyzer module that reads the Trigger Info 
//
// $Id:  $
// $Author:  $
// $Date:  $
//
// Original author G. Pezzullo
//

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
#include "BFieldGeom/inc/BFieldManager.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"

//Conditions
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"

//Dataproducts
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloTrigSeed.hh"
#include "RecoDataProducts/inc/HelixSeed.hh"
#include "RecoDataProducts/inc/KalSeed.hh"
#include "RecoDataProducts/inc/TriggerInfo.hh"

//ROOT
#include "TH1F.h"
#include "TH2F.h"

#include <cmath>
// #include <iostream>
#include <string>
// #include <map>
#include <vector>



namespace mu2e {


  class ReadTriggerInfo : public art::EDAnalyzer {

  public:

    enum {
      kNTrigInfo     = 20,
      kNTrackTrig    = 10,
      kNTrackTrigVar = 6,
      kNHelixTrig    = 10,
      kNHelixTrigVar = 7,
      kNCaloCalib    = 5,
      kNCaloCalibVar = 5,
      kNCaloOnly     = 5,
      kNCaloOnlyVar  = 5      
    };

    struct  trigInfo_ {
      int           counts;
      int           exclusive_counts;
      std::string   label;
      
      trigInfo_ ():counts(0), exclusive_counts(0){}
    };

    
    explicit ReadTriggerInfo(fhicl::ParameterSet const& pset);
    virtual ~ReadTriggerInfo() { }

    virtual void beginJob();
    virtual void endJob();
    virtual void endSubRun(const art::SubRun& sr);

    // This is called for each event.
    virtual void analyze(const art::Event& e);
    virtual void beginRun(const art::Run & run);

    void   findTrigIndex        (std::vector<trigInfo_> Vec, std::string ModuleLabel, int &Index);
    void   fillTrackTrigInfo    (int TrkTrigIndex, const KalSeed*  KSeed);
    void   fillHelixTrigInfo    (int HelTrigIndex, const HelixSeed*HSeed);
    void   fillCaloTrigSeedInfo (int CTrigSeedIndex, const CaloTrigSeed*HCl);
    void   fillCaloCalibTrigInfo(int ClCalibIndex, const CaloCluster*  HCl);

    void   findCorrelatedEvents (std::vector<string>& VecLabels, double &NCorrelated);
    void   evalTriggerRate      ();

  private:

    int      _diagLevel;
    size_t   _nMaxTrig;
    
    int      _nTrackTrig;
    int      _nCaloTrig;
    int      _nCaloCalibTrig;

    double   _duty_cycle;

    int      _nProcess;
    int      _numEvents;
    
    double   _bz0;
    
    std::vector<trigInfo_>    _trigAll;	     
    std::vector<trigInfo_>    _trigFinal;    
    std::vector<trigInfo_>    _trigCaloOnly; 
    std::vector<trigInfo_>    _trigCaloCalib;
    std::vector<trigInfo_>    _trigTrack;    
    std::vector<trigInfo_>    _trigHelix;    
    std::vector<trigInfo_>    _trigEvtPS;    
    
    TH1F *_hTrigInfo  [kNTrigInfo];
    TH2F *_h2DTrigInfo[kNTrigInfo];

    TH1F *_hTrkInfo [kNTrackTrig][kNTrackTrigVar];
    TH1F *_hHelInfo [kNHelixTrig][kNHelixTrigVar];

    TH1F *_hCaloCalibInfo[kNCaloCalib][kNCaloCalibVar];
    TH1F *_hCaloOnlyInfo [kNCaloOnly][kNCaloOnlyVar];
    
    TH1F *_hTrigBDW[kNTrigInfo];
  };


  ReadTriggerInfo::ReadTriggerInfo(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset), 
    _diagLevel     (pset.get<int>   ("diagLevel", 0)),
    _nMaxTrig      (pset.get<size_t>("nFilters", 50)),
    _nTrackTrig    (pset.get<size_t>("nTrackTriggers", 4)),
    _nCaloTrig     (pset.get<size_t>("nCaloTriggers", 4)),
    _nCaloCalibTrig(pset.get<size_t>("nCaloCalibTriggers", 4)),
    _duty_cycle    (pset.get<float> ("dutyCycle", 1.)),
    _nProcess(0), 
    _numEvents(0)
  {
    _trigAll.      resize(_nMaxTrig);	     
    _trigFinal.    resize(_nMaxTrig);    
    _trigCaloOnly. resize(_nMaxTrig); 
    _trigCaloCalib.resize(_nMaxTrig);
    _trigTrack.    resize(_nMaxTrig);    
    _trigHelix.    resize(_nMaxTrig);    
    _trigEvtPS.    resize(_nMaxTrig);    
  }
  
  void ReadTriggerInfo::beginJob(){

    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory trigInfoDir = tfs->mkdir("trigInfo");

    _hTrigInfo[0]   = trigInfoDir.make<TH1F>("hTrigInfo_global"    , "Global Trigger rejection"                   , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));       
    _hTrigInfo[1]   = trigInfoDir.make<TH1F>("hTrigInfo_track"     , "Calo-only Triggers rejection"               , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));       
    _hTrigInfo[2]   = trigInfoDir.make<TH1F>("hTrigInfo_calo"      , "Track Triggers rejection"                   , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));       
    _hTrigInfo[3]   = trigInfoDir.make<TH1F>("hTrigInfo_evtPS"     , "Event prescaler Trigger bits distribution"  , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));       
    _hTrigInfo[4]   = trigInfoDir.make<TH1F>("hTrigInfo_helix"     , "HelixSeed Triggers rejection"               , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));       
    _hTrigInfo[5]   = trigInfoDir.make<TH1F>("hTrigInfo_caloCalib" , "Calo Calibration rejection"                 , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));       
    _hTrigInfo[6]   = trigInfoDir.make<TH1F>("hTrigInfo_final"     , "Global Trigger rejection of the paths"      , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));       

    _hTrigInfo[10]  = trigInfoDir.make<TH1F>("hTrigInfo_unique_all", "Events found only by each Trig path"        , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));       
    _hTrigInfo[11]  = trigInfoDir.make<TH1F>("hTrigInfo_unique"    , "Events found only by each Trig path"        , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));       

    _h2DTrigInfo[0] = trigInfoDir.make<TH2F>("h2DTrigInfo_map_all" , "Trigger correlation map from all filters"   , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5), (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));       
    _h2DTrigInfo[1] = trigInfoDir.make<TH2F>("h2DTrigInfo_map"     , "Trigger correlation map"                    , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5), (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));   

    for (int i=0; i<_nCaloTrig; ++i){
      art::TFileDirectory caloInfoDir = tfs->mkdir(Form("caloOnly_%i",i));
      _hCaloOnlyInfo[i][0] = caloInfoDir.make<TH1F>(Form("hEPeak_%i"   , i), "peak energy; E[MeV]"        , 400, 0, 200);
      _hCaloOnlyInfo[i][1] = caloInfoDir.make<TH1F>(Form("hR1Max1_%i"   , i), "ring1 max; ring1max [MeV]" , 400, 0, 200);
      _hCaloOnlyInfo[i][2] = caloInfoDir.make<TH1F>(Form("hR1Max2_%i"   , i), "ring1 max; ring1max2 [MeV]", 400, 0, 200);
    
    }
    
    for (int i=0; i<_nCaloCalibTrig; ++i){
      art::TFileDirectory caloCalibInfoDir = tfs->mkdir(Form("caloCalib_%i",i));
      _hCaloCalibInfo[i][0] = caloCalibInfoDir.make<TH1F>(Form("hE_%i"   , i), "Cluster energy; E[MeV]", 800, 0, 800);
      _hCaloCalibInfo[i][1] = caloCalibInfoDir.make<TH1F>(Form("hN_%i"   , i), "Cluster size; nCrystalHits", 101, -0.5, 100.5);
   
    }
    
    for (int i=0; i<_nTrackTrig; ++i){
      art::TFileDirectory trkInfoDir  = tfs->mkdir(Form("trk_%i",i));
      _hTrkInfo[i][0] = trkInfoDir.make<TH1F>(Form("hP_%i"    , i), "Track Momentum; p[MeV/c]", 400, 0, 200);
      _hTrkInfo[i][1] = trkInfoDir.make<TH1F>(Form("hPt_%i"   , i), "Track Pt; p_{t} [MeV/c]", 400, 0, 200);
      _hTrkInfo[i][2] = trkInfoDir.make<TH1F>(Form("hNSh_%i"  , i), "N-StrawHits; nStrawHits", 101, -0.5, 100.5);
      _hTrkInfo[i][3] = trkInfoDir.make<TH1F>(Form("hD0_%i"   , i), "Track impact parameter; d0 [mm]", 801, -400.5, 400.5);
      _hTrkInfo[i][4] = trkInfoDir.make<TH1F>(Form("hChi2d_%i", i), "Track #chi^{2}/ndof;#chi^{2}/ndof", 100, 0, 50);
      _hTrkInfo[i][5] = trkInfoDir.make<TH1F>(Form("hClE_%i"  , i), "calorimeter Cluster energy; E [MeV]", 240, 0, 120);
    }
    
    for (int i=0; i<_nTrackTrig; ++i){
      art::TFileDirectory helinfoDir  = tfs->mkdir(Form("helix_%i",i));
      _hHelInfo[i][0] = helinfoDir.make<TH1F>(Form("hP_%i"    , i), "Helix Momentum; p[MeV/c]", 400, 0, 200);
      _hHelInfo[i][1] = helinfoDir.make<TH1F>(Form("hPt_%i"   , i), "Helix Pt; p_{t} [MeV/c]", 400, 0, 200);
      _hHelInfo[i][2] = helinfoDir.make<TH1F>(Form("hNSh_%i"  , i), "N-StrawHits; nStrawHits", 101, -0.5, 100.5);
      _hHelInfo[i][3] = helinfoDir.make<TH1F>(Form("hD0_%i"   , i), "Helix impact parameter; d0 [mm]", 801, -400.5, 400.5);
      _hHelInfo[i][4] = helinfoDir.make<TH1F>(Form("hChi2dXY_%i"  , i), "Helix #chi^{2}_{xy}/ndof;#chi^{2}_{xy}/ndof"      , 100, 0, 50);
      _hHelInfo[i][5] = helinfoDir.make<TH1F>(Form("hChi2dZPhi_%i", i), "Helix #chi^{2}_{z#phi}/ndof;#chi^{2}_{z#phi}/ndof", 100, 0, 50);
      _hHelInfo[i][6] = helinfoDir.make<TH1F>(Form("hClE_%i"  , i), "calorimeter Cluster energy; E [MeV]", 240, 0, 120);
    }
    
    art::TFileDirectory trigBDWDir = tfs->mkdir("trigBDW");

    _hTrigBDW[0]   = trigBDWDir.make<TH1F>("hTrigBDW_global"    , "Trigger bandwidth; ; rate [Hz]"                   , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));       
    _hTrigBDW[1]   = trigBDWDir.make<TH1F>("hTrigBDW_cumulative", "Cumulative Trigger bandwidth; ; rate [Hz]"        , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));       
    
  }

  void ReadTriggerInfo::endJob(){

    int    indexTrigInfo11(0);

    //fill the histograms
    for (size_t i=0; i<_trigAll.size(); ++i ){
      _hTrigInfo  [0]->GetXaxis()->SetBinLabel(i+1, _trigAll[i].label.c_str());
      _h2DTrigInfo[0]->GetXaxis()->SetBinLabel(i+1, _trigAll[i].label.c_str());

      if (_trigAll[i].counts > 0) {
	_hTrigInfo[0]->SetBinContent(i+1, _nProcess/_trigAll[i].counts);
	for (size_t j=0; j<_trigAll.size(); ++j ){
	  _h2DTrigInfo[0]->GetYaxis()->SetBinLabel(j+1, _trigAll[j].label.c_str());
	}
      }

      _hTrigInfo[1]->GetXaxis()->SetBinLabel(i+1, _trigTrack[i].label.c_str());
      if (_trigTrack[i].counts > 0) _hTrigInfo[1]->SetBinContent(i+1, _nProcess/_trigTrack[i].counts);

      _hTrigInfo[2]->GetXaxis()->SetBinLabel(i+1, _trigCaloOnly[i].label.c_str());
      if (_trigCaloOnly[i].counts > 0) _hTrigInfo[2]->SetBinContent(i+1, _nProcess/_trigCaloOnly[i].counts);

      _hTrigInfo[3]->GetXaxis()->SetBinLabel(i+1, _trigEvtPS[i].label.c_str());
      if (_trigEvtPS[i].counts > 0) _hTrigInfo[3]->SetBinContent(i+1, _trigEvtPS[i].counts);

      _hTrigInfo[4]->GetXaxis()->SetBinLabel(i+1, _trigHelix[i].label.c_str());
      if (_trigHelix[i].counts > 0) _hTrigInfo[4]->SetBinContent(i+1, _trigHelix[i].counts);

      _hTrigInfo[5]->GetXaxis()->SetBinLabel(i+1, _trigCaloCalib[i].label.c_str());
      if (_trigCaloCalib[i].counts > 0) _hTrigInfo[5]->SetBinContent(i+1, _trigCaloCalib[i].counts);

      if (_trigFinal[i].counts > 0) {
	_hTrigInfo  [6]->GetXaxis()->SetBinLabel(i+1, _trigFinal[i].label.c_str());
	_hTrigInfo  [6]->SetBinContent(i+1, _nProcess/_trigFinal[i].counts);
      }

      //fill  the histograms that shows how many events were found exclusively by each trigger path
      _hTrigInfo  [10]->GetXaxis()->SetBinLabel(i+1, _trigAll[i].label.c_str());
      double    content_trigInfo11 = _hTrigInfo [10]->GetBinContent(i+1);
      if (content_trigInfo11>0){
	_hTrigInfo  [11]->GetXaxis()->SetBinLabel(indexTrigInfo11 +1, _trigAll[i].label.c_str());
	_hTrigInfo  [11]->SetBinContent(indexTrigInfo11 +1, content_trigInfo11);
	++indexTrigInfo11;
      }

    }

    //now let's filter the 2D correlation histogram with only those that actually triggered at least one event
    int                nbinsx = _h2DTrigInfo[0]->GetNbinsX();
    int                nbinsy = _h2DTrigInfo[0]->GetNbinsY();
    std::vector<int>   binsToSkip;

    for (int i=0; i<nbinsx; ++i){
      bool used(false);

      for (int j=0; j<nbinsy; ++j){
	if (_h2DTrigInfo[0]->GetBinContent(i+1, j+1) > 0) {
	  used = true;
	  break;
	}
      }
      if (!used) binsToSkip.push_back(i);
    }

    int   index_x(0);
    for (int i=0; i<nbinsx; ++i){
      int    counts = std::count(binsToSkip.begin(), binsToSkip.end(), i);
      if (counts >= 1)       continue;
      //set the label
      _h2DTrigInfo[1]->GetXaxis()->SetBinLabel(index_x+1, _trigAll[i].label.c_str());

      int    index_y(0);

      for (int j=0; j<nbinsy; ++j){
	counts = std::count(binsToSkip.begin(), binsToSkip.end(), j);
	if (counts >= 1)       continue;
	double  content =  _h2DTrigInfo[0]->GetBinContent(i+1, j+1);
	_h2DTrigInfo[1]->SetBinContent(index_x+1, index_y+1, content);
	
	//set the label
	if (index_x == 0){
	  _h2DTrigInfo[1]->GetYaxis()->SetBinLabel(index_y+1, _trigAll[j].label.c_str());
	}

	++index_y;
      }
      ++index_x;
    }
    
    // now evaluate the bandwidth
    // NOTE: "evalTriggerrate" re-order the vectors _trigFinal
    evalTriggerRate();
  }
  
  //--------------------------------------------------------------------------------  
  void   ReadTriggerInfo::evalTriggerRate        (){
    //order the array with the filter used at the end of each path
    std::sort(_trigFinal.begin(), _trigFinal.end(), [](const auto a, const auto b) {return a.counts < b.counts; });
    
    ConditionsHandle<AcceleratorParams> accPar("ignored");
    double    mbtime         = accPar->deBuncherPeriod;
    double    mean_mb_rate   = 1./(mbtime/CLHEP::s)*_duty_cycle;

    bool      isFirst(true);
    int       index(0);
   
    std::vector<string>    labels_by_rate;

    for (size_t i=0; i< _trigFinal.size(); ++i){
      double  nEvents = (double)_trigFinal[i].counts;
      if ( nEvents <= 1e-3)                 continue;

      labels_by_rate.push_back(_trigFinal[i].label);

      double  eff   = nEvents/(double)_nProcess;
      double  rate  = mean_mb_rate*eff;
      _hTrigBDW[0]->GetXaxis()->SetBinLabel(index+1, _trigFinal[i].label.c_str());
      _hTrigBDW[1]->GetXaxis()->SetBinLabel(index+1, _trigFinal[i].label.c_str());
      _hTrigBDW[0]->SetBinContent(index+1, rate);
      
      if (isFirst) {
      	_hTrigBDW[1]->SetBinContent(index+1, rate);
      	//	rate_ref  = rate;
      	isFirst = false;
      }else{
	double    nCorrelated(0);
	findCorrelatedEvents(labels_by_rate, nCorrelated);

	rate = _hTrigBDW[1]->GetBinContent(index) + (nEvents-nCorrelated)/(double)_nProcess*mean_mb_rate;
      	_hTrigBDW[1]->SetBinContent(index+1, rate);
      }

      ++index;
    }

  }



  void   ReadTriggerInfo::findCorrelatedEvents(std::vector<string>& VecLabels, double &NCorrelated){
    
    NCorrelated = 0;
    
    const char* label_ref = VecLabels.at(VecLabels.size()-1).c_str();

    //    char* label(0);
    
    int        nLabels = VecLabels.size() -1;
    int        nbins   = _h2DTrigInfo[1]->GetNbinsX();
    for(int i=0; i<nbins; ++i){
      const char* label = _h2DTrigInfo[1]->GetXaxis()->GetBinLabel(i+1);
      if (std::strcmp(label_ref, label) != 0)      continue;
      
      for (int k=0; k<nLabels; ++k){
	label_ref = VecLabels.at(k).c_str();
	for (int j=0; j<nbins; ++j){
	  if (j == i)      break;
	  label =   _h2DTrigInfo[1]->GetYaxis()->GetBinLabel(j+1);
	  if (std::strcmp(label_ref, label) != 0)        continue;
	  NCorrelated += _h2DTrigInfo[1]->GetBinContent(i+1, j+1);
	}
      }
      break;
    }
    
    
    
  }

  //================================================================
  void   ReadTriggerInfo::beginRun(const art::Run & run){
    // get bfield
    GeomHandle<BFieldManager> bfmgr;
    GeomHandle<DetectorSystem> det;
    CLHEP::Hep3Vector vpoint_mu2e = det->toMu2e(CLHEP::Hep3Vector(0.0,0.0,0.0));
    _bz0 = bfmgr->getBField(vpoint_mu2e).z();
  }

  void ReadTriggerInfo::endSubRun(const art::SubRun& sr){}

  void ReadTriggerInfo::analyze(const art::Event& event) {

    ++_nProcess;
    
    std::vector<art::Handle<TriggerInfo> > hTrigInfoVec;
    event.getManyByType(hTrigInfoVec);
    
    art::Handle<TriggerInfo>       hTrigInfo;
    TriggerFlag                    prescalerFlag       = TriggerFlag::prescaleRandom;
    TriggerFlag                    trackFlag           = TriggerFlag::track;
    TriggerFlag                    helixFlag           = TriggerFlag::helix;
    TriggerFlag                    caloFlag            = TriggerFlag::caloCluster;
    TriggerFlag                    caloCalibFlag       = TriggerFlag::caloCalib;
    TriggerFlag                    caloTrigSeedFlag    = TriggerFlag::caloTrigSeed;
    TriggerFlag                    caloOrTrackFlag     = trackFlag; caloOrTrackFlag.merge(caloFlag); caloOrTrackFlag.merge(caloCalibFlag); caloOrTrackFlag.merge(caloTrigSeedFlag);// caloOrTrackFlag.merge(helixFlag);
    
    std::vector<int>   trigFlagAll_index, trigFlag_index;
    
    for (size_t i=0; i<hTrigInfoVec.size(); ++i){
      hTrigInfo = hTrigInfoVec.at(i);
      if (!hTrigInfo.isValid())         continue;
      const TriggerInfo* trigInfo  = hTrigInfo.product();
      const TriggerFlag  flag      = trigInfo->triggerBits();

      std::string    moduleLabel   = hTrigInfo.provenance()->moduleLabel();
      int            index_all(0);         
      int            index(0);         
      
      //fill the Global Trigger bits info
      findTrigIndex(_trigAll, moduleLabel, index_all);
      _trigAll[index_all].label  = moduleLabel;
      //      if ( flag.hasAnyProperty(caloOrTrackFlag)){ 
	//      }

      findTrigIndex(_trigFinal, moduleLabel, index);
      if ( flag.hasAnyProperty(caloOrTrackFlag)){ 
	_trigFinal[index].label  = moduleLabel;
	_trigFinal[index].counts = _trigFinal[index].counts + 1;
	_trigAll[index_all].counts = _trigAll[index_all].counts + 1;
	trigFlagAll_index.push_back(index_all);
      }
      //fill the Calo-Only Trigger bits info
      findTrigIndex(_trigCaloOnly, moduleLabel, index);
      if ( flag.hasAnyProperty(caloFlag) || flag.hasAnyProperty(caloTrigSeedFlag)){ 
	_trigCaloOnly[index].label  = moduleLabel;
	_trigCaloOnly[index].counts = _trigCaloOnly[index].counts + 1;
	const CaloTrigSeed*clseed = trigInfo->caloTrigSeed().get();
	if(clseed) fillCaloTrigSeedInfo(index, clseed);
	trigFlag_index.push_back(index_all);
      }

      //fill the CaloCalib Trigger bits info
      findTrigIndex(_trigCaloCalib, moduleLabel, index);
      if ( flag.hasAnyProperty(caloCalibFlag)){ 
	_trigCaloCalib[index].label  = moduleLabel;
	_trigCaloCalib[index].counts = _trigCaloCalib[index].counts + 1;
	const CaloCluster*cluster = trigInfo->caloCluster().get();
	if(cluster) fillCaloCalibTrigInfo(index, cluster);
	trigFlag_index.push_back(index_all);
      }
      
      //fill the Track Trigger bits info
      findTrigIndex(_trigTrack, moduleLabel, index);
      if ( flag.hasAnyProperty(trackFlag)){ 
	_trigTrack[index].label  = moduleLabel;
	_trigTrack[index].counts = _trigTrack[index].counts + 1;
	const KalSeed*kseed = trigInfo->track().get();
	if(kseed) fillTrackTrigInfo(index, kseed);
	trigFlag_index.push_back(index_all);
      }
       //fill the Helix Trigger bits info
      findTrigIndex(_trigHelix, moduleLabel, index);
      if ( flag.hasAnyProperty(helixFlag)){ 
	_trigHelix[index].label  = moduleLabel;
	_trigHelix[index].counts = _trigHelix[index].counts + 1;
	const HelixSeed*hseed = trigInfo->helix().get();
	if(hseed) fillHelixTrigInfo(index, hseed);
      }
      
      //fill the Event-Prescaler Trigger bits info
      findTrigIndex(_trigEvtPS, moduleLabel, index);
      if ( flag.hasAnyProperty(prescalerFlag)){ 
	_trigEvtPS[index].label  = moduleLabel;
	_trigEvtPS[index].counts = _trigEvtPS[index].counts + 1;
      }
      
      
    }//end loop over the TriggerInfo Handles

    //now fill the correlation matrix
    for (size_t i=0; i<trigFlagAll_index.size(); ++i){
      for (size_t j=0; j<trigFlagAll_index.size(); ++j){
	_h2DTrigInfo[0]->Fill(trigFlagAll_index.at(i), trigFlagAll_index.at(j));
      }
    }
    
    if (trigFlagAll_index.size() == 1) _hTrigInfo[10]->Fill(trigFlagAll_index.at(0));

  }
  
  void   ReadTriggerInfo::findTrigIndex(std::vector<trigInfo_> Vec, std::string ModuleLabel, int &Index){
    //reset the index value
    Index = 0;
    for (size_t i=0; i<Vec.size(); ++i){
      if (Vec[i].label == ModuleLabel) { 
	Index = i;
	break;
      }else if (Vec[i].label != ""){
	Index = i+1;
      }
    }
  }

  void   ReadTriggerInfo::fillTrackTrigInfo(int TrkTrigIndex, const KalSeed*KSeed){
    int                nsh = (int)KSeed->hits().size();
    KalSegment const& fseg = KSeed->segments().front();

    double     ndof  = std::max(1.0,nsh - 5.0);
    double     p     = fseg.mom();
    double     chi2d = KSeed->chisquared()/ndof;
    double     pt    = p*std::cos(std::atan(fseg.helix().tanDip()));
    double     d0    = fseg.helix().d0();
    double     clE(-1.);
    if (KSeed->caloCluster()) clE = KSeed->caloCluster()->energyDep();

    _hTrkInfo[TrkTrigIndex][0]->Fill(p);
    _hTrkInfo[TrkTrigIndex][1]->Fill(pt);
    _hTrkInfo[TrkTrigIndex][2]->Fill(nsh);
    _hTrkInfo[TrkTrigIndex][3]->Fill(d0);
    _hTrkInfo[TrkTrigIndex][4]->Fill(chi2d);
    _hTrkInfo[TrkTrigIndex][5]->Fill(clE);
  }

  
  void   ReadTriggerInfo::fillHelixTrigInfo(int HelTrigIndex, const HelixSeed*HSeed){
    int        nsh       = (int)HSeed->hits().size();
    float      mm2MeV    = (3./10.)*_bz0;

    double     p         = HSeed->helix().momentum()*mm2MeV;
    double     chi2dZPhi = HSeed->helix().chi2dZPhi();
    double     chi2dXY   = HSeed->helix().chi2dXY();
    double     pt        = HSeed->helix().radius()*mm2MeV;
    double     d0        = HSeed->helix().rcent() - HSeed->helix().radius();
    double     clE(-1.);
    if (HSeed->caloCluster()) clE = HSeed->caloCluster()->energyDep();

    _hHelInfo[HelTrigIndex][0]->Fill(p);
    _hHelInfo[HelTrigIndex][1]->Fill(pt);
    _hHelInfo[HelTrigIndex][2]->Fill(nsh);
    _hHelInfo[HelTrigIndex][3]->Fill(d0);
    _hHelInfo[HelTrigIndex][4]->Fill(chi2dXY);
    _hHelInfo[HelTrigIndex][5]->Fill(chi2dZPhi);
    _hHelInfo[HelTrigIndex][6]->Fill(clE);
  }

  void   ReadTriggerInfo::fillCaloCalibTrigInfo(int ClCalibIndex, const CaloCluster*HCl){
    int        clsize    = HCl->size();
    double     energy    = HCl->energyDep();

    _hCaloCalibInfo[ClCalibIndex][0]->Fill(energy);
    _hCaloCalibInfo[ClCalibIndex][1]->Fill(clsize);
  }

   void   ReadTriggerInfo::fillCaloTrigSeedInfo(int Index, const CaloTrigSeed*HCl){
    _hCaloOnlyInfo[Index][0]->Fill(HCl->epeak());
    _hCaloOnlyInfo[Index][1]->Fill(HCl->ring1max());
    _hCaloOnlyInfo[Index][2]->Fill(HCl->ring1max2());
  }

  
  
}  

DEFINE_ART_MODULE(mu2e::ReadTriggerInfo);
