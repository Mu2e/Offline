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
      kNTrigInfo     = 10,
      kNTrackTrig    = 10,
      kNTrackTrigVar = 6,
      kNHelixTrig    = 10,
      kNHelixTrigVar = 7,
      kNCaloCalib    = 5,
      kNCaloCalibVar = 5,
      kNCaloOnly     = 5,
      kNCaloOnlyVar  = 5      
    };
    
    explicit ReadTriggerInfo(fhicl::ParameterSet const& pset);
    virtual ~ReadTriggerInfo() { }

    virtual void beginJob();
    virtual void endJob();
    virtual void endSubRun(const art::SubRun& sr);

    // This is called for each event.
    virtual void analyze(const art::Event& e);
    virtual void beginRun(const art::Run & run);

    void   findTrigIndex        (std::vector<std::string> Vec, std::string ModuleLabel, int &Index);
    void   fillTrackTrigInfo    (int TrkTrigIndex, const KalSeed*  KSeed);
    void   fillHelixTrigInfo    (int HelTrigIndex, const HelixSeed*HSeed);
    void   fillCaloTrigSeedInfo (int CTrigSeedIndex, const CaloTrigSeed*HCl);
    void   fillCaloCalibTrigInfo(int ClCalibIndex, const CaloCluster*  HCl);

  private:

    int      _diagLevel;
    size_t   _nMaxTrig;
    
    int      _nTrackTrig;
    int      _nCaloTrig;
    int      _nCaloCalibTrig;

    int      _nProcess;
    int      _numEvents;
    
    double   _bz0;
    
    std::vector<std::string>  _trigLabels;
    std::vector<int>          _trigCounts;
    
    std::vector<std::string>  _trigLabelsCaloOnly;
    std::vector<int>          _trigCountsCaloOnly;
    
    std::vector<std::string>  _trigLabelsCaloCalib;
    std::vector<int>          _trigCountsCaloCalib;
    
    std::vector<std::string>  _trigLabelsTrack;
    std::vector<int>          _trigCountsTrack;
    
    std::vector<std::string>  _trigLabelsHelix;
    std::vector<int>          _trigCountsHelix;
    
    std::vector<std::string>  _trigLabelsEvtPS;
    std::vector<int>          _trigCountsEvtPS;
    
    TH1F *_hTrigInfo  [kNTrigInfo];
    TH2F *_h2DTrigInfo[kNTrigInfo];

    TH1F *_hTrkInfo [kNTrackTrig][kNTrackTrigVar];
    TH1F *_hHelInfo [kNHelixTrig][kNHelixTrigVar];

    TH1F *_hCaloCalibInfo[kNCaloCalib][kNCaloCalibVar];
    TH1F *_hCaloOnlyInfo [kNCaloOnly][kNCaloOnlyVar];
  };


  ReadTriggerInfo::ReadTriggerInfo(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset), 
    _diagLevel     (pset.get<int>("diagLevel", 0)),
    _nMaxTrig      (pset.get<size_t>("nFilters", 50)),
    _nTrackTrig    (pset.get<size_t>("nTrackTriggers", 4)),
    _nCaloTrig     (pset.get<size_t>("nCaloTriggers", 4)),
    _nCaloCalibTrig(pset.get<size_t>("nCaloCalibTriggers", 4)),
    _nProcess(0), 
    _numEvents(0)
  {
    _trigLabels.resize(_nMaxTrig);
    _trigCounts.resize(_nMaxTrig);

    _trigLabelsCaloOnly.resize(_nMaxTrig);
    _trigCountsCaloOnly.resize(_nMaxTrig);

    _trigLabelsCaloCalib.resize(_nMaxTrig);
    _trigCountsCaloCalib.resize(_nMaxTrig);

    _trigLabelsTrack.resize(_nMaxTrig);
    _trigCountsTrack.resize(_nMaxTrig);

    _trigLabelsHelix.resize(_nMaxTrig);
    _trigCountsHelix.resize(_nMaxTrig);

    _trigLabelsEvtPS.resize(_nMaxTrig);
    _trigCountsEvtPS.resize(_nMaxTrig);
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

    _h2DTrigInfo[0] = trigInfoDir.make<TH2F>("h2DTrigInfo_map"     , "Trigger correlation map"                    , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5), (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));       

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
    

  }

  void ReadTriggerInfo::endJob(){
    //set the labels of the histrograms
    for (size_t i=0; i<_trigCounts.size(); ++i ){
      _hTrigInfo  [0]->GetXaxis()->SetBinLabel(i+1, _trigLabels[i].c_str());
      _h2DTrigInfo[0]->GetXaxis()->SetBinLabel(i+1, _trigLabels[i].c_str());
      if (_trigCounts[i] > 0) {
	_hTrigInfo[0]->SetBinContent(i+1, _nProcess/_trigCounts[i]);
	for (size_t j=0; j<_trigCounts.size(); ++j ){
	  _h2DTrigInfo[0]->GetYaxis()->SetBinLabel(j+1, _trigLabels[j].c_str());
	  if (_trigCounts[j] > 0) _h2DTrigInfo[0]->SetBinContent(i+1, j+1, _trigCounts[j]);
	}
      }

      _hTrigInfo[1]->GetXaxis()->SetBinLabel(i+1, _trigLabelsTrack[i].c_str());
      if (_trigCountsTrack[i] > 0) _hTrigInfo[1]->SetBinContent(i+1, _nProcess/_trigCountsTrack[i]);

      _hTrigInfo[2]->GetXaxis()->SetBinLabel(i+1, _trigLabelsCaloOnly[i].c_str());
      if (_trigCountsCaloOnly[i] > 0) _hTrigInfo[2]->SetBinContent(i+1, _nProcess/_trigCountsCaloOnly[i]);

      _hTrigInfo[3]->GetXaxis()->SetBinLabel(i+1, _trigLabelsEvtPS[i].c_str());
      if (_trigCountsEvtPS[i] > 0) _hTrigInfo[3]->SetBinContent(i+1, _trigCountsEvtPS[i]);

      _hTrigInfo[4]->GetXaxis()->SetBinLabel(i+1, _trigLabelsHelix[i].c_str());
      if (_trigCountsHelix[i] > 0) _hTrigInfo[4]->SetBinContent(i+1, _trigCountsHelix[i]);

      _hTrigInfo[5]->GetXaxis()->SetBinLabel(i+1, _trigLabelsCaloCalib[i].c_str());
      if (_trigCountsCaloCalib[i] > 0) _hTrigInfo[5]->SetBinContent(i+1, _trigCountsCaloCalib[i]);
    }
  }

  //================================================================
  void ReadTriggerInfo::beginRun(const art::Run & run){
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
    TriggerFlag                    caloOrTrackFlag     = trackFlag; caloOrTrackFlag.merge(caloFlag); caloOrTrackFlag.merge(caloCalibFlag); caloOrTrackFlag.merge(caloTrigSeedFlag); caloOrTrackFlag.merge(helixFlag);
    
    for (size_t i=0; i<hTrigInfoVec.size(); ++i){
      hTrigInfo = hTrigInfoVec.at(i);
      if (!hTrigInfo.isValid())         continue;
      const TriggerInfo* trigInfo  = hTrigInfo.product();
      const TriggerFlag  flag      = trigInfo->triggerBits();

      std::string    moduleLabel   = hTrigInfo.provenance()->moduleLabel();
      int            index(0);         
      
      //fill the Global Trigger bits info
      findTrigIndex(_trigLabels, moduleLabel, index);
      _trigLabels[index] = moduleLabel;
      if ( flag.hasAnyProperty(caloOrTrackFlag)){ 
	_trigCounts[index] = _trigCounts[index] + 1;
      }
      
      //fill the Calo-Only Trigger bits info
      findTrigIndex(_trigLabelsCaloOnly, moduleLabel, index);
      if ( flag.hasAnyProperty(caloFlag) || flag.hasAnyProperty(caloTrigSeedFlag)){ 
	_trigLabelsCaloOnly[index] = moduleLabel;
	_trigCountsCaloOnly[index] = _trigCountsCaloOnly[index] + 1;
	const CaloTrigSeed*clseed = trigInfo->caloTrigSeed().get();
	if(clseed) fillCaloTrigSeedInfo(index, clseed);
      }

      //fill the CaloCalib Trigger bits info
      findTrigIndex(_trigLabelsCaloCalib, moduleLabel, index);
      if ( flag.hasAnyProperty(caloCalibFlag)){ 
	_trigLabelsCaloCalib[index] = moduleLabel;
	_trigCountsCaloCalib[index] = _trigCountsCaloCalib[index] + 1;
	const CaloCluster*cluster = trigInfo->caloCluster().get();
	if(cluster) fillCaloCalibTrigInfo(index, cluster);
      }
      
      //fill the Track Trigger bits info
      findTrigIndex(_trigLabelsTrack, moduleLabel, index);
      if ( flag.hasAnyProperty(trackFlag)){ 
	_trigLabelsTrack[index] = moduleLabel;
	_trigCountsTrack[index] = _trigCountsTrack[index] + 1;
	const KalSeed*kseed = trigInfo->track().get();
	if(kseed) fillTrackTrigInfo(index, kseed);
      }
       //fill the Helix Trigger bits info
      findTrigIndex(_trigLabelsHelix, moduleLabel, index);
      if ( flag.hasAnyProperty(helixFlag)){ 
	_trigLabelsHelix[index] = moduleLabel;
	_trigCountsHelix[index] = _trigCountsHelix[index] + 1;
	const HelixSeed*hseed = trigInfo->helix().get();
	if(hseed) fillHelixTrigInfo(index, hseed);
      }
      
      //fill the Event-Prescaler Trigger bits info
      findTrigIndex(_trigLabelsEvtPS, moduleLabel, index);
      if ( flag.hasAnyProperty(prescalerFlag)){ 
	_trigLabelsEvtPS[index] = moduleLabel;
	_trigCountsEvtPS[index] = _trigCountsEvtPS[index] + 1;
      }
      
      
    }//end loop over the TriggerInfo Handles
  }
  
  void   ReadTriggerInfo::findTrigIndex(std::vector<std::string> Vec, std::string ModuleLabel, int &Index){
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
