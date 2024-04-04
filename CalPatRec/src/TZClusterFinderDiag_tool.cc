
#include "fhiclcpp/ParameterSet.h"

#include "Offline/CalPatRec/inc/TZClusterFinder_types.hh"

#include "Offline/Mu2eUtilities/inc/McUtilsToolBase.hh"
#include "Offline/Mu2eUtilities/inc/ModuleHistToolBase.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/IntensityInfoTimeCluster.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/ProtonBunchIntensity.hh"

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art/Framework/Principal/Event.h"

#include "TH1.h"
#include "TProfile.h"


namespace mu2e {

  using namespace TZClusterFinderTypes;

  class TZClusterFinderDiag : public mu2e::ModuleHistToolBase {

  public:

    enum {
      kNEventHistsSets = 1,
      kNTimeClusterHistsSets = 1
    };

    struct EventHists {
      TH1F*     nClusters;
      TH1F*     nProtonsPred;
      TH1F*     nProtonsMC;
      TH1F*     nProtPredMinusMC;
      TProfile* nProtPredvsMC;
      TProfile* nProtMCvsPOT;
      TProfile* nProtPredvsPOT;
    };

    struct TimeClusterHists {
      TH1F* clusterSlopes;
      TH1F* cHitsInClusters;
      TH1F* Chi2DOFHist;
    };

    struct Hists {
      EventHists*       _eventHists[kNEventHistsSets];
      TimeClusterHists* _timeClusterHists[kNTimeClusterHistsSets];
    };

  protected:
    Hists                            _hist;
    Data_t*                          _data;
    std::vector<mcSimIDs>            _mcHits;
    int                              _mcTruth;
    int                              _simIDThresh;
    std::unique_ptr<McUtilsToolBase> _mcUtils;
    const mu2e::SimParticle*         _simParticle;
    float                            _nPOT;

  public:

    TZClusterFinderDiag(const fhicl::Table<mu2e::TZClusterFinderTypes::Config>& config);
    ~TZClusterFinderDiag();

  private:

    int  bookEventHistograms      (EventHists*        Hist, art::TFileDirectory* Dir);
    int  bookTimeClusterHistograms(TimeClusterHists*  Hist, art::TFileDirectory* Dir);

    int  fillEventHistograms      (EventHists*       Hist, Data_t* Data);
    int  fillTimeClusterHistograms(TimeClusterHists* Hist, Data_t* Data, int loopIndex);


    virtual int bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) override ;
    virtual int fillHistograms(void* Data, int Mode = -1) override ;
  };

//-----------------------------------------------------------------------------
  TZClusterFinderDiag::TZClusterFinderDiag(const fhicl::Table<mu2e::TZClusterFinderTypes::Config>& config):
    _mcTruth               (config().mcTruth()               ),
    _simIDThresh           (config().simIDThresh()           )
  {
    if (_mcTruth != 0) _mcUtils = art::make_tool<McUtilsToolBase>(config().mcUtils,"mcUtils");
    else               _mcUtils = std::make_unique<McUtilsToolBase>();
  }

//-----------------------------------------------------------------------------
  TZClusterFinderDiag::~TZClusterFinderDiag() {}

//-----------------------------------------------------------------------------
  int TZClusterFinderDiag::bookEventHistograms(EventHists* Hist, art::TFileDirectory* Dir) {
    Hist->nClusters = Dir->make<TH1F>("nClusters" , "number of clusters", 50, 0.0, 50.0);
    Hist->nProtonsPred = Dir->make<TH1F>("nProtonsPred" , "number of protons clusters counted (>= 15 straw hits)", 30, 0.0, 30.0 );
    Hist->nProtonsMC = Dir->make<TH1F>("nProtonsMC" , "number of MC truth protons (>= 15 straw hits)", 30, 0.0, 30.0 );
    Hist->nProtPredMinusMC = Dir->make<TH1F>("nProtPredMinusMC" , "number of protons predicted minus truth (>= 15 straw hits)", 20, -10.0, 10.0 );
    Hist->nProtPredvsMC = Dir->make<TProfile>("nProtPredvsMC" , "profile number of protons predicted vs truth (>= 15 straw hits)", 30, 0, 30, 0, 30, "i" );
    Hist->nProtMCvsPOT = Dir->make<TProfile>("nProtMCvsPOT" , "profile MC number of protons vs nPOT", 15, 1e6, 120e6, 0, 30, "i" );
    Hist->nProtPredvsPOT = Dir->make<TProfile>("nProtPredvsPOT" , "profile number of protons predicted vs nPOT", 15, 1e6, 120e6, 0, 30, "i" );

    return 0;
  }

//-----------------------------------------------------------------------------
  int TZClusterFinderDiag::bookTimeClusterHistograms(TimeClusterHists* Hist, art::TFileDirectory* Dir) {
    Hist->clusterSlopes = Dir->make<TH1F>("clusterSlopes" , "cluster slopes", 100, -0.06, 0.06 );
    Hist->cHitsInClusters = Dir->make<TH1F>("cHitsInClusters" , "number of combo hits in cluster", 100, 0., 100. );
    Hist->Chi2DOFHist = Dir->make<TH1F>("Chi2DOFHist" , "line fit reduced chi-squared", 1000, -2.0, 50. );
    return 0;
  }

//-----------------------------------------------------------------------------
  int TZClusterFinderDiag::bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) {
    char folder_name[20];
    TH1::AddDirectory(0);

//-----------------------------------------------------------------------------
// book event histograms
//-----------------------------------------------------------------------------
    for (int i=0; i<kNEventHistsSets; i++) {
      sprintf(folder_name,"evt_%i",i);
      art::TFileDirectory tfdir = Tfs->mkdir(folder_name);
      _hist._eventHists[i] = new EventHists;
      bookEventHistograms(_hist._eventHists[i],&tfdir);
    }

//-----------------------------------------------------------------------------
// book time cluster histograms
//-----------------------------------------------------------------------------
    for (int i=0; i<kNTimeClusterHistsSets; i++) {
      sprintf(folder_name,"tcl_%i",i);
      art::TFileDirectory tfdir = Tfs->mkdir(folder_name);
      _hist._timeClusterHists[i] = new TimeClusterHists;
      bookTimeClusterHistograms(_hist._timeClusterHists[i],&tfdir);
    }

    return 0;

  }

//-----------------------------------------------------------------------------
  int TZClusterFinderDiag::fillEventHistograms(EventHists* Hist, Data_t* Data) {

    Hist->nClusters->Fill(Data->_nTZClusters);

    // ---------------------
    // make proton predicted plots
    // ---------------------
    int nProtonsCounted = (int)Data->_iiTC->nProtonTCs();
    Hist->nProtonsPred->Fill(nProtonsCounted);

    // ---------------------
    // make proton truth plots
    // ---------------------
    if ( _mcTruth != 0 ) {
      int nProtonsTruth = 0;
      for (size_t i=0; i<_mcHits.size(); i++) {
        if (_mcHits[i].pdgID != 2212) {continue;}
        if (_mcHits[i].nHits >= _simIDThresh) { nProtonsTruth++; }
      }
      int diff = nProtonsCounted - nProtonsTruth;
      Hist->nProtonsMC->Fill(nProtonsTruth);
      Hist->nProtPredMinusMC->Fill(diff);
      Hist->nProtPredvsMC->Fill(nProtonsTruth, nProtonsCounted, 1);
      Hist->nProtMCvsPOT->Fill(_nPOT, nProtonsTruth);
      Hist->nProtPredvsPOT->Fill(_nPOT, nProtonsCounted);
    }

    return 0;
  }

//-----------------------------------------------------------------------------
  int TZClusterFinderDiag::fillTimeClusterHistograms(TimeClusterHists* Hist, Data_t* Data, int loopIndex) {

    // fill slopes of clusters
    float slope = Data->lineSlope.at(loopIndex);
    Hist->clusterSlopes->Fill(slope);

    // fill number of combo hits in clusters
    int numHits = (int)Data->_tcColl->at(loopIndex)._strawHitIdxs.size();
    Hist->cHitsInClusters->Fill(numHits);

    // fill reduced chi-squared of lines fit to clusters
    float reducedChiSquare = Data->chi2DOF.at(loopIndex);
    Hist->Chi2DOFHist->Fill(reducedChiSquare);

    return 0;
  }

//-----------------------------------------------------------------------------
// Mode is not used here
//-----------------------------------------------------------------------------
  int TZClusterFinderDiag::fillHistograms(void* Data, int Mode) {
    _data = (Data_t*) Data;

    // get PBI info
    _nPOT  = -1.;
    art::Handle<ProtonBunchIntensity> evtWeightH;
    _data->_event->getByLabel("PBISim", evtWeightH);
    if (evtWeightH.isValid()){
      _nPOT  = (float)evtWeightH->intensity();
    }

    // fill mcSimIDs data members simID, nHits, and pdgID
    if ( _mcTruth != 0 ) {
      for (size_t i=0; i<_data->_chColl2->size(); ++i){
        std::vector<StrawDigiIndex> shids;
        _data->_chColl2->fillStrawDigiIndices(i,shids);
        int SimID = _mcUtils->strawHitSimId(_data->_event,shids[0]);
        _simParticle = _mcUtils->getSimParticle(_data->_event,shids[0]);
        int _pdgID = _simParticle->pdgId();
        if (i==0) {
          mcSimIDs sim_info;
          sim_info.simID = SimID;
          sim_info.pdgID = _pdgID;
          sim_info.nHits = 1;
          _mcHits.push_back(sim_info);
          continue;
        }
        int alreadyStored = 0;
        for (size_t j=0; j<_mcHits.size(); j++) {
          if (_mcHits[j].simID == SimID) {
            alreadyStored = 1;
            _mcHits[j].nHits++;
            break;
          }
        }
        if (alreadyStored == 1) {continue;}
        else {
          mcSimIDs sim_info;
          sim_info.simID = SimID;
          sim_info.pdgID = _pdgID;
          sim_info.nHits = 1;
          _mcHits.push_back(sim_info);
        }
      }
    }

//-----------------------------------------------------------------------------
// fill event histograms
//-----------------------------------------------------------------------------
    fillEventHistograms(_hist._eventHists[0],_data);

//-----------------------------------------------------------------------------
// fill time cluster histograms
//-----------------------------------------------------------------------------
    for (int i=0; i<(int)_data->_tcColl->size(); i++) {
      fillTimeClusterHistograms(_hist._timeClusterHists[0],_data, i);
    }

    // clear _mcHits before next event
    _mcHits.clear();

    return 0;
  }

  DEFINE_ART_CLASS_TOOL(TZClusterFinderDiag)

}
