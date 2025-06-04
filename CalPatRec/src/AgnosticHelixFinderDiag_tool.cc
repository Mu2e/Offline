#include "fhiclcpp/ParameterSet.h"

#include "Offline/CalPatRec/inc/AgnosticHelixFinder_types.hh"

#include "Offline/Mu2eUtilities/inc/ModuleHistToolBase.hh"

#include "art/Framework/Principal/Event.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"

#include "TH1.h"
#include "TProfile.h"

namespace mu2e {

  using namespace AgnosticHelixFinderTypes;

  class AgnosticHelixFinderDiag : public mu2e::ModuleHistToolBase {

  public:
    enum { kNEventHistsSets = 1,
           kNTimeClusterHistsSets = 1,
           kNHelixSeedHistsSets = 1,
           kNLineInfoHistsSets = 1,
           kNLineSegmentHistsSets = 1
    };

    struct EventHists {
      TH1F* nHelices;
      TH1F* nTimeClusters;
    };

    struct TimeClusterHists {
      TH1F* nHelicesPerTC;
      TH1F* nComboHitsPerTC;
      TH1F* nStrawHitsPerTC;
    };

    struct HelixSeedHists {
      TH1F* eDepAvgPerHS;
    };

    struct LineInfoHists {
      TH1F* nHitsRatioPerHS;
    };

    struct LineSegmentHists {
      TH1F* chi2dof;
      TH1F* maxHitGap;
    };

    struct Hists {
      EventHists* _eventHists[kNEventHistsSets];
      TimeClusterHists* _timeClusterHists[kNTimeClusterHistsSets];
      HelixSeedHists* _helixSeedHists[kNHelixSeedHistsSets];
      LineInfoHists* _lineInfoHists[kNLineInfoHistsSets];
      LineSegmentHists* _lineSegmentHists[kNLineSegmentHistsSets];
    };

  protected:
    Hists _hist;
    diagInfo* _data;

  public:
    AgnosticHelixFinderDiag(const fhicl::Table<mu2e::AgnosticHelixFinderTypes::Config>& config);
    ~AgnosticHelixFinderDiag();

  private:
    int bookEventHistograms(EventHists* Hist, art::TFileDirectory* Dir);
    int bookTimeClusterHistograms(TimeClusterHists* Hist, art::TFileDirectory* Dir);
    int bookHelixSeedHistograms(HelixSeedHists* Hist, art::TFileDirectory* Dir);
    int bookLineInfoHistograms(LineInfoHists* Hist, art::TFileDirectory* Dir);
    int bookLineSegmentHistograms(LineSegmentHists* Hist, art::TFileDirectory* Dir);

    int fillEventHistograms(EventHists* Hist, diagInfo* Data);
    int fillTimeClusterHistograms(TimeClusterHists* Hist, diagInfo* Data, int loopIndex);
    int fillHelixSeedHistograms(HelixSeedHists* Hist, diagInfo* Data, int loopIndex);
    int fillLineInfoHistograms(LineInfoHists* Hist, diagInfo* Data, int loopIndex);
    int fillLineSegmentHistograms(LineSegmentHists* Hist, diagInfo* Data, int loopIndex);

    virtual int bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) override;
    virtual int fillHistograms(void* Data, int Mode = -1) override;
  };

  //-----------------------------------------------------------------------------
  AgnosticHelixFinderDiag::AgnosticHelixFinderDiag(const fhicl::Table<mu2e::AgnosticHelixFinderTypes::Config>& config) {}

  //-----------------------------------------------------------------------------
  AgnosticHelixFinderDiag::~AgnosticHelixFinderDiag() {}

  //-----------------------------------------------------------------------------
  int AgnosticHelixFinderDiag::bookEventHistograms(EventHists* Hist, art::TFileDirectory* Dir) {

    Hist->nHelices = Dir->make<TH1F>("nHelices", "number of helices found per event", 30, 0.0, 30.0);
    Hist->nTimeClusters =
      Dir->make<TH1F>("nTimeClusters", "number of time clusters per event", 100, 0.0, 100.0);

    return 0;
  }

  //-----------------------------------------------------------------------------
  int AgnosticHelixFinderDiag::bookTimeClusterHistograms(TimeClusterHists* Hist, art::TFileDirectory* Dir) {

    Hist->nHelicesPerTC =
      Dir->make<TH1F>("nHelicesPerTC", "number of helices found per TC", 30, 0, 30.0);
    Hist->nComboHitsPerTC =
      Dir->make<TH1F>("nComboHitsPerTC", "number of combo hits per TC", 200, 0, 200.0);
    Hist->nStrawHitsPerTC =
      Dir->make<TH1F>("nStrawHitsPerTC", "number of straw hits per TC", 300, 0, 300.0);

    return 0;
  }

  //-----------------------------------------------------------------------------
  int AgnosticHelixFinderDiag::bookHelixSeedHistograms(HelixSeedHists* Hist, art::TFileDirectory* Dir) {

    Hist->eDepAvgPerHS =
      Dir->make<TH1F>("eDepAvgPerHS", "average combo hit eDep per helix seed", 1000, 0, 0.01);

    return 0;
  }

  //-----------------------------------------------------------------------------
  int AgnosticHelixFinderDiag::bookLineInfoHistograms(LineInfoHists* Hist, art::TFileDirectory* Dir) {

    Hist->nHitsRatioPerHS =
      Dir->make<TH1F>("nHitsRatioPerHS", "ratio of circle fit hits to phi fit hits", 1000, 0, 5.0);

    return 0;
  }

  //-----------------------------------------------------------------------------
  int AgnosticHelixFinderDiag::bookLineSegmentHistograms(LineSegmentHists* Hist, art::TFileDirectory* Dir) {

    Hist->chi2dof = Dir->make<TH1F>("chi2dof", "chi2dof of line segments", 300, 0, 30.0);
    Hist->maxHitGap = Dir->make<TH1F>("maxHitGap", "max dZ between adjacent points", 600, 0, 600.0);

    return 0;
  }

  //-----------------------------------------------------------------------------
  int AgnosticHelixFinderDiag::bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) {
    char folder_name[20];
    TH1::AddDirectory(0);

    //-----------------------------------------------------------------------------
    // book event histograms
    //-----------------------------------------------------------------------------
    for (int i = 0; i < kNEventHistsSets; i++) {
      sprintf(folder_name, "evt_%i", i);
      art::TFileDirectory tfdir = Tfs->mkdir(folder_name);
      _hist._eventHists[i] = new EventHists;
      bookEventHistograms(_hist._eventHists[i], &tfdir);
    }

    //-----------------------------------------------------------------------------
    // book time cluster histograms
    //-----------------------------------------------------------------------------
    for (int i = 0; i < kNTimeClusterHistsSets; i++) {
      sprintf(folder_name, "tcl_%i", i);
      art::TFileDirectory tfdir = Tfs->mkdir(folder_name);
      _hist._timeClusterHists[i] = new TimeClusterHists;
      bookTimeClusterHistograms(_hist._timeClusterHists[i], &tfdir);
    }

    //-----------------------------------------------------------------------------
    // book helix seed histograms
    //-----------------------------------------------------------------------------
    for (int i = 0; i < kNHelixSeedHistsSets; i++) {
      sprintf(folder_name, "hs_%i", i);
      art::TFileDirectory tfdir = Tfs->mkdir(folder_name);
      _hist._helixSeedHists[i] = new HelixSeedHists;
      bookHelixSeedHistograms(_hist._helixSeedHists[i], &tfdir);
    }

    //-----------------------------------------------------------------------------
    // book lineInfo histograms
    //-----------------------------------------------------------------------------
    for (int i = 0; i < kNLineInfoHistsSets; i++) {
      sprintf(folder_name, "li_%i", i);
      art::TFileDirectory tfdir = Tfs->mkdir(folder_name);
      _hist._lineInfoHists[i] = new LineInfoHists;
      bookLineInfoHistograms(_hist._lineInfoHists[i], &tfdir);
    }

    //-----------------------------------------------------------------------------
    // book line segment histograms
    //-----------------------------------------------------------------------------
    for (int i = 0; i < kNLineSegmentHistsSets; i++) {
      sprintf(folder_name, "ls_%i", i);
      art::TFileDirectory tfdir = Tfs->mkdir(folder_name);
      _hist._lineSegmentHists[i] = new LineSegmentHists;
      bookLineSegmentHistograms(_hist._lineSegmentHists[i], &tfdir);
    }

    return 0;
  }

  //-----------------------------------------------------------------------------
  int AgnosticHelixFinderDiag::fillEventHistograms(EventHists* Hist, diagInfo* Data) {

    Hist->nHelices->Fill(Data->nHelices);
    Hist->nTimeClusters->Fill(Data->nTimeClusters);

    return 0;
  }

  //-----------------------------------------------------------------------------
  int AgnosticHelixFinderDiag::fillTimeClusterHistograms(TimeClusterHists* Hist, diagInfo* Data,
                                                         int loopIndex) {

    // fill per tc info
    Hist->nHelicesPerTC->Fill(Data->timeClusterData.at(loopIndex).nHelices);
    Hist->nComboHitsPerTC->Fill(Data->timeClusterData.at(loopIndex).nComboHits);
    Hist->nStrawHitsPerTC->Fill(Data->timeClusterData.at(loopIndex).nStrawHits);

    return 0;
  }

  //-----------------------------------------------------------------------------
  int AgnosticHelixFinderDiag::fillHelixSeedHistograms(HelixSeedHists* Hist, diagInfo* Data,
                                                         int loopIndex) {

    // fill per hs info
    Hist->eDepAvgPerHS->Fill(Data->helixSeedData.at(loopIndex).eDepAvg);

    return 0;
  }

  //-----------------------------------------------------------------------------
  int AgnosticHelixFinderDiag::fillLineInfoHistograms(LineInfoHists* Hist, diagInfo* Data,
                                                      int loopIndex) {

    // fill per line info
    Hist->nHitsRatioPerHS->Fill(Data->lineInfoData.at(loopIndex).nHitsRatio);

    return 0;
  }

  //-----------------------------------------------------------------------------
  int AgnosticHelixFinderDiag::fillLineSegmentHistograms(LineSegmentHists* Hist, diagInfo* Data,
                                                         int loopIndex) {

    // fill per line segment info
    Hist->chi2dof->Fill(Data->lineSegmentData.at(loopIndex).chi2dof);
    Hist->maxHitGap->Fill(Data->lineSegmentData.at(loopIndex).maxHitGap);

    return 0;
  }

  //-----------------------------------------------------------------------------
  // Mode is not used here
  //-----------------------------------------------------------------------------
  int AgnosticHelixFinderDiag::fillHistograms(void* Data, int Mode) {

    _data = (diagInfo*)Data;

    //-----------------------------------------------------------------------------
    // fill event histograms
    //-----------------------------------------------------------------------------
    fillEventHistograms(_hist._eventHists[0], _data);

    //-----------------------------------------------------------------------------
    // fill time cluster histograms
    //-----------------------------------------------------------------------------
    for (int i = 0; i < (int)_data->timeClusterData.size(); i++) {
      fillTimeClusterHistograms(_hist._timeClusterHists[0], _data, i);
    }

    //-----------------------------------------------------------------------------
    // fill helix seed histograms
    //-----------------------------------------------------------------------------
    for (int i = 0; i < (int)_data->helixSeedData.size(); i++) {
      if (_data->helixSeedData.at(i).eDepAvg > 1e-6) {
        fillHelixSeedHistograms(_hist._helixSeedHists[0], _data, i);
      }
    }

    //-----------------------------------------------------------------------------
    // fill lineInfo seed histograms
    //-----------------------------------------------------------------------------
    for (int i = 0; i < (int)_data->lineInfoData.size(); i++) {
      fillLineInfoHistograms(_hist._lineInfoHists[0], _data, i);
    }


    //-----------------------------------------------------------------------------
    // fill line segment histograms
    //-----------------------------------------------------------------------------
    for (int i = 0; i < (int)_data->lineSegmentData.size(); i++) {
      fillLineSegmentHistograms(_hist._lineSegmentHists[0], _data, i);
    }

    return 0;
  }

  DEFINE_ART_CLASS_TOOL(AgnosticHelixFinderDiag)

} // namespace mu2e
