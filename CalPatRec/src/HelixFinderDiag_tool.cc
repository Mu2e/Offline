
#include "fhiclcpp/ParameterSet.h"

#include "Offline/CalPatRec/inc/HelixFinder_types.hh"

#include "Offline/Mu2eUtilities/inc/ModuleHistToolBase.hh"

#include "art/Framework/Principal/Event.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"

#include "TH1.h"
#include "TProfile.h"

namespace mu2e {

using namespace HelixFinderTypes;

class HelixFinderDiag : public mu2e::ModuleHistToolBase {

public:
  enum { kNEventHistsSets = 1, kNTimeClusterHistsSets = 2 };

  struct EventHists {
    TH1F* moduleTime;
    TH1F* nHelices;
    TH1F* nTimeClusters;
    TProfile* moduleTimeVSnTimeClusters;
  };

  struct TimeClusterHists {
    TH1F* nHelicesPerTC;
    TH1F* nComboHitsPerTC;
    TH1F* nStrawHitsPerTC;
    TH1F* timePerTC;
    TProfile* timePerTCVSnComboHitsPerTC;
  };

  struct Hists {
    EventHists* _eventHists[kNEventHistsSets];
    TimeClusterHists* _timeClusterHists[kNTimeClusterHistsSets];
  };

protected:
  Hists _hist;
  diagInfo* _data;

public:
  HelixFinderDiag(const fhicl::Table<mu2e::HelixFinderTypes::Config>& config);
  ~HelixFinderDiag();

private:
  int bookEventHistograms(EventHists* Hist, art::TFileDirectory* Dir);
  int bookTimeClusterHistograms(TimeClusterHists* Hist, art::TFileDirectory* Dir);

  int fillEventHistograms(EventHists* Hist, diagInfo* Data);
  int fillTimeClusterHistograms(TimeClusterHists* Hist, diagInfo* Data, int loopIndex);

  virtual int bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) override;
  virtual int fillHistograms(void* Data, int Mode = -1) override;
};

//-----------------------------------------------------------------------------
HelixFinderDiag::HelixFinderDiag(const fhicl::Table<mu2e::HelixFinderTypes::Config>& config) {}

//-----------------------------------------------------------------------------
HelixFinderDiag::~HelixFinderDiag() {}

//-----------------------------------------------------------------------------
int HelixFinderDiag::bookEventHistograms(EventHists* Hist, art::TFileDirectory* Dir) {

  Hist->moduleTime =
      Dir->make<TH1F>("moduleTime", "time (ms) per event spent at module", 10000, 0.0, 500.0);
  Hist->nHelices = Dir->make<TH1F>("nHelices", "number of helices found per event", 30, 0.0, 30.0);
  Hist->nTimeClusters =
      Dir->make<TH1F>("nTimeClusters", "number of time clusters per event", 100, 0.0, 100.0);
  Hist->moduleTimeVSnTimeClusters = Dir->make<TProfile>(
      "moduleTimeVSnTimeClusters", "moduleTimeVSnTimeClusters", 100, 0.0, 100.0, 0.0, 500.0, "i");

  return 0;
}

//-----------------------------------------------------------------------------
int HelixFinderDiag::bookTimeClusterHistograms(TimeClusterHists* Hist, art::TFileDirectory* Dir) {

  Hist->nHelicesPerTC =
      Dir->make<TH1F>("nHelicesPerTC", "number of helices found per TC", 30, 0, 30.0);
  Hist->nComboHitsPerTC =
      Dir->make<TH1F>("nComboHitsPerTC", "number of combo hits per TC", 200, 0, 200.0);
  Hist->nStrawHitsPerTC =
      Dir->make<TH1F>("nStrawHitsPerTC", "number of straw hits per TC", 300, 0, 300.0);
  Hist->timePerTC =
      Dir->make<TH1F>("timePerTC", "time (ms) searching for helix per TC", 50000, 0.0, 25.0);
  Hist->timePerTCVSnComboHitsPerTC = Dir->make<TProfile>(
      "timePerTCVSnComboHitsPerTC", "timePerTCVSnComboHitsPerTC", 200, 0.0, 200.0, 0.0, 25.0, "i");

  return 0;
}

//-----------------------------------------------------------------------------
int HelixFinderDiag::bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) {
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

  return 0;
}

//-----------------------------------------------------------------------------
int HelixFinderDiag::fillEventHistograms(EventHists* Hist, diagInfo* Data) {

  Hist->moduleTime->Fill(Data->moduleTime);
  Hist->nHelices->Fill(Data->nHelices);
  Hist->nTimeClusters->Fill(Data->nTimeClusters);
  Hist->moduleTimeVSnTimeClusters->Fill(Data->nTimeClusters, Data->moduleTime, 1);

  return 0;
}

//-----------------------------------------------------------------------------
int HelixFinderDiag::fillTimeClusterHistograms(TimeClusterHists* Hist, diagInfo* Data,
                                               int loopIndex) {

  // fill per tc info
  Hist->nHelicesPerTC->Fill(Data->timeClusterData.at(loopIndex).nHelices);
  Hist->nComboHitsPerTC->Fill(Data->timeClusterData.at(loopIndex).nComboHits);
  Hist->nStrawHitsPerTC->Fill(Data->timeClusterData.at(loopIndex).nStrawHits);
  Hist->timePerTC->Fill(Data->timeClusterData.at(loopIndex).time);
  Hist->timePerTCVSnComboHitsPerTC->Fill(Data->timeClusterData.at(loopIndex).nComboHits,
                                         Data->timeClusterData.at(loopIndex).time, 1);

  return 0;
}

//-----------------------------------------------------------------------------
// Mode is not used here
//-----------------------------------------------------------------------------
int HelixFinderDiag::fillHistograms(void* Data, int Mode) {

  _data = (diagInfo*)Data;

  //-----------------------------------------------------------------------------
  // fill event histograms
  //-----------------------------------------------------------------------------
  fillEventHistograms(_hist._eventHists[0], _data);

  //-----------------------------------------------------------------------------
  // fill time cluster histograms
  //-----------------------------------------------------------------------------
  for (int i = 0; i < (int)_data->timeClusterData.size(); i++) {
    // fill folder 0 always
    fillTimeClusterHistograms(_hist._timeClusterHists[0], _data, i);
    // fill folder 1 only if TC takes more than 3.0 ms to process
    if (_data->timeClusterData.at(i).time > 3.0) {
      fillTimeClusterHistograms(_hist._timeClusterHists[1], _data, i);
    }
  }

  return 0;
}

DEFINE_ART_CLASS_TOOL(HelixFinderDiag)

} // namespace mu2e
