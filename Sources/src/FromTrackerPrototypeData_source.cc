#include <iostream>
#include <fstream>
#include <boost/utility.hpp>
#include <cassert>
#include <set>
#include <string>

#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "fhiclcpp/types/Name.h"
#include "art/Framework/IO/Sources/Source.h"
#include "art/Framework/Core/InputSourceMacros.h"
#include "art/Framework/IO/Sources/SourceHelper.h"
#include "art/Framework/Principal/RunPrincipal.h"
#include "art/Framework/Principal/SubRunPrincipal.h"
#include "art/Framework/Principal/EventPrincipal.h"
#include "art/Framework/IO/Sources/put_product_in_principal.h"
#include "canvas/Persistency/Provenance/Timestamp.h"
#include "canvas/Persistency/Provenance/RunID.h"
#include "canvas/Persistency/Provenance/SubRunID.h"
#include "canvas/Persistency/Provenance/EventID.h"
#include "canvas/Persistency/Provenance/BranchType.h"
#include "canvas/Persistency/Provenance/ProductID.h"
#include "canvas/Persistency/Provenance/canonicalProductName.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "Offline/DataProducts/inc/TrkTypes.hh"
#include "Offline/DataProducts/inc/StrawId.hh"

#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"

#include "Offline/SeedService/inc/SeedService.hh"

using CLHEP::Hep3Vector;
using CLHEP::HepLorentzVector;

using namespace std;

namespace mu2e {
  struct Config
  {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    fhicl::Sequence<std::string> inputFiles{Name("fileNames"),Comment("Input flat TTree Root file")};
    fhicl::Atom<unsigned int> runNumber{Name("runNumber"), Comment("Run number")};
    fhicl::Atom<unsigned int> maxEvents{Name("maxEvents"), Comment("Max number of events"), 0};

    // These are used by art and are required.
    fhicl::Atom<std::string> module_label{Name("module_label"), Comment("Art module label"), ""};
    fhicl::Atom<std::string> module_type{Name("module_type"), Comment("Art module type"), ""};
  };
  typedef fhicl::WrappedTable<Config> Parameters;


  //================================================================
  class TrackerPrototypeDataDetail : private boost::noncopyable {
    std::string myModuleLabel_;
    art::SourceHelper const& pm_;
    unsigned runNumber_; // from ParSet
    art::SubRunID lastSubRunID_;
    std::set<art::SubRunID> seenSRIDs_;

    std::string currentFileName_;

    TFile *currentFile_ = nullptr;
    TTree *tree_ = nullptr;

    TBranch *bsamples;
    int _panel, _channel;
    int _run;
    uint32_t _ewm;
    uint32_t _tdcCal, _tdcHV;
    uint16_t _totCal, _totHV;
    uint16_t _pmp;
    std::vector<uint16_t> *_samples = 0;
    int entryIndex_;

    mu2e::StrawDigiCollection digis_;
    mu2e::StrawDigiADCWaveformCollection digiadcs_;

    unsigned currentSubRunNumber_; // from file
    unsigned currentEventNumber_;
    unsigned maxEvents_;

    bool first_call;
    int printAtEvent;
    int currentEvent;

    // A helper function used to manage the principals.
    // This is boilerplate that does not change if you change the data products.
    void managePrincipals ( int runNumber,
        int subRunNumber,
        int eventNumber,
        art::RunPrincipal*&    outR,
        art::SubRunPrincipal*& outSR,
        art::EventPrincipal*&  outE);

    public:
    TrackerPrototypeDataDetail(const Parameters &conf,
        art::ProductRegistryHelper &,
        const art::SourceHelper &);

    void readFile(std::string const& filename, art::FileBlock*& fb);

    bool readNext(art::RunPrincipal* const& inR,
        art::SubRunPrincipal* const& inSR,
        art::RunPrincipal*& outR,
        art::SubRunPrincipal*& outSR,
        art::EventPrincipal*& outE);

    void closeCurrentFile();
  };

  //----------------------------------------------------------------
  TrackerPrototypeDataDetail::TrackerPrototypeDataDetail(const Parameters& conf,
      art::ProductRegistryHelper& rh,
      const art::SourceHelper& pm)
    : myModuleLabel_("FromTrackerPrototypeData")
      , pm_(pm)
      , runNumber_(conf().runNumber())
      , currentSubRunNumber_(-1U)
      , currentEventNumber_(0)
      , maxEvents_(conf().maxEvents())
  {
    if(!art::RunID(runNumber_).isValid()) {
      throw cet::exception("BADCONFIG", " FromTrackerPrototypeData: ")
        << " fhicl::ParameterSet specifies an invalid runNumber = "<<runNumber_<<"\n";
    }

    rh.reconstitutes<mu2e::StrawDigiCollection,art::InEvent>(myModuleLabel_);
    rh.reconstitutes<mu2e::StrawDigiADCWaveformCollection,art::InEvent>(myModuleLabel_);

    currentFileName_ = conf().inputFiles()[0];

    currentFile_ = new TFile(currentFileName_.c_str());
    tree_ = dynamic_cast<TTree*>(currentFile_->Get("T"));

    tree_->SetBranchAddress("run",&_run);
    tree_->SetBranchAddress("panel",&_panel);
    tree_->SetBranchAddress("channel",&_channel);
    tree_->SetBranchAddress("ewm",&_ewm);
    tree_->SetBranchAddress("tdcCal",&_tdcCal);
    tree_->SetBranchAddress("tdcHV",&_tdcHV);
    tree_->SetBranchAddress("totCal",&_totCal);
    tree_->SetBranchAddress("totHV",&_totHV);
    tree_->SetBranchAddress("pmp",&_pmp);

    bsamples = 0;
    tree_->SetBranchAddress("samples",&_samples,&bsamples);

    entryIndex_ = 0;

    tree_->GetEntry(0);
    currentSubRunNumber_ = _run;
    currentEventNumber_ = _ewm;

    first_call = true;
    printAtEvent = 0;
    currentEvent = 0;

  }


  //----------------------------------------------------------------
  void TrackerPrototypeDataDetail::readFile(const std::string& filename, art::FileBlock*& fb) {

    currentFileName_ = filename;

    currentFile_ = new TFile(currentFileName_.c_str());
    tree_ = dynamic_cast<TTree*>(currentFile_->Get("T"));

    tree_->SetBranchAddress("run",&_run);
    tree_->SetBranchAddress("panel",&_panel);
    tree_->SetBranchAddress("channel",&_channel);
    tree_->SetBranchAddress("ewm",&_ewm);
    tree_->SetBranchAddress("tdcCal",&_tdcCal);
    tree_->SetBranchAddress("tdcHV",&_tdcHV);
    tree_->SetBranchAddress("totCal",&_totCal);
    tree_->SetBranchAddress("totHV",&_totHV);
    tree_->SetBranchAddress("pmp",&_pmp);

    bsamples = 0;
    tree_->SetBranchAddress("samples",&_samples,&bsamples);

    entryIndex_ = 0;

    tree_->GetEntry(0);
    currentSubRunNumber_ = _run;
    currentEventNumber_ = _ewm;

    fb = new art::FileBlock(art::FileFormatVersion(1, "TrackerPrototypeDataInput"), currentFileName_);
  }

  //----------------------------------------------------------------
  void TrackerPrototypeDataDetail::closeCurrentFile() {
    currentFileName_ = "";
    currentFile_->Close();
  }

  //----------------------------------------------------------------
  bool TrackerPrototypeDataDetail::readNext(art::RunPrincipal* const& inR,
      art::SubRunPrincipal* const& inSR,
      art::RunPrincipal*& outR,
      art::SubRunPrincipal*& outSR,
      art::EventPrincipal*& outE)
  {
    if (maxEvents_ > 0 && static_cast<unsigned int>(currentEvent) >= maxEvents_)
      return false;
    while (true) {
      if (entryIndex_ >= tree_->GetEntries()){
        return false;
      }

      tree_->GetEntry(entryIndex_);
      if (static_cast<unsigned>(_run) != currentSubRunNumber_ || _ewm != currentEventNumber_){
        managePrincipals(runNumber_, currentSubRunNumber_, currentEventNumber_, outR, outSR, outE);
        std::unique_ptr<mu2e::StrawDigiCollection> updigis(new mu2e::StrawDigiCollection);
        for (size_t i=0;i<digis_.size();i++){
          updigis->push_back(digis_[i]);
        }
        std::unique_ptr<mu2e::StrawDigiADCWaveformCollection> updigiadcs(new mu2e::StrawDigiADCWaveformCollection);
        for (size_t i=0;i<digiadcs_.size();i++){
          updigiadcs->push_back(digiadcs_[i]);
        }
        art::put_product_in_principal(std::move(updigis), *outE, myModuleLabel_);
        art::put_product_in_principal(std::move(updigiadcs), *outE, myModuleLabel_);
        digis_.clear();
        digiadcs_.clear();
        currentSubRunNumber_ = _run;
        currentEventNumber_ = _ewm;
        return true;
      }

      mu2e::TrkTypes::TDCValues tdcs = {_tdcCal, _tdcHV};
      mu2e::TrkTypes::TOTValues tots = {_totCal, _totHV};
      mu2e::TrkTypes::ADCValue firmwarepmp = _pmp;
      mu2e::TrkTypes::ADCWaveform adc;

      Long64_t tentry = tree_->LoadTree(entryIndex_);
      bsamples->GetEntry(tentry);
      for (size_t j=0;j<_samples->size();j++){
        adc.push_back(_samples->at(j));
      }
      mu2e::StrawId sid(0,_panel,_channel);
      auto sd = mu2e::StrawDigi(sid,tdcs,tots,firmwarepmp);
      auto sda = mu2e::StrawDigiADCWaveform(adc);
      digis_.push_back(sd);
      digiadcs_.push_back(sda);
      entryIndex_++;

      if (entryIndex_ == tree_->GetEntries()){
        managePrincipals(runNumber_, currentSubRunNumber_, currentEventNumber_, outR, outSR, outE);
        std::unique_ptr<mu2e::StrawDigiCollection> updigis(new mu2e::StrawDigiCollection);
        for (size_t i=0;i<digis_.size();i++){
          updigis->push_back(digis_[i]);
        }
        std::unique_ptr<mu2e::StrawDigiADCWaveformCollection> updigiadcs(new mu2e::StrawDigiADCWaveformCollection);
        for (size_t i=0;i<digiadcs_.size();i++){
          updigiadcs->push_back(digiadcs_[i]);
        }
        art::put_product_in_principal(std::move(updigis), *outE, myModuleLabel_);
        art::put_product_in_principal(std::move(updigiadcs), *outE, myModuleLabel_);
        digis_.clear();
        digiadcs_.clear();
        std::cout << "DONE!" << std::endl;
        return true;
      }
    }
    return true;
  } // readNext()


  // Each time that we encounter a new run, a new subRun or a new event, we need to make a new principal
  // of the appropriate type.  This code does not need to change as the number and type of data products changes.
  void TrackerPrototypeDataDetail::managePrincipals ( int runNumber,
      int subRunNumber,
      int eventNumber,
      art::RunPrincipal*&    outR,
      art::SubRunPrincipal*& outSR,
      art::EventPrincipal*&  outE){

    art::Timestamp ts;

    if (currentEvent == printAtEvent){
      std::cout << "Event " << currentEvent << std::endl;
      printAtEvent = (printAtEvent+1)*2-1;
    }
    currentEvent++;


    if (first_call){
      outR = pm_.makeRunPrincipal(runNumber, ts);
      first_call = false;
    }
    art::SubRunID newID(runNumber, subRunNumber);

    if(newID != lastSubRunID_) {
      std::cout << "Subrun " << subRunNumber << std::endl;
      // art takes ownership of the object pointed to by outSR and will delete it at the appropriate time.
      outSR = pm_.makeSubRunPrincipal(runNumber,
          subRunNumber,
          ts);

    }
    lastSubRunID_ = newID;

    // art takes ownership of the object pointed to by outE and will delete it at the appropriate time.
    outE = pm_.makeEventPrincipal(runNumber, subRunNumber, eventNumber, ts, false);

  } // managePrincipals()
  //----------------------------------------------------------------

} // namespace mu2e

typedef art::Source<mu2e::TrackerPrototypeDataDetail> FromTrackerPrototypeData;
DEFINE_ART_INPUT_SOURCE(FromTrackerPrototypeData)
