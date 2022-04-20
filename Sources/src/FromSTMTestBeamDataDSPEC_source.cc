#include <iostream>
#include <fstream>
#include <boost/utility.hpp>
#include <cassert>
#include <set>
#include <string>
#include <vector>

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
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "Offline/RecoDataProducts/inc/STMTestBeamPacketDefns.hh"
#include "Offline/RecoDataProducts/inc/STMDigi.hh"


using namespace std;

namespace mu2e {
  struct Config
  {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    fhicl::Sequence<std::string> inputFiles{Name("fileNames"),Comment("Input binary file")};
    fhicl::Atom<unsigned int> runNumber{Name("runNumber"), Comment("Run number")};
    fhicl::Atom<unsigned int> maxEvents{Name("maxEvents"), Comment("Max number of events")};

    // These are used by art and are required.
    fhicl::Atom<std::string> module_label{Name("module_label"), Comment("Art module label"), ""};
    fhicl::Atom<std::string> module_type{Name("module_type"), Comment("Art module type"), ""};
  };
  typedef fhicl::WrappedTable<Config> Parameters;


  //================================================================
  class STMTestBeamDataDSPECDetail : private boost::noncopyable {
    std::string myModuleLabel_;
    art::SourceHelper const& pm_;
    unsigned runNumber_; // from ParSet
    art::SubRunID lastSubRunID_;
    std::set<art::SubRunID> seenSRIDs_;
    
    unsigned currentFileNumber_;
    std::string currentFileName_;
    std::ifstream* currentFile_ = nullptr;
      
    mu2e::STMDigiCollection digis_;

    unsigned currentSubRunNumber_;
    unsigned currentEventNumber_;
    unsigned maxEvents_;
    unsigned nBins_; // from file

    int printAtEvent;

    // A helper function used to manage the principals.
    // This is boilerplate that does not change if you change the data products.
    void managePrincipals ( int runNumber,
        int subRunNumber,
        int eventNumber,
        art::RunPrincipal*&    outR,
        art::SubRunPrincipal*& outSR,
        art::EventPrincipal*&  outE);

    public:
    STMTestBeamDataDSPECDetail(const Parameters &conf,
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
  STMTestBeamDataDSPECDetail::STMTestBeamDataDSPECDetail(const Parameters& conf,
      art::ProductRegistryHelper& rh,
      const art::SourceHelper& pm)
    : myModuleLabel_("FromSTMTestBeamDataDSPEC")
      , pm_(pm)
      , runNumber_(conf().runNumber())
    , currentFileNumber_(0)
      , currentSubRunNumber_(-1U)
      , currentEventNumber_(0)
      , maxEvents_(conf().maxEvents())
    , nBins_(0)
  {
    if(!art::RunID(runNumber_).isValid()) {
      throw cet::exception("BADCONFIG", " FromSTMTestBeamDataDSPEC: ")
        << " fhicl::ParameterSet specifies an invalid runNumber = "<<runNumber_<<"\n";
    }

    rh.reconstitutes<mu2e::STMDigiCollection,art::InEvent>(myModuleLabel_);

    currentSubRunNumber_ = 0;
    currentEventNumber_ = 0;

    printAtEvent = 0;

  }


  //----------------------------------------------------------------
  void STMTestBeamDataDSPECDetail::readFile(const std::string& filename, art::FileBlock*& fb) {

    currentFileName_ = filename;

    currentFile_ = new std::ifstream(currentFileName_.c_str(), std::ios::in);
    if (!currentFile_->is_open()) {
      throw cet::exception("MakeSTMDigisFromBin") << "A problem opening DSPEC file " << currentFileName_ << std::endl;
    }
    // Read through files until we get to the "$DATA:" keywoerd
    std::string line;
    bool data_found = false;
    while (!data_found) {
      std::getline(*currentFile_, line);
      if (line.find("$DATA:") != std::string::npos) {
	data_found = true;
	break;
      }
    }
    // next line gives the number of data lines
    std::getline(*currentFile_, line, ' '); // first word is 0
    std::getline(*currentFile_, line, ' '); // second word is number of data lines
    nBins_ = std::stoi(line);

    ++currentSubRunNumber_;
    currentEventNumber_ = 1;
    fb = new art::FileBlock(art::FileFormatVersion(1, "STMTestBeamDataDSPECInput"), currentFileName_);
  }

  //----------------------------------------------------------------
  void STMTestBeamDataDSPECDetail::closeCurrentFile() {
    currentFileName_ = "";
    currentFile_->close();
    delete currentFile_;
    currentFile_ = nullptr;
  }

  //----------------------------------------------------------------
  bool STMTestBeamDataDSPECDetail::readNext(art::RunPrincipal* const& inR,
      art::SubRunPrincipal* const& inSR,
      art::RunPrincipal*& outR,
      art::SubRunPrincipal*& outSR,
      art::EventPrincipal*& outE)
  {
    if (currentFile_->eof()) {
      return false;
    }
    managePrincipals(runNumber_, currentSubRunNumber_, currentEventNumber_, outR, outSR, outE);
    std::unique_ptr<mu2e::STMDigiCollection> outputSTMDigis(new mu2e::STMDigiCollection);

    std::string line;
    for (unsigned i_bin = 0; i_bin < nBins_; ++i_bin) {
      std::getline(*currentFile_, line);
      int n_counts = std::stoi(line);
      std::vector<int16_t> adcs;
      adcs.push_back(i_bin); // just a single ADC here
      for (int i_count = 0; i_count < n_counts; ++i_count) {
	STMDigi stm_digi(0, 0, 0, 0, 0, 0, adcs);
	outputSTMDigis->push_back(stm_digi);
      }
    }
    std::getline(*currentFile_, line);

    art::put_product_in_principal(std::move(outputSTMDigis), *outE, myModuleLabel_);

    ++currentEventNumber_;

    return true;
  } // readNext()


  // Each time that we encounter a new run, a new subRun or a new event, we need to make a new principal
  // of the appropriate type.  This code does not need to change as the number and type of data products changes.
  void STMTestBeamDataDSPECDetail::managePrincipals ( int runNumber,
      int subRunNumber,
      int eventNumber,
      art::RunPrincipal*&    outR,
      art::SubRunPrincipal*& outSR,
      art::EventPrincipal*&  outE){

    art::Timestamp ts;

    //    if (eventNumber == printAtEvent){
    //      printAtEvent = (printAtEvent+1)*2-1;
    //    }
    std::cout << "AE: run, subrun, event = " << runNumber << ", " << subRunNumber << ", " << eventNumber << std::endl;

    art::SubRunID newID(runNumber, subRunNumber);

    if(newID != lastSubRunID_) {
      outR = pm_.makeRunPrincipal(runNumber, ts);
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

typedef art::Source<mu2e::STMTestBeamDataDSPECDetail> FromSTMTestBeamDataDSPEC;
DEFINE_ART_INPUT_SOURCE(FromSTMTestBeamDataDSPEC);
