#include <iostream>
#include <fstream>
#include <boost/utility.hpp>
#include <cassert>
#include <set>
#include <string>

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

#include "Offline/RecoDataProducts/inc/STMTestBeamBinaryPacket.hh"
#include "Offline/RecoDataProducts/inc/STMDigi.hh"


using namespace std;

namespace mu2e {
  struct Config
  {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    fhicl::Sequence<std::string> inputFiles{Name("fileNames"),Comment("Input binary file")};
    fhicl::Atom<unsigned int> runNumber{Name("runNumber"), Comment("Run number")};
    fhicl::Atom<unsigned int> maxEvents{Name("maxEvents"), Comment("Max number of events"), 0};

    // These are used by art and are required.
    fhicl::Atom<std::string> module_label{Name("module_label"), Comment("Art module label"), ""};
    fhicl::Atom<std::string> module_type{Name("module_type"), Comment("Art module type"), ""};
  };
  typedef fhicl::WrappedTable<Config> Parameters;


  //================================================================
  class STMTestBeamDataDetail : private boost::noncopyable {
    std::string myModuleLabel_;
    art::SourceHelper const& pm_;
    unsigned runNumber_; // from ParSet
    art::SubRunID lastSubRunID_;
    std::set<art::SubRunID> seenSRIDs_;
    
    unsigned currentFileNumber_;
    std::string currentFileName_;
    std::ifstream* currentFile_ = nullptr;
      
    mu2e::STMDigiCollection digis_;

    unsigned currentSubRunNumber_; // from file
    unsigned currentEventNumber_;
    unsigned maxEvents_;

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
    STMTestBeamDataDetail(const Parameters &conf,
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
  STMTestBeamDataDetail::STMTestBeamDataDetail(const Parameters& conf,
      art::ProductRegistryHelper& rh,
      const art::SourceHelper& pm)
    : myModuleLabel_("FromSTMTestBeamData")
      , pm_(pm)
      , runNumber_(conf().runNumber())
    , currentFileNumber_(0)
      , currentSubRunNumber_(-1U)
      , currentEventNumber_(0)
      , maxEvents_(conf().maxEvents())
  {
    if(!art::RunID(runNumber_).isValid()) {
      throw cet::exception("BADCONFIG", " FromSTMTestBeamData: ")
        << " fhicl::ParameterSet specifies an invalid runNumber = "<<runNumber_<<"\n";
    }

    rh.reconstitutes<mu2e::STMDigiCollection,art::InEvent>(myModuleLabel_);

    currentSubRunNumber_ = 0;
    currentEventNumber_ = 0;

    printAtEvent = 0;

  }


  //----------------------------------------------------------------
  void STMTestBeamDataDetail::readFile(const std::string& filename, art::FileBlock*& fb) {

    currentFileName_ = filename;

    currentFile_ = new std::ifstream(currentFileName_.c_str(), std::ios::in | std::ios::binary);
    if (!currentFile_->is_open()) {
      throw cet::exception("MakeSTMDigisFromBin") << "A problem opening binary file " << currentFileName_ << std::endl;
    }
    ++currentSubRunNumber_;
    currentEventNumber_ = 1;
    fb = new art::FileBlock(art::FileFormatVersion(1, "STMTestBeamDataInput"), currentFileName_);
  }

  //----------------------------------------------------------------
  void STMTestBeamDataDetail::closeCurrentFile() {
    currentFileName_ = "";
    currentFile_->close();
    delete currentFile_;
    currentFile_ = nullptr;
  }

  //----------------------------------------------------------------
  bool STMTestBeamDataDetail::readNext(art::RunPrincipal* const& inR,
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
    STMTestBeamBinaryPacket stmTestBeamPacket[1];
    while (currentFile_->read((char *) &stmTestBeamPacket[0], sizeof(STMTestBeamBinaryPacket))) {
      std::cout << stmTestBeamPacket[0] << std::endl;
      //      if (!currentFile_) {
      //	throw cet::exception("MakeSTMDigisFromBin") << "A problem reading binary file " << currentFileName_ << std::endl;
      //      }
      STMDigi stm_digi(stmTestBeamPacket[0].trigTimeOffset, stmTestBeamPacket[0].ADC0);
      outputSTMDigis->push_back(stm_digi);
    }
    art::put_product_in_principal(std::move(outputSTMDigis), *outE, myModuleLabel_);

    return true;
  } // readNext()


  // Each time that we encounter a new run, a new subRun or a new event, we need to make a new principal
  // of the appropriate type.  This code does not need to change as the number and type of data products changes.
  void STMTestBeamDataDetail::managePrincipals ( int runNumber,
      int subRunNumber,
      int eventNumber,
      art::RunPrincipal*&    outR,
      art::SubRunPrincipal*& outSR,
      art::EventPrincipal*&  outE){

    art::Timestamp ts;

    if (eventNumber == printAtEvent){
      std::cout << "Event " << eventNumber << std::endl;
      printAtEvent = (printAtEvent+1)*2-1;
    }

    art::SubRunID newID(runNumber, subRunNumber);

    if(newID != lastSubRunID_) {
      outR = pm_.makeRunPrincipal(runNumber, ts);
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

typedef art::Source<mu2e::STMTestBeamDataDetail> FromSTMTestBeamData;
DEFINE_ART_INPUT_SOURCE(FromSTMTestBeamData);
