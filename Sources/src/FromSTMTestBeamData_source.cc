//
// Source module to read binary files from the STM test beam
// at gELBE collected in April 2022
//
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

#include "Offline/Sources/inc/STMTestBeamHeaders.hh"
#include "Offline/Sources/inc/STMTestBeamFileNameDecoder.hh"
#include "Offline/DataProducts/inc/STMTestBeamEventInfo.hh"
#include "Offline/RecoDataProducts/inc/STMWaveformDigi.hh"


using namespace std;

namespace mu2e {
  struct Config
  {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    fhicl::Sequence<std::string> inputFiles{Name("fileNames"),Comment("Input binary file")};
    fhicl::Atom<unsigned int> maxEvents{Name("maxEvents"), Comment("Max number of events")};
    fhicl::Atom<int> verbosityLevel{Name("verbosityLevel"), Comment("Verbosity level")};

    // These are used by art and are required.
    fhicl::Atom<std::string> module_label{Name("module_label"), Comment("Art module label"), ""};
    fhicl::Atom<std::string> module_type{Name("module_type"), Comment("Art module type"), ""};
  };
  typedef fhicl::WrappedTable<Config> Parameters;


  //================================================================
  class STMTestBeamDataDetail : private boost::noncopyable {
    std::string myModuleType_;
    art::SourceHelper const& pm_;
    unsigned runNumber_ = 0; // from file name
    art::SubRunID lastSubRunID_;
    std::set<art::SubRunID> seenSRIDs_;

    unsigned currentFileNumber_ = 0;
    std::string currentFileName_ = "";
    std::unique_ptr<std::ifstream> currentFile_ = nullptr;

    unsigned currentSubRunNumber_ = -1U; // from file
    unsigned currentEventNumber_ = 0;
    unsigned maxEvents_ = 0;
    int verbosityLevel_ = 0;
    mu2e::STMChannel channel_; // the channel (HPGe or LaBr) will be different for each file
    int binaryFileVersion_ = 0;

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
    : myModuleType_(conf().module_type())
    , pm_(pm)
    , maxEvents_(conf().maxEvents())
    , verbosityLevel_(conf().verbosityLevel())
  {
    rh.reconstitutes<mu2e::STMWaveformDigiCollection,art::InEvent>(myModuleType_, "HPGe");
    rh.reconstitutes<mu2e::STMWaveformDigiCollection,art::InEvent>(myModuleType_, "LaBr");
    rh.reconstitutes<mu2e::STMTestBeamEventInfo,art::InEvent>(myModuleType_);
  }


  //----------------------------------------------------------------
  void STMTestBeamDataDetail::readFile(const std::string& filename, art::FileBlock*& fb) {

    // Open the input file
    currentFileName_ = filename;
    currentFile_ = std::make_unique<std::ifstream>(currentFileName_.c_str(), std::ios::in | std::ios::binary);
    if (!currentFile_->is_open()) {
      throw cet::exception("FromSTMTestBeamData") << "A problem opening binary file " << currentFileName_ << std::endl;
    }

    mu2e::STMTestBeam::FileNameDecoder decoder(currentFileName_);
    channel_ = decoder.extractSTMChannel();
    runNumber_ = decoder.extractRunNumber();
    binaryFileVersion_ = decoder.extractBinaryFileVersion(runNumber_);
    currentSubRunNumber_ = decoder.extractSubRunNumber();
    currentEventNumber_ = 1;

    fb = new art::FileBlock(art::FileFormatVersion(1, "STMTestBeamDataInput"), currentFileName_); // art takes ownership
  }

  //----------------------------------------------------------------
  void STMTestBeamDataDetail::closeCurrentFile() {
    currentFileName_ = "";
    currentFile_.reset();
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

    // Create the STMWaveformDigiCollection that we will write to the art event
    // NB although we are creating one collection for each detector, in the test beam only one
    //    detector was connected at a time so one of these collections will be empty
    std::unique_ptr<mu2e::STMWaveformDigiCollection> outputHPGeWaveformDigis(new mu2e::STMWaveformDigiCollection());
    std::unique_ptr<mu2e::STMWaveformDigiCollection> outputLaBrWaveformDigis(new mu2e::STMWaveformDigiCollection());

    // Read the trigger header
    STMTestBeam::TriggerHeader trigger_header;
    if(currentFile_->read(reinterpret_cast<char *>(&trigger_header), sizeof(STMTestBeam::TriggerHeader))) {
      managePrincipals(runNumber_, currentSubRunNumber_, currentEventNumber_, outR, outSR, outE);

      if(verbosityLevel_ > 0) {
        std::cout << trigger_header << std::endl;
      }
      if (!trigger_header.checkFixedHeader()) {
        throw cet::exception("FromSTMTestBeamData") << "Fixed header word (0x" << std::setfill('0') << std::setw(8) << std::hex << (trigger_header.getFixedHeader()) << ") != 0xDEADBEEF" << std::endl;
      }

      // binary file version 2 now contains the unix time stamp
      if (binaryFileVersion_ == 2) {
        uint16_t unixtime[4];
        currentFile_->read(reinterpret_cast<char *>(&unixtime[0]), sizeof(unixtime));
        uint64_t time = (uint64_t(unixtime[3]) << 48) | (uint64_t(unixtime[2]) << 32) | (uint64_t(unixtime[1]) << 16) | (uint64_t(unixtime[0]));
        if (verbosityLevel_ > 0) {
          using time_point = std::chrono::system_clock::time_point;
          time_point header_timepoint(std::chrono::duration_cast<time_point::duration>(std::chrono::milliseconds(time)));
          std::time_t header_t = std::chrono::system_clock::to_time_t(header_timepoint);
          std::cout << "Unix Time: " << std::put_time(std::gmtime(&header_t), "%c %Z") << std::endl;
        }
      }

      // Get event information for this trigger and put in the event
      std::unique_ptr<mu2e::STMTestBeamEventInfo> outputEvtInfo(new mu2e::STMTestBeamEventInfo(trigger_header.getTriggerMode(), trigger_header.getTriggerTime()));
      art::put_product_in_principal(std::move(outputEvtInfo), *outE, myModuleType_);

      // Get the number of slices in this trigger
      int n_slices = trigger_header.getNSlices();
      for (int i_slice = 0; i_slice < n_slices; ++i_slice) {
        STMTestBeam::SliceHeader slice_header[1];
        // Read the slice header
        if(currentFile_->read(reinterpret_cast<char *>(&slice_header[0]), sizeof(STMTestBeam::SliceHeader))) {
          if(verbosityLevel_ > 0) {
            std::cout << slice_header[0] << std::endl;
          }

          // Get the number of ADC samples in this slice
          unsigned long int n_adc_samples = slice_header[0].getNADC();

          // Read the ADC samples into a vector for later
          std::vector<int16_t> adcs;
          adcs.reserve(n_adc_samples);
          int16_t adc[1];
          for (unsigned long int i_adc_sample = 0; i_adc_sample < n_adc_samples; ++i_adc_sample) {
            currentFile_->read(reinterpret_cast<char *>(&adc[0]), sizeof(int16_t));
            adcs.push_back(adc[0]);
          }

          // Create the STMWaveformDigi and put it in the event
          STMWaveformDigi stm_waveform(trigger_header.getTriggerOffset(), adcs);
          if (channel_ == STMChannel::HPGe) {
            outputHPGeWaveformDigis->push_back(stm_waveform);
          }
          else if (channel_ == STMChannel::LaBr) {
            outputLaBrWaveformDigis->push_back(stm_waveform);
          }
          else {
            throw cet::exception("FromSTMTestBeamData") << "Trying to create a waveform with an invalid STMChannel (" << channel_ << ")" << std::endl;
          }
        }
        else { return false; }
      }
    }
    else { return false; }

    art::put_product_in_principal(std::move(outputHPGeWaveformDigis), *outE, myModuleType_, "HPGe");
    art::put_product_in_principal(std::move(outputLaBrWaveformDigis), *outE, myModuleType_, "LaBr");

    ++currentEventNumber_;

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

typedef art::Source<mu2e::STMTestBeamDataDetail> FromSTMTestBeamData;
DEFINE_ART_INPUT_SOURCE(FromSTMTestBeamData)
