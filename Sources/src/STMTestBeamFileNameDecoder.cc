//
// Functions for FileNameDecoder
//

#include "Offline/Sources/inc/STMTestBeamFileNameDecoder.hh"
#include "canvas/Persistency/Provenance/RunID.h"
#include "cetlib_except/exception.h"

namespace mu2e {

  namespace STMTestBeam {
    FileNameDecoder::FileNameDecoder(std::string filename) {
      // Extract the tokens from the filename
      std::string delimeter = ".";
      size_t currentPos = 0;
      size_t previousPos = currentPos;
      while ( (currentPos = filename.find(delimeter, previousPos)) != std::string::npos) {
        _tokens.push_back(filename.substr(previousPos, currentPos-previousPos));
        previousPos = currentPos+1;
      }
      _tokens.push_back(filename.substr(previousPos, currentPos-previousPos)); // add the final token which gets missed in the above loop

      unsigned int n_expected_fields = 6;
      if (_tokens.size() != n_expected_fields) {
        throw cet::exception("FromSTMTestBeamData") << "Number of fields in filename (" << _tokens.size() << ") is not the same number we expected (" << n_expected_fields << "). Filename = " << filename << std::endl;
      }
    }

    mu2e::STMChannel FileNameDecoder::extractSTMChannel() {
      // Get the STM Channel from the configuration field
      std::string configuration = _tokens.at(3);
      if (configuration.find("HPGe") != std::string::npos) {
        return mu2e::STMChannel::HPGe;
      }
      else if (configuration.find("LaBr") != std::string::npos) {
        return mu2e::STMChannel::LaBr;
      }
      else {
        throw cet::exception("FromSTMTestBeamData") << "Cannot determine the channel from the configuration field (" << configuration << "). This should contain either \"HPGe\" or \"LaBr\"" << std::endl;
      }
    }

    unsigned FileNameDecoder::extractRunNumber() {
      // Get the run and subrun numbers from the sequencer
      std::string sequencer = _tokens.at(4);
      std::string runNo = sequencer.substr(0, 6); // sequencer always has 6 characters for run
      unsigned runNumber = std::stoi(runNo);

      if(!art::RunID(runNumber).isValid()) {
        throw cet::exception("BADCONFIG", " FromSTMTestBeamData: ")
          << " fhicl::ParameterSet specifies an invalid runNumber = "<<runNumber<<"\n";
      }

      // We need to restrict the number to our reserved run range
      // as of writing, this is defined in the Mu2e wiki
      if(runNumber < 101000 || runNumber > 101999) {
        throw cet::exception("FromSTMTestBeamData") << "Run number (" << runNumber << ") is outside of our reserved run range (101000 -- 101999)" << std::endl;
      }
      return runNumber;
    }

    unsigned FileNameDecoder::extractSubRunNumber() {
      // Get the run and subrun numbers from the sequencer
      std::string sequencer = _tokens.at(4);
      std::string subrunNo = sequencer.substr(7, 7+8); // sequencer always has 6 characters for run followed by an underscore and then 8 characters for the sub run
      return std::stoi(subrunNo);
    }

    int FileNameDecoder::extractBinaryFileVersion(unsigned runNumber) {
      // There are two different versions of the binary file:
      //  - v2 contains a unix timestamp after the trigger header
      //  - v1 does not
      // For the time being, just hardcode which runs belong to which (ideally would be in a DB?)
      if (runNumber == 101001) {
        return 2;
      }
      else if (runNumber >= 101002 && runNumber <= 101014) {
        return 1;
      }
      else {
        return 2;
      }
    }
  }
}
