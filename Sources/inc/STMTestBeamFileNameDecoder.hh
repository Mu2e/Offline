#ifndef Sources_STMTestBeamFileNameDecoder_hh
#define Sources_STMTestBeamFileNameDecoder_hh
//
// Class to decode the STM test beam data file names
//

#include "Offline/DataProducts/inc/STMChannel.hh"

namespace mu2e {

  namespace STMTestBeam {
    class FileNameDecoder {
    public:
      FileNameDecoder(std::string filename);

      mu2e::STMChannel extractSTMChannel();
      unsigned extractRunNumber();
      unsigned extractSubRunNumber();
      int extractBinaryFileVersion(unsigned runNumber);

    private:
      std::vector<std::string> _tokens;
    };

  }
}
#endif
