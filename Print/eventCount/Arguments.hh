#ifndef ROOTtools_eventCount_Arguments_hh
#define ROOTtools_eventCount_Arguments_hh
//
// Parse and validate the argument list for the artProductSizes executable.
//
// See the implementation of usage() for the documentation.
//

#include <string>
#include <vector>

namespace mu2e {

  struct Arguments{

    // Style of printout.
    enum PrintStyle { minimal, full, events, subruns, sam};

    Arguments ( int argc, char**argv );

    std::vector<std::string> fileNames;

    PrintStyle style = minimal;

  private:
    void usage() const;

  };

}

#endif /* ROOTtools_eventCount_Arguments_hh */
