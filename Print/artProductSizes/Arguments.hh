#ifndef ROOTtools_artProductSizes_Arguments_hh
#define ROOTtools_artProductSizes_Arguments_hh
//
// Parse and validate the argument list for the artProductSizes executable.
//
// See the implementation of usage() for the documentation.
//

#include <iostream>
#include <string>
#include <vector>

namespace mu2e {

  struct Arguments{

    Arguments ( int argc, char**argv );

    double minimumFractionDefault;
    double minimumFraction;
    std::vector<std::string> fileNames;

  private:
    double toDouble( char const*, double low, double high);
    void usage() const;

  };

}

#endif /* ROOTtools_artProductSizes_Arguments_hh */
