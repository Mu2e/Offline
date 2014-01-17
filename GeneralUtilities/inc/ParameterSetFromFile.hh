#ifndef GeneralUtilities_ParameterSetFromFile_h
#define GeneralUtilities_ParameterSetFromFile_h

//
// Given a filename, look for that file in FHICL_FILE_PATH.
// If found, open it, interpret it as a FHiCL file and turn
// it into a fhicl::ParameterSet object.
//

#include "fhiclcpp/ParameterSet.h"

#include <string>
#include <iosfwd>

namespace mu2e {

  class ParameterSetFromFile {

public:
    ParameterSetFromFile( std::string const& fileName );

    fhicl::ParameterSet const& pSet() const { return _pSet; }

    void printKeys( std::ostream& ) const;

private:

    std::string         _fileName;
    fhicl::ParameterSet _pSet;

  };

}

#endif /* GeneralUtilities_ParameterSetFromFile_h */
