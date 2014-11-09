#ifndef ROOTtools_eventCount_FileInfo_hh
#define ROOTtools_eventCount_FileInfo_hh
//
// Find number of runs, subruns and events in an art format ROOT file.
//
// See the implementation of usage() for the documentation.
//

#include <iosfwd>
#include <string>
#include <vector>

namespace mu2e {

  struct FileInfo{

    FileInfo ( std::string const& name );

    std::string   filename;
    unsigned long events     = 0;
    unsigned long runs       = 0;
    unsigned long subRuns    = 0;
    bool          openable   = false;
    bool          hasEvents  = false;
    bool          hasRuns    = false;
    bool          hasSubRuns = false;

    void fullPrint    ( std::ostream& ) const;
    void minimalPrint ( std::ostream& ) const;

    bool allOK () const {
      return openable && hasEvents && hasRuns && hasSubRuns;
    }

    std::string status() const{
      return (allOK()) ? "OK " : "BAD";
    }
  };

}

#endif /* ROOTtools_eventCount_FileInfo_hh */
