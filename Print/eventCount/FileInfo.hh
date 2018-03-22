#ifndef ROOTtools_eventCount_FileInfo_hh
#define ROOTtools_eventCount_FileInfo_hh
//
// Find number of runs, subruns and events in an art format ROOT file.
//
// See the implementation of usage() for the documentation.
//
//#include "canvas/Persistency/Provenance/SubRunID.h"
#include "canvas/Persistency/Provenance/EventID.h"
//#include "art/Framework/Principal/Event.h"
//#include "art/Framework/Principal/Run.h"
//#include "art/Framework/Principal/SubRun.h"

#include "TFile.h"
#include <iosfwd>
#include <string>
#include <vector>

namespace mu2e {

  struct FileInfo{

    FileInfo ( std::string const& name, int level=1 );
 
    std::string   filename;
    unsigned long events     = 0;
    unsigned long runs       = 0;
    unsigned long subRuns    = 0;
    bool          openable   = false;
    bool          hasEvents  = false;
    bool          hasRuns    = false;
    bool          hasSubRuns = false;

    std::vector<art::EventID>  eventList;
    std::vector<art::SubRunID> subrunList;

    void fullPrint    ( std::ostream& ) const;
    void minimalPrint ( std::ostream& ) const;
    void eventPrint   ( std::ostream& ) const;
    void subrunPrint  ( std::ostream& ) const;
    void samPrint     ( std::ostream& ) const;

    bool allOK () const {
      return openable && hasEvents && hasRuns && hasSubRuns;
    }

    std::string status() const{
      return (allOK()) ? "OK " : "BAD";
    }

    void treeInfo ( std::string const& treeName,   // input argument
		    TFile* file,                   // input argument
		    bool&              exists,     // output argument
		    unsigned long&     nEntries    // output argument
		    );

    void makeVectors(TFile* file);
  };

}

#endif /* ROOTtools_eventCount_FileInfo_hh */
