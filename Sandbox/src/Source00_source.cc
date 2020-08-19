// A very simple source module.
//
// Original author: Rob Kutschke, 2014
// Adapted from Sources/src/FromEMFMARSFileWeighted_source.cc
//

#include <iostream>
#include <fstream>
#include <boost/utility.hpp>

#include "fhiclcpp/types/Name.h"
#include "art/Framework/IO/Sources/Source.h"
#include "art/Framework/Core/InputSourceMacros.h"
#include "art/Framework/IO/Sources/SourceHelper.h"
#include "canvas/Persistency/Provenance/Timestamp.h"
#include "canvas/Persistency/Provenance/RunID.h"
#include "canvas/Persistency/Provenance/SubRunID.h"
#include "canvas/Persistency/Provenance/EventID.h"
#include "art/Framework/Principal/RunPrincipal.h"
#include "art/Framework/Principal/SubRunPrincipal.h"
#include "art/Framework/Principal/EventPrincipal.h"
#include "art/Framework/IO/Sources/put_product_in_principal.h"

namespace mu2e {

  class Source00Detail : private boost::noncopyable {

  public:

    // The signatures and return types of the c'tor and the three member functions
    // are prescribed by the art::Source class template.
    Source00Detail(const fhicl::ParameterSet&,
                         art::ProductRegistryHelper&,
                   const art::SourceHelper&);

    void readFile(std::string const& filename,
                  art::FileBlock*& fb);

    bool readNext(art::RunPrincipal*    const& inR,
                  art::SubRunPrincipal* const& inSR,
                  art::RunPrincipal*&          outR,
                  art::SubRunPrincipal*&       outSR,
                  art::EventPrincipal*&        outE);

    void closeCurrentFile();

  private:

    // Control the level of printout.
    int verbosity_;

    // The data products will appear in the art::Event with the following module label:
    std::string myModuleLabel_;

    // Helper object needed by the member function managePrincipals.
    art::SourceHelper const& pm_;

    // Used to identify boundaries between runs and subruns.
    art::SubRunID lastSubRunID_;

    // The currently active file.
    std::string currentFileName_;
    std::ifstream currentFile_;

    // A helper function used to manage the principals.
    // This is boilerplate that does not change if you change the data products.
    void managePrincipals ( int runNumber,
                            int subRunNumber,
                            int eventNumber,
                            art::RunPrincipal*&    outR,
                            art::SubRunPrincipal*& outSR,
                            art::EventPrincipal*&  outE);

  };

  Source00Detail::Source00Detail(const fhicl::ParameterSet&        pset,
                                       art::ProductRegistryHelper& rh,
                                 const art::SourceHelper&        pm)
    : verbosity_( pset.get<int>("detail.verbosity",0))
    , myModuleLabel_("Source00")
    , pm_(pm)
    , lastSubRunID_()
    , currentFileName_()
    , currentFile_()
  {
    // This is the analog of the call to produces<T> in a producer module.
    // Same pattern is used for art::InRun and art::InSubRun objects.
    rh.reconstitutes<int,art::InEvent>(myModuleLabel_);

    if ( verbosity_ > 0 ){
      std::cout << "Source00Detail: constructor" << std::endl;
    }
  }

  //----------------------------------------------------------------
  void Source00Detail::readFile(const std::string& filename, art::FileBlock*& fb) {

    if ( verbosity_ > 1 ){
      std::cout << "Source00Detail: open file " << filename << std::endl;
    }

    // Open the input file.
    currentFileName_    = filename;
    currentFile_.open(filename.c_str());

    // A FileBlock is the object art uses to maintain state information about the file.
    // art takes ownership of this pointer and will call delete on it at the appropriate time.
    fb = new art::FileBlock(art::FileFormatVersion(1, "Source00"), currentFileName_);

  }

  void Source00Detail::closeCurrentFile() {

    if ( verbosity_ > 1 ){
      std::cout << "Source00Detail: close file " << currentFileName_ << std::endl;
    }

    currentFileName_ = "";
    currentFile_.close();
  }

  bool Source00Detail::readNext(art::RunPrincipal*    const& inR,
                                art::SubRunPrincipal* const& inSR,
                                art::RunPrincipal*&          outR,
                                art::SubRunPrincipal*&       outSR,
                                art::EventPrincipal*&        outE)
  {

    // Read one line from the input file.
    int runNumber, subRunNumber, eventNumber, magic;
    currentFile_ >> runNumber >> subRunNumber >> eventNumber >> magic;

    if ( verbosity_ > 2 ){
      std::cout << "Source00Detail: readNext: "
                << runNumber    << " "
                << subRunNumber << " "
                << eventNumber  << " "
                << magic
                << std::endl;
    }

    // Tell art that we have reached of this file; they may be more files to come.
    if ( !currentFile_ ) return false;

    // Create principals as needed.
    managePrincipals( runNumber, subRunNumber, eventNumber, outR, outSR, outE );

    // Create an empty data product and then fill it with information read from the input file.
    std::unique_ptr<int> product(new int);
    *product = magic;

    // This is the analog of event.put in a producer module.
    art::put_product_in_principal(std::move(product), *outE, myModuleLabel_);

    // Tell art that this is NOT the end of file.
    return true;

  } // readNext()

  // Each time that we encounter a new run, a new subRun or a new event, we need to make a new principal
  // of the appropriate type.  This code does not need to change as the number and type of data products changes.
  void Source00Detail::managePrincipals ( int runNumber,
                                          int subRunNumber,
                                          int eventNumber,
                                          art::RunPrincipal*&    outR,
                                          art::SubRunPrincipal*& outSR,
                                          art::EventPrincipal*&  outE){

    art::Timestamp ts;

    art::SubRunID newID(runNumber, subRunNumber);
    if(newID.runID() != lastSubRunID_.runID()) {
      // art takes ownership of the object pointed to by outR and will delete it at the appropriate time.
      outR = pm_.makeRunPrincipal(runNumber, ts);
      if ( verbosity_ > 1 ){
        std::cout << "Source00Detail: making run principal: "
                  << newID.runID() << " "
                  << lastSubRunID_.runID()
                  << std::endl;
      }
    }

    if(newID != lastSubRunID_) {
      // art takes ownership of the object pointed to by outSR and will delete it at the appropriate time.
      outSR = pm_.makeSubRunPrincipal(runNumber,
                                      subRunNumber,
                                      ts);
      if ( verbosity_ > 1 ){
        std::cout << "Source00Detail: making subRun principal: "
                  << newID << " "
                  << lastSubRunID_
                  << std::endl;
      }
    }
    lastSubRunID_ = newID;

    // art takes ownership of the object pointed to by outE and will delete it at the appropriate time.
    outE = pm_.makeEventPrincipal(runNumber, subRunNumber, eventNumber, ts, false);

  } // managePrincipals()

} // namespace mu2e

typedef art::Source<mu2e::Source00Detail> Source00;
DEFINE_ART_INPUT_SOURCE(Source00);
